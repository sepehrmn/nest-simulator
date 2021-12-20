/*
 *  matco_synapse.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef MATCO_SYNAPSE_H
#define MATCO_SYNAPSE_H

// Includes from nestkernel:
#include "connection.h"

namespace nest
{

/* BeginUserDocs: synapse, plasticity

Short description
+++++++++++++++++

Synapse type for LTD/LTP based on firing rates

Description
+++++++++++

matco_synapse is a connection to create synapses with basic plasticity
following [1] and [2].

.. warning::

   This synaptic plasticity rule does not take
   :doc:`precise spike timing <simulations_with_precise_spike_times>` into
   account. When calculating the weight update, the precise spike time part
   of the timestamp is ignored.

Parameters
++++++++++

=== ======  =====================================
**Individual properties**
-------------------------------------------------
theta      real    Threshold
=== ======  =====================================

Remarks:

The common properties can only be set by SetDefaults and apply to all
synapses of the model.

References
++++++++++

.. [1] Tomasello et. al. 2018 
.. [2] 

Transmits
+++++++++

SpikeEvent

See also
++++++++

iaf_matco_2018

EndUserDocs */

template < typename targetidentifierT >
class matco_synapse : public Connection< targetidentifierT >
{
public:
  typedef CommonSynapseProperties CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  matco_synapse();

  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  matco_synapse( const matco_synapse& ) = default;

  // Explicitly declare all methods inherited from the dependent base
  // ConnectionBase. This avoids explicit name prefixes in all places these
  // functions are used. Since ConnectionBase depends on the template parameter,
  // they are not automatically found in the base class.
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_delay;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;

  /**
   * Default Destructor.
   */
  virtual ~matco_synapse()
  {
  }

  /**
   * Get all properties of this connection and put them into a dictionary.
   */
  virtual void get_status( DictionaryDatum& d ) const;

  /**
   * Set properties of this connection from the values given in dictionary.
   */
  virtual void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param cp Common properties to all synapses (empty).
   */
  void send( Event& e, thread t, const CommonSynapseProperties& cp );

  class ConnTestDummyNode : public ConnTestDummyNodeBase
  {
  public:
    // Ensure proper overriding of overloaded virtual functions.
    // Return values from functions are ignored.
    using ConnTestDummyNodeBase::handles_test_event;
    port
    handles_test_event( SpikeEvent&, rport )
    {
      return invalid_port_;
    }
  };

  void
  check_connection( Node& s, Node& t, rport receptor_type, const CommonPropertiesType& )
  {
    ConnTestDummyNode dummy_target;
    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );
  }

  //! allows efficient initialization from ConnectorModel::add_connection()
  void
  set_weight( double w )
  {
    weight_ = w;
  }

private:
  double weight_; //!< Synaptic weight

  double t_lastspike_; //!< Time point of last spike emitted

  double omega_E; 

  std::vector< double > firing_rates_;
  std::vector< double > membrane_potentials_;
  std::vector< double > deltas_;
  std::vector< double > weights_;
  std::vector< double > plasticity_flags_;

};


/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param p The port under which this connection is stored in the Connector.
 */
template < typename targetidentifierT >
inline void
matco_synapse< targetidentifierT >::send( Event& e, thread t, const CommonSynapseProperties& )
{
  // propagation t_lastspike -> t_spike, t_lastspike_ = 0 initially, p_ = 1
  const double t_spike = e.get_stamp().get_ms();
  const double h = t_spike - t_lastspike_;

  // send the spike to the target
  Node* target = get_target( t );


  DictionaryDatum d = new Dictionary();
  //( *d )[ nest::names::V_m ] = 0.0;

  target->get_status(d);

  double learning_rate = 0.0008;
  const double pre_th = 0.05;
  const double post_th = 0.15;
  const double tau_favg= 30.;


  double V_m = ( *d )[ names::V_m ];
  bool phi = ( *d )[ names::phi ];
  
  //V_m = target-> get_V_m_()

  omega_E += (-omega_E + phi) / tau_favg;

  int plasticity_type = 9;

   if ( tau_favg >= pre_th &&  V_m >= post_th )
   {
      learning_rate = learning_rate;
      plasticity_type = 0;
   } 

  else if (tau_favg >= pre_th && ((pre_th <= V_m ) && (V_m < post_th) ))
  {
      learning_rate = -learning_rate;
      plasticity_type = 1;
  }

  else if (tau_favg < pre_th && V_m >= post_th  )
  {
      learning_rate = -learning_rate;
      plasticity_type = 2;
  }
   
  else
  {
      learning_rate = 0;
      plasticity_type = 3;
  }

  weight_ += weight_ * learning_rate;
  
  membrane_potentials_.push_back( V_m );
  plasticity_flags_.push_back( plasticity_type );
  firing_rates_.push_back( omega_E );
  weights_.push_back( weight_ );
  deltas_.push_back( learning_rate );
  plasticity_flags_.push_back( plasticity_type );

  e.set_receiver( *get_target( t ) );
  e.set_weight( weight_);
  e.set_delay_steps( get_delay_steps() );
  e.set_rport( get_rport() );
  e();

  t_lastspike_ = t_spike;
}

template < typename targetidentifierT >
matco_synapse< targetidentifierT >::matco_synapse()
  : ConnectionBase()
  , weight_( 1.0 )
  , t_lastspike_( 0.0 )
{
}

template < typename targetidentifierT >
void
matco_synapse< targetidentifierT >::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );

  ( *d )[ names::y ] = firing_rates_;
  ( *d )[ names::y1 ] = membrane_potentials_;
  ( *d )[ names::y2 ] = deltas_;
  ( *d )[ names::y_0 ] = weights_;
  ( *d )[ names::y_1 ] = plasticity_flags_;

  def< long >( d, names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
matco_synapse< targetidentifierT >::set_status( const DictionaryDatum& d, ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );

  updateValue< double >( d, names::weight, weight_ );
}

} // namespace

#endif // MATCO_SYNAPSE_H
