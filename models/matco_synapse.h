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

// Includes from models:
#include "updater_device.h"

// Includes from nestkernel:
#include "connection.h"
#include "iaf_matco_2018.h"

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

class matcoCommonProperties : public CommonSynapseProperties
{
public:
  /**
   * Default constructor.
   * Sets all property values to defaults.
   */
  matcoCommonProperties();

  /**
   * Get all properties and put them into a dictionary.
   */
  void get_status( DictionaryDatum& d ) const;

  /**
   * Set properties from the values given in dictionary.
   */
  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  Node* get_node();

  long get_ut_node_id() const;

  updater_device* ut_;
};

inline long
matcoCommonProperties::get_ut_node_id() const
{
  if ( ut_ != 0 )
  {
    return ut_->get_node_id();
  }
  else
  {
    return -1;
  }
}

template < typename targetidentifierT >
class matco_synapse : public Connection< targetidentifierT >
{
public:
  typedef matcoCommonProperties CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;

  //SecondaryEvent* get_secondary_event();

  static constexpr ConnectionModelProperties properties = ConnectionModelProperties::HAS_DELAY
  | ConnectionModelProperties::IS_PRIMARY | ConnectionModelProperties::SUPPORTS_HPC
  | ConnectionModelProperties::SUPPORTS_LBL;

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
      return invalid_port;
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

  void force_update_weight( thread t,
    double t_trig,
    const matcoCommonProperties& cp );

private:
  double weight_; //!< Synaptic weight

  double t_lastspike_; //!< Time point of last spike emitted

  double omega_E_; 

  double tau_;

  bool phi_;

  double theta_;
  double theta_minus_;
  double theta_plus_;

  std::vector< double > firing_rates_;
  std::vector< double > membrane_potentials_;
  std::vector< double > deltas_;
  std::vector< double > weights_;
  std::vector< double > plasticity_flags_;

  //std::vector< double > tmp_filter_times( 0 );
  //S_.plasticity_flags_.swap( tmp_filter_times );

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
  phi_ = 1;
  // propagation t_lastspike -> t_spike, t_lastspike_ = 0 initially, p_ = 1
  const double t_spike = e.get_stamp().get_ms();
  const double h = t_spike - t_lastspike_;

  e.set_weight( weight_ );
  e.set_delay_steps( get_delay_steps() );
  e.set_receiver( *get_target( t ) );
  e.set_rport( get_rport() );
  e();

  t_lastspike_ = t_spike;
}



template < typename targetidentifierT >
matco_synapse< targetidentifierT >::matco_synapse()
  : ConnectionBase()
  , weight_( 1.0 )
  , t_lastspike_( 0.0 )
  , omega_E_(0.0)
  , tau_(30.)
  , phi_(0)
  , theta_(0.05)
  , theta_minus_(0.14)
  , theta_plus_(0.15)
{
}

template < typename targetidentifierT >
constexpr ConnectionModelProperties matco_synapse< targetidentifierT >::properties;

template < typename targetidentifierT >
void
matco_synapse< targetidentifierT >::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );

  ( *d )[ names::rate ] = firing_rates_;
  ( *d )[ names::V_m ] = membrane_potentials_;
  ( *d )[ names::delta ] = deltas_;
  ( *d )[ names::weight ] = weights_;
  ( *d )[ names::type_id ] = plasticity_flags_;

  ( *d )[ names::tau ] = tau_;

  def< long >( d, names::size_of, sizeof( *this ) );
}


  class ConnTestDummyNode : public ConnTestDummyNodeBase
  {
  public:
    // Ensure proper overriding of overloaded virtual functions.
    // Return values from functions are ignored.
    using ConnTestDummyNodeBase::handles_test_event;
    port
    handles_test_event( SpikeEvent&, rport )
    {
      return invalid_port;
    }
  };

template < typename targetidentifierT >
void
matco_synapse< targetidentifierT >::set_status( const DictionaryDatum& d, ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );

  updateValue< double >( d, names::weight, weight_ );
  updateValue< double >( d, names::tau, tau_ );
  updateValue< double >( d, names::theta, theta_ );
  updateValue< double >( d, names::theta_minus, theta_minus_ );
  updateValue< double >( d, names::theta_plus, theta_plus_ );
}

template < typename targetidentifierT >
inline void
matco_synapse< targetidentifierT >::force_update_weight( thread t,
  const double t_trig,
  const matcoCommonProperties& cp )
{
  

  iaf_matco_2018* target = reinterpret_cast<iaf_matco_2018*>(get_target( t ));

  double learning_rate = 0.0008;
  
  double V_m = target->get_V_m();

  omega_E_ += (-omega_E_ + phi_) / tau_;
  if (phi_ == 1)
  {
    phi_ = 0;
  }
  

  int plasticity_type = 9;
   
   // LTP
   if ( omega_E_>= theta_ &&  V_m >= theta_plus_ )
   {
      learning_rate = learning_rate;
      plasticity_type = 0;
   } 

  // LTD (homosynaptic)
  else if (omega_E_ >= theta_ && ((theta_minus_ <= V_m) && (V_m < theta_plus_)))
  {
      learning_rate = -learning_rate;
      plasticity_type = 1;
  }

  // LTD (heterosynaptic)
  else if (omega_E_ < theta_ && V_m >= theta_plus_)
  {
      learning_rate = -learning_rate;
      plasticity_type = 2;
  }
  
  // No plasticity
  else
  {
      learning_rate = 0;
      plasticity_type = 3;
  }

  weight_ += weight_ * learning_rate;
  
  membrane_potentials_.push_back( V_m );
  plasticity_flags_.push_back( plasticity_type );
  firing_rates_.push_back( omega_E_ );
  weights_.push_back( weight_ );
  deltas_.push_back( learning_rate );

}

// template < typename targetidentifierT >
// SecondaryEvent*
// matco_synapse< targetidentifierT >::get_secondary_event()
// {
//   return new matco_synapse();
// }

 } // namespace

#endif // MATCO_SYNAPSE_H
