///*
// *  kay_phillips_connection.h
// *
// *  This file is part of NEST.
// *
// *  Copyright (C) 2004 The NEST Initiative
// *
// *  NEST is free software: you can redistribute it and/or modify
// *  it under the terms of the GNU General Public License as published by
// *  the Free Software Foundation, either version 2 of the License, or
// *  (at your option) any later version.
// *
// *  NEST is distributed in the hope that it will be useful,
// *  but WITHOUT ANY WARRANTY; without even the implied warranty of
// *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// *  GNU General Public License for more details.
// *
// *  You should have received a copy of the GNU General Public License
// *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
// *
// */

#ifndef KAY_PHILLIPS_CONNECTION_H
#define KAY_PHILLIPS_CONNECTION_H

// Includes from nestkernel:
#include "connection.h"

/* BeginDocumentation
 Name: KAY_PHILLIPS synapse - Synapse for use with binary_kp_1994.

 Description:
 This synapse should be used with the binary_kp_1994 model.

 Parameters:
    The following parameters can be set in the status dictionary:
    eta     double - learning rate parameter

 References:
  [1] J Kay and WA Phillips (1994). Technical Report CCCN-15.

 Sends: SpikeEvent

 FirstVersion: January 2017
 Author: Sepehr Mahmoudian
 SeeAlso: binary_kp_1994
*/

namespace nest
{

   template < typename targetidentifierT >
   class KP1994Connection : public Connection< targetidentifierT >
   {
   public:
       typedef CommonSynapseProperties CommonPropertiesType;
       typedef Connection< targetidentifierT > ConnectionBase;

       /**
        * Default Constructor.
        * Sets default values for all parameters. Needed by GenericConnectorModel.
        */
       KP1994Connection():
       ConnectionBase(),
       weight_(0);
       {
       }

       /**
        * Copy constructor.
        * Needs to be defined properly in order for GenericConnector to work.
        */
       kay_phillips_connection( const kay_phillips_connection& );

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
       virtual ~kay_phillips_connection()
       {
       }

       /**
        * Get all properties of this connection and put them into a dictionary.
        */
       void get_status( DictionaryDatum& d ) const;

       /**
        * Set properties of this connection from the values given in dictionary.
        */
       void set_status( const DictionaryDatum& d, ConnectorModel& cm );

       /**
        * Send an event to the receiver of this connection.
        * \param e The event to send
        * \param t_lastspike Point in time of last spike sent.
        * \param cp Common properties to all synapses (empty).
        */
       void send( Event& e,
                  thread t,
                  double t_lastspike,
                  const CommonSynapseProperties& cp );

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

       /*
        * This function calls check_connection on the sender and checks if the
        * receiver accepts the event type and receptor type requested by the sender.
        * Node::check_connection() will either confirm the receiver port by returning
        * true or false if the connection should be ignored.
        * We have to override the base class' implementation, since for STDP
        * connections we have to call register_stdp_pl_connection on the target
        * neuron to inform the Archiver to collect spikes for this connection.
        * Further, the STDP dopamine synapse requires a volume transmitter to be set
        * before any simulation is performed. Checking this satisfies ticket #926.
        *
        * \param s The source node
        * \param r The target node
        * \param receptor_type The ID of the requested receptor type
        * \param t_lastspike last spike produced by presynaptic neuron (in ms)
        */
       void
       check_connection( Node& s,
                         Node& t,
                         rport receptor_type,
                         double,
                         const CommonPropertiesType& )
       {
           ConnTestDummyNode dummy_target;
           ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );
       }

       void
       send( Event& e, thread t, double, const CommonSynapseProperties& )
       {
         e.set_weight( weight_ );
         e.set_delay( get_delay_steps() );
         e.set_receiver( *get_target( t ) );
         e.set_rport( get_rport() );
         e();
       }

       //! allows efficient initialization from ConnectorModel::add_connection()
       void
       set_weight( double w )
       {
           weight_ = w;
       }

   private:
       double weight_; //!< synpatic weight

   };



   template < typename targetidentifierT >
   KP1994Connection< targetidentifierT >::KP1994Connection()
           : ConnectionBase()
           , weight_( 1.0 )
   {
   }

   template < typename targetidentifierT >
   void
   KP1994Connection< targetidentifierT >::get_status( DictionaryDatum& d ) const
   {
       ConnectionBase::get_status( d );
       def< double >( d, names::weight, weight_ );
       def< long >( d, names::size_of, sizeof( *this ) );
   }

   template < typename targetidentifierT >
   void
   KP1994Connection< targetidentifierT >::set_status( const DictionaryDatum& d,
                                                  ConnectorModel& cm )
   {
       ConnectionBase::set_status( d, cm );
       updateValue< double >( d, names::weight, weight_ );
   }

   /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param p The port under which this connection is stored in the Connector.
   * \param t_lastspike Time point of last spike emitted
   */
   void
   send( Event& e, thread t, double, const CommonSynapseProperties& )
   {
     e.set_weight( weight_ );
     e.set_delay( get_delay_steps() );
     e.set_receiver( *get_target( t ) );
     e.set_rport( get_rport() );
     e();
   }

} // namespace

#endif // KAY_PHILLIPS_CONNECTION_H
