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
//
//#ifndef KAY_PHILLIPS_CONNECTION_H
//#define KAY_PHILLIPS_CONNECTION_H
//
//// Includes from nestkernel:
//#include "connection.h"
//
///* BeginDocumentation
//  Name: KAY_PHILLIPS synapse - Synapse with depression after Hill & Tononi (2005).
//
//  Description:
//  This synapse implements the depression model described in [1].
//  Updating weights does not take place here. this is needed so that "syn_spec"
//
//  Parameters:
//     The following parameters can be set in the status dictionary:
//     tauP     double - synaptic vesicle pool recovery time constant [ms]
//
//  References:
//   [1] J Kay and WA Phillips (1994). Technical Report CCCN-15.
//
//  Sends: SpikeEvent
//
//  FirstVersion: January 2017
//  Author: Sepehr Mahmoudian
//  SeeAlso: kay_phillips_neuron
//*/
//
//namespace nest
//{
//
//    template < typename targetidentifierT >
//    class KP1994Connection : public Connection< targetidentifierT >
//    {
//    public:
//        typedef CommonSynapseProperties CommonPropertiesType;
//        typedef Connection< targetidentifierT > ConnectionBase;
//
//        /**
//         * Default Constructor.
//         * Sets default values for all parameters. Needed by GenericConnectorModel.
//         */
//        kay_philips_connection();
//
//        /**
//         * Copy constructor.
//         * Needs to be defined properly in order for GenericConnector to work.
//         */
//        kay_philips_connection( const kay_philips_connection& );
//
//        // Explicitly declare all methods inherited from the dependent base
//        // ConnectionBase. This avoids explicit name prefixes in all places these
//        // functions are used. Since ConnectionBase depends on the template parameter,
//        // they are not automatically found in the base class.
//        using ConnectionBase::get_delay_steps;
//        using ConnectionBase::get_delay;
//        using ConnectionBase::get_rport;
//        using ConnectionBase::get_target;
//
//        /**
//         * Default Destructor.
//         */
//        virtual ~kay_philips_connection()
//        {
//        }
//
//        /**
//         * Get all properties of this connection and put them into a dictionary.
//         */
//        virtual void get_status( DictionaryDatum& d ) const;
//
//        /**
//         * Set properties of this connection from the values given in dictionary.
//         */
//        virtual void set_status( const DictionaryDatum& d, ConnectorModel& cm );
//
//        /**
//         * Send an event to the receiver of this connection.
//         * \param e The event to send
//         * \param t_lastspike Point in time of last spike sent.
//         * \param cp Common properties to all synapses (empty).
//         */
//        void send( Event& e,
//                   thread t,
//                   double t_lastspike,
//                   const CommonSynapseProperties& cp );
//
//        class ConnTestDummyNode : public ConnTestDummyNodeBase
//        {
//        public:
//            // Ensure proper overriding of overloaded virtual functions.
//            // Return values from functions are ignored.
//            using ConnTestDummyNodeBase::handles_test_event;
//            port
//            handles_test_event( SpikeEvent&, rport )
//            {
//                return invalid_port_;
//            }
//        };
//
//        void
//        check_connection( Node& s,
//                          Node& t,
//                          rport receptor_type,
//                          double,
//                          const CommonPropertiesType& )
//        {
//            ConnTestDummyNode dummy_target;
//            ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );
//        }
//
//        //! allows efficient initialization from ConnectorModel::add_connection()
//        void
//        set_weight( double w )
//        {
//            weight_ = w;
//        }
//
//    private:
//        double weight_; //!< synpatic weight
//
//    };
//
//
///**
// * Send an event to the receiver of this connection.
// * \param e The event to send
// * \param p The port under which this connection is stored in the Connector.
// * \param t_lastspike Time point of last spike emitted
// */
//    template < typename targetidentifierT >
//    inline void
//    HTConnection< targetidentifierT >::send( Event& e,
//                                             thread t,
//                                             double t_lastspike,
//                                             const CommonSynapseProperties& )
//    {
//        double h = e.get_stamp().get_ms() - t_lastspike;
//        Node* target = get_target( t );
//        // t_lastspike_ = 0 initially
//
//        // propagation t_lastspike -> t_spike, t_lastspike_ = 0 initially, p_ = 1
//        p_ = 1 - ( 1 - p_ ) * std::exp( -h / tau_P_ );
//
//        // send the spike to the target
//        e.set_receiver( *target );
//        e.set_weight( weight_ * p_ );
//        e.set_delay( get_delay_steps() );
//        e.set_rport( get_rport() );
//        e();
//
//        // reduce pool after spike is sent
//        p_ *= ( 1 - delta_P_ );
//    }
//
//    template < typename targetidentifierT >
//    HTConnection< targetidentifierT >::HTConnection()
//            : ConnectionBase()
//            , weight_( 1.0 )
//            , tau_P_( 50.0 )
//            , delta_P_( 0.2 )
//            , p_( 1.0 )
//    {
//    }
//
//    template < typename targetidentifierT >
//    KP1994Connection< targetidentifierT >::KP1994Connection( const KP1994Connection& rhs )
//            : ConnectionBase( rhs )
//            , weight_( rhs.weight_ )
//            , tau_P_( rhs.tau_P_ )
//            , delta_P_( rhs.delta_P_ )
//            , p_( rhs.p_ )
//    {
//    }
//
//    template < typename targetidentifierT >
//    void
//    KP1994Connection< targetidentifierT >::get_status( DictionaryDatum& d ) const
//    {
//        ConnectionBase::get_status( d );
//        def< double >( d, names::weight, weight_ );
//        def< long >( d, names::size_of, sizeof( *this ) );
//    }
//
//    template < typename targetidentifierT >
//    void
//    KP1994Connection< targetidentifierT >::set_status( const DictionaryDatum& d,
//                                                   ConnectorModel& cm )
//    {
//        ConnectionBase::set_status( d, cm );
//
//        updateValue< double >( d, names::weight, weight_ );
//        updateValue< double >( d, names::contextual_field, tau_P_ );
//        updateValue< double >( d, names::receptive_field, tau_P_ );
//    }
//
//} // namespace
//
//#endif // KAY_PHILLIPS_CONNECTION_H
