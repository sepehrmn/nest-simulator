/*
 *  bpid_kp_2017.h
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

/* BeginDocumentation
  Name: bpid_kp_2017 - binary stochastic neuron introduced in [1].
  Description:
  The neuron model described in [1][2]. In [3] the modulatory attributes of this
  model are highlighted.
  References:
   [1] J Kay (1994) Technical Report, Biomathematics and Statistics Scotland
   [2] J Kay and WA Phillips (1994). Technical Report CCCN-15.
   [3] D Smyth, J Kay, WA Phillips (1996). Network: Computation in Neural Systems.
   [4] Partial information decomposition as a unified approach to the specification of neural goal functions
Michael Wibral, Viola Priesemann, Jim W. Kay, Joseph T. Lizier, William A. Phillips

  Note that the receptive_field and contextual_field include w_0 and v_0
  respectively, and that are used to calculate the output probability.
  Add w_0 to receptive_field and v_0 to contextual_field to get the value of
  only the neural input added.

  // by default integration_type is additive and k1 and k2 are set to 1.0
  // if this is changed to modulatory, k1 and k2 are changed to 0.5 and 2.0,
  // respectively, unless specified otherwise. If the integration type of the neuron is changed,
  // and k1/k2 are not specified, they are set to the default values for that integration type.

  FirstVersion: November 2018
  Author: Sepehr Mahmoudian
  SeeAlso: bpid_kp_connection
*/

// Includes from librandom:
#include "exp_randomdev.h"

// Includes from libnestutil:
#include "logging.h"

// Includes from nestkernel:
#include "conn_parameter.h"
#include "exceptions.h"
#include "kernel_manager.h"

#include "archiving_node.h"
#include "connection.h"
#include "event.h"
#include "nest_types.h"
#include "recordables_map.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"

#ifndef bpid_kp_2017_H
#define bpid_kp_2017_H

namespace nest
{
    class bpid_kp_2017 : public Archiving_Node
    {

    public:
        bpid_kp_2017();
        bpid_kp_2017( const bpid_kp_2017& );

        /**
         * Import sets of overloaded virtual functions.
         * @see Technical Issues / Virtual Functions: Overriding, Overloading,
         * and Hiding
         */
        using Node::handle;
        using Node::handles_test_event;

        port send_test_event( Node&, rport, synindex, bool );

        void handle( SpikeEvent& );
        void handle( DataLoggingRequest& );

        port handles_test_event( SpikeEvent&, rport );
        port handles_test_event( DataLoggingRequest&, rport );

        void get_status( DictionaryDatum& ) const;
        void set_status( const DictionaryDatum& );

    private:

      enum SynapseTypes
      {
        SPIKE_RECEPTOR = 0, // for spike_detector
        RF,
        CF
      };

      enum IntegrationTypes
      {
        ADDITIVE = 0,
        MODULATORY
      };

    void init_state_( const Node& proto );
    void init_buffers_();
    void calibrate();
    void update( const Time&, const long, const long );

    void learn_();

    // The next two classes need to be friends to access the State_ class/member
    friend class RecordablesMap< bpid_kp_2017 >;
    friend class UniversalDataLogger< bpid_kp_2017 >;

    // Mapping of recordables names to access functions
    static RecordablesMap< bpid_kp_2017 > recordablesMap_;

    struct Parameters_
    {


      // k1 parameter of the activation function
      double k1_;
      // k2 parameter of the activation function
      double k2_;
      // integration type (e.g., additive or modulatory)
      int integration_type_;

      // phis
      std::vector< double > phis_;

      double alpha_;
      double beta_;
      double eta_;

      // Update interval
      double interval_;

      Parameters_();
      void get( DictionaryDatum& ) const;
      void set( const DictionaryDatum& );
    };

    struct State_
    {
      int update_counter_;

      double E;  // Average output probability for joint r and c
      double E_c_;  // Average output probability conditioned on c
      double E_r_;  // Average output probability conditioned on r

      // Output probability
      double theta_; // Output probability

      // Integrated receptive field bias
      double w_0_;
      // Integrated contextual field bias
      double v_0_;

      // integrated receptive field. B_.spikes_rf_ - w_0
      double receptive_field_;
      // integrated contextual field. B_.spikes_cf_ - v_0
      double contextual_field_;

      // Default initialization
      State_();
      void get( DictionaryDatum& ) const;
      void set( const DictionaryDatum& );
    };

    struct Buffers_
    {
      // Buffers for storing RF and CF spikes
      RingBuffer spikes_rf_;
      RingBuffer spikes_cf_;

      // Logger for all analog data
      UniversalDataLogger< bpid_kp_2017 > logger_;

      Buffers_( bpid_kp_2017& );
      Buffers_( const Buffers_&, bpid_kp_2017& );
    };

    struct Variables_
    {
      int r_counter = 0;
      int c_counter = 0;

      librandom::RngPtr rng_; // pointer to thread specific random generator
    };

    // Access functions for UniversalDataLogger -------------------------------

    double
    get_receptive_field_() const
    {
      return S_.receptive_field_;
    }

    double
    get_contextual_field_() const
    {
      return S_.contextual_field_;
    }

    double
    get_theta_() const
    {
      return S_.theta_;
    }

    // ----------------------------------------------------------------

    Parameters_ P_;
    State_ S_;
    Buffers_ B_;
    Variables_ V_;

};

   inline port
   bpid_kp_2017::handles_test_event( SpikeEvent&, rport receptor_type )
   {
     if ( !( receptor_type == SPIKE_RECEPTOR
          || receptor_type == RF
          || receptor_type == CF ) )
     {
       throw UnknownReceptorType( receptor_type, get_name() );
     }

     else
     {

       std::string rt_str = std::to_string(receptor_type);

       if ( receptor_type == RF )
       {
         rt_str.append(std::to_string(V_.r_counter));
         V_.r_counter += 1;
       }

       else
       {
         rt_str.append(std::to_string(V_.c_counter));
         V_.c_counter += 1;
       }

       receptor_type = std::stoi(rt_str);

     }

     return receptor_type;
   }

   inline port
   bpid_kp_2017::send_test_event( Node& target,
     rport receptor_type,
     synindex,
     bool )
   {
     SpikeEvent e;
     e.set_sender( *this );
     return target.handles_test_event( e, receptor_type );
   }

   inline port
   bpid_kp_2017::handles_test_event( DataLoggingRequest& dlr, rport receptor_type )
   {
     if ( receptor_type != 0 )
     {
       throw UnknownReceptorType( receptor_type, get_name() );
     }
     return B_.logger_.connect_logging_device( dlr, recordablesMap_ );
   }


} //namespace

#endif //bpid_kp_2017_H
