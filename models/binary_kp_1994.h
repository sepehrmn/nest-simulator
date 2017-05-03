/*
 *  binary_kp_1994.h
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
  Name: binary_kp_1994 - binary stochastic neuron introduced in [1].
  Description:
  The neuron model described in [1][2]. In [3] the modulatory attributes of this
  model are highlighted.
  References:
   [1] J Kay (1994) Technical Report, Biomathematics and Statistics Scotland
   [2] J Kay and WA Phillips (1994). Technical Report CCCN-15.
   [3] D Smyth, J Kay, WA Phillips (1996). Network: Computation in Neural Systems.
  FirstVersion: May 2017
  Author: Sepehr Mahmoudian
  SeeAlso: kp_1994_connection
*/

#include "archiving_node.h"
#include "connection.h"
#include "event.h"
#include "nest_types.h"
#include "recordables_map.h"
//#include "universal_data_logger.h"

#ifndef BINARY_KP_1994_H
#define BINARY_KP_1994_H

namespace nest
{
    class binary_kp_1994 : public Archiving_Node
    {

    public:
        binary_kp_1994();
        binary_kp_1994( const binary_kp_1994& );

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
        SPIKE_RECEPTOR = 0,
        RF,
        CF
      };

      void init_state_( const Node& proto );
      void init_buffers_();
      void calibrate();
      void update( const Time&, const long, const long );

      // The next two classes need to be friends to access the State_ class/member
      friend class RecordablesMap< iaf_psc_exp >;
      friend class UniversalDataLogger< iaf_psc_exp >;

      // Mapping of recordables names to access functions
      static RecordablesMap< binary_kp_1994 > recordablesMap_;

    };

    struct Parameters_
    {

      // k1 parameter of the activation function
      double k1_;
      // k2 parameter of the activation function
      double k2_;

      Parameters_();
    };

    struct State_
    {
      // Output probability
      double Theta_;

      // Integrated receptive field bias
      double w_0;
      // Integrated contextual field bias
      double v_0;

      // Default initialization
      State_();
    };

    struct Buffers_
    {
      // Buffers for storing RF and CF spikes
      RingBuffer spikes_rf_;
      RingBuffer spikes_cf_;

      // Logger for all analog data
      UniversalDataLogger< iaf_psc_exp > logger_;

      Buffers_( binary_kp_1994& );
      Buffers_( const Buffers_&, binary_kp_1994& );
    };

    struct Variables_
    {
      // integrated receptive field. B_.spikes_rf_ - w_0
      double receptive_field_;
      // integrated contextual field. B_.spikes_cf_ - v_0
      double contextual_field_;

    };

    // Access functions for UniversalDataLogger -------------------------------

    inline double
    get_weighted_spikes_ex_() const
    {
      return V_.weighted_spikes_ex_;
    }

    inline double
    get_weighted_spikes_in_() const
    {
      return V_.weighted_spikes_in_;
    }
   }

   inline port
   binary_kp_1994::handles_test_event( SpikeEvent&, rport receptor_type )
   {
     if ( !( receptor_type == SPIKE_RECEPTOR
          || receptor_type == RF
          || receptor_type == CF ) )
     {
       throw UnknownReceptorType( receptor_type, get_name() );
     }
     else
       return receptor_type;
   }

   
} //namespace

#endif //BINARY_KP_1994_H
