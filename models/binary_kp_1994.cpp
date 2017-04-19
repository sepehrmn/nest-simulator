
/*
 *  binary_kp_1994.cpp
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

#include "binary_kp_1994.h"

/* ----------------------------------------------------------------
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap< nest::binary_kp_1994 >
  nest::binary_kp_1994::recordablesMap_;

namespace nest
{

  // Override the create() method with one call to RecordablesMap::insert_()
  // for each quantity to be recorded.
  template <>
  void
  RecordablesMap< iaf_psc_exp >::create()
  {
    // use standard names whereever you can for consistency!
    insert_( names::V_m, &iaf_psc_exp::get_V_m_ );
    insert_( names::weighted_spikes_ex, &iaf_psc_exp::get_weighted_spikes_ex_ );
    insert_( names::weighted_spikes_in, &iaf_psc_exp::get_weighted_spikes_in_ );
    insert_( names::I_syn_ex, &iaf_psc_exp::get_I_syn_ex_ );
    insert_( names::I_syn_in, &iaf_psc_exp::get_I_syn_in_ );
  }
  }

  /* ----------------------------------------------------------------
   * Default constructors defining default parameters and state
   * ---------------------------------------------------------------- */

  nest::iaf_psc_exp::Parameters_::Parameters_()
    : Tau_( 10.0 )             // in ms
    , C_( 250.0 )              // in pF
    , t_ref_( 2.0 )            // in ms
    , E_L_( -70.0 )            // in mV
    , I_e_( 0.0 )              // in pA
    , Theta_( -55.0 - E_L_ )   // relative E_L_
    , V_reset_( -70.0 - E_L_ ) // in mV
    , tau_ex_( 2.0 )           // in ms
    , tau_in_( 2.0 )           // in ms
  {
  }

  nest::iaf_psc_exp::State_::State_()
    : i_0_( 0.0 )
    , i_syn_ex_( 0.0 )
    , i_syn_in_( 0.0 )
    , V_m_( 0.0 )
    , r_ref_( 0 )
  {
  }

/*
 * template specialization must be placed in namespace
 *
 * Override the create() method with one call to RecordablesMap::insert_()
 * for each quantity to be recorded.
 */
 template <>
 void
 RecordablesMap< binary_kp_1994 >::create()
 {
   // use standard names wherever you can for consistency!
   insert_( names::I_syn_in, &amat2_psc_exp::get_I_syn_in_ );
 }


 /* ----------------------------------------------------------------
  * Parameter and state extractions and manipulation functions
  * ---------------------------------------------------------------- */
  
// receptive_field;
// contextual_field
void
binary_kp_1994::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  S_.get( d );
  Archiving_Node::get_status( d );

  DictionaryDatum receptor_types = new Dictionary();

  ( *receptor_types )[ "receptive_field" ] = RF;
  ( *receptor_types )[ "contextual_field" ] = CF;

  ( *d )[ "receptor_types" ] = receptor_types;
  ( *d )[ names::recordables ] = recordablesMap_.get_list();
}

void
nest::amat2_psc_exp::update( Time const& origin,
  const long from,
  const long to )
{
  assert(
    to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
  assert( from < to );

  // evolve from timestep 'from' to timestep 'to' with steps of h each
  for ( long lag = from; lag < to; ++lag )
  {
        SpikeEvent se;
        kernel().event_delivery_manager.send( *this, se, lag );
      }
    }
    else
    {
      --S_.r_;
    } // neuron is totally refractory (cannot generate spikes)

    // log state data
    B_.logger_.record_data( origin.get_steps() + lag );
  }
}
