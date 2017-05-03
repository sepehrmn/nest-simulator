
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
    insert_( names::weighted_spikes_ex, &binary_kp_1994::get_weighted_spikes_ex_ );
    insert_( names::weighted_spikes_ex, &kkk::get_weighted_spikes_ex_ );
  }


  /* ----------------------------------------------------------------
   * Default constructors defining default parameters and state
   * ---------------------------------------------------------------- */

  nest::iaf_psc_exp::Parameters_::Parameters_()
    : Tau_( 10.0 )             // in ms
  {
  }

  nest::iaf_psc_exp::State_::State_()
    : i_0_( 0.0 )
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
   insert_( names::receptive_field, &binary_kp_1994::get_receptive_field_ );
   insert_( names::contextual_field, &binary_kp_1994::get_contextual_field_ );
   contextual_field
 }


 /* ----------------------------------------------------------------
  * Parameter and state extractions and manipulation functions
  * ---------------------------------------------------------------- */

void
binary_kp_1994::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  S_.get( d );
  Archiving_Node::get_status( d );

  DictionaryDatum receptor_types = new Dictionary();

  // For easy assignment, returning the enum value for RF and CF
  ( *receptor_types )[ "receptive_field" ] = RF;
  ( *receptor_types )[ "contextual_field" ] = CF;
  ( *d )[ "receptor_types" ] = receptor_types;

  ( *d )[ names::recordables ] = recordablesMap_.get_list();
}

void
binary_kp_1994::update( Time const& origin,
  const long from,
  const long to )
{
  assert(
    to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
  assert( from < to );

    // log state data
    B_.logger_.record_data( origin.get_steps() + lag );
 }

}
