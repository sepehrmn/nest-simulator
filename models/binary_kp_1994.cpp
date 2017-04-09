
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

namespace nest
{

    /*
 * template specialization must be placed in namespace
 *
 * Override the create() method with one call to RecordablesMap::insert_()
 * for each quantity to be recorded.
 */
template <>
void
RecordablesMap< aeif_psc_exp >::create()
{
  // use standard names whereever you can for consistency!
  insert_(
    names::V_m, &aeif_psc_exp::get_y_elem_< aeif_psc_exp::State_::V_M > );
}

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