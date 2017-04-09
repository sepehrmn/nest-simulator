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
  The neuron model described in [1]. In [2] the modulatory attributes of this
  model are highlighted.
  References:
   [1] J Kay and WA Phillips (1994). Technical Report CCCN-15.
   [2] D Smyth, J Kay, WA Phillips (1996). Network: Computation in Neural Systems.
  FirstVersion: April 2017
  Author: Sepehr Mahmoudian
  SeeAlso: kp_1994_connections
*/

#ifndef BINARY_KP_1994_H
#define BINARY_KP_1994_H

namespace nest
{
    class binary_kp_1994 : public Node
    {
        //using Node::handle;
        using Node::handles_test_event;

        void get( DictionaryDatum& ) const; //!< Store current values in dictionary
        void set( const DictionaryDatum& ); //!< Set values from dicitonary
    };

inline port
binary_kp_1994::handles_test_event(SpikeEvent &e, rport receptor_type)
{
    if (receptor_type != 0) {
        throw UnknownReceptorType(receptor_type, get_name());
    }

    S_.node_gids_.push_back(e.get_sender().get_gid());
    return S_.node_gids_.size() - 1; // -1 because we want rports to start from 0
}

} //namespace

#endif //BINARY_KP_1994_H
