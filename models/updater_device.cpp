/*
 *  volume_transmitter.cpp
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

#include "updater_device.h"

// C++ includes:
#include <numeric>

// Includes from nestkernel:
#include "connector_base.h"
#include "exceptions.h"
#include "kernel_manager.h"
#include "spikecounter.h"

// Includes from libnestutil:
#include "dict_util.h"

// Includes from sli:
#include "arraydatum.h"
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"

/* ----------------------------------------------------------------
 * Default constructor defining default parameters
 * ---------------------------------------------------------------- */

nest::updater_device::Parameters_::Parameters_()
  : deliver_interval_( 1 ) // in steps of mindelay
{
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
nest::updater_device::Parameters_::get( DictionaryDatum& d ) const
{
  def< long >( d, names::deliver_interval, deliver_interval_ );
}

void ::nest::updater_device::Parameters_::set( const DictionaryDatum& d, Node* node )
{
  updateValueParam< long >( d, names::deliver_interval, deliver_interval_, node );
}

/* ----------------------------------------------------------------
 * Default and copy constructor for volume transmitter
 * ---------------------------------------------------------------- */

nest::updater_device::updater_device()
  : Node()
  , P_()
  , local_device_id_( 0 )
{
}

nest::updater_device::updater_device( const updater_device& n )
  : Node( n )
  , P_( n.P_ )
  , local_device_id_( n.local_device_id_ )
{
}

void
nest::updater_device::init_buffers_()
{
}

void
nest::updater_device::pre_run_hook()
{
}

void
nest::updater_device::update( const Time&, const long from, const long to )
{
  
  for ( long lag = from; lag < to; ++lag )
  {
  // all spikes stored in spikecounter_ are delivered to the target synapses
    if ( ( kernel().simulation_manager.get_slice_origin().get_steps() + to )
      % ( P_.deliver_interval_ * kernel().connection_manager.get_min_delay() )
    == 0 )
   {
      double t_trig = Time( Time::step( kernel().simulation_manager.get_slice_origin().get_steps() + to ) ).get_ms();
      kernel().connection_manager.force_update_weight( get_node_id(), t_trig );
    }
  }
}

void
nest::updater_device::handle( SpikeEvent& e )
{
}
