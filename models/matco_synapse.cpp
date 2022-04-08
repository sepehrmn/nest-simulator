/*
 *  stdp_dopamine_synapse.cpp
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

#include "matco_synapse.h"

// Includes from nestkernel:
#include "common_synapse_properties.h"
#include "connector_model.h"
#include "event.h"
#include "kernel_manager.h"

// Includes from sli:
#include "dictdatum.h"

namespace nest
{
//
// Implementation of class STDPDopaCommonProperties.
//

matcoCommonProperties::matcoCommonProperties()
  : CommonSynapseProperties()
  , ut_( 0 )
{
}

void
matcoCommonProperties::get_status( DictionaryDatum& d ) const
{
  CommonSynapseProperties::get_status( d );
  if ( ut_ != 0 )
  {
    def< long >( d, names::U, ut_->get_node_id() );
  }
  else
  {
    def< long >( d, names::U, -1 );
  }

}

void
matcoCommonProperties::set_status( const DictionaryDatum& d, ConnectorModel& cm )
{
  CommonSynapseProperties::set_status( d, cm );

  long utnode_id;
  if ( updateValue< long >( d, names::U, utnode_id ) )
  {
    const thread tid = kernel().vp_manager.get_thread_id();
    Node* ut = kernel().node_manager.get_node_or_proxy( utnode_id, tid );
    ut_ = dynamic_cast< updater_device* >( ut );
    if ( ut_ == 0 )
    {
      throw BadProperty( "Must be updater device" );
    }
  }
}

Node*
matcoCommonProperties::get_node()
{
  if ( ut_ == 0 )
  {
    throw BadProperty( "No weight updater assigned to synapse." );
  }
  else
  {
    return ut_;
  }
}

} // of namespace nest
