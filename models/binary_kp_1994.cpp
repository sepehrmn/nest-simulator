
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

// Includes from libnestutil:
#include "numerics.h"

// Includes from nestkernel:
#include "exceptions.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"

/* ----------------------------------------------------------------
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap< nest::binary_kp_1994 >
  nest::binary_kp_1994::recordablesMap_;

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
  RecordablesMap< binary_kp_1994 >::create()
  {
    // use standard names whereever you can for consistency!
    insert_( names::receptive_field, &binary_kp_1994::get_receptive_field_ );
    insert_( names::contextual_field, &binary_kp_1994::get_contextual_field_ );
    insert_( names::theta, &binary_kp_1994::get_theta_ );
  }
} //namespace


  /* ----------------------------------------------------------------
   * Default and copy constructor for node
   * ---------------------------------------------------------------- */

   nest::binary_kp_1994::binary_kp_1994()
     : Archiving_Node()
     , P_()
     , S_()
     , B_( *this )
   {
     recordablesMap_.create();
   }

   nest::binary_kp_1994::binary_kp_1994( const binary_kp_1994& n )
     : Archiving_Node( n )
     , P_( n.P_ )
     , S_( n.S_ )
     , B_( n.B_, *this )
   {

   }

   /* ----------------------------------------------------------------
    * Default constructors defining default parameters and state
    * ---------------------------------------------------------------- */

   nest::binary_kp_1994::Parameters_::Parameters_()
     : k1_( 0.5 )
     , k2_( 2.0 )
     , k3_( 0.0 )
   {
   }

   nest::binary_kp_1994::State_::State_()
     : theta_( 0.0 )
     ,  w_0_ ( 0.0 )
     ,  v_0_ ( 0.0 )
   {
   }

 /* ----------------------------------------------------------------
  * Parameter and state extractions and manipulation functions
  * ---------------------------------------------------------------- */

void
nest::binary_kp_1994::get_status( DictionaryDatum& d ) const
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

  ( *d )[ names::receptive_field ] = S_.receptive_field_;
  ( *d )[ names::contextual_field ] = S_.contextual_field_;

  ( *d )[ names::recordables ] = recordablesMap_.get_list();
}

void
nest::binary_kp_1994::set_status( const DictionaryDatum& d )
{
  P_.set( d);
  P_.set( d);
}

// Parameters

void
nest::binary_kp_1994::Parameters_::get( DictionaryDatum& d ) const
{
  def< double >( d, names::k1, k1_ );
  def< double >( d, names::k2, k2_ );
  def< double >( d, names::k3, k3_ );
}

void
nest::binary_kp_1994::Parameters_::set( const DictionaryDatum& d )
{
  updateValue< double >( d, names::k1, k1_ );
  updateValue< double >( d, names::k2, k2_ );
  updateValue< double >( d, names::k3, k3_ );
}

// State

void
nest::binary_kp_1994::State_::get( DictionaryDatum& d) const
{
  def< double >( d, names::theta, theta_);
  def< double >( d, names::w_0, w_0_ );
  def< double >( d, names::v_0, v_0_ );
  def< double >( d, names::receptive_field, receptive_field_ );
  def< double >( d, names::contextual_field, contextual_field_ );
}

void
nest::binary_kp_1994::State_::set( const DictionaryDatum& d )
{
  updateValue< double >( d, names::theta, theta_ );
  updateValue< double >( d, names::w_0, w_0_ );
  updateValue< double >( d, names::v_0, v_0_ );
  updateValue< double >( d, names::receptive_field, receptive_field_ );
  updateValue< double >( d, names::contextual_field, contextual_field_ );
}

// Buffers

nest::binary_kp_1994::Buffers_::Buffers_( binary_kp_1994& n )
  : logger_( n )
{
}

nest::binary_kp_1994::Buffers_::Buffers_( const Buffers_&, binary_kp_1994& n )
  : logger_( n )
{
}

/* ----------------------------------------------------------------
* Node initialization functions
* ---------------------------------------------------------------- */

void
nest::binary_kp_1994::init_state_( const Node& proto )
{
  const binary_kp_1994& pr = downcast< binary_kp_1994 >( proto );
  S_ = pr.S_;
}

void
nest::binary_kp_1994::init_buffers_()
{
  B_.spikes_rf_.clear(); // includes resize
  B_.spikes_cf_.clear(); // includes resize
  B_.logger_.reset();
  Archiving_Node::clear_history();
}

void
nest::binary_kp_1994::calibrate()
{
  B_.logger_.init();
}

/* ----------------------------------------------------------------
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void
nest::binary_kp_1994::update( const Time& origin, const long from, const long to )
{
  assert(
    to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
  assert( from < to );

  for ( long lag = from; lag < to; ++lag )
  {
    // log state data
    B_.logger_.record_data( origin.get_steps() + lag );
  }
}

void
nest::binary_kp_1994::handle( SpikeEvent& e )
{
  assert( e.get_delay() > 0 );

  if ( e.get_rport() == RF )
  {
    B_.spikes_rf_.add_value( e.get_rel_delivery_steps(
                               kernel().simulation_manager.get_slice_origin() ),
      e.get_weight() * e.get_multiplicity() );
  }
  else
  {
    B_.spikes_cf_.add_value( e.get_rel_delivery_steps(
                               kernel().simulation_manager.get_slice_origin() ),
      e.get_weight() * e.get_multiplicity() );
  }
}

void
nest::binary_kp_1994::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}
