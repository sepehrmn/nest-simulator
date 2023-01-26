/*
 *  iaf_matco_2018.cpp
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

#include "iaf_matco_2018.h"

// C++ includes:
#include <limits>

// Includes from libnestutil:
#include "dict_util.h"
#include "numerics.h"

// Includes from nestkernel:
#include "event_delivery_manager_impl.h"
#include "exceptions.h"
#include "kernel_manager.h"
#include "ring_buffer_impl.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"

/* ----------------------------------------------------------------
 * Recordables map  e
 * ---------------------------------------------------------------- */

nest::RecordablesMap< nest::iaf_matco_2018 > nest::iaf_matco_2018::recordablesMap_;

namespace nest
{
// Override the create() method with one call to RecordablesMap::insert_()
// for each quantity to be recorded.
template <>
void
RecordablesMap< iaf_matco_2018 >::create()
{
  // use standard names whereever you can for consistency!
  insert_( names::V_m, &iaf_matco_2018::get_V_m );
  insert_( names::I_syn_ex, &iaf_matco_2018::get_I_syn_ex_ );
  insert_( names::I_syn_in, &iaf_matco_2018::get_I_syn_in_ );
}
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

nest::iaf_matco_2018::Parameters_::Parameters_()
  : Tau_( 10.0 )             // in ms (timesteps?)
  , k1_( 0.01 )
  , I_e_( 0.0 )              // in pA
  , Theta_( 0.18)              // spiking threshold
  , tau_ex_( 2.5 )
  , tau_in_( 5.0 )
  , alpha_( 7.0 )
{
}

nest::iaf_matco_2018::State_::State_()
  : i_syn_ex_( 0.0 )
  , i_syn_in_( 0.0 )
  , V_m_( 0.0 )
  , omega_( 0.0 )
  , phi_( false )
{
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
nest::iaf_matco_2018::Parameters_::get( DictionaryDatum& d ) const
{
  def< double >( d, names::I_e, I_e_ );
  def< double >( d, names::theta, Theta_); // threshold value
  def< double >( d, names::tau, Tau_ );
  def< double >( d, names::tau_syn_ex, tau_ex_ );
  def< double >( d, names::tau_syn_in, tau_in_ );
}

double
nest::iaf_matco_2018::Parameters_::set( const DictionaryDatum& d, Node* node )
{
  updateValueParam< double >( d, names::theta, Theta_, node );
  updateValueParam< double >( d, names::I_e, I_e_, node );
  updateValueParam< double >( d, names::tau, Tau_, node );
  updateValueParam< double >( d, names::tau_syn_ex, tau_ex_, node );
  updateValueParam< double >( d, names::tau_syn_in, tau_in_, node );

  if ( Tau_ <= 0 || tau_ex_ <= 0 || tau_in_ <= 0 )
  {
    throw BadProperty( "Membrane and synapse time constants must be strictly positive." );
  }

  return 5.;
}

void
nest::iaf_matco_2018::State_::get( DictionaryDatum& d, const Parameters_& p ) const
{
  def< double >( d, names::V_m, V_m_); // Membrane potential
  def< double >( d, names::omega, omega_); //
  def< bool >( d, names::phi, phi_); 
}

void
nest::iaf_matco_2018::State_::set( const DictionaryDatum& d, const Parameters_& p, double delta_EL, Node* node )
{
  updateValueParam< double >( d, names::V_m, V_m_, node ); 
}

nest::iaf_matco_2018::Buffers_::Buffers_( iaf_matco_2018& n )
  : logger_( n )
{
}

nest::iaf_matco_2018::Buffers_::Buffers_( const Buffers_&, iaf_matco_2018& n )
  : logger_( n )
{
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

nest::iaf_matco_2018::iaf_matco_2018()
  : ArchivingNode()
  , P_()
  , S_()
  , B_( *this )
{
  recordablesMap_.create();
}

nest::iaf_matco_2018::iaf_matco_2018( const iaf_matco_2018& n )
  : ArchivingNode( n )
  , P_( n.P_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
nest::iaf_matco_2018::init_buffers_()
{
  B_.spike_exc_.clear(); // includes resize
  B_.spike_inh_.clear(); // includes resize

  B_.input_buffer_.clear(); // includes resize
  B_.logger_.reset();
  ArchivingNode::clear_history();
}

void
nest::iaf_matco_2018::pre_run_hook()
{
  // ensures initialization in case mm connected after Simulate
  B_.logger_.init();


  //V_.PSConInit_E = 1.0 * numerics::e / P_.tau_synE;
  //V_.PSConInit_I = 1.0 * numerics::e / P_.tau_synI;

  const double h = Time::get_resolution().get_ms();

  V_.rng_ = get_vp_specific_rng( get_thread() );
}

void
nest::iaf_matco_2018::update( const Time& origin, const long from, const long to )
{
  assert( to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
  assert( from < to );

  const double h = Time::get_resolution().get_ms();


  // evolve from timestep 'from' to timestep 'to' with steps of h each
  for ( long lag = from; lag < to; ++lag )
  {
   
    //double t_trig = Time( Time::step( kernel().simulation_manager.get_slice_origin().get_steps() + to ) ).get_ms();
    //kernel().connection_manager.force_update_weight( get_node_id(), lag );
    
    // get read access to the correct input-buffer slot
    const index input_buffer_slot = kernel().event_delivery_manager.get_modulo( lag );
    auto& input = B_.input_buffer_.get_values_all_channels( input_buffer_slot );

    // the spikes arriving at T+1 have an immediate effect on the state of the
    // neuron

    // add incoming spike
    //S_.omega_ += B_.spike_exc_.get_value( lag ) / P_.tau_ex_;
    //S_.omega_ += B_.spike_inh_.get_value( lag ) / P_.tau_in_;

    //TODO: fix for inhibitory
    S_.V_m_ += (-S_.V_m_ + P_.k1_ * (B_.spike_exc_.get_value( lag ) + B_.spike_inh_.get_value( lag ))) / P_.tau_ex_;

    if ( ((S_.V_m_) - P_.alpha_ * S_.omega_) > P_.Theta_) //  threshold crossing
    {
      set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );

      SpikeEvent se;
      kernel().event_delivery_manager.send( *this, se, lag );
      S_.phi_ = true;
    }

    else
    {
      S_.phi_ = false;
    }

    S_.omega_ += (-S_.omega_ + S_.phi_) / P_.Tau_;

    // reset all values in the currently processed input-buffer slot
    B_.input_buffer_.reset_values_all_channels( input_buffer_slot );

    // log state data
    B_.logger_.record_data( origin.get_steps() + lag );
  }
}

void
nest::iaf_matco_2018::handle( SpikeEvent& e )
{
  assert( e.get_delay_steps() > 0 );

  if ( e.get_weight() > 0.0 )
  {
    B_.spike_exc_.add_value( e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ),
      e.get_weight() * e.get_multiplicity() );
  }

  else
  {
    B_.spike_inh_.add_value( e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ),
      -e.get_weight() * e.get_multiplicity() );
  }

}

void
nest::iaf_matco_2018::handle( CurrentEvent& e )
{
  assert( e.get_delay_steps() > 0 );

  const double c = e.get_current();
  const double w = e.get_weight();

  const index input_buffer_slot = kernel().event_delivery_manager.get_modulo(
    e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ) );

  // add weighted current; HEP 2002-10-04
  if ( 0 == e.get_rport() )
  {
    B_.input_buffer_.add_value( input_buffer_slot, Buffers_::I0, w * c );
  }
  if ( 1 == e.get_rport() )
  {
    B_.input_buffer_.add_value( input_buffer_slot, Buffers_::I1, w * c );
  }
}

void
nest::iaf_matco_2018::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}
