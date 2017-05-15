
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

  /*
   * template specialization must be placed in namespace
   *
   * Override the create() method with one call to RecordablesMap::insert_()
   * for each quantity to be recorded.
   */
  template <>
  void
  RecordablesMap< iaf_psc_exp >::create()
  {
    // use standard names whereever you can for consistency!
    insert_( names::receptive_field, &binary_kp_1994::get_receptive_field_ );
    insert_( names::contextual_field, &binary_kp_1994::get_contextual_field_ );
    insert_( names::contextual_field, &binary_kp_1994::get_theta_ );
  }


  /* ----------------------------------------------------------------
   * Default constructors defining default parameters and state
   * ---------------------------------------------------------------- */


  nest::iaf_psc_exp::Parameters_::Parameters_()
    : k1_( 0.5 ),
      k2_( 2.0 ),
      k3_( 0.0 )
  {
  }

  nest::iaf_psc_exp::State_::State_()
    : theta_( 0.0 ),
      w_0 ( 0.0 ),
      v_0 ( 0.0 )
  {
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

  ( *d )[ names::receptive_field ] = V_.receptive_field_;
  ( *d )[ names::contextual_field ] = V_.contextual_field_;
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

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
nest::binary_kp_1994::init_state_( const Node& proto )
{
  const iaf_psc_exp& pr = downcast< iaf_psc_exp >( proto );
  S_ = pr.S_;
}

void
nest::binary_kp_1994::init_buffers_()
{
  B_.spikes_ex_.clear(); // includes resize
  B_.spikes_in_.clear(); // includes resize
  B_.currents_.clear();  // includes resize
  B_.logger_.reset();
  Archiving_Node::clear_history();
}

void
nest::binary_kp_1994::calibrate()
{
  B_.currents_.resize( 2 );
  // ensures initialization in case mm connected after Simulate
  B_.logger_.init();

  const double h = Time::get_resolution().get_ms();

  // numbering of state vaiables: i_0 = 0, i_syn_ = 1, V_m_ = 2

  // commented out propagators: forward Euler
  // needed to exactly reproduce Tsodyks network

  // these P are independent
  V_.P11ex_ = std::exp( -h / P_.tau_ex_ );
  // P11ex_ = 1.0-h/tau_ex_;

  V_.P11in_ = std::exp( -h / P_.tau_in_ );
  // P11in_ = 1.0-h/tau_in_;

  V_.P22_ = std::exp( -h / P_.Tau_ );
  // P22_ = 1.0-h/Tau_;

  // these are determined according to a numeric stability criterion
  V_.P21ex_ = propagator_32( P_.tau_ex_, P_.Tau_, P_.C_, h );
  V_.P21in_ = propagator_32( P_.tau_in_, P_.Tau_, P_.C_, h );

  // P21ex_ = h/C_;
  // P21in_ = h/C_;

  V_.P20_ = P_.Tau_ / P_.C_ * ( 1.0 - V_.P22_ );
  // P20_ = h/C_;

  // TauR specifies the length of the absolute refractory period as
  // a double in ms. The grid based iaf_psc_exp can only handle refractory
  // periods that are integer multiples of the computation step size (h).
  // To ensure consistency with the overall simulation scheme such conversion
  // should be carried out via objects of class nest::Time. The conversion
  // requires 2 steps:
  //     1. A time object r is constructed defining  representation of
  //        TauR in tics. This representation is then converted to computation
  //        time steps again by a strategy defined by class nest::Time.
  //     2. The refractory time in units of steps is read out get_steps(), a
  //        member function of class nest::Time.
  //
  // Choosing a TauR that is not an integer multiple of the computation time
  // step h will leed to accurate (up to the resolution h) and self-consistent
  // results. However, a neuron model capable of operating with real valued
  // spike time may exhibit a different effective refractory time.

  V_.RefractoryCounts_ = Time( Time::ms( P_.t_ref_ ) ).get_steps();
  // since t_ref_ >= 0, this can only fail in error
  assert( V_.RefractoryCounts_ >= 0 );
}

void
nest::binary_kp_1994::update( const Time& origin, const long from, const long to )
{
  assert(
    to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
  assert( from < to );

  // evolve from timestep 'from' to timestep 'to' with steps of h each
  for ( long lag = from; lag < to; ++lag )
  {
    if ( S_.r_ref_ == 0 ) // neuron not refractory, so evolve V
    {
      S_.V_m_ = S_.V_m_ * V_.P22_ + S_.i_syn_ex_ * V_.P21ex_
        + S_.i_syn_in_ * V_.P21in_ + ( P_.I_e_ + S_.i_0_ ) * V_.P20_;
    }
    else
    {
      --S_.r_ref_;
    } // neuron is absolute refractory

    // exponential decaying PSCs
    S_.i_syn_ex_ *= V_.P11ex_;
    S_.i_syn_in_ *= V_.P11in_;

    // add evolution of presynaptic input current
    S_.i_syn_ex_ += ( 1. - V_.P11ex_ ) * S_.i_1_;

    // the spikes arriving at T+1 have an immediate effect on the state of the
    // neuron

    V_.weighted_spikes_ex_ = B_.spikes_ex_.get_value( lag );
    V_.weighted_spikes_in_ = B_.spikes_in_.get_value( lag );

    S_.i_syn_ex_ += V_.weighted_spikes_ex_;
    S_.i_syn_in_ += V_.weighted_spikes_in_;

    if ( S_.V_m_ >= P_.Theta_ ) // threshold crossing
    {
      S_.r_ref_ = V_.RefractoryCounts_;
      S_.V_m_ = P_.V_reset_;

      set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );

      SpikeEvent se;
      kernel().event_delivery_manager.send( *this, se, lag );
    }

    // set new input current
    S_.i_0_ = B_.currents_[ 0 ].get_value( lag );
    S_.i_1_ = B_.currents_[ 1 ].get_value( lag );

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
