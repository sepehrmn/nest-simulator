
/*
 *  bpid_kp_2017.cpp
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

#include "bpid_kp_2017.h"

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

nest::RecordablesMap< nest::bpid_kp_2017 >
  nest::bpid_kp_2017::recordablesMap_;

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
  RecordablesMap< bpid_kp_2017 >::create()
  {
    // use standard names whereever you can for consistency!
    insert_( names::receptive_field, &bpid_kp_2017::get_receptive_field_ );
    insert_( names::contextual_field, &bpid_kp_2017::get_contextual_field_ );
    insert_( names::theta, &bpid_kp_2017::get_theta_ );
  }
} //namespace


  /* ----------------------------------------------------------------
   * Default and copy constructor for node
   * ---------------------------------------------------------------- */

   nest::bpid_kp_2017::bpid_kp_2017()
     : Archiving_Node()
     , P_()
     , S_()
     , B_( *this )
   {
     recordablesMap_.create();
   }

   nest::bpid_kp_2017::bpid_kp_2017( const bpid_kp_2017& n )
     : Archiving_Node( n )
     , P_( n.P_ )
     , S_( n.S_ )
     , B_( n.B_, *this )
   {

   }

   /* ----------------------------------------------------------------
    * Default constructors defining default parameters and state
    * ---------------------------------------------------------------- */

   nest::bpid_kp_2017::Parameters_::Parameters_()
     : k1_( 1.0 )
     , k2_( 1.0 )
     , integration_type_(ADDITIVE)
     , phis_()
     , alpha_( 0.01 )
     , beta_( 0.01 )
     , eta_( 0.01 )
     , interval_( 1.0 ) // ms
   {
   }

   nest::bpid_kp_2017::State_::State_()
     : E( 0.01 )
     , E_c_( 0.01 )
     , E_r_( 0.01 )
     , theta_( 0.0 )
     ,  w_0_ ( 0.0 )
     ,  v_0_ ( 0.0 )
     , receptive_field_ ( 0.0 )
     , contextual_field_ ( 0.0 )
   {
   }

 /* ----------------------------------------------------------------
  * Parameter and state extractions and manipulation functions
  * ---------------------------------------------------------------- */

void
nest::bpid_kp_2017::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  S_.get( d );
  Archiving_Node::get_status( d );

  DictionaryDatum receptor_types = new Dictionary();

  // For easy assignment, returning the enum value for RF and CF
  ( *receptor_types )[ "receptive_field" ] = RF;
  ( *receptor_types )[ "contextual_field" ] = CF;
  ( *d )[ "receptor_types" ] = receptor_types;

  DictionaryDatum integration_types = new Dictionary();
  ( *integration_types )[ "additive" ] = ADDITIVE;
  ( *integration_types )[ "modulatory" ] = MODULATORY;
  ( *d )[ "integration_types" ] = integration_types;

  ( *d )[ names::recordables ] = recordablesMap_.get_list();

  ( *d )[ names::receptive_field ] = S_.receptive_field_;
  ( *d )[ names::contextual_field ] = S_.contextual_field_;
}

void
nest::bpid_kp_2017::set_status( const DictionaryDatum& d )
{
  P_.set( d );
  S_.set( d );
}

// Parameters

void
nest::bpid_kp_2017::Parameters_::get( DictionaryDatum& d ) const
{
  def< double >( d, names::alpha, alpha_ );
  def< double >( d, names::beta, beta_ );
  def< double >( d, names::eta, eta_ );
  def< double >( d, names::k1, k1_ );
  def< double >( d, names::k2, k2_ );
  def< long >( d, names::integration_type, integration_type_ ); //TODO fix from long. int does not work
  def< double >( d, names::interval, interval_ );
  def< std::vector<double> >( d, names::phis, phis_ );
}

void
nest::bpid_kp_2017::Parameters_::set( const DictionaryDatum& d )
{
  updateValue< double >( d, names::alpha, alpha_ );
  updateValue< double >( d, names::beta, beta_ );
  updateValue< double >( d, names::eta, eta_ );
  updateValue< double >( d, names::k1, k1_ );
  updateValue< double >( d, names::k2, k2_ );
  updateValue< long >( d, names::integration_type, integration_type_ ); //TODO fix from long. int does not work
  updateValue< double >( d, names::interval, interval_ );

  if ( !d->known( names::k1 ) or !d->known( names::k2 )  )
  {
    if (integration_type_ == ADDITIVE)
    {
      k1_ = 1.0;
      k2_ = 1.0;
    }

    else if (integration_type_ == MODULATORY)
    {
      k1_ = 0.5;
      k2_ = 2.0;
    }
  }

  updateValue< std::vector< double > >( d, names::phis, phis_ );

  // P_.phis_ =
  // getValue< std::vector< double > >( d, names::phis );

}

// State

void
nest::bpid_kp_2017::State_::get( DictionaryDatum& d) const
{
  def< double >( d, names::theta, theta_);
  def< double >( d, names::w_0, w_0_ );
  def< double >( d, names::v_0, v_0_ );
  def< double >( d, names::receptive_field, receptive_field_ );
  def< double >( d, names::contextual_field, contextual_field_ );
}

void
nest::bpid_kp_2017::State_::set( const DictionaryDatum& d )
{
  updateValue< double >( d, names::theta, theta_ );
  updateValue< double >( d, names::w_0, w_0_ );
  updateValue< double >( d, names::v_0, v_0_ );
  updateValue< double >( d, names::receptive_field, receptive_field_ );
  updateValue< double >( d, names::contextual_field, contextual_field_ );
}

// Buffers

nest::bpid_kp_2017::Buffers_::Buffers_( bpid_kp_2017& n )
  : logger_( n )
{
}

nest::bpid_kp_2017::Buffers_::Buffers_( const Buffers_&, bpid_kp_2017& n )
  : logger_( n )
{
}

/* ----------------------------------------------------------------
* Node initialization functions
* ---------------------------------------------------------------- */

void
nest::bpid_kp_2017::init_state_( const Node& proto )
{
  const bpid_kp_2017& pr = downcast< bpid_kp_2017 >( proto );
  S_ = pr.S_;
}

void
nest::bpid_kp_2017::init_buffers_()
{
  B_.spikes_rf_.clear(); // includes resize
  B_.spikes_cf_.clear(); // includes resize
  B_.logger_.reset();
  Archiving_Node::clear_history();
}

void
nest::bpid_kp_2017::calibrate()
{
  // ensures initialization in case mm connected after Simulate
  B_.logger_.init();
  V_.rng_ = kernel().rng_manager.get_rng( get_thread() );
}

/* ----------------------------------------------------------------
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

 void
 nest::bpid_kp_2017::learn_( )
 {

   //if P_.interval_ == V_.current_interval_:


 }

void
nest::bpid_kp_2017::update( const Time& origin, const long from, const long to )
{
  assert(
    to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
  assert( from < to );

  for ( long lag = from; lag < to; ++lag )
  {

    // update the input current
    // the buffer for incoming spikes for every time step contains the
    // difference
    // of the total input h with respect to the previous step, so sum them up
    S_.receptive_field_ = B_.spikes_rf_.get_value( lag ) - S_.w_0_;

    // update the input current
    // the buffer for incoming spikes for every time step contains the
    // difference
    // of the total input h with respect to the previous step, so sum them up
    S_.contextual_field_ = B_.spikes_cf_.get_value( lag ) - S_.v_0_;

    double activ_val = 0.0;

    if (P_.integration_type_ == ADDITIVE)
    {
      activ_val = P_.k1_ * S_.receptive_field_ + P_.k2_ * S_.contextual_field_;
    }

    else if (P_.integration_type_ == MODULATORY)
    {
      activ_val = S_.receptive_field_ *
                ( P_.k1_ + ( 1 - P_.k1_ ) * exp( P_.k2_ * S_.receptive_field_ * S_.contextual_field_ ) );
    }

    // Floating point comparison according to ''The art of computer programming by Knuth''-
    // https://stackoverflow.com/questions/17333/what-is-the-most-effective-way-for-float-and-double-comparison
    // NOTE: infinity and NAN are ignores the if and else if.

    assert(!std::isnan(activ_val) && !std::isinf(activ_val));

    if ((activ_val - 20.) > ( (fabs(activ_val) < fabs(20.) ? fabs(20.) : fabs(activ_val)) * pow(10,-6)))
    {
      S_.theta_  = 1.0;
    }

    else if ((-20 - activ_val) > ( (fabs(-20) < fabs(activ_val) ? fabs(activ_val) : fabs(-20.)) * pow(10,-6)))
    {
      S_.theta_  = 0.0;
    }

    else
    {
      S_.theta_ = 1 / (1 + exp(-activ_val));
    }

    // The if statement below makes sure that the neuron only fires once in a
    // discrete timestep (defined by the network min_delay)
    if ( from == 0 && lag == to-1 )
    {
      // Training
      learn_();

      // threshold crossing
      if ( V_.rng_->drand() < S_.theta_ )
      {
        SpikeEvent se;
        kernel().event_delivery_manager.send( *this, se, lag );
      }

    } // of if (update now)

    // log state data
    B_.logger_.record_data( origin.get_steps() + lag );
  }
}

void
nest::bpid_kp_2017::handle( SpikeEvent& e )
{
  assert( e.get_delay() > 0 );

  // String of the rport
  std::string str = std::to_string(e.get_rport());

  // the receptor number
  std::string num_string = str.substr(0, str.length() - 1);
  int rt_num = std::stoi (num_string);
  
  // Type of the receptor (last digit)
  std::string type_string = str.substr(str.length() - 1, str.length());
  int rt_type = std::stoi (type_string);

  if ( rt_type == RF )
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
nest::bpid_kp_2017::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}
