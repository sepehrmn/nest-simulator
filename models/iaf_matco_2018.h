/*
 *  iaf_matco_2018.h
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

#ifndef IAF_MATCO_2018_H
#define IAF_MATCO_2018_H

// Includes from nestkernel:
#include "archiving_node.h"
#include "connection.h"
#include "event.h"
#include "nest_types.h"
#include "recordables_map.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"

namespace nest
{

/* BeginUserDocs: neuron, integrate-and-fire, current-based

Short description
+++++++++++++++++

integrate-and-fire neuron model with rate-based adaptation

Description
+++++++++++

iaf_matco_2018 is an implementation of an integrate-and-fire model
according to [1]_.

The threshold crossing leads to a spike and an adaptation

Remarks:

The present implementation uses individual variables for the
components of the state vector and the non-zero matrix elements of
the propagator. Because the propagator is a lower triangular matrix,
no full matrix multiplication needs to be carried out and the
computation can be done "in place", i.e. no temporary state vector
object is required.

The template support of recent C++ compilers enables a more succinct
formulation without loss of runtime performance already at minimal
optimization levels. A future version of iaf_matco_2018 will probably
address the problem of efficient usage of appropriate vector and
matrix objects.

.. note::

  If `tau_m` is very close to `tau_syn_ex` or `tau_syn_in`, the model
  will numerically behave as if `tau_m` is equal to `tau_syn_ex` or
  `tau_syn_in`, respectively, to avoid numerical instabilities.

  For implementation details see the
  `IAF_neurons_singularity <../model_details/IAF_neurons_singularity.ipynb>`_ notebook.

iaf_matco_2018 can handle current input in two ways: Current input
through receptor_type 0 are handled as stepwise constant current
input as in other iaf models, i.e., this current directly enters
the membrane potential equation. Current input through
receptor_type 1, in contrast, is filtered through an exponential
kernel with the time constant of the excitatory synapse,
tau_syn_ex. For an example application, see [4]_.

For conversion between postsynaptic potentials (PSPs) and PSCs,
please refer to the ``postsynaptic_potential_to_current`` function in
:doc:`PyNEST Microcircuit: Helper Functions <../auto_examples/Potjans_2014/helpers>`.

Parameters
++++++++++

The following parameters can be set in the status dictionary.

===========  =======  ========================================================
 tau_m        ms      Membrane time constant
 tau_V_th     
 V_m          mV      Membrane potential in mV
 V_th         mV      Spike threshold in mV
 I_e          pA      Constant input current
 t_spike      ms      Point in time of last spike
===========  =======  ========================================================

References
++++++++++

.. [1] Tsodyks M, Uziel A, Markram H (2000). Synchrony generation in recurrent
       networks with frequency-dependent synapses. The Journal of Neuroscience,
       20,RC50:1-5. URL: https://infoscience.epfl.ch/record/183402
.. [2] Rotter S,  Diesmann M (1999). Exact simulation of
       time-invariant linear systems with applications to neuronal
       modeling. Biologial Cybernetics 81:381-402.
       DOI: https://doi.org/10.1007/s004220050570
.. [3] Diesmann M, Gewaltig M-O, Rotter S, & Aertsen A (2001). State
       space analysis of synchronous spiking in cortical neural
       networks. Neurocomputing 38-40:565-571.
       DOI: https://doi.org/10.1016/S0925-2312(01)00409-X
.. [4] Schuecker J, Diesmann M, Helias M (2015). Modulated escape from a
       metastable state driven by colored noise. Physical Review E 92:052119
       DOI: https://doi.org/10.1103/PhysRevE.92.052119

Sends
+++++

SpikeEvent

Receives
++++++++

SpikeEvent, CurrentEvent, DataLoggingRequest

See also
++++++++

matco_synapse

EndUserDocs */

class iaf_matco_2018 : public ArchivingNode
{

public:
  iaf_matco_2018();
  iaf_matco_2018( const iaf_matco_2018& );

  /**
   * Import sets of overloaded virtual functions.
   * @see Technical Issues / Virtual Functions: Overriding, Overloading, and
   * Hiding
   */
  using Node::handle;
  using Node::handles_test_event;

  port send_test_event( Node&, rport, synindex, bool ) override;

  void handle( SpikeEvent& ) override;
  void handle( CurrentEvent& ) override;
  void handle( DataLoggingRequest& ) override;

  port handles_test_event( SpikeEvent&, rport ) override;
  port handles_test_event( CurrentEvent&, rport ) override;
  port handles_test_event( DataLoggingRequest&, rport ) override;

  void get_status( DictionaryDatum& ) const override;
  void set_status( const DictionaryDatum& ) override;

  inline double
  get_V_m() const
  {
    return S_.V_m_;
  }

  inline double
  get_phi() const
  {
    return S_.phi_;
  }

private:
  void init_buffers_() override;
  void pre_run_hook() override;

  void update( const Time&, const long, const long ) override;

  // intensity function
  //bool phi_() const;

  // The next two classes need to be friends to access the State_ class/member
  friend class RecordablesMap< iaf_matco_2018 >;
  friend class UniversalDataLogger< iaf_matco_2018 >;

  // ----------------------------------------------------------------

  /**
   * Independent parameters of the model.
   */
  struct Parameters_
  {
    /** adaptation time constant in ms. */
    double Tau_;

    /** Membrane time constant in ms. */
    double k1_;

    /** External current in pA */
    double I_e_;

    /** Absolute threshold value **/
    double Theta_;

    /** Time constant of excitatory synaptic current in ms. */
    double tau_ex_;

    /** Time constant of inhibitory synaptic current in ms. */
    double tau_in_;

    double alpha_;

    Parameters_(); //!< Sets default parameter values

    void get( DictionaryDatum& ) const; //!< Store current values in dictionary

    /** Set values from dictionary.
     * @returns Change in reversal potential E_L, to be passed to State_::set()
     */
    double set( const DictionaryDatum&, Node* node );
  };

  // ----------------------------------------------------------------

  /**
   * State variables of the model.
   */
  struct State_
  {
    // state variables
    double i_syn_ex_; //!< postsynaptic current for exc. inputs, variable 1
    double i_syn_in_; //!< postsynaptic current for inh. inputs, variable 1
    double V_m_;      //!< membrane potential, variable 2

    /** Firing rate */
    double omega_;
    bool phi_;


    State_(); //!< Default initialization

    void get( DictionaryDatum&, const Parameters_& ) const;

    /** Set values from dictionary.
     * @param dictionary to take data from
     * @param current parameters
     * @param Change in reversal potential E_L specified by this dict
     */
    void set( const DictionaryDatum&, const Parameters_&, const double, Node* );
  };

  // ----------------------------------------------------------------

  /**
   * Buffers of the model.
   */
  struct Buffers_
  {
    Buffers_( iaf_matco_2018& );
    Buffers_( const Buffers_&, iaf_matco_2018& );

    /** buffers and sums up incoming spikes/currents */
    RingBuffer spike_exc_;
    RingBuffer spike_inh_;

    //! Indices for access to different channels of input_buffer_
    enum
    {
      SYN_IN = 0,
      SYN_EX,
      I0,
      I1,
      NUM_INPUT_CHANNELS
    };

    /** buffers and sums up incoming spikes/currents */
    MultiChannelInputBuffer< NUM_INPUT_CHANNELS > input_buffer_;

    //! Logger for all analog data
    UniversalDataLogger< iaf_matco_2018 > logger_;
  };

  // ----------------------------------------------------------------

  /**
   * Internal variables of the model.
   */
  struct Variables_
  {
    /** Amplitude of the synaptic current.
        This value is chosen such that a postsynaptic potential with
        weight one has an amplitude of 1 mV.
        @note mog - I assume this, not checked.
    */
    //    double PSCInitialValue_;

    double P20_;

    double weighted_spikes_ex_;
    double weighted_spikes_in_;


    RngPtr rng_; //!< random number generator of my own thread
  };

  // Access functions for UniversalDataLogger -------------------------------

  //! Read out the membrane potential
  //inline double
  //get_V_m_() const
  //{
  //  return S_.V_m_;
  //}



  inline double
  get_I_syn_ex_() const
  {
    return S_.i_syn_ex_;
  }

  inline double
  get_I_syn_in_() const
  {
    return S_.i_syn_in_;
  }

  // ----------------------------------------------------------------

  /**
   * @defgroup iaf_matco_2018
   * Instances of private data structures for the different types
   * of data pertaining to the model.
   * @note The order of definitions is important for speed.
   * @{
   */
  Parameters_ P_;
  State_ S_;
  Variables_ V_;
  Buffers_ B_;
  /** @} */

  //! Mapping of recordables names to access functions
  static RecordablesMap< iaf_matco_2018 > recordablesMap_;
};


inline port
nest::iaf_matco_2018::send_test_event( Node& target, rport receptor_type, synindex, bool )
{
  SpikeEvent e;
  e.set_sender( *this );
  return target.handles_test_event( e, receptor_type );
}

inline port
iaf_matco_2018::handles_test_event( SpikeEvent&, rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}

inline port
iaf_matco_2018::handles_test_event( CurrentEvent&, rport receptor_type )
{
  if ( receptor_type == 0 )
  {
    return 0;
  }
  else if ( receptor_type == 1 )
  {
    return 1;
  }
  else
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
}

inline port
iaf_matco_2018::handles_test_event( DataLoggingRequest& dlr, rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return B_.logger_.connect_logging_device( dlr, recordablesMap_ );
}

inline void
iaf_matco_2018::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  S_.get( d, P_ );
  ArchivingNode::get_status( d );

  ( *d )[ names::recordables ] = recordablesMap_.get_list();
}

inline void
iaf_matco_2018::set_status( const DictionaryDatum& d )
{
  Parameters_ ptmp = P_;                       // temporary copy in case of errors
  const double delta_EL = ptmp.set( d, this ); // throws if BadProperty
  State_ stmp = S_;                            // temporary copy in case of errors
  stmp.set( d, ptmp, delta_EL, this );         // throws if BadProperty

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  ArchivingNode::set_status( d );

  // if we get here, temporaries contain consistent set of properties
  P_ = ptmp;
  S_ = stmp;
}

} // namespace

#endif // iaf_matco_2018
