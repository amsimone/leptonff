! Installed ISO_C_BINDING interface; C headers define the ABI.
module leptonff
  use iso_c_binding
  implicit none
  integer(c_int),parameter :: LEPTONFF_PROCESS_W_PLUS=24, LEPTONFF_PROCESS_W_MINUS=-24, LEPTONFF_PROCESS_Z=23
  integer(c_int),parameter :: LEPTONFF_PHOTONIC_LL=1, LEPTONFF_PHOTONIC_NLL=2, LEPTONFF_LEPTONIC_NLL=3
  integer(c_int),parameter :: LEPTONFF_QED_HVP=4
  integer(c_int),parameter :: LEPTONFF_MUON_MINUS=0, LEPTONFF_MUON_PLUS=1, LEPTONFF_PHOTON=2
  integer(c_int),parameter :: LEPTONFF_FF_ORDER1=1, LEPTONFF_FF_COUPLED_NLL=2
  integer(c_int),parameter :: LEPTONFF_ALPHA_ON_SHELL=0, LEPTONFF_ALPHA_MSBAR=1
  integer(c_int),parameter :: LEPTONFF_ALPHA_FIXED=0, LEPTONFF_ALPHA_ONE_LOOP=1, LEPTONFF_ALPHA_TWO_LOOP=2
  integer(c_int),parameter :: LEPTONFF_ONE_FAMILY=0, LEPTONFF_LEPTON_THRESHOLDS=1
  type,bind(C) :: leptonff_fragmentation_config
    integer(c_size_t) :: struct_size
    integer(c_int) :: abi_version,model,alpha_scheme,running,flavours
    real(c_double) :: alpha,alpha_reference_scale,muon_mass,initial_scale,final_scale
    integer(c_size_t) :: grid_bins
    real(c_double) :: log_span,absolute_tolerance,relative_tolerance
  end type
  interface
    subroutine leptonff_fragmentation_config_default(c) bind(C)
      import leptonff_fragmentation_config
      type(leptonff_fragmentation_config),intent(out) :: c
    end subroutine
    integer(c_int) function leptonff_fragmentation_create(c,context) bind(C)
      import leptonff_fragmentation_config,c_ptr,c_int
      type(leptonff_fragmentation_config),intent(in) :: c
      type(c_ptr),intent(out) :: context
    end function
    subroutine leptonff_fragmentation_destroy(context) bind(C)
      import c_ptr
      type(c_ptr),value :: context
    end subroutine
    integer(c_int) function leptonff_fragmentation_moment(context,parent,observed,n,result) bind(C)
      import c_ptr,c_int,c_double
      type(c_ptr),value :: context
      integer(c_int),value :: parent,observed
      real(c_double),value :: n
      real(c_double),intent(out) :: result
    end function
    integer(c_int) function leptonff_fragmentation_integral(context,parent,observed,lower,upper,result) bind(C)
      import c_ptr,c_int,c_double
      type(c_ptr),value :: context
      integer(c_int),value :: parent,observed
      real(c_double),value :: lower,upper
      real(c_double),intent(out) :: result
    end function
    integer(c_int) function leptonff_fragmentation_convolve(context,parent,observed,x,f,data,result) bind(C)
      import c_ptr,c_funptr,c_int,c_double
      type(c_ptr),value :: context,data
      type(c_funptr),value :: f
      integer(c_int),value :: parent,observed
      real(c_double),value :: x
      real(c_double),intent(out) :: result
    end function
    function leptonff_last_error() result(message) bind(C)
      import c_ptr
      type(c_ptr) :: message
    end function
    function leptonff_version() result(message) bind(C)
      import c_ptr
      type(c_ptr) :: message
    end function
  end interface
  type,bind(C) :: leptonff_decay_config
    integer(c_size_t) :: struct_size
    integer(c_int) :: abi_version,process,model
    real(c_double) :: alpha,muon_mass,mu_frag
    integer(c_size_t) :: grid_bins
    real(c_double) :: log_span
    integer(c_int) :: hard_samples
    real(c_double) :: soft_cutoff,collinear_cutoff,minimum_z
    real(c_double) :: muon_pt_min,muon_eta_max,neutrino_pt_min,transverse_mass_min
    real(c_double) :: mass_min,mass_max,hadronic_matching_scale
  end type
  type,bind(C) :: leptonff_prediction
    real(c_double) :: born,order1,matched,order2,higher
  end type
  type,bind(C) :: leptonff_observable
    type(c_ptr) :: observable,edges
    integer(c_size_t) :: edge_count
  end type
  type,bind(C) :: leptonff_decay_diagnostics
    integer(c_int) :: abi_version,photonic_unitary_conditioning
    real(c_double) :: a_initial,source_number_moment,raw_inclusive_factor,wrap_bound,tail_bound
  end type
  interface
    subroutine leptonff_decay_config_default(c) bind(C)
      import leptonff_decay_config
      type(leptonff_decay_config),intent(out) :: c
    end subroutine
    integer(c_int) function leptonff_decay_create(c,observables,count,context) bind(C)
      import leptonff_decay_config,c_ptr,c_size_t,c_int
      type(leptonff_decay_config),intent(in) :: c
      type(c_ptr),value :: observables
      integer(c_size_t),value :: count
      type(c_ptr),intent(out) :: context
    end function
    integer(c_size_t) function leptonff_decay_output_count(context) bind(C)
      import c_ptr,c_size_t
      type(c_ptr),value :: context
    end function
    integer(c_int) function leptonff_decay_evaluate(context,first,second,seed,outputs,capacity) bind(C)
      import c_ptr,c_double,c_int,c_int64_t,c_size_t,leptonff_prediction
      type(c_ptr),value :: context
      real(c_double),intent(in) :: first(4),second(4)
      integer(c_int64_t),value :: seed
      type(leptonff_prediction),intent(out) :: outputs(*)
      integer(c_size_t),value :: capacity
    end function
    subroutine leptonff_decay_destroy(context) bind(C)
      import c_ptr
      type(c_ptr),value :: context
    end subroutine
    integer(c_int) function leptonff_decay_get_diagnostics(context,diagnostics) bind(C)
      import c_ptr,c_int,leptonff_decay_diagnostics
      type(c_ptr),value :: context
      type(leptonff_decay_diagnostics),intent(out) :: diagnostics
    end function
    integer(c_int) function leptonff_decay_response(context,parent,seed,ncos,nphi,outputs,capacity) bind(C)
      import c_ptr,c_double,c_int,c_int64_t,c_size_t
      type(c_ptr),value :: context
      real(c_double),intent(in) :: parent(4)
      integer(c_int64_t),value :: seed
      integer(c_int),value :: ncos,nphi
      real(c_double),intent(out) :: outputs(*)
      integer(c_size_t),value :: capacity
    end function
  end interface
end module leptonff
