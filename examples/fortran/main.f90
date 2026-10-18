program decay_example
  use iso_c_binding
  implicit none
  type,bind(C) :: config_type
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
  type,bind(C) :: prediction
    real(c_double) :: born,order1,matched,order2,higher
  end type
  interface
    subroutine defaults(c) bind(C,name='leptonff_decay_config_default')
      import config_type
      type(config_type),intent(out) :: c
    end subroutine
    integer(c_int) function create(c,observables,count,context) bind(C,name='leptonff_decay_create')
      import config_type,c_ptr,c_size_t,c_int
      type(config_type),intent(in) :: c
      type(c_ptr),value :: observables
      integer(c_size_t),value :: count
      type(c_ptr),intent(out) :: context
    end function
    integer(c_size_t) function output_count(context) bind(C,name='leptonff_decay_output_count')
      import c_ptr,c_size_t
      type(c_ptr),value :: context
    end function
    integer(c_int) function evaluate(context,first,second,seed,outputs,capacity) bind(C,name='leptonff_decay_evaluate')
      import c_ptr,c_double,c_int,c_int64_t,c_size_t,prediction
      type(c_ptr),value :: context
      real(c_double),intent(in) :: first(4),second(4)
      integer(c_int64_t),value :: seed
      type(prediction),intent(out) :: outputs(*)
      integer(c_size_t),value :: capacity
    end function
    subroutine destroy(context) bind(C,name='leptonff_decay_destroy')
      import c_ptr
      type(c_ptr),value :: context
    end subroutine
  end interface
  type(config_type) :: config
  type(c_ptr) :: context
  type(prediction),allocatable :: outputs(:)
  integer(c_size_t) :: count
  integer(c_int) :: status
  real(c_double) :: mu(4),nu(4)
  call defaults(config)
  config%grid_bins=128;config%hard_samples=4
  if(config%struct_size/=c_sizeof(config))error stop 'C ABI layout mismatch'
  status=create(config,c_null_ptr,0_c_size_t,context)
  if(status/=0)error stop 'Decay context creation failed'
  count=output_count(context);allocate(outputs(count))
  mu=[50._c_double,40._c_double,0._c_double,30._c_double]
  nu=[50._c_double,-40._c_double,0._c_double,-30._c_double]
  status=evaluate(context,mu,nu,42_c_int64_t,outputs,count)
  if(status/=0)error stop 'Decay evaluation failed'
  print *, 'Born acceptance, matched response, bins:',outputs(1)%born,outputs(1)%matched,count-1
  call destroy(context)
end program
