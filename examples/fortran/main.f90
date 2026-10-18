program decay_example
  use leptonff
  implicit none
  type(leptonff_decay_config) :: config
  type(c_ptr) :: context
  type(leptonff_fragmentation_config) :: ff_config
  type(leptonff_decay_diagnostics) :: diagnostics
  type(c_ptr) :: ff_context
  real(c_double) :: moment
  type(leptonff_prediction),allocatable :: outputs(:)
  type(leptonff_prediction),allocatable :: batch(:)
  integer(c_size_t) :: count
  integer(c_int) :: status
  real(c_double) :: mu(4),nu(4)
  real(c_double) :: first(4,2),second(4,2)
  integer(c_int64_t) :: seeds(2)
  call leptonff_decay_config_default(config)
  config%process=LEPTONFF_PROCESS_W_PLUS
  config%grid_bins=128;config%hard_samples=4
  if(config%struct_size/=c_sizeof(config))error stop 'C ABI layout mismatch'
  status=leptonff_decay_create(config,c_null_ptr,0_c_size_t,context)
  if(status/=0)error stop 'Decay context creation failed'
  status=leptonff_decay_get_diagnostics(context,diagnostics)
  if(status/=0)error stop 'Decay diagnostics failed'
  if(abs(diagnostics%source_number_moment+.719311605774509_c_double)>1.e-12_c_double) &
    error stop 'Photonic source convention mismatch'
  count=leptonff_decay_output_count(context);allocate(outputs(count))
  mu=[50._c_double,40._c_double,0._c_double,30._c_double]
  nu=[50._c_double,-40._c_double,0._c_double,-30._c_double]
  status=leptonff_decay_evaluate(context,mu,nu,42_c_int64_t,outputs,count)
  if(status/=0)error stop 'Decay evaluation failed'
  print *, 'Born acceptance, matched response, bins:',outputs(1)%born,outputs(1)%matched,count-1
  first(:,1)=mu;first(:,2)=mu;second(:,1)=nu;second(:,2)=nu;seeds=42_c_int64_t
  allocate(batch(2*count))
  status=leptonff_decay_evaluate_batch(context,first,second,seeds,2_c_size_t,2_c_int,batch,2*count)
  if(status/=0)error stop 'Batch evaluation failed'
  if(batch(1)%matched/=outputs(1)%matched .or. batch(count+1)%matched/=outputs(1)%matched) &
    error stop 'Batch/scalar result mismatch'
  call leptonff_decay_destroy(context)
  call leptonff_fragmentation_config_default(ff_config)
  if(ff_config%struct_size/=c_sizeof(ff_config))error stop 'Fragmentation ABI layout mismatch'
  status=leptonff_fragmentation_create(ff_config,ff_context)
  if(status/=0)error stop 'Fragmentation context creation failed'
  status=leptonff_fragmentation_moment(ff_context,LEPTONFF_MUON_MINUS,LEPTONFF_MUON_MINUS,2._c_double,moment)
  if(status/=0)error stop 'Fragmentation moment failed'
  if(abs(moment-.9768684454148703_c_double)>1.e-10_c_double)error stop 'Fragmentation reference mismatch'
  print *, 'N=2 fragmentation moment:',moment
  call leptonff_fragmentation_destroy(ff_context)
end program
