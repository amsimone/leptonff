! SPDX-License-Identifier: MIT
! Native API client. All amplitudes, phase space and FKS subtraction are McMule's.
program leptonff_mcmule_reference
  use mcmule
  use load_user, only: bin_kind
  use vegas_m, only: combine_quant
  use global_def, only: pcounter
  use integrands, only: integration_ps => ps
  use leptonff_mcmule_maps, only: fsr3, fsr4, fsr4_direct
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer, parameter :: obs_count=14
  real(kind=prec), target :: lower(obs_count)=-10001._prec, upper(obs_count)=1.01_prec
  real(kind=prec) :: values(obs_count,0:2), errors(obs_count,0:2)
  real(kind=prec) :: parent_Q, support_mass, xi, mu_factor, soft, coll, nts
  real(kind=prec) :: beta_i,beta_f,born,physical_alpha
  real(kind=prec), parameter :: moment_indices(5)=[1._prec,2._prec,3._prec,5._prec,10._prec]
  real(kind=prec), parameter :: cuts(9)=[.25_prec,.5_prec,.75_prec,.9_prec,.95_prec,.98_prec,.99_prec,.995_prec,.999_prec]
  integer :: calls, iterations, warm_calls, warm_iterations, seed, unit, i, j
  character(len=900) :: output,arg,mode,piece,angular_map='native'
  if(command_argument_count()/=15 .and. command_argument_count()/=16) stop 'Usage: output piece mode Q minmass calls iterations warmcalls warmiters seed xi mu_factor softcut collcut nts [native|fsr]'
  call get_command_argument(1,output)
  call get_command_argument(2,piece)
  call get_command_argument(3,mode)
  call get_command_argument(4,arg);read(arg,*)parent_Q
  call get_command_argument(5,arg);read(arg,*)support_mass
  call get_command_argument(6,arg);read(arg,*)calls
  call get_command_argument(7,arg);read(arg,*)iterations
  call get_command_argument(8,arg);read(arg,*)warm_calls
  call get_command_argument(9,arg);read(arg,*)warm_iterations
  call get_command_argument(10,arg);read(arg,*)seed
  call get_command_argument(11,arg);read(arg,*)xi
  call get_command_argument(12,arg);read(arg,*)mu_factor
  call get_command_argument(13,arg);read(arg,*)soft
  call get_command_argument(14,arg);read(arg,*)coll
  call get_command_argument(15,arg);read(arg,*)nts
  if(command_argument_count()==16)call get_command_argument(16,angular_map)
  if(trim(angular_map)/='native' .and. trim(angular_map)/='fsr' .and. trim(angular_map)/='fsr-direct')stop 'Invalid angular importance map'
  if(parent_Q<=1 .or. support_mass<0 .or. support_mass>=parent_Q .or. calls<2 .or. iterations<2 .or. &
     warm_calls<2 .or. warm_iterations<1 .or. seed<=0 .or. xi<=0 .or. xi>1 .or. mu_factor<=0 .or. &
     soft<=0 .or. coll<=0 .or. nts<=soft) stop 'Invalid reference configuration'
  if(trim(mode)/='shape' .and. trim(mode)/='inclusive') stop 'Invalid measurement mode'
  select case(trim(piece))
  case('ee2mm0','ee2mmFMM','ee2mmRMM','ee2mmRFMMMM','ee2mmRRMMMM')
  case default
    stop 'Only pure photonic muon-line pieces are supported'
  end select
  if(trim(mode)=='shape' .and. (trim(piece)=='ee2mm0' .or. trim(piece)=='ee2mmFMM')) &
    stop 'Born/no-real-photon shape is identically zero; do not run VEGAS on it'
  open(newunit=unit,file=trim(output),status='new',action='write');close(unit)
  call set_observable(obs_count,1,lower,upper,measurement,user_initialisation=initialize,user_integration=configure_mapping)
  bin_kind=1
  call runmcmule(warm_calls,warm_iterations,calls,iterations,seed,xi,xi,trim(piece),tfilename='checkpoint.vegas')
  call combine_quant(iterations,values,errors)
  if(any(.not.ieee_is_finite(values)) .or. any(.not.ieee_is_finite(errors))) stop 'Nonfinite reference result'
  if(any(values(:,0)/=0) .or. any(values(:,2)/=0)) stop 'Observable outside the one-bin integration range'
  beta_i=sqrt(1-4*Me**2/scms);beta_f=sqrt(1-4*Mm**2/scms)
  born=4*pi/(3*scms)*beta_f/beta_i*(1+2*Me**2/scms)*(1+2*Mm**2/scms)
  physical_alpha=1/137.035999084_prec
  open(newunit=unit,file=trim(output),status='replace',action='write')
  write(unit,'(a)') '{"schema":"leptonff-mcmule-fixed-order-v1",'
  write(unit,'(a)') '"piece":"'//trim(piece)//'","measurement":"'//trim(mode)//'",'
  write(unit,'(a)') '"angular_map":"'//trim(angular_map)//'",'
  write(unit,'(a,es25.17e3,a)') '"Q":',parent_Q,','
  write(unit,'(a,es25.17e3,a)') '"mass":',Mm/1000,','
  write(unit,'(a,es25.17e3,a)') '"electron_mass":',Me/1000,','
  write(unit,'(a,es25.17e3,a)') '"support_min_mass":',support_mass,','
  write(unit,'(a,es25.17e3,a)') '"physical_alpha":',physical_alpha,','
  write(unit,'(a,es25.17e3,a)') '"internal_alpha":',alpha,','
  write(unit,'(a,es25.17e3,a)') '"born_alpha_one_MeV_minus2":',born,','
  write(unit,'(a,es25.17e3,a)') '"xi":',xi,','
  write(unit,'(a,es25.17e3,a)') '"mu_factor":',mu_factor,','
  write(unit,'(a,es25.17e3,a)') '"softcut":',soft,','
  write(unit,'(a,es25.17e3,a)') '"collcut":',coll,','
  write(unit,'(a,es25.17e3,a)') '"nts_switch":',ntsSwitch,','
  write(unit,'(a,i0,a)') '"calls_per_iteration":',calls,','
  write(unit,'(a,i0,a)') '"iterations":',iterations,','
  write(unit,'(a,i0,a)') '"warmup_calls_per_iteration":',warm_calls,','
  write(unit,'(a,i0,a)') '"warmup_iterations":',warm_iterations,','
  write(unit,'(a,i0,a)') '"seed":',seed,','
  write(unit,'(a,i0,a)') '"bad_phase_space_points":',pcounter(2),','
  write(unit,'(a)') '"current":"pure vector","fsr":true,"isr":false,"ifi":false,"pairs":false,"weak":false,"running_alpha":false,'
  write(unit,'(a)') '"coupling":"fixed on-shell; upstream alpha=1, restore powers externally","observables":['
  do i=1,obs_count
    if(i<=5)then
      write(unit,'(a,es25.17e3,a)') '{"observable":"moment","value":',moment_indices(i),','
    else
      write(unit,'(a,es25.17e3,a)') '{"observable":"cumulative","value":',cuts(i-5),','
    endif
    write(unit,'(a,es25.17e3,a)') '"integral":',values(i,1),','
    write(unit,'(a,es25.17e3,a)') '"mc_error":',errors(i,1),'}'
    if(i<obs_count)write(unit,'(a)') ','
  enddo
  write(unit,'(a)') ']}'
  close(unit)
contains
  subroutine initialize
    call initflavour('mu-e',(1000*parent_Q)**2)
    musq=(1000*parent_Q*mu_factor)**2
    softcut=soft;collcut=coll;ntsSwitch=nts
    ! No soft photon exists in the virtual two-body contribution. Activating
    ! its near-to-soft real-emission branch would dereference a missing term.
    if(trim(piece)=='ee2mmFMM' .or. trim(piece)=='ee2mm0')ntsSwitch=-1
    nel=0;nmu=0;ntau=0;nhad=0
    pol1=0;pol2=0;pol3=0;pol4=0
    do j=1,obs_count
      write(names(j),'(a,i2.2)') 'obs',j
    enddo
  end subroutine

  subroutine configure_mapping(x,ndim)
    integer::ndim
    real(kind=prec)::x(ndim)
    if(trim(angular_map)=='native')return
    select case(trim(piece))
    case('ee2mmRMM','ee2mmRFMMMM')
      integration_ps=>fsr3
    case('ee2mmRRMMMM')
      integration_ps=>fsr4
      if(trim(angular_map)=='fsr-direct')integration_ps=>fsr4_direct
    end select
  end subroutine

  function measurement(p1,p2,p3,p4,p5,p6,p7) result(var)
    real(kind=prec),intent(in)::p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
    real(kind=prec)::var(obs_count),P(4),K(4),loss,z,power,total
    integer::k_index,n
    P=p1+p2
    if(trim(mode)=='inclusive')then
      z=sq(p3+p4)/sq(P)
      pass_cut=z>=(support_mass/parent_Q)**2
      userweight=1;var=1;return
    endif
    K=p5
    if(trim(piece)=='ee2mmRRMMMM')K=K+p6
    ! Photon loss is evaluated without subtracting two nearly equal masses.
    loss=(s(P,K)-sq(K))/sq(P)
    if(loss < -1.e-12_prec .or. loss>1+1.e-9_prec) stop 'Invalid radiative mass loss'
    loss=max(0._prec,min(1._prec,loss));z=1-loss
    pass_cut=.true.
    if(z<(support_mass/parent_Q)**2)pass_cut=.false.
    if(abs(sq(p3+p4)/sq(P)-z)>1.e-8_prec)stop 'Muon/photon observable definitions disagree'
    ! Positive energy-loss weight guides adaptation. Multiplying it by var
    ! restores F-1, separately on every real and FKS counterevent.
    userweight=loss
    if(loss==0)pass_cut=.false.
    do k_index=1,5
      n=int(moment_indices(k_index))-1
      power=1;total=0
      do while(n>0)
        total=total+power;power=power*z;n=n-1
      enddo
      var(k_index)=-total
    enddo
    do k_index=1,9
      var(k_index+5)=0
      if(z<cuts(k_index) .and. loss>0)var(k_index+5)=-1/loss
    enddo
  end function
end program
