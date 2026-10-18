! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2019-2022, respective authors of MCFM.
! Copyright (C) 2026 Simone Amoroso.
! Replacement analysis hook for MCFM 10.3.1. The host's bookplot routine
! accumulates real/dipole entries in temporary bins and squares their SUM.
! Evaluate each FF observable completely before passing its signed weight.
subroutine nplotter_user(pjet,wt,wt2,nd)
  use types
  use iso_c_binding, only: c_int,c_double
  implicit none
  include 'mxpart.f'
  include 'nplot.f'
  include 'nproc.f'
  include 'ptilde.f'
  include 'first.f'
  integer :: nplotmax
  common/nplotmax/nplotmax
  real(dp),intent(in) :: pjet(mxpart,4),wt,wt2
  integer,intent(in) :: nd
  interface
    integer(c_int) function evaluate(proc,p3,p4,ref3,ref4,r,s) bind(C,name='leptonff_mcfm_evaluate_all')
      import c_int,c_double
      integer(c_int),value :: proc
      real(c_double),intent(in) :: p3(4),p4(4),ref3(4),ref4(4)
      real(c_double),intent(out) :: r(33),s(8,17)
    end function
  end interface
  real(c_double) :: p3(4),p4(4),ref3(4),ref4(4),r(33),s(8,17)
  real(dp) :: weight,eta,pt
  integer :: n,tag,i,j,k,bin,charged,code,status
  character(len=16) :: photonic,production_only
  character(len=6) :: component
  character(len=8) :: observable
  character(len=20) :: spectrum_label
  character(len=2) :: number
  character(len=6),parameter :: components(8)=[character(len=6) :: &
    'born','nlo','match','order2','higher','psborn','photos','psdiff']
  character(len=20),parameter :: labels(33)=[character(len=20) :: &
    'ff_born','ff_expanded','ff_resummed','ff_matched','ff_shift','ff_beyond_order1','ff_order2', &
    'photos_born','photos_post','photos_shift','photos_minus_ff','photos_radiated', &
    'hard_projected','hard_geometry','hard_shift','decay_order1','decay_matched','photos_minus_decay','photos_minus_order1', &
    'boundary_dressing','projected_dressing','geometry_dressing','dressing_shift', &
    'dressed_matched','photos_minus_dressed','dressed_minus_order1', &
    'coupled_projected','coupled_geometry','coupled_shift','coupled_matched','coupled_order2', &
    'coupled_minus_order1','photos_minus_coupled']
  r=0._dp;s=0._dp;eta=1000._dp;pt=0._dp
  if(first)then
    tag=1
  else
    tag=2
    p3=pjet(3,:);p4=pjet(4,:)
    ref3=p3;ref4=p4
    if(nd>0)then
      ref3=ptilde(0,3,:);ref4=ptilde(0,4,:)
    endif
    call leptonff_archive_add(nproc,p3,p4,ref3,ref4,wt)
    call get_environment_variable('LEPTONFF_PRODUCTION_ONLY',production_only,status=status)
    if(status/=0 .or. trim(production_only)/='1')then
      code=evaluate(nproc,p3,p4,ref3,ref4,r,s)
      if(code/=0)error stop 'leptonff MCFM evaluation failed'
    endif
    charged=3
    if(nproc==1)charged=4
    pt=sqrt(pjet(charged,1)**2+pjet(charged,2)**2)
    if(pt>0._dp)eta=asinh(pjet(charged,3)/pt)
  endif
  n=nextnplot-1
  do i=1,33
    n=n+1
    weight=wt*r(i)
    call bookplot(n,tag,trim(labels(i)),.5_dp,weight,weight**2,0._dp,1._dp,1._dp,'lin')
    ! Eta is unchanged by longitudinal FF rescaling. For Z this is the first
    ! parent-direction entry per selected pair, not the both-muon histogram.
    if(i<=12)then
    n=n+1
    call bookplot(n,tag,trim(labels(i))//'_eta',eta,weight,weight**2,-2.4_dp,2.4_dp,.4_dp,'lin')
    endif
  enddo
  ! Separate one-bin histograms retain the host's real/dipole grouping for
  ! EVERY physical bin, including empty and overflow bins. No Born eta proxy.
  call get_environment_variable('LEPTONFF_PHOTONIC_LL',photonic,status=status)
  do k=1,17
    if(nproc==31)then
      if(k<=8)then
        observable='mass';bin=k
      elseif(k<=13)then
        observable='minpt';bin=k-8
      else
        observable='maxeta';bin=k-13
      endif
    else
      if(k<=6)then
        observable='pt';bin=k
      elseif(k<=13)then
        observable='mt';bin=k-6
      else
        observable='abseta';bin=k-13
      endif
    endif
    write(number,'(i2.2)')bin
    do j=1,8
      component=components(j)
      if(status==0 .and. trim(photonic)=='1' .and. j==4)component='ff2'
      spectrum_label='sp_'//trim(observable)//number//'_'//trim(component)
      n=n+1
      weight=wt*s(j,k)
      call bookplot(n,tag,trim(spectrum_label),.5_dp,weight,weight**2,0._dp,1._dp,1._dp,'lin')
    enddo
  enddo
  if(.not.first)call leptonff_cov_add(wt,r,s)
  ! The dispatcher resets nextnplot on EVERY call, not only during booking.
  nextnplot=n+1
  if(first)then
    first=.false.
    nplotmax=n
  endif
end subroutine

! Optional read-only diagnostic called by an isolated, instrumented storage
! object immediately BEFORE MCFM resets its raw iteration histograms.
! Preserve real/dipole grouping and include zero-valued iterations.
subroutine leptonff_iteration_audit(stage,part,ips,iteration,ncall)
  use types
  use iso_fortran_env, only: int64
  use MCFMStorage, only: masterStorage
  implicit none
  integer,intent(in) :: stage,part,ips,iteration
  integer(int64),intent(in) :: ncall
  integer :: nplotmax
  common/nplotmax/nplotmax
  integer :: k,unit
  real(dp) :: value,variance,calls
  character(len=:),allocatable :: title
  open(newunit=unit,file='iteration_rates.csv',status='unknown',position='append',action='write')
  calls=real(ncall,dp)
  do k=1,nplotmax
    associate(h=>masterStorage(part,ips)%histCentral%histos(k))
      if(.not.h%initialized())cycle
      if(h%nbins/=1)cycle
      title=h%title
      if(index(title,'ff_')/=1 .and. index(title,'hard_')/=1 .and. &
         index(title,'decay_')/=1 .and. index(title,'photos_')/=1 .and. &
         index(title,'boundary_')/=1 .and. index(title,'projected_')/=1 .and. &
         index(title,'geometry_')/=1 .and. index(title,'dressing_')/=1 .and. &
         index(title,'dressed_')/=1 .and. index(title,'coupled_')/=1 .and. index(title,'sp_')/=1)cycle
      value=h%xx(1)
      variance=(calls*h%xxsq(1)-value**2)/(calls-1)
      write(unit,'(4(i0,","),i0,",",a,2(",",es25.17))') &
        stage,part,ips,iteration,ncall,title,value,variance
    end associate
  enddo
  close(unit)
  call leptonff_cov_flush(stage,part,ips,iteration,ncall)
end subroutine
