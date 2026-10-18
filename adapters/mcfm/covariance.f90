! Optional event covariance for the single-thread, central-PDF NLO W/Z pilot.
! begin/end surround ONE integrator call. Real/dipole sums are accepted only
! at the host's actual temporary-histogram commit site, after all dipoles.
module leptonff_covariance
  use types
  implicit none
  integer,parameter :: dimension=144
  real(dp),save :: event(dimension)=0._dp,sums(dimension)=0._dp,products(dimension,dimension)=0._dp
  logical,save :: committed=.false.,enabled=.false.,initialized=.false.,photonic=.false.
contains
  subroutine configure()
    character(len=16) :: value
    integer :: status
    if(initialized)return
    call get_environment_variable('LEPTONFF_SPECTRA',value,status=status)
    enabled=status==0 .and. trim(value)=='1'
    if(enabled)then
      call get_environment_variable('OMP_NUM_THREADS',value,status=status)
      if(status/=0 .or. trim(value)/='1')error stop 'Spectrum covariance requires one thread'
    endif
    call get_environment_variable('LEPTONFF_PHOTONIC_LL',value,status=status)
    photonic=status==0 .and. trim(value)=='1'
    initialized=.true.
  end subroutine
end module

subroutine leptonff_cov_begin()
  use leptonff_covariance
  call configure()
  call leptonff_archive_begin()
  committed=.false.
  if(.not.enabled)return
  event=0._dp
end subroutine

subroutine leptonff_cov_add(weight,r,s)
  use leptonff_covariance
  implicit none
  real(dp),intent(in) :: weight,r(33),s(8,17)
  integer :: indices(8)
  integer,parameter :: coupled_indices(8)=[1,16,30,31,32,8,9,33]
  integer,parameter :: photonic_indices(8)=[1,16,24,7,26,8,9,25]
  integer :: i,j
  if(.not.enabled)return
  indices=coupled_indices
  if(photonic)indices=photonic_indices
  do i=1,8
    event(i)=event(i)+weight*r(indices(i))
  enddo
  do j=1,17
    do i=1,8
      event(8+8*(j-1)+i)=event(8+8*(j-1)+i)+weight*s(i,j)
    enddo
  enddo
end subroutine

subroutine leptonff_cov_real_commit()
  use leptonff_covariance
  committed=.true.
end subroutine

subroutine leptonff_cov_end()
  use leptonff_covariance
  implicit none
  include 'kpart.f'
  integer :: i,j,a,b,active(dimension),count,accepted
  accepted=1
  if(kpart==kreal .and. .not.committed)accepted=0
  call leptonff_archive_end(accepted)
  if(.not.enabled)return
  if(kpart==kreal .and. .not.committed)return
  sums=sums+event
  count=0
  do i=1,dimension
    if(event(i)==0._dp)cycle
    count=count+1;active(count)=i
  enddo
  do b=1,count
    j=active(b)
    do a=1,b
      i=active(a)
      products(i,j)=products(i,j)+event(i)*event(j)
    enddo
  enddo
end subroutine

subroutine leptonff_cov_flush(stage,part,ips,iteration,ncall)
  use leptonff_covariance
  use MCFMStorage, only: masterStorage
  use iso_fortran_env, only: int64
  implicit none
  integer,intent(in) :: stage,part,ips,iteration
  integer(int64),intent(in) :: ncall
  integer :: nplotmax
  common/nplotmax/nplotmax
  integer :: i,j,k,unit,idx,count
  character(len=32) :: labels(dimension)
  character(len=20) :: rates(8)
  character(len=20),parameter :: coupled_rates(8)=[character(len=20) :: &
    'ff_born','decay_order1','coupled_matched','coupled_order2','coupled_minus_order1', &
    'photos_born','photos_post','photos_minus_coupled']
  real(dp) :: calls,covariance,scale
  call leptonff_archive_flush(stage,part,ips,iteration,ncall)
  if(.not.enabled)return
  if(stage<1)error stop 'Event covariance excludes warm-up stages'
  rates=coupled_rates
  if(photonic)rates=[character(len=20) :: 'ff_born','decay_order1','dressed_matched','ff_order2', &
    'dressed_minus_order1','photos_born','photos_post','photos_minus_dressed']
  labels='';count=0
  do k=1,nplotmax
    associate(h=>masterStorage(part,ips)%histCentral%histos(k))
      if(.not.h%initialized())cycle
      if(h%nbins/=1)cycle
      idx=0
      do i=1,8
        if(h%title==trim(rates(i)))idx=i
      enddo
      if(index(h%title,'sp_')==1)then
        count=count+1;idx=8+count
      endif
      if(idx==0)cycle
      if(idx>dimension)error stop 'Too many covariance channels'
      labels(idx)=h%title
      scale=max(1._dp,abs(h%xx(1)),abs(sums(idx)))
      if(abs(h%xx(1)-sums(idx))>1e-8_dp*scale)error stop 'Event covariance sum differs from host'
      scale=max(1._dp,abs(h%xxsq(1)),abs(products(idx,idx)))
      if(abs(h%xxsq(1)-products(idx,idx))>1e-8_dp*scale)error stop 'Event covariance diagonal differs from host'
    end associate
  enddo
  if(any(labels==''))error stop 'Missing covariance channels'
  open(newunit=unit,file='event_covariance.csv',status='unknown',position='append',action='write')
  calls=real(ncall,dp)
  do j=1,dimension
    do i=1,j
      covariance=(calls*products(i,j)-sums(i)*sums(j))/(calls-1)
      write(unit,'(4(i0,","),i0,",",a,",",a,",",es25.17)') &
        stage,part,ips,iteration,ncall,trim(labels(i)),trim(labels(j)),covariance
    enddo
  enddo
  close(unit)
  sums=0._dp;products=0._dp
end subroutine
