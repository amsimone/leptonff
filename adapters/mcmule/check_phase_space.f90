! SPDX-License-Identifier: MIT
program check_phase_space
  use mcmule
  use global_def, only: initparameters
  use phase_space, only: psx3_fks,psx4_fkss
  use leptonff_mcmule_maps, only: fsr3,fsr4,fsr4_direct
  implicit none
  integer::n,i,j,k,seed_size
  integer,allocatable::seed_values(:)
  real(kind=prec)::x(8),p1(4),p2(4),q1(4),q2(4),q3(4),q4(4),w,v,total,total2,mean,error,Q
  real(kind=prec)::shell,closure,difference(4)
  character(len=100)::arg
  n=1000000;Q=3
  if(command_argument_count()>0)then
    call get_command_argument(1,arg);read(arg,*)n
  endif
  if(command_argument_count()>1)then
    call get_command_argument(2,arg);read(arg,*)Q
  endif
  call initparameters
  call initflavour('mu-e',(1000*Q)**2)
  call random_seed(size=seed_size);allocate(seed_values(seed_size))
  seed_values=731999;call random_seed(put=seed_values)
  write(*,'(a)') '['
  do k=1,5
    total=0;total2=0;shell=0;closure=0
    do i=1,n
      call random_number(x)
      select case(k)
      case(1)
        call psx3_fks(x,p1,Me,p2,Me,q1,Mm,q2,Mm,q3,w)
      case(2)
        call fsr3(x,p1,Me,p2,Me,q1,Mm,q2,Mm,q3,w)
      case(3)
        call psx4_fkss(x,p1,Me,p2,Me,q1,Mm,q2,Mm,q3,q4,w)
      case(4)
        call fsr4(x,p1,Me,p2,Me,q1,Mm,q2,Mm,q3,q4,w)
      case(5)
        call fsr4_direct(x,p1,Me,p2,Me,q1,Mm,q2,Mm,q3,q4,w)
      end select
      v=w*x(1)
      if(k>=3)v=v*x(2)
      total=total+v;total2=total2+v*v
      if(w==0)cycle
      difference=p1+p2-q1-q2-q3
      if(k>=3)difference=difference-q4
      closure=max(closure,maxval(abs(difference))/sqrt(scms))
      shell=max(shell,abs(sq(q1)-Mm*Mm)/scms,abs(sq(q2)-Mm*Mm)/scms,abs(sq(q3))/scms)
      if(k>=3)shell=max(shell,abs(sq(q4))/scms)
    enddo
    mean=total/n;error=sqrt(max(0._prec,total2/n-mean*mean)/(n-1))
    write(*,'(a,i0,a,es25.17e3,a,es25.17e3,a,es25.17e3,a,es25.17e3,a)') &
      '{"map_index":',k,',"volume":',mean,',"mc_error":',error,',"mass_shell_residual":',shell, &
      ',"momentum_residual":',closure,'}'
    if(k<5)write(*,'(a)') ','
  enddo
  write(*,'(a)') ']'
end program
