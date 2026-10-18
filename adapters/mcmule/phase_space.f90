! SPDX-License-Identifier: MIT
! Importance coordinates only. Phase-space weights and soft configurations
! come from McMule; amplitudes and FKS operations are untouched.
module leptonff_mcmule_maps
  use global_def, only: prec, pi, scms, Mm
  use phase_space, only: psx3_fks, psx4_fkss, psx4_35_36_fkss, ksofts, youts
  use functions, only: make_mlm
  implicit none
contains
  subroutine logarithmic_angle(u,delta,v,jac)
    real(kind=prec),intent(in)::u,delta
    real(kind=prec),intent(out)::v,jac
    real(kind=prec)::t,h,d,sign
    if(u<.5_prec)then
      t=2*u;sign=1
    else
      t=2*u-1;sign=-1
    endif
    h=log(1+1/delta)
    d=delta*(exp(h*t)-1)
    v=(1+sign*(1-d))/2
    jac=delta*h*exp(h*t)
  end subroutine

  function rotation(u,phi) result(R)
    real(kind=prec),intent(in)::u,phi
    real(kind=prec)::R(4,4),c,s,a,b
    c=2*u-1;s=sqrt(max(0._prec,1-c*c));a=cos(2*pi*phi);b=sin(2*pi*phi)
    R=0
    R(1,1)=a*c;R(1,2)=-b;R(1,3)=a*s
    R(2,1)=b*c;R(2,2)=a;R(2,3)=b*s
    R(3,1)=-s;R(3,3)=c;R(4,4)=1
  end function

  subroutine fsr3(ra,p1,m1,p2,m2,q1,m3,q2,m4,q3,weight)
    real(kind=prec),intent(in)::ra(5),m1,m2,m3,m4
    real(kind=prec),intent(out)::p1(4),p2(4),q1(4),q2(4),q3(4),weight
    real(kind=prec)::x(5),R(4,4),jac,beta,delta,pairmass2
    x=ra
    pairmass2=scms*(1-ra(1))
    jac=1
    if(pairmass2>4*m3*m3 .and. m3==m4)then
      beta=sqrt(1-4*m3*m3/pairmass2)
      if(beta>1.e-8_prec)then
        delta=4*m3*m3/(pairmass2*beta*(1+beta))
        call logarithmic_angle(ra(5),delta,x(5),jac)
      endif
    endif
    ! Generate about the photon axis, then restore its full solid angle.
    x(2)=1;x(3)=0
    call psx3_fks(x,p1,m1,p2,m2,q1,m3,q2,m4,q3,weight)
    R=rotation(ra(2),ra(3))
    q1=matmul(R,q1);q2=matmul(R,q2);q3=matmul(R,q3)
    ksofts(1)=make_mlm(matmul(R,ksofts(1)%momentum))
    youts(1)=2*ra(2)-1
    weight=weight*jac
  end subroutine

  subroutine fsr4(ra,p1,m1,p2,m2,q1,m3,q2,m4,q3,q4,weight)
    real(kind=prec),intent(in)::ra(8),m1,m2,m3,m4
    real(kind=prec),intent(out)::p1(4),p2(4),q1(4),q2(4),q3(4),q4(4),weight
    real(kind=prec)::x(8),R(4,4),jac1,jac2,delta
    x=ra
    ! Fix the pair-rest-frame muon axis, sample both photon directions about
    ! it, then integrate the missing global solid angle by a common rotation.
    delta=2*m3*m3/scms
    call logarithmic_angle(ra(3),delta,x(3),jac1)
    call logarithmic_angle(ra(4),delta,x(4),jac2)
    x(7)=0;x(8)=1
    call psx4_fkss(x,p1,m1,p2,m2,q1,m3,q2,m4,q3,q4,weight)
    R=rotation(ra(8),ra(7))
    q1=matmul(R,q1);q2=matmul(R,q2);q3=matmul(R,q3);q4=matmul(R,q4)
    ksofts(1)=make_mlm(matmul(R,ksofts(1)%momentum))
    ksofts(2)=make_mlm(matmul(R,ksofts(2)%momentum))
    youts(1)=ksofts(1)%momentum(3)/ksofts(1)%momentum(4)
    youts(2)=ksofts(2)%momentum(3)/ksofts(2)%momentum(4)
    weight=weight*jac1*jac2
  end subroutine

  subroutine fsr4_direct(ra,p1,m1,p2,m2,q1,m3,q2,m4,q3,q4,weight)
    real(kind=prec),intent(in)::ra(8),m1,m2,m3,m4
    real(kind=prec),intent(out)::p1(4),p2(4),q1(4),q2(4),q3(4),q4(4),weight
    real(kind=prec)::x(8),jac1,jac2,delta,e1,e2
    integer::root
    x=ra;delta=2*m3*m3/scms
    call logarithmic_angle(ra(3),delta,x(3),jac1)
    call logarithmic_angle(ra(4),delta,x(4),jac2)
    root=1
    if(ra(8)>=.5_prec)root=2
    x(8)=2*ra(8)-(root-1)
    ! This upstream map uses its two mass arguments for both incoming and
    ! outgoing lines. Its phase-space Jacobian depends on the OUTGOING masses.
    ! Supply those masses, then restore the actual incoming electron momenta.
    ! Both kinematic roots must be sampled, with inverse probability two.
    call psx4_35_36_fkss(x,p1,m3,p2,m4,q1,q2,q3,q4,weight,root)
    e1=(scms+m1*m1-m2*m2)/(2*sqrt(scms));e2=sqrt(scms)-e1
    p1=[0._prec,0._prec,sqrt(e1*e1-m1*m1),e1]
    p2=[0._prec,0._prec,-sqrt(e2*e2-m2*m2),e2]
    youts(1)=ksofts(1)%momentum(3)/ksofts(1)%momentum(4)
    youts(2)=ksofts(2)%momentum(3)/ksofts(2)%momentum(4)
    weight=weight*2*jac1*jac2
  end subroutine
end module
