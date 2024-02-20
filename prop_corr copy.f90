program dynamics

implicit none
integer :: nst, inist, natom, ndim, n1, n2, i, j, k, n, nr, im, branch
integer ::  nc, mult, nfr, natomf
real (kind=8), allocatable, dimension (:) :: ElPhase
real (kind=8) CompTime, TimeStep, time, E0, E1, Mau, val, E_tot, dr
complex(8)    :: ii
real (kind=8), allocatable, dimension (:) :: R0,R1,P0,P1,V0,V1,F0,F1,Fb,M,Es0,Es1,Esb
real (kind=8), allocatable, dimension (:,:) :: Fs0,Fs1,Fsb
real (kind=8), allocatable, dimension (:,:,:) :: Cs0,Cs1,Csb
complex (kind=8), allocatable, dimension (:) :: A0,A1,Ab
complex (kind=8), allocatable, dimension (:,:) :: HE_0,HE_1,HE_b
character (len=1), allocatable, dimension (:) :: atom
character (len=50) f_inp
ii=(0.,1.)
Mau=1822.887
call getarg(1,f_inp)

open(1,file=trim(f_inp)//'.0')
read(1,*) natom, nst
read(1,*) branch
read(1,*) time, TimeStep
ndim=natom*3

!Allocation
allocate(r0(ndim),p0(ndim),r1(ndim),p1(ndim), atom(natom))

! Read initial point data
read(1,*)nc, mult
do i=1,natom
  read(1,*) atom(i),r0(i*3-2:i*3)
end do
read(1,*)
do i=1,natom
  read(1,*) p0(i*3-2:i*3)
end do

!Allocation
allocate (A0(nst), A1(nst), Ab(nst))
allocate (Es0(nst), Es1(nst),Esb(nst),Fs0(ndim,nst),Fs1(ndim,nst),Fsb(ndim,nst))
allocate (Cs0(ndim,nst,nst),Cs1(ndim,nst,nst),Csb(ndim,nst,nst))
allocate (HE_0(nst,nst), HE_1(nst,nst),HE_b(nst,nst))
allocate (F0(ndim),F1(ndim),Fb(ndim),V0(ndim),V1(ndim))
allocate (M(ndim),ElPhase(nst))

! Read amplitudes
read(1,*)
do i=1,nst
  read(1,*) A0(i)
end do

! Set masses 
do i=1,natom
  if (atom(i)=='C') then 
    M(i*3-2:i*3)=12*Mau
  else if (atom(i)=='N') then 
    M(i*3-2:i*3)=14*Mau
  else if (atom(i)=='H') then 
    M(i*3-2:i*3)=Mau
  else if (atom(i)=='D') then 
    M(i*3-2:i*3)=2*Mau
  else if (atom(i)=='F') then
    M(i*3-2:i*3)=19*Mau
  else
    Print*,' Atom ',atom(i), ' is not suppoted'
    stop
  end if
end do

! Read potential energies
read(1,*)
do n=1,nst
  read(1,*) val, i  
  Es0(i)=val
end do

! Read forces
read(1,*)
nr=nst*ndim
do n=1,nr
  read(1,*) val, i,j
  Fs0(j,i)=val
end do

! Read couplings
read(1,*)
nr=ndim*nst*(nst-1)/2
do n=1,nr
  read(1,*) val,i,k,j
  Cs0(j,i,k)=val
  Cs0(j,k,i)=-Cs0(j,i,k)
end do

close(1)

!Read preliminary final point data
open(1,file=trim(f_inp)//'.p')
read(1,*) natom, nst
read(1,*) branch
read(1,*) time, TimeStep

read(1,*)
do i=1,natom
  read(1,*) atom(i),r1(i*3-2:i*3)
end do
read(1,*)
do i=1,natom
  read(1,*) p1(i*3-2:i*3)
end do

read(1,*)
do i=1,nst
  read(1,*) A1(i)
end do

read(1,*)
do n=1,nst
  read(1,*) val, i  
  Es1(i)=val
end do

read(1,*)
nr=nst*ndim
do n=1,nr
    read(1,*) val, i,j
    Fs1(j,i)=val
end do

read(1,*)
nr=ndim*nst*(nst-1)/2
do n=1,nr
      read(1,*) val,i,k,j
      Cs1(j,i,k)=val
      Cs1(j,k,i)=-Cs1(j,i,k)
end do

close(1)

!Electronic phase
ElPhase(1)=1.
do j=2,nst
   val=sum(Cs0(:,1,j)*Cs1(:,1,j))/sqrt(sum(Cs0(:,1,j)**2)*sum(Cs1(:,1,j)**2))
   ElPhase(j)=sign(1.,val)
   if(abs(val).lt.0.5.and.abs(A1(j)).ge.0.35) print'(a30,i3,a17,f8.4)','!! Warring: the sign for state',j,' is not reliable!',val 
end do

do i=2,nst
  Cs1(:,i,:)=Cs1(:,i,:)*ElPhase(i)
  Cs1(:,:,i)=Cs1(:,:,i)*ElPhase(i)
end do
!print'(7f5.1)',ElPhase

! Velocities and Ehrenfest force
V0=P0/M
F0=CompForceEhr(A0,Fs0,Es0,Cs0,nst,ndim)
E0=sum(Es0*(abs(A0)**2))+sum(P0*V0)/2.
V1=P1/M
F1=CompForceEhr(A1,Fs1,Es1,Cs1,nst,ndim)
E1=sum(Es1*(abs(A1)**2))+sum(P1*V1)/2.




!Electronic Hamiltonian at initial point
do n1=1,nst
  HE_0(n1,n1)=Es0(n1)+77.67785291
  do n2=n1+1,nst
    HE_0(n1,n2)=-ii*sum(V0*Cs0(:,n1,n2))
    HE_0(n2,n1)=-HE_0(n1,n2)
  end do
end do




!Electronic Hamiltonian at preliminary point
  do n1=1,nst
    HE_1(n1,n1)=Es1(n1)+77.67785291
    do n2=n1+1,nst
      HE_1(n1,n2)=-ii*sum(V1*Cs1(:,n1,n2))
      HE_1(n2,n1)=-HE_1(n1,n2)
    end do
  end do


!Propagation
 Ab=matmul( magnus2( -ii*HE_0, -ii*HE_0, TimeStep/20), A0 )
  Esb=0.05*Es1+0.95*Es0
  Fsb=0.05*Fs1+0.95*Fs0  
  Csb=0.05*Cs1+0.95*Cs0    
  F1 = CompForceEhr(Ab,Fsb,Esb,Csb,nst,ndim)/10.
  do im=1,9
    HE_b=(im*HE_1+(10-im)*HE_0)*0.1
    Esb=(0.1*im+0.05)*Es1+(0.95-im*0.1)*Es0
    Fsb=(0.1*im+0.05)*Fs1+(0.95-im*0.1)*Fs0
    Csb=(0.1*im+0.05)*Cs1+(0.95-im*0.1)*Cs0
    A1 = matmul( magnus2( -ii*HE_b, -ii*HE_b, TimeStep/10.), Ab )
    Ab=A1
    Fb=CompForceEhr(A1,Fsb,Esb,Csb,nst,ndim)
    F1 = F1+Fb/10.
  end do
  A1 = matmul( magnus2( -ii*HE_1, -ii*HE_1, TimeStep/20.), Ab )   

  P1=P0+TimeStep * F1



!Fragmentation
nfr=0
natomf=natom
do i=1,natom
val=sqrt(sum(f0(3*i-2:3*i)**2))
if(val.lt.1.e-5) then
   nfr=i
   mult=mult-1
   if (mult==2) mult=4
   natomf=natom-1
   open(1,file=trim(f_inp)//'.diss',position='append')
   write(1,*) time, i
   close(1)
   exit
end if
end do


!Write final point data
open(2,file=trim(f_inp)//'.1')
write(2,'(2i5)') natomf, nst
write(2,'(i10)') branch
write(2,'(2f15.6)') time, TimeStep

! Write geometry and momenta
write(2,'(2i3)')nc,mult
do i=1,natom
  if (i.ne.nfr) write(2,'(a1,3f22.16)') atom(i),r1(i*3-2:i*3)
end do
write(2,*)
do i=1,natom
  if (i.ne.nfr) write(2,'(3e25.16)') p1(i*3-2:i*3)
end do

! Write amplitudes
write(2,*)
do i=1,nst
  write(2,*) A1(i)
end do
write(2,*)

!Write electronic structure data
do i=1,nst
  write(2,'(e25.16,i8)') Es1(i), i
end do

write(2,*)
do i=1,nst
  n=0
  do j=1,ndim
    if (j==nfr*3.or.j==nfr*3-1.or.j==nfr*3-2) cycle
    n=n+1
    write(2,'(e25.16,2i8)') Fs1(j,i), i,n
  end do
end do

write(2,*)
do i=1,nst-1
  do k= i+1,nst
    n=0
    do j=1,ndim
    if (j==nfr*3.or.j==nfr*3-1.or.j==nfr*3-2) cycle
    n=n+1
      write(2,'(e25.16,3i8)') Cs1(j,i,k), i,k,n
    end do
  end do
end do

close(2)

!Calculate and write the total energy
E_tot=sum(Es1*abs(A1)**2)+0.5*sum((p1**2)/M)
dr=sqrt(sum((r1-r0)**2))
open(unit=2,file=trim(f_inp)//'.te')
write(2,*) E_tot
write(2,*) dot_product(F1,(r1-r0))/dr,-(sum(Es1*abs(A1)**2)-sum(Es0*abs(A0)**2))/dr
do i=1,nst
   write(2,*) i, 0.5*dot_product((Fs1(:,i)+Fs0(:,i)),(r1-r0))/dr, -(Es1(i)-Es0(i))/dr
enddo
close(2)


deallocate (A0, A1, Ab)
deallocate (Es0,Es1,Esb,Fs0,Fs1,Fsb)
deallocate (Cs0,Cs1,Csb)
deallocate (HE_0,HE_1,HE_b)
deallocate (F0,F1,Fb,V0,V1)
deallocate (M,ElPhase)
deallocate(r0,p0,r1,p1,atom)




contains 

      Function CompForceEhr(A,F,E,C,nst,ndim) result(ForceVector)
      implicit none
      real (kind=8) :: ForceVector(ndim)
      real (kind=8) :: F(ndim,nst)
      real(kind=8)::C(ndim,nst,nst)
      real (kind=8) :: f1(ndim),f2(ndim)
      real (kind=8) :: E(nst),ae
      complex(8)    :: a(nst)
      integer  :: ndim,nst,i,j
      F1=0.d0
      F2=0.d0
      do i=1,nst
        F1=F1+F(:,i)*cdabs(A(i))**2
      end do
      do i=1,nst
        do j=i+1,nst
          ae=2.d0*dreal(dconjg(A(i))*A(j))*(E(i)-E(j))
          F2=F2+ae*C(:,i,j)
        end do
      end do
      ForceVector=F1+F2
      return
      end function CompForceEhr








function magnus2( H0, H1, dt ) result( magH )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! -
complex*16,intent(in) :: H0(:,:), H1(:,:)
real*8 :: dt
complex*16,dimension(size(H0,1),size(H0,2)) :: magH

complex*16,dimension(size(H0,1),size(H0,2)) :: a0, a1, a2, W1, Htr    !derivatives at t/2
complex*16 :: Hav
integer :: ndim, i

ndim = size(H0,1)

! Calculate the average
Hav = 0.d0
do i = 1, ndim
   Hav = Hav + H0(i,i) + H1(i,i)
enddo
Hav = Hav / dble( 2*ndim )

! THe trace matrix
Htr = 0.d0
do i = 1, ndim
   Htr(i,i) = Hav
enddo

a0 =        (H1+H0)/2.d0 - Htr

W1 = dt * a0

!magH = matrix_exponential( W1 )*exp(Hav*dt)
 magH = exp_pade( W1)*exp(Hav*dt)

return
end function magnus2

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! -
function exp_pade( A, t_in ) result ( expA )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! -

complex*16,intent(in) :: A(:,:)
real*8,    optional   :: t_in
complex*16            :: expA(size(A,1),size(A,2))

integer :: n

integer :: ideg, m , ipiv(size(A,1)), &
           iexph, & ! locates the start of expA inthe work array wsp
           lwsp, ldh, ns, iflag
complex*16 :: wsp( 4*size(A,1)**2 + 10 +1 )
real*8 ::  t

t = 1.d0
if( present(t_in) ) t = t_in

ideg = 6
m    = size(A,1)
ldh  = m
lwsp = size(wsp)

if( maxval(abs(A)) < 1.d-7 )then
   expA = 0.
   do n = 1, m
      expA(n,n) = (1.d0,0.d0)
   enddo
else
   call ZGPADM( ideg, m, t, A, ldh, wsp, lwsp, ipiv, iexph, ns, iflag )
   expA = reshape( wsp(iexph:iexph+m**2-1), (/m,m/) )
endif


return
end function exp_pade




end



