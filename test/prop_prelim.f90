program dynamics

implicit none
integer :: nst, inist, natom, ndim, n1, n2, i, j, k,n, nr, im, branch
integer :: nc, mult
real (kind=8), allocatable, dimension (:) :: ElPhase
real (kind=8) CompTime, TimeStep, time, E0, E1, Mau,  val
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

! Read geometry
read(1,*) nc, mult
do i=1,natom
  read(1,*) atom(i),r0(i*3-2:i*3)
end do
read(1,*)
do i=1,natom
  read(1,*) p0(i*3-2:i*3)
end do

!Allocation
allocate (A0(nst), A1(nst), Ab(nst))
allocate (Es0(nst),Fs0(ndim,nst))
allocate (Cs0(ndim,nst,nst))
allocate (HE_0(nst,nst))
allocate (F0(ndim),V0(ndim),V1(ndim))
allocate (M(ndim))

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
  else if (atom(i)=='O') then
    M(i*3-2:i*3)=16*Mau
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
!do i=1,nst
!  do j=1,ndim
    read(1,*) val, i,j
    Fs0(j,i)=val
!  end do
end do

! Read couplings
read(1,*)
nr=ndim*nst*(nst-1)/2
do n=1,nr
!do i=1,nst-1
!  do k=i+1,nst
!    do j=1,ndim
      read(1,*) val,i,k,j
      Cs0(j,i,k)=val
      Cs0(j,k,i)=-Cs0(j,i,k)
!    end do
!  end do
end do

close(1)


!Electronic structure
!call ElStr(R0,Es0,Fs0,Cs0,nst,ndim)
!Electronic phase
!ElPhase=1.

! Velocities and Ehrenfest force
V0=P0/M
F0=CompForceEhr(A0,Fs0,Es0,Cs0,nst,ndim)
E0=sum(Es0*(abs(A0)**2))+sum(P0*V0)/2.

!Electronic Hamiltonian
do n1=1,nst
  HE_0(n1,n1)=Es0(n1)+77.67785291
  do n2=n1+1,nst
    HE_0(n1,n2)=-ii*sum(V0*Cs0(:,n1,n2))
    HE_0(n2,n1)=-HE_0(n1,n2)
  end do
end do


!Begin propagation

  Ab=matmul( magnus2( -ii*HE_0, -ii*HE_0, TimeStep/20), A0 )
  F0 = CompForceEhr(A0,Fs0,Es0,Cs0,nst,ndim)/10.
!print*,abs(ab)**2
  do im=1,9
    A1 = matmul( magnus2( -ii*HE_0, -ii*HE_0, TimeStep/10), Ab )
    Ab=A1
!print*,abs(ab)**2
    F0 = F0+CompForceEhr(Ab,Fs0,Es0,Cs0,nst,ndim)/10.
  end do
!  A1 = matmul( magnus2( -ii*HE_0, -ii*HE_0, TimeStep/20), Ab )


R1 = R0 + TimeStep*V0 + TimeStep**2/2.d0 * F0/M
P1 =               P0 + TimeStep * F0
V1=P1/M

write(6,*) "Mass", M
write(6,*)
write(6,*) "Scaled forces",  F0/M
write(6,*) 
write(6,*) "Increment scaled", TimeStep**2/2.d0
write(6,*)
write(6,*) "Changed Forces: ", TimeStep**2/2.0d0 * F0/M
write(6,*)

time=time+TimeStep

!Write preliminary final point data
open(2,file=trim(f_inp)//'.p')
write(2,'(2i5)') natom, nst
write(2,'(i10)') branch
write(2,'(2f15.6)') time, TimeStep

! Write geometry and momenta
write(2,'(2i3)') nc, mult
do i=1,natom
  write(2,'(a1,3f22.16)') atom(i),r1(i*3-2:i*3)
end do
write(2,*)
do i=1,natom
  write(2,'(3e25.16)') p1(i*3-2:i*3)
end do

! Write amplitudes
write(2,*)
do i=1,nst
  write(2,*) A1(i)
end do
write(2,*)
close(2)

deallocate (A0, A1, Es0,Fs0,Cs0,HE_0,F0,V0,V1, M)
deallocate(r0,p0,r1,p1, atom)

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
! write(6,*) "H0: ", H0
! write(6,*) "H1: ", H1
! write(6,*) "Hav: ", Hav
! write(6,*) "Htr: ", Htr
! write(6,*) "********************"
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



