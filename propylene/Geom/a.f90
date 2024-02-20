real, dimension(9,500) :: x,y,z,px,py,pz
character(len=1), dimension(9) :: a
a(1)='C'
a(2:4)='H'
a(5)='C'
a(6)='H'
a(7)='C'
a(8)='H'
a(9)='H'


open(1,file='X')
do i=1,9
read(1,*)x(i,:)
end do
close(1)
open(1,file='Y')
do i=1,9
read(1,*)y(i,:)
end do
close(1)
open(1,file='Z')
do i=1,9
read(1,*)z(i,:)
end do
close(1)
open(1,file='Px')
do i=1,9
read(1,*)px(i,:)
end do
close(1)
open(1,file='Py')
do i=1,9
read(1,*)py(i,:)
end do
close(1)
open(1,file='Pz')
do i=1,9
read(1,*)pz(i,:)
end do
close(1)

open(1,file='Geom')
do i=1,500
write(1,*) 9
do j=1,9
write(1,'(a1,3f17.10)') a(j), x(j,i), y(j,i), z(j,i)
enddo
write(1,*) 'momenta'
do j=1,9
write(1,*) px(j,i), py(j,i), pz(j,i)
enddo
enddo
close(1) 
end
