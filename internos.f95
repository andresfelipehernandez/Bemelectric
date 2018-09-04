integer ni
real lx,ly

nx=10
ny=10
ni=nx*ny
lx=linmax(1)-limin(1)
ly=linmax(2)-limin(2)
lx=lx/(nx+1)
ly=ly/(ny+1)


c1=0
do i=1,nx
  do j=1,ny
    c1=c1+1
    cordi(c1,1)=limin(1)+lx*i      !llenar las tablas
    cordi(c1,2)=limin(2)+ly*j
  enddo
enddo

open(10,file='cordi.txt',status='unknown')
do i,ni
  write(10,*)cordi(1,1) cordi(1,2)
  