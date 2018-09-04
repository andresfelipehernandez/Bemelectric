subroutine assembly

use vars

allocate (ma(nn,nn))
allocate (b(nn,1))
allocate (bcc(nn))
allocate (bc(nn))

open(20,file='Condfront.txt',status='unknown')

!Ciclo en el que se asignan las siguientes condiciones de frontera al rectángulo:
!T=0 sobre el lado derecho. Sobre este lado todos los nodos tienen primera coordenda y1=8.35526
!T=1 sobre el lado izquierdo. Sobre este lado todos los nodos tienen primera coordenda y1=-9.73684
!Q=0 sobre el lado inferior. Sobre este lado todos los nodos tienen segunda coordenda y2=-1.94079
!Q=0 sobre el lado superior. Sobre este lado todos los nodos tienen segunda coordenda y2=6.77632
!Sobre las esquinas se asigna condición de temperatura.

do j=1,nn
  if((coord(j,1)>4).and.(coord(j,1)<6)) then !con 1 es x
    bcc(j)=1
    bc(j)=8
    write(20,*)'Nodo numero',j,coord(j,1),coord(j,2),bcc(j),bc(j)
elseif ((coord(j,2)>4).and.(coord(j,2)<6)) then 
    bcc(j)=1
    bc(j)=12
    write(20,*)'Nodo numero',j,coord(j,1),coord(j,2),bcc(j),bc(j)
elseif  then 
    bcc(j)=2
    bc(j)=0
    write(20,*)'Nodo numero',j,coord(j,1),coord(j,2),bcc(j),bc(j)
endif
!do j=1,nn
 ! if(coord(j,1)==minval(coord(:,1)))then
  !  bcc(j)=1
   ! bc(j)=1
    !write(20,*)'Nodo numero',j,coord(j,1),coord(j,2),bcc(j),bc(j),'Izquierda'
  !else if (coord(j,1)==maxval(coord(:,1))) then
   ! bcc(j)=1
    !bc(j)=0
    !write(20,*)'Nodo numero',j,coord(j,1),coord(j,2),bcc(j),bc(j),'Derecha'
  !else if ((coord(j,2)==minval(coord(:,2))).and.(coord(j,1)>minval(coord(:,1))).and.(coord(j,1)<maxval(coord(:,1)))) then
   ! bcc(j)=2
    !bc(j)=0
    !write(20,*)'Nodo numero',j,coord(j,1),coord(j,2),bcc(j),bc(j),'Abajo'
  !else if ((coord(j,2)==maxval(coord(:,2))).and.(coord(j,1)>minval(coord(:,1))).and.(coord(j,1)<maxval(coord(:,1)))) then
   ! bcc(j)=2
    !bc(j)=0
    !write(20,*)'Nodo numero',j,coord(j,1),coord(j,2),bcc(j),bc(j),'Arriba'
  !endif
  
enddo

close(20)
open(15,file='matrizA.txt',status='unknown')

!Ciclo de creación de la matriz A,ma, y del vector B,b, para el sistema AX=B

b=0.0d0
ma=0.0d0

do i=1,nn
  if (bcc(i)==1) then
    ma(:,i)=intt(:,i)
    b(1:nn,1)=b(1:nn,1)-intq(1:nn,i)*bc(i)
    write(15,*)'Columna',i,ma(:,i)
  elseif (bcc(i)==2) then
    ma(:,i)=intq(:,i)
    b(1:nn,1)=b(1:nn,1)-intt(1:nn,i)*bc(i)
    write(15,*)'Columna',i,ma(:,i)
    endif
enddo
write(15,*)'vector',b(1:nn,1)

close (15)

return
end

    
subroutine disassembly

use vars

allocate (ste(nn),sca(nn))

do i=1,nn
  if (bcc(i)==1) then !Temperatura
    ste(i)=bc(i)
    sca(i)=b(i,1)
  else if (bcc(i)==2) then ! Calor
    ste(i)= b(i,1)
    sca(i)=bc(i)
  endif
enddo

open(16,file='datos.txt',status='unknown')
do i=1,nn
write(16,*)i,ste(i),sca(i)
enddo

close (16)

return
end

