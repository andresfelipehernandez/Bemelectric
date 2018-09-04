subroutine integralt(npg)

!En esta subrutina se calcula la segunda integral sobre la frontera que aparecen en
!la ecuación integral básica, la correspondiente a la solución fundamental y el calor.

use vars

real*8 x(1,2),y(3,2),phi(3),dphi(3),psi(npg),ypsi(2),dypsi(2),xb(2),r,jac,no(2),ts,pi,wg(npg)

pi=3.1415926535897932384626433832795

psi(1)=-0.960289856497537
psi(2)=-0.796666477413627
psi(3)=-0.525532409916329
psi(4)=-0.183434642495650
psi(5)=0.183434642495650
psi(6)=0.525532409916329
psi(7)=0.796666477413627
psi(8)=0.960289856497536

wg(1)=0.101228536290376
wg(2)=0.222381034453375
wg(3)=0.313706645877887
wg(4)=0.362683783378362
wg(5)=0.362683783378362
wg(6)=0.313706645877887
wg(7)=0.222381034453374
wg(8)=0.101228536290376    

allocate(intt(nn,nn))       !La matriz intt almacena los valores de la integral para cada x (x es un nodo tomado sobre la frontera)

intt=0

do j=1,nn
  x(1,1:2)=coord(j,1:2)                                           
  !write(*,*)'gamma numero',j,x(1,1:2)
  do k=1,ne                                           !Este ciclo resuelve la primera matriz para un valor de x fijo; para ello se suman las integrales resultantes en cada elemento, k, sobre el borde
    !write(*,*)'elemento numero',k                        
    do ii=1,3                                         !Este ciclo lee las coordenadas de los tres nodos que conforman el elemento k: y(1,1:2),y(2,1:2) y y(3,1:2)
      y(ii,1:2)=coord(conect(k,ii),1:2)
      !write(*,*)k,ii,y(ii,1:2)
    enddo
    
    do i=1,npg
      phi(1)=-0.5*psi(i)*(1-psi(i))
      phi(2)=(1+psi(i))*(1-psi(i))
      phi(3)=0.5*psi(i)*(1+psi(i))
      dphi(1)=psi(i)-0.5
      dphi(2)=-2*psi(i)
      dphi(3)=psi(i)+0.5
      !write(*,*)phi(1),phi(2),phi(3),dphi(1),dphi(2),dphi(3)
      
      ypsi=0.0 ; dypsi=0.0
      
      !En este ciclo se genera la función paramétrica que describe cada elemento en términos de los puntos de Gauss.
      
      do ii=1,3                                             
        ypsi(1:2)=ypsi(1:2)+y(ii,1:2)*phi(ii)
        dypsi(1:2)=dypsi(1:2)+y(ii,1:2)*dphi(ii)
      enddo
      
      !Cálculo de los factores necesarios para hallar la integral.
        
      xb(1:2)=x(1,1:2)-ypsi(1:2)
      r=dsqrt(xb(1)**2+xb(2)**2)
      jac=dsqrt(dypsi(1)**2+dypsi(2)**2)
     ! no(1)=dypsi(2)/jac
     ! no(2)=-dypsi(1)/jac
      ts=-0.5*(dlog(r))/(pi)
      !write(*,*)ts,jac
      
      do ii=1,3
        intt(j,conect(k,ii))=intt(j,conect(k,ii))+ts*jac*wg(i)*phi(ii)
      enddo
      !write(*,*)intt(j,conect(k,ii))
    enddo
  enddo
  
open(14,file='matrizt.txt',status='unknown')
write(14,*)'Gamma numero',j,intt(j,:)
  
enddo

close(14)

return
end

