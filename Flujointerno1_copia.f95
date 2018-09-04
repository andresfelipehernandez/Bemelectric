subroutine fluxint1(npg)

!En esta subrutina se calcula una de las dos integrales sobre la frontera que aparecen en
!la ecuación integral básica, la correspondiente a la temperatura y la derivada normal de la solución fundamental.

use vars

real*8 x(1,2),y(3,2),phi(3),dphi(3),psi(npg),ypsi(2),dypsi(2),xb(2),r,jac,no(2),dqs,pi,wg(npg),tpsi,tp(3)
allocate(stefluxi(ni))
stefluxi=0.0

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

do j=1,ni
  x(1,1:2)=cordi(j,1:2)                                           
  !write(*,*)'gamma numero',j,x(1,1:2)
  do k=1,ne                                           !Este ciclo resuelve la primera matriz para un valor de x fijo; para ello se suman las integrales resultantes en cada elemento, k, sobre el borde
    !write(*,*)'elemento numero',k                        
    do ii=1,3                                         !Este ciclo lee las coordenadas de los tres nodos que conforman el elemento k: y(1,1:2),y(2,1:2) y y(3,1:2)
      y(ii,1:2)=coord(conect(k,ii),1:2)
      tp(ii)=ste(conect(k,ii))
      !write(*,*)k,ii,y(ii,1:2)
    end do
    
    do i=1,npg
      phi(1)=-0.5*psi(i)*(1-psi(i))
      phi(2)=(1+psi(i))*(1-psi(i))
      phi(3)=0.5*psi(i)*(1+psi(i))
      dphi(1)=psi(i)-0.5
      dphi(2)=-2*psi(i)
      dphi(3)=psi(i)+0.5
      !write(*,*)phi(1),phi(2),phi(3),dphi(1),dphi(2),dphi(3)
      
      ypsi=0.0 ; dypsi=0.0 ; tpsi=0.0
            
      !En este ciclo se genera la función paramétrica que describe cada elemento en términos de los puntos de Gauss.
      
      do ii=1,3                                             
        ypsi(1:2)=ypsi(1:2)+y(ii,1:2)*phi(ii)
        tpsi=tpsi+tp(ii)*phi(ii)
        dypsi(1:2)=dypsi(1:2)+y(ii,1:2)*dphi(ii)
      enddo
      
      !Cálculo de los factores necesarios para hallar la integral.

      xb(1:2)=x(1,1:2)-ypsi(1:2)                      
      r=dsqrt(xb(1)**2+xb(2)**2)
      jac=dsqrt(dypsi(1)**2+dypsi(2)**2)
      no(1)=dypsi(2)/jac
      no(2)=-dypsi(1)/jac
      dqs=0.5*(no(1)-(2*xb(1)*(xb(1)*no(1)+xb(2)*no(2)))/(r**2))/(pi*r**2)
      !write(*,*)qs,jac
      stefluxi(j)=stefluxi(j)+dqs*tpsi*jac*wg(i)
      
      !write(*,*)intq(j,conect(k,ii))
    enddo
  enddo
  
!open(12,file='matrizqi.txt',status='unknown')
!write(12,*)'Gamma numero',j,intq(j,:)
 
enddo

!close(12)

return
end