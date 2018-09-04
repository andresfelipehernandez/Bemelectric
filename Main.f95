module vars
integer a1,nn,ne,npg,ni
real*8,allocatable :: coord(:,:),bc(:),intq(:,:),intt(:,:),ma(:,:),b(:,:),ste(:),sca(:),cordi(:,:),stei(:),scai(:),T(:)
real*8,allocatable :: stefluxi(:),scafluxi(:),Q(:,:)
integer,allocatable :: conect(:,:),bcc(:)
endmodule

!Programa

program Ingresodedatos

use vars
call inputdata
npg=8
call integralq(npg)
call integralt(npg)
call assembly
call SLNPD
call disassembly 
call integralqi(npg)
call integralti(npg)
call fluxint1(npg)
call fluxint2(npg)
stop
end