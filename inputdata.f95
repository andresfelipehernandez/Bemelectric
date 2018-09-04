subroutine inputdata

!En esta subrutina se realizan dos procedimientos:
!1. Se lee el archivo generado por el software GID al enmallar la frontera de una región y se extrae de él la información concerniente a: coordenadas de nodos y nodos que conforman cada elemento.
!2. Se imprime la información leída en un nuevo archivo, que en este caso se denomina "Nodos y elementos"

use vars !Modulo definido en el programa


character*11 texto     !Variable definida con 11 caracteres y que se encaragara de recorrer el archivo de entrada por sus primeras líneas en las que aparacen letras y no números

open(10,file='Frontera.txt',status='unknown')! Abre el archivo donde se enciuentran las coordenadas de los nodos y los elementos 

texto='1234567890&'
a1=1

write(*,*)'Ingrese numero de nodos' ! Para hacer el recorrido en el archivo de texto.
read(*,*)nn
ne=nn/2!por ser cuadráticos
    
write(*,*)'El numero de elementos en la frontera es',ne

allocate(coord(nn,2))     !El número de nodos generados al enmallar la frontera (y por ende el número de elementos) determina el número de filas en los arreglos coord(:,:) y conect(:,:) 
allocate(conect(ne,3))    !esto obliga a reservar o separar espacio en memoria para cuando se definan estas variables. En el caso de elementos cuadráticos tres nodos son ingresados por cada elemento.

!Recorrido por las primeras líneas del archivo

do while (a1==1)
  read(10,*)texto
  if (texto=='Coordinates')a1=0
enddo


do i=1,nn
  read(10,*)dummy,coord(i,1),coord(i,2),dummy !Lee las coordenadas de cada nodo, donde esta Dummy lo ignora
enddo

open(11,file='Nodos y elementos.txt',status='unknown') ! Abre el archivo Nodos y elementos

do i=1,nn 
  write(11,*)i,coord(i,1),coord(i,2) !Imprime las coordenadas en el archivo nodos y elementos 
enddo

!Se repite el trabajo anterior con la información de los elementos.

a1=1

do while(a1==1)
  read(10,*)texto
  !write(*,*)texto
  if(texto=='Elements   ')a1=0
enddo    

do i=1,ne
  read(10,*)dummy,conect(i,1),conect(i,3),conect(i,2)
enddo
  
do i=1,ne
  write(11,*)i,conect(i,1),conect(i,2),conect(i,3)
enddo

close(10)

close(11)

!deallocate(coord,conect)

!Condiciones de frontera: 
!Temperatura igual a 1 en el lado izquierdo (incluyendo las esquinas) y temperatura igual a cero en el lado derecho (incluyendo las esquinas)
!Calor igual a cero en los dos lados restantes (sin incluir las esquinas)


open(40,file='malla interior.txt',status='unknown')

texto='1234567890&'
a1=1
!n=0
  
!Ingreso de información: Número de nodos(nn) y Tipo de elementos:Lineales o cuadráticos. Cálculo del número de elementos(ne)

write(*,*)'Ingrese número de nodos en el interior'
read(*,*)ni
   
allocate(cordi(ni,2))


do while (a1==1)
  read(40,*)texto
  if (texto=='Coordinates')a1=0
enddo

do i=1,ni
  read(40,*)dummy,cordi(i,1),cordi(i,2),dummy
enddo

open(41,file='Nodos al interior.txt',status='unknown')
  
!El siguiente ciclo imprime las coordenadas de los nodos
!en el archivo Nodos al interior

do i=1,ni 
  write(41,*)i,cordi(i,1),cordi(i,2)
enddo

close(40)
close(41)

return
end
