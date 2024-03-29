PROGRAM three
        IMPLICIT NONE
        INTEGER, PARAMETER :: n=508330, bines=12
        INTEGER :: i,j,k,bin
        REAL :: a,smass_min, smass_max, smas, me, w 
        REAl, DIMENSION(n) :: smass, met, wp
        REAL, DIMENSION(bines) :: med, stdv
        real :: g(3)
!***************************************************************************
        open(10,file='field.dat',status='old')
        k = 0
        smass_min = 999
        smass_max = -999
        do i=1,n
                read(10,*) a,a,a,a,a,a,a,a,smas,a,a,a,me,w
                if (smas < 8.) cycle 
                if (me < -90.) cycle
                        !print*, smas 
                        k = k + 1
                        smass(k) = smas
                        met(k)   = me
                        wp(k)    = w
                        if (smas >= smass_max) smass_max = smas
                        if (smas <= smass_min) smass_min = smas
                
        enddo
        close(10)
!*************************************************************************
        print*, smass_max, smass_min
        call medias(smass,met,wp,k,bines,smass_max,smass_min,med,stdv)
                
        do i=1,bines
                print*, med(i), stdv(i)
        enddo

ENDPROGRAM three

SUBROUTINE medias(x,y,k,bines,mx,mn,med,stdv)
        implicit none
        integer :: i,j,k,n,bines,bin
        real, dimension(k) :: x,y
        integer, dimension(bines) :: y2
        real,dimension(bines) :: med,y1,stdv
        real :: mx, mn, abin
        
        y1  = 0
        y2  = 0
        med = 0
        abin = (mx-mn)/real(bines)
        do i = 1, k
                bin = int(( x(i) - mn)/abin) + 1
                if (bin == bines+1) bin=bines
                y1(bin) = y1(bin) + y(i)
                y2(bin) = y2(bin) + 1
        enddo
        do i=1,bines
                med(i) = y1(i)/real(y2(i))
        enddo
        do i=1,k
                bin = int(( x(i) - mn)/abin) + 1
                if (bin == bines+1) bin=bines
                stdv(bin) = stdv(bin) + (x(i)-med(bin))**2/real(y2(bin)) 
        enddo
ENDSUBROUTINE

SUBROUTINE mediasW(x,y,w,n,k,bines,mx,mn,med,stdv)
        implicit none
        integer :: i,j,k,n,bines,bin
        real, dimension(n) :: x,y,w
        real,dimension(bines) :: med,y1, y2, stdv
        real :: mx, mn, abin
        
        y1  = 0
        y2  = 0
        med = 0
        abin = (mx-mn)/real(bines)
        do i = 1, k
                bin = int(( x(i) - mn)/abin) + 1
                if (bin == bines+1) bin=bines
                y1(bin) = y1(bin) + y(i)*w(i)
                y2(bin) = y2(bin) + w(i)
        enddo
        do i=1,bines
                med(i) = y1(i)/y2(i)
        enddo
        do i=1,k
                bin = int(( x(i) - mn)/abin) + 1
                if (bin == bines+1) bin=bines
                stdv(bin) = stdv(bin) + w(i)*(x(i)-med(bin))**2/y2(bin)
        enddo
ENDSUBROUTINE

!SUBROUTINE mediada(x,y,n,k,bines)
!!        implicit none
!        integer :: i,n,k,bines
!        real, dimension(n) :: x, y
!        real, dimension(k) :: x0, y0
        
!        do i=1,k
!                x0(i) = x(i)
!                y0(i) = y(i)
!        enddo      
!        call orden(k,x0)
        
        
        
!        abin = (mx-mn)/real(bines)
!        do i = 1, k
!                bin = int(( x(i) - mn)/abin) + 1
!                if (bin == bines+1) bin=bines
!                y1(bin) = y1(bin) + y(i)
!                y2(bin) = y2(bin) + 1
!:w        enddo

         

!ENDSUBROUTINE 

SUBROUTINE ORDEN(NELEM,ARREG)
! -----------------------------------------------------
!ORDENACION POR BURBUJA ("buble sort") de un arreglo
! unidimensional, de menor a mayor.
!
! NELEM = Número de elementos del arreglo
! ARREG = Arreglo unidimensional a ordenar
! -----------------------------------------------------
IMPLICIT NONE
INTEGER :: NELEM
REAL :: ARREG(*)
!-----------------------------------------------------
INTEGER:: I,J
REAL:: AUX
!-----------------------------------------------------
IF (NELEM.LT.2) RETURN
DO I=1,NELEM-1
DO J=1,NELEM-I
IF (ARREG(J).GT.ARREG(J+1)) THEN
AUX = ARREG(J)
ARREG(J) = ARREG(J+1)
ARREG(J+1) = AUX
ENDIF
ENDDO
ENDDO
RETURN
END
SUBROUTINE sort2(arr,slave)
USE nrtype; USE nrutil, ONLY : assert_eq
USE nr, ONLY : indexx
IMPLICIT NONE
REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr,slave
!Sorts an array arr into ascending order using Quicksort, while making the corresponding
!rearrangement of the same-size array slave . The sorting and rearrangement are performed
!by means of an index array.
INTEGER(I4B) :: ndum
INTEGER(I4B), DIMENSION(size(arr)) :: index
ndum=assert_eq(size(arr),size(slave),'sort2')
call indexx(arr,index)
!Make the index array.
arr=arr(index)
!Sort arr.
slave=slave(index)
!Rearrange slave.
END SUBROUTINE sort2
