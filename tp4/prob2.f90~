!gfortran -fbacktrace -fbounds-check program.f90
program main
       use numerical
       implicit none
       integer,parameter::k2=selected_real_kind(10,20)
       integer,parameter::k1=selected_real_kind(6,10)
       integer :: i, j,k, cumulos, miembros
       integer :: nlines_g, nlines_c, c  
       real(kind=k2), allocatable, dimension(:) :: ar, dec
       real(kind=k1) :: distancia
       real(kind=k2) :: a_cum, d_cum
       integer, allocatable :: inde(:)
       integer, allocatable :: tipo(:)
       real(kind=k2) :: a, d, a2, d2,z , dist, theta, zc
       real :: c_luz, H,  pi  , radio    
       
       c_luz = 3.e5
       H = 70.
       pi = acos(-1.)
       
       open(10,file='centros.dat',status='unknown')
       open(11,file='galaxias.dat',status='unknown')
       open(15,file='distcentros.dat',status='unknown') 
       
       call count_lines(10,nlines_c)
       call count_lines(11,nlines_g)
      !--------------- leo lineas que no uso
       read(10,*)          
       read(11,*) 
      !---------------
       do i=1,nlines_c - 1 !nlines_c-1
        read(10,*) cumulos, miembros, zc, a_cum, d_cum
                allocate(ar(miembros),dec(miembros),tipo(miembros))
                
                do j=1,miembros
                        read(11,*)  c,c,ar(j),dec(j),tipo(j)     
                enddo  
                do j=1,miembros
                        a=ar(j)
                        d=dec(j)
                        
                        call distanciaytheta_rad(a_cum, d_cum, zc, a, d, zc, dist, theta)

                        distancia = tan(theta*(pi/180.))*c_luz*zc/H
                        write(15,*) distancia, tipo(j)
                enddo 
                deallocate(ar,dec,tipo)
       enddo

       close(10) 
       close(11)
       close(15)
endprogram 

!.................................................................
        subroutine count_lines(nunit,nlines)
        implicit none
        integer, intent(in)            :: nunit
        integer, intent(out)           :: nlines
        integer                        :: check
        nlines = 0
        do
        read(nunit,*,iostat=check)
        if(check /= 0)exit
         nlines = nlines + 1
        end do
        rewind(nunit)
        end subroutine count_lines
!.................................................................

      subroutine distanciaytheta_rad(alf1,del1,z1,alf2,del2,z2,d,theta)

!*************** VARIABLES
      integer,parameter::k2=selected_real_kind(10,20)
      real(kind=k2):: ALF1,DEL1,ALF2,DEL2,Z1,Z2,z1t,z2t
      real(kind=k2):: alft(2),DELt(2)
      real(kind=k2):: D
      real(kind=k2):: D1,P1,THETA


      z1t=z1
      z2t=z2
!*************** CALCULO DE LA DISTANCIA MUTUA
      Itt=1
      Jtt=2
      ALFt(Itt)=ALF1
      DELt(Itt)=DEL1
      ALFt(Jtt)=ALF2
      DELt(Jtt)=DEL2

      pi=3.141592654
      fc=pi/180.
      ALFt(1)=ALFt(1)*fc
      DELt(1)=DELt(1)*fc
      ALFt(2)=ALFt(2)*fc
      DELt(2)=DELt(2)*fc


        P1=ABS(alft(Itt)-alft(Jtt))
!        IF (P1<=180.) GOTO 410
        IF (P1<=pi) GOTO 410
!        P1=360.-P1
        P1=2.*pi-P1
 410    if (abs(SIN(Delt(Jtt))*SIN(Delt(Itt))+COS(Delt(Jtt))*COS(Delt(Itt))*COS(P1))>1.) then
        write (*,*) 'argu',SIN(Delt(Jtt))*SIN(Delt(Itt))+COS(Delt(Jtt))*COS(Delt(Itt))*COS(P1)
        d1=-1.
        goto 224
        end if
        D1=ACOS(SIN(Delt(Jtt))*SIN(Delt(Itt))+COS(Delt(Jtt))*COS(Delt(Itt))*COS(P1))   
  224   IF (D1<0.) write(*,*) 'd1<0 en j,i=',jtt,itt

        THETA=D1


!************************ DISTANCIAS PROYECTADAS
!      IF(IZF==0) THEN
!      WRITE(*,*) '!!! EN EL DISTANCIATHETA LOS Z SON IGUALES'
!      IZF=1
!      END IF
!      Z1t=10000.
!      Z2t=10000.

!      z1t=0.
!      z2t=0.

      IF (Z1t==0.AND.Z2t==0) GOTO 987

!****** H=100
!      Z1t=Z1t/100.
!      Z2t=Z2t/100.
!****** H=70
      Z1t=Z1t/70.
      Z2t=Z2t/70.
      D1=(Z1t**2.+Z2t**2.-2.*Z1t*Z2t*COS(D1))**.5
 987  D=D1
!      WRITE(*,*) '-------',D
      THETA=THETA/fc
      RETURN
      END
