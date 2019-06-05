program calculadora_de_distancias
        implicit none
        integer i
        real z, d_luminosity, d_angular, zi(6006), add(6006), ld(6006), cd, cmsa
        
! ingreso redshift por terminal
        write(*,*) 'ingrese z'
        read(*,*) z
        
        open(11,file='cosmology.dat',status='unknown',action='read')
!armo los vectores con todos los datos de la cosmology
        do i=1,6000
                read(11,*) zi(i),cd,add(i),ld(i),cmsa
        enddo
         close(11)
!aca busco entre que valores de la tabla esta mi z, y realizo una interpolacion lineal.        
        do i=2,6000
                if ( zi(i) > z .and. zi(i-1) <= z) then ! .and. z < zi(i) ) then 
                       d_luminosity=ld(i-1) + (ld(i)-ld(i-1))*(z-zi(i-1))/(zi(i)-zi(i-1))
                       d_angular=add(i-1) + (add(i)-add(i-1))*(z-zi(i-1))/(zi(i)-zi(i-1))
                       go to 4                  
                endif
        enddo
4        write(*,*) 'distancia luminosidad',d_luminosity,'Mpc'
         write(*,*) 'distancia diametro angular',d_angular,'Mpc'
endprogram
