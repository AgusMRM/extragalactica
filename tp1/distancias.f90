program calculadora_distancias
        implicit none
        integer j,i
        real z, d_luminosity, d_angular, zi(6000), add(6000), ld(6000), cd, cmsa
        
        open(10,file='tabla_redshift.dat',status='old',action='read')
        open(11,file='cosmology.dat',status='old',action='read')
        open(12,file='distancias.dat',status='unknown')

        do i=1,6000
                read(11,*) zi(i),cd,add(i),ld(i),cmsa
        enddo
         close(11)

         do j=1,50000
                read(10,*) z
                 do i=2,6000
                   if ( zi(i) > z .and. zi(i-1) <= z) then ! .and. z < zi(i) ) then 
                          d_luminosity=ld(i-1) + (ld(i)-ld(i-1))*(z-zi(i-1))/(zi(i)-zi(i-1))
                          d_angular=add(i-1) + (add(i)-add(i-1))*(z-zi(i-1))/(zi(i)-zi(i-1))
                          go to 4                  
                   endif
                 enddo
         4 write(12,*) d_luminosity, d_angular, z
         enddo
          close(10)
          close(12)
endprogram


