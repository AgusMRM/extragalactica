program cinco
        implicit none
        integer i
        real ra,dec,z,bestobjid,pmu,pmg,pmr,pmi,pmz,mmu,mmg,mmr,mmi,mmz,eu
        real eg,er,ei,ez,pR50u,pR50g,pR50r,pR50i,pR50z,pR90u,pR90g,pR90r,pR90i,pR90z,fracDevr,veldisp  
        real d_luminosity, d_angular, Mpu, Mpg, Mpr, Mpi, Mpz, Mmu_abs, Mmg_abs, MMr_abs, Mmz_abs, Mmi_abs 
        real ur, gr
! leo la tabla de datos
        open(7,file='tabla.dat',status='old',action='read')
        open(8,file='distancias.dat',status='old',action='read')
        open(9,file='output.dat',status='unknown',action='write')

        do i=1,50000
              read(7,*) ra,dec,z,bestobjid,pmu,pmg,pmr,pmi,pmz,mmu,mmg,mmr,mmi,mmz,eu,eg,er,ei,ez &
                        ,pR50u,pR50g,pR50r,pR50i,pR50z,pR90u,pR90g,pR90r,pR90i,pR90z,fracDevr,veldisp  
              read(8,*) d_luminosity, d_angular
              if (pmr > 14.5 .and. pmr < 17.77) then
 !------------------------------CALCULO MAGNITUDES ABSOLUTAS---------------------------             
              Mpu=(pmu-eu) - 25. -5.*log10(d_luminosity)
              Mpg=(pmg-eg) - 25. -5.*log10(d_luminosity)
              Mpr=(pmr-er) - 25. -5.*log10(d_luminosity)
              Mpi=(pmi-ei) - 25. -5.*log10(d_luminosity)
              Mpz=(pmz-ez) - 25. -5.*log10(d_luminosity)
              
              Mmu_abs=(pmu-eu) - 25. -5.*log10(d_luminosity)
              Mmg_abs=(mmu-eg) - 25. -5.*log10(d_luminosity)
              Mmr_abs=(mmu-er) - 25. -5.*log10(d_luminosity)
              Mmi_abs=(mmu-ei) - 25. -5.*log10(d_luminosity)
              Mmz_abs=(mmu-ez) - 25. -5.*log10(d_luminosity)
!--------------------------------------------------------------------------------------
!-------------------------------CALCULO COLORES----------------------------------------
              ur=(mmu-eu)-(mmr-er)
              gr=(mmg-mmr)-(mmr-er)
!--------------------------------------------------------------------------------------
              write(9,*) z,Mpu, Mpg, Mpr, Mpi, Mpz, Mmu_abs, Mmg_abs, MMr_abs, Mmz_abs, Mmi_abs, ur, gr
              endif
        enddo

         close(7)
         close(8)
         close(9) 


endprogram
