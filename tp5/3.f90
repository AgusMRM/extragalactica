PROGRAM three
        IMPLICIT NONE
        INTEGER, PARAMETER :: n=508330, bines=12
        INTEGER :: i,j,k,bin1, bin2
        REAL, DIMENSION (n) :: smass, met, wp
        REAL :: a, smas, me, smass_min, smass_max, met_min, met_max,abin, abin2, w
        INTEGER, DIMENSION(bines) :: num
        REAL, DIMENSION (bines) :: val_smass,val_met, media_smass, media_met, stdv_smass, stdv_met, num2
     !************************************************************   
        open(10,file='field.dat',status='old')
        k = 0
        smass_min = 999
        smass_max = -999
        met_min = 999
        met_max = -999

        do i=1,n
                read(10,*) a,a,a,a,a,a,a,a,smas,a,a,a,me,w
                if (smas < 8.) cycle 
                if (me < -90.) cycle 
                        k = k + 1
                        smass(k) = smas
                        met(k)   = me
                        wp(k)    = w
                        !if (smas > smass_max) smass_max = smas
                        !if (smas < smass_min) smass_min = smas
                        if (me >= met_max) met_max = me
                        if (me <= met_min) met_min = me
                
        enddo
        close(10)
        print*, k, met_min, met_max
        smass_max=11.5
        smass_min=8.
     !************************************************************
        abin  = (smass_max - smass_min)/real(bines)
        abin2 = (met_max - met_min)/real(bines)
        num   = 0
        num2  = 0
        val_smass  = 0
        val_met = 0
        do i=1,k
                if (smass(i) <= smass_min .or. smass(i) >= smass_max) cycle
               bin1       = int(((smass(i) - smass_min))/abin) + 1
                !bin2       = int(((met(i) - met_min))/abin2) + 1
                !val_smass(bin1)  = val_smass(bin1) + smass(i)
                val_met(bin1)  = val_met(bin1) + met(i)
                num(bin1)  = num(bin1) + 1.
                !num2(bin1)  = num2(bin1) + wp(i)
        enddo
        
     !************************************************************    
        open(11,file='medias2.dat',status='unknown')
        do i=1,bines
                !media_smass(i) = val_smass(i)/num(i)
                media_met(i) = val_met(i)/(num2(i))
        enddo
        do i=1,k
                if (smass(i) <= smass_min .or. smass(i) >= smass_max) cycle
                bin1             = int(((smass(i) - smass_min))/abin) + 1
                !stdv_met(bin1) = stdv_met(bin1) + wp(i)*(media_met(bin1)-met(i))**2/(num2(bin1))                  
                stdv_met(bin1) = stdv_met(bin1) + (media_met(bin1)-met(i))**2/(num2(bin1))                  

        enddo
        do i=1,bines
                
                write(11,*) i, media_met(i), (i-0.5)*abin+smass_min, sqrt(stdv_met(i))
        enddo
        close(11)
ENDPROGRAM three