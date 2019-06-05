PROGRAM prog
        ImPlIcIt NoNe
        INTEGER :: i, n, nbin, bin
        
        REAL, ALLOCATABLE, DIMENSION(:) :: distance, galaxy, ldist
        REAL, ALLOCATABLE, DIMENSION(:) :: dE, dS0, dS, dx, d
        REAL :: dmin, dmax, abin, suma
        
        dmin=0.
        dmax=-99999
        n=5725
        nbin=15.

        ALLOCATE(distance(n),galaxy(n),ldist(n))
        ALLOCATE(dE(nbin),dS0(nbin), dS(nbin), dx(nbin), d(nbin))
        
        open(10,file='distcentros.dat',status='old')
        open(11,file='datosgraf2.at',status='unknown')
        do i=1,n
                read(10,*) distance(i), galaxy(i)
                ldist(i) = distance(i) 
                if (ldist(i)>=dmax) dmax=ldist(i)
               ! if (lden(i)<=dmin) dmin=lden(i)
        enddo
        print*, dmin, dmax
!-------------------------------------------------------------
        abin=(dmax-dmin)/float(nbin)
        dE=0
        dS0=0
        dS=0
        d=0
        do i=1,n
                bin = int((ldist(i)-dmin)/abin) + 1
                d(bin) = d(bin) + 1
                if (galaxy(i)==1) dE(bin)=dE(bin)+1
                if (galaxy(i)==2) dS0(bin)=dS0(bin)+1
                if (galaxy(i)==3 .or. galaxy(i) == 4) dS(bin)=dS(bin)+1
        enddo
        do i=1,nbin
                suma = dE(i) + dS0(i) + dS(i) 
                write(11,*) dE(i)/suma, dS0(i)/suma, dS(i)/suma, dmin+(i-.5)*abin
        enddo
        close(10)
        close(11)

        deallocate(distance,galaxy,ldist,de,ds0,ds,dx,d)  
        

ENDPROGRAM
