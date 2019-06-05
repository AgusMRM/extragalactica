program l
        implicit none
        real, allocatable :: xm(:)
        integer i, bin, nbin
        real mmin, mmax, mag, abin, peso,mp


        mmin=-23.
        mmax=-16.
        
        nbin=10
        abin=(mmax-mmin)/float(nbin)
        allocate(xm(nbin))

        xm=0
        
        open(10,file='datos.dat',status='unknown')
        open(11,file='hist.dat',status='unknown')
        
        do i=1, 539439
                read(10,*) mag, peso
                
                bin=int((mag-mmin)/abin) +1
                
                xm(bin)=xm(bin)+ peso
        enddo

        do i=1,nbin
                write(11,*) xm(i)/abin, i*abin + mmin
        enddo
        close(10)
        close(11) 

endprogram
