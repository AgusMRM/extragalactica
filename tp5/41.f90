PROGRAM three
        IMPLICIT NONE
        INTEGER, PARAMETER :: n=508330, bines=10
        INTEGER :: i,j,k,bin,l,u,h
        REAL :: a,smass_min, smass_max, smas, ssfr, w, median, q_25, q_75 ,abin,c
        REAl, DIMENSION(n) :: smass1, met1, wp1, smass2, met2, wp2,ssfr1,ssfr2
        REAL,DIMENSION(n) :: x,y,wi,cum 
        REAL, DIMENSION(bines) :: med, stdv,xm,mdn
        real :: m(5),w0(5),cum0(5)
        real, allocatable :: nose1(:), nose2(:), nose3(:)
!***************************************************************************
        open(10,file='field.dat',status='old')
        k = 0
        h = 0
        do i=1,n
                read(10,*) a,a,a,a,a,a,a,c,smas,a,ssfr,a,a,w
                if (smas < 8.) cycle 
                if (ssfr < -20.) cycle
                if (c > 2.86) then
                        k = k + 1
                        smass1(k) = smas
                        ssfr1(k)   = ssfr
                        wp1(k)    = w
                else
                        h = h + 1
                        smass2(h) = smas
                        ssfr2(h)   = ssfr
                        wp2(h)    = w
                endif        
                
        enddo
        print*, h,k
        close(10)
!*************************************************************************
        smass_max = 11.5
        smass_min = 8.
      !  call medias(smass,met,k,bines,smass_max,smass_min,xm,med,stdv)
      !  open(22,file='mediasmain.dat',status='unknown')
      !  do i=1,bines
      !          write(22,*) med(i), stdv(i), xm(i)
      !  enddo
      !  close(22)
!************************************************************************
       ! do i =1,5
       !         print*, 'ingrese valores y pesos'
       !         read(*,*) m(i), w0(i)
       ! enddo
        !call wmediana(5,m,w0,median,q_25,q_75,cum0)
       
        abin=(smass_max-smass_min)/real(bines)
        open(14,file='41E.dat',status='unknown')
        l = 0
        do i=1,bines
                do j=1,k
                        !if (smass(j) < 8.) cycle
                        !if(met(j)z-90)
                        if (smass1(j) >= (i*abin + smass_min) .or. smass1(j) < ((i-1)*abin + smass_min))  cycle
                        l = l + 1
                        x(l) = smass1(j)
                        y(l) = ssfr1(j)
                        wi(l) = 1. !wp(j)
                enddo

                allocate(nose1(l),nose2(l),nose3(l))
                do u=1,l
                        nose1(u) = x(u)
                        nose2(u) = y(u)
                        nose3(u) = wi(u) 
                enddo
                call sort2(l,nose2,nose3)
                call wmediana(l,nose2,nose3,median,q_25,q_75,cum)
                write(14,*) i,(i-.5)*abin + smass_min, median, median-q_25, q_75-median
                l = 0
                
                deallocate(nose1,nose2,nose3)
        enddo
        close(14)
        
        abin=(smass_max-smass_min)/real(bines)
        open(14,file='41L.dat',status='unknown')
        l = 0
        do i=1,bines
                do j=1,h
                        !if (smass(j) < 8.) cycle
                        !if(met(j)z-90)
                        if (smass2(j) >= (i*abin + smass_min) .or. smass2(j) < ((i-1)*abin + smass_min))  cycle
                        l = l + 1
                        x(l) = smass2(j)
                        y(l) = ssfr2(j)
                        wi(l) = 1. !wp(j)
                enddo

                allocate(nose1(l),nose2(l),nose3(l))
                do u=1,l
                        nose1(u) = x(u)
                        nose2(u) = y(u)
                        nose3(u) = wi(u) 
                enddo
                call sort2(l,nose2,nose3)
                call wmediana(l,nose2,nose3,median,q_25,q_75,cum)
                write(14,*) i,(i-.5)*abin + smass_min, median, median-q_25, q_75-median
                l = 0
                
                deallocate(nose1,nose2,nose3)
        enddo
        close(14)

ENDPROGRAM three

SUBROUTINE medias(x,y,k,bines,mx,mn,xm,med,stdv)
        implicit none
        integer :: i,j,k,n,bines,bin
        real, dimension(k) :: x,y
        integer, dimension(bines) :: y2
        real,dimension(bines) :: med,y1,stdv,xm
        real :: mx, mn, abin
        
         
        y1  = 0
        y2  = 0
        med = 0
        abin = (mx-mn)/real(bines)
        stdv = 0
        do i = 1, k
                if (x(i) < mn .or. x(i) > mx) cycle
                bin = int(( x(i) - mn)/abin) + 1
                if (bin == bines+1) bin=bines
                y1(bin) = y1(bin) + y(i)
                y2(bin) = y2(bin) + 1
        enddo
        do i=1,bines
                med(i) = y1(i)/real(y2(i))
                xm(i)  = (i-0.5)*abin + mn
        enddo
        do i=1,k
                if (x(i) < mn .or. x(i) > mx) cycle
                bin = int(( x(i) - mn)/abin) + 1
                if (bin == bines+1) bin=bines
                if (bin < 1 ) bin=1
                stdv(bin) = stdv(bin) + (y(i)-med(bin))**2 
        enddo
        do i=1,bines
                stdv(i) = sqrt(stdv(i)/real(y2(i)))
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

subroutine wmediana(n,x,w,median,q_25,q_75,cum)
        implicit none
        integer :: i,k,n
        real :: median, q_25, q_75, q25, q75, q50
        real,dimension(n) :: cum, w, x
        
        cum = 0
        cum(1) = w(1)
        do i=2,n
                cum(i) = cum(i-1) + w(i)
        enddo
        q50 = cum(n)/2.
        q25 = cum(n)/4.
        q75 = cum(n)*3./4.

        call locate(cum,n,q50,k)
        median = (x(k+1) + x(k))/2.
       
        call locate(cum,n,q25,k)
        q_25 = (x(k+1) + x(k))/2.

        call locate(cum,n,q75,k)
        q_75 = (x(k+1) + x(k))/2.

endsubroutine wmediana

SUBROUTINE locate(xx,n,x,j)
           INTEGER j,n
           REAL x,xx(n)
           INTEGER jl,jm,ju
           jl=0 
           ju=n+1  
           10 if(ju-jl.gt.1)then
           jm=(ju+jl)/2 
           if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
           jl=jm
           else
           ju=jm
           endif
           goto 10 
           endif
           if(x.eq.xx(1))then
           j=1
           else if(x.eq.xx(n))then
           j=n-1
           else
           j=jl
           endif
           return
   END

SUBROUTINE sort2(n,arr,brr)
      INTEGER n,M,NSTACK
      REAL arr(n),brr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,b,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          b=brr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
11        continue
          i=0
2         arr(i+1)=a
          brr(i+1)=b
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        temp=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
          temp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
          temp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
          temp=brr(l+1)
          brr(l+1)=brr(l)
          brr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
        b=brr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        temp=brr(i)
        brr(i)=brr(j)
        brr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        brr(l)=brr(j)
        brr(j)=b
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort2'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
