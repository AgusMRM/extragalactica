module cosm
        implicit none
        integer i, j
        integer, parameter :: nlist=6000, ngalaxy=544449
        real, dimension(ngalaxy) :: z,m,r,k,mc,e
        real, dimension(nlist) :: zlist, cd, ad, ld, dllist
        real :: dl, da, Mag,x
endmodule

program calc
        use cosm
        implicit none
        real :: a,d,m1,e1,k1, distance, dlum
        real, parameter :: rmin=14.5, rmax=17.77
        real :: dlmin, dlmax, zmin, zmax, vol, func
        external func, distance, dlum
        open(9,file='mgs.dat',status='old',action='read')
        open(10,file='cosmology.dat',status='old',action='read')
        open(11,file='datosprueba.dat',status='unknown',action='write')


!-----------------------armo vector cosmologico..............
        do i=1,6
        read(10,*)
        enddo
        do i=1,nlist
                read(10,*) zlist(i), cd(i), ad(i), dllist(i)       
        enddo
!---------------------cargo datos y corrijo magnitudes.....................
        do i=1,ngalaxy
                read(9,*) a,d,z(i),m1,m1,m(i),m1,m1,e1,e1,e(i),e1,e1,k1,k1,k1,k1,k1,k1,k1,k(i)
                mc(i)=m(i)-e(i)

                dl=distance(z(i))*(1.+z(i))
                da=distance(z(i))/(1.+z(i))

                Mag = mc(i) - 25. -5. * log10(dl) + .1 + k(i)   ! el .1 es de la correcion ab del program anterior..
               
                if (Mag <= -23. .or. Mag >= -16.) cycle
                dlmin=10**(((-rmin+Mag)+25.)/(-5.))
                dlmax=10**(((-rmax+Mag)+25.)/(-5.))
               
                zmin=dlum(dlmin)
                zmax=dlum(dlmax)
                 
                call qromb(func, zmin, zmax, vol) 
                write(11,*) Mag, 1./vol

        enddo
        close(9)
        close(10)
        close(11)     
      

endprogram

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

function distance(z0)
        use cosm
        implicit none
        real :: dist, z0, distance
        call locate(zlist,nlist,z0,j)
        distance= cd(j) + (cd(j+1)-cd(j))*(z0-zlist(j))/(zlist(j+1)-zlist(j))
        return
endfunction   
       
function dlum(dlc)
        use cosm
        implicit none
        real dlc, dlum 
         
        call locate(dllist, nlist, dlc, j)
        dlum= zlist(j) + ((zlist(j+1)-zlist(j))/(dllist(j+1)-dllist(j)))*(dlc-dllist(j))
        return
endfunction

function func(x)
        implicit none
        real num, den, func,x, distance
        num=distance(x)
        den=((.3*1+x)**3+(1.-0.3-0.7)*(1.+x)**2+.7)**.5
        func=num/den
endfunction

      SUBROUTINE qromb(func,a,b,ss)
      INTEGER JMAX,JMAXP,K,KM
      REAL a,b,func,ss,EPS
      EXTERNAL func
      PARAMETER (EPS=1.e-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
!CU    USES polint,trapzd
      INTEGER j
      REAL dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call trapzd(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)
11    continue
      pause 'too many steps in qromb'
      END 

      include 'polint.f'
      include 'trapzd.f'
