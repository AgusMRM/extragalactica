 module cosm
        implicit none
        integer i, j
        integer, parameter :: nlist=6000, ngalaxy=544449
        real, dimension(ngalaxy) :: z,m,r,k,mc,e
        real, dimension(nlist) :: zlist, cd, ad, ld, dllist
        real :: dl, da, Mag,x
endmodule

PROGRAM QSO
        use cosm
        IMPLICIT NONE
        REAL, ALLOCATABLE :: xm0(:), xm1(:), xm2(:)
        INTEGER bin, nbin, num
        REAL :: z, magi, Mi, vol, zmin, zmax, dlmin, dlmax
        REAL :: imin, imax, dlum, func

        nbin=10

        ALLOCATE(xm0(nbin))
        ALLOCATE(xm1(nbin))
        ALLOCATE(xm2(nbin))
               open(12,file='cosmology.dat',status='old',action='read')

!-----------------------armo vector cosmologico..............
        do i=1,6
        read(10,*)
        enddo
        do i=1,nlist
                read(12,*) zlist(i), cd(i), ad(i), dllist(i)       
        enddo
        close(12)
!----- hago una seleccion de magnitudes i para los quasar------------

        open(10,file='qso.dat',status='unknown')
        open(11,file='seleccionqso.dat',status='unknown')
        read(10,*) 
        num=0
        do i=1,38060
                read(10,*) z, magi, Mi
                if (magi>= 15. .and. magi<=19.) then
                        dlmin=10**(((-imin+Mi)+25.)/(-5.))
                        dlmax=10**(((-imax+Mi)+25.)/(-5.))

                        zmin=dlum(dlmin)
                        zmax=dlum(dlmax)

                        call qromb(func, zmin, zmax, vol)
                        num=num + 1
                        write(11,*) Mi, 1./vol
                endif
        enddo
        close(10)
        close(11)
!-------------------------------------------------------------------
!----- calculo los pesos--------------------------------------------         


       

        
ENDPROGRAM

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
function func(x)
        implicit none
        real num, den, func,x, distance
        num=distance(x)
        den=((.3*1+x)**3+(1.-0.3-0.7)*(1.+x)**2+.7)**.5
        func=num/den
endfunction
function dlum(dlc)
        use cosm
        implicit none
        real dlc, dlum 
         
        call locate(dllist, nlist, dlc, j)
        dlum= zlist(j) + ((zlist(j+1)-zlist(j))/(dllist(j+1)-dllist(j)))*(dlc-dllist(j))
        return
endfunction

function distance(z0)
        use cosm
        implicit none
        real :: dist, z0, distance
        call locate(zlist,nlist,z0,j)
        distance= cd(j) + (cd(j+1)-cd(j))*(z0-zlist(j))/(zlist(j+1)-zlist(j))
        return
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
