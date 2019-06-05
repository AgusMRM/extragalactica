module cosmo
implicit none
integer, parameter :: nlist=6000
real, dimension(nlist) :: zlist,ad,ld,cd
endmodule cosmo

program distancias
        use cosmo
        implicit none
        integer :: i, j
        real :: z, distance,dl,da

        open(10,file='cosmology.dat',status='old',action='read')
        do i=1,6
        read(10,*) 
        enddo
        do i=1,nlist
        read(10,*) zlist(i),cd(i),ad(i),ld(i)
        enddo
        close(10)
        z=.3
        dl=distance(z)*(1.+z)
        da=distance(z)/(1.+z)
        print*, dl, da
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

function distance(z)
        use cosmo
        implicit none
        integer :: j
        real :: distance, z 
        call locate(zlist,nlist,z,j)
                      
        distance = cd(j)+(cd(j+1)-cd(j))*(z-zlist(j))/(zlist(j+1)-zlist(j))

        return 
endfunction
