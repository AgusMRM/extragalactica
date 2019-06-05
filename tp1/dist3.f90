
module cosmo
        implicit none
        integer, parameter :: nlist=6000, ngalaxy=50000
        
        real, dimension(nlist) :: zlist, ad, ld, cd
        character :: bestobjid
        
endmodule

program dist
        use cosmo
        implicit none
        real :: z, distance, dl, da
        real :: pmu,pmg,pmr,pmi,pmz,mmu,mmg,mmr,mmi,mmz,eu,eg,er,ei,ez
        real :: ra, dec, s
        real, dimension(5) :: ab, m_ab, Ext, m_ap
        real, dimension(5) :: magnitud, p50, p90
        integer :: i,j,k
        data ab/.036,.012,.1,.028,.4/
        

        open(10,file='cosmology.dat', status='old', action='read')
        do i=1,6
        read(10,*)
        enddo

        do i=1,nlist
        read(10,*) zlist(i),cd(i),ad(i),ld(i)
        enddo
        close(10)
        
        open(11,file='tabla.dat',status='old',action='read')
        open(12,file='resultados.dat',status='unknown')
        read(11,*)
        do i=1,ngalaxy
                read(11,*) ra,dec,z, bestobjid, m_ap(1),m_ap(2), m_ap(3), m_ap(4), m_ap(5), mmu, mmg, mmr, mmi,mmz, &
                           ext(1),ext(2),ext(3),ext(4),ext(5), p50(1),p50(2),p50(3),p50(4),p50(5),p90(1),p90(2),p90(3), &
                           p90(4),p90(5)
                dl=distance(z)*(1.+z)
                da=distance(z)/(1.+z)
                
                call magnitude(m_ap,ext,dl,magnitud) 
                                 
                 
                write(12,*) magnitud
                
        enddo
        close(11)
        close(12)
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
        real ::distance, z
        call locate(zlist,nlist,z,j)
        distance = cd(j) + (cd(j+1)-cd(j))*(z-zlist(j))/(zlist(j+1)-zlist(j))

        return
endfunction
subroutine Magnitude(m_ap,ext,dl,magnitud)
        use cosmo
        implicit none
        real, dimension(5) :: magnitud, m_ap,ext, ab
        real dl
        integer k
        do k=1,5
         magnitud(k)=(m_ap(k)-ext(k))-25.-5*log10(dl)+ab(k)
                
                
        enddo
endsubroutine
