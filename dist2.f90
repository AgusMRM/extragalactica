module cosmo
        implicit none
        integer, parameter :: nlist=6000, ngalaxy=50000
        real, dimension(nlist) :: zlist, ad, ld, cd
        character(len=18) :: bestobjid

endmodule

program dist
        use cosmo
        implicit none
        integer :: i, j, k
        real :: z, distance, dl, da, l50, l90,pi
        real :: pmu,pmg,pmr,pmi,pmz,mmu,mmg,mmr,mmi,mmz,eu,eg,er,ei,ez
        real :: ra, dec, s, c, fracdev, veldisp, ur, gr, surf
        real, dimension(5) :: ab, m_ab, Ext, m_app, m_apm, magnitud_m, magnitud_p, p50, p90
        data ab/.036,.012,.1,.028,.4/
        open(10,file='cosmology.dat', status='old', action='read')
        do i=1,6
        read(10,*)
        enddo
        
        pi=acos(-1.)
        do i=1,nlist
        read(10,*) zlist(i),cd(i),ad(i),ld(i)
        enddo
        close(10)
        
        open(11,file='tabla.dat',status='old',action='read')
        open(12,file='resultados.dat',status='unknown')
        read(11,*)
        do i=1,ngalaxy
                read(11,*) ra, dec, z, bestobjid, (m_app(k), k=1,5), (m_apm(k), k=1,5), &
                         (ext(k), k=1,5), (p50(k), k=1,5), (p90(k), k=1,5), fracDev, veldisp 
                 do k=1,5
                        m_app(k) = m_app(k) - ext(k)
                        m_apm(k) = m_apm(k) - ext(k)

                enddo

                if (m_app(3) < 14.5 .or. m_app(3) > 17.77) cycle
                if (p50(3) < 1.5) cycle
                dl=distance(z)*(1.+z)
                da=distance(z)/(1.+z)
                
                do k=1,5   
                magnitud_p(k) = m_app(k) - 25. -5.*log10(dl) +ab(k)
                magnitud_m(k) = m_apm(k) - 25. -5.*log10(dl) + ab(k)
                
                enddo
                C=p90(3)/p50(3)
                UR=magnitud_m(1)-magnitud_m(3)
                GR=magnitud_m(2)-magnitud_m(3)
                l50=da*p50(3)*pi/(180.*60.*60.)*1000.
                l90=da*p90(3)*pi/(180.*60.*60.) *1000.
                
                surf=m_app(3)+ab(3)+2.5*log10(p50(3)**2)+ 0.26  
                write(12,120) ra,dec,z ,magnitud_p,magnitud_m , C , UR, GR , l50, l90 , fracdev , veldisp, surf  
        enddo
        close(11)
        close(12)
120 format(2(f8.4,1x),f6.4,1x, 10(f7.3,1x), 1x, 3(f6.3,1x), 2(f6.3,1x), f5.3, 1x, f8.2, 1x, f7.3)
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

function magnitud(m_ap, ext, dl, ab)
        implicit none
        real, dimension(5) :: m_ap, ext, ab, magnitud
        real :: dl
        magnitud=(m_ap-ext) - 25. -5.*log10(dl) + ab
        return        
endfunction
