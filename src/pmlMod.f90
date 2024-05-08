module pmlMod
    use constantsMod

    implicit none

    contains

    subroutine dampingProfile(ddx,ddz,alphax,alphaz,kx,kz,vx,vz,vp,xmin,xmax,zmin,zmax,pml_delta,rc,amax,kmax,iv,pmlloc)
        real(kind=CUSTOM_REAL), dimension(:) :: vx,vz
        real(kind=CUSTOM_REAL), dimension(:) :: ddx,ddz,alphax,alphaz,kx,kz
        real(kind=CUSTOM_REAL), dimension(Np) :: ar,br
        real(kind=CUSTOM_REAL) :: pml_delta,vp,rc,amax,kmax, xmin, xmax, zmin, zmax
        integer, dimension(2) :: pmlloc
        integer, dimension(NP) :: iv

        if ( pmlloc(1) == -1) then !xmin
            ar=sqrt((xmin+pml_delta)**2)
            br=sqrt(vx(iv)**2)
            br=abs(ar-br)
            ddx(iv) = -3.0 / (2*pml_delta) * 5000 * alog(real(rc))*( br/pml_delta )**2
            !ddx(iv) = -40.0 / (3*pml_delta * 4000) * alog(real(rc))*( br/pml_delta )**2
            alphax(iv) =amax*(1-br/pml_delta)
            kx(iv)=1+(kmax-1)*(br/pml_delta)**2
        else if ( pmlloc(1) == 1 ) then !xmax
            ar=sqrt((xmax-pml_delta)**2)
            br=sqrt(vx(iv)**2)
            br=abs(ar-br)
            ddx(iv) = -3.0 / (2*pml_delta) * 5000 * alog(real(rc))*( br/pml_delta )**2
            !ddx(iv) = -40.0 / (3*pml_delta * 4000) * alog(real(rc))*( br/pml_delta )**2
            alphax(iv) = amax*(1- br/pml_delta )
            kx(iv)=1+(kmax-1)*(br/pml_delta)**2
        endif
        if ( pmlloc(2) == -1 ) then !zmin
            ar=sqrt((zmin+pml_delta)**2)
            br=sqrt(vz(iv)**2)
            br=abs(ar-br)
            ddz(iv) = -3.0 / (2*pml_delta) * 5000 * alog(real(rc))*( br/pml_delta )**2
            !ddz(iv) = -40.0 / (3*pml_delta * 4000) * alog(real(rc))*( br/pml_delta )**2
            alphaz(iv) = amax*(1- br/pml_delta )
            kz(iv)=1+(kmax-1)*(br/pml_delta)**2
        else if ( pmlloc(2) == 1 ) then !zmax
            ar=sqrt((zmax-pml_delta)**2)
            br=sqrt(vz(iv)**2)
            br=abs(ar-br)
            ddz(iv) = -3.0 / (2*pml_delta) * 5000 * alog(real(rc))*( br/pml_delta )**2
            !ddz(iv) = -40.0 / (3*pml_delta * 4000) * alog(real(rc))*( br/pml_delta )**2
            alphaz(iv) = amax*(1- br/pml_delta )
            kz(iv)=1+(kmax-1)*( br/pml_delta )**2
        endif
    end subroutine dampingProfile
end module pmlMod
