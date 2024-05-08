module simplexMod
    ! Evaluate 2D orthonormal polynomial on simplex at (a,b) oft order (i,j)
    use constantsMod
    use jacobiMod

    implicit none

    contains

    subroutine simplex2DP(Pr,a,b,i,j)
        real(kind=CUSTOM_REAL), dimension(:), intent(out) :: Pr
        real(kind=CUSTOM_REAL), dimension(:), intent(in):: a,b
        real(kind=CUSTOM_REAL), dimension(size(a)) :: h1,h2

        integer, intent(in) :: i,j

        call jacobiP(h1,a,zero,zero,i)
        call jacobiP(h2,b,two*i+one,zero,j)

        Pr(:) = sqrt(2.0)*h1*h2*(1.0-b)**i
    end subroutine simplex2DP

    subroutine gradSimplex2DP(dPdr,dPds,a,b,id,jd)
        ! computes the differentiation matrices Dr and Ds
        real(kind=CUSTOM_REAL), dimension(:), intent(in) ::a,b
        integer, intent(in) :: id,jd
        real(kind=CUSTOM_REAL), dimension(size(a)), intent(out) :: dPdr,dPds
        real(kind=CUSTOM_REAL), dimension(size(a)) :: fa,dfa
        real(kind=CUSTOM_REAL), dimension(size(b)) :: gb,dgb
        real(kind=CUSTOM_REAL), dimension(size(a)) :: temp

        call jacobiP(fa,a,zero,zero,id)
        call gradJacobiP(dfa,a,zero,zero,id)
        call jacobiP(gb,b,two*id+one,zero,jd)
        call gradJacobiP(dgb,b,two*id+one,zero,jd)

        ! r-deriverate
        !d/dr = da/dr d/da + db/dr d/db = (2/(1-s)) d/da = (2/(1-b)) d/da

        dPdr(:)=dfa(:)*gb(:)
        if (id>0) then
            dPdr(:)=dPdr(:)*((0.5*(1.0-b(:)))**(id-1.0))
        end if

        ! s-deriverate
        !d/ds = ((1+a)/2)/((1-b)/s) d/da + d/db

        dPds(:) = dfa(:)*(gb(:)*(0.5*(1.0+a(:))))
        if (id>0) then
            dPds(:)=dPds(:)*((0.5*(1.0-b(:)))**(id-1.0))
        end if

        temp(:) = dgb(:)*((0.5*(1.0-b(:)))**id)
        if (id>0) then
            temp(:)=temp(:)-0.5*id*gb(:)*((0.5*(1.0-b(:)))**(id-1.0))
        end if
        dPds(:) = dPds(:)+fa(:)*temp(:)

        !normalize
        dPdr(:) = 2.0**(id+0.5)*dPdr(:)
        dPds(:) = 2.0**(id+0.5)*dPds(:)
    end subroutine gradSimplex2DP
end module simplexMod
