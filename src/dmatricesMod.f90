module dmatricesMod
    ! module to calculate the differentiation matrices on the simplex at (r,s)
    use constantsMod
    use vandermondeMod
    use errorMessage

    implicit none

    contains

    subroutine dmatrices2d(dr,ds,r,s,v, errmsg)
        type(error_message) :: errmsg
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: r,s
        real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: v
        real(kind=CUSTOM_REAL), dimension(size(v(:,1)),Np) :: dr,ds,w,vr,vs
        character(len=11) :: myname = "dmatrices2d"

        call addTrace(errmsg, myname)

        call gradVdm2D(vr,vs,r,s)
        call invVdm2D(v,w,0, errmsg)

        dr=matmul(vr,w)
        ds=matmul(vs,w)
    end subroutine dmatrices2d
end module dmatricesMod
