module geometricFactorsMod
    ! module to compute the metric elements for local mappings of the element
    use constantsMod
    use errorMessage
    implicit none

    contains

    subroutine geometricFactors2d(rx,sx,ry,sy,J,x,y,Dr,Ds, errmsg)
        !Input
        type(error_message) :: errmsg
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: x,y
        real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: Dr,Ds
        !Output
        real(kind=CUSTOM_REAL), dimension(size(x)) :: rx,sx,ry,sy,J
        !Local
        real(kind=CUSTOM_REAL), dimension(size(x)) :: xr,xs,yr,ys
        integer :: i
        character(len=18) :: myname = 'geometricFactors2d'

        xr = matmul(Dr,x)
        xs = matmul(Ds,x)
        yr = matmul(Dr,y)
        ys = matmul(Ds,y)

        J = -xs*yr + xr*ys

        do i=1,size(J)
            if (J(i) .le. 0.0) then
                call add(errmsg, 2, "Jacobian negative or zero!", myname)
                if (.level.errmsg == 2) then; call print(errmsg, .true.); stop; endif
            end if
        end do

        rx = ys/J
        sx = -yr/J
        ry = -xs/J
        sy = xr/J
    end subroutine geometricFactors2d
end module geometricFactorsMod
