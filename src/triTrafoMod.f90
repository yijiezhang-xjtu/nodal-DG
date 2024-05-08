module triTrafoMod
    ! module to transform between different tri representations
    use constantsMod

    implicit none

    contains

    subroutine rsToAb(r,s,a,b)
        ! transform from reference triange into reference quad
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: r,s
        real(kind=CUSTOM_REAL), dimension(:), intent(out) :: a,b
        integer :: i

        do i=1,size(r)
            if (abs(s(i) - 1.) > EPS) then
                a(i) = 2.0*(1.0+r(i))/(1.0-s(i))-1.0
            else
                a(i) = -1.0
            end if
        end do
        b(:)=s(:)
    end subroutine rsToAb

    subroutine xyToRs(x,y,r,s)
        ! transform from equilateral tri into reference tri
        real(kind=CUSTOM_REAL), dimension(:), intent(out) :: r,s
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: x,y
        real(kind=CUSTOM_REAL), dimension(Np) :: L1,L2,L3

        L1(:) = (sqrt(3.0)*y(:)+1.0)/3.0
        L2(:) = (-3.0*x(:) - sqrt(3.0) * y(:) + 2.0)/6.0
        L3(:) = ( 3.0*x(:) - sqrt(3.0) * y(:) + 2.0)/6.0
        r(:) = -L2(:) + L3(:) - L1(:)
        s(:) = -L2(:) - L3(:) + L1(:)
    end subroutine xyToRs
end module triTrafoMod
