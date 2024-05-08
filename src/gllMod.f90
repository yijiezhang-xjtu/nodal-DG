module gllMod
    use constantsMod
    ! module to generate gll points
    implicit none

    contains

    function getGll()
        ! function to get gll points for order NGLL-1
        real(kind=CUSTOM_REAL), dimension(NGLL) :: x,GetGll
        real(kind=CUSTOM_REAL), dimension(NGLL,NGLL):: V
        real(kind=CUSTOM_REAL) :: xold
        integer :: i,k,p

        p = NGLL-1

        do i=1,NGLL
            x(i) = cos(pi*(i-1)/p)
        end do

        ! Newton-Raphson iteration
        do i=1,NGLL
            xold = 2.
            do while (abs(x(i)-xold) > eps)
                xold = x(i)
                V(i,1) = 1.
                V(i,2) = x(i)
                do k=2,p
                    V(i,k+1) = ((2.*k-1.)*x(i)*V(i,k)-(k-1.)*V(i,k-1))/k
                end do
                x(i) = xold-(x(i)*V(i,NGLL)-V(i,p))/(NGLL*V(i,NGLL))
            end do
        end do
        GetGLL=x
    end function getGll
end module gllMod
