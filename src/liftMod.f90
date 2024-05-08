module liftMod
    ! module to compute the lifting from surf to vol for dg
    use constantsMod
    use vandermondeMod
    use errorMessage
    implicit none

    contains

    subroutine lift2d(lift,fmask,v,r,s,errmsg)
        type(error_message) :: errmsg
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: r,s
        real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: v
        integer, dimension(NGLL,3), intent(in) :: fmask
        real(kind=CUSTOM_REAL), dimension(Np,3*NGLL), intent(out) :: lift
        real(kind=CUSTOM_REAL), dimension(Np,3*NGLL) :: temp
        real(kind=CUSTOM_REAL), dimension(Np,3*NGLL) :: emat
        real(kind=CUSTOM_REAL), dimension(NGLL) :: faceR
        real(kind=CUSTOM_REAL), dimension(NGLL,NGLL) :: v1d, masse
        integer :: i,j,k,is


        character(len=6) :: myname = "lift2d"

        call addTrace(errmsg, myname)

        emat = 0.0
        lift = 0.0

        do is=1,3
            if (is<3) then
                faceR(:) = r(fmask(:,is))
            else
                faceR(:) = s(fmask(:,is))
            endif
            call vdm1d(v1d,faceR)

            masse=matmul(v1d,transpose(v1d))
            call invert(masse, errmsg)

            do i=1,NGLL
                k=1
                do j=(is-1)*NGLL+1,is*NGLL
                    emat(fmask(i,is),j)=masse(i,k)
                    k=k+1
                end do
            end do
        enddo

        ! build lift matrix
        temp = matmul(transpose(v),emat)
        lift= matmul(v,temp)
    end subroutine lift2d
end module liftMod
