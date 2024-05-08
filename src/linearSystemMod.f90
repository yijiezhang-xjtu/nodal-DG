module linearSystemMod
    ! module to solve linear systems
    use constantsMod
    use errorMessage

    implicit none

    contains

    subroutine solveLinearSystemQR(A,b,x,nrow,ncol, errmsg, typesize)
        ! solve a overdetermined linear system with QR methode (lapack SGELS) A*x=b
        type(error_message) :: errmsg
        real(kind=CUSTOM_REAL), dimension(:,:) :: A
        real(kind=CUSTOM_REAL), dimension(:) :: x
        real(kind=CUSTOM_REAL), dimension(:) :: b
        integer :: nrow,ncol
        real(kind=CUSTOM_REAL), dimension(nrow,ncol) :: A_tmp
        integer :: LDB,LWORK,INFO, typesize
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: x_tmp
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: WORK

        character(len=19) :: myname = "solveLinearSystemQR"

        call addTrace(errmsg, myname)

        A_tmp = A

        LDB = max(nrow,ncol)
        allocate(x_tmp(LDB))
        x_tmp(:) = 0.
        x_tmp(1:nrow) = b(:)
        allocate(WORK(1));
        LWORK = -1

        if (typesize == SIZE_REAL) then
            call SGELS('N', nrow, ncol, 1, A_tmp, nrow, x_tmp, LDB, WORK, LWORK, INFO)
        else
            call DGELS('N', nrow, ncol, 1, A_tmp, nrow, x_tmp, LDB, WORK, LWORK, INFO)
        endif
        if(INFO/=0) then
            call add(errmsg, 2, "Error in solveLinearSystemQR: LWORK", myname)
            if (.level.errmsg == 2) then; call print(errmsg, .true.); stop; endif
        end if

        LWORK = int(WORK(1))
        deallocate(WORK)
        allocate(WORK(LWORK))

        if (typesize == SIZE_REAL) then
            call SGELS('N', nrow, ncol, 1, A_tmp, nrow, x_tmp, LDB, WORK, LWORK, INFO)
        else
            call DGELS('N', nrow, ncol, 1, A_tmp, nrow, x_tmp, LDB, WORK, LWORK, INFO)
        endif
        if(INFO/=0) then
            call add(errmsg, 2, "Error in solveLinearSystemQR: Solve", myname)
            if (.level.errmsg == 2) then; call print(errmsg, .true.); stop; endif
        end if
        x=x_tmp(1:ncol)
    end subroutine solveLinearSystemQR
end module linearSystemMod
