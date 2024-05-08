module materialsMod

    use parameterMod
    use constantsMod
    use errorMessage
    use linearSystemMod

    implicit none

    interface setupMaterials
        module procedure setupMaterialsElastic
        !module procedure setupMaterialsPoroelastic
    end interface setupMaterials

    type materialVar
        real(kind=custom_real) :: vp         !P-wave velocity
        real(kind=custom_real) :: vs         !S-wave velocity
        real(kind=custom_real) :: vpu        !unrelaxed P-wave velocity
        real(kind=custom_real) :: vsu        !unrelaxed S-wave velocity
        real(kind=custom_real) :: lambda     !1st Lame parameter
        real(kind=custom_real) :: mu         !2nd Lame parameter
        real(kind=custom_real) :: rho        !Density
        real(kind=custom_real) :: imp_vp     !Impedance with respect to vp
        real(kind=custom_real) :: imp_vs     !Impedance with respect to vs
        real(kind=custom_real) :: qp         !Qualityfactor for P-wave speed
        real(kind=custom_real) :: qs         !Qualityfactor for S-wave speed
        real(kind=CUSTOM_REAL), dimension(nMB) :: ylambda   !Anelastic coefficient regarding the 1st Lame parameter
        real(kind=CUSTOM_REAL), dimension(nMB) :: ymu       !Anelastic coefficient regarding the 2st Lame parameter
        real(kind=CUSTOM_REAL), dimension(nMB) :: omegaL    !Relaxation frequencies of the Maxwellbodys
    end type materialVar

    type materialIndizes
        integer, dimension(:), allocatable :: type           !Index of the material-type that is used in each element
        integer, dimension(:), allocatable :: pml            !Index on wether an element belongs to a pml or not
    end type materialIndizes
    contains

    subroutine setupMaterialsElastic(par, mat, matInd, nelem, iregion, errmsg)
        !This subroutine initializes the material variables for the different cases (elastic, attenuation)
        !input
        type(parameterVar) :: par
        type(error_message) :: errmsg
        integer :: nelem
        integer, dimension(:), allocatable :: iregion
        !in/out
        type(materialIndizes) :: matInd
        type(materialVar), dimension(:), allocatable :: mat
        !local
        character(len=14) :: myname = "setupMaterials"

        call addTrace(errmsg, myname)

        ! materials values
        call readMaterialProperties(par, mat, 19, iregion, trim('mesh/matprop'), errmsg)

        ! materials array
        call readMaterialFile(par, matInd, nelem, 19, trim('mesh/mat'), errmsg)
    end subroutine setupMaterialsElastic

    subroutine readMaterialProperties(par, mat, lu, iregion, filename, errmsg)
        !Subroutine to read properties of the materials from the cubit/trelis-generated file
        !input
        type(parameterVar) :: par
        type (error_message) :: errmsg
        integer, intent(in) :: lu
        character(len=*), intent(in) :: filename
        !in/output
        type(materialVar), dimension(:), allocatable :: mat
        integer, dimension(:), allocatable :: iregion              !Some parameter for the external model. Not sure if is needed anymore... TM TM
        !local
        integer :: i
        integer :: dummy
        integer :: ios
        character(len=22) :: myname = "readMaterialProperties"

        call addTrace(errmsg, myname)

        open(unit=lu,file=trim(filename), status='old', iostat = ios)
        if (ios /= 0) then
            call add(errmsg,2,'could not open: '//trim(filename),myname)
            call print(errmsg)
            stop
        endif

        read(lu,*) par%matn
        allocate(mat(par%matn))
        allocate(iregion(par%matn))

        do i = 1, par%matn
           !vp,vs,rho,qp,qs
           read(lu,*) iregion(i), dummy, mat(i)%vp, mat(i)%vs, mat(i)%rho, mat(i)%qp ,mat(i)%qs
           !compute mu and lambda and the impedances
           mat(i)%mu = mat(i)%vs**2 * mat(i)%rho
           mat(i)%lambda = mat(i)%vp**2 * mat(i)%rho - 2*mat(i)%mu
           mat(i)%imp_vp = mat(i)%vp * mat(i)%rho
           mat(i)%imp_vs = mat(i)%vs * mat(i)%rho
           mat(i)%vpu = mat(i)%vp
           mat(i)%vsu = mat(i)%vs
        end do
        close(lu)
        deallocate(iregion)
    end subroutine readMaterialProperties

    

    !Subroutine to read the material index file. Preferable this file should be in materials mod.
    subroutine readMaterialFile(par, matInd, nelem, lu, filename, errmsg)
        !input
        type(parameterVar) :: par
        type(error_message) :: errmsg
        integer :: lu
        integer :: nelem
        character(len=*) :: filename
        !in/out
        type(materialIndizes) :: matInd
        !loacl
        integer :: i
        integer :: dummy
        integer :: ios
        character(len=16) :: myname = "readMaterialFile"

        call addTrace(errmsg, myname)

        if (par%log) write(*,'(a80)') "|                              read material file                              |"
        allocate(matInd%type(nelem))
        allocate(matInd%pml(nelem))
        !filename=trim('mesh/mat')
        open(unit=lu,file=trim(filename), status='old', iostat = ios)
        if (ios /= 0) then
            call add(errmsg,2,'could not open: '//trim(filename),myname)
            call print(errmsg)
            stop
        endif
        do i=1,nelem
           ! materials
           read(lu,*) dummy,matInd%type(i),matInd%pml(i)   !read Material index -> which material is used in that element, and the index identifiing pmls
        end do
        close(lu)
    end subroutine readMaterialFile

    subroutine calcAPAM(A,AP,AM,vmax,vmin)
        !input
        real(kind=custom_real), dimension(:,:), intent(in) :: A
        !output
        real(kind=custom_real), dimension(:,:), intent(out) :: AP, AM
        real(kind=custom_real), intent(out) :: vmax, vmin
        !local variables
        type(error_message) :: errmsg
        character(len=8) :: myname = 'calcAPAM'
        character(len=100) :: errstr
        real(kind=custom_real), dimension(:,:), allocatable :: Awork
        real(kind=custom_real), dimension(:), allocatable :: WR, WI    !contain real and imaginary part of eigenvalues, respectively
        real(kind=custom_real), dimension(:,:), allocatable :: VR      !matrix containing the eigenvectors as columns
        real(kind=custom_real), dimension(:,:), allocatable :: LambdaP, LambdaM  !matrix containing the eigenvalues on the diagonal elements
        real(kind=custom_real), dimension(:,:), allocatable :: VL,invVR
        real(kind=custom_real), dimension(:), allocatable :: work
        integer :: N, lwork, info, i
        integer, dimension(:), allocatable :: ipiv
        integer, dimension(2) :: shapeA

        shapeA=shape(A)
        N=shapeA(1)
        allocate(Awork(N,N))
        allocate(ipiv(N))
        allocate(WR(N))
        allocate(WI(N))
        allocate(VR(N,N))
        allocate(LambdaP(N,N))
        allocate(LambdaM(N,N))
        allocate(invVR(N,N))
        Awork = A

        !calculate eigenvectors and eigenvalues
        if (custom_real==size_double) then
            !do workspace query
            allocate(work(1)); lwork = -1
            !SYNTAX: dgeev( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )   (LAPACK: http://www.netlib.org/lapack/double/dgeev.f)
            call dgeev('N','V', N, Awork, N, WR, WI, VL, N, VR, N, work, lwork, info)
            if (info/=0) then
                write(errstr,*) 'Workspace query failed: LAPACK routine DGEEV returned INFO = ',info
                call add(errmsg,2,trim(errstr),myname)
                stop
            endif
            lwork = int(work(1))
            deallocate(work); allocate(work(lwork))
            !compute eigenvalues and eigenvectors
            call dgeev('N','V', N, Awork, N, WR, WI, VL, N, VR, N, work, lwork, info)
            if (info/=0) then
                write(errstr,*) 'Computation of eigenvalues and -vectors failed: LAPACK routine DGEEV returned INFO = ',info
                call add(errmsg,2,trim(errstr),myname)
                stop
            endif
            deallocate(work)
        elseif (custom_real==size_real) then
            !do workspace query
            allocate(work(1)); lwork = -1
            !SYNTAX: sgeev( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )   (LAPACK: http://www.netlib.org/lapack/double/dgeev.f)
            call sgeev('N','V', N, Awork, N, WR, WI, VL, N, VR, N, work, lwork, info)
            if (info/=0) then
                write(errstr,*) 'Workspace query failed: LAPACK routine SGEEV returned INFO = ',info
                call add(errmsg,2,trim(errstr),myname)
                stop
            endif
            lwork = int(work(1))
            deallocate(work); allocate(work(lwork))
            !compute eigenvalues and eigenvectors
            call sgeev('N','V', N, Awork, N, WR, WI, VL, N, VR, N, work, lwork, info)
            if (info/=0) then
                write(errstr,*) 'Computation of eigenvalues and -vectors failed: LAPACK routine SGEEV returned INFO = ',info
                call add(errmsg,2,trim(errstr),myname)
                stop
            endif
            deallocate(work)
        else
            call add(errmsg,2,'CUSTOM_REAL is neither SIZE_REAL nor SIZE_DOUBLE',myname)
            stop
        endif

        LambdaP = 0.
        LambdaM = 0.
        do i=1,N
            if (WR(i) > 0) then
                LambdaP(i,i)=WR(i)
            else
                LambdaM(i,i)=WR(i)
            endif
        enddo

        invVR = VR
        if (custom_real==size_double) then
            call dgetrf(N, N, invVR, N, ipiv, info)
            !do workspace query
            allocate(work(1)); lwork = -1
            call dgetri(N, invVR, N, ipiv, work, lwork, info)
            lwork = int(work(1))
            deallocate(work); allocate(work(lwork))
            call dgetri(N, invVR, N, ipiv, work, lwork, info)
            deallocate(work)
        elseif (custom_real==size_real) then
            call sgetrf(N, N, invVR, N, ipiv, info)
            !do workspace query
            allocate(work(1)); lwork = -1
            call sgetri(N, invVR, N, ipiv, work, lwork, info)
            lwork = int(work(1))
            deallocate(work); allocate(work(lwork))
            call sgetri(N, invVR, N, ipiv, work, lwork, info)
            deallocate(work)
        else
            call add(errmsg,2,'CUSTOM_REAL is neither SIZE_REAL nor SIZE_DOUBLE',myname)
            stop
        endif

        AP = 2.*matmul(VR,matmul(LambdaP,invVR))
        AM = 2.*matmul(VR,matmul(LambdaM,invVR))
        vmax = maxval(WR)
        vmin = 300000.
        do i=1,size(WR)
            if (WR(i) > 0) vmin = min(vmin,WR(i))
        enddo

        if (allocated(Awork)) deallocate(Awork)
        if (allocated(ipiv)) deallocate(ipiv)
        if (allocated(WR)) deallocate(WR)
        if (allocated(WI)) deallocate(WI)
        if (allocated(VR)) deallocate(VR)
        if (allocated(LambdaP)) deallocate(LambdaP)
        if (allocated(LambdaM)) deallocate(LambdaM)
        if (allocated(invVR)) deallocate(invVR)
    end subroutine calcAPAM

end module materialsMod
