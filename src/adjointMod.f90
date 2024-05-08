module adjointMod
    !use constantsMod
    use meshMod
    use parameterMod
    use mpiMod
    use timeSeries
    use dtMod
    use collectMovieMod
    implicit none

    interface makeTimeSeries
        module procedure makeTimeSeriesDP
        module procedure makeTimeSeriesSP
    end interface
    interface cubicInterpolation
        module procedure cubicInterpolationArray
        module procedure cubicInterpolationSnglVal
    end interface


    type :: invVar
        real(kind=CUSTOM_REAL) :: upperfreq                                         ! cutoff frequency for low-pass filtering
        integer :: inv_type                                                         ! flag for the search type of the search direction (1: Steepest descent, 2: Conjugate gradient, 3: BFGS)
        real(kind=CUSTOM_REAL), dimension(3) :: stepsize                            ! Length of test step sizes to be used
        real(kind=CUSTOM_REAL) :: min_step, max_step, min_step_BFGS, max_step_BFGS  ! minimum and maximum values for test step sizes during search for best step size
        real(kind=CUSTOM_REAL) :: min_vp, max_vp, min_vs, max_vs                    ! Minimum and maximum valocity values for the reconstructed models
    end type invVar
contains

    subroutine readInvfile(this,myrank,par, errmsg)
    ! Read the file "data/invpar" with information regarding the inversion
        implicit none
        type(error_message) :: errmsg
        type (InvVar) :: this
        type(parameterVar):: par
        integer ::myrank

        !local variables
        character(len=80) :: filename
        integer :: ier

        filename=trim('data/invpar')
        open(unit=19, file=trim(filename), status='old', action='read', iostat=ier)
        if (myrank==0) then
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a26, a14, a12, a28)") "|                          ","Begin reading ", filename, "...                        |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
        endif
        ! Cycle to read the parameters form the parameter file. The order of appearence of the parameters is not important

        call readFloatPar(this%upperfreq, "upperfreq", filename, 0, errmsg)
        call readIntPar(this%inv_type, "inv_type", filename, 0, errmsg)
        call readFloatPar(this%stepsize(1), "step_1", filename, 0, errmsg)
        call readFloatPar(this%stepsize(2), "step_2", filename, 0, errmsg)
        call readFloatPar(this%stepsize(3), "step_3", filename, 0, errmsg)
        call readFloatPar(this%min_step, "min_step", filename, 0, errmsg)
        call readFloatPar(this%max_step, "max_step", filename, 0, errmsg)
        call readFloatPar(this%min_step_BFGS, "min_step_BFGS", filename, 0, errmsg)
        call readFloatPar(this%max_step_BFGS, "max_step_BFGS", filename, 0, errmsg)
        call readFloatPar(this%min_vp, "min_vp", filename, 0, errmsg)
        call readFloatPar(this%max_vp, "max_vp", filename, 0, errmsg)
        call readFloatPar(this%min_vs, "min_vs", filename, 0, errmsg)
        call readFloatPar(this%max_vs, "max_vs", filename, 0, errmsg)


        ! print information read from file to screen
        if (par%log.and.myrank==0) then
            write (*,"(a40, f10.1, a30)") "|                    Highest frequency: ", this%upperfreq, "                             |"
            write (*,"(a40, i10, a30)")   "|                     Inversion method: ", this%inv_type,  "                             |"
            write (*,"(a40, f10.3, a30)") "|                            Stepsizes: ", this%stepsize(1),"                             |"
            write (*,"(a40, f10.3, a30)") "|                                       ", this%stepsize(2),"                             |"
            write (*,"(a40, f10.3, a30)") "|                                       ", this%stepsize(3),"                             |"
            write (*,"(a40, f10.1, a30)") "|               Minimum test step size: ", this%min_step, "                             |"
            write (*,"(a40, f10.1, a30)") "|               Maximum test step size: ", this%max_step, "                             |"
            write (*,"(a40, f10.1, a30)") "|               Minimum BFGS step size: ", this%min_step_BFGS, "                             |"
            write (*,"(a40, f10.1, a30)") "|               Maximum BFGS step size: ", this%max_step_BFGS, "                             |"
            write (*,"(a40, f10.1, a30)") "|                  Minimal vp velocity: ", this%min_vp,  "                             |"
            write (*,"(a40, f10.1, a30)") "|                  Maximal vp velocity: ", this%max_vp,  "                             |"
            write (*,"(a40, f10.1, a30)") "|                  Minimal vs velocity: ", this%min_vs,  "                             |"
            write (*,"(a40, f10.1, a30)") "|                  Maximal vs velocity: ", this%max_vs,  "                             |"
        endif
        close(19)
    end subroutine readInvfile

    subroutine LowPassZeroPhaseFilter(stf,nt,deltat,upperfreq, tapern, frac)
    ! Lowpass filter a timeseries to a given cutoff frequency
        type(time_Series) :: timeseries, temp                       ! Use time_series type from timeSeriesMod
        real(kind=CUSTOM_REAL), dimension(:) :: stf                 ! source time function to be filtered
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: stf_temp   ! temp array for stf
        real(kind=CUSTOM_REAL) :: deltat, upperfreq                 ! deltat: time sampling, upperfreq: cutoff frequency
        integer :: nt,j, frac                                       ! nt: length of stf, frac: fraction of stf to be tapered on both ends
        logical :: tapern                                           ! flag for tappering

        allocate(stf_temp(size(stf)))

        call makeTimeSeries(timeseries,nt,deltat,stf)               ! Create TimeSeries to apply following functions
        if (tapern) timeseries=hanningTaperTimeSeries(timeseries,real(nt/frac*deltat))  ! Tapering
        timeseries=lowPassButterworthRecursiveTimeSeries(timeseries,real(upperfreq),4)  ! Low-pass filter of degree 4
        temp=ReverseTimeSeries(timeseries)                                              ! Reverse time series to get zero-time-shift filter
        temp=lowPassButterworthRecursiveTimeSeries(temp,real(upperfreq),4)              ! Filter again
        timeseries=ReverseTimeSeries(temp)                                              ! Reverse time series to old order
        do j=1, nt                                                                      ! Extract information from timeSeries type to old stf array
            stf_temp(j)=getSampleTimeSeries(timeseries,j)
        enddo
        stf=stf_temp
        deallocate(stf_temp)
    end subroutine LowPassZeroPhaseFilter

    subroutine cubicInterpolationArray(t_old, y_old, t_new, y_new)
        ! do an interpolation based on four points and a cubic polynomial. Finds the interpolated value y_new at times t_new.
        ! t_old and t_new need to be ordered at equally spaced!
        real(kind=8), dimension(:) :: t_old, t_new              ! time samples
        real(kind=8), dimension(:) :: y_old, y_new              ! values at time samples
        real(kind=8) ::a,b,c,d,y12,t23,t34,t24,y23,t12,t13,y34,t1t1,t1t2,t2t3,t3t3,t2t2,t3t4,t4t4, t1pt2    ! interpolation variables
        integer :: i, ind, length_old, length_new               ! counter
        length_old=size(t_old)                                  ! just get the size of old and new array
        length_new=size(t_new)
        do i=1,length_new                                       ! for all values of the new array
            if (i <= length_old) then                           ! Need this if for next if
                if (abs(t_new(i)-t_old(i)) < 1.E-6) then        ! If both time values are extremly close to each other, just copy the y value
                    y_new(i)=y_old(i)
                else
                    ! Find the position of the t_old value closest but smaller as currrently used t_new value
                    ! Found index is always the second of four points for the cubic interpolation
                    call findindex(real(t_old,kind=CUSTOM_REAL),real(t_new(i),kind=CUSTOM_REAL),ind)
                    if (ind==1) then    ! Reset index if the interpolation is close to the ends of the old array
                        ind=2
                    else if (length_old-ind==1) then
                        ind=length_old-2
                    else if (length_old-ind==0) then
                        ind=length_old-2
                    endif
                    ! Get parameters for interpolation
                    y12=y_old(ind-1)-y_old(ind)
                    t23=t_old(ind)-t_old(ind+1)
                    t34=t_old(ind+1)-t_old(ind+2)
                    t24=t_old(ind)-t_old(ind+2)
                    y23=y_old(ind)-y_old(ind+1)
                    t12=t_old(ind-1)-t_old(ind)
                    t13=t_old(ind-1)-t_old(ind+1)
                    y34=y_old(ind+1)-y_old(ind+2)
                    t1t1=t_old(ind-1)*t_old(ind-1)
                    t1t2=t_old(ind-1)*t_old(ind)
                    t2t3=t_old(ind)*t_old(ind+1)
                    t3t3=t_old(ind+1)*t_old(ind+1)
                    t2t2=t_old(ind)*t_old(ind)
                    t3t4=t_old(ind+1)*t_old(ind+2)
                    t4t4=t_old(ind+2)*t_old(ind+2)
                    t1pt2=t_old(ind-1)+t_old(ind)
                    ! get parameters of cubic polynomial
                    a=(y12*t23*t34*t24-y23*t12*t34*(t24+t13)+y34*t23*t13*t12)/(t12*t23*t34)/&
                        ((t1t1+t1t2-t2t3-t3t3)*t24-(t2t2+t2t3-t3t4-t4t4)*t13)
                    b=(y23*t34-y34*t23)/(t23*t34*t24)-a*(t2t2+t2t3-t3t4-t4t4)/t24
                    c=y12/t12-a*(t1t1+t1t2+t2t2)-b*t1pt2
                    d=y_old(ind)-a*t_old(ind)**3-b*t_old(ind)**2-c*t_old(ind)
                    ! Get new value
                    y_new(i)=a*t_new(i)**3+b*t_new(i)**2+c*t_new(i)+d
                endif
            else
                ! Similar as above
                call findindex(real(t_old,kind=CUSTOM_REAL),real(t_new(i),kind=CUSTOM_REAL),ind)
                if (ind==1) then
                    ind=2
                else if (length_old-ind==1) then
                    ind=length_old-2
                else if (length_old-ind==0) then
                    ind=length_old-2
                endif
                y12=y_old(ind-1)-y_old(ind)
                t23=t_old(ind)-t_old(ind+1)
                t34=t_old(ind+1)-t_old(ind+2)
                t24=t_old(ind)-t_old(ind+2)
                y23=y_old(ind)-y_old(ind+1)
                t12=t_old(ind-1)-t_old(ind)
                t13=t_old(ind-1)-t_old(ind+1)
                y34=y_old(ind+1)-y_old(ind+2)
                t1t1=t_old(ind-1)*t_old(ind-1)
                t1t2=t_old(ind-1)*t_old(ind)
                t2t3=t_old(ind)*t_old(ind+1)
                t3t3=t_old(ind+1)*t_old(ind+1)
                t2t2=t_old(ind)*t_old(ind)
                t3t4=t_old(ind+1)*t_old(ind+2)
                t4t4=t_old(ind+2)*t_old(ind+2)
                t1pt2=t_old(ind-1)+t_old(ind)
                a=(y12*t23*t34*t24-y23*t12*t34*(t24+t13)+y34*t23*t13*t12)/(t12*t23*t34)/&
                    ((t1t1+t1t2-t2t3-t3t3)*t24-(t2t2+t2t3-t3t4-t4t4)*t13)
                b=(y23*t34-y34*t23)/(t23*t34*t24)-a*(t2t2+t2t3-t3t4-t4t4)/t24
                c=y12/t12-a*(t1t1+t1t2+t2t2)-b*t1pt2
                d=y_old(ind)-a*t_old(ind)**3-b*t_old(ind)**2-c*t_old(ind)
                y_new(i)=a*t_new(i)**3+b*t_new(i)**2+c*t_new(i)+d
            endif
        enddo
    end subroutine cubicInterpolationArray

    subroutine cubicInterpolationSnglVal(t_old, y_old, t_new, y_new)
    ! Does the same as cubicInterpolationArray for just one value of t_new
        real(kind=8), dimension(:) :: t_old, y_old
        real(kind=8) :: t_new
        real(kind=CUSTOM_REAL) :: y_new
        real(kind=8) ::a,b,c,d,y12,t23,t34,t24,y23,t12,t13,y34,t1t1,t1t2,t2t3,t3t3,t2t2,t3t4,t4t4, t1pt2
        integer :: ind, length_old
        length_old=size(t_old)

        call findindex(real(t_old,kind=CUSTOM_REAL),real(t_new,kind=CUSTOM_REAL),ind)
        if (ind==1) then
            ind=2
        else if (length_old-ind==1) then
            ind=length_old-2
        else if (length_old-ind==0) then
            ind=length_old-2
        endif
        y12=y_old(ind-1)-y_old(ind)
        t23=t_old(ind)-t_old(ind+1)
        t34=t_old(ind+1)-t_old(ind+2)
        t24=t_old(ind)-t_old(ind+2)
        y23=y_old(ind)-y_old(ind+1)
        t12=t_old(ind-1)-t_old(ind)
        t13=t_old(ind-1)-t_old(ind+1)
        y34=y_old(ind+1)-y_old(ind+2)
        t1t1=t_old(ind-1)*t_old(ind-1)
        t1t2=t_old(ind-1)*t_old(ind)
        t2t3=t_old(ind)*t_old(ind+1)
        t3t3=t_old(ind+1)*t_old(ind+1)
        t2t2=t_old(ind)*t_old(ind)
        t3t4=t_old(ind+1)*t_old(ind+2)
        t4t4=t_old(ind+2)*t_old(ind+2)
        t1pt2=t_old(ind-1)+t_old(ind)
        a=(y12*t23*t34*t24-y23*t12*t34*(t24+t13)+y34*t23*t13*t12)/(t12*t23*t34)/&
            ((t1t1+t1t2-t2t3-t3t3)*t24-(t2t2+t2t3-t3t4-t4t4)*t13)
        b=(y23*t34-y34*t23)/(t23*t34*t24)-a*(t2t2+t2t3-t3t4-t4t4)/t24
        c=y12/t12-a*(t1t1+t1t2+t2t2)-b*t1pt2
        d=y_old(ind)-a*t_old(ind)**3-b*t_old(ind)**2-c*t_old(ind)
        y_new=real(a*t_new**3+b*t_new**2+c*t_new+d, kind=CUSTOM_REAL)
    end subroutine cubicInterpolationSnglVal

    subroutine findindex(t_old,t_new,ind)
    ! finds the index (ind) of a value t_new within the array t_old. The found index corresponds to the closest but smaller value in t_old compared to t_new
        real(kind=CUSTOM_REAL), dimension(:) :: t_old
        real(kind=CUSTOM_REAL) :: t_new
        integer :: ind, left, right, middle, size_t

        ind=-1
        left=1
        size_t=size(t_old)
        right=size_t
        middle=(right+left)/2
        if (t_new >= t_old(right)) then     ! check if t_new is larger than the largest t_old value
            left=right
            ind=right   ! set index to largest value of t_old
            ! if t_new is only one time step larger as the largest t_old value it can still be used
            if (t_new > 1.0001*(2*t_old(right)-t_old(right-1))) then ! factor 1.0001 to handle bugs caused by numeric accuracy of numbers
                write(*,*) 't_new value outside of t_old range', t_new, t_old(right), 2*t_old(right)-t_old(right-1)
                call stop_mpi()
            endif
        endif
        if (t_new<t_old(1)) then    ! check if t_new value is too small
            write(*,*) 't_new value smaller than smallest value of t_old array', t_new, t_old(1)
            call stop_mpi()
        endif
        do while (right-left > 1)   ! divide search array in halfs until value is found
            middle=(right+left)/2

            if (t_old(middle) > t_new) then
                right=middle
            else if (t_old(middle) < t_new) then
                left=middle
            else
                ind=middle
                exit
            endif
        enddo
        if (ind/=middle) ind=left   ! set ind to the found index value
    end subroutine findindex


    subroutine makeTimeSeriesDP(timeseries, nt, deltat, stf)
    ! Interface function to create a time series from timeSeriesMod for double precision
        type(Time_Series) :: timeseries
        real(kind=8) :: deltat
        real(kind=8), dimension(:) :: stf
        integer :: nt

        call createDPFromDataTimeSeries(timeseries,nt,dble(0.),real(deltat),real(stf))
    end subroutine

    subroutine makeTimeSeriesSP(timeseries, nt, deltat, stf)
    ! Interface function to create a time series from timeSeriesMod for single precision
        type(Time_Series) :: timeseries
        real(kind=4) :: deltat
        real(kind=4), dimension(:) :: stf
        integer :: nt

        call createSPFromDataTimeSeries(timeseries,nt,0.,deltat,stf)
    end subroutine

!   Some functions in experimental state for smoothing
!
!    subroutine MPI_Send_Elem(arr_send, arr_rec, mesh, myrank)
!        real(kind=CUSTOM_REAL), dimension(:) :: arr_send
!        real(kind=CUSTOM_REAL), dimension(:,:) :: arr_rec
!        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: send
!        type(MeshVar) :: mesh
!        integer :: i,j, dest, tag, req, reqrec, myrank
!
!        allocate(send(mesh%mpi_ne,mesh%mpi_nn))
!        do i=1,mesh%mpi_nn
!            do j=1,mesh%mpi_ne ! loop over interface elements
!                if ( mesh%mpi_connection(i,j,1) >0) then
!                    send(j,i)= arr_send(mesh%mpi_connection(i,j,1))   ! sende Werte
!                end if
!            end do
!        end do
!        tag=0
!        do i=1,mesh%mpi_nn
!            dest=mesh%mpi_neighbor(i)-1
!            call isendV_real(send(:,i),mesh%mpi_ne,dest,tag,req,CUSTOM_REAL)
!            call irecV_real(arr_rec(:,i),mesh%mpi_ne,dest,tag,reqrec,CUSTOM_REAL)
!            call wait_req(req)
!            call wait_req(reqrec)
!        end do
!        deallocate(send)
!    end subroutine MPI_Send_Elem
!
!    subroutine get_modelmisfit(myrank, mesh, inv, rec_vp, rec_vs, chi_model)
!        integer :: myrank
!        type(meshVar) :: mesh
!        type(InvVar) :: inv
!        integer :: i, j
!        real(kind=CUSTOM_REAL), dimension(:,:) :: rec_vp, rec_vs
!        real(kind=CUSTOM_REAL) :: chi_model_sngl, chi_model, u
!
!        chi_model_sngl=0.
!        do i=1,mesh%nelem
!            u=mesh%vp(i)*mesh%smooth_A(i,1)
!            if (myrank==0) write(*,*) i, 1, u, mesh%vp(i), mesh%smooth_A(i,:)
!            do j=2,4
!                if (mesh%smooth_A(i,j)>0) then
!                    u = u + mesh%vp(mesh%smooth_A(i,j))
!                    if (myrank==0) write(*,*) i, j, u, mesh%vp(mesh%smooth_A(i,j)), mesh%smooth_A(i,j), 'normal'
!                endif
!                if (mesh%smooth_A(i,j)<0) then
!                    u = u + rec_vp(-mesh%smooth_A(i,j),mesh%mpi_interface(4,j-1,i))
!                    if (myrank==0) write(*,*) i, j, u, rec_vp(-mesh%smooth_A(i,j),mesh%mpi_interface(4,j-1,i)), mesh%smooth_A(i,j), 'MPI'
!                endif
!            enddo
!            chi_model_sngl=chi_model_sngl + u**2
!            u=mesh%vs(i)*mesh%smooth_A(i,1)
!            do j=2,4
!                if (mesh%smooth_A(i,j)>0) u = u + mesh%vs(mesh%smooth_A(i,j))
!                if (mesh%smooth_A(i,j)<0) u = u + rec_vs(-mesh%smooth_A(i,j),mesh%mpi_interface(4,j-1,i))
!            enddo
!            chi_model_sngl=chi_model_sngl + u**2
!        enddo
!        call sum_real_all(chi_model_sngl,chi_model,CUSTOM_REAL)
!    end subroutine get_modelmisfit
!
!
!    subroutine get_modelsmoothing(myrank, mesh, inv, rec_vp, rec_vs, model_cp, model_cs)
!        integer :: myrank
!        type(meshVar) :: mesh
!        type(InvVar) :: inv
!        real(kind=CUSTOM_REAL), dimension(:,:) :: rec_vp, rec_vs
!        real(kind=CUSTOM_REAL), dimension(:) :: model_cp, model_cs
!        real(kind=CUSTOM_REAL), dimension(:), allocatable :: u
!        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: rec_u
!        integer :: i,j
!
!        allocate(u(mesh%nelem))
!        allocate(rec_u(mesh%mpi_ne,mesh%mpi_nn))
!
!        do i=1,mesh%nelem   ! u=Am
!            u(i)=mesh%vp(i)*mesh%smooth_A(i,1)
!            do j=2,4
!                if (mesh%smooth_A(i,j)>0) u(i) = u(i) + mesh%vp(mesh%smooth_A(i,j))
!                if (mesh%smooth_A(i,j)<0) u(i) = u(i) + rec_vp(-mesh%smooth_A(i,j),mesh%mpi_interface(4,j-1,i))
!            enddo
!        enddo
!        call MPI_Send_Elem(u, rec_u, mesh, myrank)  ! MPI for u
!        do i=1,mesh%nelem   ! cp=Au
!            model_cp(i)=u(i)*mesh%smooth_A(i,1)
!            do j=2,4
!                if (mesh%smooth_A(i,j)>0) model_cp(i) = model_cp(i) + u(mesh%smooth_A(i,j))
!                if (mesh%smooth_A(i,j)<0) model_cp(i) = model_cp(i) + rec_u(-mesh%smooth_A(i,j),mesh%mpi_interface(4,j-1,i))
!            enddo
!        enddo
!        !same for vs
!        do i=1,mesh%nelem   ! u=Am
!            u(i)=mesh%vs(i)*mesh%smooth_A(i,1)
!            do j=2,4
!                if (mesh%smooth_A(i,j)>0) u(i) = u(i) + mesh%vs(mesh%smooth_A(i,j))
!                if (mesh%smooth_A(i,j)<0) u(i) = u(i) + rec_vs(-mesh%smooth_A(i,j),mesh%mpi_interface(4,j-1,i))
!            enddo
!        enddo
!        call MPI_Send_Elem(u, rec_u, mesh, myrank)  ! MPI for u
!        do i=1,mesh%nelem   ! cs=Au
!            model_cs(i)=u(i)*mesh%smooth_A(i,1)
!            do j=2,4
!                if (mesh%smooth_A(i,j)>0) model_cs(i) = model_cs(i) + u(mesh%smooth_A(i,j))
!                if (mesh%smooth_A(i,j)<0) model_cs(i) = model_cs(i) + rec_u(-mesh%smooth_A(i,j),mesh%mpi_interface(4,j-1,i))
!            enddo
!        enddo
!        deallocate(u, rec_u)
!    end subroutine get_modelsmoothing

end module adjointMod
