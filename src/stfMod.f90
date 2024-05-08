module stfMod

    use constantsMod
    use fileunitMod
    use errorMessage

    implicit none

    contains

    function stfRicker(t,f0,t0,factor)
        implicit none
        real(kind=custom_real) :: stfRicker
        real(kind=custom_real) :: f0,t0,factor,t
        real(kind=custom_real) :: aval
        aval = pi*pi*f0*f0
        stfRicker =  factor * (1.-2.*aval*(t-t0)**2.) * exp(-aval*(t-t0)**2.)
    end function stfRicker

    function stfDiffRicker(t,f0,t0,factor)
        implicit none
        real(kind=custom_real) :: stfDiffRicker
        real(kind=custom_real) :: f0,t0,factor,t
        real(kind=custom_real) :: aval
        aval = pi*pi*f0*f0
        stfDiffRicker =  2*aval*factor*(t-t0)*exp(-aval*(t-t0)**2.)*(2*aval*(t-t0)**2-3)
    end function stfDiffRicker

    function stfGauss(t,f0,t0,factor)
        implicit none
        real(kind=custom_real) :: stfGauss
        real(kind=custom_real) :: f0,t0,factor,t
        real(kind=custom_real) :: aval
        aval = pi*pi*f0*f0
        stfGauss = -factor * sqrt(pi) * f0 * exp(-aval*(t-t0)**2.)
    end function stfGauss

    function stfDiffGauss(t,f0,t0,factor)
        implicit none
        real(kind=custom_real) :: stfDiffGauss
        real(kind=custom_real) :: f0,t0,factor,t
        real(kind=custom_real) :: aval
        aval = pi*pi*f0*f0
        stfDiffGauss = factor * sqrt(pi) * f0 * 2. * aval * (t-t0) * exp(-aval*(t-t0)**2.)
    end function stfDiffGauss

    function stfSin3(t,f0,t0,factor)
        implicit none
        real(kind=custom_real) :: stfSin3
        real(kind=custom_real) :: f0,t0,factor,t
        if ((0 < t-t0) .and. (t-t0 < (1.0/f0))) then
            stfSin3 = factor * sin(pi*(t-t0)*f0)**3
        else
            stfSin3 = 0.0
        end if
    end function stfSin3

    function stfDiffSin3(t,f0,t0,factor)
        implicit none
        real(kind=custom_real) :: stfDiffSin3
        real(kind=custom_real) :: f0,t0,factor,t
        if ((0 < t-t0) .and. (t-t0 < (1.0/f0))) then
            stfDiffSin3 = factor * 3 * sin(pi*(t-t0)*f0)**2 * cos(pi*(t-t0)*f0)**2
        else
            stfDiffSin3 = 0.0
        end if
    end function stfDiffSin3

    subroutine stfExternal(dnew,tnew,dt,nt,filter,w,nl,filename,differentiate, errmsg)
        implicit none
        type(error_message) :: errmsg
        real(kind=CUSTOM_REAL) :: dt
        integer :: nt,w,nl
        integer :: differentiate,n
        character(len=*) :: filename
        logical :: filter,file_exists
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: data,time
        real(kind=CUSTOM_REAL), dimension(nt) :: dnew,tnew,dtmp
        real(kind=CUSTOM_REAL) :: rtemp, dt_old,t0,x0
        integer :: io
        integer :: ns,i,j,x0_f
        integer :: fu

        character(len=11) :: myname = "stfExternal"

        call addTrace(errmsg, myname)

        !open sourcefile
        inquire(file=trim(filename), exist=file_exists)
        if (.not.file_exists) then
            call add(errmsg,2,'File '//trim(filename)//' does not exist!',myname)
            call print(errmsg)
            stop
        end if

        ns = 0
        io=0

        fu = getFileUnit(97)

        open(unit=fu,file=trim(filename),action = 'read')
        do while (io >= 0)
            read(fu,*,iostat=io) rtemp,rtemp
            ns = ns+1
        end do
        close(fu)
        ns = ns-1


        allocate(data(ns),time(ns))
        fu = getFileUnit(97)

        open(unit=fu,file=trim(filename),action = 'read')
        do i=1,ns
            read(fu,*) time(i),data(i)
        end do
        close(fu)

        dt_old = time(2)-time(1)
        t0 = time(1)

        dnew(:) = 0.0
        tnew(:) = 0.0
        do i=1,nt
            tnew(i)=(i-1.)*dt
            x0 = ((i-1.)*dt) / dt_old +1
            x0_f = int(x0)
            if (x0_f<=nl) then
                do j=1,i+nl
                    dnew(i) = dnew(i)+data(j)*lanczos((x0-j),nl)
                end do
            else if (x0 > (ns - nl)) then
                dnew(i) = data(ns)
            else
                do j=(x0_f-nl),x0_f+nl
                    dnew(i) = dnew(i)+data(j)*lanczos((x0-j),nl)
                end do
            end if
        end do

        if (filter) then
            dnew=movFilter(dnew,w,nt)
        end if

        select case (differentiate)
            case (0) !no differentiation
                continue
            case (1) !differentiate once
                n = size(dnew)
                dtmp(1) = (dnew(2) - dnew(1))/dt
                do i = 2, (n-1)
                    dtmp(i) = (dnew(i+1) - dnew(i-1))/(2*dt)
                end do
                dtmp(n) = (dnew(n) - dnew(n-1))/dt
                dnew = dtmp
        end select

        tnew=tnew+t0
        deallocate(data,time)
    end subroutine stfExternal

    function lanczos(x,l)
        real(kind=CUSTOM_REAL) :: x, lanczos
        integer :: l

        if (x < epsilon(x)) then
            lanczos=1.0
        else if (abs(x) < l) then
            lanczos = sin(pi*x)/(pi*x) * sin((pi*x)/l)/((pi*x)/l)
        else
            lanczos= 0.0
        end if
    end function lanczos

    function movFilter(data,w,nt)
        integer :: nt,w
        real(kind=CUSTOM_REAL) , dimension(nt) :: data
        real(kind=CUSTOM_REAL) , dimension(nt) :: movFilter
        integer :: i,j,n,m

        real(kind=CUSTOM_REAL) :: temp

        do i=1,nt
            if (i<=w) then
                n=1
                m=i+w
            else if ( i >= nt-w) then
                n=i-w
                m=nt-1
            else
                n=i-w
                m=i+w
            end if
            temp=0.0
            do j=n,m
                temp=temp+data(j)
            end do
            movFilter(i) = temp / ((2*w)+1)
        end do
    end function movFilter
end module stfMod
