program solver
    ! Program to solve the 2D (poro)elastic wave equation with the discontinous galerkin methode in 2D on triangular mesh
    use parameterMod
    use meshMod
    use timeloopMod
    use sourceReceiverMod
    use errorMessage
    use adjointMod ! TBU
    use collectMovieMod ! TBU

    implicit none

    type(error_message) :: errmsg
    type (parameterVar) :: par
    type (meshVar) :: mesh_change
    type (srcVar) :: src
    type(movie_parameter) :: movie
    integer :: world_size
	integer :: myrank
    integer :: t1, t2, count_max, rate
    integer :: n_src_tmp, n_src
	real(kind=CUSTOM_REAL) :: tges
    logical :: file_exists
    logical, dimension(:), allocatable :: act_src
    character(len=80) :: filename
    character(len=6) :: myname = "solver"
    character(len=7) :: inv_key='forward'

    !Create new errormessagetrace
    call new(errmsg,myname)

    ! start MPI
    call init_mpi()
    call comm_size(world_size)
    call comm_rank(myrank)

    ! read Parfile
	call readParfile(par, movie, myrank, errmsg)

    write(filename,"('/meshVar',i6.6)") myrank+1    ! read in mesh in separate variable to change mesh later
    filename=trim(outpath)//trim(filename)
    inquire(file=trim(filename), exist=file_exists)
    if (file_exists) then
        call readMeshVar(mesh_change,filename, errmsg)
    else
        write(*,*) "error in databases, files not existing"
        call add(errmsg, 2, "Error in databases, file "//trim(filename)//" does not exist!", myname, filename)
        call print(errmsg)
        call stop_mpi()
    end if

    if (.level.errmsg == 2) then
        call print(errmsg)
        stop
    endif

        ! get total number of sources
    n_src_tmp=0
    if (mesh_change%has_src) then
        write(filename,"('srcVar',i6.6)") myrank+1
        inquire(file=trim(outpath)//trim(filename), exist=file_exists)
        if (file_exists) then
            call readSrcVar(src,trim(outpath)//filename)
        else
            call add(errmsg, 2, "File does not exist!", myname, filename)
            call print(errmsg)
            call stop_mpi()
        end if
        n_src_tmp=src%nsrc
        call deallocSrcVar(src)
    endif
    call sum_int_all(n_src_tmp, n_src)

    if (par%nproc > 1) call sync_mpi()

    if (par%log .and. myrank == 0) then
        write(*,'(a80)') "|------------------------------------------------------------------------------|"
        write(*,'(a80)') "|                                start solver                                  |"
    end if

    t1 = 0
    t2 = 0
    if (par%log.and.myrank==0) call system_clock(t1, rate, count_max)
    if (par%log.and.myrank==0) tges=0.
    allocate(act_src(n_src))
    act_src=.true.

    call sync_mpi()
    call timeloop2d(par, inv_key, 0, 1, act_src, movie, myrank, errmsg)
    if (par%log .and. myrank == 0)  write(*,'(a80)') "|------------------------------------------------------------------------------|"
    deallocate(act_src)

    if (par%log .and. myrank == 0) then
        write(*,'(a80)') "|                                 end solver                                   |"
        write(*,'(a80)') "|------------------------------------------------------------------------------|"
    end if

    ! stop MPI
    call finalize_mpi()
end program solver
