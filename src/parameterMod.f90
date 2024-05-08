module parameterMod

    ! module to read a parameter file to set up the simulation
    use constantsMod
    use fileParameterMod
    use errorMessage

    implicit none

    type :: parameterVar
        !private
        character(len=80) :: title              !Title of the simulation
        logical :: log = .true.                 !Write information onto screen

        !external model
        character(len=255) :: externalfilename  !File containing the external model
        !character(len=255) :: extmatpropfilename !File containing the external matprop
        logical :: extvel                       !read external model?
        !logical :: extmatprop                   !read external matprop file?

        integer :: nproc                        !Number of processors for the calculation

        !Materials
        integer :: matn                         !Number of different materials

        !Timeintegration
        integer :: timeint                      !Type of timeintegration (1:euler 2:rk2 3:rk3)
        integer :: nt                           !Number of timesteps
        logical :: autont                       !Calculate number of timesteps automatically?
        logical :: autodt                       !Calculate timestep automatically?
        real(kind=CUSTOM_REAL) :: dt            !timestep (needs to be set if autodt = .false.)
        real(kind=CUSTOM_REAL) :: t_total       !total simulated time
        real(kind=CUSTOM_REAL) :: cfl           !Courantnumber (cfl value for dt)
        real(kind=CUSTOM_REAL) :: simt0         !Starting time of simulation

        !Absorbing boundaries / PML parameters
        logical :: set_pml                      !pml
        real(kind=CUSTOM_REAL) :: pml_delta !pml
        real(kind=CUSTOM_REAL) :: pml_rc, pml_kmax, pml_afac !pml
        logical :: use_trigger                  !use sta_lta trigger for energy monitoring
        integer :: avg_window1, avg_window2     !lta and sta windows
        real(kind=CUSTOM_REAL) :: sta_lta_trigger !threshold for trigger

        !Sources
        logical :: shift_sources

        !receiver
        logical :: global_rec_angle             !use globally or locally defined rotation angles
        real(kind=CUSTOM_REAL) :: rec_angle     !rotate receivers
        integer :: subsampling_factor           !subsampling of seismograms by this factor

        !Seismograms
        logical :: div                          ! Enables the seperate calculation of the radial component of the seismogram
        logical :: curl                         ! Enables the seperate calculation of the tangential component of the seismogram
        logical :: autoshift                    ! Enables the automatic shift of the seismogram by 1.2/f0
        real(kind=CUSTOM_REAL) :: plott0        ! Offset for the timeaxis for the seismogram

    end type parameterVar

    type movie_parameter
        !A series of switches to select which files are created for plotting
        integer :: frame        !Number of timesteps for each frame of the movie
        logical :: movie        !Select if any moviefiles are created
        logical :: trimesh      !Create files with average in each element
        logical :: points       !Create files with data for each point
        logical :: displacement !Plot displacement field
        logical :: velocity     !Plot velocity field
        logical :: stress       !Plot stress field
    end type

    contains

    subroutine readParfile(this, movie, myrank, errmsg)
        implicit none
        type(error_message) :: errmsg
        type (parameterVar) :: this
        !type(lsi_parameter) :: lsipar
        type(movie_parameter) :: movie
        integer, intent(in) :: myrank

        !local variables
        character(len=80) :: filename
        character(len=11) :: myname = "readParfile"

        integer :: ier
        logical :: log = .true.

        call addTrace(errmsg,myname)

        filename=trim('data/parfile')
        open(unit=19, file=trim(filename), status='old', action='read', iostat=ier)
        if (ier /= 0) then
            call add(errmsg,2,'Could not open file.' ,myname, filename)
            call print(errmsg)
            stop
        endif

        ! Cycle to read the parameters form the parameter file. The order of appearence of the parameters is not important
        call readStringPar(this%title, "title", filename, 0, errmsg)
        call readStringPar(this%externalfilename, "externalfilename", filename, 0, errmsg)
        call readIntPar(this%nproc, "nproc", filename, 0, errmsg)
        call readIntPar(this%timeint, "timeint", filename, 0, errmsg)
        call readIntPar(movie%frame, "frame", filename, 0, errmsg)
        call readIntPar(this%avg_window1, "avg_window1", filename, 0, errmsg)
        call readIntPar(this%avg_window2, "avg_window2", filename, 0, errmsg)
        call readIntPar(this%subsampling_factor, "subsampling_factor", filename, 0, errmsg)
        call readFloatPar(this%cfl, "cfl", filename, 0, errmsg)
        call readFloatPar(this%pml_delta, "pml_delta", filename, 0, errmsg)
        call readFloatPar(this%pml_rc,"pml_rc", filename, 0, errmsg)
        call readFloatPar(this%pml_kmax, "pml_kmax", filename, 0, errmsg)
        call readFloatPar(this%pml_afac,"pml_afac", filename, 0, errmsg)
        call readFloatPar(this%sta_lta_trigger,"sta_lta_trigger", filename, 0, errmsg)
        call readFloatPar(this%dt,"dt", filename, 0, errmsg)
        call readFloatPar(this%simt0,"simt0", filename, 0, errmsg)
        call readLogicalPar(this%log, "log", filename, 0, errmsg)
        call readLogicalPar(this%extvel, "extvel", filename, 0, errmsg)
        call readLogicalPar(this%set_pml, "set_pml", filename, 0, errmsg)
        call readLogicalPar(this%global_rec_angle, "global_rec_angle", filename, 0, errmsg)
        if (this%global_rec_angle) then
            call readFloatPar(this%rec_angle,"rec_angle", filename, 0, errmsg)
        end if
        call readLogicalPar(this%autont, "autont", filename, 0, errmsg)
        if (this%autont) then
            call readFloatPar(this%t_total,"t_total", filename, 0, errmsg)
        else
            call readIntPar(this%nt, "nt", filename, 0, errmsg)
        end if
        call readLogicalPar(this%autodt, "autodt", filename, 0, errmsg)
        call readLogicalPar(this%use_trigger, "use_trigger", filename, 0, errmsg)
        call readLogicalPar(movie%movie, "movie", filename, 0, errmsg)
        call readLogicalPar(movie%displacement, "save_movie_displacement", filename, 0, errmsg)
        call readLogicalPar(movie%velocity, "save_movie_velocity", filename, 0, errmsg)
        call readLogicalPar(movie%stress, "save_movie_stress", filename, 0, errmsg)
        call readLogicalPar(movie%trimesh, "save_movie_trimesh", filename, 0, errmsg)
        call readLogicalPar(movie%points, "save_movie_points", filename, 0, errmsg)
        call readLogicalPar(this%div, "div", filename, 0, errmsg)
        call readLogicalPar(this%curl, "curl", filename, 0, errmsg)
        call readLogicalPar(this%autoshift, "autoshift", filename, 0, errmsg)
        if (.not. this%autoshift) then
            call readFloatPar(this%plott0,"plott0", filename, 0, errmsg)
        end if
        call readLogicalPar(this%shift_sources, "shift_sources", filename, 0, errmsg)
        close(19)

        !test if certain parameters have been read/entered correctly
		call ErrorMessages(this, movie, myname, errmsg, filename)

        log = this%log

        if (log .and. myrank == 0) then
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, a10, a30)")   "|              Title of the simulation: ", this%title, "                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                               Global parameters                              |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, i3,  a37)")   "|           Precision for real numbers: ", CUSTOM_REAL*8, " bits                               |"
            write (*,"(a40, i2,  a38)")   "|                                Order: ", ORDER, "                                     |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                               Basic parameters                               |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, l10, a30)")   "|                           Create log: ", this%log, "                             |"
            write (*,"(a40, i10, a30)")   "|                 Number of processors: ", this%nproc, "                             |"
            write (*,"(a40, l10, a30)")   "|        Load external velocitiy model: ", this%extvel, "                             |"
            if (this%extvel) write (*,*) trim(this%externalfilename)
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                               Movie parameters                               |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, l10, a30)")   "|                           Save movie: ", movie%movie, "                             |"
            if (movie%movie) then
                write (*,"(a40, i10, a30)")   "|   Number of timesteps for movieframe: ", movie%frame, "                             |"
                write (*,"(a40, l10, a30)")   "|           Create TriMesh-movie-files: ", movie%trimesh, "                             |"
                write (*,"(a40, l10, a30)")   "|            Create Points-movie-files: ", movie%points, "                             |"
                write (*,"(a40, l10, a30)")   "|                  Plot velocity field: ", movie%velocity, "                             |"
                write (*,"(a40, l10, a30)")   "|              Plot displacement field: ", movie%displacement, "                             |"
                write (*,"(a40, l10, a30)")   "|                    Plot stress field: ", movie%stress, "                             |"
            end if
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                         Parameters for timeintegration                       |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, i10, a30)")   "|                      Timeintegartion: ", this%timeint, "                             |"
            write (*,"(a40, l10, a30)")   "| Calculate number of timesteps autom.: ", this%autont, "                             |"
            if (this%autont) then
                write (*,"(a40, es10.3, a30)") "|                 Total simulated time: ", this%t_total, "                             |"
            else
                write (*,"(a40, i10, a30)")   "|                  Number of timesteps: ", this%nt, "                             |"
            end if
            write (*,"(a40, l10, a30)")   "|                               autodt: ", this%autodt, "                             |"
            if (.not. this%autodt) then
                write (*,"(a40, e10.3, a30)") "|                                   dt: ", this%dt, "                             |"
            endif
            write (*,"(a40, f10.1, a30)") "|                            cfl value: ", this%cfl, "                             |"
            write (*,"(a40, f10.1, a30)") "|          Starting time of simulation: ", this%simt0, "                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                                 PML parameters                               |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, l10, a30)")   "|                          pml enabled: ", this%set_pml, "                             |"
            write (*,"(a40, f10.1, a30)") &
              "|                        pml thickness: ", this%pml_delta, "                             |"
            write (*,"(a40, f10.3, a30)") "|           pml reflection coefficient: ", this%pml_rc, "                             |"
            write (*,"(a40, f10.1, a30)") &
              "|                             pml kmax: ", this%pml_kmax, "                             |"
            write (*,"(a40, f10.1, a30)") &
              "|                      factor for amax: ", this%pml_afac, "                             |"
            write (*,"(a40, i10, a30)")   "|                           lta window: ", this%avg_window1, "                             |"
            write (*,"(a40, i10, a30)")   "|                           sta window: ", this%avg_window2, "                             |"
            write (*,"(a40, f10.1, a30)") &
              "|                    sta_lta threshold: ", this%sta_lta_trigger, "                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                          Parameters for receivers                            |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, l10, a30)") "|                     global_rec_angle: ", this%global_rec_angle, "                             |"
            if (this%global_rec_angle) then
                write (*,"(a40, f10.1, a30)") "|                            rec_angle: ", this%rec_angle, "                             |"
            end if
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                           Parameters for sources                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, l10, a30)") "|                        shift_sources: ", this%shift_sources, "                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                       Parameters regarding seismograms                       |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, i10, a30)")   "|                   subsampling_factor: ", this%subsampling_factor, "                             |"
            write (*,"(a40, l10, a30)")   "|  Auto-shift the seismogram by 1.2/f0: ", this%autoshift, "                             |"
            if (.not. this%autoshift) then
                write (*,"(a40, f10.7, a30)") "|                               Offset: ", this%plott0, "                             |"
            end if
            write (*,"(a40, l10, a30)")   "|                 Calculate divergence: ", this%div, "                             |"
            write (*,"(a40, l10, a30)")   "|                       Calculate curl: ", this%curl, "                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
        end if
    end subroutine readParfile

    subroutine ErrorMessages(par, movie, myname, errmsg, filename)
        ! This subroutine produces error messages if certain parameters are set in the wrong way.

        !input:
        type(parameterVar) :: par
        type(error_message) :: errmsg
        type(movie_parameter) :: movie
        character(len=*) :: myname, filename

        if (movie%movie) then
            if (movie%frame < 0) call add(errmsg, 2, "Number of time steps per movie frame needs to be positive.", myname, filename)
            if (.not. par%autont) then
                if (movie%frame > par%nt) call add(errmsg, 2, "Number of time steps per movie frame less than the total number of time steps.", myname, filename)
            end if
            if (movie%frame == 0) call add(errmsg, 1, "If the number of time steps per movie frame is 0 no movie-files are created.", myname, filename)
        end if
        if (par%nproc <= 0) call add(errmsg, 2, "Number of cores needs to be positive.", myname, filename)
        if (par%nproc == 1) call add(errmsg, 1, "Single core calculations are not advised.", myname, filename)
        if (par%cfl <= 0.0) call add(errmsg, 2, "Parameter 'cfl' has to be positive. Recommended range is 0.0 < cfl <= 1.0", myname, filename)
        if (par%cfl > 1.0) call add(errmsg, 1, "Parameter 'cfl' might be too high. Recommended range is 0.0 < cfl <=0 1.0", myname, filename)
        if (.not. par%autodt .and. par%dt < epsilon(par%dt)) call add(errmsg, 2, "Choose dt /= 0, or autodt = .true.!", myname, filename )
        if( par%timeint < 1 .or. par%timeint > 4) call add(errmsg, 2, "Parameter to select time stepping variant is out of range. Select 1, 2, 3 or 4.", myname, filename)
        if (par%subsampling_factor <= 0) call add(errmsg, 2, "Parameter 'subsampling_factor' must be positive.", myname, filename)
    end subroutine
end module parameterMod
