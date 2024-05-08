program movie_pro
    use parameterMod
    use collectMovieMod
    use errorMessage

    type (parameterVar) :: par
    type (movie_parameter) :: movie
    type (error_message) :: errmsg
    integer :: myrank = 0
    character(len=13) :: myname = "program_movie"

    !Create new error message for this program
    call new(errmsg, myname)

    !call writeLogo()

    ! read Parfile
	call readParfile(par, movie, myrank, errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif

    if (par%log) then
       write (*, "(a80)") "|------------------------------------------------------------------------------|"
       write (*, "(a80)") "|                                start movie                                   |"
       write (*, "(a80)") "|------------------------------------------------------------------------------|"
    end if

    call collectMovie(par, movie, 1, errmsg)
    
    if (par%log) then
        write(*,'(a80)') "|------------------------------------------------------------------------------|"
        write(*,'(a80)') "|                                  end movie                                   |"
        write(*,'(a80)') "|------------------------------------------------------------------------------|"
    end if
end program
