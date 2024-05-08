program mesher
    !program to convert the given mesh files to build a MPI version of the code
    use parameterMod
    use meshMod
    use plotMod
    use errorMessage

    implicit none

    type(error_message) :: errmsg
    type (parameterVar) :: par
    type (meshVar) :: mesh
    type(movie_parameter) :: movie
    character(len=80) :: filename
    integer :: myrank = 0
    character(len=6) :: myname = "mesher"

    !Create new errormessagetrace
    call new(errmsg,myname)

    ! read Parfile
	call readParfile(par, movie, myrank, errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif

    ! log?
    if (par%log) then
        write(*,'(a80)') "|------------------------------------------------------------------------------|"
        write(*,'(a80)') "|                               start mesher                                   |"
        write(*,'(a80)') "|------------------------------------------------------------------------------|"
    end if

	call createRegularMesh(mesh, par, errmsg)

    if (par%log) write(*,'(a80)') "|                                  plot mesh                                   |"
    call triangulation_order3_plot ( trim(outpath)//"mesh.ps", mesh%ncoord, dble(mesh%coord), mesh%nelem, mesh%elem, 2, 2 )
    filename = "pointsmesh.txt"
    call plotPoints2d(mesh%vx,mesh%vz,trim(outpath)//filename)
    if (.level.errmsg == 1) call print(errmsg)
    call deallocMeshvar(mesh)
    
    if (par%log) then
        write(*,'(a80)') "|------------------------------------------------------------------------------|"
        write(*,'(a80)') "|                                end mesher                                    |"
        write(*,'(a80)') "|------------------------------------------------------------------------------|"
    end if
end program mesher
