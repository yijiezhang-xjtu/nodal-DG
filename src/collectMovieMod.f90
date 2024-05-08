module collectMovieMod
    use parameterMod
    use meshMod
    use vtkMod
    use errorMessage
    use progressbarMod
    use plotMod
    use constantsMod

    implicit none

contains

    subroutine collectMovie(par, movie, inv_step, errmsg)
        !input
        type(parameterVar) :: par
        type(movie_parameter) :: movie
        type(error_message) :: errmsg
        !local
        type(meshVar),dimension(:), allocatable :: db
        type(srcVar) :: src
        integer :: iproc, inv_step, nsrc, isrc
        integer :: nframe,c,d,e,i
        integer :: n_src
        integer, dimension(:), allocatable :: ndata,ndata2,ndata3
        logical :: make_progress
        ! file
        logical :: file_exists
        character(len=80) ::filename
        character(len=7) ::  srcstring
        character(len=15) :: logstring
        ! data
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: uplot
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: ux
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: uz
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: vplot
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: vel_x
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: vel_z
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: sigmaxx
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: sigmazz
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: sigmaxz
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: pressure
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: xyzplot
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: vpplot
        !real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: tempdata
        integer :: glob_nelem, glob_nglob, glob_ncoord
        real(kind=CUSTOM_REAL) , dimension(:,:), allocatable :: glob_coord, glob_coord2
        integer , dimension(:,:), allocatable :: glob_elem

        !For the errormessage:
        character(len=12) :: myname = "collectMovie"

        !Add module to the trace of the errormessage
        call addTrace(errmsg, myname)

        if (par%log) then
            write (*, "(a80)") "|                                generate movie                                |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
        end if
        ! load database

        nframe=par%nt/movie%frame

        allocate(db(par%nproc))
        allocate(ndata(par%nproc))
        allocate(ndata2(par%nproc))
        allocate(ndata3(par%nproc))

        do iproc=1,par%nproc
            write(filename,"('meshVar',i6.6)") iproc
            filename=trim(outpath)//trim(filename)
            inquire(file=trim(filename), exist=file_exists)
            if (file_exists) then
                call readMeshVar(db(iproc),filename, errmsg)

            else
                call add(errmsg, 2, "Error in the database. File "//trim(filename)//" does not exist.", myname, filename)
                if (.level.errmsg == 2) then
                    call print(errmsg)
                    stop
                endif
            end if
        end do

        ! get pointers
        glob_nelem=0
        glob_nglob=0
        glob_ncoord=0
        do iproc=1,par%nproc
            glob_nelem=glob_nelem+db(iproc)%nelem
            glob_nglob=glob_nglob+db(iproc)%nglob
            glob_ncoord=glob_ncoord+db(iproc)%ncoord ! only guess for number of coords
            ndata(iproc) = db(iproc)%nglob
            ndata2(iproc) = db(iproc)%ncoord
            ndata3(iproc) = db(iproc)%nelem
        end do

        allocate(xyzplot(3,glob_nglob))

        ! read coords
        d=1
        e=ndata(1)
        do iproc=1,par%nproc
            xyzplot(1,d:e) = db(iproc)%vx(:)
            xyzplot(2,d:e) = 0.0
            xyzplot(3,d:e) = db(iproc)%vz(:)
            if (iproc < par%nproc) then
                d=d+ndata(iproc)
                e=e+ndata(iproc+1)
            end if
        end do

        !build global mesh
        allocate(glob_elem(3,glob_nelem))
        allocate(glob_coord(3,glob_ncoord))

        !build global elements
        c=1
        do iproc=1,par%nproc
            do i=1,db(iproc)%nelem
                glob_elem(1,c)=db(iproc)%loc2glob_nodes(db(iproc)%elem(1,i))
                glob_elem(2,c)=db(iproc)%loc2glob_nodes(db(iproc)%elem(2,i))
                glob_elem(3,c)=db(iproc)%loc2glob_nodes(db(iproc)%elem(3,i))
                c=c+1
            end do
        end do

        !build global coordinates
        glob_coord=-1e17
        do iproc=1,par%nproc
            do i=1,db(iproc)%ncoord
                glob_coord(1,db(iproc)%loc2glob_nodes(i)) = db(iproc)%coord(1,i)
                glob_coord(2,db(iproc)%loc2glob_nodes(i)) = 0.0
                glob_coord(3,db(iproc)%loc2glob_nodes(i)) = db(iproc)%coord(2,i)
            end do
        end do

        c=0
        do i=1,glob_ncoord
            if (glob_coord(1,i)>-1e17) then
                c=c+1
            end if
        end do
        allocate(glob_coord2(3,c))
        c=0
        do i=1,glob_ncoord
            if (glob_coord(1,i)>-1e17) then
                c=c+1
                glob_coord2(:,c)=glob_coord(:,i)
            end if
        end do

        n_src=0

        do iproc=1, par%nproc
            write(filename,"('srcVar',i6.6)") iproc
            inquire(file=trim(outpath)//trim(filename), exist=file_exists)
            if (file_exists) then
                call readSrcVar(src,trim(outpath)//filename)
                n_src=n_src + src%nsrc
                call deallocSrcVar(src)
            end if
        enddo
        nsrc=1
        make_progress = .true.

        do isrc=1, nsrc
            srcstring=''
            logstring='               '
            if (movie%displacement) then
                if (par%log) write (*,"(a29,a15,a36)") "| Creating Displacement files",logstring,"                                   |"
                allocate(uplot(glob_nglob,nframe))
                allocate(ux(glob_nglob,nframe))
                allocate(uz(glob_nglob,nframe))
                call readplotbin(par%nproc,movie%frame,par%nt,errmsg,uplot,outpath,"moviedata_normU"//srcstring,glob_nglob,nframe,ndata)
                call readplotbin(par%nproc,movie%frame,par%nt,errmsg,ux,outpath,"moviedata_ux"//srcstring,glob_nglob,nframe,ndata)
                call readplotbin(par%nproc,movie%frame,par%nt,errmsg,uz,outpath,"moviedata_uz"//srcstring,glob_nglob,nframe,ndata)
                if (movie%trimesh) then
                    if (par%log) write (*,"(a37,a15,a28)") "| Creating trimesh displacement files",logstring,"                           |"
                    call generateTriMeshFiles(db, uplot, outpath, "norm_u"//srcstring, nframe, par%nproc, glob_elem, glob_coord2, make_progress)
                    call generateTriMeshFiles(db, ux, outpath, "ux"//srcstring, nframe, par%nproc, glob_elem, glob_coord2, make_progress)
                    call generateTriMeshFiles(db, uz, outpath, "uz"//srcstring, nframe, par%nproc, glob_elem, glob_coord2, make_progress)
                end if
                if (movie%points) then
                    if (par%log) write (*,"(a36,a15,a29)") "| Creating points displacement files",logstring,"                            |"
                    call generatePointFiles(uplot, outpath, "norm_u"//srcstring, xyzplot, nframe, make_progress)
                    call generatePointFiles(ux, outpath, "ux"//srcstring, xyzplot, nframe, make_progress)
                    call generatePointFiles(uz, outpath, "uz"//srcstring, xyzplot, nframe, make_progress)
                end if
                deallocate(uplot)
                deallocate(ux)
                deallocate(uz)
            end if

            if (movie%velocity) then
                if (par%log) write (*,"(a25,a15,a40)") "| Creating Velocity files",logstring,"                                       |"
                allocate(vplot(glob_nglob,nframe))
                allocate(vel_x(glob_nglob,nframe))
                allocate(vel_z(glob_nglob,nframe))
                call readplotbin(par%nproc,movie%frame,par%nt,errmsg,vplot,outpath,"moviedata_normV"//srcstring,glob_nglob,nframe,ndata)
                call readplotbin(par%nproc,movie%frame,par%nt,errmsg,vel_x,outpath,"moviedata_vx"//srcstring,glob_nglob,nframe,ndata)
                call readplotbin(par%nproc,movie%frame,par%nt,errmsg,vel_z,outpath,"moviedata_vz"//srcstring,glob_nglob,nframe,ndata)
                if (movie%trimesh) then
                    if (par%log) write (*,"(a33,a15,a32)") "| Creating trimesh Velocity files",logstring,"                               |"
                    call generateTriMeshFiles(db, vplot, outpath, "norm_v"//srcstring, nframe, par%nproc, glob_elem, glob_coord2, make_progress)
                    call generateTriMeshFiles(db, vel_x, outpath, "vx"//srcstring, nframe, par%nproc, glob_elem, glob_coord2, make_progress)
                    call generateTriMeshFiles(db, vel_z, outpath, "vz"//srcstring, nframe, par%nproc, glob_elem, glob_coord2, make_progress)
                end if
                if (movie%points) then
                    if (par%log) write (*,"(a32,a15,a33)") "| Creating points Velocity files",logstring,"                                |"
                    call generatePointFiles(vplot, outpath, "norm_v"//srcstring, xyzplot, nframe, make_progress)
                    call generatePointFiles(vel_x, outpath, "vx"//srcstring, xyzplot, nframe, make_progress)
                    call generatePointFiles(vel_z, outpath, "vz"//srcstring, xyzplot, nframe, make_progress)
                end if
                deallocate(vplot)
                deallocate(vel_x)
                deallocate(vel_z)
            end if

            if (movie%stress) then
                if (par%log) write (*,"(a23,a15,a42)") "| Creating Stress files",logstring,"                                         |"
                allocate(sigmaxx(glob_nglob,nframe))
                allocate(sigmazz(glob_nglob,nframe))
                allocate(sigmaxz(glob_nglob,nframe))
                call readplotbin(par%nproc,movie%frame,par%nt,errmsg,sigmaxx,outpath,"moviedata_sigmaXX"//srcstring,glob_nglob,nframe,ndata)
                call readplotbin(par%nproc,movie%frame,par%nt,errmsg,sigmazz,outpath,"moviedata_sigmaZZ"//srcstring,glob_nglob,nframe,ndata)
                call readplotbin(par%nproc,movie%frame,par%nt,errmsg,sigmaxz,outpath,"moviedata_sigmaXZ"//srcstring,glob_nglob,nframe,ndata)
                if (movie%trimesh) then
                    if (par%log) write (*,"(a31,a15,a34)") "| Creating trimesh stress files",logstring,"                                 |"
                    call generateTriMeshFiles(db, sigmaxx, outpath, "sigmaxx"//srcstring, nframe, par%nproc, glob_elem, glob_coord2, make_progress)
                    call generateTriMeshFiles(db, sigmazz, outpath, "sigmazz"//srcstring, nframe, par%nproc, glob_elem, glob_coord2, make_progress)
                    call generateTriMeshFiles(db, sigmaxz, outpath, "sigmaxz"//srcstring, nframe, par%nproc, glob_elem, glob_coord2, make_progress)
                end if
                if (movie%points) then
                    if (par%log) write (*,"(a30,a15,a35)") "| Creating points stress files",logstring,"                                  |"
                    call generatePointFiles(sigmaxx, outpath, "sigmaxx"//srcstring, xyzplot, nframe, make_progress)
                    call generatePointFiles(sigmazz, outpath, "sigmazz"//srcstring, xyzplot, nframe, make_progress)
                    call generatePointFiles(sigmaxz, outpath, "sigmaxz"//srcstring, xyzplot, nframe, make_progress)
                end if
                deallocate(sigmaxx)
                deallocate(sigmazz)
                deallocate(sigmaxz)
            end if
        enddo

        if (allocated(db)) deallocate(db)
        if (allocated(ndata)) deallocate(ndata)
        if (allocated(ndata2)) deallocate(ndata2)
        if (allocated(ndata3)) deallocate(ndata3)
        if (allocated(xyzplot)) deallocate(xyzplot)
        if (allocated(glob_elem)) deallocate(glob_elem)
        if (allocated(glob_coord)) deallocate(glob_coord)
    end subroutine collectMovie

    subroutine readplotbin(nproc,nstep,nt,errmsg,dat,path,name1,glob_nglob,nframe,ndata)
        !input
        integer :: nproc
        integer :: nt
        integer :: glob_nglob
        integer :: nframe
        integer :: nstep
        type(error_message) :: errmsg
        character(len=*) :: name1, path

        real(kind=CUSTOM_REAL), dimension(glob_nglob,nframe):: dat
        integer, dimension(:):: ndata
        !local
        integer :: iproc, c, it,e,d

        d=1
        e=ndata(1)

        do iproc=1,nproc
            c=1
            do it=nstep,nt,nstep
                if (name1 == 'moviedata_uz') then
                    call readPlotBinaries(d, e, c, dat, path, name1, iproc, it, errmsg)
                else
                    call readPlotBinaries(c, d, e, dat, path, name1, iproc, it, errmsg)
                end if
                c=c+1
            end do
            if (iproc < nproc) then
                d=d+ndata(iproc)
                e=e+ndata(iproc+1)
            end if
        end do
    end subroutine readplotbin

    subroutine readplotbinElem(nproc,iter,errmsg,dat,path,name1,glob_nelem,ndata)
        !input
        integer :: nproc
        integer :: iter
        integer :: glob_nelem
        type(error_message) :: errmsg
        character(len=*) :: name1, path

        real(kind=CUSTOM_REAL), dimension(glob_nelem):: dat
        integer, dimension(:):: ndata
        !local
        integer :: iproc,e,d

        d=1
        e=ndata(1)

        do iproc=1,nproc
            call readPlotBinariesInv(d, e, dat, path, name1, iproc, iter, errmsg)
            if (iproc < nproc) then
                d=d+ndata(iproc)
                e=e+ndata(iproc+1)
            end if
        end do
    end subroutine readplotbinElem
end module collectMovieMod
