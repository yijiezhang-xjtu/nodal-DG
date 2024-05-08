module timeloopMod
    ! module to calculate the timeloop
    use constantsMod
    use meshMod
    use parameterMod
    use waveMod
    use stfMod
    use sourceReceiverMod
    use matrixMod
    use plotMod
    use mpiMod
    use fileunitMod
    use pmlMod
    use errorMessage
    use calcFluxMod
    use timestampMod
    use derMod
    use adjointMod
    use, intrinsic :: iso_fortran_env

    implicit none

    contains

	subroutine timeloop2d(par, inv_key, iter_step, run_number, act_src_temp, movie, myrank, errmsg)
        type(error_message) :: errmsg
        type(parameterVar) :: par
        type(meshVar) :: mesh
        type(srcVar) ::src
        type(recVar) ::rec
        type(movie_parameter) :: movie
        ! time variables
        real(kind=CUSTOM_REAL) :: dt
        real(kind=CUSTOM_REAL) :: f0,f0tmp
        real(kind=CUSTOM_REAL), pointer, dimension(:) :: t0
        real(kind=CUSTOM_REAL) :: t0max, t0maxtmp
        real(kind=CUSTOM_REAL) :: time
        real(kind=CUSTOM_REAL) :: simt0
        real(kind=CUSTOM_REAL) :: fcrit,avg_energy1,avg_energy2,sta_lta
        integer :: timecrit
        logical :: pmlcrit = .true.
        ! free
        real(kind=custom_real), dimension(:,:), allocatable :: free
        ! indices
        integer :: i,j,r, ind, counter
        integer :: iglob
        integer, dimension(Np) :: iv
        integer :: ie, it, is
        integer :: ier
        ! RK
        integer :: irk, nrk, wrk
        real(kind=custom_real), dimension(:,:), allocatable :: resQ, resU, resQ2, resQ3, resQ4 !Runge-Kutta residual for 4th order rk
        ! rec variabels
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: recInt,recTemp
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: plotv,plotw,plotux,plotuz,plotax,plotaz,plotp1,plotp2,plotv1x,plotv1z,plotv2x,plotv2z,plot_r,plot_t
        real(kind=CUSTOM_REAL), dimension(1) :: r_v,s_v
        real(kind=CUSTOM_REAL) :: ux_temp, uz_temp, vx_temp, vz_temp, ax_temp, az_temp, angle_temp
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: energy,all_energy, energy_kin, energy_pot, np_zeros
        integer, dimension(:), pointer :: recelemv
        integer :: nsamples, isample
        ! source
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: srcInt,srcTemp
        real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: plotstf
        real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: plotDiffstf
        real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: srcArray
        real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: srcArrayM
        integer, dimension(:), pointer :: srcelemv
        real(kind=CUSTOM_REAL), dimension(2,2) :: M,Ms,Ma
        real(kind=CUSTOM_REAL) :: rk_time, stf_val, stf_val_2, stf_val_diff
        integer :: n_src_tmp, n_src
        ! PML
        logical, dimension(:), allocatable :: pmlcheck
        integer, dimension(:,:), allocatable :: pmlloc
        real(kind=CUSTOM_REAL) :: amax
        real(kind=CUSTOM_REAL) :: xmean, zmean
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: ddx,ddz,alphax,alphaz,kx,kz
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: fprime, gprime
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: resFprime, resGprime, resFprime2, resGprime2, resFprime3, resGprime3, resFprime4, resGprime4  !runge Kutte 4th order residual storage
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: fprimen, fprimem, gprimen, gprimem

        ! fields
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: ux,uz,ax,az, uplot, vplot, v1plot, v2plot
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: rQ,Q,Qn,Qm, rQm,e,u
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: ftemp,gtemp,htemp,qtemp
        real(kind=CUSTOM_REAL), dimension(Np) :: dFdr,dFds,dGdr,dGds
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: flux
        real(kind=CUSTOM_REAL), dimension(Np,Np) :: invmass, vdmTinv, mass
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: stf, stf_old
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: anelasticvar
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: elasticfluxvar
        real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: APA, aAPA
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: T, invT, VT, VTfree, aVT, aVTfree, at

        !output
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: div
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: curl

        !logical :: movie
        logical :: file_exists
        character(len=80) ::filename
        integer :: time_shift_movie

        ! usefulls
        real(kind=CUSTOM_REAL) :: onehalf = 1./2.
        real(kind=CUSTOM_REAL) :: onethree = 1./3.
        real(kind=CUSTOM_REAL) :: twothree = 2./3.
        real(kind=CUSTOM_REAL) :: onefor = 1./4.
        real(kind=CUSTOM_REAL) :: threefor = 3./4.
		real(kind=CUSTOM_REAL) :: onesix = 1./6.
        real(kind=CUSTOM_REAL) :: dummy
        character(len=80) :: int_string

        ! timer
        type(timestampVar) :: timestamp
        real(kind = custom_real) :: localtime
        ! MPI
        integer :: myrank
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: q_send
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: q_rec
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: qi
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: qi_test
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: rq_send
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: rq_rec
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rqi
        integer :: c,k
        integer :: dest
        integer :: req, req_r, tag, tag_r
        integer ,dimension(:), allocatable:: req1
        real(kind=CUSTOM_REAL) :: maxv
        real(kind=CUSTOM_REAL) :: maxu

        integer :: dimens
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: a_sigma_tmp
        character(len=7) :: inv_key
		integer :: iter_step, run_number, forward_steps, i_src
        logical, dimension(:) :: act_src_temp
        logical, dimension(:), allocatable :: act_src
        character(len=7) :: srcstring

        ! error message
        character(len=10) :: myname = "timeloop2d"

        call addTrace(errmsg, myname)

        write(filename,"('meshVar',i6.6)") myrank+1
        inquire(file=trim(outpath)//trim(filename), exist=file_exists)
        if (file_exists) then
            call readMeshVar(mesh,trim(outpath)//filename, errmsg)
        else
            call add(errmsg, 2, "Error in databases, file "//trim(filename)//" does not exist!", myname, filename)
            call print(errmsg)
            call stop_mpi()
        end if

        if (par%autodt) then
            dt = mesh%dtfactor*par%cfl
        else
            dt = par%dt
        endif
        if (par%autont) then
            par%nt = int(par%t_total/dt) + 1
        endif

        if (mesh%has_src) then
            ! load sources
            write(filename,"('srcVar',i6.6)") myrank+1
            inquire(file=trim(outpath)//trim(filename), exist=file_exists)
            if (file_exists) then
                call readSrcVar(src,trim(outpath)//filename)
            else
                call add(errmsg, 2, "File does not exist!", myname, filename)
                call print(errmsg)
                call stop_mpi()
            end if
            f0tmp = maxval(src%srcf0)
        else
            f0tmp = 0.0
        endif
        ! get the maximum f0
        call maxval_real_all(f0tmp,f0,CUSTOM_REAL)

        ! load receivers
        if (mesh%has_rec) then
            write(filename,"('recVar',i6.6)") myrank+1
            inquire(file=trim(outpath)//trim(filename), exist=file_exists)
            if (file_exists) then
                call readRecVar(rec,trim(outpath)//filename)
            else
                call add(errmsg, 2, "Error in recVar, file does not exist!", myname, filename)
                call print(errmsg)
                call stop_mpi()
            end if
        end if

        call sync_mpi()

        ! allocate fields
        dimens = 5
        allocate(ux(mesh%nglob),uz(mesh%nglob),ax(mesh%nglob),az(mesh%nglob),uplot(mesh%nglob),vplot(mesh%nglob),v1plot(mesh%nglob),v2plot(mesh%nglob))
        allocate(rQ(mesh%nglob,dimens),Q(mesh%nglob,dimens),Qn(mesh%nglob,dimens),Qm(mesh%nglob,dimens), rQm(mesh%nglob, dimens), resQ(mesh%nglob,dimens), resU(mesh%nglob, 2), u(mesh%nglob, 2))
		allocate(resQ2(mesh%nglob,dimens), resQ3(mesh%nglob,dimens), resQ4(mesh%nglob,dimens))
        allocate(ftemp(Np,dimens),gtemp(Np,dimens),htemp(Np,dimens),qtemp(Np,dimens))
        allocate(e(mesh%nglob,3))
        allocate(flux(3*NGLL,dimens))
        allocate(stf(par%nt))
        allocate(energy(par%nt),energy_kin(par%nt), energy_pot(par%nt))
        allocate(all_energy(par%nt), a_sigma_tmp(Np, 3))
        allocate(act_src(size(act_src_temp)))
        act_src=act_src_temp
        energy     = 0.0
        energy_kin = 0.0
        energy_pot = 0.0
        all_energy = 0.0

        if (mesh%has_src) then
            allocate(srcInt(Np,src%nsrc),srcTemp(Np,src%nsrc))
            allocate(t0(src%nsrc))
            if (inv_key == 'forward') then
                allocate(srcArray(Np,2,src%nsrc))
            end if
            allocate(srcArrayM(Np,3,src%nsrc))
            allocate(plotstf(par%nt,2,src%nsrc))
            allocate(plotDiffstf(par%nt,2,src%nsrc))
        end if
        if (mesh%has_rec) then
            allocate(recInt(Np,rec%nrec),recTemp(Np,rec%nrec))
            nsamples = int(par%nt/par%subsampling_factor)
            allocate(plotv(rec%nrec,nsamples), plotw(rec%nrec,nsamples))
            allocate(plotux(rec%nrec,nsamples), plotuz(rec%nrec,nsamples))
            allocate(plotax(rec%nrec,nsamples), plotaz(rec%nrec,nsamples))
            allocate(plot_r(rec%nrec,nsamples), plot_t(rec%nrec,nsamples))
        end if

        ! mpi
        allocate(q_send(NGLL*mesh%mpi_ne*dimens,mesh%mpi_nn))
        allocate(q_rec(NGLL*mesh%mpi_ne*dimens,mesh%mpi_nn))
        allocate(qi(NGLL,dimens,mesh%mpi_ne,mesh%mpi_nn))
        allocate(qi_test(NGLL,dimens,mesh%mpi_ne,mesh%mpi_nn))
        allocate(req1(mesh%mpi_nnmax))
        allocate(rq_send(NGLL*mesh%mpi_ne*dimens,mesh%mpi_nn))
        allocate(rq_rec(NGLL*mesh%mpi_ne*dimens,mesh%mpi_nn))
        allocate(rqi(NGLL,dimens,mesh%mpi_ne,mesh%mpi_nn))
        ! pml
        allocate(pmlcheck(mesh%nelem))
        pmlcheck=.false.
        if(par%set_pml) then
            allocate(fprime(mesh%nglob,dimens))
            allocate(gprime(mesh%nglob,dimens))
            allocate(resFprime(mesh%nglob,dimens),resFprime2(mesh%nglob,dimens),resFprime3(mesh%nglob,dimens),resFprime4(mesh%nglob,dimens))
            allocate(resGprime(mesh%nglob,dimens),resGprime2(mesh%nglob,dimens),resGprime3(mesh%nglob,dimens),resGprime4(mesh%nglob,dimens))
            allocate(fprimen(mesh%nglob,dimens),fprimem(mesh%nglob,dimens))
            allocate(gprimen(mesh%nglob,dimens),gprimem(mesh%nglob,dimens))
            allocate(np_zeros(mesh%nglob))
        end if
        if (mesh%has_rec) then
            plotux = 0.0
            plotuz = 0.0
            plotv  = 0.0
            plotw  = 0.0
            plotax = 0.0
            plotaz = 0.0
        end if

        q      = 1e-24!eps
        e      = 1e-24!eps
        rQ     = 1e-24!eps
        Qn     = 1e-24!eps
        resQ   = 1e-24!eps
		resQ2  = 1e-24!eps
		resQ3  = 1e-24!eps
		resQ4  = 1e-24!eps

        ux     = 0.0
        uz     = 0.0
        ax     = 0.0
        az     = 0.0
        uplot  = 0.0
        vplot  = 0.0
        v1plot = 0.0
        v2plot = 0.0

        forward_steps = 0
        time_shift_movie = 0
        timecrit = 0

        if (par%use_trigger) then
            fcrit    = 5/f0
            timecrit = int(((1.2/f0)*3)/dt)
        endif

        nrk = 0
        wrk = 0
        if (par%timeint == 1) then
            write(int_string, '(a80)') "|                        use runge kutta fourth order                          |"
            nrk = 4 ! rk4
            wrk = 1
        else if (par%timeint == 2) then
            write(int_string, '(a80)') "|                       use tvd runge kutta third order                        |"
            nrk = 3 ! rk3
            wrk = 2
        else if (par%timeint == 3) then
            write(int_string, '(a80)') "|                       use les-runge kutta fourth order                       |"
            nrk = 5 ! les-rk4
            wrk = 3
        end if
        if (par%log.and.myrank==0) write(*,'(a80)') "|------------------------------------------------------------------------------|"
        if (par%log.and.myrank==0) write(*,'(a80)') int_string

        if(par%set_pml) then
            if (par%log.and.myrank==0) then
                write(*,"(a80)") "|------------------------------------------------------------------------------|"
                write(*,'(a80)') "|                        using NPML boundary conditions                        |"
                write(*,'(a40, f5.2, a35)') "|                                 xmin: ", mesh%pmlxmin, "                                  |"
                write(*,'(a40, f5.2, a35)') "|                                 xmax: ", mesh%pmlxmax, "                                  |"
                write(*,'(a40, f5.2, a35)') "|                                 zmin: ", mesh%pmlzmin, "                                  |"
                write(*,'(a40, f5.2, a35)') "|                                 zmax: ", mesh%pmlzmax, "                                  |"
                write(*,'(a40, f5.2, a35)') "|                            pml delta: ", par%pml_delta, "                                  |"
                write(*,'(a40, f5.2, a35)') "|                                   rc: ", par%pml_rc, "                                  |"
                write(*,'(a40, f5.2, a35)') "|                                 kamx: ", par%pml_kmax, "                                  |"
                write(*,'(a40, f5.2, a35)') "|                                 afac: ", par%pml_afac, "                                  |"
                write(*,'(a37,f8.2,a4,f8.2,a23)') '|                 x axis ranges from ', mesh%pmlxmin, ' to ', mesh%pmlxmax, ' |'
                write(*,'(a37,f8.2,a4,f8.2,a23)') '|                 z axis ranges from ', mesh%pmlzmin, ' to ', mesh%pmlzmax, ' |'
            end if
            amax=par%pml_afac*pi*f0

            allocate(ddx(mesh%nglob), ddz(mesh%nglob))
            allocate(alphax(mesh%nglob), alphaz(mesh%nglob))
            allocate(kx(mesh%nglob),kz(mesh%nglob))
            allocate(pmlloc(mesh%nelem,2))
            ddx=0
            ddz=0
            kx=1
            kz=1
            alphax=0.
            alphaz=0.
            pmlloc=0
            pmlcheck=.false.
            do ie=1, mesh%nelem
                iv=mesh%ibool(:,ie)
                if (mesh%pml(ie)>0) then
                    pmlcheck(ie)=.true.      ! flag element as PML element
                    xmean=0.
                    zmean=0.
                    do j=1,NP               ! get center of element to determine its position
                        xmean=xmean+mesh%vx(iv(j))/Np
                        zmean=zmean+mesh%vz(iv(j))/Np
                    enddo
                    if (xmean<mesh%pmlxmin+par%pml_delta) pmlloc(ie,1)=-1     ! define to which PML (x,z) an element belongs (possibly multiple)
                    if (xmean>mesh%pmlxmax-par%pml_delta) pmlloc(ie,1)=1
                    if (zmean<mesh%pmlzmin+par%pml_delta) pmlloc(ie,2)=-1
                    if (zmean>mesh%pmlzmax+par%pml_delta) pmlloc(ie,2)=1
                    call dampingProfile(ddx,ddz,alphax,alphaz,kx,kz,mesh%vx,mesh%vz,mesh%vp(ie),mesh%pmlxmin,mesh%pmlxmax,mesh%pmlzmin,mesh%pmlzmax,par%pml_delta,par%pml_rc,amax,par%pml_kmax,iv,pmlloc(ie,:))    ! get damping profiles
                endif
            enddo
            fprime = 0.
            gprime = 0.
            if (wrk == 3) then !les-rk4
                resFprime = 0.
                resGprime = 0.
            end if
			if (wrk == 1) then !rk4
                resFprime = 0.
                resGprime = 0.
				resFprime2 = 0.
                resGprime2 = 0.
				resFprime3 = 0.
                resGprime3 = 0.
				resFprime4 = 0.
                resGprime4 = 0.
            end if
        end if

        ! free surface conditions
        allocate(free(dimens,dimens))
        free=0.0
        free(1,1)=-2
        free(2,2)=0
        free(3,3)=-2
        free(4,4)=0
        free(5,5)=0

        ftemp=0.
        gtemp=0.

        allocate(elasticfluxvar(mesh%nelem,4))
        allocate(APA(mesh%nelem,dimens,dimens))
        allocate(T(mesh%nelem,3,dimens,dimens), invT(mesh%nelem,3,dimens,dimens))
        allocate(VT(mesh%nelem,3,dimens,dimens), VTfree(mesh%nelem,3,dimens,dimens))

        do ie=1,mesh%nelem
            elasticfluxvar(ie,1) = -(mesh%lambda(ie)+2*mesh%mu(ie))
            elasticfluxvar(ie,2) = -mesh%lambda(ie)
            elasticfluxvar(ie,3) = -mesh%mu(ie)
            elasticfluxvar(ie,4) = -1.0/mesh%rho(ie)
            APA(ie,:,:) = getAPA(mesh%vpu(ie),mesh%vsu(ie),mesh%rho(ie),mesh%lambda(ie),mesh%mu(ie))
            do is=1,3
                T(ie,is,:,:) =getT(mesh%nx(is*NGLL,ie),mesh%nz(is*NGLL,ie),dimens)
                invT(ie,is,:,:) =getinvT(mesh%nx(is*NGLL,ie),mesh%nz(is*NGLL,ie),dimens)
                VT(ie,is,:,:) = matmul(T(ie,is,:,:),matmul(APA(ie,:,:),invT(ie,is,:,:)))
                VTfree(ie,is,:,:) = matmul(T(ie,is,:,:),matmul(matmul(APA(ie,:,:),free),invT(ie,is,:,:)))
            enddo
        enddo

        vdmTinv=mesh%vdm
        vdmTinv=transpose(vdmTinv)
        call invert(vdmTinv, errmsg)
        allocate(srcelemV(mesh%nelem))
        allocate(recelemV(mesh%nelem))
        srcelemv=0
        recelemv=0

        ms = 0.
        if (inv_key=='forward') then
            if (mesh%has_src) then
                do i=1,src%nsrc
                    r_v(1) = src%srcrs(1,i)
                    s_v(1) = src%srcrs(2,i)
                    call vdm2D(srcTemp(:,i),r_v,s_v)
                    srcInt(:,i)=matmul(vdmTinv,srcTemp(:,i))
                    t0(i)=1.2/src%srcf0(i)
                    srcelemv(src%srcelem(i)) = i

                    if (src%srctype(i) == 1) then  ! momenttensor
                        M(1,1)=src%srcm(1,i) !Mxx
                        M(1,2)=src%srcm(3,i) !Mxz
                        M(2,1)=src%srcm(3,i) !Mzx
                        M(2,2)=src%srcm(2,i) !Mzz
                        do j=1,2
                            do k=1,2
                                ms(j,k)=0.5*(M(j,k)+M(k,j))
                                ma(j,k)=0.5*(M(j,k)-M(k,j))
                            end do
                        end do
                    end if
                end do
                t0maxtmp = maxval(t0)
            else
                t0maxtmp=0.
            end if
        endif

        call sync_mpi()
        call maxval_real_all(t0maxtmp,t0max,CUSTOM_REAL)

        if (.not. par%shift_sources) then
            if (mesh%has_src) then
                do i=1,src%nsrc
                    t0(i)=0.
                enddo
            endif

            if (-t0max < par%simt0) then
                simt0 = -t0max
                if (myrank == 0) then
                    call add(errmsg, 1, "Source(s) started earlier than chosen simt0, so we changed simt0!", myname)
                end if
                if (par%log .and. myrank == 0) then
                    write(*,'(a80)') "|------------------------------------------------------------------------------|"
                    write(*,'(a80)') "|        Source(s) start earlier than chosen simt0, so we change simt0!        |"
                    write(*,'(a40, f10.7, a30)') "|                       original simt0: ", par%simt0, "                             |"
                    write(*,'(a40, f10.7, a30)') "|                            new simt0: ", simt0, "                             |"
                end if
            else
                simt0 = par%simt0
            endif
        else
            simt0 = par%simt0
        endif

        if (par%autoshift) then
            if (par%shift_sources) then
                ! This shifts the seismograms so that t=0 actually is the maximum of the wavelet with the smallest frequency (Ricker, Gaussian).
                par%plott0 = t0max ! 1.2/f0
            else
                par%plott0 = 0.
            end if
            ! else par%plott0 remains the same as provided in the parfile
        end if

        ! set up rec interpolation
        if (mesh%has_rec) then
            do i=1,rec%nrec
                r_v(1) = rec%recrs(1,i)
                s_v(1) = rec%recrs(2,i)
                call vdm2D(recTemp(:,i),r_v,s_v)
                recInt(:,i)=matmul(vdmTinv,recTemp(:,i))
                recelemv(rec%recelem(i)) = i
            end do
        end if

        ! choose source time function
        if (inv_key=='forward') then
            if (mesh%has_src) then
                do i=1,src%nsrc
                    if (act_src(src%srcnr(i))) then
                        select case (src%srcstf(i))
                            case (1) !GAUSS
                                do it=1,par%nt
                                    time = (float(it)-1.)*dt + simt0
                                    select case (src%srctype(i))
                                        case (0)  ! single force
                                            plotstf(it,1,i) = time
                                            plotstf(it,2,i) = -stfGauss(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                                        case (1)  ! moment tensor
                                            plotDiffstf(it,2,i) = -stfDiffGauss(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                                            plotDiffstf(it,1,i) = time
                                    end select
                                enddo
                            case (2) !RICKER
                                do it=1,par%nt
                                    time = (float(it)-1.)*dt + simt0
                                    select case (src%srctype(i))
                                        case (0)  ! single force
                                            plotstf(it,1,i) = time
                                            plotstf(it,2,i) = -stfRicker(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                                        case (1)  ! moment tensor
                                            plotDiffstf(it,1,i) = time
                                            plotDiffstf(it,2,i) = -stfDiffRicker(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                                    end select
                                enddo
                            case (3) !SIN^3
                                do it=1,par%nt
                                    time = (float(it)-1.)*dt + simt0
                                    select case (src%srctype(i))
                                        case (0)  ! single force
                                            plotstf(it,1,i) = time
                                            plotstf(it,2,i) = -stfSin3(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                                        case (1)  ! moment tensor
                                            plotDiffstf(it,1,i) = time
                                            plotDiffstf(it,2,i) = -stfDiffSin3(time,src%srcf0(i),t0(i)+src%delay(i),src%srcfactor(i))
                                    end select
                                enddo
                            case (4) !EXTERNAL
                                filename=src%extwavelet(i)
                                select case (src%srctype(i))
                                    case (0)  ! single force
                                        call stfExternal(plotstf(:,2,i),plotstf(:,1,i),dt,par%nt,.true.,6,3,trim(filename),0, errmsg)
                                    case (1)  ! moment tensor
                                        call stfExternal(plotDiffstf(:,2,i),plotDiffstf(:,1,i),dt,par%nt,.true.,6,3,trim(filename),1, errmsg)
                                end select
                            case default
                                call add(errmsg, 2, "Chose a valid source time function. Available functions are listed in the 'source' parameter file.", myname, "data/source")
                                call print(errmsg)
                                call stop_mpi()
                        end select

                        if (run_number == 1) then
                            ! write stf
                            select case (src%srctype(i))
                                case (0)  ! single force
                                    write(filename,"('stf',i6.6)") src%srcnr(i)
                                    
                                    open(unit=27,file=trim(outpath)//trim(filename),status='unknown')
                                    do it=1,par%nt
                                        write(27,*) plotstf(it,1,i)-par%plott0,plotstf(it,2,i)
                                    end do
                                    close(27)
                                case (1)  ! moment tensor
                                    write(filename,"('stfdiff',i6.6)") src%srcnr(i)
                                    
                                    open(unit=27,file=trim(outpath)//trim(filename),status='unknown')
                                    do it=1,par%nt
                                        write(27,*) plotDiffstf(it,1,i)-par%plott0,plotDiffstf(it,2,i)
                                    end do
                                    close(27)
                            end select
                        end if
                    end if
                end do
            end if
        endif

        qm = q
        tag = 0
        q_send = 0.
        do i=1,mesh%mpi_nn
            do ie=1,mesh%mpi_ne ! loop over interface elements
                do k=1,dimens
                    do j=1,NGLL
                        if ( mesh%mpi_connection(i,ie,1) >0) then
                            q_send((ie-1)*dimens*NGLL + (k-1)*NGLL + j,i) = &
                                qm(mesh%ibt(j,mesh%mpi_connection(i,ie,2),mesh%mpi_connection(i,ie,1)),k)
                        end if
                    end do
                end do
            end do
        end do ! all interfaces

        !send and rec
        do i=1,mesh%mpi_nn
            dest=mesh%mpi_neighbor(i)-1
            call isendV_real(q_send(:,i),(mesh%mpi_ne*dimens*NGLL),dest,tag,req,CUSTOM_REAL)
            call irecV_real(q_rec(:,i),(mesh%mpi_ne*dimens*NGLL),dest,tag,req_r,CUSTOM_REAL)
            call wait_req(req)
            call wait_req(req_r)
        end do

        !unpack mpi buffers
        do i=1,mesh%mpi_nn
            c = 1
            do ie=1,mesh%mpi_ne
                do k=1,dimens
                    do j=1,NGLL
                        if ( mesh%mpi_connection(i,ie,1) > 0) then
                            qi(j,k,ie,i) = q_rec(c,i)
                        end if
                        c = c+1
                    end do
                end do
            end do
        end do


        call sync_mpi()

        invmass=matmul(mesh%vdm,transpose(mesh%vdm))
        mass=invmass
        call invert(mass, errmsg)

        n_src_tmp=0
        if (mesh%has_src) then
            n_src_tmp=src%nsrc
        endif
        call sum_int_all(n_src_tmp, n_src)

        if (par%log.and.myrank==0) then
            call initialTimestamp(timestamp)
            write(*,"(a80)") "|------------------------------------------------------------------------------|"
            write(*,"(a80)") "|                             starting timeloop...                             |"
            write(*, "(a40, es12.5, a28)") "|                            Time step: ", dt, "                           |"
            write(*, "(a40, i12, a28)")    "|                 Number of time steps: ", par%nt, "                           |"
            write(*, "(a40, es12.5, a28)") "|                 Total simulated time: ", dt*par%nt, "                           |"
			write(*,"(a80)") "|------------------------------------------------------------------------------|"
        end if

        do it=1,par%nt
            stf(it) = (float(it)-1.) * dt + simt0
        enddo
        
        pmlcrit=.true.

        ! ---------------------------------------------------------------------------------------------
        ! ------------------------------------timeloop-------------------------------------------------
        ! ---------------------------------------------------------------------------------------------

        do it=1,par%nt ! timeloop
            if (myrank == 0) then
                localtime = (float(it)-1.) * dt + simt0
            end if
            rk_time=0.
			if (inv_key=='forward') then
                if (mesh%has_rec .and. mod(it,par%subsampling_factor) == 0) then
                    isample = int(it/par%subsampling_factor)
                    do r=1,rec%nrec
                        iv = mesh%ibool(:,rec%recelem(r))
                        ! Calculate divergence and curl of the velocity component in order to create separate seismograms for
                        ! Radial and tangential component
                        if (par%curl) then
                            call curl2d(q(iv,4),q(iv,5),mesh%rx(iv),mesh%sx(iv),mesh%rz(iv),mesh%sz(iv),mesh%Dr,mesh%Ds,curl)
                            do j = 1, Np
                                plot_t(r, isample) = plot_t(r, isample) + curl(j,2)*recint(j,r)
                            end do
                        end if
                        if (par%div) then
                            call div2d(div,q(iv,4),q(iv,5),mesh%Dr,mesh%Ds,mesh%rx(iv),mesh%sx(iv),mesh%rz(iv),mesh%sz(iv))
                            do j = 1, Np
                                plot_r(r, isample) = plot_r(r, isample) + div(j)*recint(j,r)
                            end do
                        end if

                        do j=1,Np
                            iglob=mesh%ibool(j,rec%recelem(r))
                            plotux(r,isample) = plotux(r,isample) + ux(iglob)  * recint(j,r)
                            plotuz(r,isample) = plotuz(r,isample) + uz(iglob)  * recint(j,r)
                            plotv(r,isample)  = plotv(r,isample)  + q(iglob,4) * recint(j,r)
                            plotw(r,isample)  = plotw(r,isample)  + q(iglob,5) * recint(j,r)
                            plotax(r,isample) = plotax(r,isample) + ax(iglob)  * recint(j,r)
                            plotaz(r,isample) = plotaz(r,isample) + az(iglob)  * recint(j,r)
                            !if (par%poroelastic .and. mesh%nfluids >= 1) then
                            !    plotp1(r,isample)  = plotp1(r,isample) + q(iglob,6) * recint(j,r)
                            !    plotv1x(r,isample) = plotv1x(r,isample) + q(iglob,7) * recint(j,r)
                            !    plotv1z(r,isample) = plotv1z(r,isample) + q(iglob,8) * recint(j,r)
                            !    if (mesh%nfluids == 2) then
                            !        plotp2(r,isample)  = plotp2(r,isample) + q(iglob,9) * recint(j,r)
                            !        plotv1x(r,isample) = plotv1x(r,isample) + q(iglob,7) * recint(j,r)
                            !        plotv1z(r,isample) = plotv1z(r,isample) + q(iglob,8) * recint(j,r)
                            !    endif
                            !endif
                        end do
                        ! rotate receivers
                        angle_temp = (rec%rec_ang(r)*PI)/180.
                        ux_temp = cos(angle_temp) * plotux(r,isample) + sin(angle_temp) * plotuz(r,isample)
                        uz_temp = -sin(angle_temp) * plotux(r,isample) + cos(angle_temp) * plotuz(r,isample)
                        vx_temp = cos(angle_temp) * plotv(r,isample) + sin(angle_temp) * plotw(r,isample)
                        vz_temp = -sin(angle_temp) * plotv(r,isample) + cos(angle_temp) * plotw(r,isample)
                        ax_temp = cos(angle_temp) * plotax(r,isample) + sin(angle_temp) * plotaz(r,isample)
                        az_temp = -sin(angle_temp) * plotax(r,isample) + cos(angle_temp) * plotaz(r,isample)
                        plotux(r,isample) = ux_temp
                        plotuz(r,isample) = uz_temp
                        plotv(r,isample) = vx_temp
                        plotw(r,isample) = vz_temp
                        plotax(r,isample) = ax_temp
                        plotaz(r,isample) = az_temp
                    end do
                end if
            endif

            tag   = 0
            tag_r = 0
            stf_val = 0.
            stf_val_diff = 0.
            do irk=1,nrk ! runge-kutta loop
                time = (float(it) - 1. + rk_time) * dt + simt0
                !if (.not.par%inversion .or. (par%inversion .and. inv_key=='forward')) then
				if (inv_key=='forward') then
                    if (mesh%has_src) then
                        do i=1, src%nsrc
                            if (act_src(src%srcnr(i))) then
                                select case (src%srctype(i))
                                    case (0)  ! single force
                                        call cubicInterpolation(dble(plotstf(:,1,i)-par%plott0),dble(plotstf(:,2,i)),dble(time-par%plott0), stf_val)
                                        srcArray(:,1,i)=srcInt(:,i)*(-sin((src%srcangle_force(i)*PI)/180.)*stf_val)/mesh%rho(src%srcelem(i))
                                        srcArray(:,2,i)=srcInt(:,i)*(+cos((src%srcangle_force(i)*PI)/180.)*stf_val)/mesh%rho(src%srcelem(i))
                                        
                                        iv=mesh%ibool(:,src%srcelem(i))
                                        srcArray(:,1,i)=matmul(invmass,srcArray(:,1,i))/mesh%jacobian(iv)
                                        srcArray(:,2,i)=matmul(invmass,srcArray(:,2,i))/mesh%jacobian(iv)
                                        
                                    case (1)  ! moment tensor
                                        call cubicInterpolation(dble(plotDiffstf(:,1,i)-par%plott0),dble(plotDiffstf(:,2,i)),dble(time-par%plott0), stf_val_diff)
                                        srcArrayM(:,1,i)=srcInt(:,i)*(ms(1,1)*stf_val_diff)
                                        srcArrayM(:,2,i)=srcInt(:,i)*(ms(2,2)*stf_val_diff)
                                        srcArrayM(:,3,i)=srcInt(:,i)*(ms(1,2)*stf_val_diff)

                                        iv=mesh%ibool(:,src%srcelem(i))
                                        srcArrayM(:,1,i)=matmul(invmass,srcArrayM(:,1,i))/mesh%jacobian(iv)
                                        srcArrayM(:,2,i)=matmul(invmass,srcArrayM(:,2,i))/mesh%jacobian(iv)
                                        srcArrayM(:,3,i)=matmul(invmass,srcArrayM(:,3,i))/mesh%jacobian(iv)
                                end select
                            endif
                        end do
                    end if
                endif

                qm=q
                rQm = rQ
                if(par%set_pml) then
                    fprimem=fprime
                    gprimem=gprime
                end if
                q_send(:,:) = 0
                rq_send(:,:) = 0
                ! mpi comunication
                ! build send buffer
                do i=1,mesh%mpi_nn
                    do ie=1,mesh%mpi_ne ! loop over interface elements
                        do k=1,dimens
                            do j=1,NGLL
                                if ( mesh%mpi_connection(i,ie,1) >0) then
                                    q_send((ie-1)*dimens*NGLL + (k-1)*NGLL + j,i) = &
                                        qm(mesh%ibt(j,mesh%mpi_connection(i,ie,2),mesh%mpi_connection(i,ie,1)),k)
                                    rq_send((ie-1)*dimens*NGLL + (k-1)*NGLL + j,i) = &
                                        rqm(mesh%ibt(j,mesh%mpi_connection(i,ie,2),mesh%mpi_connection(i,ie,1)),k)
                                end if
                            end do
                        end do
                    end do
                end do ! all interfaces

                ! send and rec
                do i=1,mesh%mpi_nn
                    dest=mesh%mpi_neighbor(i)-1
                    call isendV_real(q_send(:,i),(mesh%mpi_ne*dimens*NGLL),dest,tag,req,CUSTOM_REAL)
                    call irecV_real(q_rec(:,i),(mesh%mpi_ne*dimens*NGLL),dest,tag,req_r,CUSTOM_REAL)
                    call wait_req(req)
                    call wait_req(req_r)
                    call isendV_real(rq_send(:,i),(mesh%mpi_ne*dimens*NGLL),dest,tag_r,req,CUSTOM_REAL)
                    call irecV_real(rq_rec(:,i),(mesh%mpi_ne*dimens*NGLL),dest,tag_r,req_r,CUSTOM_REAL)
                    call wait_req(req)
                    call wait_req(req_r)
                end do

                ! unpack mpi buffers
                do i=1,mesh%mpi_nn
                    c = 1
                    do ie=1,mesh%mpi_ne
                        do k=1,dimens
                            do j=1,NGLL
                                if ( mesh%mpi_connection(i,ie,1) > 0) then
                                    qi(j,k,ie,i) = q_rec(c,i)
                                    rqi(j,k,ie,i) = rq_rec(c,i)
                                end if
                                c = c+1
                            end do
                        end do
                    end do
                end do

                do ie=1,mesh%nelem ! element loop
                    iv=mesh%ibool(:,ie)
                    call calcFlux(ie, iv, mesh, elasticfluxvar(ie,:), free, qm, qi, rQ, ftemp, gtemp, fprimem, gprimem, pmlcheck, kx, kz)
                    if (inv_key=='forward') then
                        if (mesh%has_src) then
                            do i=1,src%nsrc
                                if (act_src(src%srcnr(i))) then
                                    if (src%srcelem(i)==ie) then      ! add source term
                                        if (src%srctype(i) == 0) then
                                            !single force
                                            rq(iv,4:5)=rq(iv,4:5)+srcArray(:,1:2,i)
                                            
                                        elseif (src%srctype(i) == 1) then
                                            rq(iv,1:3)=rq(iv,1:3)+srcArrayM(:,1:3,i)
                                        endif
                                    endif
                                endif
                            enddo
                        endif
                    else
                        if (mesh%has_rec) then
                            do i=1,rec%nrec
                                if (rec%recelem(i)==ie) then      ! add source term
                                    rq(iv,4:5)=rq(iv,4:5)+srcArray(:,1:2,i)
                                endif
                            enddo
                        endif
                    endif

                    ! do runge kutta
                    ! RK 4
                    if (wrk==1) then
                        if (irk==1) then
							!yn
                            qn(iv,:)=q(iv,:)
                            if(par%set_pml) then
                                fprimen(iv,:)=fprime(iv,:)
                                gprimen(iv,:)=gprime(iv,:)
                            end if
							!k1
							resQ(iv, :) = dt*rq(iv,:)
							!p1
							q(iv, :) = q(iv, :) + onehalf*resQ(iv, :)
							!pml
							if (pmlcheck(ie)) then
								do i=1,dimens
									resFprime(iv,i) = dt*(-alphax(iv)*fprime(iv,i) - ddx(iv)/kx(iv)*(fprime(iv,i)+ftemp(:,i)))
									resGprime(iv,i) = dt*(-alphaz(iv)*gprime(iv,i) - ddz(iv)/kz(iv)*(gprime(iv,i)+gtemp(:,i)))
									fprime(iv,i)    = fprime(iv,i) + onehalf*resFprime(iv,i)
									gprime(iv,i)    = gprime(iv,i) + onehalf*resGprime(iv,i)
								end do
							end if
							
                            rk_time=0.5
                        end if !rk==1
						
                        if (irk==2) then
							!k2
							resQ2(iv, :) = dt*rq(iv,:)
							!p2
							q(iv, :) = qn(iv, :) + onehalf*resQ2(iv, :)
							!pml
							if (pmlcheck(ie)) then
								do i=1,dimens
									resFprime2(iv,i) = dt*(-alphax(iv)*fprime(iv,i) - ddx(iv)/kx(iv)*(fprime(iv,i)+ftemp(:,i)))
									resGprime2(iv,i) = dt*(-alphaz(iv)*gprime(iv,i) - ddz(iv)/kz(iv)*(gprime(iv,i)+gtemp(:,i)))
									fprime(iv,i)    = fprime(iv,i) + onehalf*resFprime2(iv,i)
									gprime(iv,i)    = gprime(iv,i) + onehalf*resGprime2(iv,i)
								end do
							end if
                            rk_time=0.5
                        end if !rk==2
                        if (irk==3) then
							!k3
							resQ3(iv, :) = dt*rq(iv,:)
							!p3
							q(iv, :) = qn(iv, :) + resQ3(iv, :)
							!pml
							if (pmlcheck(ie)) then
								do i=1,dimens
									resFprime3(iv,i) = dt*(-alphax(iv)*fprime(iv,i) - ddx(iv)/kx(iv)*(fprime(iv,i)+ftemp(:,i)))
									resGprime3(iv,i) = dt*(-alphaz(iv)*gprime(iv,i) - ddz(iv)/kz(iv)*(gprime(iv,i)+gtemp(:,i)))
									fprime(iv,i)    = fprime(iv,i) + resFprime3(iv,i)
									gprime(iv,i)    = gprime(iv,i) + resGprime3(iv,i)
								end do
							end if
							rk_time=1.0
                        end if !rk==3
						if (irk==4) then
							!k4
							resQ4(iv, :) = dt*rq(iv,:)
							!p4
							q(iv, :) = qn(iv, :) + onesix*(resQ(iv, :)+2.0*resQ2(iv, :)+2.0*resQ3(iv, :)+resQ4(iv, :))
							!pml
							if (pmlcheck(ie)) then
								do i=1,dimens
									resFprime4(iv,i) = dt*(-alphax(iv)*fprime(iv,i) - ddx(iv)/kx(iv)*(fprime(iv,i)+ftemp(:,i)))
									resGprime4(iv,i) = dt*(-alphaz(iv)*gprime(iv,i) - ddz(iv)/kz(iv)*(gprime(iv,i)+gtemp(:,i)))
									fprime(iv,i)    = fprime(iv,i) + onesix*(resFprime(iv,i)+2.0*resFprime2(iv,i)+2.0*resFprime3(iv,i)+resFprime4(iv,i))
									gprime(iv,i)    = gprime(iv,i) + onesix*(resGprime(iv,i)+2.0*resGprime2(iv,i)+2.0*resGprime3(iv,i)+resGprime4(iv,i))
								end do
							end if
                        end if !rk==4
						!tvd-rk3
                    else if (wrk==2) then
                        if (irk==1) then
                            qn(iv,:)=q(iv,:)
                            if(par%set_pml) then
                                fprimen(iv,:)=fprime(iv,:)
                                gprimen(iv,:)=gprime(iv,:)
                            end if
                            q(iv,:) = q(iv,:) + dt*rq(iv,:)
                            !PML
                            if (pmlcheck(ie)) then
                                do i=1,dimens
                                    fprime(iv,i) = fprime(iv,i) + dt*(-alphax(iv)*fprime(iv,i) - ddx(iv)/kx(iv)*(fprime(iv,i)+ftemp(:,i)))
                                    gprime(iv,i) = gprime(iv,i) + dt*(-alphaz(iv)*gprime(iv,i) - ddz(iv)/kz(iv)*(gprime(iv,i)+gtemp(:,i)))
                                end do
                            end if
                            rk_time=1.
                        end if !rk==1
                        if (irk==2) then
                            q(iv,:) = threefor*qn(iv,:) + onefor*(q(iv,:) + dt*rq(iv,:))
                            !PML
                            if (pmlcheck(ie)) then
                                do i=1,dimens
                                    fprime(iv,i) = threefor*fprimen(iv,i) + onefor*(fprime(iv,i) + dt*(-alphax(iv)*fprime(iv,i) - ddx(iv)/kx(iv)*(fprime(iv,i)+ftemp(:,i))))
                                    gprime(iv,i) = threefor*gprimen(iv,i) + onefor*(gprime(iv,i) + dt*(-alphaz(iv)*gprime(iv,i) - ddz(iv)/kz(iv)*(gprime(iv,i)+gtemp(:,i))))
                                end do
                            end if
                            rk_time=0.5
                        end if !rk==2
                        if (irk==3) then
                            q(iv,:) = onethree*qn(iv,:) + twothree*(q(iv,:) + dt*rq(iv,:))
                            !PML
                            if (pmlcheck(ie)) then
                                do i=1,dimens
                                    fprime(iv,i) = onethree*fprimen(iv,i) + twothree*(fprime(iv,i) + dt*(-alphax(iv)*fprime(iv,i) - ddx(iv)/kx(iv)*(fprime(iv,i)+ftemp(:,i))))
                                    gprime(iv,i) = onethree*gprimen(iv,i) + twothree*(gprime(iv,i) + dt*(-alphaz(iv)*gprime(iv,i) - ddz(iv)/kz(iv)*(gprime(iv,i)+gtemp(:,i))))
                                end do
                            end if
                        end if !rk==3
                    ! les-rk4
                    else if (wrk == 3) then
                        resQ(iv, :) = rk4a(irk)*resQ(iv, :) + dt*rq(iv,:)
                        q(iv, :) = q(iv, :) + rk4b(irk)*resQ(iv, :)
                        !pml
                        if (pmlcheck(ie)) then
                            do i=1,dimens
                                resFprime(iv,i) = rk4a(irk)*resFprime(iv,i) + dt*(-alphax(iv)*fprime(iv,i) - ddx(iv)/kx(iv)*(fprime(iv,i)+ftemp(:,i)))
                                resGprime(iv,i) = rk4a(irk)*resGprime(iv,i) + dt*(-alphaz(iv)*gprime(iv,i) - ddz(iv)/kz(iv)*(gprime(iv,i)+gtemp(:,i)))
                                fprime(iv,i)    = fprime(iv,i) + rk4b(irk)*resFprime(iv,i)
                                gprime(iv,i)    = gprime(iv,i) + rk4b(irk)*resGprime(iv,i)
                            end do
                        end if
                        if (irk < nrk) rk_time=rk4c(irk+1)
					end if
                    if (irk == nrk) then
                        !integrate fields with euler, for plotting
                        !This is done inside of the rk-loop so it is done for each iteration of the rk-loop.
                        !In case of rk2 it is done twice and in case of rk4 it would be done 4 times.
                        !Displacement
                        ux(iv) = ux(iv) + dt*q(iv,4)
                        uz(iv) = uz(iv) + dt*q(iv,5)
                        uplot(iv) = sqrt(ux(iv)**2+uz(iv)**2)
                        if (movie%movie) then
                            if (movie%velocity) then
                                !Calculate the norm of the velocity
                                vplot(iv) = sqrt(q(iv,4)**2+q(iv,5)**2)
                            endif
                        endif

                        !Acceleration
                        ax(iv) = rq(iv,4)
                        az(iv) = rq(iv,5)
                        ! Strain
                        a_sigma_tmp(:,1)=q(iv,1)
                        a_sigma_tmp(:,2)=q(iv,2)
                        a_sigma_tmp(:,3)=q(iv,3)
                        e(iv,1) = (-mesh%lambda(ie) * a_sigma_tmp(:,2) + (2 * mesh%mu(ie) + mesh%lambda(ie)) * a_sigma_tmp(:,1))/ (4*mesh%mu(ie)**2 + 4 * mesh%lambda(ie)*mesh%mu(ie))
                        e(iv,2) = ((2 * mesh%mu(ie) + mesh%lambda(ie)) * a_sigma_tmp(:,2) - mesh%lambda(ie)*a_sigma_tmp(:,1))/(4*mesh%mu(ie)**2 + 4*mesh%lambda(ie)*mesh%mu(ie))
                        e(iv,3) = a_sigma_tmp(:,3)/(2 * mesh%mu(ie))
                    endif
                end do !element loop
                call sync_mpi()
            end do ! RK-LOOP

            do ie=1,mesh%nelem
                iv=mesh%ibool(:,ie)
                do i=1,Np
                    energy_kin(it) = energy_kin(it) + mesh%rho(ie) * (q(iv(i),4)**2 + q(iv(i),5)**2)
                    energy_pot(it) = energy_pot(it) + 0.5 * (q(iv(i),1)*e(iv(i),1) + q(iv(i),2)*e(iv(i),2) + q(iv(i),3)*e(iv(i),3))
                end do
            end do
            energy(it) = energy_kin(it) + energy_pot(it)
            call sum_real_all(energy(it),all_energy(it), CUSTOM_REAL)

            ! check if energy grows to much, its done with sta/lta trigger to switch off the pml if necessary
            if(par%use_trigger) then
                if(par%set_pml) then
                    if(pmlcrit) then
                        if( (it>timecrit) ) then
                            avg_energy1=0
                            do i=1,par%avg_window1
                                avg_energy1=avg_energy1+all_energy(it-i)
                            end do
                            avg_energy1=avg_energy1/par%avg_window1

                            avg_energy2=0
                            do i=1,par%avg_window2
                                avg_energy2=avg_energy2+all_energy(it-i)
                            end do
                            avg_energy2=avg_energy2/par%avg_window2

                            sta_lta=avg_energy2/avg_energy1

                            if( abs(sta_lta -1) > par%sta_lta_trigger) then
                                write(*,*)
                                write(*,*) " !!!!!!!!!! warning, energy grows to much and pml is switched off !!!!!!"
                                write(*,*)
                                pmlcrit=.false.
                                pmlcheck(:)=.false.
                            end if
                        end if
                    end if
                end if
            end if
			if (inv_key=='forward') then
                if (mod(it,movie%frame) == 0 .or. it == par%nt) then
                    call maxval_real(maxval(uplot),maxu,CUSTOM_REAL)
                    call maxval_real(maxval(vplot),maxv,CUSTOM_REAL)
                    if (movie%movie) then
                        if (movie%velocity) then
                            call writePlotBinaries(vplot, outpath,"moviedata_normV", myrank, it)
                            call writePlotBinaries(q(:,4),outpath, "moviedata_vx", myrank, it)
                            call writePlotBinaries(q(:,5), outpath,"moviedata_vz", myrank, it)
                        end if
                        if (movie%stress) then
                            call writePlotBinaries(q(:,1), outpath,"moviedata_sigmaXX", myrank, it)
                            call writePlotBinaries(q(:,2), outpath,"moviedata_sigmaZZ", myrank, it)
                            call writePlotBinaries(q(:,3), outpath,"moviedata_sigmaXZ", myrank, it)
                        end if
                        if (movie%displacement) then
                            call writePlotBinaries(uplot, outpath,"moviedata_normU", myrank, it)
                            call writePlotBinaries(ux, outpath,"moviedata_ux", myrank, it)
                            call writePlotBinaries(uz, outpath,"moviedata_uz", myrank, it)
                        end if
                    endif
                    if (myrank == 0) then
                        write(filename,"('energytemp')")
                        call plotPoints2D(stf(:),all_energy,trim(outpath)//filename)
                        call outputstamp(par, localtime, it, timestamp, maxu, maxv)
                    end if
                endif
            else
                if (inv_key=='forward') then
                    if (mod(it,movie%frame) == 0 .or. it == par%nt) then
                        write(srcstring,"('_src',i3.3)") i_src
                        call maxval_real(maxval(uplot),maxu,CUSTOM_REAL)
                        call maxval_real(maxval(vplot),maxv,CUSTOM_REAL)
                        if (run_number==1) then
                            if (movie%movie) then
                                if (movie%velocity) then
                                    call writePlotBinaries(vplot, outpath,"moviedata_normV"//srcstring, myrank, it)
                                    call writePlotBinaries(q(:,4),outpath, "moviedata_vx"//srcstring, myrank, it)
                                    call writePlotBinaries(q(:,5), outpath,"moviedata_vz"//srcstring, myrank, it)
                                end if
                                if (movie%stress) then
                                    call writePlotBinaries(q(:,1), outpath,"moviedata_sigmaXX"//srcstring, myrank, it)
                                    call writePlotBinaries(q(:,2), outpath,"moviedata_sigmaZZ"//srcstring, myrank, it)
                                    call writePlotBinaries(q(:,3), outpath,"moviedata_sigmaXZ"//srcstring, myrank, it)
                                end if
                                if (movie%displacement) then
                                    call writePlotBinaries(uplot, outpath,"moviedata_normU"//srcstring, myrank, it)
                                    call writePlotBinaries(ux, outpath,"moviedata_ux"//srcstring, myrank, it)
                                    call writePlotBinaries(uz, outpath,"moviedata_uz"//srcstring, myrank, it)
                                end if
                            endif
                        endif
                    endif
                else
                    if (mod(it+movie%frame-time_shift_movie,movie%frame) == 0) then
                        write(srcstring,"('_src',i3.3)") i_src
                        call maxval_real(maxval(uplot),maxu,CUSTOM_REAL)
                        call maxval_real(maxval(vplot),maxv,CUSTOM_REAL)
                        if (movie%movie) then
                            if (movie%velocity) then
                                call writePlotBinaries(vplot, outpath,"moviedata_normV_adj"//srcstring, myrank, par%nt-it+movie%frame)
                                call writePlotBinaries(q(:,4),outpath,"moviedata_vx_adj"//srcstring, myrank, par%nt-it+movie%frame)
                                call writePlotBinaries(q(:,5),outpath,"moviedata_vz_adj"//srcstring, myrank, par%nt-it+movie%frame)
                            end if
                            if (movie%stress) then
                                call writePlotBinaries(q(:,1),outpath,"moviedata_sigmaXX_adj"//srcstring, myrank, par%nt-it+movie%frame)
                                call writePlotBinaries(q(:,2),outpath, "moviedata_sigmaZZ_adj"//srcstring, myrank, par%nt-it+movie%frame)
                                call writePlotBinaries(q(:,3),outpath, "moviedata_sigmaXZ_adj"//srcstring, myrank, par%nt-it+movie%frame)
                            end if
                            if (movie%displacement) then
                                call writePlotBinaries(uplot,outpath,"moviedata_normU_adj"//srcstring, myrank, par%nt-it+movie%frame)
                                call writePlotBinaries(ux, outpath,"moviedata_ux_adj"//srcstring, myrank, par%nt-it+movie%frame)
                                call writePlotBinaries(uz, outpath,"moviedata_uz_adj"//srcstring, myrank, par%nt-it+movie%frame)
                            end if
                        endif
                    endif
                endif
            endif
        end do! timeloop

        ! save seismos
        if (mesh%has_rec) then
            do r=1,rec%nrec
                write(filename,"('seismo.x.',i7.7,'.sdu')") rec%recnr(r)
                call plotPoints2D(stf(par%subsampling_factor::par%subsampling_factor)-par%plott0,plotux(r,:),trim(outpath)//filename)
                write(filename,"('seismo.z.',i7.7,'.sdu')") rec%recnr(r)
                call plotPoints2D(stf(par%subsampling_factor::par%subsampling_factor)-par%plott0,plotuz(r,:),trim(outpath)//filename)
                write(filename,"('seismo.x.',i7.7,'.sdv')") rec%recnr(r)
                call plotPoints2D(stf(par%subsampling_factor::par%subsampling_factor)-par%plott0,plotv(r,:),trim(outpath)//filename)
                write(filename,"('seismo.z.',i7.7,'.sdv')") rec%recnr(r)
                call plotPoints2D(stf(par%subsampling_factor::par%subsampling_factor)-par%plott0,plotw(r,:),trim(outpath)//filename)
                write(filename,"('seismo.x.',i7.7,'.sda')") rec%recnr(r)
                call plotPoints2D(stf(par%subsampling_factor::par%subsampling_factor)-par%plott0,plotax(r,:),trim(outpath)//filename)
                write(filename,"('seismo.z.',i7.7,'.sda')") rec%recnr(r)
                call plotPoints2D(stf(par%subsampling_factor::par%subsampling_factor)-par%plott0,plotaz(r,:),trim(outpath)//filename)
                if (par%div) then
                    write(filename,"('seismo.r.',i7.7,'.sdv')") rec%recnr(r)
                    call plotPoints2D(stf(par%subsampling_factor::par%subsampling_factor)-par%plott0,plot_r(r,:),trim(outpath)//filename)
                end if
                if (par%curl) then
                    write(filename,"('seismo.t.',i7.7,'.sdv')") rec%recnr(r)
                    call plotPoints2D(stf(par%subsampling_factor::par%subsampling_factor)-par%plott0,plot_t(r,:),trim(outpath)//filename)
                end if
                do i=1,size(act_src)
                    if (act_src(i)) i_src=i
                enddo
            end do
        end if

        if (myrank==0) then
            write(filename,"('energy')")
            call plotPoints2D(stf(:)-par%plott0,all_energy,trim(outpath)//filename)
        end if

		if (par%log.and.myrank== 0) then
            write(*,"(a80)") "|                              timeloop finished                               |"
            write(*,"(a80)") "|------------------------------------------------------------------------------|"
        end if
        close(27)

        deallocate(ux,uz,ax,az,uplot,vplot,v1plot,v2plot)
        deallocate(rQ,Q,Qn,Qm,e,resQ, resU, u)
		deallocate(resQ2,resQ3,resQ4)
        deallocate(energy, all_energy, energy_kin, energy_pot, a_sigma_tmp)
        deallocate(stf)
        deallocate(elasticfluxvar, APA, T, invT, VT, VTfree)

        deallocate(srcelemV)
        deallocate(recelemV)

        if (mesh%has_src) then
            deallocate(srcInt,srcTemp)
            deallocate(t0)
            if (allocated(srcArray)) deallocate(srcArray)
            deallocate(srcArrayM)
            deallocate(plotstf)
            deallocate(plotDiffstf)
        end if
        if (mesh%has_rec) then
            deallocate(recInt,recTemp)
            deallocate(plotv,plotw,plotux,plotuz,plotax,plotaz)
        end if

        deallocate(q_send)
        deallocate(q_rec)
        deallocate(qi)
        deallocate(qi_test)
        deallocate(req1)
        deallocate(pmlcheck)

        if(par%set_pml) then
            deallocate(fprime)
            deallocate(gprime)
            deallocate(fprimen,fprimem)
            deallocate(gprimen,gprimem)
            deallocate(ddx, ddz)
            deallocate(alphax, alphaz)
            deallocate(kx,kz)
            deallocate(np_zeros)
            deallocate(resFprime, resGprime)
        end if
        call deallocMeshvar(mesh)

    end subroutine timeloop2d
end module timeloopMod
