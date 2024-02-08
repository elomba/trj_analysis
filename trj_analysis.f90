program trj_analysis
    !
    !  This program performs an analysis of a netcdf trajectory file from LAMMPS.
    !  At this stage it computes:
    !     - Pair distribution functions
    !     - Static structure factors
    !     - Cluster analysis (average cluster profiles, cluster-cluster rdf's and sq's)
    !                         size and radii distributions, and a trajectory file with the
    !                         evolution of cluster  com's 
    !     - Dynamics (position and velocity correlation functions)
    !   The input is provided as a set of namelist (see attached example) 
    !
    !   Programmed in NVIDIA CUDA Fortran
    !
    !   A. Diaz-Pozuelo & E. Lomba, Madrid February 2024 
    !
    !
    use mod_precision
    use mod_common
    use mod_input
    use mod_nc
    use mod_nc_conf
    use mod_cells
    use mod_clusters
    use mod_sq
    use mod_rdf
    use mod_densprof
    use mod_log
    use mod_thermo
    use mod_dyn, only : dyn_init, dyn_clear, rtcorr, print_rtcor
    use cudafor
    implicit none

    integer :: io = 0, ioerr, istat, ncid_in
    integer :: argc, ncstart, i
    logical :: selectall = .false., clnfound = .true. 
    real :: t0 = 0, t1 = 0

    call cpu_time(time_cpu_start)
    call cpu_time(cpu0)

    ! Get CUDA properties
    istat = cudaSetDevice(0)
    istat = cudaGetDeviceProperties(gpu_properties, 0)
    shmsize = gpu_properties%sharedMemPerBlock
    maxthread = gpu_properties%maxThreadsPerBlock
    istat = cudaEventCreate(startEvent)
    istat = cudaEventCreate(stopEvent)

    ! Command line arguments control
    argc = command_argument_count()
    if (argc /= 1) then
        stop 'analisys: You must specify only ONE argument: the name of the &
                & input file\r\nexample: analisys.exe input_file.nml'
    end if
    call get_command_argument(1, input_filename)
    ! Load JSON input file & init log system
    call read_input_file()
    call log_init()
    ! Check that maximum number of threads is not surpassed
    if (nthread > maxthread/8) then
        nthread = maxthread/8
        write(*,'("** Warning: number of threads reset to",I3)')nthread
    endif
    ! Selected specie's types pre netcdf load
    if (sp_types_selected(1) == 0) then
        selectall = .true.
    else
        wtypes = sp_types_selected
    end if

    ! First load from netcdf input file. Load global attributes
    call check(nf90_open(path=trj_input_file, mode=NF90_WRITE, ncid=ncid_in), ioerr)
    call read_nc_cfg(ncid_in, 1, io, io_log_file)

    ! Set number of species from netcdf file or from selected species from JSON
    if (selectall) then
        nsp = ntypes ! from netcdf
    else
        nsp = size(sp_types_selected) ! from JSON
    end if
    ! If number of species is not equal to size of list's properties, print error
    if (nsp /= size(sp_labels) .and. nsp /= size(mat) .and. nsp /= size(bsc)) then
        print *, 'ERROR: nsp not equal to size of sp_labels/sp_atomic_weight/sp_scattering'
        print *, 'Check your JSON input file!'
        stop
    end if
    nit = nsp*(nsp + 1)/2
    ! Set number of configurations to read from netcdf & from and to number of configurations
    if (ncfs_from_to(1) == 0) then
        ncfs_from_to(1) = nconf_i
        ncfs_from_to(2) = 1
        ncfs_from_to(3) = nconf_i
    end if
    nconf = ncfs_from_to(1)
    ncfs_from_to(3) = ncfs_from_to(3) + 1

    ! Keytrj set & control (0: only positions ; 1: positions & velocities)
    if (allocated(v)) then
        keytrj = 1
    else
        keytrj = 0
    end if

    ! Modules to run
    clnfound = .true.
    if (rdf_sq_cl_dyn_conf(1) == .true.) run_rdf = .true.
    run_sq = .false. 
    if (rdf_sq_cl_dyn_conf(2) == .true.) run_sq = .true.
    if (rdf_sq_cl_dyn_conf(3) == .true.) then
        !
        ! Cluster analysis needs rdf's and s(q)'s to be computed
        ! This is also modified in input.f90
        !
        run_clusters = .true.
        run_rdf = .true.
        run_sq = .true.
        clnfound = .false.
    end if
    if (rdf_sq_cl_dyn_conf(4) == .true.) run_dyn = .true. 
   ! Deactivate use of cell lists if clusters not analysed 
    if (clnfound) then
        rcl = -1.
        use_cell = .false.
    endif

    ! Init common variables & print log_01
    call common_init(natoms, ndim, nthread, idir, conf(4)%units,conf(4)%scale, nsp)
    ! call log_01()

    if (idir > 0) then
        call prof_init()
    end if
    ! Simulation begin, every configuration selected
    do i = 1, ncfs_from_to(1)
        ! In the first configuration init modules
        if (i == 1) then
            if (use_cell) call cells_init_pre_nc_read(nmol)
            if (run_clusters) call clusters_init(nmol)
            if (run_dyn) call dyn_init()
        end if

        ! Read i configuration from netcdf input file
        call cpu_time(t0)
        ncstart = ncfs_from_to(2) + (i - 1)*(ncfs_from_to(3) - ncfs_from_to(2))/ncfs_from_to(1)
        write (io_log_file, '(/,A,I0)'), 'Reading configuration number: ', ncstart
        call read_nc_cfg(ncid_in, ncstart, io, io_log_file)
        ! Pre configuration analysis data transformations, corrections and format
        call trans_ncdfinput()
        call cpu_time(t1)
        tread = tread + t1 - t0

        ! In each configuration run modules
        if (run_clusters) then
            if (use_cell) then
                if (i == 1) call cells_init_post_nc_read()
                call cells_build()
            end if
        end if
        if (use_cell) call cells_reset_struct() ! Put inside above if ?? if use_cell == false, it's run
        if (run_rdf) then
            if (i==1) call RDF_init(nsp)
        endif
        if (run_sq) then
            if (i == 1) then
                call sq_init(nmol, nsp, nbcuda)
                if (run_clusters) call clusters_sq_init()
            end if
        end if
        !done(:) = 0 !!! WARN, try to move inside function

        ! Transfer data to GPU
        call transfer_cpu_gpu(ndim)

        ! BFS cluster search
        if (run_clusters) then
            call cluster_search()
        end if
        ! Thermodynamics calculus
        call thermo_calc(i)

        ! Run RDF
        if (run_rdf) call RDFcomp(Nmol, i, nbcuda, nthread)

        ! Compute density profile along idir direction
        if (idir > 0) then
            call profile_comp(nthread, ndim, idir, pwall, deltar)
        end if

        ! Compute SQ
        if (run_sq) call SQcalc()
        ! Compute cluster properties
        if (run_clusters) call cluster_analysis(i)
        ! Compute dynamics
        if (run_dyn) then
            call rtcorr(i)
        end if 

        ! Print periodic output
        call print_output(i)
    end do

    ! Normalize density profiles computed along the non-periodic dimension
    if (idir > 0) then
        call normdenspr(nconf)
    end if

    ! Print out S(Q)'s
    if (run_sq) then
        call sq_transfer_gpu_cpu()
        call printSQ(Nmol)
    end if

    do i = 1, nsp
        write (*, '(" ** ",i6," atoms of type ",i2)') ntype(i), i
    end do

    ! Print partial pair distribution functions
    if (run_rdf) then
         call printrdf(rcl, lsmax)
    end if
    if (run_dyn) then
        call print_rtcor()
    endif
    if (run_clusters) then
        call print_clusinfo(nqmin, Nmol)
    end if

    ! Cleaning house
    if (run_sq) call sq_clear()
    if (run_rdf) call rdf_clear()
    if (run_clusters) call clusters_clear()
    if (use_cell) call cells_clear()
    if (run_dyn) call dyn_clear()
    if (idir > 0) call prof_clear()
    call common_clear()
    call log_clear()
    call input_clear()
    call cpu_time(time_cpu_stop)

    ! Print simulation time
    time_total = time_total + (time_cpu_stop - time_cpu_start)
    write (*, '(/,A,F15.7,A,/)') '**** Total time: ', time_total, ' seconds'
end program trj_analysis
