program trj_analysis
    !
    !  This program performs an analysis of a netcdf trajectory file from LAMMPS.
    !  At this stage it computes:
    !     - Pair distribution functions
    !     - Static structure factors
    !     - Cluster analysis (average cluster profiles, cluster-cluster rdf's and sq's)
    !                         size and radii distributions, and a trajectory file with the
    !                         evolution of cluster  com's 
    !     - Dynamics (position and velocity correlation functions, 
    !                 intermediate scattering functions, frequency dependent S(q,w)'s, complete and self)
    !                 Only total quantities are computed in mixtures. Use the selection option if
    !                 single component quantities are needed
    !     - Kinetic energy (if velocities present in the trajectory file)
    !     - Potential energ (if compute pe/peratom pressent in the trajectory file)
    !     - Stress tensor (pressure)  (if compute stress/atom pressent in the trajectory file)
    !       LAMMPS table file with RSQ tabulation -IMPORTANT!!-)
    !     - The code allows for selection of specific components in mixtures
    !   The input is provided as a set of namelist (see attached example) 
    !
    !   Restrictions: Only orthogonal cells are allowed, particle numbers must remain constant all along
    !                 the simulation
    !
    !   Programmed in NVIDIA CUDA Fortran
    !
    !   A. Diaz-Pozuelo & E. Lomba, Madrid/Santiago de Compostela, fall 2024 
    !
    !   Important notice: in LAMMPS script the following computes must be included in the dump
    !                     in order to compute potential energies and pressures
    !             
    !       compute stress all stress/atom NULL
    !       compute ener all pe/atom
    !       dump trj1 all netcdf ${Ndump} run.nc  id type x y z vx vy vz c_stress[*] c_ener
    !    In the first INPUT namelist optional character variables "ener_name" and "press_name" refer to the
    !    names of the computes
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

    ! Get CUDA properties from device 0 (can be set from environmente variable CUDA_VISIBLE_DEVICES)

    istat = cudaSetDevice(0)
    if (istat == 0) then
        istat = cudaGetDeviceProperties(gpu_properties, 0)
    else
        write(*,"('*** Unrecoverable error: no GPU available !!')")
        stop
    end if
    shmsize = gpu_properties%sharedMemPerBlock
    maxthread = gpu_properties%maxThreadsPerBlock
    ! Initialize CUDA timing
    istat = cudaEventCreate(startEvent)
    istat = cudaEventCreate(stopEvent)
    ! Print program header
    write(*,"(/80('*')/'*',78(' '),'*')")
    write(*,"('*    Program trj_analysis: analyzing LAMMPS trajectory in NETCDF format',t80,'*')")
    write(*,"('*',t80,'*')")
    write(*,"('*    Using GPU with CUDA nvfortran/nvcc >= 11.6',t80,'*')")
    write(*,"('*',t80,'*')")
    write(*,"('*    Version 0.2.30 October 2024',,t80,'*')")
    write(*,"('*',78(' '),'*'/80('*')/)")
    call printDevPropShort(gpu_properties, 0)
    ! Command line arguments control
    argc = command_argument_count()
    if (argc /= 1) then
        stop '!!! Error: You must specify only ONE argument: the name of the &
                & input file\r\nexample: trj_analisys.exe input_file.nml'
    end if
    call get_command_argument(1, input_filename)
    ! Load namelist input file & init log system
    call read_input_file()
    call log_init()
    ! Check that maximum number of threads is not surpassed
    if (nthread > maxthread/8) then
        nthread = maxthread/8
        write(*,'("** Warning: number of threads reset to",I3)')nthread
    endif

    ! First load from netcdf input file. Load global attributes
    call check(nf90_open(path=trj_input_file, mode=NF90_WRITE, ncid=ncid_in), ioerr)
    if (ioerr .ne. 0) then
        stop("** UNRECOVERABLE ERROR: cannot open NetCDF trajectory file !")
    endif 
    !
    ! Read header and details of the NETCDF trajectory files: check consistency with input data
    !
    call read_nc_cfg(ncid_in, 1, io, io_log_file)
    ! Set number of species from netcdf file or from selected species from namelist
    if (nsp > ntypes) then 
        print *, ' ERROR: number of species in input file is',nsp,' larger than that in netcdf file ',ntypes
        STOP
    else if (nsp == ntypes) Then
        selectall = .true.
        nmol = natoms
    else if (nsp < ntypes) then 
        allocate(wtypes(nsp))
        wtypes = sp_types_selected        
        call reset_nmol(nmol)
    end if 
    ! Number of different interactions
    nit = nsp*(nsp + 1)/2
    ! Set number of configurations to read from netcdf & from and to number of configurations
    if (ncfs_from_to(1) == 0) then
        ncfs_from_to(1) = nconf_i
        ncfs_from_to(2) = 1
        ncfs_from_to(3) = nconf_i
    end if
    nconf = ncfs_from_to(1)
    ncfs_from_to(3) = ncfs_from_to(3) + 1

    ! Modules to run
    clnfound = .true.
    ! Radial ditribution function
    if (rdf_sq_cl_dyn_sqw_conf(1) == .true.) run_rdf = .true.
    run_sq = .false. 
    ! Static structure factors
    if (rdf_sq_cl_dyn_sqw_conf(2) == .true.) run_sq = .true.
    ! Cluster analysis
    if (rdf_sq_cl_dyn_sqw_conf(3) == .true.) then
        !
        ! Cluster analysis needs rdf's and s(q)'s to be computed
        ! This is also modified in input.f90
        !
        run_clusters = .true.
        run_rdf = .true.
        run_sq = .true.
        clnfound = .false.
    end if
    ! Dynamic correlations
    if (rdf_sq_cl_dyn_sqw_conf(4) == .true.) run_dyn = .true. 
   ! Deactivate use of cell lists if clusters not analysed 
    if (clnfound) then
        rcl = -1.
        use_cell = .false.
    endif
    ! Activate analysis of dynamic structure factor 
    if (rdf_sq_cl_dyn_sqw_conf(5) == .true.) then
        run_sqw = .true.
        run_dyn = .true.
        run_sq = .true.
    endif
    ! Confined system (density profile analysis in the direction of confinement)
    if (rdf_sq_cl_dyn_sqw_conf(6) == .true.) run_rdf = .true. 

    ! Init common variables & print outs
    call common_init(nmol, ndim, nthread, idir, conf(4)%units,conf(4)%scale, nsp)
    
    ! Initialize profiles : idir defines the direction of confinement (x,y,z->1,2,3)
    if (idir > 0) then
        call prof_init()
    end if
    ! Analysis begins every configuration selected
    do i = 1, ncfs_from_to(1)
        ! In the first configuration init modules
        if (i == 1) then
            if (use_cell) call cells_init_pre_nc_read(nmol)
            if (run_clusters) call clusters_init(nmol)
            if (run_dyn) call dyn_init()
        end if
        ! Read i configuration from netcdf input file
        call cpu_time(t0)
        ! Jumps configurations to be read
        ncstart = ncfs_from_to(2) + (i - 1)*(ncfs_from_to(3) - ncfs_from_to(2))/ncfs_from_to(1)
        ! Read-in a full configuration
        call read_nc_cfg(ncid_in, ncstart, io, io_log_file)
        ! Correct of end of file reached
        if (io<0) then
            ncfs_from_to(1)=i-1
            Exit
        endif 
        ! Pre configuration analysis data transformations, corrections and format
        if (ntypes == nsp) then
            ! All species selected
            call trans_ncdfinput()
        else
            ! Only some species selected from the configuration
            call select_ncdfinput()
        endif
        call cpu_time(t1)
        tread = tread + t1 - t0
        ! Over each configuration run selected modules
        if (run_clusters) then
            ! Linked cells for cluster analysis
            if (use_cell) then
                if (i == 1) call cells_init_post_nc_read()
                call cells_build()
            end if
        end if
        if (use_cell) call cells_reset_struct() ! Put inside above if ?? if use_cell == false, it's run
        !
        if (i==1) then
            if (run_rdf) call RDF_init(nsp)
            if (run_sq) call sq_init(nmol, nsp, nbcuda)
            if (run_clusters) call clusters_sq_init()
        end if

        ! Transfer data to GPU
        call transfer_cpu_gpu(ndim)

        ! BFS cluster search
        if (run_clusters) call cluster_search()
        
        ! kinetic energy (if velocities available) 
        if (ex_vel) call thermo_kin(i, ndim)
        ! Run RDF
        if (run_rdf) call RDFcomp(Nmol, i, nbcuda, nthread)

        ! Compute density profile along idir direction
        if (idir > 0) then
            call profile_comp(nthread, ndim, idir, pwall, deltar)
        end if
        ! Compute SQ
        if (run_sq) call SQcalc()
        !
        ! Compute cluster properties
        if (run_clusters) call cluster_analysis(i)
        
        ! Compute potential energy
        if (run_thermo) call poteng(i, natoms)

        ! Compute pressure
        if (ex_stress) call stress_calc(i, ndim, natoms)
     
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
        if (run_thermo) then
            call printPotEngCl()
            call printPotEngClCl()
        end if
    end if

    ! Cleaning house

    if (run_sq) call sq_clear()
    if (run_rdf) call rdf_clear()
    if (run_clusters) call clusters_clear()
    if (use_cell) call cells_clear()
    if (run_dyn) call dyn_clear()
    if (idir > 0) call prof_clear()
    if (run_thermo) call thermo_clear()
    call common_clear()
    call log_clear()
    call input_clear()
    call cpu_time(time_cpu_stop)

    ! Print simulation time
    time_total = time_total + (time_cpu_stop - time_cpu_start)
    write (*, '(/,A,F15.7,A,/)') '**** Total time: ', time_total, ' seconds'
end program trj_analysis
