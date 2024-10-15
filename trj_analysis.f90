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
   ! use mod_densprof
    use mod_log
    use mod_thermo
    use mod_dyn, only : dyn_init, dyn_clear, rtcorr, print_rtcor
    use mod_util, only : gpu_and_header, clean_memory, init_modules, reformat_input_conf, &
                         basic_init, print_results
    use cudafor
    implicit none

    integer :: io = 0, ioerr, istat, ncid_in
    integer :: argc, ncstart, i
    logical :: selectall = .false., clnfound = .true., first_configuration=.true. 
    real :: t0 = 0, t1 = 0
    ! Command line arguments control
    argc = command_argument_count()
    if (argc /= 1) then
        stop '!!! Error: You must specify only ONE argument: the name of the &
                & input file\r\nexample: trj_analisys.exe input_file.nml'
    end if
    call get_command_argument(1, input_filename)
    !
    ! initialize timers
    !
    call cpu_time(time_cpu_start)
    call cpu_time(cpu0)

    ! Get CUDA properties from device 0 (can be set from environmente variable CUDA_VISIBLE_DEVICES)
    call gpu_and_header(startEvent,stopEvent)
   
  
    ! Load namelist input file & init log system
    call read_input_file()
    call log_init()
   
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
    
    ! Analysis begins from first configuration selected
    do i = 1, ncfs_from_to(1)
        ! In the first configuration basic initialization 
        if (first_configuration) call basic_init(use_cell,run_clusters,run_dyn,confined,nmol)
        ! Read i configuration from netcdf input file
        call cpu_time(t0)
        ! Jumps configurations to be read: ncstart controls starting conf in netcdf file
        ncstart = ncfs_from_to(2) + (i - 1)*(ncfs_from_to(3) - ncfs_from_to(2))/ncfs_from_to(1)
        ! Read-in a full configuration
        call read_nc_cfg(ncid_in, ncstart, io, io_log_file)
        ! Pre configuration analysis data transformations, corrections and format
        call reformat_input_conf(io,ncfs_from_to(1),i,ntypes,nsp)
        ! Exit the loop when EOF reached 
        if (io<0) exit
        call cpu_time(t1)
        ! Accumulate i/o time
        tread = tread + t1 - t0
        ! Over each configuration run selected modules (first initialize)
        if (first_configuration) call init_modules(use_cell, run_rdf, run_sq, run_clusters, nsp, nmol, nbcuda)
        ! Linked cells for cluster analysis
        if (use_cell) then
            call cells_build()
            call cells_reset_struct()
        end if
        !
        ! Transfer data to GPU
        call transfer_cpu_gpu(ndim)

        ! BFS cluster search
        if (run_clusters) call cluster_search()
        
        ! kinetic energy (if velocities available) 
        if (ex_vel) call thermo_kin(i, ndim)
        ! Run RDF
        if (run_rdf) call RDFcomp(Nmol, i, nbcuda, nthread)

        ! Compute density profile along idir direction
        if (confined) then
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
        if (run_dyn) call rtcorr(i)

        ! Print periodic output
        call print_output(i)
        first_configuration = .false.
    end do

    ! Normalize density profiles computed along the non-periodic dimension
    if (confined) then
        call normdenspr(nconf)
    end if

    ! Programme printouts
    call print_results(run_sq,  run_rdf, run_dyn, run_clusters, run_thermo, &
                        ntype, nsp, lsmax, nmol, nqmin, rcl)

    ! Cleaning house

    call clean_memory(run_sq,run_rdf,run_clusters,run_thermo,use_cell,run_dyn,confined)
    call cpu_time(time_cpu_stop)

    ! Print simulation time
    time_total = time_total + (time_cpu_stop - time_cpu_start)
    write (*, '(/,A,F15.7,A,/)') '**** Total time: ', time_total, ' seconds'
end program trj_analysis
