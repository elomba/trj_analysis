!===============================================================================
! Program: trj_analysis
!===============================================================================
! Purpose:
!   GPU-accelerated analysis of LAMMPS NetCDF trajectory files for comprehensive
!   structural, thermodynamic, and dynamic characterization of molecular systems.
!   Source: https://www.lammps.org/
!
! Analysis Capabilities:
!   STRUCTURAL:
!     - Radial distribution functions (RDF) g_αβ(r) with species resolution
!     - Static structure factors S(Q) with adaptive Q-sampling
!     - Steinhardt orientational order parameters Q_l (2D and 3D)
!     - Density and charge density profiles for confined systems
!     - Cluster-cluster correlations (RDF and S(Q))
!
!   CLUSTER ANALYSIS:
!     - DBSCAN-based cluster identification
!     - Cluster size and radii distributions
!     - Average cluster density profiles
!     - Cluster shape descriptors (sphericity, cylindricity)
!     - Center-of-mass trajectory generation (centers.lammpstrj)
!     - Internal energy distributions (per particle and total)
!
!   DYNAMICS:
!     - Mean squared displacement (MSD) and velocity autocorrelation
!     - Intermediate scattering functions F(Q,t) and Fs(Q,t)
!     - Dynamic structure factors S(Q,ω) and Ss(Q,ω)
!     - Multi-buffer time correlation averaging with configurable origins
!     - Shear viscosity from stress tensor autocorrelation
!
!   THERMODYNAMICS:
!     - Kinetic energy and temperature (from velocities)
!     - Potential energy (from LAMMPS pe/atom compute)
!     - Pressure tensor (from LAMMPS stress/atom compute)
!     - Energy histograms with configurable binning
!
! Input Configuration:
!   Input file format: Fortran namelists (see README.md and examples/)
!   
!   Key namelists:
!     /INPUT/       - General parameters, module selection, file paths
!     /INPUT_SP/    - Species selection (types, labels, masses)
!     /INPUT_RDF/   - RDF grid and cutoff parameters
!     /INPUT_SQ/    - Structure factor Q-range and scattering lengths
!     /INPUT_CL/    - Cluster analysis thresholds and binning
!     /INPUT_ORD/   - Order parameter settings and neighbor criteria
!     /INPUT_CONF/  - Confinement geometry (wall positions)
!     /INPUT_DYN/   - Dynamics buffers, time limits, and windowing
!     /INPUT_SQW/   - Dynamic S(Q,ω) Q-values and time ranges
!
!   LAMMPS dump requirements for full functionality:
!     compute stress all stress/atom NULL
!     compute ener all pe/atom
!     dump trj1 all netcdf ${Ndump} run.nc id type x y z vx vy vz q c_stress[*] c_ener
!
!   Minimum requirements (positions only):
!     dump trj1 all netcdf ${Ndump} run.nc id type x y z
!
! Restrictions:
!   - NetCDF trajectory format required (LAMMPS dump netcdf)
!   - Orthogonal simulation cells only (no triclinic)
!   - Constant particle number (NpT allowed with minor S(Q) errors)
!   - Atomic-level analysis (molecular internal DOF not considered)
!   - LAMMPS units: "real" or "lj" (automatic conversion)
!
! Output Files:
!   THERMODYNAMICS: thermo_run.dat
!   STRUCTURE:      gmixsim.dat, sq.dat, sqmix.dat, sqcl.dat, order.dat
!   DYNAMICS:       dyn.dat, fkt.dat, fskt.dat, sqw.dat, viscor.dat, dynw.dat
!   CLUSTERS:       rhoprof.dat, radii.dat, clustdistr.dat, distUcl_N.dat,
!                   distUcltot.dat, clusevol.dat, fshape.dat, ordprof_clust.dat,
!                   ordprof_clcum.dat, order_per_cl.dat, centers.lammpstrj
!
! Usage:
!   ./trj_analysis.exe input.nml GPU_device_number (optional)
!
! Units:
!   Output: LAMMPS "real" units (time in ps) or "lj" units (reduced)
!
! Authors:
!   A. Díaz-Pozuelo & E. Lomba (optimized DBSCAN contributed by R. Lomba)
!   CSIC Madrid / USC Santiago de Compostela
!   February 2026
!
! Implementation:
!   NVIDIA CUDA Fortran with GPU acceleration
!===============================================================================
program trj_analysis
   
    use mod_precision
    use mod_common
    use mod_input
    use mod_nc
    use mod_nc_conf
    use mod_cells
    use mod_clusters
    use mod_sq
    use mod_rdf
    use mod_log
    use mod_thermo
    use mod_dyn, only : dyn_init, dyn_clear, rtcorr, print_rtcor
    use mod_util, only : gpu_and_header, clean_memory, init_modules, reformat_input_conf, &
                         basic_init, print_results, select_species, reset_confs, form_dependencies
    use cudafor
    implicit none

    integer :: io = 0, ioerr, istat, ncid_in
    integer :: argc, ncstart, i, devnum
    logical :: first_configuration=.true., nstepi0 = .false. 
    real :: t0 = 0, t1 = 0
    character(len=2) :: deviceNumber
    ! Command line arguments control
    argc = command_argument_count()
    if (argc < 1) then
        write(*,"('!!! Error: You must specify at least ONE argument: the name of the ')")
        write(*,"('    input file e.g.: trj_analisys.exe input_file.nml')")
        write(*,"('    Second argument (optional) is the GPU device number to use, (default 0)')")
        stop 
    end if
    ! Get input filename and optional GPU device number
    call get_command_argument(1, input_filename)
    if (argc ==2) then
        call get_command_argument(2, deviceNumber)
        read(deviceNumber,*) devNum
    else
        devNum = 0
    end if
    !
    ! initialize timers
    !
    call cpu_time(time_cpu_start)
    call cpu_time(cpu0)

    ! Get CUDA properties from device devNum 
    call gpu_and_header(startEvent,stopEvent,devNum)
   
  
    ! Load namelist input file & init log system
    call read_input_file()
    call log_init()
    !
    ! Open netcdf trajectory file, get the file identificator ncid_in
    call check(nf90_open(path=trj_input_file, mode=NF90_WRITE, ncid=ncid_in), ioerr)
    if (ioerr .ne. 0) then
        stop("** UNRECOVERABLE ERROR: cannot open NetCDF trajectory file !")
    endif 
    !
    ! First load from netcdf input file. Load global attributes of the simulation trajectory
    ! Read header and details of the NETCDF trajectory files: check consistency with input data
    !
    call read_nc_cfg(ncid_in, 1, io, io_log_file)

    ! Set number of species from netcdf file or from selected species from namelist
    call select_species(nsp, ntypes, nmol, natoms)
    if (nstep == 0) then
        write(*,"(' !!*** Warning: initial configuration for step 0, skipping ...  ')")
        write(*,"(' !!*** Change ncfs_from_to to n 2 m in input file to avoid this message'/)") 
        nstepi0 = .true.
    end if
    ! Reset configurations to read: if first arg=0 all confs in file are read
    call reset_confs(nconf_i,nconf)
    ! Modules to run: check for dependencies 
    !
    call form_dependencies()
    ! Init common variables & print outs
    call common_init(nmol, ndim, nthread, idir, conf(4)%units,conf(4)%scale, nsp)
    ! Analysis begins from first configuration selected
    do i = 1, ncfs_from_to(1)
        ! In the first configuration basic initialization 
        if (first_configuration) call basic_init(use_cell,run_clusters,run_dyn,confined,nmol)
        ! Read i-th configuration from netcdf input file
        call cpu_time(t0)
        ! Jumps configurations to be read: ncstart controls starting conf in netcdf file
        if (nstepi0 .and. ncfs_from_to(2)==1) then
            ncstart = ncfs_from_to(2) + (i - 1)*(ncfs_from_to(3) - ncfs_from_to(2))/ncfs_from_to(1)+1
        else
            ncstart = ncfs_from_to(2) + (i - 1)*(ncfs_from_to(3) - ncfs_from_to(2))/ncfs_from_to(1)
        endif
        ! Read-in a full configuration
        call read_nc_cfg(ncid_in, ncstart, io, io_log_file)
        ! Pre configuration analysis data transformations, corrections and format
        call reformat_input_conf(io,ncfs_from_to(1),i,ntypes,nsp)
        ! Exit the loop when EOF reached 
        if (io<0) exit
        ! Over each configuration run selected modules (first initialize)
        ! Accumulate i/o time
        call cpu_time(t1)
        if (first_configuration) call init_modules(use_cell, run_rdf, run_sq, run_clusters, run_order, nsp, nmol, nbcuda)
        ! For next iteration then ...
        first_configuration = .false.
        if (nstep == 0) then
            cycle
        endif
        tread = tread + t1 - t0
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
        if (confined) call profile_comp(nthread, ndim, idir, pwall, deltar)
        ! Compute SQ
        if (run_sq) call SQcalc()
         
        ! Compute cluster properties
        if (run_clusters) call cluster_analysis(i)
        ! Compute potential energy
        if (run_thermo) call poteng(i, natoms)

        ! Compute pressure
        if (ex_stress) call stress_calc(i, ndim, natoms)
     
        ! Compute dynamics
        if (run_dyn) call rtcorr(i)

        ! Compute order parameter
        if (run_order) call compute_order(nmol, ndim, rcl, sidel)
        ! Print periodic output
        call print_output(i)
    end do
    ! Normalize density profiles computed along the non-periodic dimension
    if (confined) call normdenspr(nconf)
    ! Programme printouts
    call print_results(run_sq,  run_rdf, run_dyn, run_clusters, run_thermo, &
                        ntype, nsp, lsmax, nmol, nqmin, rcl)
    ! Cleaning house
    write(*, '(/" **** Cleaning memory ...")')
    call clean_memory(run_sq,run_rdf,run_clusters,run_thermo,use_cell,run_dyn,run_order,confined)
    write(*, '(" **** Memory cleaned...")')
    call cpu_time(time_cpu_stop)

    ! Print simulation time
    time_total = time_total + (time_cpu_stop - time_cpu_start)
    write (*, '(/,A,F15.7,A,/)') '**** Total time: ', time_total, ' seconds'
end program trj_analysis
