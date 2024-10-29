program trj_analysis
    !
    !  This program performs an analysis of a netcdf trajectory file from LAMMPS
    !  molecular dynamics package (https://www.lammps.org/).
    !
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
    !
    !   The input is provided as a set of namelist (see attached example) 
    !
    !   Restrictions: 
    !                1. Only orthogonal cells are allowed
    !                2. Particle numbers must remain constant all along the trajectory
    !                3. Molecules are analyzed in terms of their atoms, 
    !                   so far no internal degrees of freedom taken into account
    !                4. LAMMPS units must be "real" (unit conversion soon to be implemented)
    !                   
    !
    ! 
    !   Important notice: in LAMMPS script the following computes must be included in the dump
    !                     in order to compute potential energies and pressures            
    !       compute stress all stress/atom NULL
    !       compute ener all pe/atom
    !       dump trj1 all netcdf ${Ndump} run.nc  id type x y z vx vy vz c_stress[*] c_ener
    !    In the first INPUT namelist optional character variables "ener_name" and "press_name" refer to the
    !    names of the computes
    !
    !    Usage: trj_analysis.exe input.nml (input file with sequence of namelists)
    !
    ! namelist /INPUT/ log_output_file, trj_input_file, ndim, nsp, nthread, &
    !       ncfs_from_to, rdf_sq_cl_dyn_sqw_conf, nqw, ener_name, press_name, 
    !       potnbins, potengmargin
    !       Name of log file, name of netcdf trajectory file, no. of dimensions (2,3), no. of species,
    !       no. of CUDA threads (default 128), no. of configurations-start-end, modules to run 
    !       (logical vars: RDF, structure factor, cluster analysis, dynamics, dynamic S(q,w),
    !       no. of Q's for dynamic analysis, name of compute for potential energy/atom in LAMMPS,
    !       name of stress/atom, no. of bins for energy histograms (def. 100), extra margins 
    !       in energy histograms (def. 0))
    ! namelist /INPUT_SP/ sp_types_selected, sp_labels, mat
    !          IDs of selected species (if nsp<ntypes in trajectory), character labels, atomic mass
    ! namelist /INPUT_RDF/ deltar, rcrdf, nrandom
    !          grid in RDF calculation, cut-off for rdf (default half box size), no. os random origins 
    !          for calculation of local number fluctuations
    ! namelist /INPUT_SQ/ qmax, qmin, bsc
    !          Max value of Q for S(Q), max value for full calculations (all Qs for 0<Q<=qmin), 
    !          scattering lengths 
    ! namelist /INPUT_CL/ rcl, dcl, jmin, minclsize, sigma
    !          Geometric clustering distance, grid for cluster distribution, minimum cluster size for analysis,
    !          Minimum cluster size to include in the trajectory of centers of mass, particle size
    ! namelist /INPUT_CONF/ idir, pwall, pwallp
    !          Direction of confinement (1,2,3->x,y,z), position of left wall, position of right wall
    ! namelist /INPUT_DYN/ nbuffer, tmax, tlimit, jump
    !           Number of buffers (time origins) for dynamic correlation analysis
    !           buffers are separated by jump configurations (def. 10), maximum time for 
    !           correlation functions (at tmax a window function is applied for FFTs), if omitted all
    !           t values are used.
    !           Averages stored over tlimit only (in ps), and then origin for calculation is shifted
    ! namelist /INPUT_SQW/ qw, tmqw
    !          values of Q to compute F(Q,t), Fs(Q,t) and S(Q,w),Ss(Q,w) maximum times for F(Q,t),
    !          if omitted tmax is used
    !
    !    OUTPUT FILES:
    !      * Thermodynamics
    !      - thermo_run.dat (Instantaneous values of thermod. quantities)
    !      * Dynamics
    !      - dyn.dat (msd, <v(t)v(0)>)
    !      - dynw.dat w, Z(w)
    !      - fkt.dat (F(Q_i,t) for nqw Qs)
    !      - fskt.dat (F_self(Q_i,t) for nqw Qs)
    !      - sqw.dat   (S(q,w), S_self(Q,w))
    !      - viscor.dat (<p_xy(t)p_xy(0)  eta(t) (shear viscosity integral))
    !
    !      * Structure
    !      - gmixsim.dat (g_ab(r), g_clcl(r) in cluster analysis on)
    !      - sq.dat  S_NN(Q), (S_cc(Q) , S_11, S_12, S_22 in binary systems)
    !      - sqcl.dat Cluster-cluster S(Q)
    !      - sqmix.dat S_ii (i<=nsp)
    !      - sqw.dat  (S(Q_i,w), S_self(Q_i,w) for nwq Qs)
    !
    !      * Cluster analysis
    !      - rhoprof.dat Average cluster density profile (only for finite clusters) 
    !      - radii.dat Distribution of cluster gyration radii
    !      - clustdistr.dat Distribution of cluster particle size
    !      - distUcl_N.dat Distribution of cluster internal energies per particle
    !      - distUcltot.dat Distribution of cluster internal energies (total)
    !      - clusevol.dat , conf no., no. of clusters, % of particles in clusters
    !      - centers.lammpstrj trajectory of clusters centers of mass (to be visualized with Ovito)
    !                          Particle no. not constant along the trajectory !!
    !
           
    !   OUTPUT Units: LAMMPS "real" units, except time (ps)
    !   Programmed in NVIDIA CUDA Fortran
    !
    !   A. Diaz-Pozuelo & E. Lomba, Madrid/Santiago de Compostela, fall 2024 
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
    use mod_log
    use mod_thermo
    use mod_dyn, only : dyn_init, dyn_clear, rtcorr, print_rtcor
    use mod_util, only : gpu_and_header, clean_memory, init_modules, reformat_input_conf, &
                         basic_init, print_results, select_species, reset_confs, form_dependencies
    use cudafor
    implicit none

    integer :: io = 0, ioerr, istat, ncid_in
    integer :: argc, ncstart, i
    logical :: first_configuration=.true. 
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

        ! Print periodic output
        call print_output(i)
        ! For next iteration then ...
        first_configuration = .false.
    end do

    ! Normalize density profiles computed along the non-periodic dimension
    if (confined) call normdenspr(nconf)
        
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
