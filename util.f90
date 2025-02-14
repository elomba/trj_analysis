module mod_util
   use mod_common, only : shmsize, maxthread,  run_thermo, ex_stress, &
      printDevPropShort, common_clear, nit, ener_name, press_name, &
      run_sq, run_sqw, run_rdf, run_clusters, run_dyn, &
      printcudaerror, ex_qc, nconf, qcharge
   use mod_densprof, only : prof_init, prof_clear
   use mod_sq, only : sq_init, printsq, sq_clear, sq_transfer_gpu_cpu
   use mod_rdf, only : rdf_init, printrdf, rdf_clear
   use mod_dyn, only : print_rtcor, dyn_clear, dyn_init
   use mod_log, only : log_clear, log_init, printPotEngCl, print_clusinfo
   use mod_clusters, only : clusters_clear, clusters_init, clusters_sq_init
   use mod_input, only : input_clear, sp_types_selected, ncfs_from_to, rdf_sq_cl_dyn_sqw_conf, rcl
   use mod_cells, only : cells_init_post_nc_read, cells_init_pre_nc_read, cells_clear, use_cell
   use mod_thermo, only : thermo_clear
   use mod_nc_conf, only : wtypes, nmconf, orgty, wtypes
contains
   subroutine gpu_and_header(startEvent,stopEvent)
      use cudafor
      use mod_input, only : nthread
      implicit None
      type(cudaEvent), intent(inout) :: startEvent, stopEvent
      type(cudaDeviceProp) :: gpu_properties
      integer :: istat
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
      ! Check that maximum number of threads is not surpassed
      if (nthread > maxthread/8) then
         nthread = maxthread/8
         write(*,'(" !! *** Warning: number of threads reset to",I3)')nthread
      endif
   end subroutine gpu_and_header

   subroutine select_species(nsp,ntypes,nmol,natoms)
      implicit none
      integer, intent(in) :: nsp, ntypes, natoms
      integer, intent(inout) :: nmol
      if (nsp > ntypes) then
         print *, ' ERROR: number of species in input file is',nsp,' larger than that in netcdf file ',ntypes
         STOP
      else if (nsp == ntypes) Then
         nmol = natoms
      else if (nsp < ntypes) then
         if (any(sp_types_selected == 0)) then
            write(*,'("*** Error: select the species to analyze")')
            write(*,'("    when less than those in trajectory: must be among ")')
            write(*,'("    sp_types_selected=",15i3)')orgty(1:ntypes)
            stop
         endif
         wtypes(1:nsp) = sp_types_selected(1:nsp)
         call reset_nmol(nmol)
         run_thermo = .false.
         ex_stress = .false.
         ener_name = "XXX"
         press_name = "XXX"
      end if
      ! Number of different interactions
      nit = nsp*(nsp + 1)/2
   end subroutine select_species

   subroutine reset_confs(nconf_i,nconf)
      implicit None
      integer, intent(in) :: nconf_i
      integer, intent(inout) :: nconf
      if (ncfs_from_to(1) == 0) then
         ncfs_from_to(1) = nconf_i
         ncfs_from_to(2) = 1
         ncfs_from_to(3) = nconf_i
      end if
      nconf = ncfs_from_to(1)
      ncfs_from_to(3) = ncfs_from_to(3) + 1
      if (nconf > nmconf) then
         write(*,"(/' !!*** Warning: resetting Nconf to ',i5,', last conf. in trajectory !'/)") nmconf
         nconf = nmconf
      endif
   end subroutine reset_confs

   subroutine basic_init(use_cell,run_clusters,run_dyn,confined,nmol)
      implicit none
      logical, intent(in) :: use_cell,run_clusters,run_dyn, confined
      integer, intent(in) :: nmol
      if (use_cell) call cells_init_pre_nc_read(nmol)
      if (run_clusters) call clusters_init(nmol)
      if (run_dyn) call dyn_init()
      if (confined) call prof_init()
   end subroutine basic_init

   subroutine form_dependencies()
      implicit none
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
      else
         ! Deactivate cells
         use_cell = .false.
      end if
      ! Dynamic correlations
      if (rdf_sq_cl_dyn_sqw_conf(4) == .true.) run_dyn = .true.
      ! Activate analysis of dynamic structure factor
      if (rdf_sq_cl_dyn_sqw_conf(5) == .true.) then
         run_sqw = .true.
         run_dyn = .true.
         run_sq = .true.
      endif
      ! Confined system (density profile analysis in the direction of confinement)
      if (rdf_sq_cl_dyn_sqw_conf(6) == .true.) run_rdf = .true.
   end subroutine form_dependencies

   subroutine init_modules(use_cell,run_rdf,run_sq,run_clusters,nsp,nmol,nbcuda)
      implicit none
      logical, intent(IN) :: use_cell,run_rdf,run_sq,run_clusters
      integer, intent(in) :: nsp, nmol, nbcuda
      if (use_cell) call cells_init_post_nc_read()
      if (run_rdf) call RDF_init(nsp)
      if (run_sq) call sq_init(nmol, nsp, nbcuda)
      if (run_clusters) call clusters_sq_init()
   end subroutine init_modules

   subroutine clean_memory(run_sq,run_rdf,run_clusters,run_thermo,use_cell,run_dyn,confined)
      implicit none
      logical, intent(IN) :: run_sq,run_rdf,run_clusters,use_cell,run_thermo,run_dyn,confined
      if (run_sq) call sq_clear()
      if (run_rdf) call rdf_clear()
      if (run_clusters) call clusters_clear()
      if (use_cell) call cells_clear()
      if (run_dyn) call dyn_clear()
      if (confined) call prof_clear()
      if (run_thermo) call thermo_clear()
      call common_clear()
      call log_clear()
      call input_clear()
   end subroutine clean_memory

   subroutine print_results(run_sq,  run_rdf, run_dyn, run_clusters, run_thermo, ntype, nsp, lsmax, nmol, nqmin, rcl)
      use mod_precision
      use mod_input, only : mat, charge, sp_labels
      implicit none
      real(myprec) :: rcl
      integer, intent(in) :: nsp, lsmax, nmol, nqmin
      integer, dimension(nsp), intent(in) :: ntype
      logical, intent(in) :: run_sq,  run_rdf, run_dyn, run_clusters, run_thermo
      integer :: i

      ! Print out S(Q)'s
      if (run_sq) then
         call sq_transfer_gpu_cpu()
         call printSQ(Nmol)
      end if

      write(*,"(/60('-'))")
      do i = 1, nsp
         if (ex_qc) then
            write (*, '(" ** ",i6," atoms of type ",i2,", in LAMMPS (",i2,"), ",f8.4,&
            & " amu, average charge ",f8.4," e (",a4,")")') &
            & ntype(i), i, wtypes(i), mat(i), charge(i), sp_labels(i)
         else
            write (*, '(" ** ",i6," atoms of type ",i2,", in LAMMPS (",i2,"), ",f8.4, &
            & " amu,  (",a4,")")') ntype(i), i, wtypes(i), mat(i), sp_labels(i)
         endif
      end do
      if (ex_qc) write(*,'("   Net charge:",f10.5," e")')sum(qcharge(1:nmol))

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
         end if
      end if

   end subroutine print_results

   subroutine reformat_input_conf(io,final_conf,current_conf,ntypes,nsp)
      implicit none
      integer, intent(inout) :: final_conf
      integer, intent(in) :: io, ntypes, nsp, current_conf
      if (io<0) then
         final_conf=current_conf-1
      else
         ! Pre configuration analysis data transformations, corrections and format
         if (ntypes == nsp) then
            ! All species selected
            call trans_ncdfinput()
         else
            ! Only some species selected from the configuration
            call select_ncdfinput()
         endif
      endif
   end subroutine reformat_input_conf
end module mod_util
