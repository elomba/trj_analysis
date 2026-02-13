module mod_util
   use mod_common, only : shmsize, maxthread,  run_thermo, ex_stress, &
      printDevPropShort, common_clear, nit, ener_name, press_name, &
      run_sq, run_sqw, run_rdf, run_clusters, run_dyn, &
      printcudaerror, ex_qc, nconf, qcharge, lsmax, maxcln
   use mod_densprof, only : prof_init, prof_clear
   use mod_sq, only : sq_init, printsq, sq_clear, sq_transfer_gpu_cpu
   use mod_rdf, only : rdf_init, printrdf, rdf_clear
   use mod_dyn, only : print_rtcor, dyn_clear, dyn_init
   use mod_log, only : log_clear, log_init, printPotEngCl, print_clusinfo, print_order
   use mod_clusters, only : clusters_clear, clusters_init, clusters_sq_init
   use mod_input, only : input_clear, sp_types_selected, ncfs_from_to, rdf_sq_cl_dyn_sqw_conf_ord, rcl, run_order
   use mod_cells, only : cells_init_post_nc_read, cells_init_pre_nc_read, cells_clear, use_cell
   use mod_thermo, only : thermo_clear
   use mod_nc_conf, only : wtypes, nmconf, orgty, wtypes
   use mod_order, only : order_init, order_clear, compute_order, norder
contains
   subroutine gpu_and_header(startEvent,stopEvent,gpudevice)
      use cudafor
      use mod_input, only : nthread
      implicit None
      type(cudaEvent), intent(inout) :: startEvent, stopEvent
      type(cudaDeviceProp) :: gpu_properties
      integer :: istat, gpudevice
      ! Get CUDA properties from device 0 (can be set from environmente variable CUDA_VISIBLE_DEVICES)
      
      istat = cudaSetDevice(gpudevice)
      if (istat == 0) then
         istat = cudaGetDeviceProperties(gpu_properties, gpudevice)
      else
         write(*,"('*** Unrecoverable error:  GPU ',i0,' is not available !!')") gpudevice
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
      write(*,"('*    Using GPU with CUDA nvfortran/nvcc >= 25.9/13.0',t80,'*')")
      write(*,"('*',t80,'*')")
      write(*,"('*    Version 1.0 January 2026',,t80,'*')")
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
      integer :: i
      if (nsp > ntypes) then
         print *, ' ERROR: number of species in input file is',nsp,' larger than that in netcdf file ',ntypes
         STOP
      else if (nsp == ntypes) Then
         nmol = natoms
      else if (nsp < ntypes) then
         write(*,'(" ** sp_types_selected: ",15i3)')sp_types_selected(1:nsp)
         write(*,'(" ** Types in trajectory=",15i3/)')orgty(1:ntypes)
         if (any(sp_types_selected == 0)) then
            write(*,'("*** Error: select the species to analyze")')
            write(*,'("    sp_types_selected MUST be defined from ")')
            write(*,'("    types in trajectory =",15i3//)')orgty(1:ntypes)
            stop
         endif
         do i = 1, nsp 
            if (any(orgty(1:ntypes) == sp_types_selected(i))) then
               continue
            else
               ! If the species type is not in the trajectory, stop
               write(   *,'(/"*** Error: species type ",i2," is not in trajectory")') sp_types_selected(i)
               write(*,'("    types in trajectory=",15i3)')orgty(1:ntypes)
               stop
            endif
         enddo
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
      if (use_cell.or.run_order) call cells_init_pre_nc_read(nmol)
      if (run_dyn) call dyn_init()
      if (confined) call prof_init()
   end subroutine basic_init

   subroutine form_dependencies()
      implicit none
      ! Radial ditribution function
      if (rdf_sq_cl_dyn_sqw_conf_ord(1) == .true.) run_rdf = .true.
      run_sq = .false.
      ! Static structure factors
      if (rdf_sq_cl_dyn_sqw_conf_ord(2) == .true.) run_sq = .true.
      ! Cluster analysis
      if (rdf_sq_cl_dyn_sqw_conf_ord(3) == .true.) then
         !
         ! Cluster analysis needs rdf's and s(q)'s to be computed
         ! This is also modified in input.f90
         !
         run_clusters = .true.
         run_rdf = .true.
      else
         ! Deactivate cells
         if (.not. rdf_sq_cl_dyn_sqw_conf_ord(7)) use_cell = .false.
      end if
      ! Dynamic correlations
      if (rdf_sq_cl_dyn_sqw_conf_ord(4) == .true.) run_dyn = .true.
      ! Activate analysis of dynamic structure factor
      if (rdf_sq_cl_dyn_sqw_conf_ord(5) == .true.) then
         run_sqw = .true.
         run_dyn = .true.
         run_sq = .true.
      endif
      ! Confined system (density profile analysis in the direction of confinement)
      if (rdf_sq_cl_dyn_sqw_conf_ord(6) == .true.) run_rdf = .true.
      ! Calculation of order parameters set to .true
      if (rdf_sq_cl_dyn_sqw_conf_ord(7) == .true.) run_order = .true.
   end subroutine form_dependencies

   subroutine init_modules(use_cell,run_rdf,run_sq,run_clusters,run_order,nsp,nmol,nbcuda)
      implicit none
      logical, intent(IN) :: use_cell,run_rdf,run_sq,run_clusters, run_order
      integer, intent(in) :: nsp, nmol, nbcuda
      if (use_cell) call cells_init_post_nc_read()
      if (run_rdf) call RDF_init(nsp)
      if (run_sq) call sq_init(nmol, nsp, nbcuda)
      if (run_clusters) call clusters_sq_init()
      if (run_order) call order_init(norder,nmol)
   end subroutine init_modules

   subroutine clean_memory(run_sq,run_rdf,run_clusters,run_thermo,use_cell,run_dyn,run_ord,confined)
      implicit none
      logical, intent(IN) :: run_sq,run_rdf,run_clusters,use_cell,run_thermo,run_dyn,run_ord,confined
      if (run_sq) then
         call sq_clear()
         print *, ' ··· sq_clear done'
      endif
      if (run_rdf) then
         call rdf_clear()
         print *, ' ··· rdf_clear done'
      endif
      if (run_clusters) Then 
         call clusters_clear()
         print *, ' ··· clusters_clear done'
      endif
      if (use_cell) then
         call cells_clear()
         print *, ' ··· cells_clear done'
      endif
      if (run_dyn) then 
         call dyn_clear()
         print *, ' ··· dyn_clear done'
      endif
      if (run_thermo) then
         call thermo_clear()
         print *, ' ··· thermo_clear done'
      endif
      if (run_ord) then
         call order_clear()
         print *, ' ··· order_clear done'
      endif
      if (confined) then 
         call prof_clear()
         print *, ' ··· prof_clear done'
      endif
      call common_clear()
      call log_clear()
      call input_clear()
   end subroutine clean_memory

   subroutine print_results(run_sq,  run_rdf, run_dyn, run_clusters, run_thermo, ntype, nsp, lsmax, nmol, nqmin, rcl)
      use mod_precision
      use mod_input, only : mat, charge, sp_labels, cl_thresh
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
            & " amu, average charge ",f9.5," e (",a4,")")') &
            & ntype(i), i, wtypes(i), mat(i), charge(i), sp_labels(i)
         else
            write (*, '(" ** ",i6," atoms of type ",i2,", in LAMMPS (",i2,"), ",f8.4, &
            & " amu,  (",a4,")")') ntype(i), i, wtypes(i), mat(i), sp_labels(i)
         endif
      end do
      if (ex_qc) write(*,'(/3x,"···   Net charge:",f12.5," e")')sum(qcharge(1:nmol))

      ! Print partial pair distribution functions
      if (run_rdf) then
         call printrdf(rcl, lsmax)
      end if
      if (run_dyn) then
         call print_rtcor()
      endif
      if (run_clusters) then
         ! Print cluster information only if number of big clusters > threshold
         if (maxcln>cl_thresh) call print_clusinfo(nqmin, Nmol)
         ! Print last cluster configuration
         call print_last_clustconf()
         if (run_thermo) then
            call printPotEngCl()
         end if
      end if
      if (run_order) then
         call print_order()
      endif
   end subroutine print_results

   subroutine print_last_clustconf()
      ! 
      ! Print last cluster configuration in LAMMPS format
      ! Atom type is set to original type + cluster size
      !
      use mod_common, only : cluster, itype, r, sidel, nstep, ex_vel, Nconf
      use mod_input,only : ndim
      use mod_nc_conf, only : org
      implicit none
      ! maxcolor is set to 32, so cluster size distinguished by color in VMD
      integer :: i, j, k, icl, id, io_lastclconf, natcl, maxcolor=32
      natcl = sum(cluster(1:maxcln)%clsize)
      open(newunit=io_lastclconf, file='last_clconf.lammpstrj', status='replace')
      write (io_lastclconf, "('ITEM: TIMESTEP'/I12/'ITEM: NUMBER OF ATOMS'/I12/'ITEM: BOX BOUNDS pp pp pp')") nstep, natcl
      write (io_lastclconf, "(2f15.7)") (org(i,1), org(i,1)+sidel(i), i=1, ndim)
      if (ndim == 2) write (io_lastclconf, "('-0.5 0.5')")
      if (ex_vel) then
         write (io_lastclconf, "('ITEM: ATOMS id type x y z vx vy vz')")
      else
         write (io_lastclconf, "('ITEM: ATOMS id type x y z')")
      end if
      icl = 0
      do i = 1, maxcln
         j = cluster(i)%clsize
         do k = 1, j
               icl = icl + 1
               id = cluster(i)%members(k)
               if (ndim == 3) then
                  write (io_lastclconf, "(I8,I4,3F15.7,3F15.7)") &
                     & icl, mod(j,maxcolor), r(1:ndim,id)+org(1:ndim,1)
               else
                  write (io_lastclconf, "(I8,I4,2F15.7,3F15.7)") &
                     & icl, mod(j,maxcolor), r(1:ndim,id)+org(1:ndim,1), 0.0
               end if
         end do
      end do
      close(io_lastclconf)
   end subroutine print_last_clustconf

   
   subroutine reformat_input_conf(io,final_conf,current_conf,ntypes,nsp)
      implicit none
      integer, intent(inout) :: final_conf
      integer, intent(in) :: io, ntypes, nsp, current_conf
      if (io<0) then
         final_conf=current_conf-1
      else
         ! Pre configuration analysis data transformations, corrections and format
            call select_ncdfinput()
      endif
   end subroutine reformat_input_conf
end module mod_util
