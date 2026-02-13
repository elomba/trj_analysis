!===============================================================================
! Module: mod_input
!===============================================================================
! Purpose:
!   Reads and processes input parameters from namelist files for trajectory
!   analysis. Handles validation and initialization of user-specified
!   simulation parameters across multiple analysis modules.
!
! Key Functionality:
!   - Reads multiple namelists from input file:
!     * INPUT:       General parameters (dimensions, species, modules to run)
!     * INPUT_SP:    Species selection and properties
!     * INPUT_RDF:   RDF calculation parameters
!     * INPUT_SQ:    Structure factor parameters
!     * INPUT_CL:    Cluster analysis parameters
!     * INPUT_CONF:  Confinement geometry parameters
!     * INPUT_DYN:   Dynamics correlation parameters
!     * INPUT_SQW:   Dynamic structure factor parameters
!     * INPUT_ORDER: Order parameter calculation settings
!   - Validates parameter consistency and ranges
!   - Sets default values for optional parameters
!
! Main Variables:
!   nthread              - CUDA threads per block (default 128)
!   ndim                 - System dimensions (2 or 3)
!   nsp                  - Number of species to analyze
!   rcl                  - Cutoff distance for clustering and nearest neighbors
!   sp_types_selected    - Array of selected species IDs
!   rdf_sq_cl_dyn_sqw_conf_ord - Logical flags for modules to execute
!
! Subroutines:
!   read_input_file()  - Main routine to read and validate all input namelists
!   input_clear()      - Deallocates input arrays
!
! Notes:
!   - Input file must contain INPUT namelist at minimum
!   - Species selection mandatory via sp_types_selected
!   - Module dependencies automatically checked
!===============================================================================
module mod_input
   use mod_precision
   use mod_common
   use cudafor
   use sorts
   implicit none
   integer, allocatable, dimension(:) :: sp_types_selected, nw
   integer, dimension(3) :: ncfs_from_to
   character(len=4), allocatable, dimension(:) :: sp_labels
   !
   ! nthread : number of threads for CUDA kernels is preset to 64 by default, beware of larger values for big systems
   !
   integer :: nthread=128, ndim, kmin=5,  idir=0, nsp, nbuffer=2, potnbins=100, nqw=0, &
            & jump=1, norder=1, nnbond=0, cl_thresh=10
   logical :: use_cell = .true., run_order = .false., print_orderp=.false.
   logical, dimension(7) :: rdf_sq_cl_dyn_sqw_conf_ord
   real(myprec) :: deltar, rcl=-1.0, dcl, qmin, qmax, rcrdf, rclcl=0.0, &
      tmax=-1, tmaxp=-1, tlimit=-1, potengmargin=0.0
   real(myprec), allocatable, dimension(:) :: mat, bsc, charge, qw, tmqw
   integer, allocatable, dimension(:) :: orderp
   character(len=128) :: input_filename, log_output_file, trj_input_file
   !
   ! Input namelists
   namelist /INPUT/ log_output_file, trj_input_file, ndim, nsp, nthread,  &
      & ncfs_from_to,  rdf_sq_cl_dyn_sqw_conf_ord, nqw, norder, ener_name, &
      & press_name, potnbins, potengmargin, rcl, periodic
   namelist /INPUT_SP/ sp_types_selected, sp_labels, mat
   namelist /INPUT_RDF/ deltar, rcrdf, nrandom
   namelist /INPUT_SQ/ qmax, qmin, bsc
   namelist /INPUT_CL/ dcl, kmin, ndrclus, cl_thresh
   namelist /INPUT_CONF/ idir
   namelist /INPUT_DYN/ nbuffer, tmax, tmaxp, tlimit, jump
   namelist /INPUT_SQW/ qw, tmqw
   namelist /INPUT_ORDER/ orderp, print_orderp, nnbond, rclcl
contains

   subroutine read_input_file()
      integer :: io_input_file, i

      ! Initialize all module flags to false
      rdf_sq_cl_dyn_sqw_conf_ord(:) = .false.
      ! Open and read main INPUT namelist
      open (newunit=io_input_file, file=input_filename, action='read')
      read (unit=io_input_file, nml=INPUT)
      if (run_clusters) then
         write(*,"(' *** Note: rcl (NN and/or connectivity distance) set to ',f8.4,' Å')" ) rcl
      endif
      ! Allocate arrays for species properties
      allocate(sp_types_selected(nsp))
      !
      sp_types_selected(:) = 0
      allocate(sp_labels(nsp))
      allocate(mat(nsp))
      allocate(charge(nsp))
      charge(:) = 0.0
      allocate(bsc(nsp))  ! Neutron scattering lengths
      bsc(:) = 1.0
      ! Read species-specific parameters
      read (unit=io_input_file, nml=INPUT_SP)
      if (sp_types_selected(1) == 0) then
         write(*,'("*** Error: atom IDs must be specified in sp_types_selected !")')
         stop
      endif
      ! Read RDF parameters if needed (RDF, cluster analysis, or confinement)
      deltar = -1.0
      if (rdf_sq_cl_dyn_sqw_conf_ord(1) == .true. .or. rdf_sq_cl_dyn_sqw_conf_ord(3) == .true. &
      & .or. rdf_sq_cl_dyn_sqw_conf_ord(6) == .true.) read (unit=io_input_file, nml=INPUT_RDF)
      ! Read structure factor parameters if needed (S(q), cluster analysis, or S(q,w))
      qmax = -1.0
      if (rdf_sq_cl_dyn_sqw_conf_ord(2) == .true. &
      & .or. rdf_sq_cl_dyn_sqw_conf_ord(3) == .true. .or. rdf_sq_cl_dyn_sqw_conf_ord(5) == .true.  ) read (unit=io_input_file, nml=INPUT_SQ)
      if (rdf_sq_cl_dyn_sqw_conf_ord(3) == .true.) then
         read (unit=io_input_file, nml=INPUT_CL)
         if (rcl <= 0.0) then
            write(*,'("*** Error: rcl (cluster distance) must be positive to compute clusters !")')
            stop
         endif
      endif
      if (rdf_sq_cl_dyn_sqw_conf_ord(4) == .true. &
      &  .or. rdf_sq_cl_dyn_sqw_conf_ord(5) == .true. ) then 
         read (unit=io_input_file, nml=INPUT_DYN)
         if (nbuffer<2) then
            write(*,'("*** Warning: nbuffer (number of buffers for time correlations) reset to 2 !")')
         endif
      endif 
      if (rdf_sq_cl_dyn_sqw_conf_ord(5) == .true.) then
         allocate(qw(nqw),nw(nqw),tmqw(nqw))
         tmqw(:) = 0.0
         read (unit=io_input_file, nml=INPUT_SQW)
      endif
      idir = 0
      if (rdf_sq_cl_dyn_sqw_conf_ord(6) == .true.) then
         read (unit=io_input_file, nml=INPUT_CONF)
         confined = .true.
      endif
      if (rdf_sq_cl_dyn_sqw_conf_ord(7) == .true.) then
         if (norder < 1) then
            write(*,'("*** Error: norder must be at least 1 !")')
            stop
         endif
         allocate(orderp(norder))
         orderp(:) = 0
         read (unit=io_input_file, nml=INPUT_ORDER)
          if (rclcl <= 0.0 .and. rdf_sq_cl_dyn_sqw_conf_ord(3) == .true.) then
            write(*,'("*** Error: rclcl must be positive !")')
            stop
         endif
         if (rcl <= 0.0) then
            write(*,'("*** Error: rcl (NN distance) must be defined to compute order parameters !")')
            stop
         endif
         call quicksort_nri(orderp)
         if (orderp(norder) > llmax.and.ndim==3) then
            write(*,'("*** Error: order parameter index in orderp exceeds llmax = ",I3,": redimension in order.cuf !")') llmax
            stop
         endif
      endif
      close (io_input_file)

   end subroutine read_input_file

   subroutine input_clear()

      deallocate(sp_types_selected)
      deallocate(sp_labels)
      deallocate(mat)
      deallocate(bsc)

   end subroutine input_clear

end module mod_input
