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
   integer, dimension(3) :: ncfs_from_to=3*0
   !
   ! nthread : number of threads for CUDA kernels is preset to 64 by default, beware of larger values for big systems
   !
   integer :: nthread=64, ndim, minPts,  idir=3, nsp, nbuffer=2, potnbins=100, nqw=0, &
            & jump=1, norder=1, nnbond=0, cl_thresh=10, nprint=10
   logical :: use_cell = .true., run_order = .false., print_orderp=.false., &
               geometry=.true.
   logical, dimension(7) :: rdf_sq_cl_dyn_sqw_conf_ord
   real(myprec) :: deltar, dcl, qmin, qmax, rcrdf, rclcl=0.0, &
      tmax=-1, tmaxp=-1, tlimit=-1, potengmargin=0.0
   ! mat and bsc at atomic mass and scattering length, respectively by species, masa and mscat correspond to individual atoms
   real(myprec), allocatable, dimension(:) :: mat, bsc, charge, qw, tmqw
   integer, allocatable, dimension(:) :: orderp
   character(len=128) :: input_filename, log_output_file, trj_input_file, system_data_file="system.data"
   !
   ! Input namelists
   namelist /INPUT/ log_output_file, trj_input_file, ndim, nsp, nthread,  &
      & ncfs_from_to,  rdf_sq_cl_dyn_sqw_conf_ord, nqw, nslice, norder, ener_name, &
      & press_name, potnbins, potengmargin, rcl, periodic, nprint, topol, &
      &  system_data_file
   namelist /INPUT_SP/ sp_types_selected, mat, rigid, nmrigid, rigid_mols 
   namelist /INPUT_RDF/ deltar, rcrdf, nrandom
   namelist /INPUT_SQ/ qmax, qmin, bsc
   namelist /INPUT_CL/ dcl, minPts, ndrclus, cl_thresh, geometry
   namelist /INPUT_CONF/ zslice, zgrid
   namelist /INPUT_DYN/ nbuffer, tmax, tmaxp, tlimit, jump
   namelist /INPUT_SQW/ qw, tmqw
   namelist /INPUT_ORDER/ orderp, print_orderp, nnbond, rclcl
contains

   subroutine read_input_file()
      integer :: io_input_file, i
      logical :: trj_file_exists=.true.

      ! Initialize all module flags to false
      rdf_sq_cl_dyn_sqw_conf_ord(:) = .false.
      ! Open and read main INPUT namelist
      rigid_mols(:) = -1
      open (newunit=io_input_file, file=input_filename, action='read')
      read (unit=io_input_file, nml=INPUT)
      if (rdf_sq_cl_dyn_sqw_conf_ord(3) == .true.) then
         minPts = 2*ndim+1
      endif
      !
      ! Rigid rigid degrees of freedom require a topology file
      !
      if (rigid.and..not.topol) Then
         write(*,'(//"*** Fatal error: topology analysis must be activated (topol=.true.) if rigid=.true. !")')
         stop
      endif
      !
      ! Allocate arrays for species properties
      allocate(sp_types_selected(nsp))
      !
      sp_types_selected(:) = 0
      allocate(mat(nsp))
      allocate(charge(nsp))
      charge(:) = 0.0
      allocate(bsc(nsp))  ! Neutron scattering lengths
      bsc(:) = 1.0
      !
      ! Read species-specific parameters
      read (unit=io_input_file, nml=INPUT_SP)
      !
      ! Check existence of trajectory file
      inquire (file=trj_input_file, exist=trj_file_exists)
      if (.not. trj_file_exists) then
         write(*,'("*** Error: trajectory file specified in trj_input_file ",A, " does not exist !")') char(27)//'[33m'//trim(trj_input_file)//char(27)//'[0m'
         stop 
      endif
      if (sp_types_selected(1) == 0) then
         write(*,'("*** Error: atom IDs must be specified in sp_types_selected !")')
         stop
      endif
      !
      ! Check for rigid molecules
      !
      if (rigid) then
         if (nmrigid==0 .or. nmrigid > nrig_max) then
            write(*, '(" *** Fatal error: number of rigid mol types must be defined and less than ",I2," !")') nrig_max
            stop
         endif
         if (any(rigid_mols(1:nmrigid) < 0)) then
            write(*, '(" *** Fatal error: define all rigid mol types !")')
            stop
         endif
      endif
      !
      ! Define defaults for configuration file parameters if not set by user
      if (ncfs_from_to(2) == 0) then
         ncfs_from_to(2) = 1
         write(*,'(" ** Warning: ncfs_from_to(2) (starting configuration index) reset to 1 !")')
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
      if (rdf_sq_cl_dyn_sqw_conf_ord(6) == .true.) then
         if (ndim /= 3) then
            write(*,'("*** Error: system must be 3D to compute confinement properties !")')
            stop
         endif
         if (nslice < 1) then
            nslice = 1
            write(*,'("*** Warning: nslice (number of slices for confinement profile) reset to 1 !")')
            ! Definition of z-slices postponed until we read first configuration to 
            ! get box size in z
            auto_zslice = .true.
         endif
         allocate(zslice(nslice),countslice(nslice))
         allocate(countpslice(nsp, nslice))
         read (unit=io_input_file, nml=INPUT_CONF)
         countslice(:) = 0
         countsliced(:) = 0
         confined = .true.
         if (rdf_sq_cl_dyn_sqw_conf_ord(1) == .true. .or. rdf_sq_cl_dyn_sqw_conf_ord(2) == .true.) then
            twoDstruc_3D = .true.
         endif
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
      deallocate(mat)
      deallocate(bsc)

   end subroutine input_clear

end module mod_input
