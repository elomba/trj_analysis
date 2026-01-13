module mod_input
   use mod_precision
   use mod_common
   use cudafor
   use sorts
   implicit none
   integer, allocatable, dimension(:) :: sp_types_selected, nw
   integer, dimension(3) :: ncfs_from_to
   character(len=4), allocatable, dimension(:) :: sp_labels
   integer :: nthread, ndim, jmin=3, minclsize, idir=0, nsp, nbuffer=2, potnbins=100, nqw=0, &
            & jump=1, norder=1, nnbond=6
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
   namelist /INPUT_CL/ dcl, jmin, minclsize, ndrclus
   namelist /INPUT_CONF/ idir
   namelist /INPUT_DYN/ nbuffer, tmax, tmaxp, tlimit, jump
   namelist /INPUT_SQW/ qw, tmqw
   namelist /INPUT_ORDER/ orderp, print_orderp, nnbond, rclcl
contains

   subroutine read_input_file()
      integer :: io_input_file, i

      rdf_sq_cl_dyn_sqw_conf_ord(:) = .false.
      open (newunit=io_input_file, file=input_filename, action='read')
      read (unit=io_input_file, nml=INPUT)
      if (rcl>0.0) then
         write(*,"(' *** Note: rcl (NN and/or connectivity distance) set to ',f8.4,' Å')") rcl
      endif
      allocate(sp_types_selected(nsp))
      !
      sp_types_selected(:) = 0
      allocate(sp_labels(nsp))
      allocate(mat(nsp))
      allocate(charge(nsp))
      charge(:) = 0.0
      allocate(bsc(nsp))
      bsc(:) = 1.0
      read (unit=io_input_file, nml=INPUT_SP)
      if (sp_types_selected(1) == 0) then
         write(*,'("*** Error: atom IDs must be specified in sp_types_selected !")')
         stop
      endif
      deltar = -1.0
      if (rdf_sq_cl_dyn_sqw_conf_ord(1) == .true. .or. rdf_sq_cl_dyn_sqw_conf_ord(3) == .true. &
      & .or. rdf_sq_cl_dyn_sqw_conf_ord(6) == .true.) read (unit=io_input_file, nml=INPUT_RDF)
      qmax = -1.0
      if (rdf_sq_cl_dyn_sqw_conf_ord(2) == .true. &
      & .or. rdf_sq_cl_dyn_sqw_conf_ord(3) == .true. .or. rdf_sq_cl_dyn_sqw_conf_ord(5) == .true.  ) read (unit=io_input_file, nml=INPUT_SQ)
      if (rdf_sq_cl_dyn_sqw_conf_ord(3) == .true.) then
         minclsize = jmin
         read (unit=io_input_file, nml=INPUT_CL)
         if (rcl <= 0.0) then
            write(*,'("*** Error: rcl (cluster distance) must be positive to compute clusters !")')
            stop
         endif
         if (jmin > 5) then
            write(*,'("*** Warning: jmin (minimum cluster size) reset to 5 !")')
            jmin = 5
         endif
         if (minclsize < jmin) then
            write(*,'("*** Warning : minclsize reset to ",I3," !")') jmin
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
