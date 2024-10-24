module mod_input
   use mod_precision
   use mod_common
   implicit none
   integer, allocatable, dimension(:) :: sp_types_selected, nw
   integer, dimension(3) :: ncfs_from_to
   character(len=2), allocatable, dimension(:) :: sp_labels
   integer :: nthread, ndim, jmin, minclsize, idir, nsp, nbuffer, potnbins=100, nqw=0, jump=10
   logical :: use_cell = .true.
   logical, dimension(6) :: rdf_sq_cl_dyn_sqw_conf
   real(myprec) :: deltar, rcl=-1.0, dcl, qmin, qmax, sigma, pwall, pwallp, rcrdf, tmax=100.0, potengmargin=0.0
   real(myprec), allocatable, dimension(:) :: mat, bsc, qw, tmqw
   character(len=128) :: input_filename, log_output_file, trj_input_file
   namelist /INPUT/ log_output_file, trj_input_file, ndim, nsp, nthread, &
      ncfs_from_to, rdf_sq_cl_dyn_sqw_conf, nqw, ener_name, press_name, potnbins, potengmargin 
   namelist /INPUT_SP/ sp_types_selected, sp_labels, mat 
   namelist /INPUT_RDF/ deltar, rcrdf, nrandom
   namelist /INPUT_SQ/ qmax, qmin, bsc
   namelist /INPUT_CL/ rcl, dcl, jmin, minclsize, sigma
   namelist /INPUT_CONF/ idir, pwall, pwallp
   namelist /INPUT_DYN/ nbuffer, tmax, jump
   namelist /INPUT_SQW/ qw, tmqw
contains

   subroutine read_input_file()
      integer :: io_input_file

      rdf_sq_cl_dyn_sqw_conf(:) = .false.
      open (newunit=io_input_file, file=input_filename, action='read')
      read (unit=io_input_file, nml=INPUT)
      allocate(sp_types_selected(nsp))
      sp_types_selected(:)=0
      allocate(sp_labels(nsp))
      allocate(mat(nsp))
      allocate(bsc(nsp))
      bsc(:) = 1.0
      read (unit=io_input_file, nml=INPUT_SP)

      deltar = -1.0
      if (rdf_sq_cl_dyn_sqw_conf(1) == .true. .or. rdf_sq_cl_dyn_sqw_conf(3) == .true. &
       & .or. rdf_sq_cl_dyn_sqw_conf(6) == .true.) read (unit=io_input_file, nml=INPUT_RDF)
      qmax = -1.0
      if (rdf_sq_cl_dyn_sqw_conf(2) == .true. &
      & .or. rdf_sq_cl_dyn_sqw_conf(3) == .true. .or. rdf_sq_cl_dyn_sqw_conf(5) == .true.  ) read (unit=io_input_file, nml=INPUT_SQ)
      rcl = -1.0
      if (rdf_sq_cl_dyn_sqw_conf(3) == .true.) read (unit=io_input_file, nml=INPUT_CL)
      if (rdf_sq_cl_dyn_sqw_conf(4) == .true. &
      &  .or. rdf_sq_cl_dyn_sqw_conf(5) == .true. ) read (unit=io_input_file, nml=INPUT_DYN)
      if (rdf_sq_cl_dyn_sqw_conf(5) == .true.) then 
         allocate(qw(nqw),nw(nqw),tmqw(nqw))
         tmqw(:) = 0.0
         read (unit=io_input_file, nml=INPUT_SQW)
      endif
      idir = 0
      if (rdf_sq_cl_dyn_sqw_conf(6) == .true.) then
         read (unit=io_input_file, nml=INPUT_CONF)
         confined = .true.
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
