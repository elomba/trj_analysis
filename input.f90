module mod_input
   use mod_precision
   use mod_common
   implicit none
   integer, allocatable, dimension(:) :: sp_types_selected
   integer, dimension(3) :: ncfs_from_to
   character(len=2), allocatable, dimension(:) :: sp_labels
   integer :: nthread, ndim, jmin, minclsize, idir, nsp, nbuffer, thermo_junk
   logical :: use_cell = .true.
   logical, dimension(6) :: rdf_sq_cl_dyn_thermo_conf
   real(myprec) :: deltar, rcl, dcl, qmin, qmax, sigma, pwall, pwallp, rcrdf, tmax
   real(myprec), allocatable, dimension(:) :: mat, bsc
   character(len=128) :: input_filename, log_output_file, trj_input_file
   namelist /INPUT/ log_output_file, trj_input_file, ndim, nsp, nthread, &
      ncfs_from_to, rdf_sq_cl_dyn_thermo_conf
   namelist /INPUT_SP/ sp_types_selected, sp_labels, mat, bsc
   namelist /INPUT_RDF/ deltar, rcrdf
   namelist /INPUT_SQ/ qmax, qmin
   namelist /INPUT_CL/ rcl, dcl, jmin, minclsize, sigma
   namelist /INPUT_CONF/ idir, pwall, pwallp
   namelist /INPUT_DYN/ nbuffer, tmax
   namelist /INPUT_THERMO/ thermo_junk
contains

   subroutine read_input_file()
      integer :: io_input_file

      rdf_sq_cl_dyn_thermo_conf(:) = .false.
      open (newunit=io_input_file, file=input_filename, action='read')
      read (unit=io_input_file, nml=INPUT)
      
      allocate(sp_types_selected(nsp))
      sp_types_selected(:)=0
      allocate(sp_labels(nsp))
      allocate(mat(nsp))
      allocate(bsc(nsp))
      read (unit=io_input_file, nml=INPUT_SP)

      deltar = -1.0
      if (rdf_sq_cl_dyn_thermo_conf(1) == .true. .or. rdf_sq_cl_dyn_thermo_conf(3) == .true.) read (unit=io_input_file, nml=INPUT_RDF)
      qmax = -1.0
      if (rdf_sq_cl_dyn_thermo_conf(2) == .true. .or. rdf_sq_cl_dyn_thermo_conf(3) == .true.) read (unit=io_input_file, nml=INPUT_SQ)
      rcl = -1.0
      if (rdf_sq_cl_dyn_thermo_conf(3) == .true.) read (unit=io_input_file, nml=INPUT_CL)
      if (rdf_sq_cl_dyn_thermo_conf(4) == .true.) read (unit=io_input_file, nml=INPUT_DYN)
      if (rdf_sq_cl_dyn_thermo_conf(5) == .true.) read (unit=io_input_file, nml=INPUT_THERMO)
      idir = 0
      if (rdf_sq_cl_dyn_thermo_conf(6) == .true.) read (unit=io_input_file, nml=INPUT_CONF)

      close (io_input_file)

   end subroutine read_input_file

   subroutine input_clear()
      
      deallocate(sp_types_selected)
      deallocate(sp_labels)
      deallocate(mat)
      deallocate(bsc)

   end subroutine input_clear

end module mod_input
