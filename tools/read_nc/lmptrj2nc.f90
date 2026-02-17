!===============================================================================
! Program: lmptrj2nc (read_nc)
!===============================================================================
! Purpose:
!   Reads LAMMPS/AMBER NetCDF trajectory files, extracts configuration data,
!   and rewrites to new NetCDF4 format. Allows selective atom type filtering
!   and demonstrates NetCDF frame-by-frame I/O operations.
!
! Functionality:
!   - Reads LAMMPS NetCDF trajectory dump files
!   - Extracts positions, velocities, types, and metadata
!   - Supports selective atom type extraction
!   - Writes filtered data to new NetCDF4 file
!   - Handles non-periodic dimensions with automatic box size adjustment
!
! Input (Interactive):
!   - Number of configurations to read
!   - Configuration range [start, end]
!   - Input NetCDF filename
!   - Atom ID for diagnostic output
!   - Atom types to select (or 'all')
!
! Box Size Handling:
!   LAMMPS NetCDF format sets box length = 0 for non-periodic dimensions.
!   This program redefines lengths using:
!     - abs(max(r_xyz) - min(r_xyz)) + 10
!     - or 2*abs(origin_xyz) depending on dimension
!
! Compilation:
!   ifort -I /path/to/netcdf/include -L/path/to/netcdf/lib lmptrj2nc.f90 -lnetcdff
!   export LD_LIBRARY_PATH=/path/to/netcdf/lib:$LD_LIBRARY_PATH
!
! Author:
!   Enrique Lomba, IQF, Madrid, November 2023
!===============================================================================
program read_nc
  use myncdf
  use configuration
  use netcdf
  !  use ncio
  implicit none
  integer :: ncid_in, ncid_out, ncfs, i, j, ncst, atom 
  integer  :: ncstart, ncmin, ncmax, ntms, ioerr, io=0
  logical :: selectall=.false.
  write(*,"(' Enter number of confs to be read : ')",advance="NO")
  read(*,*) ncfs
  print *, " Input starting and final configuration to read in :"
  read(*,*) ncmin, ncmax
  ncmax = ncmax+1
  print *, " File name "
  read(*,*) fname
  Print *, " ID of atom to print out "
  read(*,*) atom
  print *, ' No. of types to select (-1 for all)'
  read(*,*) nwty
  if (nwty > 0) then
     allocate(wtypes(nwty))
     print *, " Types to select "
     read(*,*) wtypes(:)
  else
     selectall=.true.
  endif
  if (index(fname,".nc") > 0) print *, "** Will attempt reading NETCDF file format"
  call check(nf90_open(path = fname, mode=NF90_WRITE, ncid = ncid_in),ioerr)
  !
  if (ioerr /= nf90_noerr) then
     stop(' ** Error opening file '//fname)
  end if
     
  fnamew = "nuevo.nc"
  ! Possible cmodes : NF90_CLOBBER, NF90_CDF5, NF90_CLASSIC_MODEL, NF90_64BIT_OFFSET, NF90_NETCDF4
!  call check(nf90_create(path = fnamew, cmode=NF90_CLOBBER, ncid =&
  call check(nf90_create(path = fnamew, cmode=NF90_NETCDF4, ncid =&
       & ncid_out),ioerr)

  write(*, "(/,' ####### Print some output from the file ',A)")fname

  do j = 1, ncfs
     ncstart = ncmin+(j-1)*(ncmax-ncmin)/ncfs
     call read_lmptrj(io)
     if (io < 0) then
        write(*,"(' ** Error: reading past end of file: exiting loop '/)")
        exit
     end if
     Write(*,"(/' ** Reading  =',i9,' atoms from configuration'&
          &,i7,' of ', i7)")natoms, j, nconf
     ! no. of atoms to read from file
     count(2) = natoms
     write(*,"(10x,'Spatial =',A)") spatial(1:3)
     write(*,"(10x,'Cell spatial =',A)") cell_spatial
     write(*,"(10x,'Cell spatial =',3(A,' '))") cell_ang
     Write(*,"(10x,'Step=',i12,' Time Step ',f6.4)")step, scale
     Write(*,"(10x,'Origin at ',3f15.7)") org
     Write(*,"(10x,'Cell lengths ',3f15.7)") cell
     Write(*,"(10x,'Cell angles ',3f15.7)") cell_a
     do i=atom, atom
        write(*,"(10x,' Frame=',i10)") ncstart
        write(*,"(10x,' atom id=',i10)") idi(i,1)
        write(*,"(10x,' atom type=',i2)") ity(i,1)
        if(ex_mol) write(*,"(10x,' molecule id=',i10)") imol(i,1)
        write(*,"(10x,' coordinates=',3f15.7)")r(1:3,i,1)
        if (ex_vel) write(*,"(10x,' velocities=',3f15.7)") v(1:3,i,1)
     enddo
     if (.not.selectall) then
        ntms = 0
        do i=1, nwty
           ntms=ntms+atypes(wtypes(i))
        end do
        call write_nc_cfg_red(ncid_out,j,ntms)
     else
        call write_nc_cfg(ncid_out,j)
        ntms = natoms
     endif
     Write(*,"(//'###### Dumping configuration no. ',i8,' with ',i7,' atoms ####')")ncstart, ntms
  enddo
  call check(nf90_close(ncid_out),ioerr)
  call check(nf90_close(ncid_in),ioerr)
end program read_nc
 
