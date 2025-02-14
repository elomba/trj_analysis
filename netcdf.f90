module mod_nc
   use mod_precision
   use mod_common, Only : ex_vel, ex_force, ex_stress, run_thermo, &
      u_p, stress, ex_mol, ex_qc, periodic, voigt, &
      ener_name, press_name
   interface
      subroutine read_nc_cfg(ncid, ncstart, io, unit)
         integer, intent(in) :: ncid, ncstart
         integer, intent(in), optional :: unit
         integer, intent(out) :: io
      end subroutine read_nc_cfg
      SUBROUTINE check(istatus, ioerr)
         INTEGER, INTENT(IN) :: istatus
         INTEGER, INTENT(OUT) :: ioerr
      end SUBROUTINE check
   end interface
   ! Define a new type to hold all details of a given configuration variable
   type config
      character(len=:), allocatable :: varname
      integer :: varid
      integer :: numdims
      integer :: numattrs
      integer :: xtype
      integer, dimension(:), allocatable :: dimlen, dimids
      real(myprec) :: scale
      character(len=:), allocatable :: units
   end type config
   ! Character variables of flexible length
   type globalat
      character(len=:), allocatable :: name
      character(len=:), allocatable :: val
   end type globalat
   type dimens
      integer :: length
      character(len=:), allocatable :: name
   end type dimens
   ! Remap NetCDF F90 variable types
   character(len=*), parameter, dimension(6) :: tipos = (/"NF90_BYTE  ", "&
   &NF90_CHAR  ", "NF90_SHORT ", "NF90_INT   ", "NF90_FLOAT ", "NF90_D&
   &OUBLE"/)
contains
   ! Routine to return NetCDF error codes
   SUBROUTINE check(istatus, ioerr)
      USE netcdf
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: istatus
      INTEGER, INTENT(OUT) :: ioerr
      IF (istatus /= nf90_noerr) THEN
         write (*, *) TRIM(ADJUSTL(nf90_strerror(istatus)))
         ioerr = istatus
      else
         ioerr = nf90_noerr
      END IF
   END SUBROUTINE check
end module mod_nc
!
! Define general variables to store configuration vales
!
module mod_nc_conf
   use mod_nc
   use netcdf
   type(config), dimension(:), allocatable :: conf
   type(dimens), dimension(:), allocatable :: mydims
   type(globalat), dimension(:), allocatable :: myglobatts
   real(myprec), dimension(:, :, :), allocatable :: r, v, fxyz
   real(double), dimension(:, :, :), allocatable :: stress_i
   real(myprec), dimension(:,:), allocatable :: u_pi, qc
   real(myprec), dimension(3, 1) :: org, cell, cell_a
   real(myprec) :: time, scale
   integer, dimension(:, :), allocatable :: idi, ity, imol
   character(len=NF90_MAX_NAME) :: nombre, nom_att
   character(len=3) :: spatial, cell_spatial
   character(len=5), dimension(1:3) :: cell_ang
   character(len=55) :: fname, fnamew, attr
   integer, dimension(nf90_max_var_dims) :: rhDimIds, rdimId
   integer :: vtype, nDims, nVars, nGlobalAtts, unlimDimID, step(1), ntypes, nmconf
   integer :: start(3), count(3), starti(2), counti(2),&
   & cstart(2), ccount(2), csstart(1), cscount(1), ln, ntm&
   &, status, varid, numatts, numdims, natoms, nconf_i, nlen, nwty
   integer, dimension(:), allocatable :: ndimid, nvarid, atypes, wtypes, orgty
   integer, dimension(:), allocatable :: counter, nct
end module mod_nc_conf


subroutine read_nc_cfg(ncid, ncstart, io, unit)
   ! Read netcdf configurations from file
   use mod_nc
   use netcdf
   use mod_nc_conf
   use mod_input, only : nsp
   implicit none
   integer, intent(in) :: ncid, ncstart
   integer, intent(in), optional :: unit
   integer, intent(out) :: io
   integer, dimension(:), allocatable :: tempty
   logical, save :: first = .true., typedefined = .false.
   integer :: i, j, tipo, iunit, k, ioerr, tmty(1)
   if (.not.present(unit)) then
      iunit = 6
   else
      iunit = unit
   end if
   ! 3 spatial coordinate
   count(1) = 3
   ! read in frame after frame
   count(3) = 1
   !
   starti(:) = 1
   starti(2) = ncstart
   counti(:) = 1
   start(1) = 1
   start(2) = 1
   start(3) = ncstart
   cstart(1) = 1
   cstart(2) = ncstart
   ccount(1) = 3
   ccount(2) = 1
   if (first) then
      Write (iunit, "(//'**** Using NetCDF library version ',A)") trim(nf90_inq_libvers())
      !
      ! Inquire dimensions and number of vars
      !
      call check(nf90_inquire(ncid, ndims, nvars, nglobalatts,&
      & unlimdimid), ioerr)
      write (iunit, "(// ' No. of dimensions ',i3,' No of vars ', i3,' No. of g&
      &lobal attibutes ', i3, ' Unlimited dimension ', i3)") ndims, nvars, nglobalatts,&
      & unlimdimid
      if (nglobalatts > 0) allocate (myglobatts(nglobalatts))
      do j = 1, nglobalatts
         call check(nf90_inq_attname(ncid, NF90_GLOBAL, j, nom_att), ioerr)
         nlen = len(trim(adjustl(nom_att)))
         allocate (character(len=nlen) :: myglobatts(j)%name)
         myglobatts(j)%name = trim(adjustl(nom_att))
         call check(nf90_inquire_attribute(ncid, NF90_GLOBAL, myglobatts(j)%name, xtype&
         &=tipo, len=nlen), ioerr)
         allocate (character(len=nlen) :: myglobatts(j)%val)
         call check(nf90_get_att(ncid, NF90_GLOBAL, nom_att, attr), ioerr)
         myglobatts(j)%val = trim(adjustl(attr))
         write (iunit, "(' - Global Attribute: ',2A25,' length',i3,' valor: ',A)") nom_att&
         &, tipos(tipo), nlen, myglobatts(j)%val
      end do
      allocate (conf(nvars))
      allocate (mydims(ndims))
      write (iunit, "(/' ** Readings dim id, name, length ',/)")
      do j = 1, ndims
         call check(nf90_inquire_dimension(ncid, j, name=nombre, len=&
         & nlen), ioerr)
         mydims(j)%length = nlen
         ln = len(trim(adjustl(nombre)))
         allocate (character(len=ln) :: mydims(j)%name)
         mydims(j)%name = trim(adjustl(nombre))
         write (iunit, "( 'dimension id',i2,' nombre ',A, ' length=',i)") j, mydims(j)%name, mydims(j)%length
         if (mydims(j)%name == "frame") nmconf = mydims(j)%length
         if (mydims(j)%name == "Voigt") voigt = mydims(j)%length
      end do
      write (iunit, "(/' ** Readings vars, type, no. dimensions, no. a&
      &ttrs ',/)")
   end if
   if (ncstart > nmconf) then
      write (*, "(' ** Warning : trying to read past last configuration no. ',i7 )") nmconf
      io = -1
      return
   end if
   do i = 1, nvars
      if (first) then
         call check(nf90_inquire_variable(ncid, i, xtype=tipo, name=nombre, ndims=&
         & numdims, natts=numatts), ioerr)
         conf(i)%xtype = tipo
         conf(i)%varid = i
         write (iunit, "('varid=',i2,' ',2A15,' ndims =',i2,' nattrs=', i2)") i&
         &, nombre, tipos(tipo), numdims, numatts
         ln = len(trim(adjustl(nombre)))
         allocate (character(len=ln) :: conf(i)%varname)
         conf(i)%varname(1:ln) = trim(adjustl(nombre))
         conf(i)%numdims = numdims
         conf(i)%numattrs = numatts
         allocate (conf(i)%dimlen(numdims))
         call check(nf90_inquire_variable(ncid, i, dimids=&
         & rhdimids(1:numdims)), ioerr)
         allocate (conf(i)%dimids(numdims))
         conf(i)%dimids(1:numdims) = rhdimids(1:numdims)
         do j = 1, numdims
            call check(nf90_inquire_dimension(ncid, rhdimids(j), len=&
            & nlen), ioerr)
            conf(i)%dimlen(j) = nlen
            write (iunit, "('  ** dim id =',i3,' length = ',i6)") rhdimids(j)&
            &, nlen
         end do
         do j = 1, numatts
            call check(nf90_inq_attname(ncid, i, j, nom_att), ioerr)
            call check(nf90_inquire_attribute(ncid, i, nom_att, xtype&
            &=tipo, len=nlen), ioerr)
            write (iunit, "('    · Attribute ',2A15,' length',i3)") nom_att&
            &, tipos(tipo), nlen
         end do
         if (numatts > 0) then
            call check(nf90_inq_attname(ncid, i, 1, nom_att), ioerr)
            call check(nf90_get_att(ncid, i, nom_att, attr), ioerr)
            ln = len(trim(adjustl(attr)))
            allocate (character(len=ln) :: conf(i)%units)
            conf(i)%units = trim(adjustl(attr))
         end if
         if (numatts == 2) then
            call check(nf90_inq_attname(ncid, i, 2, nom_att), ioerr)
            call check(nf90_get_att(ncid, i, nom_att, scale), ioerr)
            conf(i)%scale = scale
         else if (numatts > 2) then
            write (*, "(' ** Error: program cannot handles more than two attr&
            &ibutes for var ',A)") nombre
            stop
         end if
      end if

      if (conf(i)%varname == "coordinates") then
         if (first) allocate (r(conf(i)%dimlen(1), conf(i)%dimlen(2), 1))
         natoms = conf(i)%dimlen(2)
         nconf_i = conf(i)%dimlen(3)
         start(3) = ncstart
         count(2) = natoms
         call check(nf90_get_var(ncid, i, r, start, count), ioerr)
         do k = 1, 3
            if (cell(k, 1) == 0) then
               periodic(k) = .false.
               cell(k, 1) = max(abs(maxval(r(k, 1:natoms, 1)) - minval(r(k&
               &, 1:natoms, 1))) + 10.0, 2*abs(org(k, 1)))
            end if
         end do
      else if (conf(i)%varname == "velocities") then
         ex_vel = .true.
         if (first) allocate (v(conf(i)%dimlen(1), conf(i)%dimlen(2), 1))
         start(3) = ncstart
         count(2) = natoms
         call check(nf90_get_var(ncid, i, v, start, count), ioerr)
      else if (conf(i)%varname == "forces") then
         ex_force = .true.
         if (first) allocate (fxyz(conf(i)%dimlen(1), conf(i)%dimlen(2), 1))
         start(3) = ncstart
         count(2) = natoms
         call check(nf90_get_var(ncid, i, fxyz, start, count), ioerr)
      else if (conf(i)%varname == "c_"//trim(adjustl(ener_name))) then
         run_thermo = .true.
         if (first) then
            allocate (u_pi(conf(i)%dimlen(1), 1))
         endif
         counti(1) = conf(i)%dimlen(1)
         call check(nf90_get_var(ncid, i, u_pi, starti, counti), ioerr)
      else if (conf(i)%varname == "c_"//trim(adjustl(press_name))) then
         ex_stress = .true.
         if (first) allocate (stress_i(conf(i)%dimlen(1), conf(i)%dimlen(2), 1))
         count(1) = conf(i)%dimlen(1)
         start(3) = ncstart
         count(2) = natoms
         call check(nf90_get_var(ncid, i, stress_i, start, count), ioerr)
      else if (conf(i)%varname == "type") then
         typedefined = .true.
         if (first) allocate (ity(conf(i)%dimlen(1), 1))
         counti(1) = conf(i)%dimlen(1)
         call check(nf90_get_var(ncid, i, ity, starti, counti), ioerr)
      else if (conf(i)%varname == "mol") then
         ex_mol = .true.
         if (first) allocate (imol(conf(i)%dimlen(1), 1))
         counti(1) = conf(i)%dimlen(1)
         call check(nf90_get_var(ncid, i, imol, starti, counti), ioerr)
      else if (conf(i)%varname == "q") then
         ex_qc = .true.
         if (first) allocate (qc(conf(i)%dimlen(1), 1))
         counti(1) = conf(i)%dimlen(1)
         call check(nf90_get_var(ncid, i, qc, starti, counti), ioerr)
      else if (conf(i)%varname == "id") then
         if (first) allocate (idi(conf(i)%dimlen(1), 1))
         counti(1) = conf(i)%dimlen(1)
         ! Reasign id in a contiguos fashion
         do j = 1, counti(1)
            idi(j, 1) = j
         end do
         call check(nf90_get_var(ncid, i, idi, starti, counti), ioerr)
      else if (conf(i)%varname == "cell_origin") then
         ccount(2) = 1
         call check(nf90_get_var(ncid, i, org, cstart, ccount), ioerr)
      else if (conf(i)%varname == "cell_lengths") then
         ccount(2) = 1
         call check(nf90_get_var(ncid, i, cell, cstart, ccount), ioerr)
      else if (conf(i)%varname == "cell_angles") then
         ccount(2) = 1
         call check(nf90_get_var(ncid, i, cell_a, cstart, ccount), ioerr)
      else if (conf(i)%varname == "cell_spatial") then
         call check(nf90_get_var(ncid, i, cell_spatial), ioerr)
      else if (conf(i)%varname == "cell_angular") then
         call check(nf90_get_var(ncid, i, cell_ang), ioerr)
      else if (conf(i)%varname == "spatial") then
         call check(nf90_get_var(ncid, i, spatial(1:3)), ioerr)
      else if (conf(i)%varname == "time") then
         csstart(1) = ncstart
         cscount(1) = 1
         call check(nf90_get_var(ncid, i, step, csstart, cscount), ioerr)
      end if
      if (ioerr /= nf90_noerr) then
         ioerr = -1
         return
      end if
   end do
   if (typedefined) then
      if (first) then
         allocate(tempty(100),wtypes(nsp))
         tempty(:) = -1
         tempty(1) = ity(1,1)
         ntypes = 1
         do i = 1, natoms
            if (.not. any(tempty(1:ntypes) == ity(i,1))) then
               ntypes = ntypes+1
               tempty(ntypes) = ity(i,1)
            endif
         end do
         write (iunit, "(/' ** Number of atoms types =',i3)") ntypes

         allocate (atypes(ntypes),orgty(ntypes))
         orgty(1:ntypes) = tempty(1:ntypes)
         if (nsp == ntypes) wtypes = orgty
         deallocate(tempty)
         atypes(:) = 0
         do i = 1, natoms
            ! Remap types
            tmty = findloc(orgty(1:ntypes),ity(i,1))
            ity(i,1) = tmty(1)
            atypes(ity(i, 1)) = atypes(ity(i, 1)) + 1
         end do
         do i = 1, ntypes
            write (iunit, "(/'       Number of atoms of types ',i2,' (',i2,')=',i8)") i, orgty(i),atypes(i)
         end do
      else
         !
         ! Atoms types must be remapped every configuration
         do i = 1, natoms
            ! Remap types
            tmty = findloc(orgty(1:ntypes),ity(i,1))
            ity(i,1) = tmty(1)
         end do
      end if
   else
      print *, " **** Warning no atom types defined in trajectory file"
   end if
   ! if (iunit .ne. 6 .and. first) close (iunit)
   first = .false.
end subroutine read_nc_cfg

subroutine reset_nmol(nmoln)
   use mod_nc_conf, only: ity_in => ity, natoms, ntypes, wtypes, orgty
   integer, intent(INOUT) :: nmoln
   integer :: i, index
   ! remap data to take into account selected atoms
   ! calculate number of atoms to consider
   index = 0
   do i = 1, natoms
      if (any(wtypes == orgty(ity_in(i,1))))then
         if(count(wtypes==orgty(ity_in(i,1)))>1) then
            write(*,"(' !!*** Error: type',i2,' appears more than onces in selection ')")ity_in(i,1)
            stop
         endif
         index=index+1
      endif
   enddo
   Nmoln = index
   if (Nmoln >0) then
      write(*,"(' !!*** Only ',i7,' atoms selected: Nmol reset ...')") nmoln
   else
      write(*,"(' !!*** Error: selected types not in netcdf trajectory file: ',8i3)")wtypes
      stop
   endif
end subroutine reset_nmol

subroutine trans_ncdfinput()
   use mod_nc_conf, only: org, cell_in => cell, r_in => r, v_in => v, f_in => fxyz, &
      ity_in => ity, nstep_in => step, natoms, ntypes, wtypes, u_pi, stress_i, qc, &
      nct, counter
   use mod_common, only: vel, r, force, cell, sidel, side, volumen, itype, bscat, tunit, &
      ntype, masa, qcharge, nstep, vector_product, nmol, ex_vel, ex_force, ex_qc, &
      tuniti, side2, stress, u_p, voigt, run_thermo, ex_stress, cntch, ncharge, chgh
   use mod_input, only: ndim, mat, bsc, rcrdf, nsp, charge
   implicit none
   integer :: i, j, it, index, ipch(1)
   logical :: pass = .true., compcharge=.true., first = .true.
   if (.not.allocated(nct)) allocate(nct(nsp))
   if (.not.allocated(counter)) allocate(counter(nsp))
   nstep = nstep_in(1)
   sidel(:) = cell_in(:, 1)
   side = Minval(sidel(1:ndim))
   ! secure rcrdf to be less that half the simulation box
   rcrdf = min(rcrdf,side/2)
   if (rcrdf > 0.0) then
      side2 = rcrdf**2
   else
      side2 = (side/2)**2
   end if
   volumen = 1.0
   do i = 1, ndim
      cell(i, i) = cell_in(i, 1)
      volumen = volumen*cell(i, i)
   end do

   !
   ! Remap coordinates so as to have first ntype(1) particles in contiguous positions, followed
   ! by ntype(2) particles and so on .. and transform from netcdf format to local vars
   !
   ntype(:) = 0
   do i = 1, natoms
      ntype(ity_in(i, 1)) = ntype(ity_in(i, 1)) + 1
   end do
   nct(1) = 0
   do i = 2, nsp
      nct(i) = nct(i - 1) + ntype(i - 1)
   end do
   counter(:) = nct(:)
   do i = 1, natoms
      it = ity_in(i, 1)
      counter(it) = counter(it) + 1
      j = counter(it)
      !
      ! NOTE: when ndim=2, z component of r_in,v_in is discarded
      !
      if (ex_vel) vel(1:ndim, j) = v_in(1:ndim, i, 1)*tunit/tuniti
      if (ex_force) force(1:ndim, j) = f_in(1:ndim, i, 1)
      if (run_thermo) u_p(j) = u_pi(i,1)
      if (ex_stress) stress(1:voigt,j)=stress_i(1:voigt,i,1)
      r(1:ndim, j) = r_in(1:ndim, i, 1)
      ! Coordinates are folded back into the simulation cell under PBC
      r(1:ndim, j) = r(1:ndim, j) - sidel(1:ndim)*int(r(1:ndim, j)&
      &/sidel(1:ndim))
      itype(j) = it
      if (ex_qc) then
         qcharge(j) = qc(j,1)
         if (compcharge) then
            if (pass) then
               cntch(:) = 0
               cntch(1) = 1
               chgh(1) = qcharge(j)
               ncharge = 1
               pass = .false.
            else
               if (any(chgh(1:ncharge)==qcharge(j))) then
                  ipch = findloc(chgh(1:ncharge),qcharge(j))
                  cntch(ipch(1)) = cntch(ipch(1))+1
               else
                  ncharge = ncharge+1
                  chgh(ncharge) = qcharge(j)
                  cntch(ncharge) = cntch(ncharge)+1
               endif
            endif
         endif
      endif
   end do
   compcharge = .false.
   do i = 1, nsp
      masa(nct(i) + 1:nct(i) + ntype(i)) = mat(i)
      bscat(nct(i) + 1:nct(i) + ntype(i)) = bscat(i)
      if(ex_qc) charge(i) = sum(qcharge(nct(i) + 1:nct(i) + ntype(i)))/ntype(i)
   end do

   compcharge = .false.


end subroutine trans_ncdfinput

subroutine select_ncdfinput()
   ! Select atoms in wtypes from netcdf file and remap coordinates in species order
   use mod_nc_conf, only: org, cell_in => cell, r_in => r, v_in => v, f_in => fxyz, &
      ity_in => ity, nstep_in => step, natoms, ntypes, wtypes, u_pi, stress_i, orgty, qc, &
      nct, counter
   use mod_common, only: vel, r, force, cell, sidel, side, volumen, itype, bscat, tunit, &
      ntype, masa, nstep, vector_product, nmol, ex_vel, ex_force, ex_qc, &
      tuniti, side2, u_p, stress, voigt, run_thermo, ex_stress, qcharge, chgh, ncharge, cntch
   use mod_input, only: ndim, mat, bsc, rcrdf, nsp, charge
   implicit none
   integer :: i, j, it(1), index, ipch(1)
   logical :: pass = .true., compcharge = .true., first=.true.
   if (.not.allocated(nct)) allocate(nct(nsp))
   if (.not.allocated(counter)) allocate(counter(nsp))
   nstep = nstep_in(1)
   sidel(:) = cell_in(:, 1)
   side = Minval(sidel(1:ndim))
   ! secure rcrdf to be less that half the simulation box
   rcrdf = min(rcrdf,side/2)
   if (rcrdf > 0.0) then
      side2 = rcrdf**2
   else
      side2 = (side/2)**2
   end if
   volumen = 1.0
   do i = 1, ndim
      cell(i, i) = cell_in(i, 1)
      volumen = volumen*cell(i, i)
   end do

   !
   ! Remap coordinates so as to have first ntype(1) particles in contiguous positions, followed
   ! by ntype(2) particles and so on .. and transform from netcdf format to local vars
   !
   ntype(:) = 0
   do i = 1, natoms
      if (any(wtypes == orgty(ity_in(i,1)))) then
         it = findloc(wtypes,orgty(ity_in(i,1)))
         ntype(it(1)) = ntype(it(1)) + 1
      endif
   end do
   nct(1) = 0
   do i = 2, nsp
      nct(i) = nct(i - 1) + ntype(i - 1)
   end do
   counter(:) = nct(:)
   do i = 1, natoms
      if (any(wtypes == orgty(ity_in(i,1)))) then
         it = findloc(wtypes,orgty(ity_in(i, 1)))
         counter(it(1)) = counter(it(1)) + 1
         j = counter(it(1))
         !
         ! NOTE: when ndim=2, z component of r_in,v_in is discarded
         !
         if (ex_vel) vel(1:ndim, j) = v_in(1:ndim, i, 1)*tunit/tuniti
         if (ex_force) force(1:ndim, j) = f_in(1:ndim, i, 1)
         if (run_thermo) u_p(j) = u_pi(i,1)
         if (ex_stress) stress(1:voigt,j)=stress_i(1:voigt,i,1)
         r(1:ndim, j) = r_in(1:ndim, i, 1)
         if (ex_qc) then
            qcharge(j) = qc(i,1)
            if (compcharge) then
               if (pass) then
                  cntch(:) = 0
                  cntch(1) = 1
                  chgh(1) = qcharge(j)
                  ncharge = 1
                  pass = .false.
               else
                  if (any(chgh(1:ncharge)==qcharge(j))) then
                     ipch = findloc(chgh(1:ncharge),qcharge(j))
                     cntch(ipch(1)) = cntch(ipch(1))+1
                  else
                     ncharge = ncharge+1
                     chgh(ncharge) = qcharge(j)
                     cntch(ncharge) = cntch(ncharge)+1
                  endif
               endif
            endif
         endif
         ! Coordinates are folded back into the simulation cell under PBC
         r(1:ndim, j) = r(1:ndim, j) - sidel(1:ndim)*int(r(1:ndim, j)&
         &/sidel(1:ndim))
         itype(j) = it(1)
      endif
   end do
   compcharge = .false.
   do i = 1, nsp
      masa(nct(i) + 1:nct(i) + ntype(i)) = mat(i)
      bscat(nct(i) + 1:nct(i) + ntype(i)) = bscat(i)
      if(ex_qc) charge(i) = sum(qcharge(nct(i) + 1:nct(i) + ntype(i)))/ntype(i)
   end do
end subroutine select_ncdfinput
