!===============================================================================
! Module: lmp_util
!===============================================================================
! Purpose:
!   Utility routines for reading LAMMPS trajectory files in various formats.
!   Supports DL_POLY HISTORY format and LAMMPS dump formats for converting
!   to NetCDF.
!
! Key Functionality:
!   - Reads DL_POLY HISTORY trajectory files
!   - Reads LAMMPS dump files (standard and custom formats)
!   - Extracts atomic positions, velocities, forces
!   - Handles different LAMMPS keytrj values (0,1,2)
!   - Processes simulation cell information
!   - Skips unwanted configurations efficiently
!
! Main Subroutines:
!   ReadDummyCfg()    - Skips configurations without storing data
!   ReadSingleCfg()   - Reads one complete configuration
!   ReadCfgH()        - Reads DL_POLY HISTORY format
!
! Supported Formats:
!   itrj = 0: DL_POLY HISTORY format
!   itrj = 1: LAMMPS dump format (positions only)
!   itrj = 2: LAMMPS dump custom format (with velocities/forces)
!
! Variables:
!   keytrj - LAMMPS trajectory type:
!            0: positions only
!            1: positions + velocities
!            2: positions + velocities + forces
!   imcon  - Periodic boundary conditions flag
!   cell   - Simulation cell matrix
!
! Notes:
!   - Used by lmptrj2nc for format conversion
!   - Handles both 2D and 3D systems
!   - Automatically detects cell periodicity
!===============================================================================
Module lmp_util
   Use comun
   Implicit none
   contains 
Subroutine ReadDummyCfg(iunit,itrj)
  Integer, Intent(IN) :: iunit, itrj
  Integer :: i, io, iatm, idum1, id
  Real (kind=8) :: dumx, dumy, dumz, rlow(ndim), rup(ndim)
  If (itrj == 0) Then
     Read(iunit,'(a8,4i10,f12.6)',iostat=io) dumy, nstep, natms,&
          & keytrj, imcon, tstep
     If (io < 0) Then
        Print *, ' Error skipping initial configurations ..'
        Stop
     End If
     If (imcon.Gt.0) Then
        Do i=1, ndim
           Read(iunit,*) cell(1:ndim,i)
        End Do
     Endif
     iatm = 0
     Do  While (iatm.Lt.natms)
        Read(iunit,"(a8,i10,2f12.6)") atnam, iatm, weight, charge
        If (keytrj.Eq.0) Then
           Read(iunit,*) dumx,dumy, dumz
        Else If (keytrj.Eq.1) Then
           Read(iunit,*) dumx,dumy, dumz
           Read(iunit,*) dumx,dumy, dumz
        Else If (keytrj.Eq.2) Then
           Read(iunit,*) dumx,dumy, dumz
           Read(iunit,*) dumx,dumy, dumz
           Read(iunit,*) dumx,dumy, dumz
        Endif
     End Do
  Else If (itrj == 1) Then
     read(iunit,*) natms
     if (ndim == 3) then
        read(iunit,'(7x,i10,10x,3f15.7)')nstep,(cell(i,i),i=1,ndim)
     else
        read(iunit,'(7x,i10,10x,3f15.7)')nstep,dumy,(cell(i,i),i=1,ndim)
     endif
     do i = 1, natms
        Read(iunit,*) idum1
     end do
  Else If (itrj == 2) Then
     Read(iunit,'(1x)')
     Read(iunit,*)nstep
     Read(iunit,'(1x)')
     Read(iunit,*)natms
     Read(iunit,'(1x)')
     cell(:,:) = 0.0
     Do i=1, ndim
        Read(iunit,*)rlow(i),rup(i)
        cell(i,i) = rup(i)-rlow(i)
     End Do

     Read(iunit,'(1x)')
     Do i=1, natms
        Read(iunit,*) idum1
     End Do
  Else If (itrj == 3) Then
     read(iunit,*) natms
     read(iunit,'(1x)')
     do i = 1, natms
        Read(iunit,*) idum1
     end do
  End If
End Subroutine ReadDummyCfg

Subroutine ReadCfg(iunit,io,iconf,pconf,istart,itrj)
  Integer, Intent(IN) :: iunit,iconf,pconf,istart,itrj
  Integer, Intent(OUT) :: io
  Integer :: i, l, j, iatm, idum1, id, ityp, imt, idi
  Real(kind=8) :: dumy, rlow(ndim), rup(ndim), rtemp(3), veltemp(3), forcetemp(3)
  logical, save :: first=.true.
  io = 0
!  print *, " Reading configuration ", itrj, natms
  If (itrj == 0) Then
     Read(iunit,'(a8,4i10,f12.6)',iostat=io) dumy, nstep, natms, keytrj, imcon, tstep
     If (io < 0) Return
     If (imcon.Gt.0) Then
        Do i=1, ndim
           Read(iunit,*) cell(1:ndim,i)
        End Do
     Endif
     iatm = 0
     l = 0
     ntype(:) = 0
     Do While (iatm.Lt.natms)
        l = l + 1
        Read(iunit,"(a8,i10,2f12.6)") atnam, iatm, weight, charge
        Do i=1, nsp
           If (atnam == atoms(i)) Then
              masa(l) = weight
              itype(l) = i
              ntype(i) = ntype(i)+1
              Exit
           End If
        End Do
        If (keytrj.Eq.0) Then
           Read(iunit,*) R(iatm,1:ndim)
        Else If (keytrj.Eq.1) Then
           Read(iunit,*) R(iatm,1:ndim)
           Read(iunit,*) vel(iatm,1:ndim)
        Else If (keytrj.Eq.2) Then
           Read(iunit,*) R(iatm,1:3)
           Read(iunit,*) vel(iatm,1:ndim)
           Read(iunit,*) force(iatm,1:ndim)
        Endif
!
     Enddo
     tmol(:) = 1
  Else If (itrj ==1) Then
     read(iunit,*,iostat=io) natms
     if (ndim == 3) then
        read(iunit,'(7x,i10,10x,3f15.7,7x,f10.7)')nstep,(cell(i,i),i=1,ndim), tstep
     else
        read(iunit,'(7x,i10,10x,3f15.7,7x,f10.7)')nstep,dumy,(cell(i,i),i=1,ndim), tstep
     endif
     ntype(:) = 0
     do i = 1, natms
        if (ndim == 3) then
           If (keytrj.Eq.0) Then
              Read(iunit,*) itype(i), R(i,1:ndim)
           Else If (keytrj.Eq.1) Then
              Read(iunit,*) itype(i), R(i,1:ndim), Vel(i,1:ndim)
           Else If (keytrj.Eq.2) Then
              Read(iunit,*) itype(i), R(i,1:ndim), Vel(i,1:ndim), force(i,ndim)
           Endif
        else
           If (keytrj.Eq.0) Then
              Read(iunit,*) itype(i), dumy, R(i,1:ndim)
           Else If (keytrj.Eq.1) Then
              Read(iunit,*) itype(i), dumy, R(i,1:ndim), dumy, Vel(i,1:ndim)
           Else If (keytrj.Eq.2) Then
              Read(iunit,*) itype(i), dumy, R(i,1:ndim), dumy, Vel(i,1:ndim), dumy, force(i,ndim)
           Endif
        endif
        masa(i) = mat(itype(i))
        ntype(itype(i)) = ntype(itype(i))+1
!
     end do
  Else If (itrj == 2) Then
     !
     ! LAMMPS trajectory files must be ordered regarding atom id's. This must be taken care of.
     ntype(:) = 0
     Read(iunit,'(1x)',iostat=io)
     Read(iunit,*)nstep
     Read(iunit,'(1x)')
     Read(iunit,*)natms
!     Print *, ' Reading LAMMPS trajectory file, sorted input ',natms,' atoms '
     Read(iunit,'(1x)')
     cell(:,:) = 0.0
     Do i=1, ndim
        Read(iunit,*)rlow(i),rup(i)
        cell(i,i) = rup(i)-rlow(i)
     End Do
     if (ndim==3) then
        Read(iunit,'(1x)')
     elseif (ndim==2) then
        Read(iunit,'(/1x)')
     endif
     id = 0
     Do i=1, natms
        id = id+1
        If (keytrj.Eq.0) Then
           Read(iunit,*) idi, ityp,  imt, Rtemp(1:3)
        Else If (keytrj.Eq.1) Then
           Read(iunit,*) idi, ityp,  imt, Rtemp(1:3), Veltemp(1:3)
        Else If (keytrj.Eq.2) Then
           Read(iunit,*) idi, ityp,  imt, Rtemp(1:3), Veltemp(1:3), forcetemp(1:3)
        Endif
        if (ndim == 3) then
           R(id,1:3) = Rtemp(1:3)
           Vel(id,1:3) = Veltemp(1:3)
           Force(id,1:3) = forcetemp(1:3)
        else
           R(id,1:2) = Rtemp(1:2)
           Vel(id,1:2) = Veltemp(1:2)
           Force(id,1:2) = forcetemp(1:2)
        endif
        !
        ! defaul lmp trajectory file uses box units (better use custom dump)
        !
        itype(id) = ityp
        masa(id) = mat(itype(id))
        bscat(id) = bsc(itype(id))

        ntype(itype(id)) = ntype(itype(id))+1
     End Do

  End If
!
! Orthogonal cells
!
  Forall (i=1:ndim) sidel(i) = cell(i,i)
  side = Minval(sidel(1:ndim))
!
! tetragonal cells
!
  If (ndim == 3) Then
     volumen = cell(1,1)*cell(2,2)*cell(3,3)
  Else
     volumen = cell(1,1)*cell(2,2)
  End If
  !
End Subroutine ReadCfg
