!===============================================================================
! Module: mod_log
!===============================================================================
! Purpose:
!   Manages output logging and reporting of simulation results during
!   trajectory analysis. Provides formatted output for thermodynamics,
!   cluster properties, and order parameters.
!
! Key Functionality:
!   - Initializes and manages log file I/O
!   - Periodic output of instantaneous properties
!   - Thermodynamic data reporting (energy, temperature, pressure)
!   - Cluster analysis statistics output
!   - Order parameter results formatting
!   - Performance monitoring (CPU time per configuration)
!   - Progress tracking during trajectory processing
!
! Main Subroutines:
!   log_init()         - Opens output log file
!   log_clear()        - Closes log file
!   print_output()     - Periodic output of instantaneous properties
!   print_clusinfo()   - Cluster analysis statistics
!   print_order()      - Order parameter results
!   printPotEngCl()    - Cluster potential energy distributions
!
! Output Files:
!   - Log file: specified in input namelist
!   - thermo_run.dat: instantaneous thermodynamic properties
!
! Notes:
!   - Output frequency controlled by nprint parameter
!   - Supports both 'real' and 'lj' unit systems
!   - Energy and pressure units automatically handled
!===============================================================================
module mod_log
   use mod_precision
   use mod_common
   use mod_input
   use mod_nc
   use mod_nc_conf
   use mod_clusters
   use mod_order, only: norder, orderp, avorder, avorder_cos, avorder_sin, &
        cluster_order_cos, cluster_order_sin, avcluster_order, avcluster_order_cos, &
        avcluster_order_sin, atomic_order_cos, atomic_order_sin, atomic_ql, rhoorderav, &
        ordercumav
   use cudafor
   use mod_thermo, only : engclus, engclpa
   implicit none
   integer :: io_log_file=777
contains
   subroutine log_init()
      use mod_input
      open (unit=io_log_file, file=log_output_file)
      call header(io_log_file)
      call header(6)
      call print_active_modules(io_log_file)
      call print_active_modules(6)
   end subroutine log_init

   subroutine log_close()
      close (io_log_file)
   end subroutine log_close
  
   subroutine print_active_modules(unit)
      integer, intent(in) :: unit
      integer :: i
      write(unit,'(//" *** Active analysis modules:")')
      if (unit == 6) then 
         write(unit,'(a,90("_"),a)') char(27)//'[33m', char(27)//'[0m'
      else
         write(unit,'(90("_"))') 
      endif
      do i = 1, 7
         if (rdf_sq_cl_dyn_sqw_conf_ord(i) == .true.) then
         select case (i)
            case (1)
               write(unit,'(" ··· Flow control: ",A," module will be executed !")') "RDF"
            case (2)
               write(unit,'(" ··· Flow control: ",A," module will be executed !")') "S(q)"
            case (3)
               write(unit,'(" ··· Flow control: ",A," module will be executed !")') "Cluster analysis"
            case (4)
               write(unit,'(" ··· Flow control: ",A," module will be executed !")') "Dynamics/Z(w)"
            case (5)
               write(unit,'(" ··· Flow control: ",A," module will be executed !")') "F(q,t)/S(q,w)"
            case (6)
               write(unit,'(" ··· Flow control: ",A," module will be executed !")') "Confinement"
            case (7)
               write(unit,'(" ··· Flow control: ",A," module will be executed !")') "Order parameters"
         end select
         endif
      end do
      if (unit == 6) then 
         write(unit,'(a,90("_"),a)') char(27)//'[33m', char(27)//'[0m'
      else
         write(unit,'(90("_"))') 
      endif
   end subroutine print_active_modules
   
   subroutine header(unit)
      implicit none
      integer, intent(in) :: unit
      if(unit == 6) write(unit,'(A)')char(27)//'[33m'
      write(unit,"(/80('*')/'*',78(' '),'*')")
      write(unit,"('*    Program trj_analysis: analyzing LAMMPS trajectory in NETCDF format',t80,'*')")
      write(unit,"('*',t80,'*')")
      write(unit,"('*    Using GPU with CUDA nvfortran/nvcc >= 25.9/13.0',t80,'*')")
      write(unit,"('*',t80,'*')")
      write(unit,"('*    Version 1.3 March, 2026',,t80,'*')")
      write(unit,"('*',78(' '),'*'/80('*')/)")
      if(unit == 6) write(unit,'(A)') char(27)//'[0m'   
   end subroutine header


   subroutine print_output(iconf)
      integer, parameter :: nther=10
      integer, intent(in) :: iconf
      logical, dimension(nther) :: mascara
      character*15, dimension(nther) :: title =(/"     T(K)","KE (Kcal/mol)","PE Kcal/mol",&
         " Pressure (bar)","      Pxx(bar)","    Pyy(bar)","     Pzz(bar)",&
         "      Pxy(bar)","    Pxz(bar)","     Pyz(bar)"/)
      character*15, dimension(nther) :: titleeV =(/"     T(K)","KE (eV)","PE (eV)",&
         " Pressure (bar)","      Pxx(bar)","    Pyy(bar)","     Pzz(bar)",&
         "      Pxy(bar)","    Pxz(bar)","     Pyz(bar)"/)
      character*15, dimension(nther) :: titlelj =(/"     T(K)","KE/epsilon","PE/epsilon",&
         " Psigma³/epsilon","      Pxxsigma³/epsilon","    Pyysigma³/epsilon","     Pzzsigma³/epsilon",&
         "      Pxysigma³/epsilon","    Pxzsigma³/epsilon","     Pyzsigma³/epsilon"/) 
      real(myprec) :: thermo_q(nther), Tfcl
      integer :: i, j
      ! Initialize thermodynamic array and select available quantities
      thermo_q(:) = 0
      mascara(1) = ex_vel   ! Temperature and kinetic energy
      mascara(2) = ex_vel
      mascara(3) = run_thermo   ! Potential energy
      mascara(4:10) = ex_stress  ! Pressure and stress components
      if (ex_vel) then
         thermo_q(1) = temperature
         thermo_q(2) = kelvinfactor*ekin*(aunit/tunit)**2/Rgas
      endif
      if (run_thermo) thermo_q(3) = epot
      if (ex_stress) then
         thermo_q(4) = pressure
         thermo_q(5:10) = pxyz(1:6)
      endif
      ! Write instantaneous thermodynamics to file
      if (sum(mascara)) then
         if (iconf == 1) then
            open(1000, file="thermo_run.dat")
            if (tunits == 'lj') then
               write(1000,"('#    Conf  ',16a15)")pack(titlelj(1:nther),mascara)
            else if (tunits == 'picosecond') then
               write(1000,"('#    Conf  ',16a15)")pack(titleeV(1:nther),mascara)
            else     
               write(1000,"('#    Conf  ',16a15)")pack(title(1:nther),mascara)
            endif
         endif
         write(1000,"(i9,12f15.5)")nstep, pack(thermo_q(1:nther),mascara)
      endif
       if (run_clusters.and.Iconf==1) then
               Write (*, "(/ ' *** Clusters >= minPts=',i3,' particles being analyzed.  '/)") minPts
         endif
      ! Print periodic progress update to console
      If (Mod(Iconf - 1, nprint) .Eq. 0) Then
         call cpu_time(cpu1)
         write (*,"(a,90('_'),a)") char(27)//'[31m', char(27)//'[0m' 
         if (tunits == 'lj') then
               Write (*, "(' ** Working on MD step no. ',i10,' time* =',f12.3,&
               & ' compute time per conf.=',f7.2,' s:'&
               & )") nstep, nstep*tstep, (cpu1 - cpu0)/nprint
         else
            ! for real and metal units, print time in ns for better readability
               Write (*, "(' ** Working on MD step no. ',i10,' time =',f10.5,&
               & ' ns, compute time per conf.=',f7.2,' s:')") nstep, nstep*tstep/1000.0, (cpu1 - cpu0)/nprint
         endif
         write (*,"(a,90('_'),a)") char(27)//'[31m', char(27)//'[0m' 
  
         cpu0 = cpu1
         if (run_thermo) then
            if (tunits == 'lj') then
               write (*, "(' ** Potential energy/epsilon=',f15.4,' Per atom=',f15.4)") epot, epotperatom
               write (*, "(' ** Average potential energy/epsilon=',f15.4, ' Per atom=',f15.4)") epotav/Iconf, epotav/(Iconf*nmol)
            else if (tunits == 'picosecond') then
               write (*, "(' ** Potential energy=',f15.4,' eV, Per atom=',f15.4,' eV')") epot, epotperatom
               write (*, "(' ** Average potential energy=',f15.4,' eV, Per atom=',f15.4,' eV')") epotav/Iconf, epotav/(Iconf*nmol)
            else
               write (*, "(' ** Potential energy=',f15.4,' Kcal/mol, Per atom=',f15.4,'Kcal/mol')") epot, epotperatom
               write (*, "(' ** Average potential energy=',f15.4,' Kcal/mol, Per atom=',f15.4,'Kcal/mol')") epotav/Iconf, epotav/(Iconf*nmol)
            endif

            if (run_clusters) then
               if (tunits == 'lj') then
                  write (*, "(' ** Average intracluster potential energy/epsilon=',f15.4,' Peratom =',f15.4)") engclus, engclpa
               else
                  write (*, "(' ** Average intracluster potential energy=',f15.4,' Kcal/mol'&
                  ,' Peratom =',f15.4,' Kcal/mol')") engclus, engclpa
               endif
            endif
         end if
         if (ex_stress) then
            if (tunits == 'lj') then
               write (*, "(' ** Pressure*sigma**3/epsilon =',f15.4,' Average =',f15.4)") pressure, pressav/Iconf 
            else
               write (*, "(' ** Pressure =',f15.4,' bar, Average =',f15.4' bar')") pressure, pressav/Iconf !kcal_a3_to_bar*pressure, kcal_a3_to_bar*pressav/Iconf
            endif
         endif
         ! Kinetic energy and temperature
         If (ex_vel) then
            if (tunits=='lj') then
               write (*, "(' ** Kinetic energy/epsilon=',f15.4,' average=',f15.4)") &
                  ekin/natoms, ecaver/(natoms*Iconf)
            else if (tunits == 'picosecond') then
               write (*, "(' ** Kinetic energy=',f15.4,' eV, average=',f15.4,' eV')") &
                  kelvinfactor*ekin*(aunit/tunit)**2/Rgas, &
                  kelvinfactor*ecaver*(aunit/tunit)**2/Rgas/Iconf
            else
               write (*, "(' ** Kinetic energy=',f15.4,' Kcal/mol, average=',f15.4,'Kcal/mol')") &
                  kelvinfactor*ekin*(aunit/tunit)**2/Rgas, kelvinfactor*ecaver*(aunit/tunit)**2/Rgas/Iconf
            
            endif 
            if (run_clusters) then
               if (tunits == 'lj') then
                  write (*, "(' ** Cluster kinetic energy/epsilon=',f15.4,' average=',f15.4)") &
                     ekincl, ekclaver/Iconf
                  write (*, "(' ** Internal cluster kinetic energy/epsilon=',f15.4,' average=',f15.4)") &
                     ekincls, ekinclsav/Iconf
               elseif (tunits == 'picosecond') then
                  write (*, "(' ** Cluster kinetic energy=',f15.4,' eV, average=',f15.4,' eV')") &
                     kelvinfactor*ekincl*(aunit/tunit)**2/Rgas, &
                     kelvinfactor*ekclaver*(aunit/tunit)**2/Rgas/Iconf
                  write (*, "(' ** Internal cluster kinetic energy=',f15.4,' eV, average=',f15.4,' eV')") &
                     kelvinfactor*ekincls*(aunit/tunit)**2/Rgas, &
                     kelvinfactor*ekinclsav*(aunit/tunit)**2/Rgas/Iconf
               else
                   write (*, "(' ** Cluster kinetic energy/epsilon=',f15.4,' Kcal/mol average=',f15.4,' Kcal/mol')") &
                     kelvinfactor*ekincl*(aunit/tunit)**2/Rgas, kelvinfactor*ekclaver*(aunit/tunit)**2/Rgas/Iconf
                  write (*, "(' ** Internal cluster kinetic energy=',f15.4,'&
                  & Kcal/mol, average=',f15.4,'Kcal&
                  &/mol')") kelvinfactor*ekincls*(aunit/tunit)**2/Rgas,&
                  & kelvinfactor*ekinclsav*(aunit/tunit)**2/Rgas/Iconf
               endif
               if (maxcln>1) then
                  Tfcl = (maxcln-1)*ndim
               else
                  Tfcl = maxcln*ndim
               endif
               if (tunits == 'lj') then
                    Write (*, "(' ** Average cluster k_b*temperature/epsilon =',f10.2&
                  &)") 2*ekclaver/(Tfcl*Iconf)
               else if (tunits == 'picosecond') then
                  Write (*, "(' ** Average cluster temperature =',f10.2&
                  &,' K')") 2*ekclaver*(aunit/tunit)**2/(Tfcl*Rgas*Iconf)
               else
                  Write (*, "(' ** Average cluster temperature =',f10.2&
                  &,' K')") 2*ekclaver*(aunit/tunit)**2/(Tfcl*Rgas*Iconf)
               endif
            end if
            if (tunits == 'lj') then
               write (*, "(' ** Average k_b*temperature/epsilon =',f10.2)") 2*ecaver/(Tfact*Iconf)
            else
               write (*, "(' ** Average temperature =',f10.2&
               &,' K')") 2*ecaver*(aunit/tunit)**2/(Tfact*Rgas*Iconf)
            endif
         endif
         if (tunits == 'lj') then
            Write (*, "(' ** Density*sigma^3=',f10.6)") natms/volumen
         else
            Write (*, "(' ** Density=',f10.6,' A^-3')") natms/volumen
         end if
         ! Cluster information
         if (run_clusters) then
            if (iconf>1.and.maxcln>cl_thresh) then
               if (geometry) then
                  write (*, "(' ** Average cluster radius',f10.3,' average &
                     &internal cluster density ',f10.7)") avradio/iconf, averdens/iconf
                  write (*, "(' ** Average cluster gyration radius',f10.3)") avrg/iconf
                  write (*, "(' ** Internal cluster density  ',f10.7)")&
                     & sum(densav(:))/maxcln
               endif
            end if
            write (*, "(' ** No. of clusters for this configuration :',i5)") maxcln
            print *, " ··Time for graph construction", tgraph/iconf
            print *, " ··Thrust time ", tthrus/iconf
            print *, " ··Time for adjacency list construction =", tadj/iconf
            print *, " ··Time for BFS cluster search =", tbfs/iconf
         end if
         if (run_rdf) print *, " ··Time for rdf ", trdf/iconf
         if (run_sq) print *, " ··Time for S(Q) ", tsQ/iconf
         if (run_order) print *, " ··Time for order parameter ", tord/iconf
         if (run_dyn) print *, " ··Time for dynamics ", tdyn/iconf
         time_gput = tthrus + tadj + tbfs + tgraph + trdf + tsQ + tord + tdyn
         print *, " ··Time config in/out  ", tread/iconf
      End if
 
   end subroutine print_output

   subroutine printPotEngCl()
      implicit none
      !
      ! Printout Potential's Energy IntraCluster
      !
      integer :: i
      open (100, file='distUcltot.dat')
      write (100, "('#     Ucl(Kcal/mol)           D(U_cl))')")
      do i = 1, potnbins
         Write (100, '(2f16.7)') (i-1)*deltapotcl+epotcl_min, epothistomixcl(i)/real(nconf)
      end do
      close (100)
      ! dstribution of intracluster configurational energy per atom
      open (100, file='distUcl_N.dat')
      write (100, "('#     U_cl(Kcal/mol)/N_cl           D(U_cl/N_cl)')")
      do i = 1, potnbins
         Write (100, '(2f16.7)') (i-1)*deltapotperatomcl+epotperatomcl_min, epotperatomhistomixcl(i)/real(nconf)
      end do
      close (100)
   end subroutine printPotEngCl

   subroutine printSQ(Nmol)
      implicit none
      !
      ! Printout S(Q)'s
      !
      real(myprec) :: x1, x2, s11, s22, s12, scc, baver, b2aver
      integer, intent(IN) :: Nmol
      integer :: i, j
      logical :: bsc_one = .true.
      character(len=128) :: fname99
      baver = sum(ntype(1:nsp)*bsc(1:nsp))/Nmol
      b2aver = sum(ntype(1:nsp)*bsc(1:nsp)**2)/Nmol
      sqf(:) = sqf(:)/baver**2
      do i = 1, nsp
         if (abs(bsc(i)-1.0) > 1.0e-4) bsc_one=.false.
      enddo
      open (100, file='sq.dat')
      open (110, file='sqmix.dat')
      if (idir>0) then
         open (120, file='sqxy.dat')
         write (120, "('#     Q        ',9x,16('S_xy(Q,',f8.3,')',8x:))") zslice(1:nslice)
         do j=1, nsp
            write(fname99,'("sqpxy_",i1,"-",i1,".dat")') j, j
            open (130+j, file=trim(adjustl(fname99))//'.dat')
            write (130+j, "('#     Q        ',9x,16('S_xy(Q,',f8.3,')',8x:))") zslice(1:nslice)
         end do
      end if
      if (nsp == 2 .and. bsc_one) then
         x1 = (real(ntype(1))/real(Nmol))
         x2 = (real(ntype(2))/real(Nmol))
         write (100, "('#           Q        S_NN(Q)          S_cc(Q)       S_11(Q)        S_22(Q)           S_12(Q)         n(Q)')")
      else
         write (100, "('#           Q       S_NN(Q)          n(Q)')")
      end if
      write (110, "('#       Q',14x,15('S_',2i1,'(Q)',9x:))") ((j, j), j=1, nsp)
      do i = 1, nqmax
         if (i*dq <= qmin .or. dq > 0.2) then
            ! For small q, print all data points. For larger q, print every 3rd point to reduce file size and noise.
            if (twoDstruc_3D) then
               ! Compute 2D structure factor in xy plane for 3D systems with confinement
               write (120, '(15f16.7)') i*dq, sqfxy(i, 1:nslice)/(Nconf*real(nq(i)))
               do j=1, nsp
                  write (130+j, '(15f16.7)') i*dq, sqfpxy(i, j, 1:nslice)/(Nconf&
                  &*real(nq(i)))
               end do
            else
               ! Compute partial structure factors and concentration-concentration
               ! structure factor for binary mixtures with bsc=1
               if (nsp == 2 .and. bsc_one) then
                  s11 = x1*sqfp(i, 1)/(ntype(1)*Nconf*real(nq(i)))
                  s22 = x2*sqfp(i, 2)/(ntype(2)*Nconf*real(nq(i)))
                  s12 = 0.5*(sqf(i)/(Nmol*Nconf*real(nq(i))) - s11 - s22)
                  scc = x2**2*s11 + x1**2*s22 - 2*x1*x2*s12
                  write (100, '(6f15.7,i12)') i*dq, sqf(i)/(Nmol*Nconf*real(nq(i)))&
                  &, scc, s11, s22, s12, nq(i)
               else
                  write (100, '(2f15.7,i12)') i*dq, sqf(i)/(Nmol*Nconf*real(nq(i))), nq(i)
               end if
               write (110, '(15f16.7)') i*dq, (sqfp(i, j)/(ntype(j)*Nconf&
               &*real(nq(i))), j=1, nsp)
            endif
         end if
      end do
      if (dq <= 0.2) then
         do i = nint(qmin/dq) + 1, nqmax - 2, 3
            if (twoDstruc_3D) then
               ! Compute 2D structure factor in xy plane for 3D systems with confinement, using 5-point moving average to reduce noise
               write (120, '(15f16.7)') i*dq, (sum(sqfxy(i - 2:i + 2,j),dim=1)/(Nconf&
               &*real(5*nq(i))), j=1, nslice)
               do j=1, nsp
                  write (130+j, '(15f16.7)') i*dq, sum(sqfpxy(i - 2:i + 2, j, 1:nslice),dim=1)/(Nconf&
                  &*real(5*nq(i)))
               end do
            else
               ! Compute partial structure factors and concentration-concentration
               ! structure factor for binary mixtures with bsc=1,
               ! using 5-point moving average to reduce noise
               if (nsp == 2 .and. bsc_one) then
                  s11 = x1*sum(sqfp(i - 2:i + 2, 1)/(ntype(1)*Nconf*real(nq(i - 2:i + 2))))/5
                  s22 = x2*sum(sqfp(i - 2:i + 2, 2)/(ntype(2)*Nconf*real(nq(i - 2:i + 2))))/5
                  s12 = 0.5*(sum(sqf(i - 2:i + 2)/(Nmol*Nconf*real(nq(i - 2:i + 2))))/5 - s11 - s22)
                  scc = x2**2*s11 + x1**2*s22 - 2*x1*x2*s12
                  write (100, '(6f15.7,i12)') i*dq, sum(sqf(i - 2:i + 2)/(Nmol*Nconf*real(nq(i - 2:i + 2))))/5&
                  &, scc, s11, s22, s12, nq(i)
               else
                  write (100, '(2f15.7,i12)') i*dq, sum(sqf(i - 2:i + 2)/(Nmol*Nconf*real(nq(i - 2:i + 2))))/5, nq(i)
               end if
               write (110, '(15f16.7)') i*dq, (sum(sqfp(i - 2:i + 2, j)/(ntype(j)*Nconf&
               &*real(nq(i - 2:i + 2))))/5, j=1, nsp)
            endif
         end do
      end if
      close (100)
      close (110)
      if (idir>0) then
         close (120)
         do j=1, nsp
            close (130+j)
         end do
      end if
   end subroutine printSQ

   subroutine printrdf(rcl, lsmax)
      !
      !   Print out pdf's
      !
      implicit none
      real(myprec), intent(IN) :: rcl
      Real(myprec) :: gmix(nspmax, nspmax), deltav, ri, xfj, xfi
      integer, intent(in) :: lsmax
      integer :: i, j, l, k, count, iunit
      character(len=128) :: fname99
      if (twoDstruc_3D) then
         count=0
         do i=1, nsp
            do j=i, nsp
               write(fname99,'("gxy_",i1,"-",i1,".dat")') i,j
               open (130+count, file=trim(adjustl(fname99))//'.dat')
               write (130+count, "('#     r        ',3x,16('g_xy(r,',f8.3,')',8x:))") zslice(1:nslice)
               count = count + 1
            end do
         end do  
      else
         if (nsp<=6) then
            fname99 = 'gmixsim.dat'
            Open (99, file=fname99)
         else
            do k=1, nsp
               write(fname99,'("gmixsim",i1,".dat")') k
               open(887+k,file=fname99)
            enddo
         endif
      endif 
      if (nrandom>0) Open (199, file="s2n.dat")
      if (run_clusters.and.maxcln>cl_thresh) then
         if (geometry) then 
            write (99, "('#       r',16x,'g_cl(r)        g_cl-cl(r)    ',5x,16('g_',2i1,'(r)',8x:))")&
            & (((j, k), k=j, nsp), j=1, nsp)
         else
            write (99, "('#       r',9x,16('g_',2i1,'(r)',8x:))")&
            & (((j, k), k=j, nsp), j=1, nsp)
         endif
         if (nrandom>0)  write (199, "('#       r',16x,'s2n(cl)       s2n')")
      else
         if (.not.twoDstruc_3D) then
         if (nsp<=6) then
            write (99, "('#       r        ',9x,16('g_',2i1,'(r)',9x:))") (((j, k), k=j, nsp), j=1, nsp)
         else
            do k=1, nsp
               write (887+k, "('#       r        ',9x,16('g_',2i1,'(r)',9x:))") ((k, j), j=1, nsp)
            enddo
         endif
      endif
         if (nrandom>0)  write (199, "('#       r',16x,'s2n')")
      end if
      Do i = 1, lsmax - 2
         ri = i*deltar
         !
         ! Compute 3d ad 2d normalizations factors
         !
         if (ndim == 3) then
            if (twoDstruc_3D) then
               ! For 2D rdf din xy plane of 3D systems with confinement, use cylindrical shell volume for normalization
               deltaV = pi*((ri + deltar/2)**2 - (ri - deltar/2)**2)*zgrid
            else
               ! For 3D systems, use spherical shell volume for normalization
               deltaV = 4*pi*((ri + deltar/2)**3 - (ri - deltar/2)**3)/3.0
            endif
         else
            ! For 2D systems, use circular shell area for normalization
            deltaV = pi*((ri + deltar/2)**2 - (ri - deltar/2)**2)
         end if
         !
         Do j = 1, nsp
            Do l = j, nsp
               if (twoDstruc_3D) then
                  ! Compute 2D rdf in xy plane for 3D systems with confinement, using appropriate normalization for cylindrical shells and accounting for slice thickness
                  gmix_xy(j,l,:) = (j/l + 1)*volumen*(zgrid/sidel(3))*histomix_xy(i, j, l,:)/(deltaV*Nconf)
                  gmix_xy(l,j,:) = gmix_xy(j,l,:)
               else
                  gmix(j, l) = (j/l + 1)*volumen*histomix(i, j, l)/(deltaV*ntype(l)*ntype(j)*Nconf)
                  gmix(l,j) = gmix(j,l)
               endif 
            End Do
         End Do
         ! Print number fluctuations to check for hyperuniformity, if requested
         if (nrandom>0) then
            if (i<lsmax) then
               if (run_clusters) then
                  Write (199, '(18f16.5)') ri,  gclr2(i)/(Nconf*nrandom)-(gclr(i)/(Nconf*nrandom))**2,&
                  & g2i(i)/(Nconf*nrandom)-(gi(i)/(Nconf*nrandom))**2
               else
                  Write (199, '(18f16.5)') ri,  g2i(i)/(Nconf*nrandom)-(gi(i)/(Nconf*nrandom))**2
               endif
            endif
         endif
         ! Print rdf's, with different formats for cluster vs. non-cluster analysis and for 2D vs. 3D systems
         if (run_clusters.and.maxcln>cl_thresh) then
            if (geometry) then
               Write (99, '(28f16.5)') ri,&
               & gclustav(i)/(deltaV*Nconf), 2*gclcl(i)/(deltaV*Nconf),(gmix(j, j:nsp), j=1, nsp)
            else
               Write (99, '(28f16.5)') ri,&
               & (gmix(j, j:nsp), j=1, nsp)
            endif
         else
            if (twoDstruc_3D) then
               iunit = 130
               do k = 1, nsp
                  do j=k, nsp
                     write (iunit, '(15f16.7)') ri, (gmix_xy(k, j, 1:nslice))
                  iunit = iunit + 1
                  end do
               end do 
            else
               if (nsp <= 6) then
                  Write (99, '(16f16.5)') ri,  (gmix(j, j:nsp), j=1, nsp)
               else
                  do k=1, nsp
                     write(887+k,'(16f16.5)') ri,  gmix(k, 1:nsp)
                  end do
               endif
            endif 
         end if
      end do
      ! Close files for rdf output, with different handling for 2D vs. 3D systems and for cluster vs. non-cluster analysis
      if (twoDstruc_3D) then 
         count = 0
         do i = 1, nsp
            do j=i, nsp
               close(130+count)
               count = count + 1
            end do
         end do
      else
         if (nsp<=6) then
            close (99)
         else
            do k=1,nsp
               close(887+k)
            end do
         endif
      endif
      if (nrandom>0) close (199)
   end subroutine printrdf

   subroutine print_clusinfo(nqmin, Nmol)
      implicit none
      integer, intent(IN) :: nqmin, Nmol
      integer :: i, ndist
      real(myprec) :: avcldens, deltaV, ri, suma, norm
      Write (*, "(' ** Average total number of particles in clusters ', f10.2)") NTclus/nconf
      Write (*, "(' ** Average total number of clusters ', I5,' larger than ',i3)") &
      nint(sum(sizedist(:))/real(nconf)), minPts
      avcldens = sum(sizedist(:)/real(nconf))/volumen
      Write (*, "(' ** Average cluster density ', f15.9)") avcldens
      open (1001, file='clustdistr.dat')
      if (geometry) then
         open (125, file='dens.dat')
         open (126, file='radii.dat')
         open (999, file='rhoprof.dat')
         write (125, "('#      rho_cl        %clusters(rho_c)   ')")
      !
      !
         do i = 1, ndrho
            write (125, "(5f15.7)") i*drho, real(densclus(i))/(sum(densclus(:))*drho)
         end do
         write (126, "('#      r                  RG_cl(r)')")
         write (999, "('#      r                  rho_cl(r)')")
         ! Average cluster density profile
         if (ndim==3) then
           write (999, '(2f15.6)') 0.0, rhoclusav(0)/(4*pi*((drclus/2)**3)/3.0*Nconf)
         else
            write (999, '(2f15.6)') 0.0, rhoclusav(0)/(pi*((drclus/2)**2)*Nconf)
         end if
         do i = 1, ndrclus
            ri = i*drclus
            if (ndim == 3) then
               deltaV = 4*pi*((ri + drclus/2)**3 - (ri - drclus/2)**3)/3.0
            else
               deltaV = pi*((ri + drclus/2)**2 - (ri - drclus/2)**2)
            end if
            if (radii(i) .ne. 0) write (126, "(2f15.7)") ri, (real(radii(i))/sum(radii(:))/drclus)
            write (999, '(2f15.6)') ri, rhoclusav(i)/(deltaV*Nconf)
         end do
      endif
      if (run_sq.and.geometry) then
         open (102, file='sqcl.dat')
         write (102, "('#      Q            S_cl-cl(Q)')")
         do i = 1, nqmin
            write (102, '(2f15.7)') i*dq, sqfcl(i)/(Nconf*real(nq(i)))
         end do
      end if
      !
      ! Normalize cluster size distribution a s*n(s)/Ntotal
      !
      write (1001, "('#   N                %clus.               rho_cl ')")
      suma = 0
      ndist = nint(real(Nmol)/real(dcl))
      norm = (0.5*(sizedist(1) + sizedist(ndist)) + sum(sizedist(2:ndist - 1)))*dcl
      do i = 1, ndist
         suma = i*sizedist(i)*dcl + suma
         if (sizedist(i) > 0) write (1001, '(4f15.7)') (i - 0.5)*dcl, real(i*sizedist(i))/(real(Nconf*dcl*Nmol)), real(sizedist(i))/norm
      end do
      close (1001)
      close (999)
      close (102)
      close (125)
      close (126)
   end subroutine print_clusinfo

   subroutine print_order()
      implicit none
      integer :: i, j, onunit
      real(myprec) :: dV
      open (newunit=onunit, file='order.dat')
      if ( run_clusters.and.maxcln>cl_thresh) then
         if (geometry) then
            if (ndim==2) then
               write (onunit, "('#   order   psi_m   psi_m_clust   ')")
            else
               write (onunit, "('#   order   Q_l(order)   Q_l_clust(order)  ')")
            end if
            do i = 1, norder
               write (onunit, '(i3,4f12.5)') orderp(i), avorder(i)/real(nconf), &
                  & avcluster_order(i)/real(nconf)
            end do
         else
            if (ndim==2) then
               write (onunit, "('#   order   psi_m   ')")
            else
               write (onunit, "('#   order   Q_l(order)   ')")
            end if
            do i = 1, norder
               write (onunit, '(i3,4f12.5)') orderp(i), avorder(i)/real(nconf)
            end do
         endif
      else
         if (ndim==2) then
            write (onunit, "('#   order  Real(psi_m)  Im(psi_m)     |psi_m|  ')")
            do i = 1, norder
               write (onunit, '(5x,i3,3f12.5,f12.5)') orderp(i), avorder_cos(i)/real(nconf), &
               & avorder_sin(i)/real(nconf), avorder(i)/real(nconf)
            end do
         else
            write (onunit, "('# order   Q_l(order)   ')")
            do i = 1, norder
               write (onunit, '(2x,i3,f12.5)') orderp(i), avorder(i)/real(nconf)
            end do
         end if
      end if 
      close (100)
      if (print_orderp) then
         open (newunit=onunit, file='order_per_mol.dat')
         if (ndim==2) then
            write (onunit, "('# mol        x           y     ',5x,15('Real(psi_m) Im(psi_m)(',i0,')',1x:))") (orderp(i), i=1, norder)
            do i = 1, nmol
               write (onunit, '(i6,15f13.5)') i, r(1:ndim,i), (atomic_order_cos(i,j), atomic_order_sin(i,j), j=1, norder)
            end do
         else
            write (onunit, "('# mol        x          y           z   ',7x,15('Q_l(',i2,')'5x:))") (orderp(i), i=1, norder)
            do i = 1, nmol
               write (onunit, '(i6,15f12.5)') i, r(1:ndim,i), atomic_ql(i,1:norder)
            end do
         end if
 
         close (onunit)
         if (geometry) then
            if (run_clusters.and.maxcln>cl_thresh) then
               open (newunit=onunit, file='order_per_clust.dat')
               if (ndim==2) then
                  write (onunit, "('# Cluster    x           y    ',5x,15('Real(psi_m) Im(psi_m)(',i0,')',1x:))")(orderp(i),i=1,norder)
                  do i = 1, maxcln
                     write (onunit, '(i5,15f13.5)') i, cluster(i)%center(1:ndim), (cluster_order_cos(i,j), cluster_order_sin(i,j), j=1, norder)
                  end do
               else
                  write (onunit, "('# Cluster    x          y           z   ',7x,15('Q_l(',i2,')'5x:))") (orderp(i), i=1, norder)
                  do i = 1, maxcln
                     write (onunit, '(i5,15f12.5)') i, cluster(i)%center(1:ndim), (cluster_ql(i,j), j=1, norder)
                  end do
               end if
               close (onunit)
               open (newunit=onunit, file='ordprof_clust.dat')
               write (onunit, "('#      r     ',15('   rho_psi',i2,'(r)':))")(orderp(i), i=1, norder)
               write (onunit, '(f10.5, 15f12.5)') 0.0, rhoorderav(0,1:norder)/(Nconf)
               do i = 1, ndrclus
                  write (onunit, '(f10.5, 15f12.5)') i*drclus, rhoorderav(i,1:norder)/(Nconf)
               end do
               close(onunit)
               open (newunit=onunit, file='ordprof_clcum.dat')
               write (onunit, "('#      r     ',15('   psi_cum',i2,'(r)':))")(orderp(i), i=1, norder)
               write (onunit, '(f10.5, 15f12.5)') 0.0, ordercumav(0,1:norder)/(Nconf)
               do i = 1, ndrclus
                  write (onunit, '(f10.5, 15f12.5)') i*drclus, ordercumav(i,1:norder)/(Nconf)
               end do
               close(onunit)
            end if
         endif
      end if   
   end subroutine print_order
end module mod_log
