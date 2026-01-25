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
   integer :: io_log_file=55
contains
   subroutine log_init()
      use mod_input
      open (unit=io_log_file, file=log_output_file)
   end subroutine log_init

   subroutine log_clear()
      close (io_log_file)
   end subroutine log_clear

   subroutine print_output(iconf)
      integer, parameter :: nther=10, nprint=10
      integer, intent(in) :: iconf
      logical, dimension(nther) :: mascara
      character*15, dimension(nther) :: title =(/"     T(K)","KE (Kcal/mol)","PE Kcal/mol",&
         " Pressure (bar)","      Pxx(bar)","    Pyy(bar)","     Pzz(bar)",&
         "      Pxy(bar)","    Pxz(bar)","     Pyz(bar)"/)
      real(myprec) :: thermo_q(nther), Tfcl
      integer :: i, j
      thermo_q(:) = 0
      mascara(1) = ex_vel
      mascara(2) = ex_vel
      mascara(3) = run_thermo
      mascara(4:10) = ex_stress
      if (ex_vel) then
         thermo_q(1) = temperature
         thermo_q(2) = kelvintokcal*ekin*(aunit/tunit)**2/Rgas
      endif
      if (run_thermo) thermo_q(3) = epot
      if (ex_stress) then
         thermo_q(4) = pressure
         thermo_q(5:10) = pxyz(1:6)
      endif
      if (sum(mascara)) then
         if (iconf == 1) then
            open(1000, file="thermo_run.dat")
            write(1000,"('#    Conf  ',16a15)")pack(title(1:nther),mascara)
         endif
         write(1000,"(i9,12f15.5)")nstep, pack(thermo_q(1:nther),mascara)
      endif

      If (Mod(Iconf - 1, nprint) .Eq. 0) Then
         call cpu_time(cpu1)
         if (tunits == 'lj') then
            if (minclsize>0) then
               Write (*, "(/' ** Working on MD step no. ',i10,' time* =',f12.3,&
               & ' cpu time per conf.=',f7.2,' s:'/&
               & ' ** Clusters >= ',i3,' particles being analyzed '/)") nstep, nstep*tstep, &
                (cpu1 - cpu0)/nprint, minclsize
            else
               Write (*, "(/' ** Working on MD step no. ',i10,' time* =',f12.3,&
               & ' cpu time per conf.=',f7.2,' s:'/&
               & /)") nstep, nstep*tstep, &
                (cpu1 - cpu0)/nprint
            endif
         else
            if (minclsize>0) then
               Write (*, "(/' ** Working on MD step no. ',i10,' time =',f10.5,&
               & ' ns, cpu time per conf.=',f7.2,' s:'/&
               & ' ** Clusters >= ',i3,' particles being analyzed '/)") nstep, nstep*tstep/1000.0, &
                (cpu1 - cpu0)/nprint, minclsize
            else
               Write (*, "(/' ** Working on MD step no. ',i10,' time =',f10.5,&
               & ' ns, cpu time per conf.=',f7.2,' s:'/&
               & /)") nstep, nstep*tstep/1000.0, &
               (cpu1 - cpu0)/nprint
            endif
         endif
         cpu0 = cpu1
         if (run_thermo) then
            if (tunits == 'lj') then
               write (*, "(' ** Potential energy/epsilon=',f15.4,' Per atom=',f15.4)") epot, epotperatom
               write (*, "(' ** Average potential energy/epsilon=',f15.4, ' Per atom=',f15.4)") epotav/Iconf, epotav/(Iconf*nmol)
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
            else
               write (*, "(' ** Kinetic energy=',f15.4,' Kcal/mol, average=',f15.4,'Kcal/mol')") &
                  kelvintokcal*ekin*(aunit/tunit)**2/Rgas, 0.00198717*ecaver*(aunit/tunit)**2/Rgas/Iconf
            
            endif 
            if (run_clusters) then
               if (tunits == 'lj') then
                  write (*, "(' ** Cluster kinetic energy/epsilon=',f15.4,' average=',f15.4)") &
                     ekincl, ekclaver/Iconf
                  write (*, "(' ** Internal cluster kinetic energy/epsilon=',f15.4,' average=',f15.4)") &
                     ekincls, ekinclsav/Iconf
               else
                   write (*, "(' ** Cluster kinetic energy/epsilon=',f15.4,' Kcal/mol average=',f15.4,' Kcal/mol')") &
                     kelvintokcal*ekincl*(aunit/tunit)**2/Rgas, kelvintokcal*ekclaver*(aunit/tunit)**2/Rgas/Iconf
                  write (*, "(' ** Internal cluster kinetic energy=',f15.4,'&
                  & Kcal/mol, average=',f15.4,'Kcal&
                  &/mol')") kelvintokcal*ekincls*(aunit/tunit)**2/Rgas,&
                  & kelvintokcal*ekinclsav*(aunit/tunit)**2/Rgas/Iconf
               endif
               if (Nu_clus>1) then
                  Tfcl = (Nu_clus-1)*ndim
               else
                  Tfcl = Nu_clus*ndim
               endif
               if (tunits == 'lj') then
                    Write (*, "(' ** Average cluster k_b*temperature/epsilon =',f10.4&
                  &)") 2*ekclaver/(Tfcl*Iconf)
               else
                  Write (*, "(' ** Average cluster temperature =',f10.4&
                  &,' K')") 2*ekclaver*(aunit/tunit)**2/(Tfcl*Rgas*Iconf)
               endif
            end if
            if (tunits == 'lj') then
               write (*, "(' ** Average k_b*temperature/epsilon =',f10.4)") 2*ecaver/(Tfact*Iconf)
            else
               write (*, "(' ** Average temperature =',f10.4&
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
            if (iconf>1) then
               write (*, "(' ** Average cluster radius',f8.3,' average &
              &cluster density ',f10.7)") avradio/iconf, averdens/iconf
            endif
            write (*, "(' ** Average cluster gyration radius',f8.3)") avrg/iconf
            write (*, "(' ** Internal cluster density  ',f10.7)")&
            & sum(densav(:))/Nu_clus
            write (*, "(' ** No. of clusters for this configuration :',i5)") Nu_clus
            print *, " ··Time for graph construction", tgraph/iconf
            print *, " ··Thrust time ", tthrus/iconf
            print *, " ··Time adjacency list construction =", tadj/iconf
            print *, " ··Time for BFS cluster search =", tbfs/iconf
         end if
         if (run_rdf) print *, " ··Time for rdf ", trdf/iconf
         if (run_sq) print *, " ··Time for S(Q) ", tsQ/iconf
         if (run_order) print *, " ··Time for order parameter ", tord/iconf
         if (run_dyn) print *, " ··Time for dynamics ", tdyn/iconf
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
      baver = sum(ntype(1:nsp)*bsc(1:nsp))/Nmol
      b2aver = sum(ntype(1:nsp)*bsc(1:nsp)**2)/Nmol
      sqf(:) = sqf(:)/baver**2
      do i = 1, nsp
         if (abs(bsc(i)-1.0) > 1.0e-4) bsc_one=.false.
      enddo
      open (100, file='sq.dat')
      open (110, file='sqmix.dat')
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
         end if
      end do
      if (dq <= 0.2) then
         do i = nint(qmin/dq) + 1, nqmax - 2, 3
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
         end do
      end if
      close (100)
      close (110)
   end subroutine printSQ

   subroutine printrdf(rcl, lsmax)
      !
      !   Print out pdf's
      !
      implicit none
      real(myprec), intent(IN) :: rcl
      Real(myprec) :: gmix(nspmax, nspmax), deltav, ri, xfj
      integer, intent(in) :: lsmax
      integer :: i, j, l, k
      if (nsp<=6) then
         fname99 = 'gmixsim.dat'
         Open (99, file=fname99)
      else
         do k=1, nsp
            write(fname99,'("gmixsim",i1,".dat")') k
            open(887+k,file=fname99)
         enddo
      endif
      if (nrandom>0) Open (199, file="s2n.dat")
      if (run_clusters) then
         write (99, "('#       r',16x,'g_cl(r)        g_cl-cl(r)        ',5x,16('g_',2i1,'(r)',8x:))")&
         & (((j, k), k=j, nsp), j=1, nsp)
         if (nrandom>0)  write (199, "('#       r',16x,'s2n(cl)       s2n')")
      else
         if (nsp<=6) then
            write (99, "('#       r        ',16x,16('g_',2i1,'(r)',9x:))") (((j, k), k=j, nsp), j=1, nsp)
         else
            do k=1, nsp
               write (887+k, "('#       r        ',9x,16('g_',2i1,'(r)',9x:))") ((k, j), j=1, nsp)
            enddo
         endif
         if (nrandom>0)  write (199, "('#       r',16x,'s2n')")
      end if

      Do i = 1, lsmax - 2
         ri = i*deltar
         !
         ! Compute 3d ad 2d normalizations factors
         !
         if (ndim == 3) then
            deltaV = 4*pi*((ri + deltar/2)**3 - (ri - deltar/2)**3)/3.0
         else
            deltaV = pi*((ri + deltar/2)**2 - (ri - deltar/2)**2)
         end if
         !
         Do j = 1, nsp

            Do l = j, nsp
               xfj = real(ntype(j), kind=8)/Real(natms, kind=8)
               gmix(j, l) = (j/l + 1)*volumen*histomix(i, j, l)/(deltaV*ntype(l)*ntype(j)*Nconf)
               gmix(l,j) = gmix(j,l)
               if (idir > 0) gmix(j, l) = densty*gmix(j, l)/rdenst

            End Do
         End Do
         if (i<lsmax) then
            if (nrandom>0) then
               if (run_clusters) then
                  Write (199, '(18f16.5)') i*deltar,  gclr2(i)/(Nconf*nrandom)-(gclr(i)/(Nconf*nrandom))**2,&
                  & g2i(i)/(Nconf*nrandom)-(gi(i)/(Nconf*nrandom))**2
               else
                  Write (199, '(18f16.5)') i*deltar,  g2i(i)/(Nconf*nrandom)-(gi(i)/(Nconf*nrandom))**2
               endif
            endif
         endif
         if (run_clusters) then
            Write (99, '(28f16.5)') i*deltar,&
            & gclustav(i)/(deltaV*Nconf), 2*gclcl(i)/(deltaV*Nconf),(gmix(j, j:nsp), j=1, nsp)
         else
            if (nsp <= 6) then
               Write (99, '(16f16.5)') i*deltar,  (gmix(j, j:nsp), j=1, nsp)
            else
               do k=1, nsp
                  write(887+k,'(16f16.5)') i*deltar,  gmix(k, 1:nsp)
               enddo
            endif
         end if
      End Do
      if (nsp<=6) then
         close (99)
      else
         do k=1,nsp
            close(887+k)
         end do
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
      nint(sum(sizedist(:))/real(nconf)), minclsize
      avcldens = sum(sizedist(:)/real(nconf))/volumen
      Write (*, "(' ** Average cluster density ', f15.9)") avcldens
      open (125, file='dens.dat')
      open (126, file='radii.dat')
      open (999, file='rhoprof.dat')
      open (1001, file='clustdistr.dat')
      write (125, "('#      rho_cl        %clusters(rho_c)   ')")
      !
      ! Note, the cluster density distribution is commputed in units reduced with and
      ! estimated particle diameters
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
      if (run_sq) then
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
      if ( run_clusters) then
         write (onunit, "('#   order   psi_m   psi_m_clust   ')")
         do i = 1, norder
            write (onunit, '(i3,4f12.5)') orderp(i), avorder(i)/real(nconf), &
            & avcluster_order(i)/real(nconf)
         end do
      else
         if (ndim==2) then
            write (onunit, "('#   order   Real(psi_m)   Im(psi_m)   |psi_m|   ')")
            do i = 1, norder
               write (onunit, '(i3,3f12.5,f12.5)') orderp(i), avorder_cos(i)/real(nconf), &
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
            write (onunit, "('# mol        x       y     ',15('Real(psi_m)     Im(psi_m)   ',i3,8x:))") (orderp(i), i=1, norder)
            do i = 1, nmol
               write (onunit, '(i6,15f12.5)') i, r(1:ndim,i), (atomic_order_cos(i,j), atomic_order_sin(i,j), j=1, norder)
            end do
         else
            write (onunit, "('# mol        x          y           z   ',7x,15('Q_l(',i2,')'5x:))") (orderp(i), i=1, norder)
            do i = 1, nmol
               write (onunit, '(i6,15f12.5)') i, r(1:ndim,i), atomic_ql(i,1:norder)
            end do
         end if
 
         close (onunit)
         if (run_clusters) then
            open (newunit=onunit, file='order_per_mol_clust.dat')
            if (ndim==2) then
               write (onunit, "('# Cluster    x       y     ',15('Real(psi_m)     Im(psi_m)   ',i3,8x:))") (orderp(i), i=1, norder)
            else
               write (onunit, "('# Cluster    x       y       z   ',15('Real(psi_m)     Im(psi_m)   ',i3,8x:))") (orderp(i), i=1, norder)
            end if
            do i = 1, nbigcl
               write (onunit, '(i5,15f12.5)') i, r(1:ndim,i), (cluster_order_cos(i,j), cluster_order_sin(i,j), j=1, norder)
            end do
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
      end if   
   end subroutine print_order
end module mod_log
