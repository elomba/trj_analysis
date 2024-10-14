module mod_log
   use mod_precision
   use mod_common
   use mod_input
   use mod_nc
   use mod_nc_conf
   use cudafor
   use mod_thermo, only : engclus, engclcl
   implicit none
   integer :: io_log_file
contains
   subroutine log_init()
      use mod_input
      open (newunit=io_log_file, file=log_output_file)
   end subroutine log_init

   subroutine log_clear()
      close (io_log_file)
   end subroutine log_clear

   subroutine print_output(iconf)
      integer, intent(in) :: iconf
      real(myprec) :: pideal_i, pideal_av
      If (Mod(Iconf - 1, 5) .Eq. 0) Then
         call cpu_time(cpu1)
         Write (*, "(/' ** Working on MD step no. ',i8,' time =',f10.5,' ns, cpu time=',f15.2&
              &/)") nstep, nstep*tstep/1000.0, cpu1 - cpu0
         cpu0 = cpu1
         if (run_thermo) then
            write (*, "(' ** Potential energy=',f15.4,' Kcal/mol, Per atom=',f15.4,'Kcal/mol')") epot, epotperatom
            write (*, "(' ** Average potential energy=',f15.4,' Kcal/mol, Per atom=',f15.4,'Kcal/mol')") epotav/Iconf, epotav/(Iconf*nmol)
            if (run_clusters) then
               write (*, "(' ** Intracluster av. potential energy=',f15.4,' Kcal/mol')") engclus
               write (*, "(' ** Intercluster potential energy=',f15.4,' Kcal/mol')") engclcl
            endif 
         end if
         if (ex_stress) then
            write (*, "(' ** Pressure =',f15.4,' bar, Avergage pressure=',f15.4' bar')") pressure, pressav/Iconf !kcal_a3_to_bar*pressure, kcal_a3_to_bar*pressav/Iconf
         endif
         ! Kinetic energy and temperature
         If (ex_vel) then
            write (*, "(' ** Kinetic energy=',f15.4,' Kcal/mol, average=',f15.4,'Kcal/mol')") &
            kelvintokcal*ekin*(aunit/tunit)**2/Rgas, 0.00198717*ecaver*(aunit/tunit)**2/Rgas/Iconf
            if (rcl > 0) then
               write (*, "(' ** Cluster kinetic energy=',f15.4,' Kcal/mol, average=',f15.4,'Kcal/mol')") &
               kelvintokcal*ekincl*(aunit/tunit)**2/Rgas, 0.00198717*ekclaver*(aunit/tunit)**2/Rgas/Iconf
               write (*, "(' ** Internal cluster kinetic energy=',f15.4,'&
                 & Kcal/mol, average=',f15.4,'Kcal&
                 &/mol')") kelvintokcal*ekincls*(aunit/tunit)**2/Rgas,&
                 & kelvintokcal*ekinclsav*(aunit/tunit)**2/Rgas/Iconf
               Tfact = nint(sum(sizedist(:))/real(Iconf))*ndim
               Write (*, "(' ** Average cluster temperature =',f10.4&
                   &,' K')") 2*ekclaver*(aunit/tunit)**2/(Tfact*Rgas*Iconf)
            end if
            Write (*, "(' ** Temperature=',f10.4&
                     &,' K average=',f10.4,'K density=',f10.6' 1&
                     &/A^3')") temperature, taver, densty
         else
            Write (*, "(' Density*=',f10.6&
                 &)") natms/volumen
         end if
         ! Cluster information 
         if (rcl > 0) then
            write (*, "(' ** Average cluster radius',f8.3,' average &
              &cluster density ',f10.7)") avradio/iconf, averdens/iconf
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

   subroutine printPotEngClCl()
      implicit none
      !
      ! Printout Potential's Energy InterCluster
      !
      integer :: i
      epotclclmean = sum(epotclcldata) / nconf
      epotclcldata = epotclcldata - epotclclmean
      epotclcldata = epotclcldata * epotclcldata
      epotclclstdv = sum(epotclcldata) / nconf
      epotclclstdv = sqrt(epotclclstdv)
      write (*, "(' ** InterCluster AVG Potential energy=',f15.4,' Kcal/mol, STD_DEV=',f15.4)") epotclclmean, epotclclstdv

   end subroutine printPotEngClCl

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
      open (101, file='sqf.dat')
      open (110, file='sqmix.dat')
      if (nsp == 2 .and. bsc_one) then
        x1 = (real(ntype(1))/real(Nmol))
        x2 = (real(ntype(2))/real(Nmol))
        write (100, "('#           Q        S_NN(Q)          S_cc(Q)       S_11(Q)        S_22(Q)           S_12(Q)         n(Q)')")
        write (101, "('#           Q        S_NN(Q)          S_cc(Q)       S_11(Q)        S_22(Q)           S_12(Q)         n(Q)')")
      else
         write (100, "('#           Q       S_NN(Q)          n(Q)')")
         write (101, "('#           Q       S_NN(Q)          n(Q)')")
      end if
      write (110, "('#       Q',14x,6('S_',2i1,'(Q)',9x:))") ((j, j), j=1, nsp)
      do i = 1, nqmax
         if (i*dq <= qmin .or. dq > 0.2) then
         write (101, '(6f15.7,i12)') i*dq, sqf(i)/(Nmol*Nconf*real(nq(i)))
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

         write (110, '(7f15.7)') i*dq, (sqfp(i, j)/(ntype(j)*Nconf&
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
            write (110, '(7f15.7)') i*dq, (sum(sqfp(i - 2:i + 2, j)/(ntype(j)*Nconf&
                                                                    &*real(nq(i - 2:i + 2))))/5, j=1, nsp)
         end do
      end if
      do i=nint(qmin/dq)+1, nqmax
        write (101, '(2f15.7,i12)') i*dq,sqf(i)/(Nmol*Nconf*real(nq(i)))
      end do
      close (101)
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
      fname99 = 'gmixsim.dat'
      Open (99, file=fname99)
      if (nrandom>0) Open (199, file="s2n.dat")
      if (rcl > 0) then
         write (99, "('#       r',16x,'g_cl(r)        g_cl-cl(r)        ',5x,16('g_',2i1,'(r)',8x:))")&
         & (((j, k), k=j, nsp), j=1, nsp)
         if (nrandom>0)  write (199, "('#       r',16x,'s2n(cl)       s2n')")
      else
         write (99, "('#       r        ',16x,16('g_',2i1,'(r)',9x:))") (((j, k), k=j, nsp), j=1, nsp)
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
               if (idir > 0) gmix(j, l) = densty*gmix(j, l)/rdenst

            End Do
         End Do
         if (i<lsmax) then
            if (nrandom>0) then
               if (rcl>0) then
                  Write (199, '(18f16.5)') i*deltar,  gclr2(i)/(Nconf*nrandom)-(gclr(i)/(Nconf*nrandom))**2,&
                  & g2i(i)/(Nconf*nrandom)-(gi(i)/(Nconf*nrandom))**2 
               else
                  Write (199, '(18f16.5)') i*deltar,  g2i(i)/(Nconf*nrandom)-(gi(i)/(Nconf*nrandom))**2 
               endif
            endif
         endif
         if (rcl > 0) then
            Write (99, '(18f16.5)') i*deltar,&
                 & gclustav(i)/(deltaV*Nconf), 2*gclcl(i)/(deltaV*Nconf),(gmix(j, j:nsp), j=1, nsp)
         else
            Write (99, '(18f16.5)') i*deltar,  (gmix(j, j:nsp), j=1, nsp)
         end if
      End Do
      close (99)
      if (nrandom>0) close (199)
   end subroutine printrdf

   subroutine print_clusinfo(nqmin, Nmol)
      implicit none
      integer, intent(IN) :: nqmin, Nmol
      integer :: i, ndist
      real(myprec) :: avcldens, deltaV, ri, suma, norm
      Write (*, "(' ** Average total number of particles in clusters ', f10.2)") NTclus/nconf
      Write (*, "(' ** Average total number of clusters ', I5)") nint(sum(sizedist(:))/real(nconf))
      avcldens = sum(sizedist(:)/real(nconf))/volumen
      Write (*, "(' ** Average cluster density ', f15.9)") avcldens
      open (125, file='dens.dat')
      open (126, file='radii.dat')
      open (999, file='rhoprof.dat')
      open (1001, file='clustdistr.dat')
      write (125, "('#      rho_cl*        rho_cl        N(rho_cl)')")
      !
      ! Note, the cluster density distribution is commputed in units reduced with and
      ! estimated particle diameters
      !
      do i = 1, ndrho
         write (125, "(5f15.7)") i*drho, i*drho/sigma**3, (real(densclus(i))/drho/Nconf)
      end do
      write (126, "('#      r                  RG_cl(r)')")
      write (999, "('#      r                  rho_cl(r)')")
      write (999, '(2f15.6)') 0.0, rhoclusav(0)/(4*pi*((deltar&
      &/2)**3)/3.0*Nconf)
      do i = 1, ndr
         ri = i*deltar
         if (ndim == 3) then
            deltaV = 4*pi*((ri + deltar/2)**3 - (ri - deltar/2)**3)/3.0
         else
            deltaV = pi*((ri + deltar/2)**2 - (ri - deltar/2)**2)
         end if
         if (radii(i) .ne. 0) write (126, "(2f15.7)") ri, (real(radii(i))/sum(radii(:))/deltar)
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
end module mod_log
