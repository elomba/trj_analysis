module mod_io
   use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
   use mod_common
   use mod_distances
   implicit none
   integer ::nsteps(100000)
contains
   subroutine read_input_file()
      integer :: io
      print *, 'Reading input file: ', input_filename
      open (newunit=io, file=input_filename, status='old', action='read')
      read (io, *) time_steps
      close (io)
      allocate (host_configurations(time_steps))
      allocate (host_distances(time_steps))
      host_distances = 1000000_sp
   end subroutine read_input_file

   subroutine read_centers_file()
      integer :: io, i, j, k, trash1, trash2, trash3
      real :: trashr
      print *, 'Reading centers file: ', centers_filename
      open (newunit=io, file=centers_filename, status='old', action='read')
      do i = 1, 5
         read (io, *)
      end do
      do i = 1, n_dim
         read (io, *) low(i), high(i)
      end do
      side = high-low
      side_div_2 = side/2.0
      rewind (unit=io)
      do i = 1, time_steps
         host_configurations(i)%id = i
         read (io, *)
         read(io,*) nsteps(i)
         read(io,*)
         read (io, *) host_configurations(i)%n_centers
         if (host_configurations(i)%n_centers < min_n_centers) &
            min_n_centers = host_configurations(i)%n_centers
         if (host_configurations(i)%n_centers > max_n_centers) &
            max_n_centers = host_configurations(i)%n_centers
         do k = 1, 5
            read (io, *)
         end do
         do j = 1, host_configurations(i)%n_centers
            read (io, *)
         end do
      end do
      rewind (unit=io)
      max_3n_centers = 500*max_n_centers
      allocate (host_trajectories(max_3n_centers, time_steps))
      do i = 1, time_steps
         do j = 1, max_3n_centers
            host_trajectories(j, i)%coordinate = (/-1.0, -1.0, -1.0/)
            host_trajectories(j, i)%velocity = (/-1.0, -1.0, -1.0/)
         end do
      end do
      allocate (clusters_time_life(max_3n_centers, time_steps))
      clusters_time_life = 0
      do i = 1, time_steps
         allocate (host_configurations(i)%centers(n_dim, host_configurations(i)%n_centers))
         allocate (host_configurations(i)%vel(n_dim, host_configurations(i)%n_centers))

         do k = 1, 9
            read (io, *)
         end do
         do j = 1, host_configurations(i)%n_centers
            read (io, *) trash1, trash2, host_configurations(i)%centers(1, j), &
               host_configurations(i)%centers(2, j), host_configurations(i)%centers(3, j),host_configurations(i)%vel(1:3, j)
         end do
      end do
      close (io)
   end subroutine read_centers_file

   subroutine write_trajectories_atl()
      integer :: io, i, j, s, id_tmp, ncid_out
      character :: fnamew*16
      print *, "Writing center's trajectories (atl): ", trajectories_filename
      n_clusters_atl = 0
      do i = 1, max_n_centers_final
         if (clusters_time_steps(i) == time_steps) then
            n_clusters_atl = n_clusters_atl + 1
         end if
      end do
    
    
      open (newunit=io, file=trajectories_filename, status='replace', action='write')
      do i = 1, time_steps
         write (io, "('ITEM: TIMESTEP'/I12/'ITEM: NUMBER OF ATOMS'/I12/'ITE&
              &M: BOX BOUNDS pp pp pp')") nsteps(i), n_clusters_atl
         do s = 1, n_dim
            write (io, "(2f15.7)") low(s), high(s)
         end do
         write (io, "('ITEM: ATOMS id type  x  y  z  vx  vy  vz')")
         id_tmp = 1
         do j = 1, max_n_centers_final
            if (clusters_time_steps(j) == time_steps) then
               write (io, '(2i10,6f15.7)') id_tmp, 1,  host_trajectories(j, i)%coordinate, host_trajectories(j, i)%velocity
               id_tmp = id_tmp + 1
            end if
         end do
      end do
      close (io)
   end subroutine write_trajectories_atl

   subroutine write_histogram()
      integer :: io, i, n_bins, bind
      print *, 'Writing histogram: ', histogram_filename
      print *, 'histogram_dt = ', histogram_dt
      n_bins = time_steps/histogram_dt
      print *, 'n_bins = ', n_bins
      allocate (histogram(n_bins))
      histogram = 0
      do i = 1, max_n_centers_final
         bind = ceiling(clusters_time_steps(i)/histogram_dt)
         histogram(bind) = histogram(bind) + 1
      end do
      open (newunit=io, file=histogram_filename, status='replace', action='write')
      do i = 1, n_bins
         write (io, '(f15.7,i6)') i*histogram_dt, histogram(i)
      end do
      close (io)
   end subroutine write_histogram

end module mod_io
