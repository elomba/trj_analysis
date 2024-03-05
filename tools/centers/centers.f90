program centers
   use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
   use mod_common
   use mod_io
   use mod_distances
   use mod_trajectories
   implicit none
   integer :: i, istat
   real :: distance_adj
   character(:), allocatable :: directory

   directory = '.'
   input_filename = directory//'/'//'Input'
   centers_filename = directory//'/'//'centers.lammpstrj'
   trajectories_filename = directory//'/'//'centers_atl.lammpstrj'
   histogram_filename = directory//'/'//'centers_histogram.dat'
   call read_input_file()
   call read_centers_file()
   print *, 'time_steps = ', time_steps
   print *, 'side = ', side
   print *, 'min_n_centers = ', min_n_centers
   print *, 'max_n_centers = ', max_n_centers
   call cpu_centers_min_distance()
   print *, 'cal_min_distance = ', min_distance
   distance_adj = 0.5
   min_distance = min_distance*distance_adj
   print *, 'distance_corrector = ', distance_adj
   print *, 'est_min_distance = ', min_distance
   call cpu_trajectories_analysis()
   print *, 'clusters_processed = ', max_n_centers_final
   call write_trajectories_atl()
   print *, 'n_clusters_atl = ', n_clusters_atl
   write (*, '(A,f7.2)') '[% cluster_atl] (cluster_atl / clusters_proc) = ', &
      (n_clusters_atl/real(max_n_centers_final))*100
   call write_histogram()

end program centers
