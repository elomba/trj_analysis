!===============================================================================
! Program: centers
!===============================================================================
! Purpose:
!   Post-processes cluster center trajectories to identify persistent
!   clusters and track their time evolution. Analyzes cluster lifetimes
!   and generates continuous trajectory files for visualization.
!
! Key Functionality:
!   - Reads cluster center trajectories from LAMMPS-format files
!   - Identifies cluster correspondences across time steps
!   - Tracks cluster time evolution and lifetimes
!   - Generates continuous trajectories for long-lived clusters
!   - Computes lifetime histograms
!
! Input Files:
!   - Input: Configuration file with parameters
!   - centers.lammpstrj: Cluster center positions from main analysis
!
! Output Files:
!   - centers_atl.lammpstrj: Filtered trajectories (clusters at-least-time)
!   - centers_histogram.dat: Cluster lifetime distribution
!
! Algorithm:
!   1. Read all cluster configurations
!   2. Compute minimum inter-cluster distance
!   3. Match clusters across consecutive timesteps
!   4. Build continuous trajectories
!   5. Filter by minimum lifetime threshold
!
! Notes:
!   - Uses minimum distance criterion for cluster matching
!   - Distance threshold auto-adjusted with safety factor
!   - Handles variable number of clusters per timestep
!===============================================================================
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

   ! Set working directory and input/output filenames
   directory = '.'
   input_filename = directory//'/'//'Input'
   centers_filename = directory//'/'//'centers.lammpstrj'
   trajectories_filename = directory//'/'//'centers_atl.lammpstrj'
   histogram_filename = directory//'/'//'centers_histogram.dat'
   ! Read input parameters
   call read_input_file()
   call read_centers_file()
   print *, 'time_steps = ', time_steps
   print *, 'side = ', side
   print *, 'min_n_centers = ', min_n_centers
   print *, 'max_n_centers = ', max_n_centers
   ! Calculate minimum distance between clusters for matching criterion
   call cpu_centers_min_distance()
   print *, 'cal_min_distance = ', min_distance
   ! Apply safety factor to distance threshold
   distance_adj = 0.5
   min_distance = min_distance*distance_adj
   print *, 'distance_corrector = ', distance_adj
   print *, 'est_min_distance = ', min_distance
   ! Build continuous trajectories by matching clusters across timesteps
   call cpu_trajectories_analysis()
   print *, 'clusters_processed = ', max_n_centers_final
   ! Write filtered trajectories to file
   call write_trajectories_atl()
   print *, 'n_clusters_atl = ', n_clusters_atl
   write (*, '(A,f7.2)') '[% cluster_atl] (cluster_atl / clusters_proc) = ', &
      (n_clusters_atl/real(max_n_centers_final))*100
   ! Write lifetime histogram
   call write_histogram()

end program centers
