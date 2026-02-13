!===============================================================================
! Module: mod_common
!===============================================================================
! Purpose:
!   Defines global data structures and parameters for cluster center
!   trajectory analysis. Provides shared access to configuration data,
!   trajectories, and analysis parameters.
!
! Key Data Structures:
!   configuration - Holds cluster centers for one timestep:
!     - id: Timestep identifier
!     - n_centers: Number of clusters
!     - centers(:,:): Cluster positions
!     - vel(:,:): Cluster velocities
!
!   center - Individual cluster properties:
!     - coordinate(3): Position
!     - velocity(3): Velocity
!
! Main Variables:
!   host_configurations() - Array of all timestep configurations
!   host_trajectories()   - Matrix of continuous cluster trajectories
!   host_distances()      - Minimum distances per timestep
!   clusters_time_life()  - Matrix tracking cluster existence
!   clusters_time_steps() - Total lifetime of each cluster
!   histogram()           - Lifetime distribution
!   side()                - Simulation box dimensions
!   min_distance          - Minimum cluster-cluster distance
!
! Parameters:
!   n_dim           - Spatial dimensions (3)
!   histogram_dt    - Histogram bin width
!
! Notes:
!   - All arrays dynamically allocated based on trajectory length
!   - Supports periodic boundary conditions
!===============================================================================
module mod_common
   use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
   implicit none
   character(:), allocatable :: input_filename
   character(:), allocatable :: centers_filename
   character(:), allocatable :: trajectories_filename
   character(:), allocatable :: histogram_filename
   integer, parameter :: n_dim = 3
   integer :: time_steps, max_3n_centers, max_n_centers_final, n_clusters_atl
   integer :: min_n_centers = 1000000, max_n_centers = 0
   integer, allocatable :: clusters_time_life(:, :)
   integer, allocatable :: clusters_time_steps(:)
   integer, allocatable :: histogram(:)
   real(sp), parameter :: histogram_dt = 5
   real(sp) :: min_distance = 1000000_sp
   real(sp) :: side(n_dim), side_div_2(n_dim),high(n_dim),low(n_dim)
   real(sp), allocatable :: host_distances(:)
   type :: configuration
      integer :: id
      integer :: n_centers
      real(sp), allocatable :: centers(:, :)
      real(sp), allocatable :: vel(:, :)
   end type
   type(configuration), allocatable :: host_configurations(:)
   type :: center
      real(sp) :: coordinate(n_dim) = 0.0_sp
      real(sp) :: velocity(n_dim) = 0.0_sp
   end type
   type(center), allocatable :: host_trajectories(:, :)
   save
contains

end module mod_common
