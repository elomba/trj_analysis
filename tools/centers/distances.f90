!===============================================================================
! Module: mod_distances
!===============================================================================
! Purpose:
!   Computes minimum distances between cluster centers across all
!   configurations. Uses periodic boundary conditions to correctly
!   handle distances in periodic systems.
!
! Key Functionality:
!   - Calculates all pairwise cluster distances
!   - Applies minimum image convention (PBC)
!   - Identifies global minimum distance
!   - Used to determine cluster matching threshold
!
! Main Subroutine:
!   cpu_centers_min_distance()
!     - Loops over all timesteps
!     - Computes all pair distances within each timestep
!     - Tracks global minimum
!     - Applies PBC using minimum image convention
!
! Algorithm:
!   For each configuration:
!     For each pair of clusters:
!       1. Compute distance vector
!       2. Apply minimum image (PBC)
!       3. Calculate Euclidean distance
!       4. Update global minimum
!
! Notes:
!   - Uses side_div_2 for efficient PBC application
!   - Distance threshold for matching = min_distance * adjustment_factor
!   - Critical for accurate cluster correspondence
!===============================================================================
module mod_distances
   use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
   use mod_common
   implicit none
contains
   subroutine cpu_centers_min_distance()
      integer :: i, j, k, s
      real(sp) :: distance
      real(sp) :: distance_arr(n_dim)
      print *, 'Proccesing minimal distances...'
      ! Loop over all timesteps and cluster pairs to find minimum separation
      do i = 1, time_steps
         do j = 1, host_configurations(i)%n_centers - 1
            do k = j + 1, host_configurations(i)%n_centers
               ! Calculate distance with PBC
               do s = 1, n_dim
                  distance_arr(s) = host_configurations(i)%centers(s, j) - host_configurations(i)%centers(s, k)
                  ! Apply minimum image convention
                  if (distance_arr(s) > side_div_2(s)) distance_arr(s) = distance_arr(s) - side(s)
                  if (distance_arr(s) < -side_div_2(s)) distance_arr(s) = distance_arr(s) + side(s)
                  distance_arr(s) = distance_arr(s)*distance_arr(s)
               end do
               distance = sqrt(sum(distance_arr))
               ! Track global minimum distance
               if (distance < min_distance) min_distance = distance
            end do
         end do
      end do
   end subroutine cpu_centers_min_distance
end module mod_distances
