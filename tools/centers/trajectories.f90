!===============================================================================
! Module: mod_trajectories
!===============================================================================
! Purpose:
!   Constructs continuous cluster trajectories by matching clusters
!   across consecutive timesteps. Analyzes cluster lifetimes and
!   persistence to identify long-lived clusters.
!
! Key Functionality:
!   - Matches clusters between consecutive configurations
!   - Builds continuous trajectories for persistent clusters
!   - Tracks cluster lifetimes
!   - Handles cluster birth and death events
!   - Identifies long-lived clusters for output
!
! Main Subroutine:
!   cpu_trajectories_analysis()
!     - Initializes trajectories from first configuration
!     - Iterates through all timesteps
!     - Matches clusters using minimum distance criterion
!     - Creates new trajectories for unmatched clusters
!     - Terminates trajectories when clusters disappear
!
! Algorithm:
!   1. Initialize first timestep trajectories
!   2. For each subsequent timestep:
!      a. Try to match each active trajectory to a cluster
!      b. Use distance < min_distance criterion
!      c. Ensure one-to-one matching
!      d. Create new trajectories for unmatched clusters
!      e. Mark disappeared clusters with negative coordinates
!
! Output:
!   clusters_time_life(:,:) - Binary matrix of cluster existence
!   clusters_time_steps(:)  - Total lifetime of each cluster
!
! Notes:
!   - Uses periodic boundary conditions for matching
!   - Prevents double-matching of clusters
!   - Negative coordinates mark terminated trajectories
!===============================================================================
module mod_trajectories
   use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
   use mod_common
   implicit none
contains
   subroutine cpu_trajectories_analysis()
      integer :: t0, t1, i, t0_center, t1_center, s
      integer, allocatable :: t0_clusters_processed(:), t1_clusters_processed(:)
      real(sp) :: distance
      real(sp) :: distance_arr(n_dim)
      print *, "Proccesing analysis of trajectories..."
      ! Initialize trajectories from first configuration
      do i = 1, host_configurations(1)%n_centers
         host_trajectories(i, 1)%coordinate = host_configurations(1)%centers(:, i)
         host_trajectories(i, 1)%velocity = host_configurations(1)%vel(:, i)
         clusters_time_life(i, 1) = 1  ! Mark cluster as alive
      end do
      max_n_centers_final = host_configurations(1)%n_centers

      ! Loop over consecutive timesteps to match clusters
      do t0 = 1, time_steps - 1
         t1 = t0 + 1
         ! Arrays to track which clusters have been matched
         allocate (t0_clusters_processed(max_n_centers_final))
         allocate (t1_clusters_processed(host_configurations(t1)%n_centers))
         t0_clusters_processed = 0
         t1_clusters_processed = 0
         ! Try to match each trajectory with a cluster in the next timestep
         do t0_center = 1, max_n_centers_final
            ! Skip terminated trajectories (negative coordinates)
            if (host_trajectories(t0_center, t0)%coordinate(1) < 0.0_sp) then
               cycle
            end if
            ! Search for matching cluster in next timestep
            do t1_center = 1, host_configurations(t1)%n_centers
               ! Calculate distance with PBC
               do s = 1, n_dim
                  distance_arr(s) = host_configurations(t1)%centers(s, t1_center) &
                                    - host_trajectories(t0_center, t0)%coordinate(s)
                  if (distance_arr(s) > side_div_2(s)) distance_arr(s) = distance_arr(s) - side(s)
                  if (distance_arr(s) < -side_div_2(s)) distance_arr(s) = distance_arr(s) + side(s)
                  distance_arr(s) = distance_arr(s)*distance_arr(s)
               end do
               distance = sqrt(sum(distance_arr))
               ! If clusters are close enough, match them
               if (distance < min_distance) then
                  ! Error check: ensure one-to-one matching
                  if (t0_clusters_processed(t0_center) == 1 .or. t1_clusters_processed(t1_center) == 1) then
                     print *, 'ERROR: center processed again!'
                     stop
                  end if
                  ! Continue this trajectory with matched cluster
                  host_trajectories(t0_center, t1)%coordinate = host_configurations(t1)%centers(:, t1_center)
                  host_trajectories(t0_center, t1)%velocity = host_configurations(t1)%vel(:, t1_center)
                  t0_clusters_processed(t0_center) = 1
                  t1_clusters_processed(t1_center) = 1
                  clusters_time_life(t0_center, t1) = 1  ! Increment lifetime
               end if
            end do
            ! Mark trajectory as terminated if no match found
            if (t0_clusters_processed(t0_center) == 0) then
               host_trajectories(t0_center, t1)%coordinate = (/-1.0, -1.0, -1.0/)
            end if
         end do
         ! Handle new clusters appearing in this timestep
         do t1_center = 1, host_configurations(t1)%n_centers
            if (t1_clusters_processed(t1_center) == 0) then
               if (max_n_centers_final < max_3n_centers) then
                  ! Create new trajectory for this cluster
                  max_n_centers_final = max_n_centers_final + 1
                  host_trajectories(max_n_centers_final, t1)%coordinate = host_configurations(t1)%centers(:, t1_center)
                  host_trajectories(max_n_centers_final, t1)%velocity = host_configurations(t1)%vel(:, t1_center)
                  clusters_time_life(max_n_centers_final, t1) = 1
               else
                  print *, 'ERROR: centers limit reached, check max_3n_centers variable!', max_3n_centers
                  stop
               end if
            end if
         end do
         deallocate (t0_clusters_processed, t1_clusters_processed)
      end do
      allocate (clusters_time_steps(max_n_centers_final))
      do i = 1, max_n_centers_final
         clusters_time_steps(i) = sum(clusters_time_life(i, :))
      end do
   end subroutine cpu_trajectories_analysis

end module mod_trajectories
