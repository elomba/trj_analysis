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
      do i = 1, time_steps
         do j = 1, host_configurations(i)%n_centers - 1
            do k = j + 1, host_configurations(i)%n_centers
               do s = 1, n_dim
                  distance_arr(s) = host_configurations(i)%centers(s, j) - host_configurations(i)%centers(s, k)
                  if (distance_arr(s) > side_div_2(s)) distance_arr(s) = distance_arr(s) - side(s)
                  if (distance_arr(s) < -side_div_2(s)) distance_arr(s) = distance_arr(s) + side(s)
                  distance_arr(s) = distance_arr(s)*distance_arr(s)
               end do
               distance = sqrt(sum(distance_arr))
               if (distance < min_distance) min_distance = distance
            end do
         end do
      end do
   end subroutine cpu_centers_min_distance
end module mod_distances
