! These various sort routines are provided as
! examples, rather than being fully optimised
! for the best possible performance, or even
! exhaustively tested. MJR 9/2019, 1/2023
!
! So far the collection is:
!
!  bubble_sort               (order N**2)
!  selection_sort            (order N**2)
!  inline_selection_sort     (order N**2)
!  insertion_sort            (order N**2)
!  binary_insertion_sort     (order N**2)
!  pair_insertion_sort       (order N**2)
!  double_insertion_sort     (order N**2)
!  shell_sort                (order N**(4/3)) (possibly)
!  odd_even_merge_sort       (order N log N log N)
!  odd_even_merge_sort_knuth (order N log N log N)
!  merge_sort               (order N log N)
!  alt_merge_sort           (order N log N)
!  heap_sort                (order N log N)
!  heap_sort_nr             (order N log N)
!  smoothsort               (order N log N)
!  quicksort                (order N log N) (usually)
!  quicksort_nr             (order N log N) (usually)
!  dpquicksort              (order N log N) (usually)
!  radixsort                (order N)
!  radixsort_opt            (order N)

! 14/10/19 and 18/10/19 indexing errors in heap sorts corrected

! 4/4/20 smooth sort and non-recursive quicksort added, and some
!        minor code tidying

! 30/12/22 add odd_even_merge_sort and odd_even_merge_sort_knuth

! 9/11/23 add radixsorts
! From  Michel Rutter's https://www.mjr19.org.uk/IT/sorts/

module sorts
  use mod_precision
  implicit none
  integer, parameter :: prec=myprec



contains


  ! This version maintains its own stack, to avoid needing to call
  ! itself recursively. By always pushing the larger "half" to the
  ! stack, and moving directly to calculate the smaller "half",
  ! it can guarantee that the stack needs no more than log_2(N)
  ! entries
  subroutine quicksort_nri(array)
    integer, intent(inout)::array(:)
    integer :: temp,pivot
    integer :: i,j,left,right,low,high
    ! If your compiler lacks storage_size(), replace
    ! storage_size(i) by 64
    integer :: stack(2,storage_size(i)),stack_ptr

    low=1
    high=size(array)
    stack_ptr=1

    do

       if (high-low.lt.50) then ! use insertion sort on small arrays
          do i=low+1,high
             temp=array(i)
             do j=i-1,low,-1
                if (array(j).le.temp) exit
                array(j+1)=array(j)
             enddo
             array(j+1)=temp
          enddo
          ! now pop from stack
          if (stack_ptr.eq.1) return
          stack_ptr=stack_ptr-1
          low=stack(1,stack_ptr)
          high=stack(2,stack_ptr)
          cycle
       endif

       ! find median of three pivot
       ! and place sentinels at first and last elements
       temp=array((low+high)/2)
       array((low+high)/2)=array(low+1)
       if (temp.gt.array(high)) then
          array(low+1)=array(high)
          array(high)=temp
       else
          array(low+1)=temp
       endif
       if (array(low).gt.array(high)) then
          temp=array(low)
          array(low)=array(high)
          array(high)=temp
       endif
       if (array(low).gt.array(low+1)) then
          temp=array(low)
          array(low)=array(low+1)
          array(low+1)=temp
       endif
       pivot=array(low+1)

       left=low+2
       right=high-1
       do
          do while(array(left).lt.pivot)
             left=left+1
          enddo
          do while(array(right).gt.pivot)
             right=right-1
          enddo
          if (left.ge.right) exit
          temp=array(left)
          array(left)=array(right)
          array(right)=temp
          left=left+1
          right=right-1
       enddo
       if (left.eq.right) left=left+1
       !          call quicksort(array(1:left-1))
       !          call quicksort(array(left:))
       if (left.lt.(low+high)/2) then
          stack(1,stack_ptr)=left
          stack(2,stack_ptr)=high
          stack_ptr=stack_ptr+1
          high=left-1
       else
          stack(1,stack_ptr)=low
          stack(2,stack_ptr)=left-1
          stack_ptr=stack_ptr+1
          low=left
       endif

    enddo
  end subroutine quicksort_nri
  
end module sorts
