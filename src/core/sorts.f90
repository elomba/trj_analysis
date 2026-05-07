!===============================================================================
! Module: sorts
!===============================================================================
! Purpose:
!   Collection of sorting algorithms for integer and real arrays.
!   Provides various implementations with different performance characteristics
!   ranging from O(N²) to O(N log N) complexity.
!
! Available Algorithms:
!   QUADRATIC O(N²):
!     - bubble_sort             - Classic bubble sort
!     - selection_sort          - Selection sort
!     - inline_selection_sort   - Optimized selection variant
!     - insertion_sort          - Insertion sort
!     - binary_insertion_sort   - Binary search insertion
!     - pair_insertion_sort     - Bidirectional insertion
!     - double_insertion_sort   - Dual-ended insertion
!
!   IMPROVED O(N^(4/3)) (possibly):
!     - shell_sort              - Shell sort with gap sequence
!
!   LOGARITHMIC O(N log N):
!     - merge_sort              - Standard merge sort
!     - alt_merge_sort          - Alternative merge implementation
!     - heap_sort               - Heap sort
!     - heap_sort_nr            - Non-recursive heap sort
!     - smoothsort              - Leonardo heap variant
!
!   O(N log N log N):
!     - odd_even_merge_sort       - Batcher's odd-even merge
!     - odd_even_merge_sort_knuth - Knuth's variant
!
!   AVERAGE O(N log N):
!     - quicksort               - Recursive quicksort
!     - quicksort_nr            - Non-recursive with stack
!     - dpquicksort             - Dual-pivot quicksort
!
!   LINEAR O(N):
!     - radixsort               - Radix sort for integers
!     - radixsort_opt           - Optimized radix sort
!
! Notes:
!   - Provided as examples, not exhaustively optimized
!   - Non-recursive variants avoid stack overflow for large arrays
!   - Radix sorts limited to integer keys
!   - Performance varies with data characteristics
!
! History:
!   - Original collection: M.J. Rutter, 9/2019
!   - 14/10/19, 18/10/19: Fixed indexing errors in heap sorts
!   - 4/4/20: Added smoothsort and non-recursive quicksort
!   - 30/12/22: Added odd-even merge sorts
!   - 9/11/23: Added radix sorts
!
! Source:
!   Michel Rutter's sorting collection
!   https://www.mjr19.org.uk/IT/sorts/
!
! Usage:
!   call quicksort(array, n)  ! In-place sort of array(1:n)
!===============================================================================
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
