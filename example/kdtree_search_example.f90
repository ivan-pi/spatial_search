program kdtree_search_example

use rbfx_kdtree_search
implicit none

real(wp), target :: x(100), y(100), points(2,100), p3d(3,100)
type(kdtree_search) :: s
integer :: idx1(10), idx2(10), idx3(10)

call random_number(x)
call random_number(y)

!
! Use x- and y-coordinate arrays
!
call s%init(x, y)
call s%find_nearest_idx(p=[0.5_wp,0.5_wp], idx=idx1)
print *, idx1

call s%destroy()

end program