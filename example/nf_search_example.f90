program nf_search_example

use rbfx_nf_search
implicit none

real(wp), target :: x(100), y(100), points(2,100), p3d(3,100)
type(nf_search) :: s
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

!
! Use single 2-by-n array of coordinates
!

points(1,:) = x
points(2,:) = y
call s%init(points)
call s%find_nearest_idx(p=[0.5_wp, 0.5_wp], idx=idx2)
print *, idx2

print *, "all indexes equal = ", all(idx1 == idx2)
call s%destroy()


!
! 3-d coordinate example
!
call random_number(p3d)
call s%init(p3d)
call s%find_nearest_idx(p=[0.5_wp, 0.5_wp, 0.5_wp], idx=idx3)
print *, idx3
call s%destroy

end program