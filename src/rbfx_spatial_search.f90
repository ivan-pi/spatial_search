!> Provides abstract interface for spatial searches
!>
!> To implement a new search you must extend the ``spatial_search`` class.
!>
module rbfx_spatial_search
implicit none

integer, parameter :: wp = kind(1.0d0)

!> Abstract class for spatial searches
type, abstract :: spatial_search
    integer :: dims !> Dimensionality of the data
contains
    procedure(find_nearest_idx_sub), deferred :: find_nearest_idx
    procedure(find_nearest_idx_dist_sub), deferred :: find_nearest_idx_dist
    procedure(find_nearest_coords_sub), deferred :: find_nearest_coords
    procedure :: destroy
end type

abstract interface
    !> Find indexes of points nearest to p
    subroutine find_nearest_idx_sub(this,p,idx)
        import spatial_search, wp
        class(spatial_search), intent(in) :: this
        real(wp), intent(in) :: p(this%dims)
        integer, intent(out), contiguous :: idx(:)
    end subroutine
    !> Find indexes and distances of points nearest to p
    subroutine find_nearest_idx_dist_sub(this,p,idx,dist)
        import spatial_search, wp
        class(spatial_search), intent(in) :: this
        real(wp), intent(in) :: p(this%dims)
        integer, intent(out), contiguous :: idx(:)
        real(wp), intent(out), contiguous :: dist(:)
    end subroutine
    !> Find shifted coordinates of points nearest p
    !> The coordinates are shifted into the local coordinate system of p.
    subroutine find_nearest_coords_sub(this,p,m,coords,dist,idx)
        import spatial_search, wp
        class(spatial_search), intent(in) :: this
        real(wp), intent(in) :: p(this%dims)
        integer, intent(in) :: m
        real(wp), intent(out) :: coords(m,this%dims)
        real(wp), intent(out), optional :: dist(m)
        integer, intent(out), optional :: idx(m)
    end subroutine
end interface

contains

    !> Destructor template method. Child classes should over-ride this method when applicable.
    subroutine destroy(this)
        class(spatial_search), intent(inout) :: this
        ! just a template method for over-riding
    end subroutine

end module
