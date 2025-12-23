
!> Provides a ``spatial_search`` child class implemented using the kdtree2_ library
!>
!> .. _nanoflann: https://github.com/ivan-pi/kdtree2
!>
module rbfx_kdtree_search
use rbfx_spatial_search, only: spatial_search, wp
use kdtree2_module, only: kdtree2, kdtree2_create, kdtree2_destroy
use kdtree2_module, only: kdtree2_result, kdtree2_n_nearest
implicit none
private

public :: kdtree_search, wp


!> kdtree2-based spatial search tree
type, extends(spatial_search) :: kdtree_search
    private
    type(kdtree2) :: tr !> Reference to a kdtree2 tree
contains
    procedure :: init =>  kdtree_create_xy
    procedure, non_overridable :: find_nearest_idx => kdtree_nearest_idx
    procedure, non_overridable :: find_nearest_idx_dist => kdtree_nearest_idx_dist
    procedure, non_overridable :: find_nearest_coords => kdtree_nearest_coords
    procedure :: destroy => kdtree_destroy
end type

contains

    subroutine kdtree_create_xy(this,x,y)
        class(kdtree_search), intent(out) :: this
        real(wp), intent(in) :: x(:), y(:)

        type(kdtree_search) :: tree

        real(wp), allocatable :: xy(:,:)

        if (size(x) /= size(y)) error stop "error: kdtree_create_xy: size mismatch of x and y"

        allocate(xy(2,size(x)))
        xy(1,:) = x
        xy(2,:) = y

        ! kdtree2 always creates a copy of the input data

        this%dims = 2
        this%tr = kdtree2_create(xy,this%dims,&
            sort=.true.,&
            rearrange=.false.,&
            bucket_size=12,&
            verbose=.true.)

    end subroutine

    subroutine kdtree_nearest_idx(this,p,idx)
        class(kdtree_search), intent(in) :: this
        real(wp), intent(in) :: p(this%dims)
        integer, intent(out), contiguous :: idx(:)

        ! Hopefully this will fit on the stack
        type(kdtree2_result) :: results(size(idx))

        call kdtree2_n_nearest(this%tr,p,size(results),results)
        idx = results%idx

    end subroutine

    subroutine kdtree_nearest_idx_dist(this,p,idx,dist)
        class(kdtree_search), intent(in) :: this
        real(wp), intent(in) :: p(this%dims)
        integer, intent(out), contiguous :: idx(:)
        real(wp), intent(out), contiguous :: dist(:)

        ! Hopefully this will fit on the stack
        type(kdtree2_result) :: results(size(idx))


        call kdtree2_n_nearest(this%tr,p,size(results),results)
        idx = results%idx
        dist = results%dis

    end subroutine

    subroutine kdtree_nearest_coords(this,p,m,coords,dist,idx)
        class(kdtree_search), intent(in) :: this
        real(wp), intent(in) :: p(this%dims)
        integer, intent(in) :: m
        real(wp), intent(out) :: coords(m,this%dims)
        real(wp), intent(out), optional :: dist(m)
        integer, intent(out), optional :: idx(m)

        ! Hopefully results can fit on the stack; if not, make
        ! this an allocatable array
        type(kdtree2_result) :: results(m)
        integer :: d

        call kdtree2_n_nearest(this%tr,p,m,results)

        ! Gather coordinates and shift into local stencil
        do d = 1, this%dims
            coords(:,d) = this%tr%input_data(d,results%idx) - p(d)
        end do

        if (present(idx)) idx = results%idx
        if (present(dist)) dist = results%dis

    end subroutine

    subroutine kdtree_destroy(this)
        class(kdtree_search), intent(inout) :: this
        this%dims = -1
        call kdtree2_destroy(this%tr)
    end subroutine

end module