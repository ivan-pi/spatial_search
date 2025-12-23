
!> Provides a ``spatial_search`` child class implemented using the kdtree2_ library
!>
!> .. _nanoflann: https://github.com/ivan-pi/kdtree2
!>
module rbfx_kdtree_search
use rbfx_spatial_search, only: spatial_search
use kdtree2_module, only: kdtree2, kdtree2_create, kdtree2_n_nearest
implicit none

!> kdtree2-based spatial search tree
type, extends(spatial_search) :: kdtree_search
    private
    type(kdtree2) :: tr !> Reference to a kdtree2 tree
contains
    procedure, non_overridable :: find_nearest_idx => kdtree_nearest_idx
    procedure, non_overridable :: find_nearest_coords => kdtree_nearest_coords
    procedure, non_overridable :: destroy => kdtree_destroy
end type

contains

    function kdtree_create_xy(x,y) result(tree)
        real(wp), intent(in) :: x(:), y(:)

        type(kdtree_search) :: tree

        real(wp), allocatable :: xy(:,:)

        if (size(x) /= size(y)) error stop "error: kdtree_create_xy: size mismatch of x and y"

        allocate(xy(2,size(x)))
        xy(1,:) = x
        xy(2,:) = y

        tree%dims = 2
        tree%tr = kdtree2_create(xy,tree%dims,&
            sort=.true.,&
            rearrange=.false.,&
            verbose=verbose)

    end function

    subroutine kdtree_nearest_idx(this,p,m,idx)
        class(kdtree_search), intent(in) :: this
        real(wp), intent(in) :: p(this%dims)
        integer, intent(in) :: m
        integer, intent(out) :: idx(m)

        ! Hopefully this will fit on the stack
        type(kdtree2_result) :: results(m)

        call kdtree2_n_nearest(this%tr,p,m,results)
        idx = results%idx

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
        call kdtree2_destroy(this%tr)
    end subroutine

end module