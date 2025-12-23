!> Provides a ``spatial_search`` child class implemented using the nanoflann_ library
!>
!> ## Notes
!> * Searching is limited to Euclidean space (|L2_Simple_Adaptor|_)
!> * Only double precision data and and default integer indices can be used
!> * Nanoflann wrappers are located in ``nf_wrapper.cpp``
!>
!> .. _nanoflann: https://jlblancoc.github.io/nanoflann/index.html
!> .. |L2_Simple_Adaptor| replace:: ``L2_Simple_Adaptor``
!> .. _L2_Simple_Adaptor: https://jlblancoc.github.io/nanoflann/structnanoflann_1_1L2__Simple__Adaptor.html
module rbfx_nf_search

use, intrinsic :: iso_c_binding, only: c_double, c_int, c_ptr, c_loc
use, intrinsic :: iso_fortran_env, only: error_unit

use rbfx_spatial_search, only: spatial_search, wp

implicit none
private

public :: nf_search, wp

enum, bind(c)
    enumerator :: nf_type_xy    = 1
    enumerator :: nf_type_xyz   = 2
    enumerator :: nf_type_array = 3
end enum

!> nanoflann-based spatial search tree
type, extends(spatial_search) :: nf_search
    private
    integer :: tree_type = 0
    type(c_ptr) :: ptr !> Reference to heap-allocated nanoflann handle

    real(c_double), pointer :: data(:,:) => null()

    real(c_double), pointer :: x(:) => null()
    real(c_double), pointer :: y(:) => null()
    real(c_double), pointer :: z(:) => null()

contains
    procedure, private :: nf_init_xy !> Initialize from x- and y-data.
    procedure, private :: nf_init_xyz !> Initialize from x-, y-, and z-data.
    procedure, private :: nf_init_array !> Initialize from an array.

    !> Generic initializer function
    generic :: init => nf_init_xy, nf_init_xyz, nf_init_array

    !> Find indices of nearest points
    procedure, non_overridable :: find_nearest_idx => nf_nearest_idx
    !> Find indices and distances of nearest points
    procedure, non_overridable :: find_nearest_idx_dist => nf_nearest_idx_dist

    !> Find shifted coordinates of nearest points
    procedure, non_overridable :: find_nearest_coords => nf_nearest_coords

    !> Destroy and free resources used by the nf_search object
    procedure :: destroy => nf_destroy
end type

! Minimal wrapper around nanoflann; see nf_wrappers.cpp
interface
    function c_nf_init(dims,n,data) bind(c)
        import c_int, c_double, c_ptr
        integer(c_int), value :: dims, n
        real(c_double), intent(in) :: data(dims,n)
        type(c_ptr) :: c_nf_init
    end function
    function c_nf_init_xy(n,x,y) bind(c)
        import c_int, c_double, c_ptr
        integer(c_int), value :: n
        real(c_double), intent(in) :: x(n), y(n)
        type(c_ptr) :: c_nf_init_xy
    end function
    function c_nf_init_xyz(n,x,y,z) bind(c)
        import c_int, c_double, c_ptr
        integer(c_int), value :: n
        real(c_double), intent(in) :: x(n), y(n), z(n)
        type(c_ptr) :: c_nf_init_xyz
    end function
    subroutine c_nf_free(ptr) bind(c)
        import c_ptr
        type(c_ptr), value :: ptr
    end subroutine
    subroutine c_nf_n_nearest(ptr,p,nn,idx,dist) bind(c)
        import c_ptr, c_double, c_int
        type(c_ptr), value :: ptr
        real(c_double), intent(in) :: p(*)
        integer(c_int), intent(in), value :: nn
        integer(c_int), intent(out) :: idx(nn)
        real(c_double), intent(out) :: dist(nn) ! squared distances
    end subroutine
end interface

contains

    !> Initialize search from contiguous array
    !> @param this The search object to be initialized
    !> @param points a contiguous d-by-n array of points. The size along the first dimension sets the dimensionality of the data.
    !>
    !> ## Notes
    !> The `points` actual argument should have either `target` or `pointer` attribute in the caller context.
    subroutine nf_init_array(this,points)
        class(nf_search), intent(out) :: this
        real(wp), intent(in), pointer, contiguous :: points(:,:)


        associate(d => size(points,1), n => size(points,2))

            if (d > n) then
                write(error_unit,*) "warning: nf_init_array: dimensionality of data is larger than number of points"
            end if

            this%dims = d
            this%tree_type = nf_type_array
            ! Use zero-based indexing
            this%data(1:d, 0:n-1) => points

            this%ptr = c_nf_init(d,n,this%data)
        end associate

    end subroutine

    !> Initialize a 2-d search from x- and y-coordinate arrays.
    !> @param this The search object to be initialized.
    !> @param x Array of x-coordinates.
    !> @param y Array of y-coordinates.
    !>
    !> ## Notes
    !> Arrays x and y should have the same size. They should have either
    !> `target` or `pointer` attribute in the caller context.
    !>
    subroutine nf_init_xy(this,x,y)
        class(nf_search), intent(out) :: this
        real(wp), intent(in), pointer, dimension(:), contiguous :: x, y

        if (size(x) /= size(y)) error stop "error: nf_init_xy: x and y have different sizes"

        associate(n => size(x))
            this%dims = 2
            this%tree_type = nf_type_xy
            this%x(0:n-1) => x
            this%y(0:n-1) => y
            this%ptr = c_nf_init_xy(n, x, y)
        end associate

    end subroutine

    !> Initialize a 3-d search from x- and y-coordinate arrays.
    !> @param this The search object to be initialized.
    !> @param x Array of x-coordinates.
    !> @param y Array of y-coordinates.
    !> @param y Array of z-coordinates.
    !>
    !> ## Notes
    !> Arrays x, y and z should be of the same size. They should have either
    !> `target` or `pointer` attribute in the caller context.
    subroutine nf_init_xyz(this,x,y,z)
        class(nf_search), intent(out) :: this
        real(wp), intent(in), pointer, contiguous :: x(:), y(:), z(:)

        if (size(x) /= size(y)) error stop "error: nf_init_xy: x and y have different sizes"
        if (size(x) /= size(z)) error stop "error: nf_init_xy: x and z have different sizes"

        this%tree_type = nf_type_xyz
        associate(n => size(x))
            this%dims = 3
            this%x(0:n-1) => x
            this%y(0:n-1) => y
            this%z(0:n-1) => z
            this%ptr = c_nf_init_xyz(n, x, y, z)
        end associate

    end subroutine

    !> Find indexes of nearest points
    subroutine nf_nearest_idx(this,p,idx)
        class(nf_search), intent(in) :: this
        real(wp), intent(in) :: p(this%dims)
        integer, intent(out), contiguous :: idx(:)

        real(wp) :: dist(size(idx))

        call c_nf_n_nearest(this%ptr,p,size(idx),idx,dist)
        idx = idx + 1

    end subroutine

    !> Find indexes and distances of nearest points
    subroutine nf_nearest_idx_dist(this,p,idx,dist)
        class(nf_search), intent(in) :: this
        real(wp), intent(in) :: p(this%dims)
        integer, intent(out), contiguous :: idx(:)
        real(wp), intent(out), contiguous :: dist(:)

        call c_nf_n_nearest(this%ptr,p,size(idx),idx,dist)
        idx = idx + 1

    end subroutine

    !> Find shifted coordinates of m-nearest points
    subroutine nf_nearest_coords(this,p,m,coords,dist,idx)
        class(nf_search), intent(in) :: this
        real(wp), intent(in) :: p(this%dims)
        integer, intent(in) :: m
        real(wp), intent(out) :: coords(m,this%dims)
        real(wp), intent(out), optional :: dist(m)
        integer, intent(out), optional :: idx(m)

        real(wp) :: dist_(m)
        integer :: idx_(m), d

        call c_nf_n_nearest(this%ptr,p,m,idx_,dist_)

        select case(this%tree_type)
        case(nf_type_xy)
            coords(:,1) = this%x(idx_) - p(1)
            coords(:,2) = this%y(idx_) - p(2)
        case(nf_type_xyz)
            coords(:,1) = this%x(idx_) - p(1)
            coords(:,2) = this%y(idx_) - p(2)
            coords(:,3) = this%z(idx_) - p(3)
        case(nf_type_array)
            do d = 1, this%dims
                coords(:,d) = this%data(d,idx_) - p(d)
            end do
        case default
            error stop "error: nf_nearest_coords: unknown tree type"
        end select

        if (present(idx)) idx = idx_ + 1
        if (present(dist)) dist = sqrt(dist_)

    end subroutine

    !> Destroy and free resources used by the nf_search object
    subroutine nf_destroy(this)
        class(nf_search), intent(inout) :: this

        this%dims = -1
        this%tree_type = 0

        call c_nf_free(this%ptr)

        nullify(this%data)
        nullify(this%x)
        nullify(this%y)
        nullify(this%z)

    end subroutine

end module
