! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_grid module

!> The grid_t type and related functions
module musica_grid

  use musica_constants,                only : musica_dk, musica_ik
  use musica_iterator,                 only : iterator_t
  use musica_string,                   only : string_t

  implicit none
  private

  public :: grid_t, grid_section_t, grid_iterator_t

  !> Sub-section of a grid
  type grid_section_t
    private
    real(kind=musica_dk) :: lower_bound_
    real(kind=musica_dk) :: upper_bound_
  contains
    procedure, private :: range_equals_range
    generic :: operator(==)  => range_equals_range
    procedure, private :: range_not_equals_range
    generic :: operator(/=)  => range_not_equals_range
    procedure :: overlap     => grid_section_overlap
    procedure :: lower_bound => grid_section_lower_bound
    procedure :: upper_bound => grid_section_upper_bound
    procedure :: range       => grid_section_range
  end type grid_section_t

  !> Iterator over a 1D grid
  type, extends(iterator_t) :: grid_iterator_t
    private
    integer(kind=musica_ik) :: current_section_ = 0
    integer(kind=musica_ik) :: last_section_    = 1
  contains
    procedure :: next           => grid_iterator_next
    procedure :: reset          => grid_iterator_reset
  end type grid_iterator_t

  !> A 1D grid
  type grid_t
    private
    type(grid_section_t), allocatable :: sections_(:)
    type(string_t)                    :: units_
  contains
    procedure :: number_of_sections
    procedure :: units
    procedure :: lower_bound
    procedure :: upper_bound
    procedure :: range => grid_range
    procedure :: section
    procedure :: property_index
    procedure :: iterator
    procedure :: overlap
    procedure, private :: grid_equals_grid
    generic :: operator(==) => grid_equals_grid
    procedure, private :: grid_not_equals_grid
    generic :: operator(/=) => grid_not_equals_grid
    !> Private constructor - should only be called by extending types
    procedure :: private_constructor_bounds
  end type grid_t

  interface grid_t
    module procedure :: constructor_bounds
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @name Functions of the grid_t type
  !!
  !! @{

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of grid_t objects from lower and upper section bounds
  function constructor_bounds( lower_bounds, upper_bounds, units )            &
      result( new_grid )

    use musica_assert,                 only : assert_msg, die_msg

    !> Grid
    type(grid_t)        :: new_grid
    !> Lower wavelength bounds for each section
    real(kind=musica_dk), intent(in) :: lower_bounds(:)
    !> Upper wavelength bounds for each section
    real(kind=musica_dk), intent(in) :: upper_bounds(:)
    !> Base grid units
    type(string_t),       intent(in) :: units

    call new_grid%private_constructor_bounds( lower_bounds, upper_bounds,     &
                                              units )

  end function constructor_bounds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of grid_t objects from lower and upper section bounds
  !!
  !! This should only be called by constructors of extending types
  subroutine private_constructor_bounds( this, lower_bounds, upper_bounds,    &
      units )

    use musica_assert,                 only : assert_msg, die_msg

    !> Grid
    class(grid_t),        intent(inout) :: this
    !> Lower wavelength bounds for each section
    real(kind=musica_dk), intent(in)    :: lower_bounds(:)
    !> Upper wavelength bounds for each section
    real(kind=musica_dk), intent(in)    :: upper_bounds(:)
    !> Base grid units
    type(string_t),       intent(in)    :: units

    call assert_msg( 584221364,                                               &
                     size( lower_bounds ) .eq. size( upper_bounds ),          &
                     "Bad grid specification" )
    allocate( this%sections_( size( lower_bounds ) ) )
    this%sections_(:)%lower_bound_ = lower_bounds(:)
    this%sections_(:)%upper_bound_ = upper_bounds(:)
    this%units_ = units

  end subroutine  private_constructor_bounds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of sections in the grid
  integer elemental function number_of_sections( this )

    !> Grid
    class(grid_t), intent(in) :: this

    number_of_sections = size( this%sections_ )

  end function number_of_sections

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the base unit for the grid
  type(string_t) function units( this )

    !> Grid
    class(grid_t), intent(in) :: this

    units = this%units_

  end function units

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the lower bound for a particular grid section
  real(kind=musica_dk) elemental function lower_bound( this, iterator )

    !> Grid
    class(grid_t),          intent(in) :: this
    !> Grid iterator pointing to section of interest
    class(grid_iterator_t), intent(in) :: iterator

    lower_bound = this%sections_( iterator%current_section_ )%lower_bound( )

  end function lower_bound

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the upper bound for a particular grid section
  real(kind=musica_dk) elemental function upper_bound( this, iterator )

    !> Grid
    class(grid_t),          intent(in) :: this
    !> Grid iterator pointing to section of interest
    class(grid_iterator_t), intent(in) :: iterator

    upper_bound = this%sections_( iterator%current_section_ )%upper_bound( )

  end function upper_bound

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the range covered by a particular grid section in the base units
  !! of the grid
  real(kind=musica_dk) elemental function grid_range( this, iterator )

    !> Grid
    class(grid_t),          intent(in) :: this
    !> Grid iterator pointing to section of interest
    class(grid_iterator_t), intent(in) :: iterator

    grid_range = this%sections_( iterator%current_section_ )%range( )

  end function grid_range

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns a grid section
  function section( this, iterator )

    !> Grid section
    type(grid_section_t)               :: section
    !> Grid
    class(grid_t),          intent(in) :: this
    !> Grid iterator pointing to section of interest
    class(grid_iterator_t), intent(in) :: iterator

    section = this%sections_( iterator%current_section_ )

  end function section

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the index on a gridded property array referenced by the iterator
  integer elemental function property_index( this, iterator )

    !> Grid
    class(grid_t),          intent(in) :: this
    !> Grid iterator pointing to section of interest
    class(grid_iterator_t), intent(in) :: iterator

    property_index = iterator%current_section_

  end function property_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns an iterator for the grid
  function iterator( this )

    !> Grid iterator
    class(grid_iterator_t), pointer    :: iterator
    !> Grid
    class(grid_t),          intent(in) :: this

    allocate( grid_iterator_t :: iterator )
    iterator%last_section_ = this%number_of_sections( )

  end function iterator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the total overlap with another grid in the grid units
  real(kind=musica_dk) elemental function overlap( a, b )

    class(grid_t), intent(in) :: a
    class(grid_t), intent(in) :: b

    integer :: i, j

    overlap = 0.0
    do i = 1, size( a%sections_ )
      do j = 1, size( b%sections_ )
        overlap = overlap + a%sections_( i )%overlap( b%sections_( j ) )
      end do
    end do

  end function overlap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compares grid_t objects for equality
  logical elemental function grid_equals_grid( a, b ) result( equals )

    class(grid_t), intent(in) :: a
    class(grid_t), intent(in) :: b

    integer :: i_section

    equals = allocated( a%sections_ ) .and. allocated( b%sections_ )          &
             .and. ( a%units_ .eq. b%units_ )
    if( .not. equals ) return
    equals = size( a%sections_ ) .eq. size( b%sections_ )
    if( .not. equals ) return
    do i_section = 1, size( a%sections_ )
      equals = equals .and.                                                   &
               a%sections_( i_section ) .eq. b%sections_( i_section )
    end do

  end function grid_equals_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compares grid_t objects for inequality
  logical elemental function grid_not_equals_grid( a, b ) result( not_equals )

    class(grid_t), intent(in) :: a
    class(grid_t), intent(in) :: b

    not_equals = .not. a .eq. b

  end function grid_not_equals_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @}

  !> @name Functions of the grid_section_t type
  !!
  !! @{

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compares grid_section_t objects for equality
  logical elemental function range_equals_range( a, b ) result( equals )

    class(grid_section_t), intent(in) :: a
    class(grid_section_t), intent(in) :: b

    equals = a%lower_bound_ .eq. b%lower_bound_ .and.                         &
             a%upper_bound_ .eq. b%upper_bound_

  end function range_equals_range

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compares grid_section_t objects for inequality
  logical elemental function range_not_equals_range( a, b )                   &
      result( not_equals )

    class(grid_section_t), intent(in) :: a
    class(grid_section_t), intent(in) :: b

    not_equals = .not. a .eq. b

  end function range_not_equals_range

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the overlap between two grid elements
  !!
  !! The grid elements are assumed to have the same base units
  real(kind=musica_dk) elemental function grid_section_overlap( a, b )        &
      result( overlap )

    class(grid_section_t), intent(in) :: a
    class(grid_section_t), intent(in) :: b

    overlap = 0.0
    if( a%lower_bound_ .le. b%upper_bound_ .and.                              &
        b%lower_bound_ .le. a%upper_bound_ ) then
      overlap = max( 0.0, min( a%upper_bound_, b%upper_bound_ ) -             &
                          max( a%lower_bound_, b%lower_bound_ ) )
    end if

  end function grid_section_overlap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the lower bound of the grid section in the base units of the grid
  real(kind=musica_dk) elemental function grid_section_lower_bound( this )

    class(grid_section_t), intent(in) :: this

    grid_section_lower_bound = this%lower_bound_

  end function grid_section_lower_bound

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the upper bound of the grid section in the base units of the grid
  real(kind=musica_dk) elemental function grid_section_upper_bound( this )

    class(grid_section_t), intent(in) :: this

    grid_section_upper_bound = this%upper_bound_

  end function grid_section_upper_bound

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the range covered by the section in the base units of the grid
  real(kind=musica_dk) elemental function grid_section_range( this )

    class(grid_section_t), intent(in) :: this

    grid_section_range = this%upper_bound_ - this%lower_bound_

  end function grid_section_range

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @}

  !> @name Functions of the grid_iterator_t type
  !!
  !! @{

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Advances the iterator
  !!
  !! Returns false if the end of the collection has been reached
  logical function grid_iterator_next( this )

    class(grid_iterator_t), intent(inout) :: this

    this%current_section_ = this%current_section_ + 1
    grid_iterator_next = this%current_section_ .le. this%last_section_

  end function grid_iterator_next

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Resets the iterator
  subroutine grid_iterator_reset( this, parent )

    class(grid_iterator_t), intent(inout)          :: this
    class(iterator_t),      intent(in),   optional :: parent

    this%current_section_ = 0

  end subroutine grid_iterator_reset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_grid
