! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_interpolator module

!> The interpolator_t type and related functions
module musica_interpolator

  use musica_constants,                only : musica_dk
  use musica_grid,                     only : grid_t, grid_iterator_t

  implicit none
  private

  public :: interpolator_t, interpolator_element_t, interpolation_strategy_i

  !> Base mapping unit for interpolation
  !!
  !! Iterators are used to identify array elements in a gridded property
  !! array. The weight is used to add a contribution to the "to" array
  !! based on one element of the "from" array:
  !! \f[
  !!   to( i ) += from( j ) * weight
  !! \f]
  type :: interpolator_element_t
    private
    class(grid_iterator_t), allocatable :: from_
    class(grid_iterator_t), allocatable :: to_
    real(kind=musica_dk) :: weight_
  end type interpolator_element_t

  interface interpolator_element_t
    procedure :: interpolator_element_constructor
  end interface

  !> Interpolator between gridded data sets
  !!
  !! An example of how to use the interpolator object, which uses the linear
  !! 1D interpolation strategy follows:
  !!
  !! \snippet test/interpolator_strategies/linear_1D.F90 Linear interpolator example
  !!
  !! Output:
  !! \code{bash}
  !!    5.0000000000000000        5.0000000000000000        10.000000000000000        10.000000000000000        15.000000000000000        15.000000000000000
  !! \endcode
  type :: interpolator_t
    private
    class(grid_t), allocatable :: from_grid_
    class(grid_t), allocatable :: to_grid_
    type(interpolator_element_t), allocatable :: map_(:)
  contains
    procedure :: interpolate
  end type interpolator_t

  interface interpolator_t
    module procedure :: constructor
  end interface

  !> Interface for interpolator strategies
  !!
  !! Specific strategies should use the properties of the two grids to
  !! generate the set of \c interface_element_t objects that define the
  !! mapping used by the interpolator.
  !!
  !! Reference for the Strategy design pattern:
  !! http://fortranwiki.org/fortran/show/Strategy+Pattern
  abstract interface
    function interpolation_strategy_i( from_grid, to_grid ) result( map )
      use musica_grid,                 only : grid_t
      import interpolator_element_t
      class(interpolator_element_t), allocatable :: map(:)
      class(grid_t), intent(in) :: from_grid
      class(grid_t), intent(in) :: to_grid
    end function interpolation_strategy_i
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @name Functions of the interpolator_t type
  !!
  !! @{

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructs an interpolator_t for a given pair of grids
  function constructor( strategy, from_grid, to_grid )                        &
      result( new_interpolator )

    type(interpolator_t)                :: new_interpolator
    procedure(interpolation_strategy_i) :: strategy
    class(grid_t), intent(in)           :: from_grid
    class(grid_t), intent(in)           :: to_grid

    new_interpolator%from_grid_ = from_grid
    new_interpolator%to_grid_   = to_grid
    new_interpolator%map_       = strategy( from_grid, to_grid )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Interpolate between grids
  subroutine interpolate( this, from, to )

    class(interpolator_t), intent(in)    :: this
    real(kind=musica_dk),  intent(in)    :: from(:)
    real(kind=musica_dk),  intent(inout) :: to(:)

    integer :: i_map

    to(:) = 0.0
    do i_map = 1, size( this%map_ )
      to( this%to_grid_%property_index( this%map_( i_map )%to_ ) ) =          &
          to( this%to_grid_%property_index( this%map_( i_map )%to_ ) ) +      &
          from( this%from_grid_%property_index( this%map_( i_map )%from_ ) ) *&
          this%map_( i_map )%weight_
    end do

  end subroutine interpolate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @}

  !> @name Functions of the interpolator_element_t type
  !!
  !! @{

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of interpolator_element_t objects
  elemental function interpolator_element_constructor( from, to, weight )     &
      result( new_element )

    type(interpolator_element_t)       :: new_element
    class(grid_iterator_t), intent(in) :: from
    class(grid_iterator_t), intent(in) :: to
    real(kind=musica_dk),   intent(in) :: weight

    new_element%from_   = from
    new_element%to_     = to
    new_element%weight_ = weight

  end function interpolator_element_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_interpolator
