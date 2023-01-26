! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!>\file
!> The musica_lookup_axis module

!> The lookup_axis_t type and related functions
module musica_lookup_axis

  use musica_constants,                only : musica_dk
  use musica_string,                   only : string_t

  implicit none
  private

  public :: lookup_axis_t

  !> A single axis for lookup tables
  type :: lookup_axis_t
    private
    type(string_t)                    :: name_
    type(string_t)                    :: dimension_name_
    real(kind=musica_dk), allocatable :: values_(:)
  contains
    procedure :: name => axis_name
    procedure :: dimension_name
    procedure :: size => axis_size
    procedure :: find
    procedure, private :: interpolate_scalar
    procedure, private :: interpolate_1D
    generic :: interpolate => interpolate_scalar, interpolate_1D
  end type lookup_axis_t

  interface lookup_axis_t
    procedure :: constructor, constructor_array
  end interface lookup_axis_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of lookup_axis_t objects
  function constructor( config ) result( new_axis )

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_io,                     only : io_t
    use musica_io_netcdf,              only : io_netcdf_t
    use musica_string,                 only : to_char

    type(lookup_axis_t) :: new_axis
    class(config_t), intent(inout) :: config

    character(len=*), parameter :: my_name = "lookup_axis_t constructor"
    class(io_t), pointer :: lookup_file
    type(string_t) :: file_path, variable_name
    type(string_t), allocatable :: dim_names(:)

    call config%get( "file path",      file_path,     my_name )
    call config%get( "variable name",  variable_name, my_name )
    new_axis%name_ = variable_name
    lookup_file => io_netcdf_t( file_path )
    dim_names = lookup_file%variable_dimensions( variable_name, my_name )
    call assert_msg( 646167484, size( dim_names ) .eq. 1,                     &
                     "Expected 1 dimension for variable '"//                  &
                     trim( variable_name%to_char( ) )//"' in file '"//        &
                     trim( file_path%to_char( ) )//"', got "//                &
                     trim( to_char( size( dim_names ) ) )//" dimensions." )
    new_axis%dimension_name_ = dim_names(1)
    call lookup_file%read( variable_name, new_axis%values_, my_name )
    deallocate( lookup_file )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructs a lookup_axis_t object from an array of axis values
  function constructor_array( axis_name, axis_values ) result( new_axis )

    type(lookup_axis_t)              :: new_axis
    class(string_t),      intent(in) :: axis_name
    real(kind=musica_dk), intent(in) :: axis_values(:)

    new_axis%name_           = axis_name
    new_axis%dimension_name_ = axis_name
    allocate( new_axis%values_( size( axis_values ) ) )
    new_axis%values_(:)      = axis_values(:)

  end function constructor_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the name of the axis
  type(string_t) function axis_name( this )

    class(lookup_axis_t), intent(in) :: this

    axis_name = this%name_

  end function axis_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the dimension name for the axis
  type(string_t) function dimension_name( this )

    class(lookup_axis_t), intent(in) :: this

    dimension_name = this%dimension_name_

  end function dimension_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of points along the axis
  integer elemental function axis_size( this )

    class(lookup_axis_t), intent(in) :: this

    axis_size = size( this%values_ )

  end function axis_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finds a value along the axis and returns the closest index lower than the
  !! value and the fractional residual relative to the distance to the next
  !! higher index.
  !!
  !! Values beyond the limits of the axis return index=1, residual=0 (lower
  !! limit) or index=max_index-1, residual=1 (upper limit)
  elemental subroutine find( this, value, index, residual )

    !> Lookup axis
    class(lookup_axis_t), intent(in)  :: this
    !> Value to locate in along the axis
    real(kind=musica_dk), intent(in)  :: value
    !> Index along the axis for found value
    !!
    !! The index corresponds to the closest value that is less than or equal
    !! to the given value.
    !!
    !! The value of index passed to this routine will be used as the first
    !! guess of the index.
    integer,            intent(inout) :: index
    !> The fractional remainder of the value relative to the distance to the
    !! next value along the axis.
    real(kind=musica_dk), intent(out) :: residual

    if( value .ge. this%values_( size( this%values_ ) ) ) then
      index = size( this%values_ ) - 1
      residual = 1.0_musica_dk
      return
    else if( value .le. this%values_( 1 ) ) then
      index = 1
      residual = 0.0_musica_dk
      return
    end if
    do while( index .gt. 1 )
      if( this%values_( index ) .le. value ) exit
      index = index - 1
    end do
    do while( index .lt. size( this%values_ ) )
      if( this%values_( index + 1 ) .gt. value ) exit
      index = index + 1
    end do
    residual = ( value - this%values_( index ) )                            &
               / ( this%values_( index + 1 ) - this%values_( index ) )

  end subroutine find

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Interpolates between two scalar values based on a given residual
  elemental subroutine interpolate_scalar( this, lower_value, upper_value,    &
      residual, results )

    class(lookup_axis_t), intent(in)  :: this
    real(kind=musica_dk), intent(in)  :: lower_value
    real(kind=musica_dk), intent(in)  :: upper_value
    real(kind=musica_dk), intent(in)  :: residual
    real(kind=musica_dk), intent(out) :: results

    results = lower_value + (upper_value - lower_value) * residual

  end subroutine interpolate_scalar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Interpolates between two 1D arrays based on a given residual
  pure subroutine interpolate_1D( this, lower_values, upper_values, residual, &
      results )

    class(lookup_axis_t), intent(in)  :: this
    real(kind=musica_dk), intent(in)  :: lower_values(:)
    real(kind=musica_dk), intent(in)  :: upper_values(:)
    real(kind=musica_dk), intent(in)  :: residual
    real(kind=musica_dk), intent(out) :: results(:)

    results(:) = lower_values(:) +                                            &
                 (upper_values(:) - lower_values(:)) * residual

  end subroutine interpolate_1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_lookup_axis
