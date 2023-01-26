! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!>\file
!> The musica_lookup_2D_axis module

!> The lookup_2D_axis_t type and related functions
module musica_lookup_2D_axis

  use musica_lookup_axis,              only : lookup_axis_t

  implicit none
  private

  public :: lookup_2D_axis_t

  !> A 2D axis for lookup tables
  type :: lookup_2D_axis_t
    private
    type(lookup_axis_t) :: axes_(2)
  contains
    procedure :: names
    procedure :: dimension_names
    procedure :: sizes
    procedure :: find
    procedure, private :: interpolate_scalar
    procedure, private :: interpolate_1D
    generic :: interpolate => interpolate_scalar, interpolate_1D
  end type lookup_2D_axis_t

  interface lookup_2D_axis_t
    procedure :: constructor, constructor_axes
  end interface lookup_2D_axis_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of 2D lookup table axis
  function constructor( config ) result( new_axis )

    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    type(lookup_2D_axis_t)               :: new_axis
    class(config_t),       intent(inout) :: config

    character(len=*), parameter :: my_name = "lookup_2D_axis_t constructor"
    type(config_t) :: axis_config
    type(string_t) :: file_path, axis_name

    call config%get( "file path", file_path, my_name )
    call config%get( "axis 1",    axis_name, my_name )
    call axis_config%empty( )
    call axis_config%add( "file path",     file_path, my_name )
    call axis_config%add( "variable name", axis_name, my_name )
    new_axis%axes_(1) = lookup_axis_t( axis_config )
    call config%get( "axis 2",    axis_name, my_name )
    call axis_config%empty( )
    call axis_config%add( "file path",     file_path, my_name )
    call axis_config%add( "variable name", axis_name, my_name )
    new_axis%axes_(2) = lookup_axis_t( axis_config )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructs a 2D lookup table axis from 2 lookup_axis_t objects
  function constructor_axes( axis1, axis2 ) result( new_axis )

    type(lookup_2D_axis_t)          :: new_axis
    type(lookup_axis_t), intent(in) :: axis1
    type(lookup_axis_t), intent(in) :: axis2

    new_axis%axes_(1) = axis1
    new_axis%axes_(2) = axis2

  end function constructor_axes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the name of each axis
  function names( this )

    use musica_string,                 only : string_t

    type(string_t)                      :: names(2)
    class(lookup_2D_axis_t), intent(in) :: this

    names(1) = this%axes_(1)%name( )
    names(2) = this%axes_(2)%name( )

  end function names

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the name of each axis
  function dimension_names( this )

    use musica_string,                 only : string_t

    type(string_t)                      :: dimension_names(2)
    class(lookup_2D_axis_t), intent(in) :: this

    dimension_names(1) = this%axes_(1)%dimension_name( )
    dimension_names(2) = this%axes_(2)%dimension_name( )

  end function dimension_names

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of points along each axis
  function sizes( this )

    integer                             :: sizes(2)
    class(lookup_2D_axis_t), intent(in) :: this

    sizes(1) = this%axes_(1)%size( )
    sizes(2) = this%axes_(2)%size( )

  end function sizes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finds values along the axes and returns the closest index lower than the
  !! value and the fractional residual relative to the distance to the next
  !! higher index for each axis.
  !!
  !! Values beyond the limits of the axis return index=1, residual=0 (lower
  !! limit) or index=max_index-1, residual=1 (upper limit)
  pure subroutine find( this, values, indices, residuals )

    use musica_constants,              only : musica_dk

    !> 2D lookup
    class(lookup_2D_axis_t), intent(in) :: this
    !> Values to locate along each axis
    real(kind=musica_dk),    intent(in) :: values(2)
    !> Index along each axis for the found values
    !!
    !! The indices correspond to the closest value that is less than or equal
    !! to the given value for each axis.
    !!!
    !! The value of the index passed to this routine will be used as the first
    !!! guess of the index.
    integer,                 intent(inout) :: indices(2)
    !> The fractional remainder of the value relative to the distance to the
    !! next value along each axis.
    real(kind=musica_dk),    intent(out)   :: residuals(2)

    integer :: i_axis

    do i_axis = 1, 2
      call this%axes_( i_axis )%find( values(    i_axis ),                    &
                                      indices(   i_axis ),                    &
                                      residuals( i_axis ) )
    end do

  end subroutine find

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Interpolates between two scalar variables based on given residuals
  pure subroutine interpolate_scalar( this, table_data, residuals, results )

    use musica_constants,              only : musica_dk

    !> 2D lookup
    class(lookup_2D_axis_t), intent(in)  :: this
    !> Lookup table data (axis 1, axis 2)
    real(kind=musica_dk),    intent(in)  :: table_data(2,2)
    !> Residuals from axis lookup (axis 1, axis 2)
    real(kind=musica_dk),    intent(in)  :: residuals(2)
    !> Interpolated results
    real(kind=musica_dk),    intent(out) :: results

    real(kind=musica_dk) :: tu, tuc, tcuc, tcu

    tu  = residuals(1) * residuals(2)
    tuc = residuals(1) - tu
    tcuc = 1.0_musica_dk - tuc - residuals(2)
    tcu  = residuals(2) - tu
    results = tcuc * table_data( 1, 1 ) +                                     &
               tuc * table_data( 2, 1 ) +                                     &
                tu * table_data( 2, 2 ) +                                     &
               tcu * table_data( 1, 2 )

  end subroutine interpolate_scalar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Interpolates between two 1D arrays based on given residuals
  pure subroutine interpolate_1D( this, table_data, residuals, results )

    use musica_constants,              only : musica_dk

    !> 2D lookup
    class(lookup_2D_axis_t), intent(in)  :: this
    !> Lookup table data (output array, axis 1, axis 2)
    real(kind=musica_dk),    intent(in)  :: table_data(:,:,:)
    !> Residuals from axis lookup (axis 1, axis 2)
    real(kind=musica_dk),    intent(in)  :: residuals(2)
    !> Interpolated results
    real(kind=musica_dk),    intent(out) :: results(:)

    real(kind=musica_dk) :: tu, tuc, tcuc, tcu

    tu  = residuals(1) * residuals(2)
    tuc = residuals(1) - tu
    tcuc = 1.0_musica_dk - tuc - residuals(2)
    tcu  = residuals(2) - tu
    results(:) = tcuc * table_data( :, 1, 1 ) +                               &
                  tuc * table_data( :, 2, 1 ) +                               &
                   tu * table_data( :, 2, 2 ) +                               &
                  tcu * table_data( :, 1, 2 )

  end subroutine interpolate_1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_lookup_2D_axis
