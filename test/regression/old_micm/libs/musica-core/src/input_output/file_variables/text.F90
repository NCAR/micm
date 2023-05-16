! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_file_variable_text module

!> The file_variable_text_t type and related functions
module musica_file_variable_text

  use musica_constants,                only : musica_dk, musica_ik
  use musica_file_dimension_range,     only : file_dimension_range_t
  use musica_file_variable,            only : file_variable_t

  implicit none
  private

  public :: file_variable_text_t

  !> A text file variable
  !!
  !! In addition to standard file_variable_t functions, a text file variable
  !! provides access to the column or row the variable is located in.
  !!
  type, extends(file_variable_t) :: file_variable_text_t
    private
    !> Row or column the variable is located in (starting at 1)
    integer(kind=musica_ik) :: id_ = -1
    !> Variable dimensions
    type(file_dimension_range_t), allocatable :: dimensions_(:)
  contains
    !> Returns the row or column the variable is located in (starting at 1)
    procedure :: id
    !> Gets the variable dimensions
    procedure :: get_dimensions
    !> Gets a sub-set of the variable data for a specified index range
    !!
    !! Data are returned after applying conversions set up during
    !! initialization.
    !!
    procedure :: get_data
    !> Outputs data to the file for a given time step
    procedure :: output
    !> Returns a flag indicating whether two file_variable_t objects refer to
    !! the same file variable
    procedure :: is_same_as
  end type file_variable_text_t

  !> Constructor
  interface file_variable_text_t
    module procedure :: constructor_name, constructor_id
  end interface file_variable_text_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates a file_variable_text_t object for an existing text file variable
  !! by name for input files, or a new file variable for output files.
  !!
  function constructor_name( file, variable_name, dimensions, config, found ) &
      result( new_obj )

    use musica_assert,                 only : die
    use musica_config,                 only : config_t
    use musica_file,                   only : file_t
    use musica_file_text,              only : file_text_t

    !> New text file variable
    class(file_variable_t), pointer :: new_obj
    !> Text file
    class(file_t), intent(inout) :: file
    !> Variable name
    character(len=*), intent(in) :: variable_name
    !> Variable dimensions
    !!
    !! Only required for adding a file variable
    !!
    type(file_dimension_range_t), intent(in), optional :: dimensions(:)
    !> Configuration describing how to match to MUSICA variables
    !!
    !! If omitted, standard matching is applied
    type(config_t), intent(inout), optional :: config
    !> Optional flag that indicates whether the variable was found in the file
    logical, intent(out), optional :: found

    integer(kind=musica_ik) :: variable_id

    call file%check_open( )
    if( present( found ) ) then
      if( file%is_file_variable( variable_name ) ) then
        found = .true.
      else
        found = .false.
        return
      end if
    end if
    select type( file )
    class is( file_text_t )
      if( file%is_input( ) .and. .not. file%is_output( ) ) then
        variable_id = file%get_variable_id( variable_name )
      else if( .not. file%is_input( ) .and. file%is_output( ) ) then
        variable_id = file%add_variable( variable_name )
      else
        call die( 113385982 )
      end if
    class default
      call die( 542896218 )
    end select
    new_obj => constructor_id( file, variable_id, config = config )

  end function constructor_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates a file_variable_text_t object for an existing text variable by id
  !!
  !! If no matching state variable is found in the MUSICA domain, a null
  !! pointer is returned.
  !!
  function constructor_id( file, variable_id, config ) result( new_obj )

    use musica_assert,                 only : die
    use musica_config,                 only : config_t
    use musica_file,                   only : file_t
    use musica_file_text,              only : file_text_t
    use musica_string,                 only : string_t

    !> New text file variable
    class(file_variable_t), pointer :: new_obj
    !> Text file
    class(file_t), intent(inout) :: file
    !> Variable row or column position (starting from 1)
    integer(kind=musica_ik), intent(in) :: variable_id
    !> Configuration describing how to match to MUSICA variables
    type(config_t), intent(inout) :: config

    character(len=*), parameter :: my_name = "text file variable constructor"
    type(string_t) :: var_name
    type(string_t) :: dimension_names(1)

    allocate( file_variable_text_t :: new_obj )
    select type( new_obj )
    class is( file_variable_text_t )
      select type( file )
      class is( file_text_t )
        var_name    = file%get_variable_name( variable_id )
        new_obj%id_ = variable_id
        new_obj%dimensions_ = file%get_variable_dimensions( variable_id )
      class default
        call die( 569766402 )
      end select
    class default
      call die( 511927843 )
    end select
    ! text files only have a single dimension, which defaults to 'time'
    call config%get( "dimension name", dimension_names(1), my_name,           &
                     default = "time" )

    call new_obj%private_constructor( config, file, var_name, dimension_names )

  end function constructor_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the row or column of the variable in the file (starting at 1)
  integer(kind=musica_ik) function id( this )

    !> Text file variable
    class(file_variable_text_t), intent(in) :: this

    id = this%id_

  end function id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gets the variable dimensions
  function get_dimensions( this ) result( dimensions )

    !> Variable dimensions
    type(file_dimension_range_t), allocatable :: dimensions(:)
    !> Text file variable
    class(file_variable_text_t), intent(in) :: this

    dimensions = this%dimensions_

  end function get_dimensions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gets a subset of the data from the file
  !!
  !! Applies any necessary conversions to the raw input data
  !!
  subroutine get_data( this, file, range, values )

    use musica_assert,                 only : die
    use musica_file,                   only : file_t
    use musica_file_text,              only : file_text_t

    !> Text file variable
    class(file_variable_text_t), intent(in) :: this
    !> Text file
    class(file_t), intent(inout) :: file
    !> Range of data to return
    !!
    !! Only one range in the array is permitted to have size > 1, so that
    !! results are always returned as a rank 1 array
    !!
    type(file_dimension_range_t), intent(in) :: range(:)
    !> Values to return
    real(kind=musica_dk), target, intent(out) :: values(:)

    select type( file )
    class is( file_text_t )
      call file%get_data( this%id_, range, values )
    class default
      call die( 980276126 )
    end select
    call this%convert_to_musica_values( values )

  end subroutine get_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Outputs data to the file for a given timestep
  !!
  !! The state_value will be converted according to the variable configuration
  !! prior to outputting it to the file.
  !!
  subroutine output( this, file, unlimited_dimension_value, indices,          &
      state_value )

    use musica_assert,                 only : die
    use musica_file,                   only : file_t
    use musica_file_text,              only : file_text_t

    !> Text file variable
    class(file_variable_text_t), intent(inout) :: this
    !> Text output file
    class(file_t), intent(inout) :: file
    !> Unlimited dimension value
    real(kind=musica_dk), intent(in) :: unlimited_dimension_value
    !> Indices of data point to output in all non-unlimited dimensions
    class(file_dimension_range_t), intent(in) :: indices(:)
    !> Value to output
    real(kind=musica_dk), intent(in) :: state_value

    real(kind=musica_dk) :: temp_val(1)

    select type( file )
      class is( file_text_t )
        temp_val(1) = state_value
        call this%convert_to_file_values( temp_val )
        call file%output( this%id_, unlimited_dimension_value, indices,       &
                          temp_val(1) )
      class default
        call die( 786618560 )
    end select

  end subroutine output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns a flag indicating whether two file_variable_t objects refer to
  !! the same file variable
  logical elemental function is_same_as( a, b )

    !> File variable a
    class(file_variable_text_t), intent(in) :: a
    !> File variable b
    class(file_variable_t), intent(in) :: b

    is_same_as = .false.
    select type( b )
    class is( file_variable_text_t )
      is_same_as = a%id_ .eq. b%id_
    end select

  end function is_same_as

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_file_variable_text
