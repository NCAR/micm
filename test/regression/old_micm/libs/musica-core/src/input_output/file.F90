! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_file module

!> The file_t type and related functions
module musica_file

  use musica_constants,                only : musica_ik, musica_dk
  use musica_string,                   only : string_t

  implicit none
  private

  public :: file_t, private_constructor, get_file_unit, free_file_unit

  !> Maximum number of IO units
  integer, parameter :: kMaxFileUnits = 200
  !> Minimum unit number
  integer, parameter :: kMinFileUnit = 11
  !> Currently used file units
  logical, save :: units_used(kMaxFileUnits) = .false.

  !> A Input/Output file
  !!
  !! Handles reading or writing to/from a file and provides information
  !! about the structure and variables in a file.
  !!
  !! \todo add example usage for file_t
  !!
  type, abstract :: file_t
    private
    !> Path to the file
    type(string_t) :: path_
    !> Indicates file is for input
    logical :: is_input_ = .false.
    !> Indicates file is for output
    logical :: is_output_ = .false.
  contains
    !> Returns the name of the file
    procedure :: name => file_name
    !> Returns whether the file is an input
    procedure :: is_input
    !> Returns whether the file is an output
    procedure :: is_output
    !> Returns the type of file as a string
    procedure(file_type), deferred :: type
    !> Returns the number of dimensions in the file
    procedure(number_of_dimensions), deferred :: number_of_dimensions
    !> Returns the number of variables in the file
    procedure(number_of_variables), deferred :: number_of_variables
    !> Returns a flag indicating whether a variable exists in the file
    procedure(is_file_variable), deferred :: is_file_variable
    !> Opens the file if it is not currently open
    procedure(check_open), deferred :: check_open
    !> Closes the file
    procedure(close), deferred :: close
    !> Prints the file properties
    procedure :: print => do_print
  end type file_t

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the type of file as a string
  type(string_t) function file_type( this )
    use musica_string,                 only : string_t
    import file_t
    !> Input/Output file
    class(file_t), intent(in) :: this
  end function file_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Opens the file if it is not open already
  subroutine check_open( this )
    import file_t
    !> Input/Output file
    class(file_t), intent(inout) :: this
  end subroutine check_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Closes the Input/Output file
  subroutine close( this )
    import file_t
    !> Input/Output file
    class(file_t), intent(inout) :: this
  end subroutine close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of dimensions in the file
  integer(kind=musica_ik) function number_of_dimensions( this )
    use musica_constants,              only : musica_ik
    import file_t
    !> Input/Output file
    class(file_t), intent(inout) :: this
  end function number_of_dimensions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of variables in the file
  integer(kind=musica_ik) function number_of_variables( this )
    use musica_constants,              only : musica_ik
    import file_t
    !> Input/Output file
    class(file_t), intent(inout) :: this
  end function number_of_variables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns a flag indicating whether a variable exists in the file
  logical function is_file_variable( this, variable_name )
    import file_t
    !> Input/Output file
    class(file_t), intent(inout) :: this
    !> Variable to look for
    character(len=*), intent(in) :: variable_name
  end function is_file_variable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the file name
  function file_name( this )

    !> File name
    type(string_t) :: file_name
    !> Input/Output file
    class(file_t), intent(in) :: this

    file_name = this%path_

  end function file_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns whether the file is an input
  logical function is_input( this )

    !> Input/Output file
    class(file_t), intent(in) :: this

    is_input = this%is_input_

  end function is_input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns whether the file is an output
  logical function is_output( this )

    !> Input/Output file
    class(file_t), intent(in) :: this

    is_output = this%is_output_

  end function is_output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Prints the file properties
  subroutine do_print( this )

    !> Input/Output file
    class(file_t), intent(in) :: this

    write(*,*) "file path: "//this%path_%to_char( )
    write(*,*) "is input:", this%is_input_
    write(*,*) "is output:", this%is_output_

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Private constructor for common data members
  !! (Should only be called by constructors of extending types)
  subroutine private_constructor( this, config, default_output_file_name )

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t

    !> Input/Output file
    class(file_t), intent(inout), pointer :: this
    !> File configuration
    type(config_t), intent(inout) :: config
    !> Default output file name
    character(len=*), intent(in), optional :: default_output_file_name

    character(len=*), parameter :: my_name = "file_t constructor"
    type(string_t) :: file_intent
    logical :: file_name_found

    ! get file path
    call config%get( "file name", this%path_, my_name,                        &
                     found = file_name_found )

    ! get file intent (input/output)
    call config%get( "intent", file_intent, my_name )
    file_intent = file_intent%to_lower( )
    if( file_intent .eq. "input" ) then
      this%is_input_ = .true.
    else if( file_intent .eq. "output" ) then
      this%is_output_ = .true.
    else if( file_intent .eq. "input/output" ) then
      this%is_input_ = .true.
      this%is_output_ = .true.
    else
      call die_msg( 165481344, "Invalid intent specified for file '"//        &
                    this%path_%to_char( )//"': "//file_intent%to_char( ) )
    end if

    if( file_name_found ) return
    if( present( default_output_file_name ) .and.                             &
        .not. this%is_input_ ) then
      this%path_ = default_output_file_name
      return
    endif
    call die_msg( 497176670, "Missing file name." )

  end subroutine private_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gets an unused file unit
  integer function get_file_unit( )

    use musica_assert,                 only : die_msg

    integer :: i
    logical :: found

    found = .false.
    do i = 1, kMaxFileUnits
      if( .not. units_used( i ) ) then
        found = .true.
        exit
      end if
    end do
    if( .not. found ) then
      call die_msg( 895680497, "Maximum number of open file units reached" )
    end if
    units_used( i ) = .true.
    get_file_unit = i + kMinFileUnit

  end function get_file_unit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Frees a file unit
  subroutine free_file_unit( unit )

    !> File unit to free
    integer, intent(in) :: unit

    units_used( unit ) = .false.

  end subroutine free_file_unit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_file
