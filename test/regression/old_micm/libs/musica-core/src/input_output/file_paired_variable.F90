! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_file_paired_variable module

!> The file_paired_variable_t type and related functions
module musica_file_paired_variable

  use musica_constants,                only : musica_dk, musica_ik
  use musica_domain_state_accessor,    only : domain_state_accessor_t
  use musica_domain_state_mutator,     only : domain_state_mutator_t
  use musica_file_variable,            only : file_variable_t
  use musica_string,                   only : string_t

  implicit none
  private

  public :: file_paired_variable_t, file_paired_variable_ptr

  !> Max length of staged data array
  integer(kind=musica_ik), parameter :: kMaxStagedData = 100

  !> Data for a single paired MUSICA <-> file variable
  type :: file_paired_variable_t
    private
    !> MUSICA variable name
    type(string_t) :: musica_name_
    !> Mutator for the variable
    class(domain_state_mutator_t), pointer :: mutator_ => null( )
    !> Accessor for the variable
    class(domain_state_accessor_t), pointer :: accessor_ => null( )
    !> Variable information
    class(file_variable_t), pointer :: variable_ => null( )
    !> Index of first staged data
    integer(kind=musica_ik) :: first_staged_index_ = 1
    !> Number of staged data
    integer(kind=musica_ik) :: number_staged_ = 0
    !> Staged data
    real(kind=musica_dk) :: staged_data_(kMaxStagedData) = -huge(1.0_musica_dk)
  contains
    !> Returns a flag indicating whether the pair includes a given file
    !! variable
    procedure :: includes_variable
    !> Gets a value from the file
    procedure :: get_file_value
    !> Sets a value in the file
    procedure :: set_file_value
    !> Gets a MUSICA domain state value
    procedure :: get_musica_value
    !> Sets a MUSICA domain state value
    procedure :: set_musica_value
    !> MUSICA variable name
    procedure :: musica_name
    !> Prints the properties of the paired variable
    procedure :: print => do_print
    !> Updates the staged data
    procedure, private :: update_staged_data
    !> Finalizes a paired variable object
    final :: finalize
  end type file_paired_variable_t

  !> Constructor
  interface file_paired_variable_t
    module procedure :: constructor
  end interface file_paired_variable_t

  !> Unique pointer to file_paired_variable_t objects
  type :: file_paired_variable_ptr
    class(file_paired_variable_t), pointer :: val_ => null( )
  contains
    !> Finalizes the pointer
    final :: file_paired_variable_ptr_finalize
  end type file_paired_variable_ptr

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates a matched MUSICA <-> File variable pair
  !!
  !! If the file variable cannot be matched to a MUSICA domain variable, a
  !! null pointer is returned.
  !!
  function constructor( file, domain, variable, create_accessor )             &
      result( new_obj )

    use musica_assert,                 only : die
    use musica_data_type,              only : kDouble
    use musica_domain,                 only : domain_t
    use musica_domain_target_cells,    only : domain_target_cells_t
    use musica_file,                   only : file_t
    use musica_file_variable,          only : file_variable_t
    use musica_property,               only : property_t

    !> New MUSICA<->File variable match
    type(file_paired_variable_t), pointer :: new_obj
    !> File to update to or from
    class(file_t), intent(inout) :: file
    !> Model domain
    class(domain_t), intent(inout) :: domain
    !> File variable
    class(file_variable_t), intent(in) :: variable
    !> Flag indicating whether to force the inclusion of an accessor
    logical, intent(in) :: create_accessor

    character(len=*), parameter :: my_name = "File updater pair constructor"
    type(string_t) :: std_units, var_name
    class(property_t), pointer :: prop
    type(domain_target_cells_t) :: all_cells

    allocate( new_obj )
    allocate( new_obj%variable_, source = variable )

    ! Find or create domain variables
    if( .not. do_match( domain, new_obj%variable_ ) ) then
      deallocate( new_obj )
      new_obj => null( )
      return
    end if

    var_name = new_obj%variable_%musica_name( )
    std_units = domain%units( var_name%to_char( ) )
    prop => property_t( my_name,                                              &
                        name = var_name%to_char( ),                           & !- state variable name
                        units = std_units%to_char( ),                         & !- MUSICA units
                        applies_to = all_cells,                               & !- target domain
                        data_type = kDouble )                                   !- data type
    if( file%is_input( ) ) then
      new_obj%mutator_ => domain%mutator( prop )
    end if
    if( file%is_output( ) .or. create_accessor ) then
      new_obj%accessor_ => domain%accessor( prop )
    end if
    new_obj%musica_name_ = var_name
    deallocate( prop )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns a flag indicating whether the pair includes a given file
  !! variable
  logical function includes_variable( this, variable )

    use musica_assert,                 only : assert

    !> Paired variable
    class(file_paired_variable_t), intent(in) :: this
    !> File variable to compare with
    class(file_variable_t), intent(in) :: variable

    call assert( 113368224, associated( this%variable_ ) )
    includes_variable = this%variable_ .eq. variable

  end function includes_variable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a value from the file
  real(kind=musica_dk) function get_file_value( this, file, index )

    use musica_file,                   only : file_t

    !> Paired variable
    class(file_paired_variable_t), intent(inout) :: this
    !> IO file
    class(file_t), intent(inout) :: file
    !> Index in the temporal dimension to update from
    integer(kind=musica_ik), intent(in) :: index

    if( index .lt. this%first_staged_index_ .or.                             &
        index .gt. ( this%first_staged_index_ + this%number_staged_ ) - 1 )  &
      call this%update_staged_data( file, index )
    get_file_value = this%staged_data_( index - this%first_staged_index_ + 1 )

  end function get_file_value

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set a value in the file
  subroutine set_file_value( this, file, time__s, value )

    use musica_assert,                 only : assert
    use musica_file,                   only : file_t
    use musica_file_dimension_range,   only : file_dimension_range_t

    !> Paired variable
    class(file_paired_variable_t), intent(inout) :: this
    !> IO file
    class(file_t), intent(inout) :: file
    !> Simulation time [s]
    real(kind=musica_dk), intent(in) :: time__s
    !> Value to set in the file
    real(kind=musica_dk), intent(in) :: value

    type(file_dimension_range_t) :: indices(0)

    call assert( 886989346, associated( this%variable_ ) )
    call this%variable_%output( file, time__s, indices, value )

  end subroutine set_file_value

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a MUSICA domain state value
  real(kind=musica_dk) function get_musica_value( this, domain_state,         &
      iterator ) result( musica_value )

    use musica_assert,                 only : assert
    use musica_domain_state,           only : domain_state_t
    use musica_domain_iterator,        only : domain_iterator_t

    !> Paired variable
    class(file_paired_variable_t), intent(in) :: this
    !> Domain state
    class(domain_state_t), intent(in) :: domain_state
    !> Domain state iterator
    class(domain_iterator_t), intent(in) :: iterator

    call assert( 660447122, associated( this%accessor_ ) )
    call domain_state%get( iterator, this%accessor_, musica_value )

  end function get_musica_value

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set a MUSICA domain state value
  subroutine set_musica_value( this, domain_state, iterator, value )

    use musica_assert,                 only : assert
    use musica_domain_state,           only : domain_state_t
    use musica_domain_iterator,        only : domain_iterator_t

    !> Paired variable
    class(file_paired_variable_t), intent(in) :: this
    !> Domain state
    class(domain_state_t), intent(inout) :: domain_state
    !> Domain state iterator
    class(domain_iterator_t), intent(in) :: iterator
    !> New value
    real(kind=musica_dk), intent(in) :: value

    call assert( 633576938, associated( this%mutator_ ) )
    call domain_state%update( iterator, this%mutator_, value )

  end subroutine set_musica_value

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the MUSICA variable name
  function musica_name( this )

    !> MUSICA variable name
    type(string_t) :: musica_name
    !> Paired variable
    class(file_paired_variable_t), intent(in) :: this

    musica_name = this%musica_name_

  end function musica_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Prints the contents of the updater
  subroutine do_print( this )

    !> Paired variable
    class(file_paired_variable_t), intent(in) :: this

    call this%variable_%print( )

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Updates the staged data to start from a given index
  subroutine update_staged_data( this, file, index )

    use musica_assert,                 only : assert, assert_msg
    use musica_file,                   only : file_t
    use musica_file_dimension_range,   only : file_dimension_range_t

    !> Paired variable
    class(file_paired_variable_t), intent(inout) :: this
    !> IO file
    class(file_t), intent(inout) :: file
    !> New starting index
    integer(kind=musica_ik), intent(in) :: index

    integer(kind=musica_ik) :: n_times
    type(file_dimension_range_t), allocatable :: var_dims(:)

    var_dims = this%variable_%get_dimensions( )
    call assert_msg( 301384964, size( var_dims ) .eq. 1, "File variables "//  &
                     "are currently restricted to one dimension." )
    n_times = min( kMaxStagedData, var_dims(1)%upper_bound( ) - index + 1 )
    call file%check_open( )
    call assert( 382334001, index .gt. 0 .and. n_times .ge. 1 )
    this%staged_data_(:) = -huge( 1.0_musica_dk )
    call var_dims(1)%set( index, index + n_times - 1 )
    call this%variable_%get_data( file, var_dims(1:1), this%staged_data_ )
    this%first_staged_index_ = index
    this%number_staged_      = n_times

  end subroutine update_staged_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalizes a paired variable object
  elemental subroutine finalize( this )

    !> Paired variable
    type(file_paired_variable_t), intent(inout) :: this

    if( associated( this%mutator_  ) ) then
      deallocate( this%mutator_  )
      this%mutator_ => null( )
    end if
    if( associated( this%accessor_ ) ) then
      deallocate( this%accessor_ )
      this%accessor_ => null( )
    end if
    if( associated( this%variable_ ) ) then
      deallocate( this%variable_ )
      this%variable_ => null( )
    end if

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalizes a unique paired variable pointer
  elemental subroutine file_paired_variable_ptr_finalize ( this )

    !> Paired variable pointer
    type(file_paired_variable_ptr), intent(inout) :: this

    if( associated( this%val_ ) ) then
      deallocate( this%val_ )
      this%val_ => null( )
    end if

  end subroutine file_paired_variable_ptr_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Attemps to find a domain state variable for a given file variable. For
  !! certain file variables a domain state variable is created.
  !!
  logical function do_match( domain, variable )

    use musica_data_type,              only : kDouble
    use musica_domain,                 only : domain_t
    use musica_domain_target_cells,    only : domain_target_cells_t
    use musica_property,               only : property_t

    !> MUSICA domain
    class(domain_t), intent(inout) :: domain
    !> File variable
    class(file_variable_t), intent(inout) :: variable

    character(len=*), parameter :: my_name = "File variable matcher"
    type(string_t) :: musica_name, musica_units
    class(property_t), pointer :: prop
    type(domain_target_cells_t) :: all_cells

    musica_name = variable%musica_name( )
    do_match = domain%is_registered( musica_name%to_char( ) )
    if( do_match ) then
      musica_units = domain%units( musica_name%to_char( ) )
      call variable%set_musica_units( musica_units%to_char( ) )
    end if

  end function do_match

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_file_paired_variable
