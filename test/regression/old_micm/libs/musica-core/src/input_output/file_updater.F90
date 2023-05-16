! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_file_updater module

!> The file_updater_t type and related functions
module musica_file_updater

  use musica_constants,                only : musica_dk, musica_ik
  use musica_file_paired_variable,     only : file_paired_variable_ptr
  use musica_string,                   only : string_t

  implicit none
  private

  public :: file_updater_t, file_updater_ptr

  !> Updater for one or more paired MUSICA <-> file variables
  !!
  !! Functions for updating MUSICA state variables from input data and
  !! updating output files from the MUSICA state.
  !!
  type :: file_updater_t
    private
    !> Name for the updater
    type(string_t) :: name_
    !> Set of paired variables
    type(file_paired_variable_ptr), allocatable :: pairs_(:)
    !> Scaling factor (for linear combinations of variables)
    real(kind=musica_dk) :: scale_factor_ = 1.0_musica_dk
    !> Working array for updating MUSICA/file variables
    real(kind=musica_dk), allocatable :: working_file_values_(:)
    !> Working array for updating MUSICA/file variables
    real(kind=musica_dk), allocatable :: working_musica_values_(:)
    !> Flag indicating whether the variable is tethered
    !! (Tethered variables are allowed to vary during the simulation, but are
    !!  reset to the input value at the end of every time step.)
    logical :: is_tethered_ = .false.
  contains
    !> Returns the name of the updater
    procedure :: name => updater_name
    !> Returns the standard MUSICA file name for the variable
    !> Indicates whether a specified variable is used in the updater
    procedure :: includes_variable
    !> Returns whether the paired variable(s) is (are) tethered
    procedure :: is_tethered
    !> Updates the state for a given index in the temporal dimension
    procedure :: update_state
    !> Returns the names of the MUSICA variables read/set by the updater
    procedure :: musica_variable_names
    !> Outputs data to the file
    procedure :: output
    !> Prints the properties of the updater
    procedure :: print => do_print
  end type file_updater_t

  !> Constructor
  interface file_updater_t
    module procedure :: constructor_single_variable
    module procedure :: constructor_linear_combination
  end interface file_updater_t

  !> Unique pointer to file_updater_t objects
  type :: file_updater_ptr
    type(file_updater_t), pointer :: val_ => null( )
  contains
    !> Finalizes the pointer
    final :: file_updater_ptr_finalize
  end type file_updater_ptr

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates an updater for a single paired MUSICA/file variable
  !!
  !! If the MUSICA domain variable cannot be found, a null pointer is
  !! returned.
  !!
  function constructor_single_variable( file, domain, variable )              &
      result( new_obj )

    use musica_domain,                 only : domain_t
    use musica_file,                   only : file_t
    use musica_file_paired_variable,   only : file_paired_variable_t
    use musica_file_variable,          only : file_variable_t

    !> New MUSICA<->File variable match
    type(file_updater_t), pointer :: new_obj
    !> File to update to or from
    class(file_t), intent(inout) :: file
    !> Model domain
    class(domain_t), intent(inout) :: domain
    !> File variable
    class(file_variable_t), intent(in) :: variable

    type(file_paired_variable_t), pointer :: pair

    new_obj => null( )
    pair => file_paired_variable_t( file, domain, variable, .false. )
    if( .not. associated( pair ) ) return
    allocate( new_obj )
    allocate( new_obj%pairs_( 1 ) )
    allocate( new_obj%working_file_values_( 1 ) )
    allocate( new_obj%working_musica_values_( 1 ) )
    new_obj%pairs_( 1 )%val_ => pair
    new_obj%is_tethered_ = variable%is_tethered( )
    new_obj%name_        = variable%musica_name( )

  end function constructor_single_variable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates an updater for a linear combination of MUSICA/file variables
  !!
  !! If any of the specified variables cannot be found in the file or the
  !! MUSICA domain, a null pointer is returned.
  !!
  function constructor_linear_combination( name, file, domain, config )       &
      result( new_obj )

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_domain,                 only : domain_t
    use musica_file,                   only : file_t
    use musica_file_paired_variable,   only : file_paired_variable_t
    use musica_file_variable,          only : file_variable_t
    use musica_file_variable_factory,  only : file_variable_builder
    use musica_iterator,               only : iterator_t

    !> New updater for MUSICA<->File linear combination of variables
    type(file_updater_t), pointer :: new_obj
    !> Name for the linear combination
    type(string_t), intent(in) :: name
    !> File to update from
    class(file_t), intent(inout) :: file
    !> Model domain
    class(domain_t), intent(inout) :: domain
    !> Linear combination configuration
    type(config_t), intent(inout) :: config

    character(len=*), parameter :: my_name =                                  &
        "File linear combination updater constructor"
    type(config_t) :: props, var_config
    type(string_t) :: var_name
    class(iterator_t), pointer :: iter
    class(file_variable_t), pointer :: var
    type(file_paired_variable_t), pointer :: pair
    integer(kind=musica_ik) :: i_var, n_vars

    call config%get( "properties", props, my_name )
    iter => props%get_iterator( )
    n_vars = 0
    do while( iter%next( ) )
      n_vars = n_vars + 1
    end do
    allocate( new_obj )
    allocate( new_obj%pairs_( n_vars ) )
    allocate( new_obj%working_file_values_( n_vars ) )
    allocate( new_obj%working_musica_values_( n_vars ) )
    new_obj%is_tethered_ = .true.
    new_obj%name_        = name
    call config%get( "scale factor", new_obj%scale_factor_, my_name,          &
                     default = new_obj%scale_factor_ )
    call iter%reset( )
    i_var = 0
    do while( iter%next( ) )
      i_var = i_var + 1
      call props%get( iter, var_config, my_name )
      var_name = props%key( iter )
      call var_config%add( "type", file%type( ), my_name )
      call var_config%add( "default variability", "tethered", my_name )
      var => file_variable_builder( var_config, file, var_name%to_char( ) )
      call assert_msg( 659534542, associated( var ), "Could not find "//      &
                       "MUSICA domain variable '"//var_name%to_char( )//"'" )
      pair => file_paired_variable_t( file, domain, var, .true. )
      if( .not. associated( pair ) ) then
        deallocate( new_obj )
        new_obj => null( )
        exit
      end if
      new_obj%pairs_( i_var )%val_ => pair
      deallocate( var )
    end do
    deallocate( iter )

  end function constructor_linear_combination

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the name of the updater
  type(string_t) function updater_name( this )

    !> File updater
    class(file_updater_t), intent(in) :: this

    updater_name = this%name_

  end function updater_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Indicates whether a specified variable is used in the updater
  logical function includes_variable( this, variable )

    use musica_assert,                 only : assert
    use musica_file_variable,          only : file_variable_t

    !> File updater
    class(file_updater_t), intent(in) :: this
    !> Variable to look for in updater
    class(file_variable_t), intent(in) :: variable

    integer(kind=musica_ik) :: i_pair

    call assert( 684123688, allocated( this%pairs_ ) )
    includes_variable = .false.
    do i_pair = 1, size( this%pairs_ )
      if( this%pairs_( i_pair )%val_%includes_variable( variable ) ) then
        includes_variable = .true.
        return
      end if
    end do

  end function includes_variable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns whether the variable(s) is (are) tethered to input conditions
  logical function is_tethered( this )

    !> File updater
    class(file_updater_t), intent(in) :: this

    is_tethered = this%is_tethered_

  end function is_tethered

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Updates a domain state for a given index in the temporal dimension
  subroutine update_state( this, file, index, iterator, state )

    use musica_assert,                 only : assert
    use musica_domain_state,           only : domain_state_t
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_file,                   only : file_t

    !> File updater
    class(file_updater_t), intent(inout) :: this
    !> File file
    class(file_t), intent(inout) :: file
    !> Index in the temporal dimension to update from
    integer(kind=musica_ik), intent(in) :: index
    !> Domain state iterator
    class(domain_iterator_t), intent(inout) :: iterator
    !> Domain state to update
    class(domain_state_t), intent(inout) :: state

    integer(kind=musica_ik) :: i_pair
    real(kind=musica_dk) :: file_total, musica_total, new_value

    call assert( 381356243, allocated( this%pairs_ ) )
    call iterator%reset( )
    do while( iterator%next( ) )
      if( size( this%pairs_ ) .eq. 1 ) then
        associate( pair => this%pairs_( 1 )%val_ )
        call pair%set_musica_value( state, iterator,                          &
                                    pair%get_file_value( file, index ) )
        end associate
        cycle
      end if
      file_total   = 0.0_musica_dk
      musica_total = 0.0_musica_dk
      do i_pair = 1, size( this%pairs_ )
        associate( pair => this%pairs_( i_pair )%val_ )
        this%working_file_values_( i_pair ) =                                 &
            pair%get_file_value( file, index ) * this%scale_factor_
        this%working_musica_values_( i_pair ) =                               &
            pair%get_musica_value( state, iterator )
        file_total   = file_total   + this%working_file_values_( i_pair )
        musica_total = musica_total + this%working_musica_values_( i_pair )
        end associate
      end do
      if( musica_total .eq. 0.0_musica_dk ) then
        do i_pair = 1, size( this%pairs_ )
          associate( pair => this%pairs_( i_pair )%val_ )
          new_value = this%working_file_values_( i_pair )
          call pair%set_musica_value( state, iterator, new_value )
          end associate
        end do
      else
        do i_pair = 1, size( this%pairs_ )
          associate( pair => this%pairs_( i_pair )%val_ )
          new_value = this%working_musica_values_( i_pair ) / musica_total    &
                      * file_total
          call pair%set_musica_value( state, iterator, new_value )
          end associate
        end do
      end if
    end do

  end subroutine update_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the names of the MUSICA variables read/set by the updater
  function musica_variable_names( this )

    use musica_assert,                 only : assert

    !> MUSICA variable names
    type(string_t), allocatable :: musica_variable_names(:)
    !> File updater
    class(file_updater_t), intent(in) :: this

    integer(kind=musica_ik) :: i_pair

    call assert( 973607705, allocated( this%pairs_ ) )
    allocate( musica_variable_names( size( this%pairs_ ) ) )
    do i_pair = 1, size( this%pairs_ )
      musica_variable_names( i_pair ) =                                       &
          this%pairs_( i_pair )%val_%musica_name( )
    end do

  end function musica_variable_names

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Outputs data to a file
  subroutine output( this, file, time__s, domain, domain_state, iterator )

    use musica_assert,                 only : assert, assert_msg
    use musica_domain,                 only : domain_t
    use musica_domain_state,           only : domain_state_t
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_file,                   only : file_t
    use musica_file_dimension_range,   only : file_dimension_range_t

    !> File updater
    class(file_updater_t), intent(inout) :: this
    !> Output file
    class(file_t), intent(inout) :: file
    !> Current simulation time [s]
    real(kind=musica_dk), intent(in) :: time__s
    !> Model domain
    class(domain_t), intent(in) :: domain
    !> Domain state
    class(domain_state_t), intent(in) :: domain_state
    !> Domain iterator
    class(domain_iterator_t), intent(inout) :: iterator

    real(kind=musica_dk) :: state_value

    call assert( 243308233, allocated( this%pairs_ ) )
    call assert( 632837520, size( this%pairs_ ) .eq. 1 )
    call iterator%reset( )
    if( iterator%next( ) ) then
      state_value = this%pairs_( 1 )%val_%get_musica_value( domain_state,     &
                                                            iterator )
      call this%pairs_( 1 )%val_%set_file_value( file, time__s, state_value )
    end if
    call assert_msg( 608265274, .not. iterator%next( ), "Output files are "// &
                     "not yet set up for multiple cells." )

  end subroutine output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Prints the contents of the updater
  subroutine do_print( this )

    !> File updater
    class(file_updater_t), intent(in) :: this

    integer(kind=musica_ik) :: i_pair

    write(*,*) "*** Updater for paired MUSICA/file variables ***"
    if( allocated( this%pairs_ ) ) then
      do i_pair = 1, size( this%pairs_ )
        call this%pairs_( i_pair )%val_%print( )
      end do
    end if
    write(*,*) "*** End file updater ***"

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalizes a unique file updater pointer
  elemental subroutine file_updater_ptr_finalize( this )

    !> File updater
    type(file_updater_ptr), intent(inout) :: this

    if( associated( this%val_ ) ) then
      deallocate( this%val_ )
      this%val_ => null( )
    end if

  end subroutine file_updater_ptr_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_file_updater
