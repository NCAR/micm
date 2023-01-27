! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_input_output_processor module

!> The input_output_processor_t type and related functions
module musica_input_output_processor

  use musica_constants,                only : musica_ik, musica_dk
  use musica_file,                     only : file_t
  use musica_file_dimension,           only : file_dimension_t
  use musica_file_updater,             only : file_updater_ptr
  use musica_domain_iterator,          only : domain_iterator_t

  implicit none
  private

  public :: input_output_processor_t, input_output_processor_ptr

  !> Input/Output processor
  !!
  !! One input/output processor should be set up for each source/destination
  !! and can be used at runtime to load model data from an external source
  !! or output model data to an external destination. All conversions to/from
  !! standard MUSICA units, interpolation, and user-specified data
  !! transformations are handled by the input/output processor.
  !!
  !! \todo add example usage for I/O processors
  !!
  type :: input_output_processor_t
    private
    !> File attributes and functions
    class(file_t), pointer :: file_ => null( )
    !> Time dimension
    class(file_dimension_t), pointer :: time_ => null( )
    !> Last time index used
    integer(kind=musica_ik) :: last_time_index_ = 1
    !> Updaters for paired linear combinations of MUSICA and I/O variables
    type(file_updater_ptr), allocatable :: linear_combination_updaters_(:)
    !> Updaters for paired MUSICA <-> I/O variables
    type(file_updater_ptr), allocatable :: single_variable_updaters_(:)
    !> Iterator over all domain cells
    class(domain_iterator_t), pointer :: iterator_ => null( )
  contains
    !> Registers a state variable for output
    procedure :: register_output_variable
    !> Gets the times corresponding to entries (for input data) [s]
    procedure :: entry_times__s
    !> Updates the model state with input data
    procedure :: update_state
    !> Returns the names of all MUSICA variables set/read by the I/O processor
    procedure :: musica_variable_names
    !> Outputs the current domain state
    procedure :: output
    !> Preprocesses input/output processor configuration
    procedure :: preprocess_input
    !> Print the input/output configuration information
    procedure :: print => do_print
    !> Loads input file variable names and units
    procedure, private :: load_input_variables
    !> Load linear combinations of variables
    procedure, private :: load_linear_combinations
    !> Returns whether a specified variable is included in the set of linear
    !! combination updaters
    procedure, private :: is_used_in_linear_combination
    !> Finalizes the input/output processor
    final :: finalize
  end type input_output_processor_t

  !> Unique pointer to input_output_processor_t objects
  type :: input_output_processor_ptr
    class(input_output_processor_t), pointer :: val_ => null( )
  contains
    !> Finalizes the pointer
    final :: input_output_processor_ptr_finalize
  end type input_output_processor_ptr

  !> Constructor
  interface input_output_processor_t
    module procedure :: constructor
  end interface input_output_processor_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates a connection to an input or output file
  !!
  !! At minimum, \c config must include a top-level key-value pair "intent",
  !! which can be either "input" or "output".
  !!
  !! A "file name" is also required for files openned for input. This is
  !! optional for output files, with the default name starting with "output"
  !! and having an appropriate file extension.
  !!
  !! Input files require the model domain object for mapping between model
  !! domain and file variables.
  !!
  function constructor( config, domain ) result( new_obj )

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_domain,                 only : domain_t
    use musica_file_dimension_factory, only : file_dimension_builder
    use musica_file_factory,           only : file_builder
    use musica_string,                 only : string_t

    !> New input/output processor
    type(input_output_processor_t), pointer :: new_obj
    !> Configuration data
    type(config_t), intent(inout) :: config
    !> Model domain
    class(domain_t), intent(inout), optional :: domain

    character(len=*), parameter :: my_name = 'I/O processor constructor'
    type(string_t) :: temp_str
    logical :: found, is_input
    type(config_t) :: linear_combos, time_config

    allocate( new_obj )

    ! set up the file
    new_obj%file_ => file_builder( config )

    ! load the variable names and units from input files
    if( new_obj%file_%is_input( ) ) then
      if( present( domain ) ) then
        call new_obj%load_input_variables( domain, config )
      else
        call die_msg( 264651720, "Input files require the model domain "//    &
                                 "during initialization for mapping." )
      end if
    else
      call time_config%empty( )
      call time_config%add( "type", new_obj%file_%type( ), my_name )
      new_obj%time_ => file_dimension_builder( time_config, new_obj%file_,    &
                                               dimension_name = "time" )
      allocate( new_obj%single_variable_updaters_( 0 ) )
      allocate( new_obj%linear_combination_updaters_( 0 ) )
    end if

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Registers a state variable for output
  subroutine register_output_variable( this, domain, domain_variable_name,    &
      units, io_variable_name )

    use musica_assert,                 only : assert, assert_msg
    use musica_config,                 only : config_t
    use musica_domain,                 only : domain_t
    use musica_file_dimension_range,   only : file_dimension_range_t
    use musica_file_updater,           only : file_updater_t
    use musica_file_variable,          only : file_variable_t
    use musica_file_variable_factory,  only : file_variable_builder
    use musica_string,                 only : string_t

    !> Input/output
    class(input_output_processor_t), intent(inout) :: this
    !> Model domain
    class(domain_t), intent(inout) :: domain
    !> Variable to register
    character(len=*), intent(in) :: domain_variable_name
    !> Units to use for the intput/output data
    !!
    !! Conversions will be handled by the input/output processor
    character(len=*), intent(in) :: units
    !> Optional custom name to use in file
    !!
    !! If not included, the standard MUSICA name will be used.
    character(len=*), intent(in), optional :: io_variable_name

    character(len=*), parameter :: my_name = "Output variable registrar"
    type(string_t) :: file_name, io_var_name
    type(config_t) :: variable_config, temp_config
    class(file_variable_t), pointer :: new_var
    type(file_updater_t), pointer :: updater
    type(file_updater_ptr), allocatable :: temp_ptrs(:)
    type(file_updater_ptr) :: new_ptr
    integer(kind=musica_ik) :: i_ptr, n_ptrs
    logical :: is_match
    ! there is currently only one dimension for output variables (time)
    type(file_dimension_range_t) :: dims(1)

    call assert( 181665946, associated( this%time_ ) )
    dims(1) = this%time_%get_range( )
    file_name = this%file_%name( )
    call assert_msg( 567596084, this%file_%is_output( ),                      &
                     "Cannot register output of '"//                          &
                     domain_variable_name//"' to input file '"//              &
                     file_name%to_char( )//"'" )
    call variable_config%add( "units", units, my_name )
    call variable_config%add( "MusicBox name", domain_variable_name, my_name )
    call temp_config%add( io_variable_name, variable_config, my_name )
    call variable_config%add( "properties", temp_config, my_name )
    call variable_config%add( "type", this%file_%type( ), my_name )
    if( present( io_variable_name ) ) then
      io_var_name = io_variable_name
    else
      io_var_name = domain_variable_name
    end if
    new_var => file_variable_builder( variable_config, this%file_,            &
                                      variable_name = io_var_name%to_char( ), &
                                      dimensions = dims )
    updater => file_updater_t( this%file_, domain, new_var )
    call assert_msg( 987906252, associated( updater ),                        &
                     "Could not find domain state variable '"//               &
                     domain_variable_name//"' to register as '"//             &
                     io_var_name%to_char( )//"' in output file '"//           &
                     file_name%to_char( )//"'" )
    new_ptr%val_ => updater
    n_ptrs = size( this%single_variable_updaters_ )
    allocate( temp_ptrs( n_ptrs ) )
    do i_ptr = 1, n_ptrs
      temp_ptrs( i_ptr )%val_ => this%single_variable_updaters_( i_ptr )%val_
      this%single_variable_updaters_( i_ptr )%val_ => null( )
    end do
    deallocate( this%single_variable_updaters_ )
    allocate( this%single_variable_updaters_( n_ptrs + 1 ) )
    do i_ptr = 1, n_ptrs
      this%single_variable_updaters_( i_ptr )%val_ => temp_ptrs( i_ptr )%val_
      temp_ptrs( i_ptr )%val_ => null( )
    end do
    this%single_variable_updaters_( n_ptrs + 1 )%val_ => new_ptr%val_
    new_ptr%val_ => null( )
    deallocate( temp_ptrs )
    deallocate( new_var )

  end subroutine register_output_variable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the times corresponding to entries (for input data) [s]
  !!
  !! These times include any adjustments specified in the configuration data
  !! and should thus correspond directly to the MUSICA simulation time.
  function entry_times__s( this )

    use musica_assert,                 only : assert

    !> Entry times [s]
    real(kind=musica_dk), allocatable :: entry_times__s(:)
    !> Input/output
    class(input_output_processor_t), intent(inout) :: this

    call assert( 579211287, associated( this%file_ ) )
    call assert( 188228761, this%file_%is_input( ) )
    call assert( 465439703, associated( this%time_ ) )
    entry_times__s = this%time_%get_values( )

  end function entry_times__s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Updates the model state with input data
  !!
  !! If a time is included, input data for the specified time (with any
  !! necessary interpolation) will be used to update domain state variables
  !! registered with the \c input_output_processor_t type during intialization.
  !!
  !! If no time is provided the first entry in the input data will be used
  !! to update the domain state (used for initial conditions).
  !!
  subroutine update_state( this, domain, domain_state, time__s, tethered_only )

    use musica_assert,                 only : assert
    use musica_domain,                 only : domain_t
    use musica_domain_state,           only : domain_state_t
    use musica_domain_target_cells,    only : domain_target_cells_t

    !> Input/output
    class(input_output_processor_t), intent(inout) :: this
    !> Model domain
    class(domain_t), intent(in) :: domain
    !> Domain state to update
    class(domain_state_t), intent(inout) :: domain_state
    !> Current simulation time [s]
    real(kind=musica_dk), intent(in), optional :: time__s
    !> Optional flag indicating whether to update for all or only tethered
    !! file variables (defaults to false)
    logical, intent(in), optional :: tethered_only

    integer(kind=musica_ik) :: i_data, i_updater
    logical :: found, l_tethered_only
    type(domain_target_cells_t) :: all_cells

    call assert( 967287335, associated( this%file_ ) )
    call assert( 682224647, this%file_%is_input( ) )
    i_data = 1
    l_tethered_only = .false.
    if( present( tethered_only ) ) l_tethered_only = tethered_only
    if( present( time__s ) ) then
      i_data = this%time_%get_index( time__s, is_exact = found,               &
                                     guess = this%last_time_index_ )
      if( .not. found ) return
    end if
    if( .not. associated( this%iterator_ ) ) then
      this%iterator_ => domain%iterator( all_cells )
    end if
    do i_updater = 1, size( this%linear_combination_updaters_ )
      associate( updater =>                                                   &
                 this%linear_combination_updaters_( i_updater )%val_ )
      if( .not. updater%is_tethered( ) .and. l_tethered_only ) cycle
      call updater%update_state( this%file_, i_data, this%iterator_,          &
                                 domain_state )
      end associate
    end do
    do i_updater = 1, size( this%single_variable_updaters_ )
      associate( updater => this%single_variable_updaters_( i_updater )%val_ )
      if( .not. updater%is_tethered( ) .and. l_tethered_only ) cycle
      call updater%update_state( this%file_, i_data, this%iterator_,          &
                                 domain_state )
      end associate
    end do

  end subroutine update_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the names of all MUSICA variables read/set by the I/O processor
  function musica_variable_names( this )

    use musica_string,                 only : string_t

    !> MUSICA variable names
    type(string_t), allocatable :: musica_variable_names(:)
    !> Input/output
    class(input_output_processor_t), intent(in) :: this

    integer(kind=musica_ik) :: i_updater
    type(string_t), allocatable :: var_names(:)

    allocate( musica_variable_names( 0 ) )
    if( allocated( this%single_variable_updaters_ ) ) then
      do i_updater = 1, size( this%single_variable_updaters_ )
        associate( updater =>                                                 &
                   this%single_variable_updaters_( i_updater )%val_ )
        var_names = updater%musica_variable_names( )
        musica_variable_names = [ musica_variable_names, var_names ]
        end associate
      end do
    end if
    if( allocated( this%linear_combination_updaters_ ) ) then
      do i_updater = 1, size( this%linear_combination_updaters_ )
        associate( updater =>                                                 &
                   this%linear_combination_updaters_( i_updater )%val_ )
        var_names = updater%musica_variable_names( )
        musica_variable_names = [ musica_variable_names, var_names ]
        end associate
      end do
    end if

  end function musica_variable_names

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Outputs the current domain state
  !!
  !! Domain state variables registered with the \c input_output_processor_t
  !! type during initialization will be output for the current simulation
  !! time after conversion to the specified output units.
  !!
  subroutine output( this, time__s, domain, domain_state )

    use musica_assert,                 only : die_msg
    use musica_domain,                 only : domain_t
    use musica_domain_state,           only : domain_state_t
    use musica_domain_target_cells,    only : domain_target_cells_t

    !> Input/output
    class(input_output_processor_t), intent(inout) :: this
    !> Current simulation time [s]
    real(kind=musica_dk), intent(in) :: time__s
    !> Model domain
    class(domain_t), intent(in) :: domain
    !> Domain state
    class(domain_state_t), intent(in) :: domain_state

    integer(kind=musica_ik) :: i_updater
    type(domain_target_cells_t) :: all_cells

    if( .not. associated( this%iterator_ ) ) then
      this%iterator_ => domain%iterator( all_cells )
    end if

    do i_updater = 1, size( this%single_variable_updaters_ )
      associate( updater => this%single_variable_updaters_( i_updater )%val_ )
      call updater%output( this%file_, time__s, domain, domain_state,         &
                           this%iterator_ )
      end associate
    end do

  end subroutine output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Preprocesses the input/output processor configuration
  subroutine preprocess_input( this, config, default_variability )

    use musica_assert,                 only : assert
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    !> Input/output
    class(input_output_processor_t), intent(in) :: this
    !> Input/output configuration data
    type(config_t), intent(out) :: config
    !> Default variability to assume for the file variables
    !! "prognosed", "tethered", or "fixed"
    character(len=*), intent(in) :: default_variability

    character(len=*), parameter :: my_name = "I/O preprocessor"
    integer(kind=musica_ik) :: i_updater, i_var
    type(string_t) :: linear_combo_name
    type(string_t), allocatable :: var_names(:)
    type(config_t) :: linear_combo, linear_combos, vars, empty, tethered
    logical :: default_tethered, found_custom_var

    call config%empty( )
    call empty%empty( )
    call tethered%empty( )
    call tethered%add( "variability", "tethered", my_name )
    default_tethered = default_variability .eq. "tethered"
    found_custom_var = .false.
    if( allocated( this%single_variable_updaters_ ) ) then
      do i_updater = 1, size( this%single_variable_updaters_ )
        associate( updater => this%single_variable_updaters_( i_updater ) )
        call assert( 994514267, associated( updater%val_ ) )
        var_names = updater%val_%musica_variable_names( )
        call assert( 915783652, size( var_names ) .eq. 1 )
        if( updater%val_%is_tethered( ) .and. .not. default_tethered ) then
          if( .not. found_custom_var ) then
            call vars%empty( )
            found_custom_var = .true.
          end if
          call vars%add( var_names( 1 )%to_char( ), tethered, my_name )
        end if
        end associate
      end do
      if( found_custom_var ) call config%add( "properties", vars, my_name )
    end if
    if( allocated( this%linear_combination_updaters_ ) ) then
      if( size( this%linear_combination_updaters_ ) .gt. 0 ) then
        call linear_combos%empty( )
        do i_updater = 1, size( this%linear_combination_updaters_ )
          associate( updater =>                                               &
                     this%linear_combination_updaters_( i_updater ) )
          call assert( 931411401, associated( updater%val_ ) )
          var_names = updater%val_%musica_variable_names( )
          call vars%empty( )
          do i_var = 1, size( var_names )
            call vars%add( var_names( i_var )%to_char( ), empty, my_name )
          end do
          call linear_combo%empty( )
          call linear_combo%add( "properties", vars, my_name )
          linear_combo_name = updater%val_%name( )
          call linear_combos%add( linear_combo_name%to_char( ), linear_combo, &
                                  my_name )
          end associate
        end do
        call config%add( "linear combinations", linear_combos, my_name )
      end if
    end if

  end subroutine preprocess_input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the input/output configuration information
  subroutine do_print( this )

    !> Input/output
    class(input_output_processor_t), intent(in) :: this

    integer(kind=musica_ik) :: i

    write(*,*) "***** Input/Output Processor Configuration *****"
    write(*,*) ""
    call this%file_%print( )
    write(*,*) ""
    call this%time_%print( )
    write(*,*) ""
    write(*,*) "---------------------"
    write(*,*) " State/Output Updaters"
    write(*,*) "---------------------"
    write(*,*) ""
    if( allocated( this%linear_combination_updaters_ ) ) then
      do i = 1, size( this%linear_combination_updaters_ )
        call this%linear_combination_updaters_( i )%val_%print( )
      end do
    end if
    if( allocated( this%single_variable_updaters_ ) ) then
      do i = 1, size( this%single_variable_updaters_ )
        call this%single_variable_updaters_( i )%val_%print( )
      end do
    end if
    write(*,*) ""
    write(*,*) "***** End Input/Output Processor Configuration *****"

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Loads input file variable names, units, and ids
  subroutine load_input_variables( this, domain, config )

    use musica_assert,                 only : assert, assert_msg
    use musica_config,                 only : config_t
    use musica_domain,                 only : domain_t
    use musica_file_dimension_factory, only : file_dimension_builder
    use musica_file_updater,           only : file_updater_t
    use musica_file_variable,          only : file_variable_ptr
    use musica_file_variable_factory,  only : file_variable_builder
    use musica_string

    !> Input/Output processor
    class(input_output_processor_t), intent(inout) :: this
    !> MUSICA domain
    class(domain_t), intent(inout) :: domain
    !> Input/output configuration
    type(config_t), intent(inout) :: config

    character(len=*), parameter :: my_name = "Input file variable loader"
    logical :: is_match, found
    type(string_t) :: file_name
    type(file_updater_t), pointer :: updater
    integer(kind=musica_ik) :: n_dims, n_vars, i_var, n_match, i_match,       &
                               i_updater
    type(file_variable_ptr), allocatable :: vars(:)
    type(config_t) :: linear_combos

    ! get linear combinations first
    call config%get( "linear combinations", linear_combos, my_name,           &
                     found = found )
    if( found ) then
      call this%load_linear_combinations( domain, linear_combos )
    else
      allocate( this%linear_combination_updaters_( 0 ) )
    end if

    call this%file_%check_open( )
    file_name = this%file_%name( )
    n_dims    = this%file_%number_of_dimensions( )
    n_vars    = this%file_%number_of_variables( )
    allocate( vars( n_vars ) )
    n_match = 0
    do i_var = 1, n_vars
      vars( i_var )%val_ =>                                                   &
          file_variable_builder( config, this%file_, variable_id = i_var )
      updater => file_updater_t( this%file_, domain, vars( i_var )%val_ )
      if( associated( updater ) .and.                                         &
          .not. this%is_used_in_linear_combination( vars( i_var )%val_ ) ) then
        deallocate( updater )
        n_match = n_match + 1
        cycle
      end if
      if( associated( updater ) ) deallocate( updater )
      if( vars( i_var )%val_%musica_name( ) .eq. "time" ) cycle
      deallocate( vars( i_var )%val_ )
      vars( i_var )%val_ => null( )
    end do
    call assert( 177225698, .not. allocated( this%single_variable_updaters_ ) )
    allocate( this%single_variable_updaters_( n_match ) )
    i_match = 0
    do i_var = 1, n_vars
      if( associated( vars( i_var )%val_ ) ) then
        if( vars( i_var )%val_%musica_name( ) .eq. "time" ) then
          call vars( i_var )%val_%set_musica_units( "s" )
          this%time_ =>                                                       &
              file_dimension_builder( config, this%file_, vars( i_var )%val_ )
          cycle
        end if
        i_match = i_match + 1
        this%single_variable_updaters_( i_match )%val_ =>                     &
            file_updater_t( this%file_, domain, vars( i_var )%val_ )
      end if
    end do
    do i_var = 1, n_vars
      if( associated( vars( i_var )%val_ ) ) deallocate( vars( i_var )%val_ )
    end do
    call assert( 177448848, i_match .eq. n_match )
    deallocate( vars )

  end subroutine load_input_variables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load linear combinations of variables
  subroutine load_linear_combinations( this, domain, config )

    use musica_assert,                 only : assert, die_msg
    use musica_config,                 only : config_t
    use musica_domain,                 only : domain_t
    use musica_file_updater,           only : file_updater_t
    use musica_iterator,               only : iterator_t
    use musica_string,                 only : string_t

    !> Input/Output processor
    class(input_output_processor_t), intent(inout) :: this
    !> MUSICA domain
    class(domain_t), intent(inout) :: domain
    !> Linear combinations configuration
    type(config_t), intent(inout) :: config

    character(len=*), parameter :: my_name =                                  &
        "Input file linear combination loader"
    class(iterator_t), pointer :: iter
    type(config_t) :: linear_combo_config
    integer(kind=musica_ik) :: i_linear_combo, n_linear_combos
    type(string_t) :: combo_name

    iter => config%get_iterator( )
    n_linear_combos = 0
    do while( iter%next( ) )
      n_linear_combos = n_linear_combos + 1
    end do
    call assert( 480140708,                                                   &
                 .not. allocated( this%linear_combination_updaters_ ) )
    allocate( this%linear_combination_updaters_( n_linear_combos ) )
    call iter%reset( )
    i_linear_combo = 0
    do while( iter%next( ) )
      i_linear_combo = i_linear_combo + 1
      associate( lc_ptr =>                                                    &
                 this%linear_combination_updaters_( i_linear_combo ) )
      call config%get( iter, linear_combo_config, my_name )
      combo_name = config%key( iter )
      lc_ptr%val_ => file_updater_t( combo_name, this%file_, domain,          &
                                     linear_combo_config )
      if( .not. associated( lc_ptr%val_ ) ) then
        call linear_combo_config%print( )
        call die_msg( 110789656, "Could not create linear combination." )
      end if
      end associate
    end do
    call assert( 101427158, i_linear_combo .eq. n_linear_combos )
    deallocate( iter )

  end subroutine load_linear_combinations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns whether a variable is used in the set of linear combination
  !! updaters
  logical function is_used_in_linear_combination( this, variable )            &
      result( is_used )

    use musica_assert,                 only : assert
    use musica_file_variable,          only : file_variable_t

    !> Input/Output processor
    class(input_output_processor_t), intent(in) :: this
    !> File variable
    class(file_variable_t), intent(in) :: variable

    integer(kind=musica_ik) :: i_updater

    is_used = .false.
    call assert( 998495082, allocated( this%linear_combination_updaters_ ) )
    do i_updater = 1, size( this%linear_combination_updaters_ )
      associate( lc => this%linear_combination_updaters_( i_updater ) )
      call assert( 428732572, associated( lc%val_ ) )
      if( lc%val_%includes_variable( variable ) ) then
        is_used = .true.
        return
      end if
      end associate
    end do

  end function is_used_in_linear_combination

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalizes an input/output processor
  subroutine finalize( this )

    !> Input/Output processor
    type(input_output_processor_t), intent(inout) :: this

    integer(kind=musica_ik) :: i_updater

    if( associated( this%file_ ) ) then
      deallocate( this%file_ )
      this%file_ => null( )
    end if
    if( associated( this%time_ ) ) then
      deallocate( this%time_ )
      this%file_ => null ( )
    end if
    if( associated( this%iterator_ ) ) then
      deallocate( this%iterator_ )
      this%iterator_ => null( )
    end if

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalizes a unique input/output pointer
  elemental subroutine input_output_processor_ptr_finalize( this )

    !> Input/output pointer
    type(input_output_processor_ptr), intent(inout) :: this

    if( associated( this%val_ ) ) then
      deallocate( this%val_ )
      this%val_ => null( )
    end if

  end subroutine input_output_processor_ptr_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_input_output_processor
