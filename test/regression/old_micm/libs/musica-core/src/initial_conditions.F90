! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_initial_conditions module

!> set_initial_conditions and related functions
module musica_initial_conditions

  use musica_constants,                only : musica_dk, musica_ik
  use musica_domain_state_mutator,     only : domain_state_mutator_t
  use musica_input_output_processor,   only : input_output_processor_ptr
  use musica_string

  implicit none
  private

  public :: initial_conditions_t

  !> Initial value for one domain state variable
  type :: initial_value_t
    private
    !> Domain state mutator
    class(domain_state_mutator_t), pointer :: mutator_ => null( )
    !> Initial value
    real(kind=musica_dk) :: value_ = -9.99999e300_musica_dk
  contains
    !> Sets the domain state variable value
    procedure :: set_state => initial_value_set_state
    !> Finalizes the object
    final :: initial_value_finalize
  end type initial_value_t

  !> Initial model state conditions
  type :: initial_conditions_t
    private
    !> Model variables affected by the initial conditions
    type(string_t), allocatable :: variable_names_(:)
    !> Input files
    type(input_output_processor_ptr), allocatable :: files_(:)
    !> Initial conditions specified in the configuration
    type(initial_value_t), allocatable :: initial_values_(:)
  contains
    !> Creates a model state loaded with the initial conditions
    procedure :: get_state
    !> Update a state for "tethered" variables from initial conditions
    procedure :: update_state
    !> Preprocesses the initial conditions for repeat model runs
    procedure :: preprocess_input
    !> Load initial conditions
    procedure, private :: set_initial_conditions
    !> Load chemical species concentrations from configuration
    procedure, private :: set_chemical_species
    !> Load environmental conditions from configuration
    procedure, private :: set_environmental_conditions
    !> Load photolysis rate constants from configuration
    procedure, private :: set_photolysis_rate_constants
  end type initial_conditions_t

  !> Constructor for initial_conditions_t objects
  interface initial_conditions_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for initial_conditions_t objects
  function constructor( config, domain ) result( new_obj )

    use musica_config,                 only : config_t
    use musica_domain,                 only : domain_t

    !> New initial conditions
    type(initial_conditions_t), pointer :: new_obj
    !> Initial conditions configuration
    type(config_t), intent(inout) :: config
    !> Model domain
    class(domain_t), intent(inout) :: domain

    allocate( new_obj )
    call new_obj%set_initial_conditions( config, domain,                      &
                                         new_obj%variable_names_ )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates a model state loaded with the initial conditions
  function get_state( this, domain ) result( new_state )

    use musica_assert,                 only : assert
    use musica_domain,                 only : domain_t
    use musica_domain_state,           only : domain_state_t
    use musica_domain_target_cells,    only : domain_target_cells_t
    use musica_domain_iterator,        only : domain_iterator_t

    !> New model state loaded with initial conditions
    class(domain_state_t), pointer :: new_state
    !> Initial conditions
    class(initial_conditions_t), intent(in) :: this
    !> Model domain
    class(domain_t), intent(in) :: domain

    class(domain_iterator_t), pointer :: iter
    integer(kind=musica_ik) :: i_file, i_value
    type(domain_target_cells_t) :: all_cells

    new_state => domain%new_state( )
    do i_file = 1, size( this%files_ )
      call assert( 765999643, associated( this%files_( i_file )%val_ ) )
      call this%files_( i_file )%val_%update_state( domain, new_state )
    end do
    iter => domain%iterator( all_cells )
    do i_value = 1, size( this%initial_values_ )
      associate( val => this%initial_values_( i_value ) )
      call iter%reset( )
      do while( iter%next( ) )
        call new_state%update( iter, val%mutator_, val%value_ )
      end do
      end associate
    end do
    deallocate( iter )

  end function get_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update a model state for variables "tethered" to initial conditions
  subroutine update_state( this, domain, state )

    use musica_assert,                 only : assert
    use musica_domain,                 only : domain_t
    use musica_domain_state,           only : domain_state_t

    !> Initial conditions
    class(initial_conditions_t), intent(in) :: this
    !> Model domain
    class(domain_t), intent(in) :: domain
    !> Model state
    class(domain_state_t), intent(inout) :: state

    integer(kind=musica_ik) :: i_file

    do i_file = 1, size( this%files_ )
      call assert( 930075706, associated( this%files_( i_file )%val_ ) )
      call this%files_( i_file )%val_%update_state( domain, state,            &
                                                    tethered_only = .true. )
    end do

  end subroutine update_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Preprocess the initial conditions for repeat model runs
  subroutine preprocess_input( this, config, domain, output_path )

    use musica_config,                 only : config_t
    use musica_domain,                 only : domain_t
    use musica_domain_state,           only : domain_state_t
    use musica_input_output_processor, only : input_output_processor_t

    !> Initial conditions
    class(initial_conditions_t), intent(in) :: this
    !> Generated initial conditions configuration
    type(config_t), intent(out) :: config
    !> Model domain
    class(domain_t), intent(inout) :: domain
    !> Output path to save data files to
    character(len=*), intent(in) :: output_path

    character(len=*), parameter :: my_name = "Initial conditions preprocessor"
    integer(kind=musica_ik) :: i_file, i_var
    type(config_t) :: config_file_opts, file_opts, temp_config
    type(string_t) :: units
    class(domain_state_t), pointer :: state
    class(input_output_processor_t), pointer :: init_cond_file

    write(*,*) "Saving initial conditions..."

    if( size( this%variable_names_ ) .eq. 0 ) return
    write(*,*) "  - saving file '"//output_path//"initial.nc'"
    call config_file_opts%empty( )
    do i_file = 1, size( this%files_ )
      call this%files_( i_file )%val_%preprocess_input( temp_config,          &
                                                        "prognosed" )
      call config_file_opts%merge_in( temp_config, my_name )
    end do
    call config%add( "initial.nc", config_file_opts, my_name )
    call file_opts%empty( )
    call file_opts%add( "intent", "output", my_name )
    call file_opts%add( "type", "netcdf", my_name )
    call file_opts%add( "file name", output_path//"initial.nc", my_name )
    init_cond_file => input_output_processor_t( file_opts )
    do i_var = 1, size( this%variable_names_ )
      associate( var_name => this%variable_names_( i_var ) )
      units = domain%units( var_name%to_char( ) )
      call init_cond_file%register_output_variable( domain,                   & !- model domain
                                                    var_name%to_char( ),      & !- variable name
                                                    units%to_char( ) )          !- units
      end associate
    end do
    state => this%get_state( domain )
    call init_cond_file%output( 0.0_musica_dk, domain, state )
    deallocate( state          )
    deallocate( init_cond_file )

    write(*,*) "...done!"

  end subroutine preprocess_input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalizes an initial_conditions_t object
  elemental subroutine initial_value_finalize( this )

    !> Initial conditions
    type(initial_value_t), intent(inout) :: this

    integer(kind=musica_ik) :: i_elem

    if( associated( this%mutator_ ) ) then
      deallocate( this%mutator_ )
      this%mutator_ => null( )
    end if

  end subroutine initial_value_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set the initial conditions in a domain state
  !!
  !! Conditions specified explicitly in the configuration data take
  !! precedence over those read from input files.
  !!
  subroutine set_initial_conditions( this, config, domain, variable_names )

    use musica_assert,                 only : assert
    use musica_config,                 only : config_t
    use musica_domain,                 only : domain_t
    use musica_input_output_processor, only : input_output_processor_t
    use musica_iterator,               only : iterator_t

    !> Initial conditions
    class(initial_conditions_t), intent(inout) :: this
    !> Initial condition configuration data
    type(config_t), intent(inout) :: config
    !> Model domain data
    class(domain_t), intent(inout) :: domain
    !> MUSICA variables for which initial conditions were set
    type(string_t), allocatable, intent(out) :: variable_names(:)

    character(len=*), parameter :: my_name = 'initial conditions'
    logical :: found
    type(config_t) :: subset, input_file_config, photolysis_config
    class(iterator_t), pointer :: iter
    type(string_t) :: temp_str
    type(string_t), allocatable :: str_array(:), musica_var_names(:)
    integer(kind=musica_ik) :: i_file, n_files

    allocate( this%initial_values_( 0 ) )
    allocate( variable_names( 0 ) )

    ! Load conditions from input files
    call config%get( "initial conditions", subset, my_name, found = found )
    if( found ) then
      iter => subset%get_iterator( )
      n_files = 0
      do while( iter%next( ) )
        n_files = n_files + 1
      end do
      allocate( this%files_( n_files ) )
      i_file = 0
      call iter%reset( )
      do while( iter%next( ) )
        i_file = i_file + 1
        call subset%get( iter, input_file_config, my_name )
        temp_str = subset%key( iter )
        str_array = temp_str%split( "." )
        temp_str = str_array( size( str_array ) )%to_lower( )
        call input_file_config%add( "type", temp_str, my_name )
        call input_file_config%add( "intent", "input", my_name )
        call input_file_config%add( "default variability", "prognosed",       &
                                    my_name )
        call input_file_config%add( "file name", subset%key( iter ), my_name )
        this%files_( i_file )%val_ =>                                         &
            input_output_processor_t( input_file_config, domain )
        call assert( 893945737, associated( this%files_( i_file )%val_ ) )
        musica_var_names = this%files_( i_file )%val_%musica_variable_names( )
        variable_names = [ variable_names, musica_var_names ]
      end do
      call assert( 274515453, i_file .eq. n_files )
      deallocate( iter )
    else
      allocate( this%files_( 0 ) )
    end if

    ! set all domain cell chemical species concentrations to specified
    ! values
    call config%get( "chemical species", subset, my_name, found = found )
    if( found ) then
      call this%set_chemical_species( subset, domain, variable_names )
    end if

    ! set all domain cell environmental conditions to specified values
    call config%get( "environmental conditions", subset, my_name,             &
                     found = found )
    if( found ) then
      call this%set_environmental_conditions( subset, domain, variable_names )
    end if

    ! set photolysis rate constants
    call config%get( "photolysis", photolysis_config, my_name, found = found )
    if( found ) then
      call this%set_photolysis_rate_constants( photolysis_config, domain,     &
                                          variable_names )
    end if

  end subroutine set_initial_conditions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set initial species concentrations for all domain cells
  subroutine set_chemical_species( this, config, domain, variable_names )

    use musica_assert,                 only : assert
    use musica_config,                 only : config_t
    use musica_data_type,              only : kDouble
    use musica_domain,                 only : domain_t
    use musica_domain_target_cells,    only : domain_target_cells_t
    use musica_iterator,               only : iterator_t
    use musica_property,               only : property_t

    !> Initial conditions
    class(initial_conditions_t), intent(inout) :: this
    !> Configuration data
    type(config_t), intent(inout) :: config
    !> Model domain data
    class(domain_t), intent(inout) :: domain
    !> MUSICA variables for which initial conditions were set
    type(string_t), allocatable, intent(inout) :: variable_names(:)

    character(len=*), parameter :: my_name = 'initial species concentrations'
    type(config_t) :: subset
    type(string_t) :: species_name
    real(kind=musica_dk) :: conc
    class(iterator_t), pointer :: species_iter
    integer(kind=musica_ik) :: i_value, n_values
    class(property_t), pointer :: prop
    type(domain_target_cells_t) :: all_cells

    call assert( 142287599, allocated( variable_names ) )
    species_iter => config%get_iterator( )
    n_values = 0
    do while( species_iter%next( ) )
      n_values = n_values + 1
    end do
    call increase_array_size( this%initial_values_, n_values )
    call species_iter%reset( )
    i_value = size( this%initial_values_ ) - n_values
    do while( species_iter%next( ) )
      i_value = i_value + 1
      species_name = "chemical_species%"//config%key( species_iter )
      call config%get( species_iter, subset, my_name )
      call subset%get( "initial value", "mol m-3", conc, my_name )
      prop => property_t( my_name,                                            &
                          name = species_name%to_char( ),                     & !- state variable name
                          units = "mol m-3",                                  & !- MUSICA units
                          data_type = kDouble,                                & !- data type
                          applies_to = all_cells )                              !- target domain
      this%initial_values_( i_value )%mutator_ => domain%mutator( prop )
      deallocate( prop )
      this%initial_values_( i_value )%value_ = conc
      variable_names = [ variable_names, species_name ]
    end do
    call assert( 332266831, i_value .eq. size( this%initial_values_ ) )
    do i_value = 1, size( this%initial_values_ )
      call assert( 204426349,                                                 &
                   associated( this%initial_values_( i_value )%mutator_ ) )
    end do

    ! clean up
    deallocate( species_iter )

  end subroutine set_chemical_species

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set environmental conditions for all domain cells
  subroutine set_environmental_conditions( this, config, domain,              &
      variable_names )

    use musica_assert,                 only : assert
    use musica_config,                 only : config_t
    use musica_data_type,              only : kDouble
    use musica_domain,                 only : domain_t
    use musica_domain_target_cells,    only : domain_target_cells_t
    use musica_iterator,               only : iterator_t
    use musica_property,               only : property_t

    !> Initial conditions
    class(initial_conditions_t), intent(inout) :: this
    !> Configuration data
    type(config_t), intent(inout) :: config
    !> Model domain data
    class(domain_t), intent(inout) :: domain
    !> MUSICA variables for which initial conditions were set
    type(string_t), allocatable, intent(inout) :: variable_names(:)

    character(len=*), parameter :: my_name =                                  &
      'initial environmental conditions'
    type(config_t) :: subset
    type(string_t) :: property_name, units
    real(musica_dk) :: property_value
    class(iterator_t), pointer :: property_iter
    integer(kind=musica_ik) :: i_value, n_values
    type(domain_target_cells_t) :: all_cells
    class(property_t), pointer :: prop

    call assert( 258704135, allocated( variable_names ) )
    property_iter => config%get_iterator( )
    n_values = 0
    do while( property_iter%next( ) )
      n_values = n_values + 1
    end do
    call increase_array_size( this%initial_values_, n_values )
    call property_iter%reset( )
    i_value = size( this%initial_values_ ) - n_values
    do while( property_iter%next( ) )
      i_value = i_value + 1
      property_name = config%key( property_iter )
      units         = domain%units( property_name%to_char( ) )
      call config%get( property_iter, subset, my_name )
      call subset%get( "initial value",  units%to_char( ), property_value,    &
                      my_name )
      prop => property_t( my_name,                                            &
                          name = property_name%to_char( ),                    & !- state variable name
                          units = units%to_char( ),                           & !- MUSICA units
                          data_type = kDouble,                                & !- data type
                          applies_to = all_cells )                              !- target domain
      this%initial_values_( i_value )%mutator_ => domain%mutator( prop )
      this%initial_values_( i_value )%value_ = property_value
      variable_names = [ variable_names, property_name ]
      deallocate( prop )
    end do
    call assert( 745759824, i_value .eq. size( this%initial_values_ ) )
    do i_value = 1, size( this%initial_values_ )
      call assert( 639812342,                                                 &
                   associated( this%initial_values_( i_value )%mutator_ ) )
    end do

    ! clean up
    deallocate( property_iter )

  end subroutine set_environmental_conditions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set photolysis rates for all domain cells
  subroutine set_photolysis_rate_constants( this, config, domain,             &
      variable_names )

    use musica_assert,                 only : assert
    use musica_config,                 only : config_t
    use musica_data_type,              only : kDouble
    use musica_domain,                 only : domain_t
    use musica_domain_target_cells,    only : domain_target_cells_t
    use musica_iterator,               only : iterator_t
    use musica_property,               only : property_t

    !> Initial conditions
    class(initial_conditions_t), intent(inout) :: this
    !> Configuration data
    type(config_t), intent(inout) :: config
    !> Model domain data
    class(domain_t), intent(inout) :: domain
    !> MUSICA variables for which initial conditions were set
    type(string_t), allocatable, intent(inout) :: variable_names(:)

    character(len=*), parameter :: my_name =                                  &
      'initial photolysis rate constants'
    type(config_t) :: subset
    type(string_t) :: photo_name, units
    real(musica_dk) :: rate_constant
    class(iterator_t), pointer :: photo_iter
    integer(kind=musica_ik) :: i_value, n_values
    class(property_t), pointer :: prop
    type(domain_target_cells_t) :: all_cells

    call assert( 122874379, allocated( variable_names ) )
    photo_iter => config%get_iterator( )
    n_values = 0
    do while( photo_iter%next( ) )
      n_values = n_values + 1
    end do
    call increase_array_size( this%initial_values_, n_values )
    call photo_iter%reset( )
    i_value = size( this%initial_values_ ) - n_values
    do while( photo_iter%next( ) )
      i_value = i_value + 1
      photo_name = "photolysis_rate_constants%"//config%key( photo_iter )
      units      = domain%units( photo_name%to_char( ) )
      call config%get( photo_iter, subset, my_name )
      call subset%get( "initial value", units%to_char( ), rate_constant,      &
                       my_name )
      prop => property_t( my_name,                                            &
                          name = photo_name%to_char( ),                       & !- state variable name
                          units = units%to_char( ),                           & !- MUSICA units
                          data_type = kDouble,                                & !- data type
                          applies_to = all_cells )                              !- target domain
      this%initial_values_( i_value )%mutator_ => domain%mutator( prop )
      this%initial_values_( i_value )%value_ = rate_constant
      variable_names = [ variable_names, photo_name ]
      deallocate( prop )
    end do
    call assert( 655047356, i_value .eq. size( this%initial_values_ ) )
    do i_value = 1, size( this%initial_values_ )
      call assert( 246924282,                                                 &
                   associated( this%initial_values_( i_value )%mutator_ ) )
    end do

    ! clean up
    deallocate( photo_iter )

  end subroutine set_photolysis_rate_constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets the domain state variable value
  subroutine initial_value_set_state( this, iterator, state )

    use musica_assert,                 only : assert
    use musica_domain_state,           only : domain_state_t
    use musica_domain_iterator,        only : domain_iterator_t

    !> Initial value
    class(initial_value_t), intent(in) :: this
    !> Domain state iterator
    class(domain_iterator_t), intent(in) :: iterator
    !> Model domain state
    class(domain_state_t), intent(inout) :: state

    call assert( 839805468, associated( this%mutator_ ) )
    call state%update( iterator, this%mutator_, this%value_ )

  end subroutine initial_value_set_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Increase the size of an initial_value_t array by a specified number of
  !! elements
  subroutine increase_array_size( array, number_of_elements )

    use musica_assert,                 only : assert

    !> Initial value array
    type(initial_value_t), allocatable, intent(inout) :: array(:)
    !> Number of elements to add to array (added to the end of the array)
    integer(kind=musica_ik), intent(in) :: number_of_elements

    integer(kind=musica_ik) :: i_elem
    type(initial_value_t), allocatable :: temp_array(:)

    if( number_of_elements .eq. 0 ) return
    call assert( 745020406, number_of_elements .gt. 0 )
    if( .not. allocated( array ) ) then
      allocate( array( number_of_elements ) )
      return
    end if
    allocate( temp_array, source = array )
    do i_elem = 1, size( array )
      array( i_elem )%mutator_ => null( )
    end do
    deallocate( array )
    allocate( array( size( temp_array ) + number_of_elements ) )
    array( :size( temp_array ) ) = temp_array(:)
    do i_elem = 1, size( temp_array )
      temp_array( i_elem )%mutator_ => null( )
    end do
    deallocate( temp_array )
    do i_elem = 1, size( array ) - number_of_elements
      call assert( 660704398, associated( array( i_elem )%mutator_ ) )
    end do
    do i_elem = size( array ) - number_of_elements + 1, size( array )
      call assert( 437973242, .not. associated( array( i_elem )%mutator_ ) )
    end do

  end subroutine increase_array_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_initial_conditions
