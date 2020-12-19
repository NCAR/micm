! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The micm_core module

!> The core_t type and related functions
module micm_core

  use micm_kinetics,                   only : kinetics_t
  use micm_ODE_solver,                 only : ODE_solver_t
  use musica_component,                only : component_t
  use musica_constants,                only : musica_dk, musica_ik
  use musica_domain_state_mutator,     only : domain_state_mutator_ptr
  use musica_domain_state_accessor,    only : domain_state_accessor_ptr,      &
                                              domain_state_accessor_t
  use musica_string,                   only : string_t


  implicit none
  private

  public :: core_t

  !> MICM core
  !!
  !! Top-level chemistry type. A core initializes the chemical scheme,
  !! solves for chemistry over given time steps, and finalizes chemistry
  !! objects.
  !!
  !! \todo ensure that MICM core is thread-safe
  type, extends(component_t) :: core_t
    private
    !> ODE solver
    class(ODE_solver_t), pointer :: ODE_solver_ => null( )
    !> Kinetics calculator
    class(kinetics_t), pointer :: kinetics_ => null( )
    !> Mutators for chemical species
    class(domain_state_mutator_ptr), pointer ::                               &
        species_mutators_(:) => null( )
    !> Accessors for chemical species
    class(domain_state_accessor_ptr), pointer ::                              &
        species_accessors_(:) => null( )
    !> Mutators for reaction rates
    class(domain_state_mutator_ptr), pointer ::                               &
        rate_mutators_(:) => null( )
    !> Accessors for photolysis rate constants
    class(domain_state_accessor_ptr), pointer ::                              &
        photolysis_rate_constant_accessors_(:) => null( )
    !> Environmental property accessors
    !! \todo move MICM environmental accessors to micm_environment module
    !! @{
    class(domain_state_accessor_t), pointer :: temperature__K_ => null ()
    class(domain_state_accessor_t), pointer :: pressure__Pa_ => null( )
    class(domain_state_accessor_t), pointer ::                                &
        number_density_air__mol_m3_ => null( )
    !> @}

    !> Working number density arrray [molec cm-3]
    real(kind=musica_dk), allocatable :: number_densities__molec_cm3_(:)
    !> Working reaction rate array [molec cm-3 s-1]
    real(kind=musica_dk), allocatable :: reaction_rates__molec_cm3_s_(:)
    !> Working photolysis rate constant array [s-1]
    real(kind=musica_dk), allocatable :: photolysis_rate_constants__s_(:)
    !> Flag indicating whether to output reaction rates
    logical :: output_reaction_rates_ = .false.
    !> Flag indicating whether to output photolysis rate constants
    logical :: output_photolysis_rate_constants_ = .false.
    !> Flag indicating whether to solve chemistry during the simulation
    logical :: solve_ = .true.
  contains
    !> Returns the name of the component
    procedure :: name => component_name
    !> Returns a description of the component purpose
    procedure :: description
    !> Solve chemistry for one or more grid cells
    procedure :: advance_state
    !> Preprocess chemistry input data
    procedure :: preprocess_input
    !> Set the initial conditions for the current time step
    procedure, private :: time_step_initialize
    !> Finalize the chemistry core
    final :: finalize
  end type core_t

  !> Constructor
  interface core_t
    module procedure constructor
  end interface core_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> MICM Core constructor
  !!
  !! Sets up chemistry objects for solving
  function constructor( config, domain, output ) result( new_obj )

    use micm_ODE_solver_factory,       only : ODE_solver_builder
    use musica_assert,                 only : assert
    use musica_config,                 only : config_t
    use musica_data_type,              only : kDouble
    use musica_domain,                 only : domain_t
    use musica_domain_target_cells,    only : domain_target_cells_t
    use musica_input_output_processor, only : input_output_processor_t
    use musica_property,               only : property_t
    use musica_property_set,           only : property_set_t
    use musica_string,                 only : string_t

    !> New MICM Core
    type(core_t), pointer :: new_obj
    !> Chemistry configuration data
    type(config_t), intent(inout) :: config
    !> Model domain
    class(domain_t), intent(inout) :: domain
    !> Output file
    class(input_output_processor_t), intent(inout) :: output

    character(len=*), parameter :: my_name = 'MICM chemistry constructor'
    integer :: i_spec, i_rxn
    type(string_t), allocatable :: species_names(:), reaction_names(:),       &
                                   photo_names(:)
    type(config_t) :: solver_opts, outputs, output_opts
    real(kind=musica_dk) :: chemistry_time_step__s
    type(property_set_t), pointer :: prop_set
    type(property_t), pointer :: prop
    logical :: found
    type(domain_target_cells_t) :: all_cells

    allocate( new_obj )

    ! Check whether to solve chemistry
    call config%get( "solve", new_obj%solve_, my_name, default = .true. )

    ! Set up the kinetics calculator
    new_obj%kinetics_ => kinetics_t( )
    call new_obj%kinetics_%species_names( species_names )
    call new_obj%kinetics_%reaction_names( reaction_names )
    call new_obj%kinetics_%photolysis_reaction_names( photo_names )

    ! Set up the solver
    call config%get( "solver", solver_opts, my_name )
    call solver_opts%add( "number of variables", size( species_names ),       &
                          my_name )
    new_obj%ODE_solver_ => ODE_solver_builder( solver_opts )

    ! Register state variables for the chemical species concentrations
    prop_set => property_set_t( )
    do i_spec = 1, size( species_names )
      prop => property_t( my_name,                                            &
                          name = species_names( i_spec )%to_char( ),          &
                          units = "mol m-3",                                  &
                          applies_to = all_cells,                             &
                          data_type = kDouble,                                &
                          default_value = 0.0_musica_dk )
      call prop_set%add( prop )
      deallocate( prop )
    end do
    call domain%register( "chemical_species", prop_set )
    deallocate( prop_set )
    new_obj%species_mutators_ => domain%mutator_set(   "chemical_species",    & !- property set name
                                                       "mol m-3",             & !- units
                                                       kDouble,               & !- data type
                                                       all_cells,             & !- variable domain
                                                       my_name )
    new_obj%species_accessors_ => domain%accessor_set( "chemical_species",    & !- property set name
                                                       "mol m-3",             & !- units
                                                       kDouble,               & !- data type
                                                       all_cells,             & !- variable domain
                                                       my_name )
    call assert( 862216278, size( new_obj%species_mutators_  ) .eq.           &
                            size( species_names ) )
    call assert( 865575051, size( new_obj%species_accessors_ ) .eq.           &
                            size( species_names ) )

    ! Register state variables for reaction rates
    prop_set => property_set_t( )
    do i_rxn = 1, size( reaction_names )
      prop => property_t( my_name,                                            &
                          name = reaction_names( i_rxn )%to_char( ),          &
                          units = "mol m-3 s-1",                              &
                          applies_to = all_cells,                             &
                          data_type = kDouble,                                &
                          default_value = 0.0_musica_dk )
      call prop_set%add( prop )
      deallocate( prop )
    end do
    call domain%register( "reaction_rates", prop_set )
    deallocate( prop_set )
    new_obj%rate_mutators_ => domain%mutator_set( "reaction_rates",           & !- property set name
                                                  "mol m-3 s-1",              & !- units
                                                  kDouble,                    & !- data type
                                                  all_cells,                  & !- variable domain
                                                  my_name )
    call assert( 746539160, size( new_obj%rate_mutators_  ) .eq.              &
                            size( reaction_names ) )

    ! Register an array of photolysis rate constants so that other model
    ! components can set photolysis rates
    prop_set => property_set_t( )
    do i_rxn = 1, size( photo_names )
      prop => property_t( my_name,                                            &
                          name = photo_names( i_rxn )%to_char( ),             &
                          units = "s-1",                                      &
                          applies_to = all_cells,                             &
                          data_type = kDouble,                                &
                          default_value = 0.0_musica_dk )
      call prop_set%add( prop )
      deallocate( prop )
    end do
    call domain%register( "photolysis_rate_constants", prop_set )
    deallocate( prop_set )
    new_obj%photolysis_rate_constant_accessors_ =>                            &
        domain%accessor_set( "photolysis_rate_constants",                     & !- property set name
                             "s-1",                                           & !- units
                             kDouble,                                         & !- data type
                             all_cells,                                       & !- variable domain
                             my_name )
    call assert( 295812541,                                                   &
                 size( new_obj%photolysis_rate_constant_accessors_ ) .eq.     &
                 size( photo_names ) )

    ! Register accessors for environmental properties
    prop => property_t( my_name, name = "temperature", units = "K",           &
                        applies_to = all_cells, data_type = kDouble )
    new_obj%temperature__K_ => domain%accessor( prop )
    deallocate( prop )
    prop => property_t( my_name, name = "pressure", units = "Pa",             &
                        applies_to = all_cells, data_type = kDouble )
    new_obj%pressure__Pa_ => domain%accessor( prop )
    deallocate( prop )
    prop => property_t( my_name, name = "number density air",                 &
                        units = "mol m-3", applies_to = all_cells,            &
                        data_type = kDouble )
    new_obj%number_density_air__mol_m3_ => domain%accessor( prop )
    deallocate( prop )

    ! Register the chemical species concentrations for output
    do i_spec = 1, size( species_names )
      call output%register_output_variable(                                   &
                            domain,                                           &
                            "chemical_species%"//                             & !- variable full name
                                species_names( i_spec )%to_char( ),           &
                            "mol m-3",                                        & !- units
                            "CONC."//species_names( i_spec )%to_char( ) )       !- output name
    end do

    ! set up optional output
    call config%get( "output", outputs, my_name, found = found )
    if( found ) then

      ! reaction rate output
      call outputs%get( "reaction rates", output_opts, my_name, found = found )
      if( found ) then
        new_obj%output_reaction_rates_ = .true.
        do i_rxn = 1, size( reaction_names )
          call output%register_output_variable(                               &
                            domain,                                           &
                            "reaction_rates%"//                               & !- variable full name
                                reaction_names( i_rxn )%to_char( ),           &
                            "mol m-3 s-1",                                    & !- units
                            "RATE."//reaction_names( i_rxn )%to_char( ) )       !- output name
        end do
      end if

      ! photolysis rate constant output
      call outputs%get( "photolysis rate constants", output_opts, my_name,    &
                        found = found )
      if( found ) then
        new_obj%output_photolysis_rate_constants_ = .true.
        do i_rxn = 1, size( photo_names )
          call output%register_output_variable(                               &
                                domain,                                       &
                                "photolysis_rate_constants%"//                & !- variable full name
                                    photo_names( i_rxn )%to_char( ),          &
                                "s-1",                                        & !- units
                                "PHOTO."//photo_names( i_rxn )%to_char( ) )
        end do
      end if
    end if

    ! Set up arrays for use during solving
    allocate( new_obj%number_densities__molec_cm3_( size( species_names ) ) )
    allocate( new_obj%reaction_rates__molec_cm3_s_( size( reaction_names ) ) )
    allocate( new_obj%photolysis_rate_constants__s_( size( photo_names ) ) )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Model component name
  type(string_t) function component_name( this )

    !> CAMP interface
    class(core_t), intent(in) :: this

    component_name = "MICM: Model-Independent Chemical Mechanisms"

  end function component_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Model component description
  type(string_t) function description( this )

    !> CAMP interface
    class(core_t), intent(in) :: this

    description = "Gas-phase chemistry solver"

  end function description

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Solve chemistry for a given number of grid cells and time step
  subroutine advance_state( this, domain_state, domain_element,               &
      current_time__s, time_step__s )

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : kAvagadro
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_domain_state,           only : domain_state_t
    use musica_string,                 only : to_char

    !> MICM chemistry
    class(core_t), intent(inout) :: this
    !> Domain state
    class(domain_state_t), intent(inout) :: domain_state
    !> Grid cell to solve
    class(domain_iterator_t), intent(in) :: domain_element
    !> Current simulation time [s]
    real(kind=musica_dk), intent(in) :: current_time__s
    !> Chemistry time step [s]
    real(kind=musica_dk), intent(in) :: time_step__s

    integer(kind=musica_ik) :: i_spec, i_rxn, error_flag

    if( .not. this%solve_ ) return

    ! Set the initial conditions for the time step
    call this%time_step_initialize( domain_state, domain_element )

    ! solve the chemistry for this time step
    call this%ODE_solver_%solve( TStart = 0.0_musica_dk,                      &
                                 TEnd   = time_step__s,                       &
                                 y      = this%number_densities__molec_cm3_,  &
                                 theKinetics = this%kinetics_,                &
                                 IErr   = error_flag )

    call assert_msg( 534725427, error_flag .eq. 0,                            &
                     "Chemistry solver failed with code "//                   &
                     to_char( error_flag ) )

    ! update the species concentrations [mol m-3]
    do i_spec = 1, size( this%species_mutators_ )
      call domain_state%update( domain_element,                               &
                                this%species_mutators_( i_spec )%val_,        &
                                this%number_densities__molec_cm3_( i_spec ) / &
                                    kAvagadro * 1.0d6 )
    end do

    ! save the reaction rates [mol m-3 s-1]
    this%reaction_rates__molec_cm3_s_ =                                       &
        this%kinetics_%reaction_rates( this%number_densities__molec_cm3_ )
    do i_rxn = 1, size( this%reaction_rates__molec_cm3_s_ )
      call domain_state%update( domain_element,                               &
                                this%rate_mutators_( i_rxn )%val_,            &
                                this%reaction_rates__molec_cm3_s_( i_rxn ) /  &
                                    kAvagadro * 1.0d6 )
    end do

  end subroutine advance_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Preprocess chemistry input data
  subroutine preprocess_input( this, config, output_path )

    use musica_assert,                 only : assert
    use musica_config,                 only : config_t

    !> MICM chemistry
    class(core_t), intent(inout) :: this
    !> Chemistry configuration
    type(config_t), intent(out) :: config
    !> Folder to save input data to
    character(len=*), intent(in) :: output_path

    character(len=*), parameter :: my_name = "MICM input preprocessor"
    type(config_t) :: solver, output, empty_config

    write(*,*) "Saving chemistry configuration..."

    call empty_config%empty( )
    call config%add( "type", "MICM", my_name )
    call config%add( "solve", this%solve_, my_name )
    call assert( 418280390, associated( this%ODE_solver_ ) )
    call this%ODE_solver_%preprocess_input( solver, output_path )
    call config%add( "solver", solver, my_name )
    call output%empty( )
    if( this%output_reaction_rates_ ) then
      call output%add( "reaction rates", empty_config, my_name )
    end if
    if( this%output_photolysis_rate_constants_ ) then
      call output%add( "photolysis rate constants", empty_config, my_name )
    end if
    if( output%number_of_children( ) .gt. 0 ) then
      call config%add( "output", output, my_name )
    end if

  end subroutine preprocess_input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set initial conditions for current time step
  !!
  !! Sets environmental conditions and calculates rate constants
  subroutine time_step_initialize( this, domain_state, cell )

    use micm_environment,              only : environment_t
    use musica_constants,              only : kAvagadro
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_domain_state,           only : domain_state_t

    !> MICM chemistry
    class(core_t), intent(inout) :: this
    !> Domain state
    class(domain_state_t), intent(inout) :: domain_state
    !> Grid cell to solve
    class(domain_iterator_t), intent(in) :: cell

    integer :: i_spec, i_photo
    type(environment_t) :: env

    ! get the current environmental conditions
    call domain_state%get( cell, this%temperature__K_, env%temperature )
    call domain_state%get( cell, this%pressure__Pa_,   env%pressure )
    call domain_state%get( cell, this%number_density_air__mol_m3_,            & ! currently in [mol m-3]
                           env%number_density_air )

    ! convert air density to non-standard units currently used in MICM
    ! [molec cm-3]
    env%number_density_air = env%number_density_air * kAvagadro * 1.0d-6

    ! get the current species concentrations [mol m-3]
    ! and convert to non-standard units currently used in MICM [molec cm-3]
    !> \todo update MICM to use [mol m-3] for number densities
    do i_spec = 1, size( this%species_accessors_ )
      call domain_state%get( cell,                                            &
                             this%species_accessors_( i_spec )%val_,          &
                             this%number_densities__molec_cm3_( i_spec ) )
      this%number_densities__molec_cm3_( i_spec ) =                           &
          this%number_densities__molec_cm3_( i_spec ) * kAvagadro * 1.0d-6
    end do

    ! update the photolysis rates
    do i_photo = 1, size( this%photolysis_rate_constant_accessors_ )
      call domain_state%get( cell,                                            &
                    this%photolysis_rate_constant_accessors_( i_photo )%val_, &
                    this%photolysis_rate_constants__s_( i_photo ) )
    end do

    ! update the kinetics for the current conditions
    call this%kinetics_%update( env, this%photolysis_rate_constants__s_ )

  end subroutine time_step_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the chemistry core
  subroutine finalize( this )

    !> MICM chemistry
    type(core_t), intent(inout) :: this

    integer :: i

    if( associated( this%ODE_solver_ ) )                                      &
      deallocate( this%ODE_solver_ )
    if( associated( this%kinetics_ ) )                                        &
      deallocate( this%kinetics_ )
    if( associated( this%species_mutators_ ) ) then
      do i = 1, size( this%species_mutators_ )
        if( associated( this%species_mutators_( i )%val_ ) ) then
          deallocate( this%species_mutators_( i )%val_ )
        end if
      end do
      deallocate( this%species_mutators_ )
    end if
    if( associated( this%species_accessors_ ) ) then
      do i = 1, size( this%species_accessors_ )
        if( associated( this%species_accessors_( i )%val_ ) ) then
          deallocate( this%species_accessors_( i )%val_ )
        end if
      end do
      deallocate( this%species_accessors_ )
    end if
    if( associated( this%rate_mutators_ ) ) then
      do i = 1, size( this%rate_mutators_ )
        if( associated( this%rate_mutators_( i )%val_ ) ) then
          deallocate( this%rate_mutators_( i )%val_ )
        end if
      end do
      deallocate( this%rate_mutators_ )
    end if
    if( associated( this%photolysis_rate_constant_accessors_ ) ) then
      do i = 1, size( this%photolysis_rate_constant_accessors_ )
        if( associated( this%photolysis_rate_constant_accessors_( i )%val_ ) )&
            then
          deallocate( this%photolysis_rate_constant_accessors_( i )%val_ )
        end if
      end do
      deallocate( this%photolysis_rate_constant_accessors_ )
    end if
    if( associated( this%temperature__K_ ) )                                  &
      deallocate( this%temperature__K_ )
    if( associated( this%pressure__Pa_ ) )                                    &
      deallocate( this%pressure__Pa_ )
    if( associated( this%number_density_air__mol_m3_ ) )                      &
      deallocate( this%number_density_air__mol_m3_ )

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_core
