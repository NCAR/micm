!> \file
!> The micm_core module

!> The core_t type and related functions
module micm_core

  use musica_domain,                   only : domain_state_mutator_ptr,       &
                                              domain_state_accessor_ptr
  use musica_string,                   only : string_t

  implicit none
  private

  public :: core_t

  !> MICM core
  !!
  !! Top-level chemistry object. The core initializes the chemical scheme,
  !! solves for chemistry over given time steps, and finalizes chemistry
  !! objects.
  type :: core_t
    private
    !> Chemical species names
    type(string_t), allocatable :: species_names_(:)
    !> Mutators for chemical species
    class(domain_state_mutator_ptr), pointer ::                               &
        species_mutators_(:) => null( )
    !> Accessors for chemical species
    class(domain_state_accessor_ptr), pointer ::                              &
        species_accessors_(:) => null( )
  contains
    !> Get the name of each chemical species in the chemistry state array
    procedure :: species_names
    !> Solve chemistry for one or more grid cells
    procedure :: solve
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

    use musica_assert,                 only : assert
    use musica_config,                 only : config_t
    use musica_domain,                 only : domain_t
    use musica_output,                 only : output_t
    use musica_string,                 only : string_t

    !> New MICM Core
    type(core_t), pointer :: new_obj
    !> Domain
    class(domain_t), intent(inout) :: domain
    !> Chemistry configuration data
    class(config_t), intent(inout) :: config
    !> Output stream
    class(output_t), intent(inout) :: output

    character(len=*), parameter :: my_name = 'MICM chemistry constructor'
    integer :: i_spec
    type(string_t), allocatable :: accessor_names(:)

    allocate( new_obj )

    ! read species from molec.json
    allocate( new_obj%species_names_( 5 ) )
    new_obj%species_names_(1) = "NO2"
    new_obj%species_names_(2) = "ISOP"
    new_obj%species_names_(3) = "O3"
    new_obj%species_names_(4) = "HCHO"
    new_obj%species_names_(5) = "NO"

    new_obj%species_mutators_ =>                                              &
      domain%register_cell_state_variable_set( "chemical_species",            & !<- variable set name
                                               "mol m-3",                     & !<- units
                                               0.0d0,                         & !<- default value
                                               new_obj%species_names_,        & !<- variable element names
                                               my_name )
    new_obj%species_accessors_ =>                                             &
      domain%cell_state_set_accessor( "chemical_species",                     & !<- variable set name
                                      "mol m-3",                              & !<- units
                                      accessor_names,                         & !<- variable element names
                                      my_name )

    call assert( 415788666, size( new_obj%species_names_ ) .eq.               &
                            size( accessor_names ) )
    do i_spec = 1, size( new_obj%species_names_ )
      call assert( 359403346, new_obj%species_names_( i_spec ) .eq.           &
                              accessor_names( i_spec ) )
      call output%register( domain,                                           &
                            "chemical_species%"//                             & !<- varible full name
                                new_obj%species_names_( i_spec )%to_char( ),  &
                            "mol m-3",                                        & !<- units
                            "CONC."//                                         & !<- output name
                                new_obj%species_names_( i_spec )%to_char( ) )
    end do

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the species names for each chemical species in the order expected for
  !! the chemistry state array
  function species_names( this )

    !> Species names
    type(string_t), allocatable :: species_names(:)
    !> MICM Core
    class(core_t), intent(in) :: this

    species_names = this%species_names_

  end function species_names

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Solve chemistry for a given number of grid cells and time step
  subroutine solve( this, domain_state, cell, time_step__s )

    use musica_constants,              only : musica_dk
    use musica_domain,                 only : domain_state_t,                 &
                                              domain_iterator_t

    !> MICM chemistry
    class(core_t), intent(inout) :: this
    !> Domain state
    class(domain_state_t), intent(inout) :: domain_state
    !> Grid cell to solve
    class(domain_iterator_t), intent(in) :: cell
    !> Chemistry time step [s]
    real(kind=musica_dk), intent(in) :: time_step__s

    write(*,*) "Solving chemistry! time step: ", time_step__s, " s"

  end subroutine solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the chemistry core
  subroutine finalize( this )

    !> MICM chemistry
    type(core_t), intent(inout) :: this

    if( associated( this%species_mutators_ ) )                                &
      deallocate( this%species_mutators_ )
    if( associated( this%species_accessors_ ) )                               &
      deallocate( this%species_accessors_ )

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_core
