!> \file
!> The micm_core module

!> The core_t type and related functions
module micm_core

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
  contains
    !> Get the name of each chemical species in the chemistry state array
    procedure :: species_names
    !> Solve chemistry for one or more grid cells
    procedure :: solve
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
  function constructor( config ) result( new_obj )

    use musica_config,                 only : config_t

    !> New MICM Core
    type(core_t) :: new_obj
    !> Chemistry configuration data
    class(config_t), intent(inout) :: config

    ! read species from molec.json
    allocate( new_obj%species_names_( 3 ) )
    new_obj%species_names_(1) = "A"
    new_obj%species_names_(2) = "B"
    new_obj%species_names_(3) = "C"

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
  subroutine solve( this, time_step__s )

    use musica_constants,              only : musica_dk

    !> MICM chemistry
    class(core_t), intent(inout) :: this
    !> Chemistry time step [s]
    real(kind=musica_dk), intent(in) :: time_step__s

    write(*,*) "Solving chemistry! time step: ", time_step__s, " s"

  end subroutine solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_core
