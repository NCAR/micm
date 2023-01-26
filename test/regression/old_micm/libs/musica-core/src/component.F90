! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_component module

!> The abstract component_t type and related functions
module musica_component

  implicit none
  private

  public :: component_t, component_ptr

  !> Model component
  !!
  !! Model components calculate diagnostics for a given model state and/or
  !! advance the model state for a given timestep.
  !!
  !! \todo add full description and example usage for component_t
  !!
  type, abstract :: component_t
  contains
    !> Returns the name of the component
    procedure(component_name), deferred :: name
    !> Returns a description of the component purpose
    procedure(description), deferred :: description
    !> Advance the model state for a given timestep
    procedure(advance_state), deferred :: advance_state
    !> Save the component configuration for future simulations
    procedure(preprocess_input), deferred :: preprocess_input
  end type component_t

  !> Unique pointer for component_t objects
  type :: component_ptr
    class(component_t), pointer :: val_
  contains
    !> Finalize the pointer
    final :: component_ptr_finalize
  end type component_ptr

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the name of the component
  type(string_t) function component_name( this )
    use musica_string,                 only : string_t
    import component_t
    !> Model component
    class(component_t), intent(in) :: this
  end function component_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns a description of the component purpose
  type(string_t) function description( this )
    use musica_string,                 only : string_t
    import component_t
    !> Model component
    class(component_t), intent(in) :: this
  end function description

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Advance the model state for a given timestep
  subroutine advance_state( this, domain_state, domain_element,               &
      current_time__s, time_step__s )
    use musica_constants,              only : musica_dk
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_domain_state,           only : domain_state_t
    import component_t
    !> Model component
    class(component_t), intent(inout) :: this
    !> Domain state
    class(domain_state_t), intent(inout) :: domain_state
    !> Domain element to advance state for
    class(domain_iterator_t), intent(in) :: domain_element
    !> Current simulation time [s]
    real(kind=musica_dk), intent(in) :: current_time__s
    !> Time step to advance state by [s]
    real(kind=musica_dk), intent(in) :: time_step__s
  end subroutine advance_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Save the component configuration for future simulations
  subroutine preprocess_input( this, config, output_path )
    use musica_config,                 only : config_t
    import component_t
    !> Model component
    class(component_t), intent(inout) :: this
    !> Model component configuration
    type(config_t), intent(out) :: config
    !> Folder to save input data to
    character(len=*), intent(in) :: output_path
  end subroutine preprocess_input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the unique pointer
  elemental subroutine component_ptr_finalize( this )

    !> Component pointer
    type(component_ptr), intent(inout) :: this

    if( associated( this%val_ ) ) then
      deallocate( this%val_ )
      this%val_ => null( )
    end if

  end subroutine component_ptr_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_component
