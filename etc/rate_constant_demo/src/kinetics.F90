! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The micm_kinetics module

!> The kinetics_t type and related functions
!!
!! \todo Consider rewriting micm_kinetics to follow musica style conventions
module micm_kinetics

  use micm_rate_constant,              only : rate_constant_ptr
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: kinetics_t

  !> Number of reactions (from the preprocessor)
  integer(kind=musica_ik), parameter :: kNumberOfReactions = 4

  !> Kinetics
  type :: kinetics_t
    !> Rate constant calculators
    type(rate_constant_ptr), allocatable :: rate_constants_(:)
    !> Current rate constant values
    real(kind=musica_dk), allocatable :: rate_constant_values_(:)
  contains
    !> Update the object for new environmental conditions
    procedure :: update_for_new_environmental_state
    !> Finalize the object
    final :: finalize
  end type kinetics_t

  !> kinetics_t constructor
  interface kinetics_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of kinetics_t objects
  function constructor( config ) result( new_obj )

    use micm_rate_constant_factory,    only : rate_constant_builder
    use musica_assert,                 only : assert
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t

    !> New kinetics
    class(kinetics_t), pointer :: new_obj
    !> Kinetics configuration data
    type(config_t), intent(inout) :: config

    character(len=*), parameter :: my_name = "Kinetics constructor"
    type(config_t) :: reaction_set, reaction_config, rate_constant_config
    class(iterator_t), pointer :: iter
    integer(kind=musica_ik) :: i_reaction

    allocate( new_obj )
    allocate( new_obj%rate_constants_(       kNumberOfReactions ) )
    allocate( new_obj%rate_constant_values_( kNumberOfReactions ) )
    call config%get( "reactions", reaction_set, my_name )
    iter => reaction_set%get_iterator( )
    i_reaction = 0
    do while( iter%next( ) )
      i_reaction = i_reaction + 1
      call reaction_set%get( iter, reaction_config, my_name )
      call reaction_config%get( "rate constant", rate_constant_config,        &
                                my_name )
      new_obj%rate_constants_( i_reaction )%val_ =>                           &
          rate_constant_builder( rate_constant_config )
      call reaction_config%finalize( )
      call rate_constant_config%finalize( )
    end do
    call reaction_set%finalize( )
    deallocate( iter )
    call assert( 618706864, i_reaction .eq. kNumberOfReactions )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update the object for new environmental conditions
  subroutine update_for_new_environmental_state( this, environment )

    use micm_environment,              only : environment_t

    !> Kinetics
    class(kinetics_t), intent(inout) :: this
    !> Environmental conditions
    class(environment_t), intent(in) :: environment

    integer(kind=musica_ik) :: i_rc

    do i_rc = 1, kNumberOfReactions
      this%rate_constant_values_( i_rc ) =                                    &
          this%rate_constants_( i_rc )%val_%calculate( environment )
      write(*,*) "Rate constant ", i_rc, ": ",                                &
                 this%rate_constant_values_( i_rc )
    end do

  end subroutine update_for_new_environmental_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the kinetics
  subroutine finalize( this )

    !> Kinetics
    type(kinetics_t), intent(inout) :: this

    integer(kind=musica_ik) :: i_rate_constant

    if( allocated( this%rate_constants_ ) ) then
      do i_rate_constant = 1, size( this%rate_constants_ )
        if( associated( this%rate_constants_( i_rate_constant )%val_ ) ) then
          deallocate( this%rate_constants_( i_rate_constant )%val_ )
        end if
      end do
      deallocate( this%rate_constants_ )
    end if

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_kinetics
