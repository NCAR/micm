! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> This musica_rate_constant_bar module

!> The rate_constant_bar_t type and related functions
module micm_rate_constant_bar

  use micm_rate_constant,              only : rate_constant_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none
  private

  public :: rate_constant_bar_t

  !> Calculator of bar-type rate constants
  type, extends(rate_constant_t) :: rate_constant_bar_t
    !> The 'A' parameter
    real(kind=musica_dk) :: A_
    !> The 'B' parameter
    real(kind=musica_dk) :: B_
    !> The 'C' parameter
    real(kind=musica_dk) :: C_
  contains
    !> Calculate the rate constant
    procedure :: calculate
  end type rate_constant_bar_t

  !> Constructor of rate_constant_bar_t objects
  interface rate_constant_bar_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of rate_constant_bar_t objects
  function constructor( config ) result( new_obj )

    use musica_assert,                 only : die
    use musica_config,                 only : config_t

    !> New bar rate constant calculator
    class(rate_constant_t), pointer :: new_obj
    !> Rate constant configuration data
    type(config_t), intent(inout) :: config

    character(len=*), parameter :: my_name = 'bar rate constant constructor'

    allocate( rate_constant_bar_t :: new_obj )
    select type( new_obj )
    class is( rate_constant_bar_t )

      ! 'A' must be specified
      call config%get( 'A', new_obj%A_, my_name )

      ! 'B' is optional with a default value of 12.0
      call config%get( 'B', new_obj%B_, my_name, default = 12.0_musica_dk )

      ! The 'C' parameter must be specified
      call config%get( 'C', new_obj%C_, my_name )

    class default
      call die( 649274138 )
    end select

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the rate constant for a given set of environmental conditions
  elemental function calculate( this, environment ) result( rate_constant )

    use micm_environment,              only : environment_t

    !> Calculated rate constant
    real(kind=musica_dk) :: rate_constant
    !> bar rate constant calculator
    class(rate_constant_bar_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in) :: environment

    ! k = A * exp( T / B ) - C * P
    rate_constant = this%A_ * exp( environment%temperature / this%B_ ) -      &
                    this%C_ * environment%pressure

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_rate_constant_bar
