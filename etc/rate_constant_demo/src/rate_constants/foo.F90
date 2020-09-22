! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> This musica_rate_constant_foo module

!> The rate_constant_foo_t type and related functions
module micm_rate_constant_foo

  use micm_rate_constant,              only : rate_constant_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none
  private

  public :: rate_constant_foo_t

  !> Calculator of foo-type rate constants
  type, extends(rate_constant_t) :: rate_constant_foo_t
    !> The 'A' parameter
    real(kind=musica_dk) :: A_
    !> The 'B' parameter
    real(kind=musica_dk) :: B_
  contains
    !> Calculate the rate constant
    procedure :: calculate
  end type rate_constant_foo_t

  !> Constructor of rate_constant_foo_t objects
  interface rate_constant_foo_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of rate_constant_foo_t objects
  function constructor( config ) result( new_obj )

    use musica_assert,                 only : die
    use musica_config,                 only : config_t

    !> New foo rate constant calculator
    class(rate_constant_t), pointer :: new_obj
    !> Rate constant configuration data
    class(config_t), intent(inout) :: config

    character(len=*), parameter :: my_name = 'foo rate constant constructor'

    allocate( rate_constant_foo_t :: new_obj )
    select type( new_obj )
    class is( rate_constant_foo_t )

      ! 'A' must be specified
      call config%get( 'A', new_obj%A_, my_name )

      ! 'B' is optional with a default value of 12.0
      call config%get( 'B', new_obj%B_, my_name, default = 12.0_musica_dk )

    class default
      call die( 121270143 )
    end select

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the rate constant for a given set of environmental conditions
  elemental function calculate( this, environment ) result( rate_constant )

    use micm_environment,              only : environment_t

    !> Calculated rate constant
    real(kind=musica_dk) :: rate_constant
    !> foo rate constant calculator
    class(rate_constant_foo_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in) :: environment

    ! k = A * T + B * P
    rate_constant = this%A_ * environment%temperature +                       &
                    this%B_ * environment%pressure

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_rate_constant_foo
