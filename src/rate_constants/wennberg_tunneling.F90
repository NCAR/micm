! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> The rate_constant_wennberg_tunneling_t type
module rate_constant_wennberg_tunneling

  implicit none
  private

  public :: rate_constant_wennberg_tunneling_t

  !> A Wennberg tunneling rate constant
  !! (Wennberg et al., Chemical Reviews 118(7):3337-3390, 2018)
  type, extends(rate_constant_t) :: rate_constant_wennberg_tunneling_t
    private
    real :: A_ = 1.0
    real :: B_ = 0.0
    real :: C_ = 0.0
  contains
    !> Returns the rate constant for a given set of conditions
    procedure :: calculate
  end type :: rate_constant_wennberg_tunneling_t

  interface rate_constant_wennberg_tunneling_t
    module procedure :: constructor
  end interface rate_constant_wennberg_tunneling_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of Wennberg tunneling rate constants
  elemental function constructor( A, B, C ) result( new_obj )

    !> New rate constant
    type(rate_constant_wennberg_tunneling_t) :: new_obj
    !> Rate constant parameters
    real, intent(in), optional :: A, B, C

    if( present( A ) ) new_obj%A_ = A
    if( present( B ) ) new_obj%B_ = B
    if( present( C ) ) new_obj%C_ = C

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the rate constant for a given set of conditions
  real elemental function calculate( this, environment )

    use environment,                   only : environment_t

    !> Reaction
    class(rate_constant_wennberg_tunneling_t), intent(in) :: this
    !> Environmental conditions
    type(environment_t), intent(in) :: environment

    associate( T => environment%temperature )
      calculate = this%A_ * exp( -this%B_ / T + this%C_ / T**3 )
    end associate

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module rate_constant_wennberg_tunneling
