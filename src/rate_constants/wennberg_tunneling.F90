! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> The rate_constant_wennberg_tunneling_t type
module micm_rate_constant_wennberg_tunneling

  use micm_rate_constant,              only : rate_constant_t
  use musica_constants,                only : musica_dk

  implicit none
  private

  public :: rate_constant_wennberg_tunneling_t

  !> A Wennberg tunneling rate constant
  !! (Wennberg et al., Chemical Reviews 118(7):3337-3390, 2018)
  type, extends(rate_constant_t) :: rate_constant_wennberg_tunneling_t
    private
    real(kind=musica_dk) :: A_ = 1.0
    real(kind=musica_dk) :: B_ = 0.0
    real(kind=musica_dk) :: C_ = 0.0
  contains
    !> Returns the rate constant for a given set of conditions
    procedure :: calculate
  end type rate_constant_wennberg_tunneling_t

  interface rate_constant_wennberg_tunneling_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of Wennberg tunneling rate constants
  function constructor( A, B, C ) result( new_obj )
    !$acc routine seq

    !> New rate constant
    type(rate_constant_wennberg_tunneling_t) :: new_obj
    !> Rate constant parameters
    real(kind=musica_dk), intent(in), optional :: A, B, C

    if( present( A ) ) new_obj%A_ = A
    if( present( B ) ) new_obj%B_ = B
    if( present( C ) ) new_obj%C_ = C

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the rate constant for a given set of conditions
  real(kind=musica_dk) function calculate( this, environment )
    !$acc routine seq

    use micm_environment,              only : environment_t

    !> Reaction
    class(rate_constant_wennberg_tunneling_t), intent(in) :: this
    !> Environmental conditions
    type(environment_t), intent(in) :: environment

    associate( T => environment%temperature )
      calculate = this%A_ * exp( -this%B_ / T + this%C_ / T**3 )
    end associate

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_rate_constant_wennberg_tunneling
