! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The micm_rate_constant module

!> The abstract rate_constant_t type and related functions
module micm_rate_constant

  use musica_constants,                only : musica_dk, musica_ik
  use micm_environment,                only : environment_t

  implicit none
  private

  public :: rate_constant_t, rate_constant_ptr

  !> Rate constant calculator
  type, abstract :: rate_constant_t
  contains
    !> Calculate the rate constant
    procedure(calculate), deferred :: calculate
  end type rate_constant_t

  !> Pointer type for building sets of rate constants
  type :: rate_constant_ptr
    class(rate_constant_t), pointer :: val_ => null( )
  end type rate_constant_ptr

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the rate constant for a given set of conditions
  elemental function calculate( this, environment ) result( rate_constant )
    use micm_environment,              only : environment_t
    use musica_constants,              only : musica_dk
    import rate_constant_t
    !> Rate constant
    real(kind=musica_dk) :: rate_constant
    !> Rate constant calculator
    class(rate_constant_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in) :: environment
  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

end module micm_rate_constant
