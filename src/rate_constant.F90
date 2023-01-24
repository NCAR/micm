! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> The abstract rate_constant_t type definition
module micm_rate_constant

  implicit none
  private

  public :: rate_constant_t

  !> A rate constant for a chemical reaction
  type, abstract :: rate_constant_t
  contains
    !> Returns the rate constant for a given set of conditions
    procedure(calculate), deferred :: calculate
  end type rate_constant_t

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the rate constant for a given set of conditions
  subroutine calculate( this, environment, rate_constant )
    use micm_environment,              only : environment_t
    use musica_constants,              only : musica_dk
    use constants,                     only : length 

    import rate_constant_t
    !> Rate constant
    class(rate_constant_t), intent(in) :: this
    !> Environmental conditions
    type(environment_t), intent(in) :: environment(length)
    !> Rate constant
    real(kind=musica_dk), intent(out) :: rate_constant(length)
  end subroutine calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

end module micm_rate_constant
