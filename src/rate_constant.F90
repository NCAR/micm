! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> The abstract rate_constant_t type definition
module rate_constant

  implicit none
  private

  public :: rate_constant_t

  !> A rate constant for a chemical reaction
  type, abstract :: rate_constant_t
  contains
    !> Returns the rate constant for a given set of conditions
    procedure :: calculate
  end type :: rate_constant_t

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the rate constant for a given set of conditions
  function calculate( this, environment )
    use environment,                   only : environment_t
    import rate_constant_t
    !> Rate constant
    class(rate_constant_t), intent(in) :: this
    !> Environmental conditions
    type(environment_t), intent(in) :: environment
  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

end module rate_constant
