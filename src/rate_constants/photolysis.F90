! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> The rate_constant_photolysis_t type
module micm_rate_constant_photolysis

  use micm_rate_constant,              only : rate_constant_t

  implicit none
  private

  public :: rate_constant_photolysis_t

  !> A photolysis rate constant
  type, extends(rate_constant_t) :: rate_constant_photolysis_t
    private
    integer :: index_ = -1
  contains
    !> Returns the rate constant for a given set of conditions
    procedure :: calculate
  end type rate_constant_photolysis_t

  interface rate_constant_photolysis_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of photolysis rate constants
  elemental function constructor( photolysis_rate_constant_index )            &
      result( new_obj )

    !> New rate constant
    type(rate_constant_photolysis_t) :: new_obj
    !> Index of the rate constant in the photolysis rate constant array
    integer, intent(in) :: photolysis_rate_constant_index

    new_obj%index_ = photolysis_rate_constant_index

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the rate constant for a given set of conditions
  real(kind=musica_dk) elemental function calculate( this, environment )

    use micm_environment,              only : environment_t
    use musica_constants,              only : musica_dk

    !> Reaction
    class(rate_constant_photolysis_t), intent(in) :: this
    !> Environmental conditions
    type(environment_t), intent(in) :: environment

    calculate = environment%photolysis_rate_constants( this%index_ )

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_rate_constant_photolysis
