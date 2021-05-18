! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> The rate_constant_arrhenius_t type
module rate_constant_arrhenius

  implicit none
  private

  public :: rate_constant_arrhenius_t

  !> An Arrhenius rate constant
  type, extends(rate_constant_t) :: rate_constant_arrhenius_t
    private
    real :: A_ = 1.0
    real :: B_ = 0.0
    real :: C_ = 0.0
    real :: D_ = 300.0
    real :: E_ = 0.0
  contains
    !> Returns the rate constant for a given set of conditions
    procedure :: calculate
  end type :: rate_constant_arrhenius_t

  interface rate_constant_arrhenius_t
    module procedure :: constructor
  end interface rate_constant_arrhenius_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of Arrhenius rate constants
  elemental function constructor( A, B, C, D, E, Ea ) result( new_obj )

    use constants,                     only : kBoltzmann

    !> New rate constant
    type(rate_constant_arrhenius_t) :: new_obj
    !> Rate constant parameters
    real, intent(in), optional :: A, B, C, D, E, Ea

    if( present( A  ) ) new_obj%A_ = A
    if( present( B  ) ) new_obj%B_ = B
    if( present( C  ) ) new_obj%C_ = C
    if( present( D  ) ) new_obj%D_ = D
    if( present( E  ) ) new_obj%E_ = E
    if( present( Ea ) ) new_obj%C_ = - Ea / kBoltzmann

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the rate constant for a given set of conditions
  real elemental function calculate( this, environment )

    use environment,                   only : environment_t

    !> Reaction
    class(rate_constant_arrhenius_t), intent(in) :: this
    !> Environmental conditions
    type(environment_t), intent(in) :: environment

    calculate = this%A_ &
      * exp( this%C_ / environment%temperature ) &
      * ( environment%temperature / this%D_ ) ** this%B_ &
      * ( 1.0 + this%E_ *  state%pressure )

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module rate_constant_arrhenius
