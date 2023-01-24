! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> The rate_constant_arrhenius_t type
module micm_rate_constant_arrhenius

  use musica_constants,                only : musica_dk
  use micm_rate_constant,              only : rate_constant_t
  use constants,                       only : length, VLEN, STREAM0

  implicit none
  private

  public :: rate_constant_arrhenius_t

  !> An Arrhenius rate constant
  type, extends(rate_constant_t) :: rate_constant_arrhenius_t
    private
    real(kind=musica_dk) :: A_ = 1.0
    real(kind=musica_dk) :: B_ = 0.0
    real(kind=musica_dk) :: C_ = 0.0
    real(kind=musica_dk) :: D_ = 300.0
    real(kind=musica_dk) :: E_ = 0.0
  contains
    !> Returns the rate constant for a given set of conditions
    procedure :: calculate
  end type rate_constant_arrhenius_t

  interface rate_constant_arrhenius_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of Arrhenius rate constants
  function constructor( A, B, C, D, E, Ea ) result( new_obj )

    use constants,                     only : kBoltzmann

    !> New rate constant
    type(rate_constant_arrhenius_t) :: new_obj
    !> Rate constant parameters
    real(kind=musica_dk), intent(in), optional :: A, B, C, D, E, Ea

    if( present( A  ) ) new_obj%A_ = A
    if( present( B  ) ) new_obj%B_ = B
    if( present( C  ) ) new_obj%C_ = C
    if( present( D  ) ) new_obj%D_ = D
    if( present( E  ) ) new_obj%E_ = E
    if( present( Ea ) ) new_obj%C_ = - Ea / kBoltzmann

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the rate constant for a given set of conditions
  subroutine calculate( this, environment, rate_constant )

    use micm_environment,              only : environment_t

    !> Reaction
    class(rate_constant_arrhenius_t), intent(in) :: this
    !> Environmental conditions
    type(environment_t), intent(in) :: environment(length)
    !> Rate constant
    real(kind=musica_dk), intent(out) :: rate_constant(length)

    ! Local variable
    integer :: i

    !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
    !$acc loop gang vector
    do i = 1, length
       rate_constant(i) = this%A_ &
           * exp( this%C_ / environment(i)%temperature ) &
           * ( environment(i)%temperature / this%D_ ) ** this%B_ &
           * ( 1.0 + this%E_ *  environment(i)%pressure )
    end do
    !$acc end parallel

  end subroutine calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_rate_constant_arrhenius
