! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> The rate_constant_wennberg_tunneling_t type
module micm_rate_constant_wennberg_tunneling

  use micm_rate_constant,              only : rate_constant_t
  use musica_constants,                only : musica_dk
  use constants,                       only : length, VLEN, STREAM0

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
  subroutine calculate( this, environment, rate_constant )

    use micm_environment,              only : environment_t

    !> Reaction
    class(rate_constant_wennberg_tunneling_t), intent(in) :: this
    !> Environmental conditions
    type(environment_t), intent(in) :: environment(length)
    !> Rate constant
    real(kind=musica_dk), intent(out) :: rate_constant(length)

    ! Local variable
    integer :: i

    !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
    !$acc loop gang vector
    do i = 1, length
       rate_constant(i) = this%A_ * exp( -this%B_ / environment(i)%temperature + &
                          this%C_ / environment(i)%temperature**3 )
    end do
    !$acc end parallel

  end subroutine calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_rate_constant_wennberg_tunneling
