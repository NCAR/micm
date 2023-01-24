! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> The rate_constant_troe_t type
module micm_rate_constant_troe

  use micm_rate_constant,              only : rate_constant_t
  use musica_constants,                only : musica_dk
  use constants,                       only : length, VLEN, STREAM0

  implicit none
  private

  public :: rate_constant_troe_t

  !> A Troe rate constant
  type, extends(rate_constant_t) :: rate_constant_troe_t
    private
    real(kind=musica_dk) :: k0_A_   = 1.0
    real(kind=musica_dk) :: k0_B_   = 0.0
    real(kind=musica_dk) :: k0_C_   = 0.0
    real(kind=musica_dk) :: kinf_A_ = 1.0
    real(kind=musica_dk) :: kinf_B_ = 0.0
    real(kind=musica_dk) :: kinf_C_ = 0.0
    real(kind=musica_dk) :: Fc_     = 0.6
    real(kind=musica_dk) :: N_      = 1.0
  contains
    !> Returns the rate constant for a given set of conditions
    procedure :: calculate
  end type rate_constant_troe_t

  interface rate_constant_troe_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of Troe rate constants
  function constructor( k0_A, k0_B, k0_C, kinf_A, kinf_B, kinf_C,   &
      Fc, N ) result( new_obj )

    !> New rate constant
    type(rate_constant_troe_t) :: new_obj
    !> Rate constant parameters
    real(kind=musica_dk), intent(in), optional :: k0_A, k0_B, k0_C, kinf_A,   &
                                                  kinf_B, kinf_C, Fc, N

    if( present( k0_A   ) ) new_obj%k0_A_   = k0_A
    if( present( k0_B   ) ) new_obj%k0_B_   = k0_B
    if( present( k0_C   ) ) new_obj%k0_C_   = k0_C
    if( present( kinf_A ) ) new_obj%kinf_A_ = kinf_A
    if( present( kinf_B ) ) new_obj%kinf_B_ = kinf_B
    if( present( kinf_C ) ) new_obj%kinf_C_ = kinf_C
    if( present( Fc     ) ) new_obj%Fc_     = Fc
    if( present( N      ) ) new_obj%N_      = N

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the rate constant for a given set of conditions
  subroutine calculate( this, environment, rate_constant )

    use micm_environment,              only : environment_t

    !> Reaction
    class(rate_constant_troe_t), intent(in) :: this
    !> Environmental conditions
    type(environment_t), intent(in) :: environment(length)
    !> Rate constant
    real(kind=musica_dk), intent(out) :: rate_constant(length)

    ! Local variable
    integer :: i
    real(kind=musica_dk) :: k0, kinf, M

    !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
    !$acc loop gang vector
    do i = 1, length
       M    = environment(i)%number_density_air
       k0   = this%k0_A_   * exp( this%k0_C_   / environment(i)%temperature )       &
              * ( environment(i)%temperature / 300.0 ) ** this%k0_B_
       kinf = this%kinf_A_ * exp( this%kinf_C_ / environment(i)%temperature )       &
              * ( environment(i)%temperature / 300.0 ) ** this%kinf_B_
       rate_constant(i) = k0 * M / ( 1.0 + k0 * M / kinf )                          &
                  * this%Fc_**( 1.0 /                                            &
                      ( 1.0 + 1.0 / this%N_ * ( log10( k0 * M / kinf ) )**2 ) )
    end do
    !$acc end parallel

  end subroutine calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_rate_constant_troe
