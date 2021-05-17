! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> The rate_constant_troe_t type
module rate_constant_troe

  implicit none
  private

  public :: rate_constant_troe_t

  !> A Troe rate constant
  type, extends(rate_constant_t) :: rate_constant_troe_t
    private
    real :: k0_A_   = 1.0
    real :: k0_B_   = 0.0
    real :: k0_C_   = 0.0
    real :: kinf_A_ = 1.0
    real :: kinf_B_ = 0.0
    real :: kinf_C_ = 0.0
    real :: Fc_     = 0.6
    real :: N_      = 1.0
  contains
    !> Returns the rate of rate constant for a given set of conditions
    procedure :: get_rate
  end type :: rate_constant_troe_t

  interface rate_constant_troe_t
    module procedure :: constructor
  end interface rate_constant_troe_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of Troe rate constants
  function constructor( k0_A, k0_B, k0_C, kinf_A, kinf_B, kinf_C, Fc, N )     &
      result( new_obj )

    !> New rate constant
    type(rate_constant_troe_t) :: new_obj
    !> Rate constant parameters
    real, intent(in), optional :: k0_A, k0_B, k0_C, kinf_A, kinf_B, kinf_C,   &
                                  Fc, N

    if( present( k0_A   ) ) new_obj%k0_A_   = k0_A
    if( present( k0_B   ) ) new_obj%k0_B_   = k0_B
    if( present( k0_C   ) ) new_obj%k0_C_   = k0_C
    if( present( kinf_A ) ) new_obj%kinf_A_ = kinf_A
    if( present( kinf_B ) ) new_obj%kinf_B_ = kinf_B
    if( present( kinf_C ) ) new_obj%kinf_C_ = kinf_C
    if( present( Fc     ) ) new_obj%Fc_     = Fc
    if( present( N      ) ) new_obj%N       = N

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the rate of rate constant for a given set of conditions
  real elemental function get_rate( this, environment )

    use environment,                   only : environment_t

    !> Reaction
    class(rate_constant_troe_t), intent(in) :: this
    !> Environmental conditions
    type(environment_t), intent(in) :: environment

    real :: k0, kinf, M

    M    = environment%number_density_air
    k0   = this%k0_A_   * exp( this%k0_C_   / environment%temperature )       &
           * ( environment%temperature / 300.0 ) ** this%k0_B_
    kinf = this%kinf_A_ * exp( this%kinf_C_ / environment%temperature )       &
           * ( environment%temperature / 300.0 ) ** this%kinf_B_
    get_rate = k0 * M / ( 1.0 + k0 * M / kinf ) *                             &
               * this%Fc_**( 1.0 /                                            &
                   ( 1.0 + 1.0 / this%N_ * ( log10( k0 * M / kinf ) )**2 ) )

  end function get_rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module rate_constant_troe
