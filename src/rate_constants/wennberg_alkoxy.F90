! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> The rate_constant_wennberg_alkoxy_t type
module micm_rate_constant_wennberg_alkoxy

  use micm_rate_constant,              only : rate_constant_t
  use musica_constants,                only : musica_dk

  implicit none
  private

  public :: rate_constant_wennberg_alkoxy_t

  !> A Wennberg NO + RO2 (alkoxy branch) rate constant
  !! (Wennberg et al., Chemical Reviews 118(7):3337-3390, 2018)
  type, extends(rate_constant_t) :: rate_constant_wennberg_alkoxy_t
    private
    real(kind=musica_dk) :: X_    = 1.0
    real(kind=musica_dk) :: Y_    = 0.0
    real(kind=musica_dk) :: a0_   = 1.0
    integer :: n_ = 0.0
  contains
    !> Returns the rate constant for a given set of conditions
    procedure :: calculate
  end type rate_constant_wennberg_alkoxy_t

  interface rate_constant_wennberg_alkoxy_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of Wennberg NO + RO2 (alkoxy branch) rate constants
  function constructor( X, Y, a0, n) result( new_obj )
    !$acc routine seq

    !> New rate constant
    type(rate_constant_wennberg_alkoxy_t) :: new_obj
    !> Rate constant parameters
    real(kind=musica_dk), intent(in), optional :: X, Y, a0
    integer,              intent(in), optional :: n

    if( present( X  ) ) new_obj%X_  = X
    if( present( Y  ) ) new_obj%Y_  = Y
    if( present( a0 ) ) new_obj%a0_ = a0
    if( present( n  ) ) new_obj%n_  = n

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the rate constant for a given set of conditions
  real(kind=musica_dk) function calculate( this, environment )
    !$acc routine seq

    use micm_environment,              only : environment_t

    !> Reaction
    class(rate_constant_wennberg_alkoxy_t), intent(in) :: this
    !> Environmental conditions
    type(environment_t), intent(in) :: environment

    real(kind=musica_dk) :: A, Z

    associate( T => environment%temperature,                                  &
               M => environment%number_density_air )
      A = calculate_A( T, M, this%n_ )
      Z = calculate_A( 293.0_musica_dk, 2.45e19_musica_dk, this%n_ )          &
            * ( 1.0 - this%a0_ ) / this%a0_
      calculate = this%X_ * exp( -this%Y_ / T ) * ( Z / ( Z + A ) )
    end associate

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate A( T, [M], n )
  real(kind=musica_dk) elemental function calculate_A( T, M, n ) result( A )

    !> Temperature [K]
    real(kind=musica_dk), intent(in) :: T
    !> Number density of air [molec cm-3]
    real(kind=musica_dk), intent(in) :: M
    !> Number of heavy atoms in peroxy radical
    integer, intent(in) :: n

    real(kind=musica_dk) :: k0M, kinf

    k0M  = 2e-22 * exp( real( n, kind=musica_dk ) ) * M
    kinf = 0.43 * ( T / 298.0 ) **( -8 )
    A = k0M / ( 1.0 + k0M / kinf )                                            &
        * 0.41**( 1.0 / ( 1.0 + ( log10( k0M / kinf ) )**2 ) )

  end function calculate_A

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_rate_constant_wennberg_alkoxy
