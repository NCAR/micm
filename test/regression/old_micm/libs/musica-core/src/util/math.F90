! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The mam_math_util module

!> Common math functions
!!
!! \todo review naming and descriptions of math functions
module musica_math

  implicit none
  private

  public :: chebyshev, weighted_chebyshev, chebyshev_function

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the value of a Chebyshev polynomial at a given point
  !!
  !! The algoritm is based on \cite William1992 page 193.
  !!
  !! If an out-of-bounds value is provided and x is out of the range
  !! \f$a<=x<=b\f$, the out-of-bounds value will be returned. If an
  !! out-of-bounds value is not provided and x is out-of-bounds, an error
  !! is thrown.
  real(kind=musica_dk) function chebyshev( a, b, c, m, x, out_of_bounds_value )

    use musica_assert,                 only : die
    use musica_constants,              only : musica_dk, musica_ik

    !> Lower bound of parameterization
    real(kind=musica_dk), intent(in) :: a
    !> Upper bound of parameterization
    real(kind=musica_dk), intent(in) :: b
    !> Chebyshev coefficients c[1...m]
    real(kind=musica_dk), intent(in) :: c(:)
    !> Number of elements of c[] to use in calculation
    integer(kind=musica_ik), intent(in) :: m
    !> Independent variable
    real(kind=musica_dk), intent(in) :: x
    !> Out-of-bounds value
    real(kind=musica_dk), intent(in), optional :: out_of_bounds_value

    integer(kind=musica_ik) :: j
    real(kind=musica_dk) :: d, dd, sv, y, y2

    if( ( x - a ) * ( x - b ) > 0._musica_dk ) then
      if( present( out_of_bounds_value ) ) then
        chebyshev = out_of_bounds_value
      else
        call die( 155206939 )
      end if
    else
      d  = 0._musica_dk
      dd = 0._musica_dk
      y  = ( 2._musica_dk * x - a - b ) / ( b - a )
      y2 = 2._musica_dk * y
      do j = m, 2, -1
        sv = d
        d  = y2 * d - dd + c( j )
        dd = sv
      end do
      chebyshev = y * d - dd + 0.5_musica_dk * c( 1 )
    end if

  end function chebyshev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculates a value from a weighted Chebyshev polynomial
  !!
  !! \todo is "weighted Chebyshev polynomial" the right term to use here?
  pure function weighted_chebyshev( number_of_coefficients, coefficients,     &
      polynomials ) result( value )

    use musica_constants,              only : musica_dk

    real(kind=musica_dk)             :: value
    integer,              intent(in) :: number_of_coefficients
    real(kind=musica_dk), intent(in) :: coefficients( number_of_coefficients )
    real(kind=musica_dk), intent(in) :: polynomials(  number_of_coefficients )

    integer :: i

    value = 0.5_musica_dk * coefficients(1)
    do i = 2, number_of_coefficients
      value = value + coefficients( i ) * polynomials( i )
    end do

  end function weighted_chebyshev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns a Chebyshev polynomial function for a given independent variable
  !!
  !! The number of polynomials returned corresponds to the size of the
  !! polynomials array
  !!
  !! \todo is "Chebyshev polynomial function" the correct name for what is
  !!       returned from math::chebyshev_function?
  pure subroutine chebyshev_function( independent_variable, polynomials )

    use musica_constants,              only : musica_dk

    !> Value to use for the independent variable in calculating the function
    real(kind=musica_dk), intent(in)  :: independent_variable
    !> Resulting function elements
    !!
    !! The number of polynomials is based on the size of the array passed to
    !! this function. The array must have at least 2 elements.
    real(kind=musica_dk), intent(out) :: polynomials(:)

    integer :: i

    associate( x => independent_variable )
      polynomials(1) = 1.0_musica_dk
      polynomials(2) = x
      do i = 3, size( polynomials )
        polynomials( i ) = 2.0_musica_dk * x * polynomials( i - 1 )           &
                           - polynomials( i - 2 )
      end do
    end associate

  end subroutine chebyshev_function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_math
