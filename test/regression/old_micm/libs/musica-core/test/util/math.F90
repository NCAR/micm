! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the musica_math module

!> Test module for the musica_math module
program test_util_math

  use musica_assert
  use musica_math

  implicit none

  character(len=256) :: failure_test_type

  if( command_argument_count( ) .eq. 0 ) then
    call test_math( )
  else if( command_argument_count( ) .eq. 1 ) then
    call get_command_argument( 1, failure_test_type )
    call failure_test( failure_test_type )
  else
    call die( 361435247 )
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Tests math functions
  subroutine test_math( )

    use musica_constants,              only : dk => musica_dk

    real(kind=dk) :: a, b, c(3), p(3)
    real(kind=dk) :: x, y

    ! test chebyshev( )
    ! \todo review Chebyshev equation and add test for in-range x value
    a = 10.0_dk
    b = 30.0_dk
    c = [ 3.0_dk, 2.5_dk, 1.2_dk ]
    x = 32.0_dk
    y = chebyshev( a, b, c, size( c ), x, out_of_bounds_value = 2.0_dk )
    call assert( 525187325, y .eq. 2.0_dk )
    x = -2.0_dk
    y = chebyshev( a, b, c, size( c ), x, out_of_bounds_value = 4.0_dk )
    call assert( 239353303, y .eq. 4.0_dk )

    ! test weighted_chebyshev( )
    c = [ 1.2_dk, 3.5_dk, 4.2_dk ]
    p = [ 4.2_dk, 1.5_dk, 6.2_dk ]
    y = 0.5_dk * c(1) + p(2) * c(2) + p(3) * c(3)
    call assert( 347172356, weighted_chebyshev( 3, c, p ) .eq. y )

    ! test chebyshev_function( )
    x = 12.5_dk
    c(1) = 1.0_dk
    c(2) = x
    c(3) = 2.0_dk * x * x - 1.0_dk
    call chebyshev_function( x, p )
    call assert( 140068949, are_equal( c, p ) )

  end subroutine test_math

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Failure tests for math functions
  subroutine failure_test( test_type )

    character(len=*), intent(in) :: test_type

    if( test_type .eq. "155206939" ) then
      call failure_test_155206939( )
    else if( test_type .eq. "155206939-2" ) then
      call failure_test_155206939_2( )
    else
      call die( 344954102 )
    end if

  end subroutine failure_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Failure test for out-of-bounds x value sent to chebyshev( )
  subroutine failure_test_155206939( )

    use musica_constants,              only : dk => musica_dk

    real(kind=dk) :: a, b, c(3), x, y

    a = 10.0_dk
    b = 20.0_dk
    c = [ 1.2_dk, 4.2_dk, 5.2_dk ]
    x = -5.0_dk
    y = chebyshev( a, b, c, size( c ), x )

  end subroutine failure_test_155206939

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Failure test for out-of-bounds x value sent to chebyshev( )
  subroutine failure_test_155206939_2( )

    use musica_constants,              only : dk => musica_dk

    real(kind=dk) :: a, b, c(3), x, y

    a = 10.0_dk
    b = 20.0_dk
    c = [ 1.2_dk, 4.2_dk, 5.2_dk ]
    x = 40.0_dk
    y = chebyshev( a, b, c, size( c ), x )

  end subroutine failure_test_155206939_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_util_math
