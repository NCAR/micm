! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the musica_interpolator_linear_1D module

!> Test module for the linaer 1D interpolator strategy
program test_interpolator_linear_1D

  implicit none

  call test_linear_strategy( )
  call example( )

contains

  subroutine test_linear_strategy( )

    use musica_assert,                 only : assert
    use musica_constants,              only : musica_dk
    use musica_grid,                   only : grid_t
    use musica_interpolator,           only : interpolator_t
    use musica_interpolator_linear_1D, only : strategy
    use musica_string,                 only : string_t

    type(interpolator_t) :: forward, backward
    type(grid_t) :: from_grid, to_grid
    type(string_t) :: units
    real(kind=musica_dk), allocatable :: lower_bounds(:), upper_bounds(:)
    real(kind=musica_dk), allocatable :: from_array(:), to_array(:)

    units = "foos"

    lower_bounds = (/ 1.0,  5.0, 20.0, 25.0 /)
    upper_bounds = (/ 5.0, 10.0, 25.0, 45.0 /)
    from_grid = grid_t( lower_bounds, upper_bounds, units )

    lower_bounds = (/ 0.0,  5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0 /)
    upper_bounds = (/ 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0 /)
    to_grid = grid_t( lower_bounds, upper_bounds, units )

    allocate( from_array( from_grid%number_of_sections( ) ) )
    allocate( to_array(   to_grid%number_of_sections( )   ) )

    forward  = interpolator_t( strategy, from_grid, to_grid   )
    backward = interpolator_t( strategy, to_grid,   from_grid )

    from_array = (/ 10.0, 60.0, 4.0, 200.0 /)
    to_array   = -9999.0

    call forward%interpolate( from_array, to_array )

    call assert( 371766980, to_array(1) .eq.  10.0 )
    call assert( 643713615, to_array(2) .eq.  60.0 )
    call assert( 694834628, to_array(3) .eq.   0.0 )
    call assert( 749314414, to_array(4) .eq.   0.0 )
    call assert( 468744699, to_array(5) .eq.  54.0 )
    call assert( 235484929, to_array(6) .eq. 100.0 )
    call assert( 793717881, to_array(7) .eq.  50.0 )
    call assert( 745554540, to_array(8) .eq.   0.0 )

    call backward%interpolate( to_array, from_array )

    call assert( 837006426, from_array(1) .eq.   8.0 )
    call assert( 141490479, from_array(2) .eq.  60.0 )
    call assert( 871333574, from_array(3) .eq.  27.0 )
    call assert( 701176670, from_array(4) .eq. 152.0 )

  end subroutine test_linear_strategy

  subroutine example( )

    use musica_assert,                 only : assert

    !! [Linear interpolator example]
    use musica_constants,              only : musica_dk
    use musica_grid,                   only : grid_t
    use musica_interpolator,           only : interpolator_t
    use musica_interpolator_linear_1D, only : strategy
    use musica_string,                 only : string_t

    type(interpolator_t) :: a
    type(grid_t) :: from_grid, to_grid
    type(string_t) :: units
    real(kind=musica_dk), allocatable :: lower_bounds(:), upper_bounds(:)
    real(kind=musica_dk) :: x(3), y(6)

    units = "foos"

    lower_bounds = (/  0.0, 10.0, 20.0 /)
    upper_bounds = (/ 10.0, 20.0, 30.0 /)
    from_grid = grid_t( lower_bounds, upper_bounds, units )

    lower_bounds = (/ 0.0,  5.0, 10.0, 15.0, 20.0, 25.0 /)
    upper_bounds = (/ 5.0, 10.0, 15.0, 20.0, 25.0, 30.0 /)
    to_grid = grid_t( lower_bounds, upper_bounds, units )

    x = (/ 10.0, 20.0, 30.0 /)

    a = interpolator_t( strategy, from_grid, to_grid )

    call a%interpolate( x, y )

    write(*,*) y
    !! [Linear interpolator example]

    call assert( 640718756, y(1) .eq.  5.0 )
    call assert( 186633364, y(2) .eq.  5.0 )
    call assert( 916476459, y(3) .eq. 10.0 )
    call assert( 746319555, y(4) .eq. 10.0 )
    call assert( 858637900, y(5) .eq. 15.0 )
    call assert( 688480996, y(6) .eq. 15.0 )

  end subroutine example

end program test_interpolator_linear_1D
