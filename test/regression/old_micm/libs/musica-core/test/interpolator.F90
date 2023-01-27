! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the musica_interpolator module

!> Test module for the musica_interpolator module
module test_interpolator

  use musica_assert
  use musica_constants,                only : musica_dk
  use musica_interpolator

  implicit none

contains

  function foo_mapper( from_grid, to_grid ) result( map )

    use musica_grid,                   only : grid_t, grid_iterator_t

    class(interpolator_element_t), allocatable :: map(:)
    class(grid_iterator_t), pointer :: iter_from, iter_to
    class(grid_t), intent(in) :: from_grid
    class(grid_t), intent(in) :: to_grid
    type(interpolator_element_t) :: local_map(3)

    iter_from => from_grid%iterator( )
    iter_to   => to_grid%iterator( )

    call assert( 306034352, iter_from%next( ) )
    call assert( 306034352, iter_to%next( ) )
    ! to(1) = from(1) * 1.0
    local_map(1) = interpolator_element_t( iter_from, iter_to, 1.0_musica_dk )

    call assert( 306034352, iter_from%next( ) )
    ! to(1) = from(2) * 0.5
    local_map(2) = interpolator_element_t( iter_from, iter_to, 0.5_musica_dk )

    call assert( 306034352, iter_to%next( ) )
    call assert( 306034352, iter_to%next( ) )
    ! to(3) = from(2) * 0.9
    local_map(3) = interpolator_element_t( iter_from, iter_to, 0.9_musica_dk )

    allocate( map, source = local_map )

    deallocate( iter_from )
    deallocate( iter_to   )

  end function foo_mapper

end module test_interpolator

program test_interpolator_driver

  implicit none

  call run_tests( )

contains

  subroutine run_tests( )

    use musica_grid,                     only : grid_t
    use musica_string,                   only : string_t
    use test_interpolator

    implicit none

    type(interpolator_t) :: interpolator
    type(grid_t) :: from_grid, to_grid
    type(string_t) :: units
    real(kind=musica_dk) :: lower(3), upper(3) ! values aren't used in test strategy
    real(kind=musica_dk) :: from_array(3), to_array(3)

    units        = "foos"
    lower        = 0.0
    upper        = 0.0
    from_grid    = grid_t( lower, upper, units )
    to_grid      = grid_t( lower, upper, units )
    interpolator = interpolator_t( foo_mapper, from_grid, to_grid )

    from_array = (/ 20.0, 400.0, 60.0 /)
    to_array   = (/ 9999.0, 9999.0, 9999.0 /)

    call interpolator%interpolate( from_array, to_array )

    call assert( 239899440, to_array(1) .eq. 220.0 )
    call assert( 736482765, to_array(2) .eq. 0.0   )
    call assert( 291020453, to_array(3) .eq. 360.0 )

  end subroutine run_tests

end program test_interpolator_driver
