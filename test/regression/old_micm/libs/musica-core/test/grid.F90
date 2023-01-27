! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the musica_grid module

!> Test module for the musica_grid module
module test_grid

  use musica_assert
  use musica_grid

  implicit none
  private

  public :: test_grid_t

  type, extends(grid_t) :: foo_grid_t
  end type foo_grid_t

  interface foo_grid_t
    procedure :: foo_grid_constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for foo_grid_t
  function foo_grid_constructor( ) result( foo )

    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t

    type(foo_grid_t) :: foo

    type(string_t) :: units
    real(kind=musica_dk) :: lower_bounds(3), upper_bounds(3)

    lower_bounds = (/ 10.0, 100.0, 200.0 /)
    upper_bounds = (/ 60.0, 150.0, 300.0 /)
    units = "foos"

    call foo%private_constructor_bounds( lower_bounds, upper_bounds, units )

  end function foo_grid_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test grid_t functionality
  subroutine test_grid_t( )

    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t

    type(grid_t) :: a, b, c
    type(foo_grid_t) :: foo

    real(kind=musica_dk), allocatable :: lower_bounds(:), upper_bounds(:)
    real(kind=musica_dk), allocatable :: compare_bounds(:)
    type(grid_section_t) :: a_sect, b_sect, foo_sect
    class(grid_iterator_t), pointer :: a_iter, b_iter, foo_iter
    type(string_t) :: units

    lower_bounds = (/ 10.0, 100.0, 200.0 /)
    upper_bounds = (/ 60.0, 150.0, 300.0 /)
    units = "foos"
    a = grid_t( lower_bounds, upper_bounds, units )
    c = grid_t( lower_bounds, upper_bounds, units )

    lower_bounds(2) = 125.0
    upper_bounds(2) = 174.0
    b = grid_t( lower_bounds, upper_bounds, units )

    deallocate( lower_bounds )
    deallocate( upper_bounds )

    foo = foo_grid_t( )

    ! test comparison operators
    call assert( 222304556, a .eq. c )
    call assert( 901478933, a .ne. b )
    call assert( 105174199, c .ne. b )
    call assert( 664860390, a .eq. foo )
    call assert( 389554982, foo .eq. c )
    call assert( 831658521, foo .ne. b )

    ! test grid property accessors
    call assert( 501446629, a%number_of_sections( )   .eq. 3 )
    call assert( 498087856, foo%number_of_sections( ) .eq. 3 )

    call assert( 149150968, c%units( )   .eq. "foos" )
    call assert( 886163904, foo%units( ) .eq. "foos" )

    a_iter   => a%iterator( )
    b_iter   => b%iterator( )
    foo_iter => foo%iterator( )

    ! grid section 1
    call assert( 200083629, a_iter%next( ) )
    call assert( 142245070, b_iter%next( ) )
    call assert( 937096565, foo_iter%next( ) )

    call assert( 807984355, a%lower_bound( a_iter )     .eq. 10.0 )
    call assert( 253849594, b%lower_bound( b_iter )     .eq. 10.0 )
    call assert( 701217440, foo%lower_bound( foo_iter ) .eq. 10.0 )

    call assert( 531060536, a%upper_bound( a_iter )     .eq. 60.0 )
    call assert( 425912032, b%upper_bound( b_iter )     .eq. 60.0 )
    call assert( 873279878, foo%upper_bound( foo_iter ) .eq. 60.0 )

    call assert( 985598223, a%range( a_iter )     .eq. 60.0 - 10.0 )
    call assert( 815441319, b%range( b_iter )     .eq. 60.0 - 10.0 )
    call assert( 645284415, foo%range( foo_iter ) .eq. 60.0 - 10.0 )

    call assert( 284078551, a%property_index( a_iter )     .eq. 1 )
    call assert( 475127511, b%property_index( b_iter )     .eq. 1 )
    call assert( 922495357, foo%property_index( foo_iter ) .eq. 1 )

    a_sect   = a%section( a_iter )
    b_sect   = b%section( b_iter )
    foo_sect = foo%section( foo_iter )

    call assert( 389679350, a_sect .eq. b_sect )
    call assert( 134813703, a_sect .eq. foo_sect )

    call assert( 388226111, a_sect%overlap( b_sect )   .eq. 50.0 )
    call assert( 864656798, a_sect%overlap( foo_sect ) .eq. 50.0 )

    call assert( 997127781, a_sect%lower_bound( )   .eq. 10.0 )
    call assert( 976975143, b_sect%lower_bound( )   .eq. 10.0 )
    call assert( 189293489, foo_sect%lower_bound( ) .eq. 10.0 )

    call assert( 919136584, a_sect%upper_bound( )   .eq. 60.0 )
    call assert( 748979680, b_sect%upper_bound( )   .eq. 60.0 )
    call assert( 296347527, foo_sect%upper_bound( ) .eq. 60.0 )

    call assert( 126190623, a_sect%range( )   .eq. 50.0 )
    call assert( 921042118, b_sect%range( )   .eq. 50.0 )
    call assert( 750885214, foo_sect%range( ) .eq. 50.0 )

    ! grid section 2
    call assert( 979594498, a_iter%next( ) )
    call assert( 526962345, b_iter%next( ) )
    call assert( 421813841, foo_iter%next( ) )

    a_sect   = a%section( a_iter )
    b_sect   = b%section( b_iter )
    foo_sect = foo%section( foo_iter )

    call assert( 427390868, a_sect   .ne. b_sect )
    call assert( 188866791, foo_sect .ne. b_sect )
    call assert( 190772325, foo_sect .eq. a_sect )

    call assert( 294467590, a_sect%overlap( b_sect )   .eq. 150.0 - 125.0 )
    call assert( 856059315, a_sect%overlap( foo_sect ) .eq. 150.0 - 100.0 )

    call assert( 629969386, a_sect%range( )   .eq. 150.0 - 100.0 )
    call assert( 284391271, b_sect%range( )   .eq. 174.0 - 125.0 )
    call assert( 346040898, foo_sect%range( ) .eq. 150.0 - 100.0 )

    ! grid section 3
    call assert( 251656937, a_iter%next( ) )
    call assert( 699024783, b_iter%next( ) )
    call assert( 246392630, foo_iter%next( ) )

    ! end of grid
    call assert( 358710975, .not. a_iter%next( ) )
    call assert( 188554071, .not. b_iter%next( ) )
    call assert( 635921917, .not. foo_iter%next( ) )

    ! check overlap function
    call assert( 773356593, a%overlap( c ) .eq. 200.0 )
    call assert( 371845453, a%overlap( b ) .eq. 175.0 )

    deallocate( a_iter   )
    deallocate( b_iter   )
    deallocate( foo_iter )

  end subroutine test_grid_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module test_grid

!> Driver for musica_grid tests
program test_grid_driver

  use test_grid

  implicit none

  call test_grid_t( )

end program test_grid_driver
