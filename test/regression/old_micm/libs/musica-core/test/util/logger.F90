! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the musica_logger module

!> Tests for the musica_logger module
program test_util_logger

  use musica_assert
  use musica_logger

  implicit none

  call test_logger_t( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Tests the logger_t functions
  subroutine test_logger_t( )

    use musica_constants,              only : dk => musica_dk
    use musica_string,                 only : string_t

    type(logger_t) :: logger
    integer :: i_step
    real(kind=dk) :: start__s, stop__s, curr__s
    type(string_t) :: line
    type(string_t), allocatable :: split_line(:)

    ! Tests logger output for short simluations
    ! excluding tests of long simulations as they would make the test too long
    start__s = 10.0_dk
    stop__s  = 110.0_dk
    open( 12, file = "output_logger.txt", status = "replace" )
    logger = logger_t( start__s, stop__s, file_unit = 12 )
    do i_step = 1, 1000
      curr__s = start__s +                                                    &
                real( i_step, kind=dk ) / 1000.0_dk * ( stop__s - start__s )
      call logger%progress( curr__s )
    end do
    close( 12 )
    open( 12, file = "output_logger.txt", status = "old" )
    read( 12, * ) line
    call assert( 176665545, line .eq. "|  Simulation  |     Simulation     | Computation | Estimated time |" )
    read( 12, * ) line
    call assert( 901244333, line .eq. "| progress [s] | time remaining [s] |  time [ms]  | remaining [ms] |" )
    do i_step = 1, 11
      read( 12, * ) line
      split_line = line%split( " ", compress = .true. )
      call assert( 227786558, size( split_line ) .eq. 4 )
      call assert( 957629653, split_line(1) .eq. ( i_step - 1 ) * 10 )
      if( i_step .lt. 11 ) then
        call assert( 169947999, split_line(2) .eq. 99 - ( i_step - 1 ) * 10 )
      else
        call assert( 462786690, split_line(2) .eq. 0 )
      end if
    end do
    close( 12 )

  end subroutine test_logger_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_util_logger
