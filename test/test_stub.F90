!> \file
!> Stub test

!> Stub test program
program test_stub

  use musica_assert,                   only : assert

  implicit none

  call a_test( )
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Stub for future tests
  subroutine a_test( )

    call assert( 788034954, .true. )

  end subroutine a_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_stub
