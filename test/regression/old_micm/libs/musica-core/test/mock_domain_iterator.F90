! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Mock domain_iterator_t class for use in tests

!> Mock domain_iterator_t class for use in tests
module test_mock_domain_iterator

  use musica_domain_iterator,          only : domain_iterator_t

  implicit none
  private

  public :: mock_domain_iterator_t

  type, extends(domain_iterator_t) :: mock_domain_iterator_t
    integer :: curr_elem_ = 0
  contains
    procedure :: next
    procedure :: reset
  end type mock_domain_iterator_t

  interface mock_domain_iterator_t
    procedure :: constructor
  end interface mock_domain_iterator_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( target_domain ) result( new_iterator )

    use musica_target,                 only : target_t

    type(mock_domain_iterator_t), pointer    :: new_iterator
    class(target_t),              intent(in) :: target_domain

    allocate( new_iterator )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function next( this )

    class(mock_domain_iterator_t), intent(inout) :: this

    if( this%curr_elem_ .ge. 3 ) next = .false.
    this%curr_elem_ = this%curr_elem_ + 1
    next = .true.

  end function next

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine reset( this, parent )

    use musica_iterator,               only : iterator_t

    class(mock_domain_iterator_t), intent(inout)          :: this
    class(iterator_t),              intent(in),  optional :: parent

    this%curr_elem_ = 0

  end subroutine reset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module test_mock_domain_iterator
