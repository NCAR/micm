! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Mock domain_state_t class for use in tests

!> Mock domain_state_t class for use in tests
module test_mock_domain_state

  use musica_domain_state,             only : domain_state_t

  implicit none
  private

  public :: mock_domain_state_t

  type, extends(domain_state_t) :: mock_domain_state_t
  contains
    procedure :: get_boolean
    procedure :: get_double
    procedure :: get_float
    procedure :: get_integer
    procedure :: update_boolean
    procedure :: update_double
    procedure :: update_float
    procedure :: update_integer
  end type mock_domain_state_t

  interface mock_domain_state_t
    procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( ) result( new_state )

    type(mock_domain_state_t), pointer :: new_state

    allocate( new_state )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_boolean( this, iterator, accessor, state_value )

    use musica_constants,              only : musica_lk
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_domain_state_accessor,  only : domain_state_accessor_t

    class(mock_domain_state_t),     intent(in)  :: this
    class(domain_iterator_t),       intent(in)  :: iterator
    class(domain_state_accessor_t), intent(in)  :: accessor
    logical(kind=musica_lk),        intent(out) :: state_value

    state_value = .false.

  end subroutine get_boolean

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_double( this, iterator, accessor, state_value )

    use musica_constants,              only : musica_dk
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_domain_state_accessor,  only : domain_state_accessor_t

    class(mock_domain_state_t),     intent(in)  :: this
    class(domain_iterator_t),       intent(in)  :: iterator
    class(domain_state_accessor_t), intent(in)  :: accessor
    real(kind=musica_dk),           intent(out) :: state_value

    state_value = 14.56_musica_dk

  end subroutine get_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_float( this, iterator, accessor, state_value )

    use musica_constants,              only : musica_rk
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_domain_state_accessor,  only : domain_state_accessor_t

    class(mock_domain_state_t),     intent(in)  :: this
    class(domain_iterator_t),       intent(in)  :: iterator
    class(domain_state_accessor_t), intent(in)  :: accessor
    real(kind=musica_rk),           intent(out) :: state_value

    state_value = 52.34_musica_rk

  end subroutine get_float

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_integer( this, iterator, accessor, state_value )

    use musica_constants,              only : musica_ik
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_domain_state_accessor,  only : domain_state_accessor_t

    class(mock_domain_state_t),     intent(in)  :: this
    class(domain_iterator_t),       intent(in)  :: iterator
    class(domain_state_accessor_t), intent(in)  :: accessor
    integer(kind=musica_ik),        intent(out) :: state_value

    state_value = 42_musica_ik

  end subroutine get_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update_boolean( this, iterator, mutator, state_value )

    use musica_constants,              only : musica_lk
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_domain_state_mutator,   only : domain_state_mutator_t

    class(mock_domain_state_t),     intent(inout) :: this
    class(domain_iterator_t),       intent(in)    :: iterator
    class(domain_state_mutator_t),  intent(in)    :: mutator
    logical(kind=musica_lk),        intent(in)    :: state_value

  end subroutine update_boolean

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update_double( this, iterator, mutator, state_value )

    use musica_constants,              only : musica_dk
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_domain_state_mutator,   only : domain_state_mutator_t

    class(mock_domain_state_t),     intent(inout) :: this
    class(domain_iterator_t),       intent(in)    :: iterator
    class(domain_state_mutator_t),  intent(in)    :: mutator
    real(kind=musica_dk),           intent(in)    :: state_value

  end subroutine update_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update_float( this, iterator, mutator, state_value )

    use musica_constants,              only : musica_rk
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_domain_state_mutator,   only : domain_state_mutator_t

    class(mock_domain_state_t),     intent(inout) :: this
    class(domain_iterator_t),       intent(in)    :: iterator
    class(domain_state_mutator_t),  intent(in)    :: mutator
    real(kind=musica_rk),           intent(in)    :: state_value

  end subroutine update_float

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update_integer( this, iterator, mutator, state_value )

    use musica_constants,              only : musica_ik
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_domain_state_mutator,   only : domain_state_mutator_t

    class(mock_domain_state_t),     intent(inout) :: this
    class(domain_iterator_t),       intent(in)    :: iterator
    class(domain_state_mutator_t),  intent(in)    :: mutator
    integer(kind=musica_ik),        intent(in)    :: state_value

  end subroutine update_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module test_mock_domain_state
