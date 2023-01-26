! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Mock domain class for use in tests

!> Mock domain class for use in tests
module test_mock_domain

  use musica_domain,                   only : domain_t

  implicit none
  private

  public :: mock_domain_t

  type, extends(domain_t) :: mock_domain_t
  contains
    procedure :: type => domain_type
    procedure :: new_state
    procedure :: iterator
    procedure :: allocate_mutator
    procedure :: allocate_accessor
  end type mock_domain_t

  interface mock_domain_t
    procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( ) result( mock_domain )

    type(mock_domain_t), pointer  :: mock_domain

    allocate( mock_domain )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function domain_type( this )

    use musica_string,                 only : string_t

    type(string_t) :: domain_type
    class(mock_domain_t), intent(in) :: this

    domain_type = "mock domain"

  end function domain_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function new_state( this )

    use musica_domain_state,           only : domain_state_t
    use test_mock_domain_state,        only : mock_domain_state_t

    class(domain_state_t), pointer    :: new_state
    class(mock_domain_t),  intent(in) :: this

    new_state => mock_domain_state_t( )

  end function new_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function iterator( this, target_domain )

    use musica_domain_iterator,        only : domain_iterator_t
    use musica_target,                 only : target_t
    use test_mock_domain_iterator,     only : mock_domain_iterator_t

    class(domain_iterator_t), pointer    :: iterator
    class(mock_domain_t),     intent(in) :: this
    class(target_t),          intent(in) :: target_domain

    iterator => mock_domain_iterator_t( target_domain )

  end function iterator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function allocate_mutator( this, target_domain, data_type, property_index ) &
      result( mutator )

    use musica_constants,              only : musica_ik
    use musica_data_type,              only : data_type_t
    use musica_domain_state_mutator,   only : domain_state_mutator_t
    use musica_target,                 only : target_t
    use test_mock_domain_state_mutator, only : mock_domain_state_mutator_t

    class(domain_state_mutator_t), pointer :: mutator
    class(mock_domain_t), intent(in) :: this
    class(target_t), intent(in) :: target_domain
    type(data_type_t), intent(in) :: data_type
    integer(kind=musica_ik), intent(in) :: property_index

    mutator => mock_domain_state_mutator_t( target_domain, data_type,         &
                                            property_index )

  end function allocate_mutator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function allocate_accessor( this, target_domain, data_type, property_index )&
      result( accessor )

    use musica_constants,              only : musica_ik
    use musica_data_type,              only : data_type_t
    use musica_domain_state_accessor,  only : domain_state_accessor_t
    use musica_target,                 only : target_t
    use test_mock_domain_state_accessor, only : mock_domain_state_accessor_t

    class(domain_state_accessor_t), pointer :: accessor
    class(mock_domain_t), intent(in) :: this
    class(target_t), intent(in) :: target_domain
    type(data_type_t), intent(in) :: data_type
    integer(kind=musica_ik), intent(in) :: property_index

    accessor => mock_domain_state_accessor_t( target_domain, data_type,       &
                                              property_index )

  end function allocate_accessor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module test_mock_domain
