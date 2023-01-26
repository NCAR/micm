! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Mock domain_state_accessor_t class for use in tests

!> Mock domain_state_accessor_t class for use in tests
module test_mock_domain_state_accessor

  use musica_domain_state_accessor,    only : domain_state_accessor_t

  implicit none
  private

  public :: mock_domain_state_accessor_t

  type, extends(domain_state_accessor_t) :: mock_domain_state_accessor_t
  end type mock_domain_state_accessor_t

  interface mock_domain_state_accessor_t
    procedure :: constructor
  end interface
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( target_domain, data_type, property_index )            &
      result( new_accessor )

    use musica_constants,              only : musica_ik
    use musica_data_type,              only : data_type_t
    use musica_target,                 only : target_t

    type(mock_domain_state_accessor_t), pointer    :: new_accessor
    class(target_t),                    intent(in) :: target_domain
    type(data_type_t),                  intent(in) :: data_type
    integer(kind=musica_ik),            intent(in) :: property_index

    allocate( new_accessor )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module test_mock_domain_state_accessor
