! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the musica_initial_conditions module

!> Tests for the musica_initial_conditions module
program test_initial_conditions

  use musica_assert
  use musica_initial_conditions

  implicit none

  call test_initial_conditions_t( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_initial_conditions_t( )

    use musica_config,                 only : config_t
    use test_mock_domain,              only : mock_domain_t

    type(config_t) :: config
    type(mock_domain_t), pointer :: domain
    type(initial_conditions_t), pointer :: init_cond

    nullify( domain )
    nullify( init_cond )
    call config%empty( )
    domain => mock_domain_t( )
    init_cond => initial_conditions_t( config, domain )

    call assert( 445981235, associated( init_cond ) )

    deallocate( init_cond )
    deallocate( domain )

  end subroutine test_initial_conditions_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_initial_conditions
