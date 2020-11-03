! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The micm_rate_constant_factory module

!> Builder of rate constant calculators
module micm_rate_constant_factory

  use micm_rate_constant,              only : rate_constant_t
  use micm_rate_constant_foo,          only : rate_constant_foo_t
  use micm_rate_constant_bar,          only : rate_constant_bar_t

  implicit none
  private

  public :: rate_constant_builder

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Builder of rate constant calculators
  !!
  !! At minimum, the \c config argument must include a top-level key-value
  !! pair "type" whose value is a valid rate constant type. Currently, these
  !! are:
  !! - "foo"
  !! - "bar"
  !!
  function rate_constant_builder( config ) result( new_rate_constant )

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    !> New rate constant calculator
    class(rate_constant_t), pointer :: new_rate_constant
    !> Rate constant configuration data
    type(config_t), intent(inout) :: config

    type(string_t) :: rate_constant_type
    character(len=*), parameter :: my_name = 'rate constant builder'

    new_rate_constant => null( )
    call config%get( 'type', rate_constant_type, my_name )
    rate_constant_type = rate_constant_type%to_lower( )

    if( rate_constant_type .eq. 'foo' ) then
      new_rate_constant => rate_constant_foo_t( config )
    else if( rate_constant_type .eq. 'bar' ) then
      new_rate_constant => rate_constant_bar_t( config )
    else
      call die_msg( 450768214, "Invalid rate constant type: '"//              &
                               rate_constant_type%to_char( )//"'" )
    end if

  end function rate_constant_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_rate_constant_factory
