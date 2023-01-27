! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the musica_lookup_axis module

!> Tests for the lookup_axis_t type
program test_lookup_axis

  use musica_assert
  use musica_lookup_axis

  implicit none

  character(len=256) :: failure_test_type

  if( command_argument_count( ) .eq. 0 ) then
    call test_lookup_axis_t( )
  else if( command_argument_count( ) .eq. 1 ) then
    call get_command_argument( 1, failure_test_type )
    call failure_test( failure_test_type )
  else
    call die( 557641360 )
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test lookup_axis_t functions
  subroutine test_lookup_axis_t( )

    use musica_constants,              only : dk => musica_dk
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    type(lookup_axis_t) :: axis
    type(config_t) :: config
    character(len=*), parameter :: my_name = "lookup_axis_t tests"
    real(kind=dk) :: residual, val
    real(kind=dk) :: vals(2)
    type(string_t) :: var_name
    integer :: index

    call config%empty( )
    call config%add( "file path",     "axis_data.nc", my_name )
    call config%add( "variable name", "foo",          my_name )
    axis = lookup_axis_t( config )

    call assert( 681915961, axis%name( )           .eq. "foo" )
    call assert( 382646847, axis%dimension_name( ) .eq. "bar" )
    call assert( 831920227, axis%size( )           .eq. 5     )

    ! data: 2.4, 4.3, 7.8, 8.0, 11.3
    index = 1
    call axis%find( 5.4_dk, index, residual )
    call assert( 257345705, index .eq. 2 )
    call assert( 185619759, almost_equal( residual, 1.1_dk / 3.5_dk ) )
    index = 2
    call axis%find( 5.4_dk, index, residual )
    call assert( 211510996, index .eq. 2 )
    call assert( 658878842, almost_equal( residual, 1.1_dk / 3.5_dk ) )
    index = 5
    call axis%find( 5.4_dk, index, residual )
    call assert( 436147686, index .eq. 2 )
    call assert( 265990782, almost_equal( residual, 1.1_dk / 3.5_dk ) )
    call axis%find( 1.2_dk, index, residual )
    call assert( 800238031, index .eq. 1 )
    call assert( 737135165, residual .eq. 0.0_dk )
    call axis%find( 12.4_dk, index, residual )
    call assert( 451301143, index .eq. 4 )
    call assert( 175995735, residual .eq. 1.0_dk )

    var_name = "foo"
    axis = lookup_axis_t( var_name,                                           &
                          [ 2.4_dk, 4.3_dk, 7.8_dk, 8.0_dk, 11.3_dk ] )

    call assert( 416260174, axis%name( )           .eq. "foo" )
    call assert( 246103270, axis%dimension_name( ) .eq. "foo" )
    call assert( 693471116, axis%size( )           .eq. 5     )

    ! data: 2.4, 4.3, 7.8, 8.0, 11.3
    index = 1
    call axis%find( 5.4_dk, index, residual )
    call assert( 523314212, index .eq. 2 )
    call assert( 418165708, almost_equal( residual, 1.1_dk / 3.5_dk ) )
    call axis%find( 1.2_dk, index, residual )
    call assert( 865533554, index .eq. 1 )
    call assert( 412901401, residual .eq. 0.0_dk )
    call axis%find( 12.4_dk, index, residual )
    call assert( 860269247, index .eq. 4 )
    call assert( 407637094, residual .eq. 1.0_dk )

    ! static interpolation functions
    call axis%interpolate( 1.2_dk, 3.2_dk, 0.5_dk, val )
    call assert( 723700073, almost_equal( val, 2.2_dk ) )
    call axis%interpolate( [ 3.1_dk, 4.3_dk ], [ 7.1_dk, 4.7_dk ], 0.75_dk, vals )
    call assert( 831206406, almost_equal( vals(1), 6.1_dk ) )
    call assert( 369951173, almost_equal( vals(2), 4.6_dk ) )

  end subroutine test_lookup_axis_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Failure tests for lookup_axis_t class
  subroutine failure_test( test_type )

    character(len=*), intent(in) :: test_type

    if( test_type .eq. "646167484" ) then
      call failure_test_646167484( )
    else
      call die( 549470575 )
    end if

  end subroutine failure_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine failure_test_646167484( )

    use musica_config,                 only : config_t

    character(len=*), parameter :: my_name = "lookup_axis_t failure test"
    type(lookup_axis_t) :: axis
    type(config_t) :: config

    call config%empty( )
    call config%add( "file path",     "axis_data.nc", my_name )
    call config%add( "variable name", "bad_foo",      my_name )
    axis = lookup_axis_t( config )

  end subroutine failure_test_646167484

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_lookup_axis
