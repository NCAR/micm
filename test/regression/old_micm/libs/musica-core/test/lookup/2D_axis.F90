! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the musica_lookup_2D_axis module

!> Tests for the lookup_2D_axis_t type
program test_lookup_2D_axis

  use musica_assert
  use musica_lookup_2D_axis

  implicit none

  call test_lookup_2D_axis_t( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test lookup_2D_axis_t functions
  subroutine test_lookup_2D_axis_t( )

    use musica_constants,              only : dk => musica_dk
    use musica_config,                 only : config_t
    use musica_lookup_axis,            only : lookup_axis_t
    use musica_string,                 only : string_t

    type(lookup_2D_axis_t) :: axis
    type(lookup_axis_t) :: axis_a, axis_b
    type(config_t) :: config
    character(len=*), parameter :: my_name = "lookup_2D_axis_t tests"
    real(kind=dk) :: residuals(2), table_data(2,2), results
    real(kind=dk) :: table_data_1D(3,2,2), results_1D(3)
    integer :: indices(2), sizes(2)
    type(string_t) :: var_name, names(2)

    call config%empty( )
    call config%add( "file path", "lookup/axis_data.nc", my_name )
    call config%add( "axis 1",    "foo",                 my_name )
    call config%add( "axis 2",    "qux",                 my_name )
    axis = lookup_2D_axis_t( config )

    names = axis%names( )
    call assert( 704047324, names(1) .eq. "foo" )
    call assert( 928684014, names(2) .eq. "qux" )

    names = axis%dimension_names( )
    call assert( 416760063, names(1) .eq. "bar" )
    call assert( 236074545, names(2) .eq. "quux" )

    sizes = axis%sizes( )
    call assert( 350298424, sizes(1) .eq. 5 )
    call assert( 286931632, sizes(2) .eq. 4 )

    ! foo = [ 2.4, 4.3, 7.8, 8.0, 11.3 ]
    ! qux = [ 5.4, 7.8, 12.2, 20.0 ]
    indices = [ 1, 1 ]
    call axis%find( [ 5.4_dk, 16.3_dk], indices, residuals )
    call assert( 321522706, indices(1) .eq. 2 )
    call assert( 169061857, indices(2) .eq. 3 )
    call assert( 544706157, almost_equal( residuals(1), 1.1_dk / 3.5_dk ) )
    call assert( 541347384, almost_equal( residuals(2), 4.1_dk / 7.8_dk ) )

    var_name = "foo"
    axis_a = lookup_axis_t( var_name,                                         &
                            [ 2.4_dk, 4.3_dk, 7.8_dk, 8.0_dk, 11.3_dk ] )
    var_name = "qux"
    axis_b = lookup_axis_t( var_name,                                         &
                            [ 5.4_dk, 7.8_dk, 12.2_dk, 20.0_dk ] )

    axis = lookup_2D_axis_t( axis_a, axis_b )

    names = axis%names( )
    call assert( 243983804, names(1) .eq. "foo" )
    call assert( 691351650, names(2) .eq. "qux" )

    names = axis%dimension_names( )
    call assert( 521194746, names(1) .eq. "foo" )
    call assert( 968562592, names(2) .eq. "qux" )

    sizes = axis%sizes( )
    call assert( 798405688, sizes(1) .eq. 5 )
    call assert( 370474259, sizes(2) .eq. 4 )

    ! foo = [ 2.4, 4.3, 7.8, 8.0, 11.3 ]
    ! qux = [ 5.4, 7.8, 12.2, 20.0 ]
    indices = [ 1, 1 ]
    call axis%find( [ 5.4_dk, 16.3_dk], indices, residuals )
    call assert( 628248784, indices(1) .eq. 2 )
    call assert( 283572859, indices(2) .eq. 3 )
    call assert( 458091880, almost_equal( residuals(1), 1.1_dk / 3.5_dk ) )
    call assert( 352943376, almost_equal( residuals(2), 4.1_dk / 7.8_dk ) )

    ! static interpolation functions
    ! interpolate scalar
    table_data(1,:) = [ 2.0_dk, 4.0_dk ]
    table_data(2,:) = [ 6.0_dk, 6.6_dk ]
    residuals       = [ 0.0_dk, 0.0_dk ]
    call axis%interpolate( table_data, residuals, results )
    call assert( 335834391, results .eq. 2.0_dk )
    residuals       = [ 1.0_dk, 0.0_dk ]
    call axis%interpolate( table_data, residuals, results )
    call assert( 747760167, results .eq. 6.0_dk )
    residuals       = [ 0.0_dk, 1.0_dk ]
    call axis%interpolate( table_data, residuals, results )
    call assert( 577603263, results .eq. 4.0_dk )
    residuals       = [ 1.0_dk, 1.0_dk ]
    call axis%interpolate( table_data, residuals, results )
    call assert( 689921608, results .eq. 6.6_dk )
    residuals       = [ 0.5_dk, 0.0_dk ]
    call axis%interpolate( table_data, residuals, results )
    call assert( 237289455, almost_equal( results, 4.0_dk ) )
    residuals       = [ 0.0_dk, 0.5_dk ]
    call axis%interpolate( table_data, residuals, results )
    call assert( 967132550, almost_equal( results, 3.0_dk ) )
    residuals       = [ 0.5_dk, 1.0_dk ]
    call axis%interpolate( table_data, residuals, results )
    call assert( 796975646, almost_equal( results, 5.3_dk ) )
    residuals       = [ 1.0_dk, 0.5_dk ]
    call axis%interpolate( table_data, residuals, results )
    call assert( 626818742, almost_equal( results, 6.3_dk ) )

    ! interpolate array
    table_data_1D(1,1,:) = [  2.0_dk,  4.0_dk ]
    table_data_1D(1,2,:) = [  6.0_dk,  6.6_dk ]
    table_data_1D(2,1,:) = [  8.4_dk,  9.4_dk ]
    table_data_1D(2,2,:) = [  2.6_dk,  5.8_dk ]
    table_data_1D(3,1,:) = [ -4.6_dk,  1.8_dk ]
    table_data_1D(3,2,:) = [ 10.4_dk, 20.8_dk ]
    residuals            = [  0.0_dk,  0.0_dk ]
    call axis%interpolate( table_data_1D, residuals, results_1D )
    call assert( 523862895, results_1D(1) .eq.  2.0_dk )
    call assert( 295867432, results_1D(2) .eq.  8.4_dk )
    call assert( 125710528, results_1D(3) .eq. -4.6_dk )
    residuals            = [  1.0_dk,  0.0_dk ]
    call axis%interpolate( table_data_1D, residuals, results_1D )
    call assert( 208539334, results_1D(1) .eq.  6.0_dk )
    call assert( 655907180, results_1D(2) .eq.  2.6_dk )
    call assert( 485750276, results_1D(3) .eq. 10.4_dk )
    residuals            = [  0.0_dk,  1.0_dk ]
    call axis%interpolate( table_data_1D, residuals, results_1D )
    call assert( 315593372, results_1D(1) .eq.  4.0_dk )
    call assert( 145436468, results_1D(2) .eq.  9.4_dk )
    call assert( 940287963, results_1D(3) .eq.  1.8_dk )
    residuals            = [  1.0_dk,  1.0_dk ]
    call axis%interpolate( table_data_1D, residuals, results_1D )
    call assert( 770131059, results_1D(1) .eq.  6.6_dk )
    call assert( 599974155, results_1D(2) .eq.  5.8_dk )
    call assert( 429817251, results_1D(3) .eq. 20.8_dk )
    residuals            = [  0.5_dk,  0.0_dk ]
    call axis%interpolate( table_data_1D, residuals, results_1D )
    call assert( 259660347, almost_equal( results_1D(1),  4.0_dk ) )
    call assert( 989503442, almost_equal( results_1D(2),  5.5_dk ) )
    call assert( 819346538, almost_equal( results_1D(3),  2.9_dk ) )
    residuals            = [  0.0_dk,  0.5_dk ]
    call axis%interpolate( table_data_1D, residuals, results_1D )
    call assert( 649189634, almost_equal( results_1D(1),  3.0_dk ) )
    call assert( 196557481, almost_equal( results_1D(2),  8.9_dk ) )
    call assert( 991408976, almost_equal( results_1D(3), -1.4_dk ) )
    residuals            = [  0.5_dk,  1.0_dk ]
    call axis%interpolate( table_data_1D, residuals, results_1D )
    call assert( 538776823, almost_equal( results_1D(1),  5.3_dk ) )
    call assert( 986144669, almost_equal( results_1D(2),  7.6_dk ) )
    call assert( 533512516, almost_equal( results_1D(3), 11.3_dk ) )
    residuals            = [  1.0_dk,  0.5_dk ]
    call axis%interpolate( table_data_1D, residuals, results_1D )
    call assert( 363355612, almost_equal( results_1D(1),  6.3_dk ) )
    call assert( 193198708, almost_equal( results_1D(2),  4.2_dk ) )
    call assert( 988050203, almost_equal( results_1D(3), 15.6_dk ) )

  end subroutine test_lookup_2D_axis_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_lookup_2D_axis
