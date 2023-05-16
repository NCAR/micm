! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the musica_property module

!> Tests for the property_t type
module test_util_property_module

  use musica_target,                   only : target_t

  implicit none

  !> Mock target
  type, extends(target_t) :: mock_target_t
  contains
    procedure :: name => target_name
    procedure :: equals_target
  end type mock_target_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the name of the target
  type(string_t) function target_name( this )

    use musica_string,                 only : string_t

    !> Target
    class(mock_target_t), intent(in) :: this

    target_name = "mock target"

  end function target_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Equality comparison
  logical function equals_target( a, b ) result( eq )

    !> Target
    class(mock_target_t), intent(in) :: a
    !> Other target
    class(target_t), intent(in) :: b

    select type( b )
    class is( mock_target_t )
      eq = .true.
    class default
      eq = .false.
    end select

  end function equals_target

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module test_util_property_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_util_property

  use musica_assert
  use musica_constants,                only : musica_dk, musica_ik, musica_lk
  use musica_property
  use test_util_property_module

  implicit none

  call test_property_t( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test property_t functionality
  subroutine test_property_t( )

    use musica_config,                 only : config_t
    use musica_data_type,              only : kBoolean, kDouble, kInteger

    character(len=*), parameter :: my_name = "property_t tests"
    type(property_t), pointer :: a, b, c, d
    type(config_t) :: prop_config
    type(mock_target_t) :: prop_target
    class(target_t), pointer :: applies_to
    real(kind=musica_dk) :: temp_double
    integer(kind=musica_ik) :: temp_int
    logical(kind=musica_lk) :: temp_bool

    ! config constructor
    call prop_config%empty( )
    call prop_config%add( "name", "prop a", my_name )
    call prop_config%add( "units", "m", my_name )
    call prop_config%add( "data type", "double", my_name )
    a => property_t( prop_config, my_name, applies_to = prop_target )

    ! standard constructor
    b => property_t( my_name,                                                 &
                     name = "my prefix%prop b",                               &
                     units = "count",                                         &
                     applies_to = prop_target,                                &
                     data_type = kInteger,                                    &
                     default_value = 7_musica_ik )

    ! constructor from existing property
    c => property_t( b, my_name,                                              &
                     name = "prop c",                                         &
                     units = "flag",                                          &
                     data_type = kBoolean,                                    &
                     default_value = .false. )

    ! assignment
    allocate( d )
    d = c

    ! equality
    call assert( 996130093, a .ne. b )
    call assert( 938291534, a .ne. c )
    call assert( 880452975, a .ne. d )
    call assert( 822614416, b .ne. c )
    call assert( 712201605, b .ne. d )
    call assert( 936838295, c .eq. d )

    ! names and prefixes
    call assert( 700219752, a%name( ) .eq. "prop a" )
    call assert( 189749040, b%name( ) .eq. "my prefix%prop b" )
    call assert( 526704075, c%name( ) .eq. "prop c" )
    call assert( 751340765, d%name( ) .eq. "prop c" )
    call assert( 418196798, a%prefix( ) .eq. "" )
    call assert( 295349839, b%prefix( ) .eq. "my prefix" )
    call assert( 914780123, c%prefix( ) .eq. "" )
    call assert( 521892063, d%prefix( ) .eq. "" )
    call assert( 746528753, a%base_name( ) .eq. "prop a" )
    call assert( 295802134, b%base_name( ) .eq. "prop b" )
    call assert( 455430424, c%base_name( ) .eq. "prop c" )
    call assert( 345017613, d%base_name( ) .eq. "prop c" )

    ! units
    call assert( 217358642, a%units( ) .eq. "m" )
    call assert( 384156773, b%units( ) .eq. "count" )
    call assert( 326318214, c%units( ) .eq. "flag" )
    call assert( 215905403, d%units( ) .eq. "flag" )

    ! data types
    call assert( 384609068, a%data_type( ) .eq. kDouble  )
    call assert( 546142892, b%data_type( ) .eq. kInteger )
    call assert( 370721681, c%data_type( ) .eq. kBoolean )
    call assert( 425201467, d%data_type( ) .eq. kBoolean )

    ! applies to
    applies_to => a%applies_to( )
    call assert( 163618274, applies_to .eq. prop_target )
    deallocate( applies_to )
    applies_to => b%applies_to( )
    call assert( 837528344, applies_to .eq. prop_target )
    deallocate( applies_to )
    applies_to => c%applies_to( )
    call assert( 949846689, applies_to .eq. prop_target )
    deallocate( applies_to )
    applies_to => d%applies_to( )
    call assert( 162165035, applies_to .eq. prop_target )
    deallocate( applies_to )

    ! default value
    temp_double = 123.654_musica_dk
    call a%get_default( temp_double ) ! no default for a, so should be unchanged
    call assert( 941675904, temp_double .eq. 123.654_musica_dk )
    temp_int = 24_musica_ik
    call b%get_default( temp_int )
    call assert( 814016933, temp_int .eq. 7_musica_ik )
    temp_bool = .true.
    call c%get_default( temp_bool )
    call assert( 412505793, .not. temp_bool )
    temp_bool = .true.
    call d%get_default( temp_bool )
    call assert( 466985579, .not. temp_bool )

    ! defined by
    call assert( 133841612, a%defined_by( ) .eq. my_name )
    call assert( 353213995, b%defined_by( ) .eq. my_name )
    call assert( 465532340, c%defined_by( ) .eq. my_name )
    call assert( 577850685, d%defined_by( ) .eq. my_name )

    ! bug fix test
    deallocate( a )
    a => property_t( my_name, name = "foo", units = "bar" )
    call assert( 975856464, a%name( )  .eq. "foo" )
    call assert( 582968404, a%units( ) .eq. "bar" )

    ! free memory
    deallocate( a )
    deallocate( b )
    deallocate( c )
    deallocate( d )

  end subroutine test_property_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_util_property
