! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the musica_domain_cell module

!> Test module for the musica_domain_cell module
program test_domain_cell

  use musica_assert
  use musica_domain_cell

  implicit none

  call test_domain_cell_t( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test domain_t functionality
  subroutine test_domain_cell_t( )

    use musica_config,                 only : config_t
    use musica_constants,              only : dk => musica_dk,                &
                                              ik => musica_ik,                &
                                              lk => musica_lk,                &
                                              rk => musica_rk
    use musica_data_type,              only : kBoolean, kDouble, kFloat,      &
                                              kInteger
    use musica_domain,                 only : domain_t
    use musica_domain_state,           only : domain_state_t
    use musica_domain_target_cells,    only : domain_target_cells_t
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_domain_state_mutator,   only : domain_state_mutator_t,         &
                                              domain_state_mutator_ptr
    use musica_domain_state_accessor,  only : domain_state_accessor_t,        &
                                              domain_state_accessor_ptr
    use musica_property,               only : property_ptr, property_t
    use musica_property_set,           only : property_set_t
    use musica_string,                 only : string_t, to_char

    class(domain_t), pointer :: domain
    type(config_t) :: config
    class(domain_state_t), pointer :: state, state2
    class(domain_iterator_t), pointer :: iter
    type(property_ptr) :: prop_ia, prop_ib, prop_fa, prop_fb
    type(property_ptr) :: prop_da, prop_db, prop_ba, prop_bb
    type(property_ptr) :: prop_base, prop_temp
    type(property_set_t) :: prop_set
    type(domain_target_cells_t) :: all_cells
    type(domain_state_mutator_ptr) :: mut_ia, mut_ib, mut_fa, mut_fb
    type(domain_state_mutator_ptr) :: mut_da, mut_db, mut_ba, mut_bb
    type(domain_state_mutator_ptr) :: mut_ia2, mut_fa2, mut_da2, mut_ba2
    type(domain_state_accessor_ptr) :: acc_ia, acc_ib, acc_fa, acc_fb
    type(domain_state_accessor_ptr) :: acc_da, acc_db, acc_ba, acc_bb
    type(domain_state_accessor_ptr) :: acc_ia2, acc_fa2, acc_da2, acc_ba2
    class(domain_state_mutator_ptr), pointer :: mut_set(:), mut_set2(:)
    class(domain_state_accessor_ptr), pointer :: acc_set(:), acc_set2(:)
    character(len=*), parameter :: my_name = "domain_cell_t tests"
    integer :: i
    type(string_t) :: temp_str
    integer(kind=ik) :: temp_int
    real(kind=rk) :: temp_float
    real(kind=dk) :: temp_double
    logical(kind=lk) :: temp_bool

    ! constructor
    ! (no config information is needed)
    call config%empty( )
    domain => domain_cell_t( config )

    ! properties to register
    prop_ia%val_ => property_t( my_name, name = "ia", units = "count",        &
                                applies_to = all_cells, data_type = kInteger, &
                                default_value = 15_ik )
    prop_ib%val_ => property_t( my_name, name = "ib", units = "count",        &
                                applies_to = all_cells, data_type = kInteger )
    prop_fa%val_ => property_t( my_name, name = "fa", units = "s",            &
                                applies_to = all_cells, data_type = kFloat,   &
                                default_value = 52.3_rk )
    prop_fb%val_ => property_t( my_name, name = "fb", units = "s",            &
                                applies_to = all_cells, data_type = kFloat )
    prop_da%val_ => property_t( my_name, name = "da", units = "K",            &
                                applies_to = all_cells, data_type = kDouble,  &
                                default_value = 93.4_dk )
    prop_db%val_ => property_t( my_name, name = "db", units = "K",            &
                                applies_to = all_cells, data_type = kDouble )
    prop_ba%val_ => property_t( my_name, name = "ba", units = "unitless",     &
                                applies_to = all_cells, data_type = kBoolean, &
                                default_value = .false. )
    prop_bb%val_ => property_t( my_name, name = "bb", units = "unitless",     &
                                applies_to = all_cells, data_type = kBoolean )

    prop_base%val_ => property_t( my_name, name = "var 1", units = "Pa",      &
                                  applies_to = all_cells, data_type = kDouble,&
                                  default_value = 0.0_dk )
    prop_set = property_set_t( )
    do i = 1, 3
      temp_str = "var "//to_char( i )
      prop_temp%val_ => property_t( prop_base%val_, my_name,                  &
                                    name = temp_str%to_char( ) )
      call assert( 293668162, associated( prop_temp%val_ ) )
      call prop_set%add( prop_temp%val_ )
      deallocate( prop_temp%val_ )
    end do

    ! register
    call domain%register( prop_ia%val_ )
    call domain%register( prop_ib%val_ )
    call domain%register( prop_fa%val_ )
    call domain%register( prop_fb%val_ )
    call domain%register( prop_da%val_ )
    call domain%register( prop_db%val_ )
    call domain%register( prop_ba%val_ )
    call domain%register( prop_bb%val_ )
    call domain%register( "my set", prop_set )

    ! lock the domain configuration
    call domain%lock( )

    ! get mutators
    mut_ia%val_   => domain%mutator( prop_ia%val_ )
    mut_ia2%val_  => domain%mutator( prop_ia%val_ )
    mut_ib%val_   => domain%mutator( prop_ib%val_ )
    mut_fa%val_   => domain%mutator( prop_fa%val_ )
    mut_fa2%val_  => domain%mutator( prop_fa%val_ )
    mut_fb%val_   => domain%mutator( prop_fb%val_ )
    mut_da%val_   => domain%mutator( prop_da%val_ )
    mut_da2%val_  => domain%mutator( prop_da%val_ )
    mut_db%val_   => domain%mutator( prop_db%val_ )
    mut_ba%val_   => domain%mutator( prop_ba%val_ )
    mut_ba2%val_  => domain%mutator( prop_ba%val_ )
    mut_bb%val_   => domain%mutator( prop_bb%val_ )
    mut_set  => domain%mutator_set( "my set", "Pa", kDouble, all_cells,       &
                                    my_name )
    mut_set2 => domain%mutator_set( "my set", "Pa", kDouble, all_cells,       &
                                    my_name )

    ! get accessors
    acc_ia%val_   => domain%accessor( prop_ia%val_ )
    acc_ia2%val_  => domain%accessor( prop_ia%val_ )
    acc_ib%val_   => domain%accessor( prop_ib%val_ )
    acc_fa%val_   => domain%accessor( prop_fa%val_ )
    acc_fa2%val_  => domain%accessor( prop_fa%val_ )
    acc_fb%val_   => domain%accessor( prop_fb%val_ )
    acc_da%val_   => domain%accessor( prop_da%val_ )
    acc_da2%val_  => domain%accessor( prop_da%val_ )
    acc_db%val_   => domain%accessor( prop_db%val_ )
    acc_ba%val_   => domain%accessor( prop_ba%val_ )
    acc_ba2%val_  => domain%accessor( prop_ba%val_ )
    acc_bb%val_   => domain%accessor( prop_bb%val_ )
    acc_set  => domain%accessor_set( "my set", "Pa", kDouble, all_cells,       &
                                    my_name )
    acc_set2 => domain%accessor_set( "my set", "Pa", kDouble, all_cells,       &
                                    my_name )

    ! get property
    prop_temp%val_ => mut_ia%val_%property( )
    call assert( 781455245, prop_temp%val_%name( ) .eq. "ia" )
    deallocate( prop_temp%val_ )
    prop_temp%val_ => mut_ib%val_%property( )
    call assert( 734432423, prop_temp%val_%name( ) .eq. "ib" )
    deallocate( prop_temp%val_ )
    prop_temp%val_ => mut_fa%val_%property( )
    call assert( 794176516, prop_temp%val_%name( ) .eq. "fa" )
    deallocate( prop_temp%val_ )
    prop_temp%val_ => mut_fb%val_%property( )
    call assert( 624019612, prop_temp%val_%name( ) .eq. "fb" )
    deallocate( prop_temp%val_ )
    prop_temp%val_ => mut_da%val_%property( )
    call assert( 118813207, prop_temp%val_%name( ) .eq. "da" )
    deallocate( prop_temp%val_ )
    prop_temp%val_ => mut_db%val_%property( )
    call assert( 231131552, prop_temp%val_%name( ) .eq. "db" )
    deallocate( prop_temp%val_ )
    prop_temp%val_ => mut_ba%val_%property( )
    call assert( 508342494, prop_temp%val_%name( ) .eq. "ba" )
    deallocate( prop_temp%val_ )
    prop_temp%val_ => mut_bb%val_%property( )
    call assert( 620660839, prop_temp%val_%name( ) .eq. "bb" )
    deallocate( prop_temp%val_ )
    prop_temp%val_ => mut_set(1)%val_%property( )
    call assert( 450503935, prop_temp%val_%name( ) .eq. "my set%var 1" )
    deallocate( prop_temp%val_ )
    prop_temp%val_ => mut_set(2)%val_%property( )
    call assert( 562822280, prop_temp%val_%name( ) .eq. "my set%var 2" )
    deallocate( prop_temp%val_ )
    prop_temp%val_ => mut_set(3)%val_%property( )
    call assert( 675140625, prop_temp%val_%name( ) .eq. "my set%var 3" )
    deallocate( prop_temp%val_ )

    ! get domain state objects
    state  => domain%new_state( )
    state2 => domain%new_state( )

    ! check units
    temp_str = domain%units( "ia" )
    call assert( 348289395, temp_str .eq. "count" )
    temp_str = domain%units( "ib" )
    call assert( 844872720, temp_str .eq. "count" )
    temp_str = domain%units( "fa" )
    call assert( 355718873, temp_str .eq. "s" )
    temp_str = domain%units( "fb" )
    call assert( 468037218, temp_str .eq. "s" )
    temp_str = domain%units( "da" )
    call assert( 297880314, temp_str .eq. "K" )
    temp_str = domain%units( "db" )
    call assert( 127723410, temp_str .eq. "K" )
    temp_str = domain%units( "ba" )
    call assert( 857566505, temp_str .eq. "unitless" )
    temp_str = domain%units( "bb" )
    call assert( 687409601, temp_str .eq. "unitless" )
    temp_str = domain%units( "my set%var 1" )
    call assert( 103047772, temp_str .eq. "Pa" )
    temp_str = domain%units( "my set%var 2" )
    call assert( 382164248, temp_str .eq. "Pa" )
    temp_str = domain%units( "my set%var 3" )
    call assert( 724383590, temp_str .eq. "Pa" )

    ! check for variables and flags
    call assert( 140444681, domain%is_registered( "ia" ) )
    call assert( 310601585, domain%is_registered( "ib" ) )
    call assert( 133274840, domain%is_registered( "fa" ) )
    call assert( 585906993, domain%is_registered( "fb" ) )
    call assert( 392097426, domain%is_registered( "da" ) )
    call assert( 169366270, domain%is_registered( "db" ) )
    call assert( 895850592, domain%is_registered( "ba" ) )
    call assert( 390644187, domain%is_registered( "bb" ) )
    call assert( 846635113, domain%is_registered( "my set%var 1" ) )
    call assert( 341428708, domain%is_registered( "my set%var 2" ) )
    call assert( 388738653, domain%is_registered( "my set%var 3" ) )
    call assert( 501056998, .not. domain%is_registered( "not there" ) )

    ! get an iterator
    iter => domain%iterator( all_cells )

    ! default values for both states
    do while( iter%next( ) )
      call state%get( iter, acc_ia%val_, temp_int )
      call assert( 203122738, temp_int .eq. 15_ik )
      call state%get( iter, acc_ia2%val_, temp_int )
      call assert( 734033211, temp_int .eq. 15_ik )
      call state%get( iter, acc_fa%val_, temp_float )
      call assert( 866504194, temp_float .eq. 52.3_rk )
      call state%get( iter, acc_fa2%val_, temp_float )
      call assert( 447457882, temp_float .eq. 52.3_rk )
      call state%get( iter, acc_da%val_, temp_double )
      call assert( 600820921, temp_double .eq. 93.4_dk )
      call state%get( iter, acc_da2%val_, temp_double )
      call assert( 430664017, temp_double .eq. 93.4_dk )
      call state%get( iter, acc_ba%val_, temp_bool )
      call assert( 260507113, temp_bool .eqv. .false. )
      call state%get( iter, acc_ba2%val_, temp_bool )
      call assert( 437833858, temp_bool .eqv. .false. )
      call state%get( iter, acc_set(1)%val_, temp_double )
      call assert( 308723537, temp_double .eq. 0.0_dk )
      call state%get( iter, acc_set(2)%val_, temp_double )
      call assert( 303459230, temp_double .eq. 0.0_dk )
      call state%get( iter, acc_set(3)%val_, temp_double )
      call assert( 363203323, temp_double .eq. 0.0_dk )
      call state%get( iter, acc_set2(1)%val_, temp_double )
      call assert( 917623318, temp_double .eq. 0.0_dk )
      call state%get( iter, acc_set2(2)%val_, temp_double )
      call assert( 359842661, temp_double .eq. 0.0_dk )
      call state%get( iter, acc_set2(3)%val_, temp_double )
      call assert( 189685757, temp_double .eq. 0.0_dk )
      call state2%get( iter, acc_ia%val_, temp_int )
      call assert( 169080824, temp_int .eq. 15_ik )
      call state2%get( iter, acc_ia2%val_, temp_int )
      call assert( 846349667, temp_int .eq. 15_ik )
      call state2%get( iter, acc_fa%val_, temp_float )
      call assert( 958668012, temp_float .eq. 52.3_rk )
      call state2%get( iter, acc_fa2%val_, temp_float )
      call assert( 170986358, temp_float .eq. 52.3_rk )
      call state2%get( iter, acc_da%val_, temp_double )
      call assert( 900829453, temp_double .eq. 93.4_dk )
      call state2%get( iter, acc_da2%val_, temp_double )
      call assert( 730672549, temp_double .eq. 93.4_dk )
      call state2%get( iter, acc_ba%val_, temp_bool )
      call assert( 560515645, temp_bool .eqv. .false. )
      call state2%get( iter, acc_ba2%val_, temp_bool )
      call assert( 390358741, temp_bool .eqv. .false. )
      call state2%get( iter, acc_set(1)%val_, temp_double )
      call assert( 785152335, temp_double .eq. 0.0_dk )
      call state2%get( iter, acc_set(2)%val_, temp_double )
      call assert( 897470680, temp_double .eq. 0.0_dk )
      call state2%get( iter, acc_set(3)%val_, temp_double )
      call assert( 727313776, temp_double .eq. 0.0_dk )
      call state2%get( iter, acc_set2(1)%val_, temp_double )
      call assert( 839632121, temp_double .eq. 0.0_dk )
      call state2%get( iter, acc_set2(2)%val_, temp_double )
      call assert( 669475217, temp_double .eq. 0.0_dk )
      call state2%get( iter, acc_set2(3)%val_, temp_double )
      call assert( 499318313, temp_double .eq. 0.0_dk )
    end do

    ! update state 2 with first set of mutators
    call iter%reset( )
    do while( iter%next( ) )
      call state2%update( iter, mut_ib%val_, 65_ik )
      call state2%update( iter, mut_fb%val_, 798.4_rk )
      call state2%update( iter, mut_db%val_, 834.9_dk )
      call state2%update( iter, mut_bb%val_, .true. )
      call state2%update( iter, mut_set(1)%val_, 458.1_dk )
      call state2%update( iter, mut_set(2)%val_, 65.23_dk )
      call state2%update( iter, mut_set(2)%val_, 95.10_dk )
    end do

    ! default values in state and updated values in state 2
    do while( iter%next( ) )
      call state%get( iter, acc_ia%val_, temp_int )
      call assert( 347434110, temp_int .eq. 15_ik )
      call state%get( iter, acc_fa%val_, temp_float )
      call assert( 689653452, temp_float .eq. 52.3_rk )
      call state%get( iter, acc_da%val_, temp_double )
      call assert( 519496548, temp_double .eq. 93.4_dk )
      call state%get( iter, acc_ba%val_, temp_bool )
      call assert( 631814893, temp_bool .eqv. .false. )
      call state%get( iter, acc_set(1)%val_, temp_double )
      call assert( 744133238, temp_double .eq. 0.0_dk )
      call state%get( iter, acc_set(2)%val_, temp_double )
      call assert( 573976334, temp_double .eq. 0.0_dk )
      call state%get( iter, acc_set(3)%val_, temp_double )
      call assert( 403819430, temp_double .eq. 0.0_dk )
      call state2%get( iter, acc_ia%val_, temp_int )
      call assert( 233662526, temp_int .eq. 15_ik )
      call state2%get( iter, acc_ib%val_, temp_int )
      call assert( 963505621, temp_int .eq. 65_ik )
      call state2%get( iter, acc_fa%val_, temp_float )
      call assert( 175823967, temp_float .eq. 52.3_rk )
      call state2%get( iter, acc_fb%val_, temp_float )
      call assert( 905667062, temp_float .eq. 798.4_rk )
      call state2%get( iter, acc_da%val_, temp_double )
      call assert( 735510158, temp_double .eq. 93.4_dk )
      call state2%get( iter, acc_db%val_, temp_double )
      call assert( 847828503, temp_double .eq. 834.9_dk )
      call state2%get( iter, acc_ba%val_, temp_bool )
      call assert( 677671599, temp_bool .eqv. .false. )
      call state2%get( iter, acc_bb%val_, temp_bool )
      call assert( 507514695, temp_bool .eqv. .true. )
      call state2%get( iter, acc_set(1)%val_, temp_double )
      call assert( 402366191, temp_double .eq. 458.1_dk )
      call state2%get( iter, acc_set(2)%val_, temp_double )
      call assert( 232209287, temp_double .eq. 65.23_dk )
      call state2%get( iter, acc_set(3)%val_, temp_double )
      call assert( 962052382, temp_double .eq. 95.10_dk )
    end do

    ! clean up memory
    deallocate( domain )
    deallocate( state )
    deallocate( state2 )
    deallocate( iter )
    deallocate( prop_ia%val_ )
    deallocate( prop_ib%val_ )
    deallocate( prop_fa%val_ )
    deallocate( prop_fb%val_ )
    deallocate( prop_da%val_ )
    deallocate( prop_db%val_ )
    deallocate( prop_ba%val_ )
    deallocate( prop_bb%val_ )
    deallocate( mut_ia%val_ )
    deallocate( mut_ib%val_ )
    deallocate( mut_fa%val_ )
    deallocate( mut_fb%val_ )
    deallocate( mut_da%val_ )
    deallocate( mut_db%val_ )
    deallocate( mut_ba%val_ )
    deallocate( mut_bb%val_ )
    deallocate( acc_ia%val_ )
    deallocate( acc_ib%val_ )
    deallocate( acc_fa%val_ )
    deallocate( acc_fb%val_ )
    deallocate( acc_da%val_ )
    deallocate( acc_db%val_ )
    deallocate( acc_ba%val_ )
    deallocate( acc_bb%val_ )
    do i = 1, size( mut_set )
      deallocate( mut_set( i )%val_ )
    end do
    deallocate( mut_set )
    do i = 1, size( mut_set2 )
      deallocate( mut_set2( i )%val_ )
    end do
    deallocate( mut_set2 )
    do i = 1, size( acc_set )
      deallocate( acc_set( i )%val_ )
    end do
    deallocate( acc_set )
    do i = 1, size( acc_set2 )
      deallocate( acc_set2( i )%val_ )
    end do
    deallocate( acc_set2 )

  end subroutine test_domain_cell_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_domain_cell
