! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The test_io_netcdf program

!> Tests for the io_netcdf_t type
program test_io_netcdf

  use musica_assert
  use musica_io_netcdf

  implicit none

  call test_io_netcdf_t( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Tests for the io_netcdf_t type
  subroutine test_io_netcdf_t( )

    use musica_constants,              only : dk => musica_dk
    use musica_io,                     only : io_t
    use musica_string,                 only : string_t

    character(len=*), parameter :: my_name = "io_netcdf_t tests"
    class(io_t), pointer :: my_file
    type(string_t) :: file_name, var_name
    real(kind=dk)               :: real0D
    real(kind=dk),  allocatable :: real1D(:)
    real(kind=dk),  allocatable :: real2D(:,:)
    real(kind=dk),  allocatable :: real3D(:,:,:)
    real(kind=dk),  allocatable :: real4D(:,:,:,:)
    integer                     :: int0D
    integer,        allocatable :: int1D(:)
    type(string_t), allocatable :: dim_names(:)

    file_name = "data/io_netcdf_test_data.nc"
    my_file => io_netcdf_t( file_name )
    call assert(312726817, associated( my_file ) )

    ! check units
    var_name = "foo"
    call assert( 600672797, my_file%variable_units( var_name, my_name )       &
                            .eq. "foobits" )

    ! scalar real
    var_name = "qux"
    call my_file%read( var_name, real0D, my_name )
    call assert( 322376438, real0D .eq. 92.37_dk )

    ! scalar int
    var_name = "quux"
    call my_file%read( var_name, int0D, my_name )
    call assert( 330171719, int0D .eq. 7 )

    ! unallocated 1D real
    var_name = "foo"
    call my_file%read( var_name, real1D, my_name )
    call assert( 441036825, allocated( real1D ) )
    call assert( 442942359, size( real1D ) .eq. 4 )
    call assert( 567529680, almost_equal( real1D( 1 ), 15.32_dk ) )
    call assert( 846646156, almost_equal( real1D( 2 ), 3.14_dk  ) )
    call assert( 394014003, almost_equal( real1D( 3 ), 26.71_dk ) )
    call assert( 558906600, almost_equal( real1D( 4 ), 19.34_dk ) )
    deallocate( real1D )

    ! pre-allocatted 1D real
    allocate( real1D( 3 ) )
    var_name = "bar"
    call my_file%read( var_name, real1D, my_name )
    call assert( 329457898, allocated( real1D ) )
    call assert( 889144089, size( real1D ) .eq. 3 )
    call assert( 155942221, almost_equal( real1D( 1 ), 51.43_dk  ) )
    call assert( 885785316, almost_equal( real1D( 2 ), 123.01_dk ) )
    call assert( 150677914, almost_equal( real1D( 3 ), 32.61_dk  ) )
    deallocate( real1D )

    ! unallocated 1D int
    var_name = "quuz"
    call my_file%read( var_name, int1D, my_name )
    call assert( 595878700, allocated( int1D ) )
    call assert( 820515390, size( int1D ) .eq. 4 )
    call assert( 145152081, int1D( 1 ) .eq. 9  )
    call assert( 199631867, int1D( 2 ) .eq. 3  )
    call assert( 364524464, int1D( 3 ) .eq. 12 )
    call assert( 811892310, int1D( 4 ) .eq. 1  )
    deallocate( int1D )

    ! pre-allocated 1D int
    allocate( int1D( 4 ) )
    var_name = "quuz"
    call my_file%read( var_name, int1D, my_name )
    call assert( 917493109, allocated( int1D ) )
    call assert( 182385707, size( int1D ) .eq. 4 )
    call assert( 977237202, int1D( 1 ) .eq. 9  )
    call assert( 807080298, int1D( 2 ) .eq. 3  )
    call assert( 971972895, int1D( 3 ) .eq. 12 )
    call assert( 519340742, int1D( 4 ) .eq. 1  )
    deallocate( int1D )

    ! unallocated 2D real
    var_name = "baz"
    call my_file%read( var_name, real2D, my_name )
    call assert( 910775563, allocated( real2D ) )
    call assert( 517887503, size( real2D, 1 ) .eq. 3 )
    call assert( 454784637, size( real2D, 2 ) .eq. 4 )
    call assert( 961896576, almost_equal( real2D( 1, 1 ), 31.2_dk      ) )
    call assert( 337654280, almost_equal( real2D( 2, 1 ), 41.3_dk      ) )
    call assert( 785022126, almost_equal( real2D( 3, 1 ), 623.34_dk    ) )
    call assert( 332389973, almost_equal( real2D( 1, 2 ), 124.24_dk    ) )
    call assert( 227241469, almost_equal( real2D( 2, 2 ), 1592.3_dk    ) )
    call assert( 674609315, almost_equal( real2D( 3, 2 ), 42.53_dk     ) )
    call assert( 221977162, almost_equal( real2D( 1, 3 ), 1.3e-7_dk    ) )
    call assert( 669345008, almost_equal( real2D( 2, 3 ), -31.6_dk     ) )
    call assert( 499188104, almost_equal( real2D( 3, 3 ), 82.3_dk      ) )
    call assert( 111564351, almost_equal( real2D( 1, 4 ), 51.64_dk     ) )
    call assert( 558932197, almost_equal( real2D( 2, 4 ), -61.7_dk     ) )
    call assert( 106300044, almost_equal( real2D( 3, 4 ), -423000.0_dk ) )
    deallocate( real2D )

    ! pre-allocated 2D real
    allocate( real2D( 3, 4 ) )
    call my_file%read( var_name, real2D, my_name )
    call assert( 301447195, allocated( real2D ) )
    call assert( 748815041, size( real2D, 1 ) .eq. 3 )
    call assert( 861133386, size( real2D, 2 ) .eq. 4 )
    call assert( 408501233, almost_equal( real2D( 1, 1 ), 31.2_dk      ) )
    call assert( 303352729, almost_equal( real2D( 2, 1 ), 41.3_dk      ) )
    call assert( 133195825, almost_equal( real2D( 3, 1 ), 623.34_dk    ) )
    call assert( 863038920, almost_equal( real2D( 1, 2 ), 124.24_dk    ) )
    call assert( 410406767, almost_equal( real2D( 2, 2 ), 1592.3_dk    ) )
    call assert( 857774613, almost_equal( real2D( 3, 2 ), 42.53_dk     ) )
    call assert( 405142460, almost_equal( real2D( 1, 3 ), 1.3e-7_dk    ) )
    call assert( 234985556, almost_equal( real2D( 2, 3 ), -31.6_dk     ) )
    call assert( 129837052, almost_equal( real2D( 3, 3 ), 82.3_dk      ) )
    call assert( 577204898, almost_equal( real2D( 1, 4 ), 51.64_dk     ) )
    call assert( 742097495, almost_equal( real2D( 2, 4 ), -61.7_dk     ) )
    call assert( 854415840, almost_equal( real2D( 3, 4 ), -423000.0_dk ) )
    deallocate( real2D )

    ! 3D unallocated variable
    var_name = "foobar"
    call my_file%read( var_name, real3D, my_name )
    call assert( 628827846, allocated( real3D ) )
    call assert( 688571939, size( real3D, 1 ) .eq. 1 )
    call assert( 230675479, size( real3D, 2 ) .eq. 3 )
    call assert( 125526975, size( real3D, 3 ) .eq. 4 )
    call assert( 850105763, almost_equal( real3D( 1, 1, 1 ), 532.123_dk  ) )
    call assert( 231414897, almost_equal( real3D( 1, 2, 1 ), 1.5e28_dk   ) )
    call assert( 343733242, almost_equal( real3D( 1, 3, 1 ), 42.5_dk     ) )
    call assert( 723638505, almost_equal( real3D( 1, 1, 2 ), 39.25_dk    ) )
    call assert( 835956850, almost_equal( real3D( 1, 2, 2 ), 4293.12_dk  ) )
    call assert( 383324697, almost_equal( real3D( 1, 3, 2 ), 9753.231_dk ) )
    call assert( 926217023, almost_equal( real3D( 1, 1, 3 ), 3.25e-19_dk ) )
    call assert( 473584870, almost_equal( real3D( 1, 2, 3 ), 4.629e10_dk ) )
    call assert( 368436366, almost_equal( real3D( 1, 3, 3 ), 7264.12_dk  ) )
    call assert( 133271062, almost_equal( real3D( 1, 1, 4 ), 8.4918e7_dk ) )
    call assert( 757965653, almost_equal( real3D( 1, 2, 4 ), 13.2_dk     ) )
    call assert( 310597807, almost_equal( real3D( 1, 3, 4 ), 8293.12_dk  ) )
    deallocate( real3D )

    ! 3D pre-allocated variable
    var_name = "foobar"
    allocate( real3D( 1, 3, 4 ) )
    call my_file%read( var_name, real3D, my_name )
    call assert( 506458779, allocated( real3D ) )
    call assert( 618777124, size( real3D, 1 ) .eq. 1 )
    call assert( 166144971, size( real3D, 2 ) .eq. 3 )
    call assert( 895988066, size( real3D, 3 ) .eq. 4 )
    call assert( 443355913, almost_equal( real3D( 1, 1, 1 ), 532.123_dk  ) )
    call assert( 338207409, almost_equal( real3D( 1, 2, 1 ), 1.5e28_dk   ) )
    call assert( 168050505, almost_equal( real3D( 1, 3, 1 ), 42.5_dk     ) )
    call assert( 615418351, almost_equal( real3D( 1, 1, 2 ), 39.25_dk    ) )
    call assert( 727736696, almost_equal( real3D( 1, 2, 2 ), 4293.12_dk  ) )
    call assert( 892629293, almost_equal( real3D( 1, 3, 2 ), 9753.231_dk ) )
    call assert( 439997140, almost_equal( real3D( 1, 1, 3 ), 3.25e-19_dk ) )
    call assert( 334848636, almost_equal( real3D( 1, 2, 3 ), 4.629e10_dk ) )
    call assert( 164691732, almost_equal( real3D( 1, 3, 3 ), 7264.12_dk  ) )
    call assert( 612059578, almost_equal( real3D( 1, 1, 4 ), 8.4918e7_dk ) )
    call assert( 441902674, almost_equal( real3D( 1, 2, 4 ), 13.2_dk     ) )
    call assert( 606795271, almost_equal( real3D( 1, 3, 4 ), 8293.12_dk  ) )
    deallocate( real3D )

    ! 4D unallocated variable
    var_name = "corge"
    call my_file%read( var_name, real4D, my_name )
    call assert( 464572470, allocated( real4D ) )
    call assert( 911940316, size( real4D, 1 ) .eq. 2 )
    call assert( 124258662, size( real4D, 2 ) .eq. 1 )
    call assert( 571626508, size( real4D, 3 ) .eq. 3 )
    call assert( 118994355, size( real4D, 4 ) .eq. 4 )
    call assert( 913845850, almost_equal( real4D( 1, 1, 1, 1 ), 532.123_dk  ) )
    call assert( 743688946, almost_equal( real4D( 2, 1, 1, 1 ), 632.123_dk  ) )
    call assert( 291056793, almost_equal( real4D( 1, 1, 2, 1 ), 1.5e28_dk   ) )
    call assert( 738424639, almost_equal( real4D( 2, 1, 2, 1 ), 2.5e28_dk   ) )
    call assert( 285792486, almost_equal( real4D( 1, 1, 3, 1 ), 42.5_dk     ) )
    call assert( 115635582, almost_equal( real4D( 2, 1, 3, 1 ), 52.5_dk     ) )
    call assert( 227953927, almost_equal( real4D( 1, 1, 1, 2 ), 39.25_dk    ) )
    call assert( 221236381, almost_equal( real4D( 2, 1, 1, 2 ), 49.25_dk    ) )
    call assert( 398563126, almost_equal( real4D( 1, 1, 2, 2 ), 4293.12_dk  ) )
    call assert( 563455723, almost_equal( real4D( 2, 1, 2, 2 ), 5293.12_dk  ) )
    call assert( 170567663, almost_equal( real4D( 1, 1, 3, 2 ), 9753.231_dk ) )
    call assert( 335460260, almost_equal( real4D( 2, 1, 3, 2 ), 1753.231_dk ) )
    call assert( 782828106, almost_equal( real4D( 1, 1, 1, 3 ), 3.25e-19_dk ) )
    call assert( 395204353, almost_equal( real4D( 2, 1, 1, 3 ), 4.25e-19_dk ) )
    call assert( 225047449, almost_equal( real4D( 1, 1, 2, 3 ), 4.629e10_dk ) )
    call assert( 389940046, almost_equal( real4D( 2, 1, 2, 3 ), 5.629e10_dk ) )
    call assert( 554832643, almost_equal( real4D( 1, 1, 3, 3 ), 7264.12_dk  ) )
    call assert( 384675739, almost_equal( real4D( 2, 1, 3, 3 ), 8264.12_dk  ) )
    call assert( 897051985, almost_equal( real4D( 1, 1, 1, 4 ), 8.4918e7_dk ) )
    call assert( 444419832, almost_equal( real4D( 2, 1, 1, 4 ), 9.4918e7_dk ) )
    call assert( 556738177, almost_equal( real4D( 1, 1, 2, 4 ), 13.2_dk     ) )
    call assert( 104106024, almost_equal( real4D( 2, 1, 2, 4 ), 23.2_dk     ) )
    call assert( 551473870, almost_equal( real4D( 1, 1, 3, 4 ), 8293.12_dk  ) )
    call assert( 446325366, almost_equal( real4D( 2, 1, 3, 4 ), 9293.12_dk  ) )
    deallocate( real4D )

    ! 4D allocated variable
    var_name = "corge"
    allocate( real4D( 2, 1, 3, 4 ) )
    call my_file%read( var_name, real4D, my_name )
    call assert( 493635311, allocated( real4D ) )
    call assert( 106011558, size( real4D, 1 ) .eq. 2 )
    call assert( 553379404, size( real4D, 2 ) .eq. 1 )
    call assert( 100747251, size( real4D, 3 ) .eq. 3 )
    call assert( 830590346, size( real4D, 4 ) .eq. 4 )
    call assert( 660433442, almost_equal( real4D( 1, 1, 1, 1 ), 532.123_dk  ) )
    call assert( 207801289, almost_equal( real4D( 2, 1, 1, 1 ), 632.123_dk  ) )
    call assert( 102652785, almost_equal( real4D( 1, 1, 2, 1 ), 1.5e28_dk   ) )
    call assert( 550020631, almost_equal( real4D( 2, 1, 2, 1 ), 2.5e28_dk   ) )
    call assert( 997388477, almost_equal( real4D( 1, 1, 3, 1 ), 42.5_dk     ) )
    call assert( 262281075, almost_equal( real4D( 2, 1, 3, 1 ), 52.5_dk     ) )
    call assert( 157132571, almost_equal( real4D( 1, 1, 1, 2 ), 39.25_dk    ) )
    call assert( 886975666, almost_equal( real4D( 2, 1, 1, 2 ), 49.25_dk    ) )
    call assert( 151868264, almost_equal( real4D( 1, 1, 2, 2 ), 4293.12_dk  ) )
    call assert( 264186609, almost_equal( real4D( 2, 1, 2, 2 ), 5293.12_dk  ) )
    call assert( 711554455, almost_equal( real4D( 1, 1, 3, 2 ), 9753.231_dk ) )
    call assert( 258922302, almost_equal( real4D( 2, 1, 3, 2 ), 1753.231_dk ) )
    call assert( 771298548, almost_equal( real4D( 1, 1, 1, 3 ), 3.25e-19_dk ) )
    call assert( 601141644, almost_equal( real4D( 2, 1, 1, 3 ), 4.25e-19_dk ) )
    call assert( 766034241, almost_equal( real4D( 1, 1, 2, 3 ), 4.629e10_dk ) )
    call assert( 313402088, almost_equal( real4D( 2, 1, 2, 3 ), 5.629e10_dk ) )
    call assert( 208253584, almost_equal( real4D( 1, 1, 3, 3 ), 7264.12_dk  ) )
    call assert( 655621430, almost_equal( real4D( 2, 1, 3, 3 ), 8264.12_dk  ) )
    call assert( 202989277, almost_equal( real4D( 1, 1, 1, 4 ), 8.4918e7_dk ) )
    call assert( 932832372, almost_equal( real4D( 2, 1, 1, 4 ), 9.4918e7_dk ) )
    call assert( 762675468, almost_equal( real4D( 1, 1, 2, 4 ), 13.2_dk     ) )
    call assert( 874993813, almost_equal( real4D( 2, 1, 2, 4 ), 23.2_dk     ) )
    call assert( 139886411, almost_equal( real4D( 1, 1, 3, 4 ), 8293.12_dk  ) )
    call assert( 317213156, almost_equal( real4D( 2, 1, 3, 4 ), 9293.12_dk  ) )
    deallocate( real4D )

    ! dimension names
    var_name = "qux"
    dim_names = my_file%variable_dimensions( var_name, my_name )
    call assert( 685336671, allocated( dim_names ) )
    call assert( 410031263, size( dim_names ) .eq. 0 )
    deallocate( dim_names )

    var_name = "corge"
    dim_names = my_file%variable_dimensions( var_name, my_name )
    call assert( 513726528, allocated( dim_names ) )
    call assert( 562942007, size( dim_names ) .eq. 4 )
    call assert( 282372292, dim_names( 1 ) .eq. "i" )
    call assert( 619327327, dim_names( 2 ) .eq. "h" )
    call assert( 166695174, dim_names( 3 ) .eq. "g" )
    call assert( 896538269, dim_names( 4 ) .eq. "f" )
    deallocate( dim_names )

    ! clean up
    deallocate( my_file )

  end subroutine test_io_netcdf_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_io_netcdf
