! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_file_variable_netcdf module

!> The file_variable_netcdf_t type and related functions
module musica_file_variable_netcdf

  use musica_constants,                only : musica_dk, musica_ik
  use musica_file_dimension_range,     only : file_dimension_range_t
  use musica_file_variable,            only : file_variable_t

  implicit none
  private

  public :: file_variable_netcdf_t

  !> A NetCDF variable
  !!
  !! In addition to standard file_variable_t functions, a NetCDF variable
  !! provides access to the NetCDF variable id.
  !!
  type, extends(file_variable_t) :: file_variable_netcdf_t
    private
    !> NetCDF variable id
    integer(kind=musica_ik) :: id_ = -1
    !> Variable dimensions
    type(file_dimension_range_t), allocatable :: dimensions_(:)
  contains
    !> Returns the NetCDF variable id
    procedure :: id
    !> Gets the variable dimensions
    procedure :: get_dimensions
    !> Gets a sub-set of the variable data for a specified index range
    !!
    !! Data are returned after applying conversions set up during
    !! initialization.
    !!
    procedure :: get_data
    !> Outputs data to the file for a given time step
    procedure :: output
    !> Returns a flag indicating whether two file_variable_t objects refer to
    !! the same file variable
    procedure :: is_same_as
    !> Prints the variable properties
    procedure :: print => do_print
  end type file_variable_netcdf_t

  !> Constructor
  interface file_variable_netcdf_t
    module procedure :: constructor_name, constructor_id
  end interface file_variable_netcdf_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates a file_variable_netcdf_t object for an existing NetCDF variable
  !> by name
  !!
  function constructor_name( file, variable_name, dimensions, config, found ) &
      result( new_obj )

    use musica_assert,                 only : assert_msg, die, die_msg
    use musica_config,                 only : config_t
    use musica_file,                   only : file_t
    use musica_file_netcdf,            only : file_netcdf_t
    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_inq_varid, nf90_def_var,   &
                                              nf90_put_att, NF90_DOUBLE

    !> New NetCDF variable
    class(file_variable_t), pointer :: new_obj
    !> NetCDF file
    class(file_t), intent(inout) :: file
    !> Variable name
    character(len=*), intent(in) :: variable_name
    !> Variable dimensions
    !!
    !! Only required for adding a file variable
    !!
    type(file_dimension_range_t), intent(in), optional :: dimensions(:)
    !> Configuration describing how to match to MUSICA variables
    !!
    !! If omitted, standard matching is applied
    type(config_t), intent(inout), optional :: config
    !> Optional flag that indicates whether the variable was found in the file
    logical, intent(out), optional :: found

    character(len=*), parameter :: my_name = "NetCDF variable constructor"
    integer(kind=musica_ik) :: variable_id
    integer(kind=musica_ik), allocatable :: dims(:)
    type(string_t) :: units
    logical :: l_found

    if( present( found ) ) then
      if( file%is_file_variable( variable_name ) ) then
        found = .true.
      else
        found = .false.
        return
      end if
    end if
    select type( file )
    class is( file_netcdf_t )
      call file%check_open( )
      if( file%is_input( ) .and. .not. file%is_output( ) ) then
        call file%check_status( 542234258,                                    &
            nf90_inq_varid( file%id( ), variable_name, variable_id ),         &
            "Error getting variable id for '"//variable_name//"'" )
      else
        if( file%is_file_variable( variable_name ) ) then
          call file%check_status( 802998523,                                  &
              nf90_inq_varid( file%id( ), variable_name, variable_id ),       &
              "Internal error" )
        else
          call assert_msg( 587613594, present( dimensions ),                  &
                           "Missing dimensions for new variable '"//          &
                           variable_name//"'" )
          call assert_msg( 586160355, present( config ),                      &
                           "Missing configuration for new variable '"//       &
                           variable_name//"'" )
          call config%get( "units", units, my_name, found = l_found )
          call assert_msg( 637281368, l_found, "No units specified for new "//&
                           "variable '"//variable_name//"'" )
          call assert_msg( 178020050, size( dimensions ) .ge. 1,              &
                           "NetCDF variable '"//variable_name//               &
                           "' must have at least 1 dimension." )
          allocate( dims( size( dimensions ) ) )
          dims(:) = dimensions(:)%id( )
          call file%check_status( 571644287,                                  &
              nf90_def_var( file%id( ), variable_name, NF90_DOUBLE,           &
                            dims, variable_id ),                              &
              "Error creating NetCDF variable '"//variable_name//"'" )
          call file%check_status( 491375293,                                  &
              nf90_put_att( file%id( ), variable_id, "units",                 &
                            units%to_char( ) ),                               &
              "Error setting units for new NetCDF variable '"//variable_name//&
              "'" )
        end if
      end if
    class default
      call die( 952578143 )
    end select
    new_obj => constructor_id( file, variable_id, config = config )

  end function constructor_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates a file_variable_netcdf_t object for an existing NetCDF variable
  !> by id
  !!
  function constructor_id( file, variable_id, config ) result( new_obj )

    use musica_assert,                 only : assert_msg, die
    use musica_config,                 only : config_t
    use musica_file,                   only : file_t
    use musica_file_netcdf,            only : file_netcdf_t
    use musica_string,                 only : string_t, to_char
    use netcdf,                        only : NF90_MAX_NAME,                  &
                                              nf90_inquire,                   &
                                              nf90_inquire_variable,          &
                                              nf90_inquire_dimension,         &
                                              nf90_inq_attname,               &
                                              nf90_get_att

    !> New NetCDF variable
    class(file_variable_t), pointer :: new_obj
    !> NetCDF file
    class(file_t), intent(inout) :: file
    !> Variable ID
    integer(kind=musica_ik), intent(in) :: variable_id
    !> Configuration describing how to match to MUSICA variables
    !!
    !! If omitted, standard matching is applied
    type(config_t), intent(inout), optional :: config

    character(len=NF90_MAX_NAME) :: name, units, dim_name
    type(string_t) :: file_name, att_name, var_name, var_units
    type(string_t), allocatable :: dimension_names(:)
    integer(kind=musica_ik) :: n_values, i_att, n_attributes, unlimited_dimid,&
                               i_dim, n_dims
    integer(kind=musica_ik), allocatable :: dimids(:)

    allocate( file_variable_netcdf_t :: new_obj )

    select type( new_obj )
    class is( file_variable_netcdf_t )
      select type( file )
      class is( file_netcdf_t )
        file_name = file%name( )
        call file%check_open( )
        call file%check_status( 404330062,                                    &
                                nf90_inquire_variable( file%id( ),            &
                                                       variable_id,           &
                                                       ndims = n_dims ),      &
                                "Error getting number of dimensions for "//   &
                                "id: "//to_char( variable_id ) )
        allocate( dimids( n_dims ) )
        call file%check_status( 206732462,                                    &
                                nf90_inquire_variable( file%id( ),            &
                                                       variable_id,           &
                                                       name = name,           &
                                                       dimids = dimids,       &
                                                       nAtts = n_attributes ),&
                            "Error getting variable information for id: "//   &
                            to_char( variable_id ) )
        var_name    = name
        var_units   = ""
        new_obj%id_ = variable_id
        allocate( new_obj%dimensions_( size( dimids ) ) )
        call file%check_status( 683946769,                                    &
                                nf90_inquire( file%id( ),                     &
                                          unlimiteddimid = unlimited_dimid ), &
                                "Error getting unlimited dimension id" )
        allocate( dimension_names( size( dimids ) ) )
        do i_dim = 1, size( dimids )
          call file%check_status( 661270149,                                  &
                              nf90_inquire_dimension( file%id( ),             &
                                                      dimids( i_dim ),        &
                                                      name = dim_name,        &
                                                      len = n_values ),       &
                              "Error getting dimensions of variable '"//      &
                              var_name%to_char( )//"'" )
          dimension_names( i_dim ) = dim_name
          new_obj%dimensions_( i_dim ) =                                      &
              file_dimension_range_t( dimids( i_dim ), 1, n_values,           &
                      is_unlimited = ( dimids( i_dim ) .eq. unlimited_dimid ) )
        end do
        units = ""
        do i_att = 1, n_attributes
          call file%check_status( 485848938,                                  &
                              nf90_inq_attname( file%id( ),                   &
                                                variable_id,                  &
                                                i_att,                        &
                                                name ),                       &
                              "Error getting attribute "//to_char( i_att )//  &
                              " name for variable '"//                        &
                              var_name%to_char( )//"'" )
          att_name = trim( name )
          att_name = att_name%to_lower( )
          if( att_name .eq. "units" .or. att_name .eq. "unit" ) then
            call file%check_status( 992960877,                                &
                                nf90_get_att( file%id( ), new_obj%id_, name,  &
                                              units ),                        &
                                "Error getting units for variable '"//        &
                                var_name%to_char( )//"'" )
            var_units = trim( units )
          end if
        end do
        call assert_msg( 738503497, var_units .ne. "",                        &
                         "No units found for variable '"//var_name%to_char( ) &
                         //"' in NetCDF file '"//file_name%to_char( )//"'" )
      class default
        call die( 185049127 )
      end select
    class default
      call die( 857053663 )
    end select

    call new_obj%private_constructor( config, file, var_name, dimension_names,&
                                      units = var_units )

  end function constructor_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the NetCDF variable index
  integer(kind=musica_ik) function id( this )

    !> NetCDF variable
    class(file_variable_netcdf_t), intent(in) :: this

    id = this%id_

  end function id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gets the variable dimensions
  function get_dimensions( this ) result( dimensions )

    !> Variable dimensions
    type(file_dimension_range_t), allocatable :: dimensions(:)
    !> NetCDF variable
    class(file_variable_netcdf_t), intent(in) :: this

    dimensions = this%dimensions_

  end function get_dimensions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gets a sub-set of the data from the file
  !!
  !! Applies any necessary conversions to the raw input data
  !!
  subroutine get_data( this, file, range, values )

    use musica_assert,                 only : assert, die, die_msg
    use musica_file,                   only : file_t
    use musica_file_netcdf,            only : file_netcdf_t
    use netcdf,                        only : nf90_get_var

    !> NetCDF variable
    class(file_variable_netcdf_t), intent(in) :: this
    !> NetCDF file
    class(file_t), intent(inout) :: file
    !> Range of data to return
    !!
    !! Only one range in the array is permitted to have size > 1, so that
    !! results are always returned as a rank 1 array
    !!
    type(file_dimension_range_t), intent(in) :: range(:)
    !> Values to return
    real(kind=musica_dk), target, intent(out) :: values(:)

    real(kind=musica_dk), pointer :: v2(:,:), v3(:,:,:), v4(:,:,:,:)
    integer(kind=musica_ik), allocatable :: n_values(:), start(:)
    integer(kind=musica_ik) :: i_val, i_dim, j_dim, temp, temp2
    logical :: found

    select type( file )
    class is( file_netcdf_t )
      allocate( start(    size( this%dimensions_ ) ) )
      allocate( n_values( size( this%dimensions_ ) ) )
      call assert( 448155683, size( this%dimensions_ ) .eq. size( range ) )
      do i_dim = 1, size( range )
        found = .false.
        do j_dim = 1, size( this%dimensions_ )
          if( range( i_dim ) .eq. this%dimensions_( j_dim ) ) then
            found = .true.
            start(    j_dim ) = range( i_dim )%lower_bound( )
            n_values( j_dim ) = range( i_dim )%upper_bound( ) -               &
                                start( j_dim ) + 1
            exit
          end if
        end do
        call assert( 205271889, found )
      end do
      temp = 1
      temp2 = 0
      do i_dim = 1, size( n_values )
        temp = temp * n_values( i_dim )
        temp2 = temp2 + n_values( i_dim )
      end do
      call assert( 394077671, temp .eq. temp2 - size( n_values ) + 1 )
      call assert( 509603293, temp .le. size( values ) )
      if( size( n_values ) .eq. 1 ) then
        call file%check_status( 448163017, nf90_get_var( file%id( ), this%id_,&
                                                         values, start,       &
                                                         n_values ),          &
                                "Error getting values for NetCDF variable" )
      else if( size( n_values ) .eq. 2 ) then
        v2( 1:n_values(1), 1:n_values(2) ) => values
        call file%check_status( 689836516, nf90_get_var( file%id( ), this%id_,&
                                                         v2, start,           &
                                                         n_values ),          &
                                "Error getting values for NetCDF variable" )
      else if( size( n_values ) .eq. 3 ) then
        v3( 1:n_values(1), 1:n_values(2), 1:n_values(3) ) => values
        call file%check_status( 184630111, nf90_get_var( file%id( ), this%id_,&
                                                          v3, start,          &
                                                          n_values ),         &
                               "Error getting values for NetCDF variable" )
      else if( size( n_values ) .eq. 4 ) then
        v4( 1:n_values(1), 1:n_values(2), 1:n_values(3), 1:n_values(4) ) =>   &
            values
        call file%check_status( 574159398, nf90_get_var( file%id( ), this%id_,&
                                                         v4, start,           &
                                                         n_values ),          &
                                "Error getting values for NetCDF variable" )
      else
        call die_msg( 468558599, "NetCDF variables above 4 dimensions are "// &
                                 "not yet supported." )
      end if
      call this%convert_to_musica_values( values )
    class default
      call die( 636979049 )
    end select

  end subroutine get_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Outputs data to the file for a given timestep
  !!
  !! The state_value will be converted according to the variable configuration
  !! prior to outputting it to the file.
  !!
  subroutine output( this, file, unlimited_dimension_value, indices,          &
      state_value )

    use musica_assert,                 only : assert, die, die_msg
    use musica_file,                   only : file_t
    use musica_file_netcdf,            only : file_netcdf_t
    use netcdf,                        only : nf90_put_var

    !> NetCDF variable
    class(file_variable_netcdf_t), intent(inout) :: this
    !> Output file
    class(file_t), intent(inout) :: file
    !> Unlimited dimension value
    real(kind=musica_dk), intent(in) :: unlimited_dimension_value
    !> Indices of data point to output in all non-unlimited dimensions
    class(file_dimension_range_t), intent(in) :: indices(:)
    !> Value to output
    real(kind=musica_dk), intent(in) :: state_value

    real(kind=musica_dk), pointer :: v2(:,:), v3(:,:,:), v4(:,:,:,:)
    integer(kind=musica_ik), allocatable :: n_values(:), start(:)
    integer(kind=musica_ik) :: i_val, i_dim, j_dim, temp
    real(kind=musica_dk), target :: file_val(1)
    logical :: found

    select type( file )
    class is( file_netcdf_t )
      call file%set_output_time__s( unlimited_dimension_value )
      allocate( start(    size( this%dimensions_ ) ) )
      allocate( n_values( size( this%dimensions_ ) ) )
      call assert( 164507326,                                                 &
                   size( this%dimensions_ ) .eq. size( indices ) + 1 )
      n_values(:) = 1
      start(1)    = file%current_time_index( )
      do i_dim = 1, size( indices )
        found = .false.
        do j_dim = 1, size( this%dimensions_ )
          if( indices( i_dim ) .eq. this%dimensions_( j_dim ) ) then
            start(    j_dim ) = indices( i_dim )%lower_bound( )
            call assert( 157075959, indices( i_dim )%lower_bound( ) .eq.      &
                                    indices( i_dim )%upper_bound( ) )
            exit
          end if
        end do
        call assert( 749183764, found )
      end do
      file_val(1) = state_value
      call this%convert_to_file_values( file_val )
      if( size( n_values ) .eq. 1 ) then
        call file%check_status( 970461681, nf90_put_var( file%id( ), this%id_,&
                                                         file_val, start,     &
                                                         n_values ),          &
                                "Error setting values for NetCDF variable" )
      else if( size( n_values ) .eq. 2 ) then
        v2( 1:n_values(1), 1:n_values(2) ) => file_val
        call file%check_status( 875676619, nf90_put_var( file%id( ), this%id_,&
                                                         v2, start,           &
                                                         n_values ),          &
                                "Error setting values for NetCDF variable" )
      else if( size( n_values ) .eq. 3 ) then
        v3( 1:n_values(1), 1:n_values(2), 1:n_values(3) ) => file_val
        call file%check_status( 928703166, nf90_put_var( file%id( ), this%id_,&
                                                         v3, start,           &
                                                         n_values ),          &
                                "Error setting values for NetCDF variable" )
      else if( size( n_values ) .eq. 4 ) then
        v4( 1:n_values(1), 1:n_values(2), 1:n_values(3), 1:n_values(4) ) =>   &
            file_val
        call file%check_status( 920080086, nf90_put_var( file%id( ), this%id_,&
                                                         v4, start,           &
                                                         n_values ),          &
                                "Error setting values for NetCDF variable" )
      else
        call die_msg( 351770815, "NetCDF variables above 4 dimensions are "// &
                                 "not yet supported." )
      end if
    class default
      call die( 368678658 )
    end select

  end subroutine output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns a flag indicating whether two file_variable_t objects refer to
  !! the same file variable
  logical elemental function is_same_as( a, b )

    !> File variable a
    class(file_variable_netcdf_t), intent(in) :: a
    !> File variable b
    class(file_variable_t), intent(in) :: b

    is_same_as = .false.
    select type( b )
    class is( file_variable_netcdf_t )
      is_same_as = a%id_ .eq. b%id_
    end select

  end function is_same_as

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Prints the properties of the variable
  subroutine do_print( this )

    use musica_string,                 only : string_t

    !> NetCDF variable
    class(file_variable_netcdf_t), intent(in) :: this

    type(string_t) :: units, musica_name, var_name

    var_name    = this%name( )
    musica_name = this%musica_name( )
    units       = this%units( )
    write(*,*) "*** Variable: "//var_name%to_char( )//" ***"
    write(*,*) "MUSICA name: "//musica_name%to_char( )
    write(*,*) "NetCDF variable id:", this%id_
    write(*,*) "units: "//units%to_char( )

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_file_variable_netcdf
