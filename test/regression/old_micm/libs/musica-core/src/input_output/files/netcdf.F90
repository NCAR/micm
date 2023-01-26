! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_file_netcdf module

!> The file_netcdf_t type and related functions
module musica_file_netcdf

  use musica_constants,                only : musica_ik, musica_dk
  use musica_file,                     only : file_t

  implicit none
  private

  public :: file_netcdf_t

  !> A NetCDF file
  type, extends(file_t) :: file_netcdf_t
    private
    !> Flag indicating whether the file is open
    logical :: is_open_ = .false.
    !> Flag indicating whether a file is in data mode
    logical :: is_in_data_mode_ = .false.
    !> NetCDF file id
    integer(kind=musica_ik) :: id_
    !> NetCDF dimension id for time (output files only)
    integer(kind=musica_ik) :: time_dimid_
    !> NetCDF variable id for time (output files only)
    integer(kind=musica_ik) :: time_varid_
    !> Last output time value [s]
    real(kind=musica_dk) :: last_output_time__s_ = -9999e300_musica_dk
    !> Last output time index
    integer(kind=musica_ik) :: last_output_time_index_ = 0
  contains
    !> Returns the type of file as a string
    procedure :: type => file_type
    !> Returns the NetCDF file id
    procedure :: id
    !> Returns the number of dimensions in the file
    procedure :: number_of_dimensions
    !> Returns the number of variables in the file
    procedure :: number_of_variables
    !> Returns a flag indicating whether a variable exists in the file
    procedure :: is_file_variable
    !> Opens the file if it is not currently open
    procedure :: check_open
    !> Enters data mode if not currently in data mode
    procedure :: check_data_mode
    !> Sets the current simulation time (output files only)
    procedure :: set_output_time__s
    !> Returns the current index in the temporal dimension
    procedure :: current_time_index
    !> Checks a returned NetCDF status code and fail if an error occurred
    procedure :: check_status
    !> Closes the file
    procedure :: close
    !> Creates a file for output
    procedure, private :: create_output_file
    !> Finalizes the file
    final :: finalize, finalize_1D_array
  end type file_netcdf_t

  !> Constructor
  interface file_netcdf_t
    module procedure :: constructor
  end interface file_netcdf_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates a file_netcdf_t object for a NetCDF file
  function constructor( config ) result( new_obj )

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_file,                   only : private_constructor

    !> New NetCDF file object
    class(file_t), pointer :: new_obj
    !> File configuration
    type(config_t), intent(inout) :: config

    allocate( file_netcdf_t :: new_obj )

    call private_constructor( new_obj, config )
    call assert_msg( 708084310,                                               &
                ( new_obj%is_input( ) .and. .not. new_obj%is_output( ) ) .or. &
                ( .not. new_obj%is_input( ) .and. new_obj%is_output( ) ),     &
                "NetCDF files must be either input or output." )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the type of file as a string
  type(string_t) function file_type( this )

    use musica_string,                 only : string_t

    !> NetCDF file
    class(file_netcdf_t), intent(in) :: this

    file_type = "netcdf"

  end function file_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the NetCDF file id
  function id( this )

    !> NetCDF file id
    integer(kind=musica_ik) :: id
    !> NetCDF file
    class(file_netcdf_t), intent(inout) :: this

    call this%check_open( )
    id = this%id_

  end function id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of dimensions in the file
  integer(kind=musica_ik) function number_of_dimensions( this )

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : musica_ik
    use musica_string,                 only : to_char, string_t
    use netcdf,                        only : nf90_inquire

    !> NetCDF file
    class(file_netcdf_t), intent(inout) :: this

    type(string_t) :: file_name

    call this%check_data_mode( )
    file_name = this%name( )
    call this%check_status( 639211977,                                        &
                        nf90_inquire( this%id_,                               &
                                      nDimensions = number_of_dimensions ),   &
                       "Error getting dimension information for NetCDF file '"&
                       //file_name%to_char( )//"'" )
    call assert_msg( 938882534, number_of_dimensions .eq. 1,                  &
                     "NetCDF files are currently only set up for one "//      &
                     "dimension of time. File '"//file_name%to_char( )//      &
                     "' has "//to_char( number_of_dimensions )//" dimensions" )

  end function number_of_dimensions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of variables in the file
  integer(kind=musica_ik) function number_of_variables( this )

    use musica_constants,              only : musica_ik
    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_inquire

    !> NetCDF file
    class(file_netcdf_t), intent(inout) :: this

    type(string_t) :: file_name

    call this%check_data_mode( )
    file_name = this%name( )
    call this%check_status( 150118878,                                        &
                       nf90_inquire( this%id_,                                &
                                     nVariables = number_of_variables ),      &
                       "Error getting variable information for NetCDF file '" &
                       //file_name%to_char( )//"'" )

  end function number_of_variables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns a flag indicating whether a variable exists in the file
  logical function is_file_variable( this, variable_name )

    use netcdf,                        only : nf90_inq_varid, NF90_NOERR

    !> NetCDF file
    class(file_netcdf_t), intent(inout) :: this
    !> Variable to look for
    character(len=*), intent(in) :: variable_name

    integer(kind=musica_ik) :: status, var_id

    call this%check_open( )
    is_file_variable = .false.
    status = nf90_inq_varid( this%id_, variable_name, var_id )
    if( status .eq. NF90_NOERR ) then
      is_file_variable = .true.
    end if

  end function is_file_variable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Opens the file if it is not open already
  subroutine check_open( this )

    use musica_assert,                 only : die_msg
    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_open, NF90_NOWRITE

    !> NetCDF file
    class(file_netcdf_t), intent(inout) :: this

    type(string_t) :: file_name

    if( this%is_open_ ) return
    file_name = this%name( )
    if( this%is_input( ) .and. .not. this%is_output( ) ) then
      call this%check_status( 172405314,                                      &
          nf90_open( file_name%to_char( ), NF90_NOWRITE, this%id_ ),          &
          "Error opening NetCDF file '"//file_name%to_char( )//"'" )
      this%is_in_data_mode_ = .true.
    else if( this%is_output( ) .and. .not. this%is_input( ) ) then
      call this%create_output_file( )
    else
      call die_msg( 154946927, "Mixed input/output NetCDF files are not "//   &
                    "supported for file '"//file_name%to_char( )//"'" )
    end if
    this%is_open_ = .true.

  end subroutine check_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Enters data mode if the file is not already in data mode
  subroutine check_data_mode( this )

    use musica_assert,                 only : assert
    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_enddef, NF90_CLOBBER

    !> NetCDF file
    class(file_netcdf_t), intent(inout) :: this

    type(string_t) :: file_name

    call this%check_open( )
    if( this%is_in_data_mode_ ) return
    call assert( 562436195, this%is_output( ) .and. .not. this%is_input( ) )
    file_name = this%name( )
    call this%check_status( 567692478, nf90_enddef( this%id_ ),               &
          "Error ending define mode for NetCDF output file '"//               &
          file_name%to_char( )//"'" )
    this%is_in_data_mode_ = .true.

  end subroutine check_data_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets the current simulation time [s] (output files only)
  subroutine set_output_time__s( this, time__s )

    use musica_assert,                 only : assert
    use netcdf,                        only : nf90_put_var

    !> NetCDF file
    class(file_netcdf_t), intent(inout) :: this
    !> Current simulation time [s]
    real(kind=musica_dk), intent(in) :: time__s

    call this%check_data_mode( )
    if( this%last_output_time__s_ .eq. time__s ) return
    call assert( 156844437, time__s .gt. this%last_output_time__s_ )
    this%last_output_time__s_    = time__s
    this%last_output_time_index_ = this%last_output_time_index_ + 1
    call this%check_status( 718436162,                                        &
        nf90_put_var( this%id_, this%time_varid_, time__s,                    &
                      start = (/ this%last_output_time_index_ /) ),           &
        "Error adding time to NetCDF file" )

  end subroutine set_output_time__s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the current index in the temporal dimension
  integer(kind=musica_ik) function current_time_index( this )

    use musica_assert,                 only : assert

    !> NetCDF file
    class(file_netcdf_t), intent(inout) :: this

    call this%check_data_mode( )
    call assert( 329646293, this%last_output_time_index_ .gt. 0 )
    current_time_index = this%last_output_time_index_

  end function current_time_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Checks a NetCDF status code and fail with a message if an error occurred
  subroutine check_status( this, code, status, error_message )

    use musica_assert,                 only : die_msg
    use musica_string,                 only : string_t
    use netcdf,                        only : NF90_NOERR, nf90_strerror

    !> NetCDF file
    class(file_netcdf_t), intent(in) :: this
    !> Unique code for the assertion
    integer(kind=musica_ik), intent(in) :: code
    !> Status code
    integer(kind=musica_ik), intent(in) :: status
    !> Error message
    character(len=*), intent(in) :: error_message

    type(string_t) :: file_name

    if( status .eq. NF90_NOERR ) return
    file_name = this%name( )
    call die_msg( code, "NetCDF file '"//file_name%to_char( )//               &
                  "': "//trim( error_message )//": "//                        &
                  trim( nf90_strerror( status ) ) )

  end subroutine check_status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Closes the NetCDF file
  subroutine close( this )

    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_close

    !> NetCDF file
    class(file_netcdf_t), intent(inout) :: this

    type(string_t) :: file_name

    if( .not. this%is_open_ ) return
    file_name = this%name( )
    call this%check_status( 633660547, nf90_close( this%id_ ),                &
                            "Error closing NetCDF file '"//                   &
                            file_name%to_char( )//"'" )
    this%is_open_                = .false.
    this%is_in_data_mode_        = .false.
    this%last_output_time__s_    = -9999e300_musica_dk
    this%last_output_time_index_ = 0

  end subroutine close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates a file for output
  !!
  !! \todo To be replaced when column data is included in I/O
  subroutine create_output_file( this )

    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_create, nf90_def_dim,      &
                                              nf90_def_var, nf90_put_att,     &
                                              NF90_CLOBBER, NF90_DOUBLE,      &
                                              NF90_UNLIMITED

    !> NetCDF file
    class(file_netcdf_t), intent(inout) :: this

    type(string_t) :: file_name

    file_name = this%name( )
    call this%check_status( 407715447,                                        &
        nf90_create( file_name%to_char( ), NF90_CLOBBER, this%id_ ),          &
        "Error opening NetCDF file '"//file_name%to_char( )//"'" )
    call this%check_status( 388625371,                                        &
        nf90_def_dim( this%id_, "time", NF90_UNLIMITED, this%time_dimid_ ),   &
        "Error creating time dimension in NetCDF file '"//                    &
        file_name%to_char( )//"'" )
    call this%check_status( 764719566,                                        &
        nf90_def_var( this%id_, "time", NF90_DOUBLE, this%time_dimid_,        &
                      this%time_varid_ ),                                     &
        "Error creating time variable in NetCDF file '"//                     &
        file_name%to_char( )//"'" )
    call this%check_status( 163555675,                                        &
        nf90_put_att( this%id_, this%time_varid_, "units", "s" ),             &
        "Error adding units attribute for time in NetCDF file '"//            &
        file_name%to_char( )//"'" )

  end subroutine create_output_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalizes a file object
  subroutine finalize( this )

    !> NetCDF file
    type(file_netcdf_t), intent(inout) :: this

    call this%close( )

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalizes an array of file objects
  subroutine finalize_1D_array( this )

    !> NetCDF file
    type(file_netcdf_t), intent(inout) :: this(:)

    integer :: i

    do i = 1, size( this )
      call this( i )%close( )
    end do

  end subroutine finalize_1D_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_file_netcdf
