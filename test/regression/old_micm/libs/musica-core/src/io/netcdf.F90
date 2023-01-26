! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_io_netcdf module

!> The io_netcdf_t type and related functions
module musica_io_netcdf

  use musica_io,                       only : io_t
  use musica_string,                   only : string_t

  implicit none
  private

  public :: io_netcdf_t

  integer, parameter :: kUnknownFileId = -9999

  !> NetCDF file reader
  type, extends(io_t) :: io_netcdf_t
    integer        :: file_id_ = kUnknownFileId
    type(string_t) :: file_name_
  contains
    !> @name Data read functions
    !! @{
    procedure :: read_0D_double
    procedure :: read_1D_double
    procedure :: read_2D_double
    procedure :: read_3D_double
    procedure :: read_4D_double
    procedure :: read_0D_int
    procedure :: read_1D_int
    !> @}
    !> Returns the dimension names for a given variable
    procedure :: variable_dimensions
    !> Returns the units for a given variable
    procedure :: variable_units
    procedure, private :: is_open
    procedure, private :: variable_id
    procedure, private :: dimension_sizes
    final :: finalize
  end type io_netcdf_t

  interface io_netcdf_t
    procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for NetCDF file readers
  function constructor( file_name ) result( new_io )

    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_open, NF90_NOWRITE

    type(io_netcdf_t), pointer    :: new_io
    type(string_t),    intent(in) :: file_name

    allocate( new_io )
    new_io%file_name_ = file_name
    call check_status( 126279520,                                             &
        nf90_open( file_name%to_char( ), NF90_NOWRITE, new_io%file_id_ ),     &
        "Error openning file '"//file_name%to_char( )//"'" )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 0D double-precision floating-pointer data
  subroutine read_0D_double( this, variable_name, container, requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : musica_dk
    use musica_string,                 only : to_char
    use netcdf,                        only : nf90_get_var

    class(io_netcdf_t),   intent(inout) :: this
    class(string_t),      intent(in)    :: variable_name
    real(kind=musica_dk), intent(out)   :: container
    character(len=*),     intent(in)    :: requestor_name

    integer :: var_id
    integer, allocatable :: dim_sizes(:)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 879207328, this%is_open( ),                              &
                     "Trying to read from unopen file: "//id_str )
    var_id    = this%variable_id( variable_name )
    dim_sizes = this%dimension_sizes( variable_name )
    call assert_msg( 712409197, size( dim_sizes ) .eq. 0,                     &
                     "Wrong number of dimensions for "//                      &
                     trim( id_str%to_char( ) )//": Expected 0 got "//         &
                     trim( to_char( size( dim_sizes ) ) ) )
    call check_status( 190408927,                                             &
                       nf90_get_var( this%file_id_, var_id, container ),      &
                       "Error getting value for "//trim( id_str%to_char( ) ) )

  end subroutine read_0D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 1D double-precision floating-pointer data
  subroutine read_1D_double( this, variable_name, container, requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : musica_dk
    use musica_string,                 only : to_char
    use netcdf,                        only : nf90_get_var

    class(io_netcdf_t),                intent(inout) :: this
    class(string_t),                   intent(in)    :: variable_name
    real(kind=musica_dk), allocatable, intent(out)   :: container(:)
    character(len=*),                  intent(in)    :: requestor_name

    integer :: var_id
    integer, allocatable :: dim_sizes(:)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 163123652, this%is_open( ),                              &
                     "Trying to read from unopen file: "//id_str )
    var_id    = this%variable_id( variable_name )
    dim_sizes = this%dimension_sizes( variable_name )
    call assert_msg( 275441997, size( dim_sizes ) .eq. 1,                     &
                     "Wrong number of dimensions for "//                      &
                     trim( id_str%to_char( ) )//": Expected 1 got "//         &
                     trim( to_char( size( dim_sizes ) ) ) )
    if( allocated( container ) ) then
      call assert_msg( 976961669, size( container ) .eq. dim_sizes(1),        &
                       "Wrong size container for "//trim( id_str%to_char( ) ) &
                       //": Expected "//trim( to_char( dim_sizes(1) ) )//     &
                       " got "//trim( to_char( size( container ) ) ) )
    else
      allocate( container( dim_sizes(1) ) )
    end if
    call check_status( 722809843,                                             &
                       nf90_get_var( this%file_id_, var_id, container ),      &
                       "Error getting value for "//trim( id_str%to_char( ) ) )

  end subroutine read_1D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 2D double-precision floating-pointer data
  subroutine read_2D_double( this, variable_name, container, requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : musica_dk
    use musica_string,                 only : to_char
    use netcdf,                        only : nf90_get_var

    class(io_netcdf_t),                intent(inout) :: this
    class(string_t),                   intent(in)    :: variable_name
    real(kind=musica_dk), allocatable, intent(out)   :: container(:,:)
    character(len=*),                  intent(in)    :: requestor_name

    integer :: var_id, i_dim
    integer, allocatable :: dim_sizes(:)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 675787021, this%is_open( ),                              &
                     "Trying to read from unopen file: "//id_str )
    var_id    = this%variable_id( variable_name )
    dim_sizes = this%dimension_sizes( variable_name )
    call assert_msg( 400481613, size( dim_sizes ) .eq. 2,                     &
                     "Wrong number of dimensions for "//                      &
                     trim( id_str%to_char( ) )//": Expected 2 got "//         &
                     trim( to_char( size( dim_sizes ) ) ) )
    if( allocated( container ) ) then
      do i_dim = 1, 2
        call assert_msg( 230324709, size( container, i_dim ) .eq.             &
                                    dim_sizes( i_dim ),                       &
                         "Wrong size container for "//                        &
                         trim( id_str%to_char( ) )//": Expected "//           &
                         trim( to_char( dim_sizes( i_dim ) ) )//     &
                         " got "//trim( to_char( size( container, i_dim ) ) ) &
                         //" for dimension "//trim( to_char( i_dim ) ) )
      end do
    else
      allocate( container( dim_sizes(1), dim_sizes(2) ) )
    end if
    call check_status( 960167804,                                             &
                       nf90_get_var( this%file_id_, var_id, container ),      &
                       "Error getting value for "//trim( id_str%to_char( ) ) )

  end subroutine read_2D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 3D double-precision floating-pointer data
  subroutine read_3D_double( this, variable_name, container, requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : musica_dk
    use musica_string,                 only : to_char
    use netcdf,                        only : nf90_get_var

    class(io_netcdf_t),                intent(inout) :: this
    class(string_t),                   intent(in)    :: variable_name
    real(kind=musica_dk), allocatable, intent(out)   :: container(:,:,:)
    character(len=*),                  intent(in)    :: requestor_name

    integer :: var_id, i_dim
    integer, allocatable :: dim_sizes(:)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 539957265, this%is_open( ),                              &
                     "Trying to read from unopen file: "//id_str )
    var_id    = this%variable_id( variable_name )
    dim_sizes = this%dimension_sizes( variable_name )
    call assert_msg( 603060131, size( dim_sizes ) .eq. 3,                     &
                     "Wrong number of dimensions for "//                      &
                     trim( id_str%to_char( ) )//": Expected 3 got "//         &
                     trim( to_char( size( dim_sizes ) ) ) )
    if( allocated( container ) ) then
      do i_dim = 1, 3
        call assert_msg( 715378476, size( container, i_dim ) .eq.             &
                                    dim_sizes( i_dim ),                       &
                         "Wrong size container for "//                        &
                         trim( id_str%to_char( ) )//": Expected "//           &
                         trim( to_char( dim_sizes( i_dim ) ) )//     &
                         " got "//trim( to_char( size( container, i_dim ) ) ) &
                         //" for dimension "//trim( to_char( i_dim ) ) )
      end do
    else
      allocate( container( dim_sizes(1), dim_sizes(2), dim_sizes(3) ) )
    end if
    call check_status( 210172071,                                             &
                       nf90_get_var( this%file_id_, var_id, container ),      &
                       "Error getting value for "//trim( id_str%to_char( ) ) )

  end subroutine read_3D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 4D double-precision floating-pointer data
  subroutine read_4D_double( this, variable_name, container, requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : musica_dk
    use musica_string,                 only : to_char
    use netcdf,                        only : nf90_get_var

    class(io_netcdf_t),                intent(inout) :: this
    class(string_t),                   intent(in)    :: variable_name
    real(kind=musica_dk), allocatable, intent(out)   :: container(:,:,:,:)
    character(len=*),                  intent(in)    :: requestor_name

    integer :: var_id, i_dim
    integer, allocatable :: dim_sizes(:)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 198190218, this%is_open( ),                              &
                     "Trying to read from unopen file: "//id_str )
    var_id    = this%variable_id( variable_name )
    dim_sizes = this%dimension_sizes( variable_name )
    call assert_msg( 650822371, size( dim_sizes ) .eq. 4,                     &
                     "Wrong number of dimensions for "//                      &
                     trim( id_str%to_char( ) )//": Expected 4 got "//         &
                     trim( to_char( size( dim_sizes ) ) ) )
    if( allocated( container ) ) then
      do i_dim = 1, 4
        call assert_msg( 820979275, size( container, i_dim ) .eq.             &
                                    dim_sizes( i_dim ),                       &
                         "Wrong size container for "//                        &
                         trim( id_str%to_char( ) )//": Expected "//           &
                         trim( to_char( dim_sizes( i_dim ) ) )//     &
                         " got "//trim( to_char( size( container, i_dim ) ) ) &
                         //" for dimension "//trim( to_char( i_dim ) ) )
      end do
    else
      allocate( container( dim_sizes(1), dim_sizes(2), dim_sizes(3),          &
                           dim_sizes(4) ) )
    end if
    call check_status( 708660930,                                             &
                       nf90_get_var( this%file_id_, var_id, container ),      &
                       "Error getting value for "//trim( id_str%to_char( ) ) )

  end subroutine read_4D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 0D integer data
  subroutine read_0D_int( this, variable_name, container, requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : musica_dk
    use musica_string,                 only : to_char
    use netcdf,                        only : nf90_get_var

    class(io_netcdf_t), intent(inout) :: this
    class(string_t),    intent(in)    :: variable_name
    integer,            intent(out)   :: container
    character(len=*),   intent(in)    :: requestor_name

    integer :: var_id
    integer, allocatable :: dim_sizes(:)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 418014896, this%is_open( ),                              &
                     "Trying to read from unopen file: "//id_str )
    var_id    = this%variable_id( variable_name )
    dim_sizes = this%dimension_sizes( variable_name )
    call assert_msg( 747800090, size( dim_sizes ) .eq. 0,                     &
                     "Wrong number of dimensions for "//                      &
                     trim( id_str%to_char( ) )//": Expected 0 got "//         &
                     trim( to_char( size( dim_sizes ) ) ) )
    call check_status( 860118435,                                             &
                       nf90_get_var( this%file_id_, var_id, container ),      &
                       "Error getting value for "//trim( id_str%to_char( ) ) )

  end subroutine read_0D_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 1D integer data
  subroutine read_1D_int( this, variable_name, container, requestor_name )

    use musica_assert,                 only : assert_msg
    use musica_constants,              only : musica_dk
    use musica_string,                 only : to_char
    use netcdf,                        only : nf90_get_var

    class(io_netcdf_t),   intent(inout) :: this
    class(string_t),      intent(in)    :: variable_name
    integer, allocatable, intent(out)   :: container(:)
    character(len=*),     intent(in)    :: requestor_name

    integer :: var_id
    integer, allocatable :: dim_sizes(:)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 121652260, this%is_open( ),                              &
                     "Trying to read from unopen file: "//id_str )
    var_id    = this%variable_id( variable_name )
    dim_sizes = this%dimension_sizes( variable_name )
    call assert_msg( 798921103, size( dim_sizes ) .eq. 1,                     &
                     "Wrong number of dimensions for "//                      &
                     trim( id_str%to_char( ) )//": Expected 1 got "//         &
                     trim( to_char( size( dim_sizes ) ) ) )
    if( allocated( container ) ) then
      call assert_msg( 346288950, size( container ) .eq. dim_sizes(1),        &
                       "Wrong size container for "//trim( id_str%to_char( ) ) &
                       //": Expected "//trim( to_char( dim_sizes(1) ) )//     &
                       " got "//trim( to_char( size( container ) ) ) )
    else
      allocate( container( dim_sizes(1) ) )
    end if
    call check_status( 458607295,                                             &
                       nf90_get_var( this%file_id_, var_id, container ),      &
                       "Error getting value for "//trim( id_str%to_char( ) ) )

  end subroutine read_1D_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the dimension names for a given variable
  function variable_dimensions( this, variable_name, requestor_name )         &
      result( dimensions )

    use musica_string,                 only : to_char
    use netcdf,                        only : NF90_MAX_NAME,                  &
                                              nf90_inquire_variable,          &
                                              nf90_inquire_dimension

    type(string_t),     allocatable :: dimensions(:)
    class(io_netcdf_t), intent(in)  :: this
    class(string_t),    intent(in)  :: variable_name
    character(len=*),   intent(in)  :: requestor_name

    integer :: var_id, i_dim, n_dims
    integer, allocatable :: dimids(:)
    type(string_t) :: id_str
    character(len=NF90_MAX_NAME) :: dim_name

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    var_id = this%variable_id( variable_name )
    call check_status( 744311319,                                             &
        nf90_inquire_variable( this%file_id_, var_id, ndims = n_dims ),       &
        "Error getting number of dimensions for "//id_str%to_char( ) )
    allocate( dimids( n_dims ) )
    call check_status( 104014576,                                             &
        nf90_inquire_variable( this%file_id_, var_id, dimids = dimids ),      &
        "Error getting dimesions for "//id_str%to_char( ) )
    allocate( dimensions( n_dims ) )
    do i_dim = 1, n_dims
      call check_status( 788714786,                                           &
          nf90_inquire_dimension( this%file_id_, dimids( i_dim ),             &
                                  name = dim_name ),&
          "Error getting dimesion size "//trim( to_char( i_dim ) )//" for "// &
          id_str%to_char( ) )
      dimensions( i_dim ) = trim( dim_name )
    end do

  end function variable_dimensions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the units for a given variable
  function variable_units( this, variable_name, requestor_name )

    use musica_string,                 only : to_char
    use netcdf,                        only : NF90_MAX_NAME,                  &
                                              nf90_get_att

    type(string_t)                 :: variable_units
    class(io_netcdf_t), intent(in) :: this
    class(string_t),    intent(in) :: variable_name
    character(len=*),   intent(in) :: requestor_name

    integer :: var_id
    type(string_t) :: id_str
    character(len=NF90_MAX_NAME) :: units

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    var_id = this%variable_id( variable_name )
    call check_status( 301987512,                                             &
                       nf90_get_att( this%file_id_, var_id, "units", units ), &
                       "Error getting units for "//trim( id_str%to_char( ) ) )
    variable_units = trim( units )

  end function variable_units

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns whether a file is open or not
  logical function is_open( this )

    class(io_netcdf_t), intent(in) :: this

    is_open = this%file_id_ .ne. kUnknownFileId

  end function is_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns a variable's id in the NetCDF file
  integer function variable_id( this, variable_name )

    use musica_assert,                 only : assert_msg
    use netcdf,                        only : nf90_inq_varid

    class(io_netcdf_t), intent(in) :: this
    class(string_t),    intent(in) :: variable_name

    call assert_msg( 249726322, this%is_open( ),                              &
                     "Trying to read from unopen file" )
    call check_status( 153462424,                                             &
                       nf90_inq_varid( this%file_id_,                         &
                                       variable_name%to_char( ),              &
                                       variable_id ), &
                       "Cannot file variable '"//                             &
                       trim( variable_name%to_char( ) )//"' in file '"//      &
                       trim( this%file_name_%to_char( ) )//"'" )

  end function variable_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the dimensions for variable in the NetCDF file
  function dimension_sizes( this, variable_name ) result( dim_sizes )

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : to_char
    use netcdf,                        only : nf90_inquire_variable,          &
                                              nf90_inquire_dimension

    integer,            allocatable :: dim_sizes(:)
    class(io_netcdf_t), intent(in)  :: this
    class(string_t),    intent(in)  :: variable_name

    integer :: var_id, n_dims, i_dim
    integer, allocatable :: dimids(:)
    type(string_t) :: id_str

    id_str = "variable '"//variable_name//"' in file '"//this%file_name_//"'"
    call assert_msg( 191887763, this%is_open( ),                              &
                     "Trying to read from unopen file: "//id_str )
    var_id = this%variable_id( variable_name )
    call check_status( 516121527,                                             &
        nf90_inquire_variable( this%file_id_, var_id, ndims = n_dims ),       &
        "Error getting number of dimensions for "//trim( id_str%to_char( ) ) )
    allocate( dimids( n_dims ) )
    call check_status( 269878960,                                             &
        nf90_inquire_variable( this%file_id_, var_id, dimids = dimids ),      &
        "Error getting dimensions for "//trim( id_str%to_char( ) ) )
    allocate( dim_sizes( n_dims ) )
    do i_dim = 1, n_dims
      call check_status( 770273353,                                           &
          nf90_inquire_dimension( this%file_id_, dimids( i_dim ),             &
                                  len = dim_sizes( i_dim ) ),                 &
          "Error getting dimension size "//trim( to_char( i_dim ) )//" for "//&
          trim( id_str%to_char( ) ) )
    end do

  end function dimension_sizes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalizes a NetCDF file reader
  subroutine finalize( this )

    use netcdf,                        only : nf90_close

    type(io_netcdf_t), intent(inout) :: this

    if( this%file_id_ .ne. kUnknownFileId ) then
      call check_status( 708311006, nf90_close( this%file_id_ ),              &
                         "Error closing file" )
    end if
    this%file_id_   = kUnknownFileId
    this%file_name_ = ""

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! @name Private NetCDF support functions
!! @{
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Checks a NetCDF status code and fail with a message if an error occurred
  subroutine check_status( code, status, error_message )

    use musica_assert,                 only : die_msg
    use netcdf,                        only : NF90_NOERR, nf90_strerror

    !> Unique code to associate with any failure
    integer,          intent(in) :: code
    !> NetCDF status code
    integer,          intent(in) :: status
    !> Error message to display on failure
    character(len=*), intent(in) :: error_message

    if( status .eq. NF90_NOERR ) return
    call die_msg( 330311277, "NetCDF error: "//trim( error_message )//": "//  &
                  trim( nf90_strerror( status ) ) )

  end subroutine check_status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @}

end module musica_io_netcdf
