! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_file_dimension_netcdf module

!> The file_dimension_netcdf_t type and related functions
module musica_file_dimension_netcdf

  use musica_constants,                only : musica_dk, musica_ik
  use musica_file_dimension,           only : file_dimension_t

  implicit none
  private

  public :: file_dimension_netcdf_t

  !> A NetCDF dimension
  type, extends(file_dimension_t) :: file_dimension_netcdf_t
    private
    !> NetCDF dimension id
    integer(kind=musica_ik) :: id_ = -1
  contains
    !> Gets the number of values associated with a dimension
    procedure, private :: number_of_values
    !> Finalizes the dimension
    final :: finalize
  end type file_dimension_netcdf_t

  !> Constructor
  interface file_dimension_netcdf_t
    module procedure :: constructor
  end interface file_dimension_netcdf_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates a file_dimension_netcdf_t object for a NetCDF dimension
  function constructor( file, variable, dimension_name ) result( new_obj )

    use musica_assert,                 only : assert, die
    use musica_file,                   only : file_t
    use musica_file_dimension_range,   only : file_dimension_range_t
    use musica_file_netcdf,            only : file_netcdf_t
    use musica_file_variable,          only : file_variable_t
    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_inq_dimid, nf90_inquire

    !> Pointer to the new NetCDF dimension object
    class(file_dimension_t), pointer :: new_obj
    !> NetCDF file
    class(file_t), intent(inout) :: file
    !> NetCDF variable associated with the dimension
    class(file_variable_t), intent(in), optional :: variable
    !> Name of the dimension
    !! (This should be provided when the dimension has indices but no values)
    character(len=*), intent(in), optional :: dimension_name

    integer(kind=musica_ik) :: unlimited_dimid
    type(string_t) :: dim_name
    type(file_dimension_range_t) :: range
    integer(kind=musica_ik) :: n_values

    if( present( variable ) ) dim_name = variable%dimension_name( 1 )
    if( present( dimension_name ) ) dim_name = dimension_name
    call assert( 855947264, dim_name .ne. "" )
    allocate( file_dimension_netcdf_t :: new_obj )
    select type( new_obj )
    class is( file_dimension_netcdf_t )
      select type( file )
      class is( file_netcdf_t )
        call file%check_open( )
        call file%check_status( 140723118,                                    &
            nf90_inq_dimid( file%id( ), dim_name%to_char( ), new_obj%id_ ),   &
            "Error finding id for dimension '"//dim_name%to_char( )//"'" )
        call file%check_status( 426925310,                                    &
            nf90_inquire( file%id( ), unlimitedDimid = unlimited_dimid ),     &
            "Error getting unlimited dimension id." )
        n_values = new_obj%number_of_values( file )
        range = file_dimension_range_t( new_obj%id_, 1, n_values,             &
                              is_unlimited = new_obj%id_ .eq. unlimited_dimid )
      class default
        call die( 345843533 )
      end select
    class default
      call die( 353013374 )
    end select
    call new_obj%private_constructor( dim_name, file, range, variable )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Prints the properties of the dimension
  subroutine do_print( this )

    use musica_string,                 only : to_char

    !> NetCDF dimension
    class(file_dimension_netcdf_t), intent(in) :: this

    write(*,*) "*** Dimension id: "//to_char( this%id_ )//" ***"
    write(*,*) "*** End Dimension ***"

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the numnber of values for the dimension from the NetCDF file
  integer(kind=musica_ik) function number_of_values( this, file )             &
      result( n_values )

    use musica_assert,                 only : die
    use musica_file,                   only : file_t
    use musica_file_netcdf,            only : file_netcdf_t
    use musica_file_variable,          only : file_variable_t
    use musica_string,                 only : string_t
    use netcdf,                        only : nf90_inquire_dimension

    !> NetCDF dimension
    class(file_dimension_netcdf_t), intent(inout) :: this
    !> NetCDF file
    class(file_t), intent(inout) :: file

    type(string_t) :: dim_name

    select type( file )
      class is( file_netcdf_t )
        dim_name = this%name( )
        call file%check_status( 649288296,                                    &
                                nf90_inquire_dimension( file%id( ),           &
                                                        this%id_,             &
                                                        len = n_values ),     &
                                "Error getting values for dimension '"//      &
                                dim_name%to_char( )//"'" )
      class default
        call die( 775240150 )
    end select

  end function number_of_values

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalizes a NetCDF dimension
  elemental subroutine finalize( this )

    !> NetCDF dimension
    type(file_dimension_netcdf_t), intent(inout) :: this

    call this%private_finalize( )

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_file_dimension_netcdf
