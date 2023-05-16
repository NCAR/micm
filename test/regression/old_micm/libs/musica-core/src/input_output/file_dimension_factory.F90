! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> File dimension factory

!> Builder of file dimension objects
module musica_file_dimension_factory

  use musica_file_dimension,           only : file_dimension_t
  use musica_file_dimension_netcdf,    only : file_dimension_netcdf_t
  use musica_file_dimension_text,      only : file_dimension_text_t

  implicit none
  private

  public :: file_dimension_builder

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Builds an file dimension object for a specific file type
  !!
  !! The \c config argument must include a top-level key-value pair "type"
  !! whose value is a valid input/output type name, Currently, these are:
  !! - "csv", "txt", or "text" (a delimited text file)
  !! - "nc" or "netcdf" (a NetCDF file)
  !!
  function file_dimension_builder( config, file, variable, dimension_name )   &
      result( new_dimension )

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_file,                   only : file_t
    use musica_file_variable,          only : file_variable_t
    use musica_string,                 only : string_t

    !> New file dimension
    class(file_dimension_t), pointer :: new_dimension
    !> Dimension configuration
    type(config_t), intent(inout) :: config
    !> Input/Output file
    class(file_t), intent(inout) :: file
    !> File variable associated with the dimension
    class(file_variable_t), intent(in), optional :: variable
    !> Name of the dimension
    !! (This should be provided when the dimension has indices but no values)
    character(len=*), intent(in), optional :: dimension_name

    type(string_t) :: file_type
    character(len=*), parameter :: my_name = 'File dimension builder'

    new_dimension => null( )
    call config%get( 'type', file_type, my_name )
    file_type = file_type%to_lower( )

    if( file_type .eq. 'txt' .or.                                             &
        file_type .eq. 'text' .or.                                            &
        file_type .eq. 'csv' ) then
      new_dimension => file_dimension_text_t( file, variable, dimension_name )
    else if( file_type .eq. 'nc' .or.                                         &
             file_type .eq. 'netcdf' ) then
      new_dimension => file_dimension_netcdf_t( file, variable,               &
                                                dimension_name )
    else
      call die_msg( 178041200, "Invalid input/output file type: '"//          &
                               file_type%to_char( )//"'" )
    end if

  end function file_dimension_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_file_dimension_factory
