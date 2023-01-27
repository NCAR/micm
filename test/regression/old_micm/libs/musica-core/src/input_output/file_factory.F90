! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Input/output file factory

!> Builder of input/output files
module musica_file_factory

  use musica_file,                     only : file_t
  use musica_file_netcdf,              only : file_netcdf_t
  use musica_file_text,                only : file_text_t

  implicit none
  private

  public :: file_builder

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Builds an input/output file object by type
  !!
  !! The \c config argument must include a top-level key-value pair "type"
  !! whose value is a valid input/output type name, Currently, these are:
  !! - "csv", "txt", or "text" (a delimited text file)
  !! - "nc" or "netcdf" (a NetCDF file)
  !!
  !! The \c config must also include a top-level key-value pair "intent",
  !! which can be either "input" or "output".
  !!
  !! A "file name" is also required for files openned for input. This is
  !! optional for output files, with the default name starting with "output"
  !! and having an appropriate file extension.
  !!
  function file_builder( config ) result( new_file )

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    !> New input/output file object
    class(file_t), pointer :: new_file
    !> Input/output file configuration
    type(config_t), intent(inout) :: config

    type(string_t) :: file_type
    character(len=*), parameter :: my_name = 'File builder'

    new_file => null( )
    call config%get( 'type', file_type, my_name )
    file_type = file_type%to_lower( )

    if( file_type .eq. 'txt' .or.                                             &
        file_type .eq. 'text' .or.                                            &
        file_type .eq. 'csv' ) then
      new_file => file_text_t( config )
    else if( file_type .eq. 'nc' .or.                                         &
             file_type .eq. 'netcdf' ) then
      new_file => file_netcdf_t( config )
    else
      call die_msg( 690482624, "Invalid input/output file type: '"//          &
                               file_type%to_char( )//"'" )
    end if

  end function file_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_file_factory
