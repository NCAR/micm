! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> File variable factory

!> Builder of file variable objects
module musica_file_variable_factory

  use musica_constants,               only : musica_dk, musica_ik
  use musica_file_variable,           only : file_variable_t
  use musica_file_variable_netcdf,    only : file_variable_netcdf_t
  use musica_file_variable_text,      only : file_variable_text_t

  implicit none
  private

  public :: file_variable_builder

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Builds an file variable object for a specific file type
  !!
  !! The \c config argument must include a top-level key-value pair "type"
  !! whose value is a valid input/output type name, Currently, these are:
  !! - "csv", "txt", or "text" (a delimited text file)
  !! - "nc" or "netcdf" (a NetCDF file)
  !!
  function file_variable_builder( config, file, variable_name, variable_id,   &
      dimensions, found ) result( new_variable )

    use musica_assert,                 only : assert_msg, die, die_msg
    use musica_config,                 only : config_t
    use musica_file,                   only : file_t
    use musica_file_dimension_range,   only : file_dimension_range_t
    use musica_string,                 only : string_t

    !> New file variable
    class(file_variable_t), pointer :: new_variable
    !> Dimension configuration
    type(config_t), intent(inout) :: config
    !> Input/Output file
    class(file_t), intent(inout) :: file
    !> Variable name
    character(len=*), intent(in), optional :: variable_name
    !> Variable ID
    integer(kind=musica_ik), intent(in), optional :: variable_id
    !> Variable dimensions
    !!
    !! Only required when adding a new file variable
    !!
    type(file_dimension_range_t), intent(in), optional :: dimensions(:)
    !> Optional flag that indicates whether the variable was found in the
    !! file
    logical, intent(out), optional :: found

    type(string_t) :: file_type
    character(len=*), parameter :: my_name = 'File variable builder'

    if( present( found ) ) found = .true.
    call assert_msg( 379618774,                                               &
                     present( variable_name ) .neqv. present( variable_id ),  &
                     "Either a variable name or a variable id must be "//     &
                     "specified to build a file variable object "//           &
                     "(but not both)" )
    new_variable => null( )
    call config%get( 'type', file_type, my_name )
    file_type = file_type%to_lower( )

    if( file_type .eq. 'txt' .or.                                             &
        file_type .eq. 'text' .or.                                            &
        file_type .eq. 'csv' ) then
      if( present( variable_name ) ) then
        new_variable => file_variable_text_t( file, variable_name, dimensions,&
                                              config, found )
      else if( present( variable_id ) ) then
        new_variable => file_variable_text_t( file, variable_id, config )
      else
        call die( 284671798 )
      end if
    else if( file_type .eq. 'nc' .or.                                         &
             file_type .eq. 'netcdf' ) then
      if( present( variable_name ) ) then
        new_variable => file_variable_netcdf_t( file, variable_name,          &
                                                dimensions, config, found )
      else if( present( variable_id ) ) then
        new_variable => file_variable_netcdf_t( file, variable_id, config )
      else
        call die( 788424964 )
      end if
    else
      call die_msg( 192668005, "Invalid input/output file type: '"//          &
                               file_type%to_char( )//"'" )
    end if

  end function file_variable_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_file_variable_factory
