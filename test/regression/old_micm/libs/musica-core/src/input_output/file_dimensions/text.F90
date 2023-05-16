! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_file_dimension_text module

!> The file_dimension_text_t type and related functions
module musica_file_dimension_text

  use musica_constants,                only : musica_dk, musica_ik
  use musica_file_dimension,           only : file_dimension_t

  implicit none
  private

  public :: file_dimension_text_t

  !> A text file dimension
  type, extends(file_dimension_t) :: file_dimension_text_t
    private
    !> Text dimension id
    integer(kind=musica_ik) :: id_ = -1
  contains
    !> Finalizes the dimension
    final :: finalize
  end type file_dimension_text_t

  !> Constructor
  interface file_dimension_text_t
    module procedure :: constructor
  end interface file_dimension_text_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates a file_dimension_text_t object for a text file dimension
  function constructor( file, variable, dimension_name ) result( new_obj )

    use musica_assert,                 only : assert, die
    use musica_file,                   only : file_t
    use musica_file_dimension_range,   only : file_dimension_range_t
    use musica_file_text,              only : file_text_t
    use musica_file_variable,          only : file_variable_t
    use musica_string,                 only : string_t

    !> Pointer to the new text file dimension object
    class(file_dimension_t), pointer :: new_obj
    !> Text file
    class(file_t), intent(inout) :: file
    !> Text file variable associated with the dimension
    class(file_variable_t), intent(in), optional :: variable
    !> Name of the dimension
    !! (This should be provided when the dimension has indices but no values)
    character(len=*), intent(in), optional :: dimension_name

    type(file_dimension_range_t) :: range
    integer(kind=musica_ik) :: n_values
    type(string_t) :: dim_name

    ! there is only a time dimension in text files
    if( present( variable ) ) dim_name = variable%dimension_name( 1 )
    if( present( dimension_name ) ) dim_name = dimension_name
    call assert( 103261211, dim_name .ne. "" )
    allocate( file_dimension_text_t :: new_obj )
    select type( new_obj )
    class is( file_dimension_text_t )
      select type( file )
      class is( file_text_t )
        new_obj%id_ = file%get_dimension_id(   dim_name    )
        n_values    = file%get_dimension_size( new_obj%id_ )
        range = file_dimension_range_t( new_obj%id_, 1, n_values,             &
                                        is_unlimited = .true. )
      class default
        call die( 361184159 )
      end select
    class default
      call die( 585820849 )
    end select
    call new_obj%private_constructor( dim_name, file, range, variable )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalizes a text file dimension
  elemental subroutine finalize( this )

    !> Text file dimension
    type(file_dimension_text_t), intent(inout) :: this

    call this%private_finalize( )

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_file_dimension_text
