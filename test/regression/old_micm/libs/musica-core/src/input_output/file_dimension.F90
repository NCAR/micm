! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_file_dimension module

!> The file_dimension_t type and related functions
module musica_file_dimension

  use musica_constants,                only : musica_dk, musica_ik
  use musica_file_dimension_range,     only : file_dimension_range_t
  use musica_file_variable,            only : file_variable_t
  use musica_string,                   only : string_t

  implicit none
  private

  public :: file_dimension_t

  !> A File dimension
  !!
  !! File dimensions are used to define the range and boundaries of
  !! file variables.
  !!
  type, abstract :: file_dimension_t
    private
    !> Dimension name
    type(string_t) :: name_
    !> Dimension size
    type(file_dimension_range_t) :: range_
    !> All dimension values present in the file
    real(kind=musica_dk), allocatable :: values_(:)
    !> File variable for this dimension
    class(file_variable_t), pointer :: variable_ => null( )
  contains
    !> Name of the dimension
    procedure :: name => get_name
    !> Gets the values for this dimension
    procedure :: get_values
    !> Gets the index for a given dimesion value
    procedure :: get_index
    !> Gets the dimension range
    procedure :: get_range
    !> Prints the properties of the dimension
    procedure :: print => do_print
    !> Private constructor
    !! (Should only be called by constructors of extending types)
    procedure :: private_constructor
    !> Finalizes the object
    !! (Should only be called by final procedures of extending types)
    procedure :: private_finalize
  end type file_dimension_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gets the name of the dimension
  function get_name( this )

    !> Dimension name
    type(string_t) :: get_name
    !> File dimension
    class(file_dimension_t), intent(in) :: this

    get_name = this%name_

  end function get_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gets the values for this dimension
  function get_values( this ) result( values )

    !> Dimension values
    real(kind=musica_dk), allocatable :: values(:)
    !> File dimension
    class(file_dimension_t), intent(in) :: this

    values = this%values_

  end function get_values

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gets the index for a given dimension value
  !!
  !! If the value does not correspond exactly to a file index, the index of
  !! closest value less the requested value is returned. The is_exact flag
  !! can be included to indicate whether an exact match was found.
  !!
  !! If the first dimension value is greater than the requested value, an
  !! index of 1 is returned.
  !!
  function get_index( this, value, is_exact, guess ) result( index )

    use musica_assert,                 only : assert
    use musica_config,                 only : config_t
    use musica_datetime,               only : datetime_t

    !> Index for the closest (without going over) value
    integer(kind=musica_ik) :: index
    !> File dimension
    class(file_dimension_t), intent(in) :: this
    !> Value to find
    real(kind=musica_dk), intent(in) :: value
    !> Flag indicating whether an exact match was found
    logical, intent(out), optional :: is_exact
    !> A guess for the index that will be used to start the search
    integer(kind=musica_ik), intent(in), optional :: guess

    integer(kind=musica_ik) :: i_val, l_guess

    if( present( is_exact ) ) is_exact = .false.
    if( this%values_( 1 ) .gt. value ) then
      index = 1
      return
    end if
    if( present( guess ) ) then
      l_guess = guess
      if( l_guess .gt. size( this%values_ ) ) l_guess = size( this%values_ )
      if( l_guess .lt. 1 ) l_guess = 1
      if( this%values_( l_guess ) .le. value ) then
        do i_val = l_guess, size( this%values_ )
          if( this%values_( i_val ) .eq. value ) then
            if( present( is_exact ) ) is_exact = .true.
            index = i_val
            return
          else if( this%values_( i_val ) .gt. value ) then
            index = i_val - 1
            return
          end if
        end do
      else
        do i_val = l_guess, 0, -1
          if( this%values_( i_val ) .eq. value ) then
            if( present( is_exact ) ) is_exact = .true.
            index = i_val
            return
          else if( this%values_( i_val ) .lt. value ) then
            index = i_val
            return
          end if
        end do
      end if
    else
      do i_val = 1, size( this%values_ )
        if( this%values_( i_val ) .eq. value ) then
          if( present( is_exact ) ) is_exact = .true.
          index = i_val
          return
        else if( this%values_( i_val ) .gt. value ) then
          index = i_val - 1
          return
        end if
      end do
    end if

  end function get_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the range for this dimension
  function get_range( this )

    use musica_assert,                 only : assert

    !> Dimension range
    type(file_dimension_range_t) :: get_range
    !> Dimension
    class(file_dimension_t), intent(in) :: this

    get_range = this%range_

  end function get_range

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Prints the properties of the dimension
  subroutine do_print( this )

    !> File dimension
    class(file_dimension_t), intent(in) :: this

    write(*,*) "*** Dimension ***"
    if( associated( this%variable_ ) ) call this%variable_%print( )
    write(*,*) "*** End Dimension ***"

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalizes the object
  !! (Should only be called by final procedures of extending types)
  elemental subroutine private_finalize( this )

    !> File dimension
    class(file_dimension_t), intent(inout) :: this

    if( associated( this%variable_ ) ) then
      deallocate( this%variable_ )
      this%variable_ => null( )
    end if

  end subroutine private_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Private constructor for common data elements
  !! (Should only be called by constructors of extending types)
  subroutine private_constructor( this, dimension_name, file, range, variable )

    use musica_assert,                 only : assert_msg
    use musica_file,                   only : file_t

    !> File dimension
    class(file_dimension_t), intent(inout) :: this
    !> Dimension name
    type(string_t), intent(in) :: dimension_name
    !> Input file
    class(file_t), intent(inout) :: file
    !> Dimension boundaries
    type(file_dimension_range_t), intent(in) :: range
    !> File variable associated with this dimension
    class(file_variable_t), intent(in), optional :: variable

    type(file_dimension_range_t) :: ranges(1)

    this%name_  = dimension_name
    this%range_ = range
    if( present( variable ) ) then
      allocate( this%variable_, source = variable )
      allocate( this%values_( range%lower_bound( ) : range%upper_bound( ) ) )
      ranges(1) = range
      call this%variable_%get_data( file, ranges, this%values_ )
    end if

  end subroutine private_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_file_dimension
