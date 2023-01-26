! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_file_dimension_range module

!> The file_dimension_range_t type and related functions
module musica_file_dimension_range

  use musica_constants,                only : musica_dk, musica_ik

  implicit none
  private

  public :: file_dimension_range_t

  !> A file dimension range
  !!
  !! File dimension ranges can be used alone or in combination to identify
  !! subsets of file variable data for input or output.
  !!
  !! \todo add example usage for file_dimension_range_t
  !!
  type :: file_dimension_range_t
    private
    !> Id of the dimension in the file
    integer(kind=musica_ik) :: dimension_id_ = -1
    !> Lower dimension index bound (inclusive)
    integer(kind=musica_ik) :: lower_bound_ = -1
    !> Upper dimension index bound (inclusive)
    integer(kind=musica_ik) :: upper_bound_ = -1
    !> Flag indicating whether the dimension is unlimited in length
    logical :: is_unlimited_ = .false.
  contains
    !> Sets the range
    procedure :: set
    !> Get the dimension id
    procedure :: id
    !> Gets the lower index bound (inclusive)
    procedure :: lower_bound
    !> Gets the upper index bound (inclusive)
    procedure :: upper_bound
    !> Returns a flag indicating whether the dimension is unlimited in length
    procedure :: is_unlimited
    !> @name Comparison operators (based only on dimension id)
    !> @{
    procedure :: equals
    generic :: operator(==) => equals
    procedure :: not_equals
    generic :: operator(/=) => not_equals
    !> @}
    !> Print the range parameters
    procedure :: print => do_print
  end type file_dimension_range_t

  !> Constructor for file dimensions ranges
  interface file_dimension_range_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for file_dimension_range_t objects
  function constructor( dimension_id, lower_bound, upper_bound, is_unlimited )&
      result( new_obj )

    !> File dimension range
    type(file_dimension_range_t) :: new_obj
    !> Dimension
    integer(kind=musica_ik), intent(in) :: dimension_id
    !> Lower index bound (inclusive)
    integer(kind=musica_ik), intent(in) :: lower_bound
    !> Upper index bound (inclusive)
    integer(kind=musica_ik), intent(in) :: upper_bound
    !> Flag indicating whether the dimension is unlimited in length
    !!
    !! Defaults to false
    logical, intent(in), optional :: is_unlimited

    new_obj%dimension_id_ = dimension_id
    new_obj%lower_bound_  = lower_bound
    new_obj%upper_bound_  = upper_bound
    if( present( is_unlimited ) ) new_obj%is_unlimited_ = is_unlimited

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets the range
  subroutine set( this, lower_bound, upper_bound )

    !> File dimension range
    class(file_dimension_range_t), intent(inout) :: this
    !> Lower index bound (inclusive)
    integer(kind=musica_ik), intent(in) :: lower_bound
    !> Upper index bound (inclusive)
    integer(kind=musica_ik), intent(in) :: upper_bound

    this%lower_bound_  = lower_bound
    this%upper_bound_  = upper_bound

  end subroutine set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the dimension id
  elemental integer(kind=musica_ik) function id( this )

    !> File dimension range
    class(file_dimension_range_t), intent(in) :: this

    id = this%dimension_id_

  end function id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gets the lower index bound (inclusive)
  elemental integer(kind=musica_ik) function lower_bound( this )

    !> File dimension range
    class(file_dimension_range_t), intent(in) :: this

    lower_bound = this%lower_bound_

  end function lower_bound

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gets the upper index bound (inclusive)
  elemental integer(kind=musica_ik) function upper_bound( this )

    !> File dimension range
    class(file_dimension_range_t), intent(in) :: this

    upper_bound = this%upper_bound_

  end function upper_bound

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns a flag indicating whether the dimension is unlimited in length
  logical function is_unlimited( this )

    !> File dimension range
    class(file_dimension_range_t), intent(in) :: this

    is_unlimited = this%is_unlimited_

  end function is_unlimited

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Equality comparison
  logical elemental function equals( a, b )

    !> Range a
    class(file_dimension_range_t), intent(in) :: a
    !> Range b
    class(file_dimension_range_t), intent(in) :: b

    equals = a%dimension_id_ .eq. b%dimension_id_

  end function equals

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Inequality comparison
  logical elemental function not_equals( a, b )

    !> Range a
    class(file_dimension_range_t), intent(in) :: a
    !> Range b
    class(file_dimension_range_t), intent(in) :: b

    not_equals = a%dimension_id_ .ne. b%dimension_id_

  end function not_equals

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Prints the range parameters
  subroutine do_print( this )

    !> File dimension range
    class(file_dimension_range_t), intent(in) :: this

    write(*,*) "*** File dimension range ***"
    write(*,*) "Dimension id: ", this%dimension_id_
    write(*,*) "Lower bound:  ", this%lower_bound_
    write(*,*) "Upper bound:  ", this%upper_bound_
    write(*,*) "Is unlimited? ", this%is_unlimited_
    write(*,*) "*** End file dimension range ***"

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_file_dimension_range
