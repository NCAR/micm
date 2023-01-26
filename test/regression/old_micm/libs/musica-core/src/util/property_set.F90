! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_property_set module

!> The property_set_t type and related functions
module musica_property_set

  use musica_constants,                only : musica_dk, musica_ik
  use musica_property,                 only : property_ptr
  use musica_string,                   only : string_t

  implicit none
  private

  public :: property_set_t

  !> A set of physical or other properties
  !!
  !! \todo add detailed description and example of property_set_t usage
  !!
  type :: property_set_t
    private
    !> Properties
    type(property_ptr), allocatable :: properties_(:)
  contains
    !> @name Assignment
    !! @{
    procedure :: assign_property_set
    generic :: assignment(=) => assign_property_set
    !> @}
    !> Returns the index of a property in the set
    procedure :: index => property_set_index
    !> Returns the number of properties in the set
    procedure :: size => property_set_size
    !> Adds a property to the set
    procedure :: add
    !> Returns a subset of the properties
    procedure :: subset
    !> @name Returns a property
    !! @{
    procedure, private :: get_by_index
    procedure, private :: get_by_name
    generic :: get => get_by_index, get_by_name
    !> @}
    !> Output the property set to a given file unit
    procedure :: output
  end type property_set_t

  !> Property set constructor
  interface property_set_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates a new property set
  function constructor( ) result( new_obj )

    !> New property set
    type(property_set_t) :: new_obj

    allocate( new_obj%properties_( 0 ) )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Assigns a property set from an existing property set
  subroutine assign_property_set( to, from )

    !> Property set to assign
    class(property_set_t), intent(inout) :: to
    !> Property set to copy
    type(property_set_t), intent(in) :: from

    integer :: i_prop

    if( allocated( to%properties_ ) ) deallocate( to%properties_ )
    if( .not. allocated( from%properties_ ) ) return
    allocate( to%properties_( size( from%properties_ ) ) )
    do i_prop = 1, size( from%properties_ )
      allocate( to%properties_( i_prop )%val_ )
      to%properties_( i_prop )%val_ = from%properties_( i_prop )%val_
    end do

  end subroutine assign_property_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the index of a property in the set
  integer(kind=musica_ik) function property_set_index( this, property )

    use musica_assert,                 only : assert, die_msg
    use musica_data_type,              only : data_type_t
    use musica_property,               only : property_t

    !> Property set
    class(property_set_t), intent(in) :: this
    !> Property to find
    class(property_t), intent(in) :: property

    integer(kind=musica_ik) :: i_prop

    property_set_index = 0
    call assert( 388666315, allocated( this%properties_ ) )
    do i_prop = 1, size( this%properties_ )
      property_set_index = property_set_index + 1
      if( this%properties_( i_prop )%val_%name( ) .eq. property%name( ) ) then
        call this%properties_( i_prop )%val_%must_equal( property )
        return
      end if
    end do
    call die_msg( 732626019, "Property set does not include property '"//     &
                  property%name( )//"' defined by '"//property%defined_by( ) )

  end function property_set_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of properties in the set
  integer(kind=musica_ik) function property_set_size( this )

    use musica_assert,                 only : assert
    use musica_data_type,              only : data_type_t

    !> Property set
    class(property_set_t), intent(in) :: this

    integer(kind=musica_ik) :: i_prop

    property_set_size = 0
    call assert( 435740314, allocated( this%properties_ ) )
    do i_prop = 1, size( this%properties_ )
      property_set_size = property_set_size + 1
    end do

  end function property_set_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adds a property to the set
  subroutine add( this, new_property )

    use musica_assert,                 only : assert
    use musica_property,               only : property_t

    !> Property set
    class(property_set_t), intent(inout) :: this
    !> Property to add
    class(property_t), intent(in) :: new_property

    integer(kind=musica_ik) :: i_prop, n_elem
    type(property_ptr), allocatable :: temp_array(:)

    call assert( 841011328, allocated( this%properties_ ) )
    n_elem = size( this%properties_ )
    allocate( temp_array( n_elem ) )
    do i_prop = 1, n_elem
      allocate( temp_array( i_prop )%val_,                                    &
                mold = this%properties_( i_prop )%val_ )
      temp_array( i_prop )%val_ = this%properties_( i_prop )%val_
      deallocate( this%properties_( i_prop )%val_ )
    end do
    deallocate( this%properties_ )
    allocate( this%properties_( n_elem + 1 ) )
    do i_prop = 1, n_elem
      allocate( this%properties_( i_prop )%val_,                              &
                mold = temp_array( i_prop )%val_ )
      this%properties_( i_prop )%val_ = temp_array( i_prop )%val_
      deallocate( temp_array( i_prop )%val_ )
    end do
    allocate( this%properties_( n_elem + 1 )%val_, mold = new_property )
    this%properties_( n_elem + 1 )%val_ = new_property

  end subroutine add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gets a property by index
  function get_by_index( this, index ) result( prop )

    use musica_assert,                 only : assert
    use musica_property,               only : property_t

    !> Property
    class(property_t), pointer :: prop
    !> Property set
    class(property_set_t), intent(in) :: this
    !> Index of the property in the set
    integer(kind=musica_ik) :: index

    call assert( 971234877, allocated( this%properties_ ) )
    call assert( 295871568, index .ge. 1 .and.                                &
                            index .le. size( this%properties_ ) )
    associate( this_prop => this%properties_( index )%val_ )
    allocate( prop, mold = this_prop )
    prop = this_prop
    end associate

  end function get_by_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gets a property by name
  function get_by_name( this, name, found ) result( prop )

    use musica_assert,                 only : assert, die_msg
    use musica_property,               only : property_t

    !> Property
    class(property_t), pointer :: prop
    !> Property set
    class(property_set_t), intent(in) :: this
    !> Property name
    character(len=*), intent(in) :: name
    !> Flag indicating whether the variable was found
    logical, intent(out), optional :: found

    integer(kind=musica_ik) :: i_prop

    if( present( found ) ) found = .false.
    prop => null( )
    call assert( 774207789, allocated( this%properties_ ) )
    do i_prop = 1, size( this%properties_ )
      if( this%properties_( i_prop )%val_%name( ) .eq. name ) then
        if( present( found ) ) found = .true.
        associate( this_prop => this%properties_( i_prop )%val_ )
        allocate( prop, mold = this_prop )
        prop = this_prop
        end associate
        return
      end if
      end do
    if( .not. present( found ) ) then
      call die_msg( 934288374, "Property '"//name//"' not found in "//        &
                    "property set." )
    end if

  end function get_by_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns a subset of the properties
  !!
  !! Optional parameters act as filters for the returned subset. Passing no
  !! optional parameters results in the full property set being returned.
  !!
  function subset( this, prefix, data_type, applies_to, defined_by )

    use musica_assert,                 only : assert
    use musica_data_type,              only : data_type_t
    use musica_target,                 only : target_t

    !> Property subset
    type(property_set_t) :: subset
    !> Property set
    class(property_set_t), intent(in) :: this
    !> Property name prefix
    character(len=*), intent(in), optional :: prefix
    !> Data type
    class(data_type_t), intent(in), optional :: data_type
    !> Model element(s) to which the properties apply
    class(target_t), intent(in), optional :: applies_to
    !> Model component that defined the properties
    character(len=*), intent(in), optional :: defined_by

    integer(kind=musica_ik) :: i_prop
    type(string_t) :: l_prefix, l_defined_by
    class(target_t), pointer :: prop_target

    if( present( prefix ) ) l_prefix = prefix
    if( present( defined_by ) ) l_defined_by = defined_by
    call assert( 518952648, allocated( this%properties_ ) )
    subset = property_set_t( )
    do i_prop = 1, size( this%properties_ )
      call assert( 240541644, associated( this%properties_( i_prop )%val_ ) )
      associate( prop => this%properties_( i_prop )%val_ )
      if( present( prefix ) ) then
        if( prop%prefix( ) .ne. l_prefix ) cycle
      end if
      if( present( data_type ) ) then
        if( prop%data_type( ) .ne. data_type ) cycle
      end if
      if( present( applies_to ) ) then
        prop_target => prop%applies_to( )
        if( associated( prop_target ) ) then
          if( prop_target .ne. applies_to ) then
            deallocate( prop_target )
            cycle
          endif
          deallocate( prop_target )
        else
          cycle
        end if
      end if
      if( present( defined_by ) ) then
        if( prop%defined_by( ) .ne. l_defined_by ) cycle
      end if
      call subset%add( prop )
      end associate
    end do

  end function subset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Output the contents of the property set to a given file unit
  subroutine output( this, file_unit )

    use musica_assert,                 only : assert
    use musica_data_type,              only : data_type_t
    use musica_string,                 only : output_table
    use musica_target,                 only : target_t

    !> Property set
    class(property_set_t), intent(in) :: this
    !> File unit
    integer(kind=musica_ik), intent(in) :: file_unit

    integer(kind=musica_ik) :: i_row
    type(string_t) :: header(5)
    type(string_t), allocatable :: table(:,:)
    type(data_type_t) :: data_type
    class(target_t), pointer :: applies_to

    call assert( 981304030, allocated( this%properties_ ) )
    allocate( table( 5, size( this%properties_ ) ) )
    header(1) = "Property"
    header(2) = "Units"
    header(3) = "Data Type"
    header(4) = "Applies To"
    header(5) = "Defined By"
    do i_row = 1, size( this%properties_ )
      associate( prop => this%properties_( i_row )%val_ )
      table( 1, i_row ) = prop%name( )
      table( 2, i_row ) = prop%units( )
      data_type = prop%data_type( )
      table( 3, i_row ) = data_type%name( )
      applies_to => prop%applies_to( )
      if( associated( applies_to ) ) then
        table( 4, i_row ) = applies_to%name( )
        deallocate( applies_to )
      else
        table( 4, i_row ) = "undefined"
      end if
      table( 5, i_row ) = prop%defined_by( )
      end associate
    end do
    call output_table( header, table, file_unit )

  end subroutine output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_property_set
