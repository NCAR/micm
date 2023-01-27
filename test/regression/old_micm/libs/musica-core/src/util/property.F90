! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_property module

!> The property_t type and related functions
module musica_property

  use musica_constants,                only : musica_dk, musica_ik,           &
                                              musica_lk, musica_rk
  use musica_data_type,                only : data_type_t, kTypeUnknown
  use musica_string,                   only : string_t
  use musica_target,                   only : target_t

  implicit none
  private

  public :: property_t, property_ptr

  !> A physical or other property
  !!
  !! \todo add detailed description and example for property_t usage
  !!
  type :: property_t
    private
    !> Name for the property
    type(string_t) :: name_
    !> Property units
    type(string_t) :: units_
    !> Property data type
    type(data_type_t) :: data_type_ = kTypeUnknown
    !> Model element(s) to which the property applies
    class(target_t), pointer :: applies_to_ => null( )
    !> Default property value
    class(*), pointer :: default_value_ => null( )
    !> Model component that defined the property
    type(string_t) :: defined_by_
  contains
    !> @name Assignment
    !! @{
    procedure :: assign_property
    generic :: assignment(=) => assign_property
    !> @}
    !> @name Equality
    !! @{
    procedure :: equals_property
    generic :: operator(==) => equals_property
    procedure :: not_equals_property
    generic :: operator(/=) => not_equals_property
    !> @}
    !> Returns the name of the property
    procedure :: name
    !> Returns the prefix of the property name, if present
    procedure :: prefix
    !> Returns the base name of the property (excluding prefixes)
    procedure :: base_name
    !> Returns the units of the property
    procedure :: units
    !> Returns the property data type
    procedure :: data_type
    !> Returns the model element(s) to which the property applies
    procedure :: applies_to
    !> Sets a given variable to the default value of the property
    procedure :: get_default
    !> Returns the name of the model component that defined the property
    procedure :: defined_by
    !> Fails with an error if two properties are not the same
    procedure :: must_equal
    !> Finalizes the object
    final :: finalize
  end type property_t

  !> Property constructor
  interface property_t
    module procedure :: constructor_config
    module procedure :: constructor
    module procedure :: constructor_property
  end interface

  !> Property pointer
  type :: property_ptr
    class(property_t), pointer :: val_ => null( )
  contains
    !> Finalize the pointer
    final :: property_ptr_finalize
  end type property_ptr

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates a new property from configuration data
  function constructor_config( config, defined_by, applies_to )               &
      result( new_obj )

    use musica_assert,                 only : die
    use musica_config,                 only : config_t

    !> New property_t object
    type(property_t), pointer :: new_obj
    !> Property configuration data
    class(config_t), intent(inout) :: config
    !> Name of calling function
    character(len=*), intent(in) :: defined_by
    !> Model element(s) to which the property applies
    class(target_t), intent(in), optional :: applies_to

    character(len=*), parameter :: my_name = "Property constructor"
    type(string_t) :: temp_str
    logical :: found

    allocate( new_obj )
    call config%get( "name", new_obj%name_, my_name )
    call config%get( "units", new_obj%units_, my_name, found = found )
    call config%get( "data type", temp_str, my_name, found = found )
    if( found ) new_obj%data_type_ = data_type_t( temp_str )
    if( present( applies_to ) ) then
      allocate( new_obj%applies_to_, source = applies_to )
    end if
    new_obj%defined_by_ = defined_by

  end function constructor_config

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates a new property
  function constructor( defined_by, name, units, applies_to, data_type,       &
      default_value ) result( new_obj )

    !> New property_t object
    class(property_t), pointer :: new_obj
    !> Name of calling function
    character(len=*), intent(in) :: defined_by
    !> Property name
    character(len=*), intent(in), optional :: name
    !> Property units
    character(len=*), intent(in), optional :: units
    !> Model element(s) to which the property applies
    class(target_t), intent(in), optional :: applies_to
    !> Property data type
    class(data_type_t), intent(in), optional :: data_type
    !> Default property value
    class(*), intent(in), optional :: default_value

    type(property_t) :: empty_prop

    !> \bug Possible bug in gfortran not picking up optical data_type?
    if( present( data_type ) ) then
      new_obj => constructor_property( empty_prop, defined_by, name = name,   &
          units = units, applies_to = applies_to, data_type = data_type,      &
          default_value = default_value )
    else
      new_obj => constructor_property( empty_prop, defined_by, name = name,   &
          units = units, applies_to = applies_to,                             &
          default_value = default_value )
    end if

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates a new property from an existing property
  function constructor_property( other, defined_by, name, units, applies_to,  &
      data_type, default_value ) result( new_obj )

    use musica_assert,                 only : die

    !> New property_t object
    class(property_t), pointer :: new_obj
    !> Property to copy
    class(property_t), intent(in) :: other
    !> Name of calling function
    character(len=*), intent(in) :: defined_by
    !> Property name
    character(len=*), intent(in), optional :: name
    !> Property units
    character(len=*), intent(in), optional :: units
    !> Model element(s) to which the property applies
    class(target_t), intent(in), optional :: applies_to
    !> Property data type
    type(data_type_t), intent(in), optional :: data_type
    !> Default property value
    class(*), intent(in), optional :: default_value

    class(property_t), pointer :: temp_obj

    allocate( temp_obj, mold = other )
    call assign_property( temp_obj, other )
    temp_obj%defined_by_ = defined_by
    if( present( name          ) ) temp_obj%name_      = name
    if( present( units         ) ) temp_obj%units_     = units
    if( present( applies_to    ) ) then
      if( associated( temp_obj%applies_to_ ) )                                &
          deallocate( temp_obj%applies_to_ )
      allocate( temp_obj%applies_to_, source = applies_to )
    end if
    if( present( data_type     ) ) temp_obj%data_type_ = data_type
    if( present( default_value ) ) then
      if( associated( temp_obj%default_value_ ) )                             &
          deallocate( temp_obj%default_value_ )
      allocate( temp_obj%default_value_, source = default_value )
    end if
    new_obj => temp_obj

  end function constructor_property

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Assigns property data from an existing property
  subroutine assign_property( to, from )

    !> Property to assign
    class(property_t), intent(inout) :: to
    !> New property data
    type(property_t), intent(in) :: from

    to%name_ = from%name_
    to%units_ = from%units_
    if( associated( to%applies_to_ ) ) deallocate( to%applies_to_ )
    to%applies_to_ => null( )
    if( associated( from%applies_to_ ) ) then
      allocate( to%applies_to_, source = from%applies_to_ )
    end if
    to%data_type_ = from%data_type_
    if( associated( to%default_value_ ) ) deallocate( to%default_value_ )
    to%default_value_ => null( )
    if( associated( from%default_value_ ) ) then
      allocate( to%default_value_, source = from%default_value_ )
    end if
    to%defined_by_ = from%defined_by_

  end subroutine assign_property

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Evalutes two properties for equality
  logical function equals_property( a, b )

    !> First property
    class(property_t), intent(in) :: a
    !> Second property
    class(property_t), intent(in) :: b

    equals_property = a%name_ .eq. b%name_
    equals_property = equals_property .and. a%units_      .eq. b%units_
    equals_property = equals_property .and. a%data_type_  .eq. b%data_type_
    if( associated( a%applies_to_ ) .and. associated( b%applies_to_ ) ) then
      equals_property = equals_property .and. a%applies_to_ .eq. b%applies_to_
    else
      equals_property = equals_property .and.                                 &
        .not. associated( a%applies_to_ ) .and.                               &
        .not. associated( b%applies_to_ )
    end if

  end function equals_property

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Evaluates two properties for inequality
  logical function not_equals_property( a, b )

    !> First property
    class(property_t), intent(in) :: a
    !> Second property
    class(property_t), intent(in) :: b

    not_equals_property = .not. a .eq. b

  end function not_equals_property

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the name of the property
  type(string_t) function name( this )

    !> Property
    class(property_t), intent(in) :: this

    name = this%name_

  end function name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the prefix of the property name, if present
  !!
  !! The name prefix indicates a sub-group to which the property belongs. If
  !! no prefix is present, and empty string is returned.
  type(string_t) function prefix( this )

    !> Property
    class(property_t), intent(in) :: this

    type(string_t), allocatable :: substrs(:)

    substrs = this%name_%split( "%" )
    if( size( substrs ) .gt. 1 ) then
      prefix = substrs( 1 )
    else
      prefix = ""
    end if

  end function prefix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the base name of a property (excluding prefixes)
  type(string_t) function base_name( this )

    !> Property
    class(property_t), intent(in) :: this

    type(string_t), allocatable :: substrs(:)

    substrs = this%name_%split( "%" )
    base_name = substrs( size( substrs ) )

  end function base_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the units of the property
  type(string_t) function units( this )

    !> Property
    class(property_t), intent(in) :: this

    units = this%units_

  end function units

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the property data type
  type(data_type_t) function data_type( this )

    !> Property
    class(property_t), intent(in) :: this

    data_type = this%data_type_

  end function data_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the model element(s) to which the property applies
  function applies_to( this )

    !> Model element(s) to which the property applies
    class(target_t), pointer :: applies_to
    !> Property
    class(property_t), intent(in) :: this

    if( associated( this%applies_to_ ) ) then
      allocate( applies_to, source = this%applies_to_ )
    else
      nullify( applies_to )
    end if

  end function applies_to

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets a given variable to the default value of the property
  !!
  !! If no default value is specified, the variable is left unchanged
  !!
  subroutine get_default( this, variable )

    use musica_assert,                 only : die, die_msg

    !> Property
    class(property_t), intent(in) :: this
    !> Variable to set to the default value
    class(*), intent(inout) :: variable

    if( .not. associated( this%default_value_ ) ) return
    select type( default_value => this%default_value_ )
    type is( integer(kind=musica_ik) )
      select type( variable )
      type is( integer(kind=musica_ik) )
        variable = default_value
      class default
        call die_msg( 297185249, "Type mismatch setting default value for "// &
                      "property '"//this%name_//"'" )
      end select
    type is( real(kind=musica_rk) )
      select type( variable )
      type is( real(kind=musica_rk) )
        variable = default_value
      class default
        call die_msg( 350211796, "Type mismatch setting default value for "// &
                      "property '"//this%name_//"'" )
      end select
    type is( real(kind=musica_dk) )
      select type( variable )
      type is( real(kind=musica_dk) )
        variable = default_value
      class default
        call die_msg( 122216333, "Type mismatch setting default value for "// &
                      "property '"//this%name_//"'" )
      end select
    type is( logical(kind=musica_lk) )
      select type( variable )
      type is( logical(kind=musica_lk) )
        variable = default_value
      class default
        call die_msg( 176696119, "Type mismatch setting default value for "// &
                      "property '"//this%name_//"'" )
      end select
    class default
      call die( 971095319 )
    end select

  end subroutine get_default

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the name of the model component that defined the property
  type(string_t) function defined_by( this )

    !> Property
    class(property_t), intent(in) :: this

    defined_by = this%defined_by_

  end function defined_by

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Fails with an error if two properties are not the same
  subroutine must_equal( this, other )

    use musica_assert,                 only : assert_msg
    use musica_string

    !> Property
    class(property_t), intent(in) :: this
    !> Propery to compare with
    class(property_t), intent(in) :: other

    call assert_msg( 916556332, this .eq. other, "Property mismatch for '"//  &
                     this%name_//"' as defined by '"//this%defined_by_//      &
                     "' and '"//other%defined_by_//"'" )

  end subroutine must_equal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize a property_t object
  elemental subroutine finalize( this )

    !> Property
    type(property_t), intent(inout) :: this

    if( associated( this%applies_to_ ) ) then
      deallocate( this%applies_to_ )
      this%applies_to_ => null( )
    end if
    if( associated( this%default_value_ ) ) then
      deallocate( this%default_value_ )
      this%default_value_ => null( )
    end if

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize a property_ptr object
  elemental subroutine property_ptr_finalize( this )

    !> Property pointer
    type(property_ptr), intent(inout) :: this

    if( associated( this%val_ ) ) then
      deallocate( this%val_ )
      this%val_ => null( )
    end if

  end subroutine property_ptr_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_property
