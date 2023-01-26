! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_domain module

!> The abstract domain_t type and related functions
module musica_domain

  use musica_constants,                only : musica_dk, musica_ik
  use musica_property_set,             only : property_set_t

  implicit none
  private

  public :: domain_t, domain_ptr

  !> A model domain of abstract structure
  !!
  !! Extending classes of domain_t define the structure of the domain and can
  !! be used to build domain state objects and related accessors/mutators.
  !!
  !! The general usage of \c domain_t objects is to:
  !! - create a domain_t object using the
  !!   \c musica_domain_factory::domain_builder function
  !! - register any needed state variables and properies using the \c domain_t
  !!   type-bound \c register \c mutator and \c accessor functions for the
  !!   domain subset you are interested in (e.g., all cells, surface cells,
  !!   columns)
  !! - use the \c domain_t type-bound \c new_state function to get a state
  !!   object to use for the domain
  !! - during solving, use the accessors and mutators registered during
  !!   initialization with the \c domain_state_t::get and
  !!   \c domain_state_t::update functions to access or modify the current
  !!   values of state variables
  !!
  !! Although the structure of the abstract domain types permits run-time
  !! registration of state parameters and variables, it is compatible with
  !! models that use a fixed set of parameters. In this case the domain
  !! registration, accessor and mutator functions would check to make sure
  !! a state variable that is requested is present in the model, and return
  !! an error or warning if they are not found.
  !!
  !! \todo develop a complete set of \c domain_t examples
  !!
  type, abstract :: domain_t
    private
    !> Flag indicating whether the domain configuration has been finalized
    logical :: is_locked_ = .false.
    !> Set of properties that define the domain state
    type(property_set_t) :: properties_
    !> Set of mutators that have been registered with the domain
    type(property_set_t) :: mutators_
    !> Set of accessors that have been registered with the domain
    type(property_set_t) :: accessors_
  contains
    !> Locks the domain configuration
    procedure :: lock
    !> Returns a flag indicating whether the domain configuration has been
    !! locked
    procedure :: is_locked
    !> @name Registers domain state properties
    !! @{
    procedure, private :: register_property
    procedure, private :: register_property_set
    generic :: register => register_property, register_property_set
    !> @}
    !> Returns the set of registered properties for the domain
    procedure :: properties
    !> Requests a mutator for a domain state property
    procedure :: mutator
    !> Requests mutators for a domain state property set
    procedure :: mutator_set
    !> Requests an accessor for a domain state property
    procedure :: accessor
    !> Requests accessors for a domain state property set
    procedure :: accessor_set
    !> Indicates whether a domain state property exists
    procedure :: is_registered
    !> Returns the units of a domain state property
    procedure :: units
    !> Outputs the registered mutators and accessors
    procedure :: output_registry
    !> Returns the domain type as a string
    procedure(domain_type), deferred :: type
    !> Creates a new state for the domain
    procedure(new_state), deferred :: new_state
    !> Returns an iterator for the domain or a supported domain subset
    procedure(iterator), deferred :: iterator
    !> Allocates a mutator for a given data type and target domain
    procedure(allocate_mutator), deferred :: allocate_mutator
    !> Allocates an accessor for a given data type and target domain
    procedure(allocate_accessor), deferred :: allocate_accessor
    !> Private constructor (should only be called by extending types)
    procedure :: private_constructor
  end type domain_t

  !> Unique pointer to domain_t objects
  type domain_ptr
    class(domain_t), pointer :: val_ => null( )
  contains
    final :: finalize
  end type domain_ptr

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the domain type as a string
  function domain_type( this )
    use musica_string,                 only : string_t
    import domain_t
    !> Domain type
    type(string_t) :: domain_type
    !> Domain
    class(domain_t), intent(in) :: this
  end function domain_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates a new domain state object
  function new_state( this )
    use musica_domain_state,           only : domain_state_t
    import domain_t
    !> New domain state
    class(domain_state_t), pointer :: new_state
    !> Domain
    class(domain_t), intent(in) :: this
  end function new_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns an iterator for the domain or a supported domain subset
  function iterator( this, target_domain )
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_target,                 only : target_t
    import domain_t
    !> New iterator
    class(domain_iterator_t), pointer :: iterator
    !> Domain
    class(domain_t), intent(in) :: this
    !> Target for the iterator
    class(target_t), intent(in) :: target_domain
  end function iterator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates a mutator for a given target domain and data type
  function allocate_mutator( this, target_domain, data_type, property_index ) &
      result( mutator )
    use musica_constants,              only : musica_ik
    use musica_data_type,              only : data_type_t
    use musica_domain_state_mutator,   only : domain_state_mutator_t
    use musica_target,                 only : target_t
    import domain_t
    !> Mutator
    class(domain_state_mutator_t), pointer :: mutator
    !> Domain
    class(domain_t), intent(in) :: this
    !> Target domain for the mutator
    class(target_t), intent(in) :: target_domain
    !> Data type for the mutatable property
    type(data_type_t), intent(in) :: data_type
    !> Index for the property among registered properties of the same
    !! target domain and data type
    integer(kind=musica_ik), intent(in) :: property_index
  end function allocate_mutator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates a accessor for a given target domain and data type
  function allocate_accessor( this, target_domain, data_type, property_index ) &
      result( accessor )
    use musica_constants,              only : musica_ik
    use musica_data_type,              only : data_type_t
    use musica_domain_state_accessor,  only : domain_state_accessor_t
    use musica_target,                 only : target_t
    import domain_t
    !> Accessor
    class(domain_state_accessor_t), pointer :: accessor
    !> Domain
    class(domain_t), intent(in) :: this
    !> Target domain for the accessor
    class(target_t), intent(in) :: target_domain
    !> Data type for the accessible property
    type(data_type_t), intent(in) :: data_type
    !> Index for the property among registered properties of the same
    !! target domain and data type
    integer(kind=musica_ik), intent(in) :: property_index
  end function allocate_accessor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Locks the domain configuration
  subroutine lock( this )

    !> Domain
    class(domain_t), intent(inout) :: this

    this%is_locked_ = .true.

  end subroutine lock

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns whether the domain configuration has been locked
  logical function is_locked( this )

    !> Domain
    class(domain_t), intent(in) :: this

    is_locked = this%is_locked_

  end function is_locked

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Registers a domain state property
  subroutine register_property( this, property )

    use musica_assert,                 only : assert, assert_msg
    use musica_property,               only : property_t
    use musica_string,                 only : string_t

    !> Domain
    class(domain_t), intent(inout) :: this
    !> Property to add to the domain registry
    class(property_t), intent(in) :: property

    logical :: found
    type(string_t) :: prop_name
    class(property_t), pointer :: existing_prop

    prop_name = property%name( )
    existing_prop => this%properties_%get( prop_name%to_char( ),              &
                                           found = found )
    if( found ) then
      call assert_msg( 460201666, existing_prop .eq. property,                &
                       "Trying to overwrite domain property '"//              &
                       trim( prop_name%to_char( ) )//"'" )
      deallocate( existing_prop )
      return
    end if
    call assert( 823250106, .not. this%is_locked_ )
    call this%properties_%add( property )

  end subroutine register_property

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Registers a domain property set
  subroutine register_property_set( this, prefix, property_set )

    use musica_assert,                 only : assert
    use musica_property,               only : property_t, property_ptr
    use musica_property_set,           only : property_set_t
    use musica_string,                 only : string_t

    !> Domain
    class(domain_t), intent(inout) :: this
    !> Prefix to attach to the property names
    character(len=*), intent(in) :: prefix
    !> Property set to add to the domain registry
    class(property_set_t), intent(in) :: property_set

    integer(kind=musica_ik) :: i_prop
    class(property_t), pointer :: prop, new_prop
    type(string_t) :: new_name, defined_by

    call assert( 995312544, .not. this%is_locked_ )
    do i_prop = 1, property_set%size( )
      prop => property_set%get( i_prop )
      new_name = prop%name( )
      if( len( prefix ) .gt. 0 ) new_name = prefix//"%"//new_name
      defined_by = prop%defined_by( )
      new_prop => property_t( prop, defined_by%to_char( ),                    &
                              name = new_name%to_char( ) )
      call this%properties_%add( new_prop )
      deallocate( new_prop )
      deallocate( prop     )
    end do

  end subroutine register_property_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the set of registered properties for the domain
  function properties( this )

    use musica_assert,                 only : assert

    !> Registered properties
    type(property_set_t) :: properties
    !> Domain
    class(domain_t), intent(in) :: this

    call assert( 203819822, this%is_locked_ )
    properties = this%properties_%subset( )

  end function properties

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Requests a mutator for a domain state property
  function mutator( this, property )

    use musica_assert,                 only : assert_msg
    use musica_domain_state_mutator,   only : domain_state_mutator_t
    use musica_property,               only : property_t
    use musica_string,                 only : string_t
    use musica_target,                 only : target_t

    !> Mutator
    class(domain_state_mutator_t), pointer :: mutator
    !> Domain
    class(domain_t), intent(inout) :: this
    !> Property to get mutator for
    class(property_t), intent(in) :: property

    type(string_t) :: prop_name, owner
    class(property_t), pointer :: prop_reg, mut_prop
    type(property_set_t) :: type_props
    class(target_t), pointer :: mutator_target
    integer(kind=musica_ik) :: prop_index

    prop_name = property%name( )
    owner     = property%defined_by( )
    prop_reg => this%properties_%get( prop_name%to_char( ) )
    call assert_msg( 422651221, property .eq. prop_reg,                       &
                     "Property mismatch registering mutator for domain "//    &
                     "property '"//prop_name%to_char( )//"'" )
    mut_prop => property_t( prop_reg, owner%to_char( ) )
    deallocate( prop_reg )
    call this%mutators_%add( mut_prop )
    mutator_target => mut_prop%applies_to( )
    type_props = this%properties_%subset( data_type = mut_prop%data_type( ),  &
                                           applies_to = mutator_target )
    prop_index = type_props%index( mut_prop )
    mutator => this%allocate_mutator( mutator_target,                         &
                                      mut_prop%data_type( ),                  &
                                      prop_index )
    call mutator%attach_property( mut_prop )
    deallocate( mutator_target )
    deallocate( mut_prop )

  end function mutator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Requests mutators for a set of domain state properties
  !!
  !! All members of the property set must have the same units and target
  !! (sub)domain. Otherwise, the mutators must be requested individually.
  !!
  function mutator_set( this, variable_name, units, data_type, applies_to,    &
      requestor )

    use musica_assert,                 only : assert
    use musica_data_type,              only : data_type_t
    use musica_domain_state_mutator,   only : domain_state_mutator_ptr
    use musica_property,               only : property_t
    use musica_string,                 only : string_t
    use musica_target,                 only : target_t

    !> Mutators for the requested state variable set
    class(domain_state_mutator_ptr), pointer :: mutator_set(:)
    !> Domain
    class(domain_t), intent(inout) :: this
    !> Name of the property set
    character(len=*), intent(in) :: variable_name
    !> Units for the property set members
    character(len=*), intent(in) :: units
    !> Data type for the property set members
    class(data_type_t), intent(in) :: data_type
    !> Model element(s) to which the properties apply
    class(target_t), intent(in) :: applies_to
    !> Name of the model component requesting the mutator
    character(len=*), intent(in) :: requestor

    class(property_t), pointer :: prop_base, prop_spec, prop_reg
    type(property_set_t) :: prop_subset
    integer(kind=musica_ik) :: num_props, i_mutator
    type(string_t) :: prop_name

    call assert( 814373165, len( trim( variable_name ) ) .gt. 0 )

    prop_subset = this%properties_%subset( prefix = variable_name )
    num_props = prop_subset%size( )
    prop_base => property_t( requestor,                                       &
                             units = units,                                   &
                             data_type = data_type,                           &
                             applies_to = applies_to )
    allocate( mutator_set( num_props ) )
    do i_mutator = 1, num_props
      prop_reg => prop_subset%get( i_mutator )
      prop_name = prop_reg%name( )
      prop_spec => property_t( prop_base, requestor,                          &
                               name = prop_name%to_char( ) )
      mutator_set( i_mutator )%val_ => this%mutator( prop_spec )
      deallocate( prop_spec )
      deallocate( prop_reg  )
    end do
    deallocate( prop_base   )

  end function mutator_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Requests a accessor for a domain state property
  function accessor( this, property )

    use musica_assert,                 only : assert_msg
    use musica_domain_state_accessor,  only : domain_state_accessor_t
    use musica_property,               only : property_t
    use musica_string,                 only : string_t
    use musica_target,                 only : target_t

    !> Accessor
    class(domain_state_accessor_t), pointer :: accessor
    !> Domain
    class(domain_t), intent(inout) :: this
    !> Property to get accessor for
    class(property_t), intent(in) :: property

    type(string_t) :: prop_name, owner
    class(property_t), pointer :: prop_reg, acc_prop
    type(property_set_t) :: type_props
    class(target_t), pointer :: accessor_target
    integer(kind=musica_ik) :: prop_index

    prop_name = property%name( )
    owner     = property%defined_by( )
    prop_reg => this%properties_%get( prop_name%to_char( ) )
    call assert_msg( 255426392, property .eq. prop_reg,                       &
                     "Property mismatch registering accessor for domain "//   &
                     "property '"//prop_name%to_char( )//"'" )
    acc_prop => property_t( prop_reg, owner%to_char( ) )
    deallocate( prop_reg )
    call this%accessors_%add( acc_prop )
    accessor_target => acc_prop%applies_to( )
    type_props = this%properties_%subset( data_type = acc_prop%data_type( ),  &
                                           applies_to = accessor_target )
    prop_index = type_props%index( acc_prop )
    accessor => this%allocate_accessor( accessor_target,                      &
                                        acc_prop%data_type( ),                &
                                        prop_index )
    call accessor%attach_property( acc_prop )
    deallocate( accessor_target )
    deallocate( acc_prop )

  end function accessor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Requests accessors for a set of domain state properties
  !!
  !! All members of the property set must have the same units and target
  !! (sub)domain. Otherwise, the accessors must be requested individually.
  !!
  function accessor_set( this, variable_name, units, data_type, applies_to,   &
      requestor )

    use musica_assert,                 only : assert
    use musica_data_type,              only : data_type_t
    use musica_domain_state_accessor,  only : domain_state_accessor_ptr
    use musica_property,               only : property_t
    use musica_string,                 only : string_t
    use musica_target,                 only : target_t

    !> Accessors for the requested state variable set
    class(domain_state_accessor_ptr), pointer :: accessor_set(:)
    !> Domain
    class(domain_t), intent(inout) :: this
    !> Name of the property set
    character(len=*), intent(in) :: variable_name
    !> Units for the property set members
    character(len=*), intent(in) :: units
    !> Data type for the property set members
    class(data_type_t), intent(in) :: data_type
    !> Model element(s) to which the properties apply
    class(target_t), intent(in) :: applies_to
    !> Name of the model component requesting the accessor
    character(len=*), intent(in) :: requestor

    class(property_t), pointer :: prop_base, prop_spec, prop_reg
    type(property_set_t) :: prop_subset
    integer(kind=musica_ik) :: num_props, i_accessor
    type(string_t) :: prop_name

    call assert( 309906178, len( trim( variable_name ) ) .gt. 0 )

    prop_subset = this%properties_%subset( prefix = variable_name )
    num_props = prop_subset%size( )
    prop_base => property_t( requestor,                                       &
                             units = units,                                   &
                             data_type = data_type,                           &
                             applies_to = applies_to )
    allocate( accessor_set( num_props ) )
    do i_accessor = 1, num_props
      prop_reg => prop_subset%get( i_accessor )
      prop_name = prop_reg%name( )
      prop_spec => property_t( prop_base, requestor,                          &
                               name = prop_name%to_char( ) )
      accessor_set( i_accessor )%val_ => this%accessor( prop_spec )
      deallocate( prop_spec )
      deallocate( prop_reg  )
    end do
    deallocate( prop_base   )

  end function accessor_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Indicates whether a domain state property has been registered
  logical function is_registered( this, property_name )

    use musica_property,               only : property_t

    !> Domain
    class(domain_t), intent(in) :: this
    !> Name of the property to look for
    character(len=*), intent(in) :: property_name

    class(property_t), pointer :: prop
    integer(kind=musica_ik) :: var_id

    prop => this%properties_%get( property_name, found = is_registered )
    if( associated( prop ) ) deallocate( prop )

  end function is_registered

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the units of a domain state property
  function units( this, property_name )

    use musica_assert,                 only : assert
    use musica_property,               only : property_t
    use musica_string,                 only : string_t

    !> Units for the property
    type(string_t) :: units
    !> Domain
    class(domain_t), intent(in) :: this
    !> Name of the registered state property
    character(len=*), intent(in) :: property_name

    class(property_t), pointer :: prop

    call assert( 425563855, len( trim( property_name ) ) .gt. 0 )
    prop => this%properties_%get( property_name )
    units = prop%units( )
    deallocate( prop )

  end function units

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Outputs the registered mutators and accessors
  subroutine output_registry( this, file_unit )

    use musica_string,                 only : string_t

    !> Domain
    class(domain_t), intent(in) :: this
    !> File unit to output to
    integer, intent(in), optional :: file_unit

    integer(kind=musica_ik) :: f

    f = 6
    if( present( file_unit ) ) f = file_unit
    write(f,*)
    write(f,*) "Registered domain properties"
    call this%properties_%output( f )
    write(f,*)
    write(f,*) "Registered mutators"
    call this%mutators_%output( f )
    write(f,*)
    write(f,*) "Registered accessors"
    call this%accessors_%output( f )

  end subroutine output_registry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Private constructor (should only be called by extending types)
  subroutine private_constructor( this )

    !> Domain
    class(domain_t), intent(inout) :: this

    this%properties_ = property_set_t( )
    this%mutators_   = property_set_t( )
    this%accessors_  = property_set_t( )

  end subroutine private_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalizes a unique domain pointer
  elemental subroutine finalize( this )

    !> Domain pointer
    type(domain_ptr), intent(inout) :: this

    if( associated( this%val_ ) ) then
      deallocate( this%val_ )
      this%val_ => null( )
    end if

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_domain
