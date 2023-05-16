! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_domain_state module

!> The abstract domain_state_t type and related functions
module musica_domain_state

  implicit none
  private

  public :: domain_state_t, domain_state_ptr

  !> Abstract domain state
  !!
  !! A domain state maintains all time-evolving conditions for a domain_t
  !! model domain.
  !!
  type, abstract :: domain_state_t
  contains
    !> @name Gets the value of a state property
    !! @{
    procedure(get_boolean), deferred :: get_boolean
    procedure(get_double),  deferred :: get_double
    procedure(get_float),   deferred :: get_float
    procedure(get_integer), deferred :: get_integer
    generic :: get => get_boolean, get_double, get_float, get_integer
    !> @}
    !> @name Updates the value of a state property
    !! @{
    procedure(update_boolean), deferred :: update_boolean
    procedure(update_double),  deferred :: update_double
    procedure(update_float),   deferred :: update_float
    procedure(update_integer), deferred :: update_integer
    generic :: update => update_boolean, update_double, update_float,         &
                         update_integer
    !> @}
  end type domain_state_t

  !> Unique pointer to domain_state_t objects
  type domain_state_ptr
    class(domain_state_t), pointer :: val_ => null( )
  contains
    final :: finalize
  end type domain_state_ptr

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gets the value of a registered state boolean property
  !!
  !! The value returned will be in the units specified when the accessor was
  !! created.
  subroutine get_boolean( this, iterator, accessor, state_value )
    use musica_constants,              only : musica_lk
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_domain_state_accessor,  only : domain_state_accessor_t
    import domain_state_t
    !> Domain state
    class(domain_state_t), intent(in) :: this
    !> Domain iterator
    class(domain_iterator_t), intent(in) :: iterator
    !> Accessor for the registered state property
    class(domain_state_accessor_t), intent(in) :: accessor
    !> Value of the property
    logical(kind=musica_lk), intent(out) :: state_value
  end subroutine get_boolean

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gets the value of a registered state double property
  !!
  !! The value returned will be in the units specified when the accessor was
  !! created.
  subroutine get_double( this, iterator, accessor, state_value )
    use musica_constants,              only : musica_dk
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_domain_state_accessor,  only : domain_state_accessor_t
    import domain_state_t
    !> Domain state
    class(domain_state_t), intent(in) :: this
    !> Domain iterator
    class(domain_iterator_t), intent(in) :: iterator
    !> Accessor for the registered state property
    class(domain_state_accessor_t), intent(in) :: accessor
    !> Value of the property
    real(kind=musica_dk), intent(out) :: state_value
  end subroutine get_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gets the value of a registered state float property
  !!
  !! The value returned will be in the units specified when the accessor was
  !! created.
  subroutine get_float( this, iterator, accessor, state_value )
    use musica_constants,              only : musica_rk
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_domain_state_accessor,  only : domain_state_accessor_t
    import domain_state_t
    !> Domain state
    class(domain_state_t), intent(in) :: this
    !> Domain iterator
    class(domain_iterator_t), intent(in) :: iterator
    !> Accessor for the registered state property
    class(domain_state_accessor_t), intent(in) :: accessor
    !> Value of the property
    real(kind=musica_rk), intent(out) :: state_value
  end subroutine get_float

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gets the value of a registered state integer property
  !!
  !! The value returned will be in the units specified when the accessor was
  !! created.
  subroutine get_integer( this, iterator, accessor, state_value )
    use musica_constants,              only : musica_ik
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_domain_state_accessor,  only : domain_state_accessor_t
    import domain_state_t
    !> Domain state
    class(domain_state_t), intent(in) :: this
    !> Domain iterator
    class(domain_iterator_t), intent(in) :: iterator
    !> Accessor for the registered state property
    class(domain_state_accessor_t), intent(in) :: accessor
    !> Value of the property
    integer(kind=musica_ik), intent(out) :: state_value
  end subroutine get_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Updates the value of a registered state boolean property
  !!
  !! The units for the value passed to this function must be the same as
  !! those specified when the mutator was created.
  subroutine update_boolean( this, iterator, mutator, state_value )
    use musica_constants,              only : musica_lk
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_domain_state_mutator,   only : domain_state_mutator_t
    import domain_state_t
    !> Domain state
    class(domain_state_t), intent(inout) :: this
    !> Domain iterator
    class(domain_iterator_t), intent(in) :: iterator
    !> Mutator for registered state property
    class(domain_state_mutator_t), intent(in) :: mutator
    !> New value
    logical(kind=musica_lk), intent(in) :: state_value
  end subroutine update_boolean

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Updates the value of a registered state double property
  !!
  !! The units for the value passed to this function must be the same as
  !! those specified when the mutator was created.
  subroutine update_double( this, iterator, mutator, state_value )
    use musica_constants,              only : musica_dk
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_domain_state_mutator,   only : domain_state_mutator_t
    import domain_state_t
    !> Domain state
    class(domain_state_t), intent(inout) :: this
    !> Domain iterator
    class(domain_iterator_t), intent(in) :: iterator
    !> Mutator for registered state property
    class(domain_state_mutator_t), intent(in) :: mutator
    !> New value
    real(kind=musica_dk), intent(in) :: state_value
  end subroutine update_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Updates the value of a registered state float property
  !!
  !! The units for the value passed to this function must be the same as
  !! those specified when the mutator was created.
  subroutine update_float( this, iterator, mutator, state_value )
    use musica_constants,              only : musica_rk
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_domain_state_mutator,   only : domain_state_mutator_t
    import domain_state_t
    !> Domain state
    class(domain_state_t), intent(inout) :: this
    !> Domain iterator
    class(domain_iterator_t), intent(in) :: iterator
    !> Mutator for registered state property
    class(domain_state_mutator_t), intent(in) :: mutator
    !> New value
    real(kind=musica_rk), intent(in) :: state_value
  end subroutine update_float

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Updates the value of a registered state integer property
  !!
  !! The units for the value passed to this function must be the same as
  !! those specified when the mutator was created.
  subroutine update_integer( this, iterator, mutator, state_value )
    use musica_constants,              only : musica_ik
    use musica_domain_iterator,        only : domain_iterator_t
    use musica_domain_state_mutator,   only : domain_state_mutator_t
    import domain_state_t
    !> Domain state
    class(domain_state_t), intent(inout) :: this
    !> Domain iterator
    class(domain_iterator_t), intent(in) :: iterator
    !> Mutator for registered state property
    class(domain_state_mutator_t), intent(in) :: mutator
    !> New value
    integer(kind=musica_ik), intent(in) :: state_value
  end subroutine update_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize a unique state pointer
  elemental subroutine finalize( this )

    !> Domain pointer
    type(domain_state_ptr), intent(inout) :: this

    if( associated( this%val_ ) ) then
      deallocate( this%val_ )
      this%val_ => null ()
    end if

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_domain_state
