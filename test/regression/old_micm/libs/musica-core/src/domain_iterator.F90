! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_domain_iterator module

!> The abstract domain_iterator_t type and related functions
module musica_domain_iterator

  use musica_constants,                only : musica_dk, musica_ik
  use musica_iterator,                 only : iterator_t

  implicit none
  private

  public :: domain_iterator_t, domain_iterator_ptr

  !> Domain iterator
  !!
  !! Domain iterators are used to access or modify state variables for
  !! sub-components of a model domain. They can be used in combination to
  !! target nested domain components (e.g., cells within a column).
  !!
  type, abstract, extends(iterator_t) :: domain_iterator_t
  end type domain_iterator_t

  !> Unique pointer to iterator_t objects
  type domain_iterator_ptr
    class(domain_iterator_t), pointer :: val_ => null( )
  contains
    final :: finalize
  end type domain_iterator_ptr

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize a unique iterator pointer
  elemental subroutine finalize( this )

    !> Domain pointer
    type(domain_iterator_ptr), intent(inout) :: this

    if( associated( this%val_ ) ) then
      deallocate( this%val_ )
      this%val_ => null( )
    end if

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_domain_iterator
