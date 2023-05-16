! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_domain_target_columns module

!> The abstract domain_domain_target_columns_t type and related functions
module musica_domain_target_columns

  use musica_constants,                only : musica_dk, musica_ik
  use musica_target,                   only : target_t

  implicit none
  private

  public :: domain_target_columns_t

  !> Target for all columns in a domain
  type, extends(target_t) :: domain_target_columns_t
  contains
    procedure :: name => target_name
    procedure :: equals_target
  end type domain_target_columns_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the name of the target
  type(string_t) function target_name( this )

    use musica_string,                 only : string_t

    !> Target
    class(domain_target_columns_t), intent(in) :: this

    target_name = "all domain columns"

  end function target_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Equality comparison
  logical function equals_target( a, b ) result( eq )

    !> Target
    class(domain_target_columns_t), intent(in) :: a
    !> Other target
    class(target_t), intent(in) :: b

    select type( b )
    class is( domain_target_columns_t )
      eq = .true.
    class default
      eq = .false.
    end select

  end function equals_target

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_domain_target_columns
