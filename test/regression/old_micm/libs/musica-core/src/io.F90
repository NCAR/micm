! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The musica_io module

!> The io_t type and related functions
module musica_io

  implicit none
  private

  public :: io_t

  !> General input/output class
  type, abstract :: io_t
  contains
    !> @name Data read functions
    !! @{
    procedure(read_0D_double), deferred :: read_0D_double
    procedure(read_1D_double), deferred :: read_1D_double
    procedure(read_2D_double), deferred :: read_2D_double
    procedure(read_3D_double), deferred :: read_3D_double
    procedure(read_4D_double), deferred :: read_4D_double
    procedure(read_0D_int),    deferred :: read_0D_int
    procedure(read_1D_int),    deferred :: read_1D_int
    generic :: read => read_0D_double, read_1D_double, read_2D_double,        &
                       read_3D_double, read_4D_double, read_0D_int,           &
                       read_1D_int
    !> @}
    !> Returns the dimension names for a given variable
    procedure(variable_dimensions), deferred :: variable_dimensions
    !> Returns the units for a given variable
    procedure(variable_units),      deferred :: variable_units
  end type io_t

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 0D double-precision floating-point data
  subroutine read_0D_double( this, variable_name, container, requestor_name )
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    import io_t
    class(io_t),          intent(inout) :: this
    class(string_t),      intent(in)    :: variable_name
    real(kind=musica_dk), intent(out)   :: container
    character(len=*),     intent(in)    :: requestor_name
  end subroutine read_0D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 1D double-precision floating-point data
  !!
  !! If \c container is unallocated, it will be allocated to the dimensions
  !! of the read variable. Otherwise, its dimensions must match those of the
  !! read variable.
  !!
  subroutine read_1D_double( this, variable_name, container, requestor_name )
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    import io_t
    class(io_t),                       intent(inout) :: this
    class(string_t),                   intent(in)    :: variable_name
    real(kind=musica_dk), allocatable, intent(out)   :: container(:)
    character(len=*),                  intent(in)    :: requestor_name
  end subroutine read_1D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 2D double-precision floating-point data
  !!
  !! If \c container is unallocated, it will be allocated to the dimensions
  !! of the read variable. Otherwise, its dimensions must match those of the
  !! read variable.
  !!
  subroutine read_2D_double( this, variable_name, container, requestor_name )
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    import io_t
    class(io_t),                       intent(inout) :: this
    class(string_t),                   intent(in)    :: variable_name
    real(kind=musica_dk), allocatable, intent(out)   :: container(:,:)
    character(len=*),                  intent(in)    :: requestor_name
  end subroutine read_2D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 3D double-precision floating-point data
  !!
  !! If \c container is unallocated, it will be allocated to the dimensions
  !! of the read variable. Otherwise, its dimensions must match those of the
  !! read variable.
  !!
  subroutine read_3D_double( this, variable_name, container, requestor_name )
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    import io_t
    class(io_t),                       intent(inout) :: this
    class(string_t),                   intent(in)    :: variable_name
    real(kind=musica_dk), allocatable, intent(out)   :: container(:,:,:)
    character(len=*),                  intent(in)    :: requestor_name
  end subroutine read_3D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 4D double-precision floating-point data
  !!
  !! If \c container is unallocated, it will be allocated to the dimensions
  !! of the read variable. Otherwise, its dimensions must match those of the
  !! read variable.
  !!
  subroutine read_4D_double( this, variable_name, container, requestor_name )
    use musica_constants,              only : musica_dk
    use musica_string,                 only : string_t
    import io_t
    class(io_t),                       intent(inout) :: this
    class(string_t),                   intent(in)    :: variable_name
    real(kind=musica_dk), allocatable, intent(out)   :: container(:,:,:,:)
    character(len=*),                  intent(in)    :: requestor_name
  end subroutine read_4D_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 0D integer data
  subroutine read_0D_int( this, variable_name, container, requestor_name )
    use musica_string,                 only : string_t
    import io_t
    class(io_t),      intent(inout) :: this
    class(string_t),  intent(in)    :: variable_name
    integer,          intent(out)   :: container
    character(len=*), intent(in)    :: requestor_name
  end subroutine read_0D_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads 1D integer data
  !!
  !! If \c container is unallocated, it will be allocated to the dimensions
  !! of the read variable. Otherwise, its dimensions must match those of the
  !! read variable.
  !!
  subroutine read_1D_int( this, variable_name, container, requestor_name )
    use musica_string,                 only : string_t
    import io_t
    class(io_t),                   intent(inout) :: this
    class(string_t),               intent(in)    :: variable_name
    integer,          allocatable, intent(out)   :: container(:)
    character(len=*),              intent(in)    :: requestor_name
  end subroutine read_1D_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the dimension names for a given variable
  function variable_dimensions( this, variable_name, requestor_name )         &
      result( dimensions )
    use musica_string,                 only : string_t
    import io_t
    type(string_t), allocatable  :: dimensions(:)
    class(io_t),      intent(in) :: this
    class(string_t),  intent(in) :: variable_name
    character(len=*), intent(in) :: requestor_name
  end function variable_dimensions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the units for a given variable
  function variable_units( this, variable_name, requestor_name )
    use musica_string,                 only : string_t
    import io_t
    type(string_t)               :: variable_units
    class(io_t),      intent(in) :: this
    class(string_t),  intent(in) :: variable_name
    character(len=*), intent(in) :: requestor_name
  end function variable_units

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

end module musica_io
