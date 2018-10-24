module const_props_mod

  implicit none

  private

  type, public :: const_props_type
     private
     character(len=16) :: name = 'UNSET'
     character(len=64) :: desc = 'UNSET'
     real :: molecwght = -huge(1.0)
     integer :: props_set = 0
   contains
     procedure :: set_name => const_props_set_name  
     procedure :: get_name => const_props_get_name
     procedure :: set_desc => const_props_set_desc  
     procedure :: get_desc => const_props_get_desc
     procedure :: set_wght => const_props_set_wght
     procedure :: get_wght => const_props_get_wght
     procedure :: print    => const_props_print
  end type const_props_type

  integer, parameter :: name_pos = 1
  integer, parameter :: desc_pos = 2
  integer, parameter :: wght_pos = 3
  
contains

  subroutine const_props_set_name(this, name)
    class(const_props_type), intent(inout) :: this
    character(len=*), intent(in) :: name

    if (btest(this%props_set, name_pos)) then
       write(*,*) 'ERROR: name already set'
       call abort()
    endif
    this%name = trim(name)
    this%props_set = ibset(this%props_set,name_pos)
  end subroutine const_props_set_name
  
  function const_props_get_name(this) result(x)
    class(const_props_type), intent(in) :: this
    character(len=16) :: x
    x =  this%name
  end function const_props_get_name

  subroutine const_props_set_desc(this, desc)
    class(const_props_type), intent(inout) :: this
    character(len=*), intent(in) :: desc

    if (btest(this%props_set, desc_pos)) then
       write(*,*) 'ERROR: desc already set'
       call abort()
    endif
    this%desc = trim(desc)
    this%props_set = ibset(this%props_set,desc_pos)
  end subroutine const_props_set_desc
  
  function const_props_get_desc(this) result(x)
    class(const_props_type), intent(in) :: this
    character(len=16) :: x
    x =  this%desc
  end function const_props_get_desc

  subroutine const_props_set_wght(this, wgt)
    class(const_props_type), intent(inout) :: this
    real, intent(in) :: wgt

    if (btest(this%props_set,wght_pos)) then
       write(*,*) 'ERROR: molec weight already set'
       call abort()
    endif
    this%molecwght = wgt
    this%props_set = ibset(this%props_set,wght_pos)
  end subroutine const_props_set_wght
  
  function const_props_get_wght(this) result(x)
    class(const_props_type), intent(in) :: this
    real :: x
    x =  this%molecwght
  end function const_props_get_wght

  subroutine const_props_print(this)
    class(const_props_type), intent(in) :: this

    write(*,'(a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    write(*,"(' const name : ',a)")    this%get_name()
    write(*,"(' const desc : ',a)")    this%get_desc()
    write(*,"(' molec wght : ',f8.4)") this%get_wght()
    write(*,'(a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

  end subroutine const_props_print
end module const_props_mod
