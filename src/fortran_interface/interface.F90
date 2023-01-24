module micm_solver_interface
  use iso_c_binding

  implicit none

  interface
    type(c_funptr) function get_solver(filepath) bind(c)
      import :: c_char, c_funptr
      character(len=1, kind=c_char), dimension(*), intent(in) :: filepath
    end function get_solver
  end interface

  interface
    subroutine solver(arg1, arg2, arg3) bind(c)
      import :: c_ptr, c_double
      real(c_double), dimension(*) :: arg1
      real(c_double), dimension(*) :: arg2
      real(c_double), dimension(*) :: arg3
    end subroutine
  end interface

end module micm_solver_interface
