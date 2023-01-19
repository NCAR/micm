program call_cpp_test
  use iso_c_binding

  implicit none

  interface
    type(c_funptr) function get_solver(filepath) bind(c)
      import :: c_char, c_funptr
      character(len=1, kind=c_char), dimension(*), intent(in) :: filepath
    end function get_solver
  end interface

  interface
    function solver(arg1, arg2, arg3)
      import :: c_ptr, c_double
      real(c_double), dimension(*) :: arg1
      real(c_double), dimension(*) :: arg2
      real(c_double), dimension(*) :: arg3
    end function
  end interface

  call run()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine run()
    type(c_funptr) :: csolver_func_pointer
    procedure(solver), pointer :: fsolver
    character(len=*, kind=c_char), parameter :: filepath = "hello.txt"

    real(c_double), dimension(:), allocatable :: arg1
    real(c_double), dimension(:), allocatable :: arg2
    real(c_double), dimension(:), allocatable :: arg3

    allocate(arg1(5))
    allocate(arg2(5))
    allocate(arg3(5))

    csolver_func_pointer = get_solver(filepath)
    call c_f_procpointer(csolver_func_pointer, fsolver)

    call fsolver(arg1, arg2, arg3)
  end subroutine run

end program call_cpp_test
