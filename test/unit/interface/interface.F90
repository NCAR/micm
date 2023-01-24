program call_cpp_test
  use micm_solver_interface, only: get_solver, solver
  use iso_c_binding

  implicit none

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
