module call_cpp_test
  use iso_c_binding

  implicit none

  interface
    type(c_funptr) function get_solver(filepath) bind(c)
      import :: c_char
      character(len=1, kind=c_char), dimension(*), intent(in) :: filepath
    end function get_solver
  end interface

  ! interface
  !   function solver(a)
  !     import :: c_ptr
  !     real(c_float), intent(in) :: a
  !     real(c_float) :: func
  !   end function
  ! end interface

  run()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine run()
    type(c_funptr) :: solver
    character(len=*, kind=c_char), parameter :: filepath = "hello.txt"

    print *, 'hello'
    solver = get_solver(filepath)
  end subroutine run

end module call_cpp_test
