module tuv_error_mod
  implicit none

contains

  subroutine tuv_error_fatal(msg)

    character(len=*), intent(in) :: msg

    write(*,*) 'ERROR: '//trim(msg)
    call exit(-1)
    
  end subroutine tuv_error_fatal
end module tuv_error_mod
