program rate_constants
  use micm_rate_constant_arrhenius, only : rate_constant_arrhenius_t
  implicit none

  type(rate_constant_arrhenius_t) :: rate

  rate = rate_constant_arrhenius_t(                                           &
    A=real(1), B=real(1), C=real(1), D=real(1), E=real(1))

end program rate_constants