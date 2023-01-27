module rate_constants
  use iso_c_binding

  implicit none

contains

  real(kind=c_double) function arrhenius_rate(temperature, pressure, a, b, c, d, e) bind(c)
    use iso_c_binding,                only : c_double
    use micm_environment,             only : environment_t
    use musica_constants,             only : musica_dk
    use micm_rate_constant_arrhenius, only : rate_constant_arrhenius_t

    type(rate_constant_arrhenius_t) :: rate
    type(environment_t)             :: env
    real(kind=c_double), value      :: temperature, pressure, a, b, c, d, e

    env%temperature = temperature
    env%pressure = pressure

    rate = rate_constant_arrhenius_t( A=a, B=b, C=c, D=d, E=e )

    arrhenius_rate = rate%calculate(env)

  end function arrhenius_rate
end module rate_constants