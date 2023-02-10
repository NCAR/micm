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

  subroutine reaction_rate_constants(temperature, pressure, pr1, pr2, pr3, rate_constants) bind(c)
    use iso_c_binding,              only : c_double
    use micm_ODE_solver_rosenbrock, only : ODE_solver_rosenbrock_t
    use musica_config,              only : config_t
    use micm_environment,           only : environment_t
    use micm_kinetics,              only : kinetics_t
    use musica_constants,           only : musica_dk, musica_ik

    real(kind=c_double), value :: temperature, pressure, pr1, pr2, pr3
    real(kind=c_double), pointer, intent(out) :: rate_constants(:)
    real(kind=musica_dk), allocatable :: f_rate_constants(:)

    type(environment_t) :: env
    class(kinetics_t), pointer :: kinetics

    env%temperature = temperature
    env%pressure = pressure
    env%photolysis_rate_constants = (/pr1, pr2, pr3/)

    kinetics => kinetics_t()
    call kinetics%update(env)

    f_rate_constants = kinetics%reaction_rate_constants()
    allocate(rate_constants(size(f_rate_constants)))
    rate_constants(:) = f_rate_constants(:)

  end subroutine reaction_rate_constants

end module rate_constants