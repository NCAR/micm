module solve_mod
  use iso_c_binding

  implicit none

contains

  subroutine solve(temperature, pressure, number_density_air, time_start,     &
    time_end, pr1, pr2, pr3, number_densities, new_number_densities) bind(c)

    use iso_c_binding,              only : c_double
    use micm_ODE_solver_rosenbrock, only : ODE_solver_rosenbrock_t
    use musica_config,              only : config_t
    use micm_environment,           only : environment_t
    use micm_kinetics,              only : kinetics_t
    use musica_constants,           only : musica_dk, musica_ik

    real(kind=c_double), value :: temperature, pressure, number_density_air, time_start, time_end, pr1, pr2, pr3
    real(kind=c_double), intent(in) :: number_densities(:)
    real(kind=c_double), pointer, intent(out) :: new_number_densities(:)

    real(kind=musica_dk), allocatable :: f_number_densities(:)

    type(config_t) :: config
    type(environment_t) :: env
    class(kinetics_t), pointer :: kinetics
    class(ODE_solver_rosenbrock_t), pointer :: solver 
    integer(kind=musica_ik) :: error_flag

    call config%from_file("regression_configs/solver.json" )

    env%temperature = temperature
    env%pressure = pressure
    env%number_density_air = number_density_air
    env%photolysis_rate_constants = (/pr1, pr2, pr3/)

    kinetics => kinetics_t()
    call kinetics%update(env)

    solver => ODE_solver_rosenbrock_t( config )

    f_number_densities = number_densities

    call solver%solve( TStart = time_start,                                   &
      TEnd = time_end,                                                        &
      y = f_number_densities,                                                 &
      theKinetics = kinetics,                                                 &
      IErr = error_flag                                                       &
    )

    allocate(new_number_densities(size(f_number_densities)))
    new_number_densities(:) = f_number_densities(:)

  end subroutine solve

end module solve_mod