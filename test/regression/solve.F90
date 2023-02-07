module solve_mod
  use iso_c_binding

  implicit none

contains

  subroutine solve() bind(c)
    use iso_c_binding,              only : c_double
    use micm_ODE_solver_rosenbrock, only : ODE_solver_rosenbrock_t
    use musica_config,              only : config_t
    use micm_environment,           only : environment_t
    use micm_kinetics,              only : kinetics_t
    use musica_constants,           only : musica_dk, musica_ik

    type(config_t) :: config
    type(environment_t) :: env
    ! type(kinetics_t) :: kinetics
    class(kinetics_t), pointer :: kinetics
    class(ODE_solver_rosenbrock_t), pointer :: solver 
    real(kind=musica_dk), allocatable :: number_densities__molec_cm3_(:)
    integer(kind=musica_ik) :: error_flag

    call config%from_file( "configs/solver.json" )


    env%temperature = 273.15
    env%pressure = 100000
    env%number_density_air = 2.7e19
    env%photolysis_rate_constants = (/1e-4, 1e-5, 1e-6/)

    kinetics => kinetics_t()
    call kinetics%update(env)

    solver => ODE_solver_rosenbrock_t( config )

    allocate( number_densities__molec_cm3_( solver%N ) )

    number_densities__molec_cm3_ = 2e15

    call solver%solve( TStart = 0.0_musica_dk,                                &
      TEnd = 1.0_musica_dk,                                                   &
      y = number_densities__molec_cm3_,                                       &
      theKinetics = kinetics,                                                 &
      IErr = error_flag                                                       &
    )

  end subroutine solve

end module solve_mod