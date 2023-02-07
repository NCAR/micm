module solve_mod
  use iso_c_binding

  implicit none

contains

  subroutine solve() bind(c)
    use iso_c_binding,              only : c_double
    use micm_ODE_solver_rosenbrock, only : ODE_solver_rosenbrock_t
    use musica_config,              only : config_t

    type(config_t) :: config
    class(ODE_solver_rosenbrock_t), pointer :: solver

    call config%from_file( "configs/solver.json" )

    solver => ODE_solver_rosenbrock_t( config )

    print *, solver%N

  end subroutine solve

end module solve_mod