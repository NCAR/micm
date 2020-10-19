! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> ODE solver factory

!> Builder of ODE solvers
module micm_ODE_solver_factory

  use micm_ODE_solver,                 only : ODE_solver_t
  use micm_ODE_solver_mozart,          only : ODE_solver_mozart_t
  use micm_ODE_solver_rosenbrock,      only : ODE_solver_rosenbrock_t

  implicit none
  private

  public :: ODE_solver_builder

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Build an ODE solver by name
  function ODE_solver_builder( config ) result( new_solver )

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    !> New domain
    class(ODE_solver_t), pointer :: new_solver
    !> Solver configuration data
    type(config_t), intent(inout) :: config

    type(string_t) :: solver_type
    character(len=*), parameter :: my_name = 'ODE solver builder'

    new_solver => null( )
    call config%get( 'type', solver_type, my_name )
    solver_type = solver_type%to_lower( )

    if( solver_type .eq. 'mozart' ) then
      new_solver => ODE_solver_mozart_t( config )
    else if( solver_type .eq. 'rosenbrock' ) then
      new_solver => ODE_solver_rosenbrock_t( config )
    else
      call die_msg( 719420755, "Invalid ODE solver type: "//                  &
                               solver_type%to_char( ) )
    end if

  end function ODE_solver_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_ODE_solver_factory
