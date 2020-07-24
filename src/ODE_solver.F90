! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The micm_ODE_solver module

!> The abstract ODE_solver_t type and related functions
module micm_ODE_solver

  use musica_constants,                only : musica_ik, musica_dk

  implicit none
  private

  public :: ODE_solver_t

  !> A general ODE solver
  !!
  !! Extending classes of ODE_solver_t can be used to solve chemical kinetics
  !! systems during the model run based on kinetics information provided by a
  !! kinetics_t object
  type, abstract :: ODE_solver_t
  private
    !> Integer parameters
    integer(musica_ik), public :: icntrl(20)
    !> Real parameters
    real(musica_dk), public :: rcntrl(20)
  contains
    !> Solve the chemical system
    procedure(solve), deferred  :: solve
  end type ODE_solver_t

interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run the solver with a provided set of kinetics data
  subroutine solve( this, Y, Tstart, Tend, T, theKinetics, Ierr )
    use micm_kinetics,               only : kinetics_t
    use musica_constants,            only : musica_dk
    import ODE_solver_t
    !> ODE solver
    class(ODE_solver_t), intent(inout) :: this
    !> The solver variables
    real(musica_dk), intent(inout) :: Y(:)
    !> Chemistry simulation start time [s]
    real(musica_dk), optional, intent(in) :: Tstart
    !> Chemistry simulation end time [s]
    real(musica_dk), optional, intent(in) :: Tend
    !> Current chemistry simulation time [s]
    real(musica_dk), optional, intent(out) :: T
    !> Kinetics information
    type(kinetics_t), optional, intent(inout) :: theKinetics
    !> Error code
    integer, intent(out) :: Ierr
  end subroutine solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end interface

end module micm_ODE_solver
