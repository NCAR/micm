module ode_solve

! integrates dy/dt = force(y)  from t=0 to t=time_step_size
! may require jacobian = d(force(y))/dy
! ignorant of anything about y, force(y)

use precision, only : r8
use chemistry_specification, only: n_prognostic_variables
use solver_specification, only: ndiv
use rosenbrock_integrator, only: Rosenbrock

implicit none

contains

  subroutine time_advance( initial_state, time_step_size, final_state, &
                           AbsTol, RelTol, ierr)
  ! returns 0 for failure, 1 for success, other integers for other data

    integer, intent(inout) :: ierr
    real(r8), intent(in)   :: initial_state(n_prognostic_variables)
    real(r8), intent(in)   :: time_step_size
    real(r8), intent(in)   :: AbsTol(n_prognostic_variables)
    real(r8), intent(in)   :: RelTol(n_prognostic_variables)
    real(r8), intent(out)  :: final_state(n_prognostic_variables)
  
    integer  :: icntrl(20), istatus(20)
    real(r8) :: state_at_current_time_step(n_prognostic_variables)
    real(r8) :: rcntrl(20), rstatus(20)

    icntrl(:) = 0 ; rcntrl(:) = 0._r8

    icntrl(1) = 1                                 ! autonomous, F depends only on Y
    icntrl(3) = 2                                 ! ros3 solver
!   rcntrl(3) = .01_r8 * time_step_size           ! initial rosenbrock step size

    state_at_current_time_step(:) = initial_state(:)

    ! using rosenbrock ros3 solver
    call Rosenbrock( n_prognostic_variables, state_at_current_time_step, &
                     0._r8, time_step_size, AbsTol, RelTol, &
                     rcntrl, icntrl, rstatus, istatus, ierr )
  
    final_state(:) = state_at_current_time_step(:)
   
  end subroutine time_advance
  
end module ode_solve

