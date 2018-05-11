module chem_solve

!----------------------------------------------------------------------
! Module which solves the chemical equations
!
! integrates dy/dt = force(y)  from t=0 to t=timeEnd
! may require jacobian = d(force(y))/dy
! ignorant of anything about y, force(y)
!----------------------------------------------------------------------

use precision, only : r8
use chemistry_specification, only: nSpecies_specified
use solver_specification, only: ndiv
use rosenbrock_integrator, only: Rosenbrock

implicit none

private
public :: chem_solve_register, chem_solve_run

integer :: nSpecies

contains

  subroutine chem_solve_register (nSpecies_loc)
    integer, intent(out) :: nSpecies_loc ! number of chemical species (NOTE -- This needs to be "protected" in cap)

    nSpecies_loc = nSpecies_specified  ! Set the number of species to pass back to the driver
    nSpecies     = nSpecies_specified  ! Set the module level nSpecies

  end subroutine chem_solve_register

  ! returns 0 for failure, 1 for success, other integers for other data
  subroutine chem_solve_run (nkReact, state_init, timeEnd, kRateConst, AbsTol, RelTol, state_final, ierr)

    integer, intent(in) :: nkReact
    real(r8), pointer, intent(in) :: state_init(:)
    real(r8), intent(in) :: timeEnd
    real(r8), pointer, intent(in) ::  kRateConst(:)   ! rates constants for each reaction
    real(r8), intent(in) :: AbsTol(:)
    real(r8), intent(in) :: RelTol(:)
    real(r8), intent(out) :: state_final(nSpecies)
    integer, intent(out) :: ierr
  
    integer  :: icntrl(20), istatus(20)
    real(r8) :: state_curr(nSpecies)
    real(r8) :: rcntrl(20), rstatus(20)
 
    real(r8) :: timeStart = 0._r8

    icntrl(:) = 0 ; rcntrl(:) = 0._r8

    icntrl(1) = 1                                 ! autonomous, F depends only on Y
    icntrl(3) = 2                                 ! ros3 solver

    state_curr(:) = state_init(:)

    ! using rosenbrock ros3 solver

    call Rosenbrock( nSpecies, state_curr, &
                     timeStart, timeEnd,  nkReact, kRateConst, AbsTol, RelTol, &
                     rcntrl, icntrl, rstatus, istatus, ierr )
    state_final(:) = state_curr(:)
   
  
  end subroutine chem_solve_run
    
  
end module chem_solve

