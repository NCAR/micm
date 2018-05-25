module chemSolve

!----------------------------------------------------------------------
! Module which solves the chemical equations
!
! integrates dy/dt = force(y)  from t=0 to t=timeEnd
! may require jacobian = d(force(y))/dy
! ignorant of anything about y, force(y)
!----------------------------------------------------------------------

use precision, only:                   r8
use chemistry_specification, only:     nSpecies_specified
use solver_specification, only:        ndiv
use rosenbrock_integrator, only:       Rosenbrock
use forcing_and_jacobian, only:        forcingParam_type

implicit none

private
public :: chemSolve_register, chemSolve_init, chemSolve_run

integer :: nSpecies

contains

  subroutine chemSolve_register (nSpecies_loc)
    integer, intent(out) :: nSpecies_loc ! number of chemical species (NOTE -- This needs to be "protected" in cap)

    nSpecies_loc = nSpecies_specified  ! Set the number of species to pass back to the driver
    nSpecies     = nSpecies_specified  ! Set the module level nSpecies

  end subroutine chemSolve_register

  subroutine chemSolve_init (absTol, relTol)
  !---------------------------------------------------------------
  ! Initialize the chemistry solver
  !---------------------------------------------------------------
    real(r8),pointer, intent(inout) :: absTol(:) ! absolute tolerance
    real(r8),pointer, intent(inout) :: relTol(:) ! relatvie tolerance

!>>FromCafe
    AbsTol(:) = 1.e-9_r8
    RelTol(:) = 1.e-4_r8
!<<FromCafe
  
  end subroutine chemSolve_init


  ! returns 0 for failure, 1 for success, other integers for other data
  subroutine chemSolve_run (nkReact, vmr_init, time_step_size, kRateConst, AbsTol, RelTol, vmr_curr, ierr)

    integer, intent(in) :: nkReact
    real(r8), pointer, intent(in) :: vmr_init(:)
    real(r8), intent(in) :: time_step_size
    real(r8), pointer, intent(in) ::  kRateConst(:)   ! rates constants for each reaction
    real(r8), intent(in) :: AbsTol(:)
    real(r8), intent(in) :: RelTol(:)
    real(r8), intent(out) :: vmr_curr(:)
    integer, intent(out) :: ierr
   
    type(forcingParam_type) :: forcingParam

    integer  :: icntrl(20), istatus(20)
    real(r8) :: rcntrl(20), rstatus(20)
 
    real(r8) :: timeStart = 0._r8
    real(r8) :: timeEnd, time

    icntrl(:) = 0 ; rcntrl(:) = 0._r8

    icntrl(1) = 1                                 ! autonomous, F depends only on Y
    icntrl(3) = 2                                 ! ros3 solver

    allocate(forcingParam%k_rateConst(nkReact))

    forcingParam%k_rateConst(:) = kRateConst(:)

    vmr_curr(:) = vmr_init(:)

    ! using rosenbrock ros3 solver
     timeEnd   = timeStart + real(ndiv,r8) * time_step_size
     time   = timeStart

    ! Advance vmr for time_step_size seconds
    do while ( time < timeEnd ) 
        call Rosenbrock( nSpecies, vmr_curr, &
                         timeStart, time_step_size,  forcingParam, AbsTol, RelTol, &
                         rcntrl, icntrl, rstatus, istatus, ierr )

        if( ierr /= 1 ) exit

        Time = Time + time_step_size
    end do
  
  end subroutine chemSolve_run
    
  
end module chemSolve

