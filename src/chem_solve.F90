module chemSolve

!----------------------------------------------------------------------
! Module which solves the chemical equations
!
! integrates dy/dt = force(y)  from t=0 to t=timeEnd
! may require jacobian = d(force(y))/dy
! ignorant of anything about y, force(y)
!----------------------------------------------------------------------

use precision,                   only: r8
use chemistry_specification,     only: nSpecies_specified
use solver_specification,        only: ndiv
use rosenbrock_integrator,       only: Rosenbrock
use forcing_and_jacobian,        only: forcingParam_type

implicit none

private
public :: chemSolve_register, chemSolve_init, chemSolve_run

integer :: nSpecies

contains

  subroutine chemSolve_register (nSpecies_loc)
  !---------------------------------------------------------------
  ! Register the chemistry constants (from the Cafe)
  !---------------------------------------------------------------

    integer, intent(out) :: nSpecies_loc ! number of chemical species (NOTE -- This needs to be "protected" in cap)

    nSpecies_loc = nSpecies_specified  ! Set the number of species to pass back to the driver
    nSpecies     = nSpecies_specified  ! Set the module level nSpecies

  end subroutine chemSolve_register

  subroutine chemSolve_init (absTol, relTol)
  use chemistry_specification, only: chemistry_init
  !---------------------------------------------------------------
  ! Initialize the chemistry solver
  !---------------------------------------------------------------
    real(r8),pointer, intent(inout) :: absTol(:) ! absolute tolerance
    real(r8),pointer, intent(inout) :: relTol(:) ! relatvie tolerance

    call chemistry_init(absTol, relTol)
  
  end subroutine chemSolve_init


  ! returns 0 for failure, 1 for success, other integers for other data
  subroutine chemSolve_run (nkReact, vmr_init, timeStepSize, kRateConst, absTol, relTol, vmr_final, ierr)
  !---------------------------------------------------------------
  ! Run the chemistry solver
  !---------------------------------------------------------------

    integer, intent(in) :: nkReact
    real(r8), pointer, intent(in) :: vmr_init(:)
    real(r8), intent(in) :: timeStepSize
    real(r8), pointer, intent(in) ::  kRateConst(:)   ! rates constants for each reaction
    real(r8), intent(in) :: absTol(:)
    real(r8), intent(in) :: relTol(:)
    real(r8), intent(out) :: vmr_final(:)             ! Final VMR
    integer, intent(out) :: ierr
   
    type(forcingParam_type) :: forcingParam
    real(r8)                :: vmr_curr(size(vmr_final))      ! Value of current VMR

    integer  :: icntrl(20), istatus(20)
    real(r8) :: rcntrl(20), rstatus(20)
 
    real(r8) :: timeStart = 0._r8
    real(r8) :: timeEnd, time

    icntrl(:) = 0 ; rcntrl(:) = 0._r8

    icntrl(1) = 1                                 ! autonomous, F depends only on Y
    icntrl(3) = 2                                 ! ros3 solver

    forcingParam%nkReact = nkReact
    allocate(forcingParam%k_rateConst(nkReact))

    ! Allocate forcingParam arrays for use internally in rosenbrock
    forcingParam%k_rateConst(:) = kRateConst(:)

    ! Start the VMR array with the initial setting
    vmr_curr(:) = vmr_init(:)

    ! Initialize the time keeping variables
    timeEnd   = timeStart + real(ndiv,r8) * timeStepSize
    time   = timeStart

    ! Advance vmr for timeStepSize seconds
    do while ( time < timeEnd ) 

       ! using rosenbrock ros3 solver
        call Rosenbrock( nSpecies, vmr_curr, &
                         timeStart, timeStepSize,  forcingParam, absTol, relTol, &
                         rcntrl, icntrl, rstatus, istatus, ierr )

        if( ierr /= 1 ) exit

        ! Status was good, proceed with next time step
        vmr_final(:) = vmr_curr(:)
        Time = Time + timeStepSize

    end do
  
    deallocate(forcingParam%k_rateConst)

  end subroutine chemSolve_run
    
  
end module chemSolve

