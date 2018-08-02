
MODULE Null_Solver

  USE precision, only : r8
  USE ODE_solver

  IMPLICIT NONE
  PUBLIC
  SAVE
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Null extension
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TYPE, extends(baseOdeSolver) :: NullSolver
    INTEGER :: N
    REAL(r8), allocatable :: absTol(:), relTol(:)
    CONTAINS
      procedure :: Initialize => NullInit
      procedure :: Run        => NullRun
  END TYPE NullSolver

CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine NullInit( this, Tstart, Tend, AbsTol, RelTol,&
                       rate_frc_jac, RCNTRL,ICNTRL,IERR )

    use rate_frc_jac_module, only : rate_frc_jac_type

    class(NullSolver)    :: this
    integer, intent(out) :: ierr
    real(r8), intent(in) :: Tstart
    real(r8), intent(in) :: Tend
    REAL(r8), INTENT(IN) :: AbsTol(:),RelTol(:)
    INTEGER,  INTENT(IN) :: ICNTRL(20)
    REAL(r8), INTENT(IN) :: RCNTRL(20)
    TYPE(rate_frc_jac_type) :: rate_frc_jac

    write(*,*) ' '
    write(*,*) 'This is the Null ODE solver package'
    write(*,*) ' '

    IERR = 1

  end subroutine NullInit

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine NullRun( this, Y, Tstart, Tend, T, &
                        rate_frc_jac, istatus, rstatus, Ierr )

      use rate_frc_jac_module, only : rate_frc_jac_type

      class(NullSolver)       :: this
      integer, intent(out)    :: Ierr
      integer, intent(inout)  :: istatus(:)
      real(r8), intent(inout) :: rstatus(:)
      real(r8), intent(inout) :: Y(:)
      real(r8), intent(out)   :: T
      real(r8), intent(in)    :: Tstart
      real(r8), intent(in)    :: Tend
      TYPE(rate_frc_jac_type) :: rate_frc_jac

      integer, parameter :: Ntotstp = 9

      ISTATUS(Ntotstp) = ISTATUS(Ntotstp) + 1
      Ierr = 1

    end subroutine NullRun

END MODULE Null_Solver
