
MODULE Half_Solver

  USE ODE_solver
  USE kinetics_module, only : kinetics_type

  IMPLICIT NONE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Half extension
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TYPE, extends(baseOdeSolver) :: HalfSolver
    CONTAINS
      procedure :: Initialize => HalfInit
      procedure :: Run => HalfRun
  END TYPE HalfSolver

CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine HalfInit( this, Tstart, Tend, AbsTol, RelTol, RCNTRL, ICNTRL, IERR )

      class(HalfSolver)       :: this
      integer, intent(out)    :: Ierr

      real(r8), optional, intent(in) :: Tstart
      real(r8), optional, intent(in) :: Tend
      REAL(r8), optional, INTENT(IN) :: AbsTol(:), RelTol(:)
      INTEGER,  optional, INTENT(IN) :: ICNTRL(:)
      REAL(r8), optional, INTENT(IN) :: RCNTRL(:)

      Ierr = 0

    end subroutine HalfInit

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine HalfRun( this, Y, Tstart, Tend, T, &
                        theKinetics, istatus, rstatus, Ierr )

      class(HalfSolver)       :: this
      integer, intent(out)    :: Ierr
      real(r8), intent(inout) :: Y(:)

      integer, optional, intent(inout)  :: istatus(:)
      real(r8), optional, intent(inout) :: rstatus(:)
      real(r8), optional, intent(out)   :: T
      real(r8), optional, intent(in)    :: Tstart
      real(r8), optional, intent(in)    :: Tend
      TYPE(kinetics_type), optional     :: theKinetics

      Y(:) = Y(:) * .5_r8
      Ierr = 0

    end subroutine HalfRun

END MODULE Half_Solver
