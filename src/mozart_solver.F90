!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Backward Euler using N-R iteration
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

MODULE Mozart_Solver

  USE ODE_solver
  USE kinetics_module, only : kinetics_type

  IMPLICIT NONE

  PUBLIC
  
!~~~>  Statistics on the work performed by the Rosenbrock method
  INTEGER, PARAMETER :: Nfun=1, Njac=2, Nstp=3, Nacc=4, &
                        Nrej=5, Ndec=6, Nsol=7, Nsng=8, &
                        Ntotstp=9, Nred=10, &
                        Ntexit=1, Nhexit=2, Nhnew = 3

  REAL(r8), PARAMETER :: ZERO = 0.0_r8, HALF = .5_r8, ONE = 1.0_r8, TEN = 10._r8
  REAL(r8), PARAMETER :: FIVEPCNT = .05_r8
  REAL(r8), PARAMETER :: ErrMin   = 1.e-10_r8
  REAL(r8), PARAMETER :: DeltaMin = 1.0E-5_r8

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Rosenbrock extension
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TYPE, extends(baseOdeSolver) :: MozartSolver
    INTEGER  :: N
    INTEGER  :: iterMax
    INTEGER  :: rejMax
    REAL(r8) :: Roundoff
    REAL(r8) :: FacRej, FacAcc
    REAL(r8) :: Hmin, Hmax, Hstart
    REAL(r8), allocatable :: AbsTol(:), RelTol(:)
    logical  :: Autonomous
    logical  :: VectorTol
    CONTAINS
      procedure :: Initialize => MozartInit
      procedure :: Run        => MozartRun
      final     :: Finalize
  END TYPE MozartSolver

CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine MozartInit( this, Tstart, Tend, AbsTol, RelTol, RCNTRL, ICNTRL, IERR )

    class(MozartSolver)  :: this
    integer, intent(out) :: ierr
    real(r8), optional, intent(in) :: Tstart
    real(r8), optional, intent(in) :: Tend
    REAL(r8), optional, INTENT(IN) :: AbsTol(:),RelTol(:)
    INTEGER,  optional, INTENT(IN) :: ICNTRL(:)
    REAL(r8), optional, INTENT(IN) :: RCNTRL(:)

    integer :: i, N

    IF( .not. present(AbsTol) .or. .not. present(RelTol) .or. &
        .not. present(ICNTRL) .or. .not. present(RCNTRL) ) THEN
      Ierr = -10
      RETURN
    ENDIF

    N = size(AbsTol)
    this%N = N

    allocate( this%AbsTol(N),this%RelTol(N) )

    IERR = 0

!~~~>  Autonomous or time dependent ODE. Default is time dependent.
    this%Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
    IF (ICNTRL(2) == 0) THEN
      this%VectorTol = .TRUE.
      this%AbsTol(:N) = AbsTol(:N)
      this%RelTol(:N) = RelTol(:N)
    ELSE
      this%VectorTol = .FALSE.
      this%AbsTol(2:N) = AbsTol(1)
      this%RelTol(2:N) = RelTol(1)
    END IF
!~~~>   The maximum number of N-R iters per sub-timestep
   IF (ICNTRL(4) == 0) THEN
      this%iterMax = 11
   ELSEIF (ICNTRL(4) > 0) THEN
      this%iterMax = ICNTRL(4)
   ELSE
      IERR = -1
      PRINT * ,'User-selected max N-R iterations : ICNTRL(4)=',ICNTRL(4)
      CALL moz_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   The maximum number of consecutive sub-step failures allowed
   IF (ICNTRL(5) == 0) THEN
      this%rejMax = 5
   ELSEIF (ICNTRL(5) > 0) THEN
      this%rejMax = ICNTRL(5)
   ELSE
      IERR = -1
      PRINT * ,'User-selected max N-R failures : ICNTRL(5)=',ICNTRL(5)
      CALL moz_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Unit roundoff (1+Roundoff>1)
   this%Roundoff = EPSILON( this%Roundoff )

!~~~>  Lower bound on the step size: (positive value)
   IF (RCNTRL(1) == ZERO) THEN
      this%Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      this%Hmin = RCNTRL(1)
   ELSE
      IERR = -3
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL moz_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Upper bound on the step size: (positive value)
   IF (RCNTRL(2) == ZERO) THEN
      this%Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      this%Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      IERR = -3
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL moz_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Starting step size: (positive value)
   IF (RCNTRL(3) == ZERO) THEN
      this%Hstart = MAX(this%Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      this%Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      IERR = -3
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL moz_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      this%FacRej = HALF
   ELSEIF (RCNTRL(6) > ZERO) THEN
      this%FacRej = MIN( RCNTRL(6),HALF )
   ELSE
      IERR = -4
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL moz_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacAcc: Factor to increase H when N-R converges
   IF (RCNTRL(7) == ZERO) THEN
      this%FacAcc = 2.0_r8
   ELSEIF (RCNTRL(7) > ZERO) THEN
      this%FacAcc = MAX( RCNTRL(7),ONE )
   ELSE
      IERR = -4
      PRINT * , 'User-selected FacAcc: RCNTRL(7)=', RCNTRL(7)
      CALL moz_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Check if tolerances are reasonable
    DO i=1,N
      IF( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= TEN*this%Roundoff) &
                              .OR. (RelTol(i) >= ONE) ) THEN
        PRINT * , ' AbsTol(',i,') = ',this%AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',this%RelTol(i)
        IERR = -5
        CALL moz_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO

    if (this%print_log_message) then
       write(*,*) ' '
       write(*,*) 'Mozart ODE solver package initialized'
       write(*,*) 'Autonomous = ',this%Autonomous
       write(*,*) 'Vectol     = ',this%VectorTol
       write(*,*) 'N,Max_NR{iter,rej} = ',this%N,this%iterMax,this%rejMax
       write(*,*) 'Hmin,Hmax,Hstart     = ',this%Hmin,this%Hmax,this%Hstart
       write(*,*) 'Fac{Rej,Acc} = ',this%FacRej,this%FacAcc
       write(*,*) 'RelTol       = ',RelTol(:)
       write(*,*) 'AbsTol       = ',AbsTol(:)
       write(*,*) ' '
    endif
 
    end subroutine MozartInit

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine MozartRun( this, Y, Tstart, Tend, T, &
                          theKinetics, istatus, rstatus, Ierr )

      class(MozartSolver)     :: this
      integer, intent(out)    :: Ierr
      integer, optional, intent(inout)  :: istatus(:)
      real(r8), optional, intent(inout) :: rstatus(:)
      real(r8), intent(inout)           :: Y(:)
      real(r8), optional, intent(out)   :: T
      real(r8), optional, intent(in)    :: Tstart
      real(r8), optional, intent(in)    :: Tend
      TYPE(kinetics_type), optional     :: theKinetics

! ~~~~ Local variables
      INTEGER  :: N
      INTEGER  :: nIter
      INTEGER  :: j
      INTEGER  :: rejCnt
      INTEGER  :: istat(20)
      REAL(r8) :: H, Hinv, Hnew
      REAL(r8) :: presentTime
      REAL(r8) :: truncError
      REAL(r8), target  :: residual(this%N)
      REAL(r8), pointer :: deltaY(:)
      REAL(r8) :: Ynew(this%N), Yerr(this%N)
      REAL(r8) :: Fcn(this%N)
      REAL(r8) :: rstat(20)
      LOGICAL  :: RejectLastH, Singular, nrConverged

!~~~>  Initial preparations
   IF( .not. present(theKinetics) .or. .not. present(Tstart) .or. &
       .not. present(Tend) ) THEN
     Ierr = -10
     RETURN
   ENDIF

   N = this%n
   presentTime = Tstart
   deltaY => residual
   RSTAT(Nhexit) = ZERO
   H = MIN( MAX(ABS(this%Hmin),ABS(this%Hstart)) , ABS(this%Hmax) )
   IF (ABS(H) <= TEN*this%Roundoff) H = DeltaMin

   Ierr   = 0
   rejCnt = 0
   rejectLastH = .FALSE.
   istat(:) = 0
   rstat(:) = ZERO

!~~~> Time loop begins below
TimeLoop: DO WHILE ( (presentTime-Tend)+this%Roundoff <= ZERO)

     IF ( ((presentTime+0.1_r8*H) == presentTime).OR.(H <= this%Roundoff) ) THEN  ! Step size too small
       Ierr = -7
       CALL moz_ErrorMsg(-7,presentTime,H,IERR)
       RETURN
     END IF

!~~~>  Limit H if necessary to avoid going beyond Tend
     H = MIN( H,ABS(Tend-presentTime) )
     Hinv = ONE/H
     Ynew(1:N) = Y(1:N)
!~~~>  Newton-Raphson iteration loop
NRloop: DO nIter = 1,this%iterMax
!~~~>   Compute the Jacobian
       ISTAT(Njac) = ISTAT(Njac) + 1
       CALL theKinetics%LinFactor( H, ONE, Ynew, Singular, istat )
       IF (Singular) THEN ! More than 5 consecutive failed decompositions
         Ierr = -8
         CALL moz_ErrorMsg(-8,presentTime,H,IERR)
         RETURN
       END IF
!~~~>   Compute the function at current time
       Fcn(:) = theKinetics%force( Ynew )
       ISTAT(Nfun) = ISTAT(Nfun) + 1
       residual(1:N) = Fcn(1:N) - (Ynew(1:N) - Y(1:N))*Hinv
!~~~>   Compute the iteration delta
       CALL theKinetics%LinSolve( deltaY )
!~~~>   Update N-R iterate
       Ynew(1:N) = Ynew(1:N) + deltaY(1:N)

       if( nIter > 1 ) then
         nrConverged = moz_NRErrorNorm( this, deltaY, Ynew, Yerr )
         if( nrConverged ) then
           EXIT NRloop
         endif
       endif
     END DO NRloop

nrHasConverged: &
     IF( nrConverged ) THEN
       Yerr(:) = theKinetics%dForcedyxForce( Fcn )
       truncError = moz_TruncErrorNorm ( this, Y, Ynew, Yerr )
       Hnew = SQRT( 2.0_r8/truncError )
       Hnew = MAX(this%Hmin,MIN(Hnew,this%Hmax))
!~~~>   Check step truncation error
acceptStep: &
!      IF( truncError <= 2.0_r8/(H*H) ) THEN
       IF( .5_r8*H*H*truncError <= (ONE + 100._r8*this%Roundoff) ) THEN
         !write(*,*) 'mozartRun: step accepted'
         Y(1:N) = Ynew(1:N)
         presentTime = presentTime + H
         IF (RejectLastH) THEN  ! No step size increase after a rejected step
           Hnew = MIN(Hnew,H)
         END IF
         RSTAT(Nhexit) = H
         RSTAT(Nhnew)  = Hnew
         RSTAT(Ntexit) = presentTime
         rejCnt = 0
         RejectLastH    = .FALSE.
         ISTAT(Nstp)    = ISTAT(Nstp) + 1
         ISTAT(Ntotstp) = ISTAT(Ntotstp) + 1
         ISTAT(Nacc)    = ISTAT(Nacc) + 1
       ELSE acceptStep
         IF( Hnew == this%Hmin .and. H == this%Hmin ) THEN
           write(*,'(a)') 'Warning: H == Hmin, can not be reduced '
           Ierr = -1
           EXIT TimeLoop
         ENDIF
         IF( RejectLastH ) then
           IF( abs(H - Hnew) <= FIVEPCNT*H ) THEN
             Hnew = .9_r8 * Hnew
             !write(*,'(''mozartRun: Hnew = '',1p,1x,g0)') Hnew
           ENDIF
         ENDIF
         IF (ISTAT(Nacc) >= 1)  ISTAT(Nrej) = ISTAT(Nrej) + 1
         RejectLastH = .TRUE.
       ENDIF acceptStep
     ELSE nrHasConverged
!~~~> N-R failed to converge
       write(*,*) 'mozartRun: Newton-Raphson fails to converge for within timestep allowed'
       write(*,'(''time, chemistry-subtimestep, iter = '',1p,2e12.4,i4)') presentTime,H,nIter
       rejCnt = rejCnt + 1
       IF( rejCnt < this%rejMax ) THEN
         Hnew = MAX( H*this%FacRej,this%Hmin )
         IF( rejCnt > 1 .and. Hnew == this%Hmin .and. H == this%Hmin ) THEN
           write(*,'(a)') 'Warning: H == Hmin, can not be reduced '
           Ierr = -1
           EXIT TimeLoop
         ENDIF
         ISTAT(Nred) = ISTAT(Nred) + 1
       ELSE
         write(*,'(a,i2)') 'Warning: Last ',this%rejMax,' sub steps failed'
         Ierr = -1
         EXIT TimeLoop
       ENDIF
     ENDIF nrHasConverged
     H = Hnew
   END DO TimeLoop

   IF( present(T) ) THEN
     T = presentTime
   ENDIF
   IF( present(istatus) ) THEN
     istatus(:) = istat(:)
   ENDIF
   IF( present(rstatus) ) THEN
     rstatus(:) = rstat(:)
   ENDIF

!   stop 'Debugging'

    end subroutine MozartRun

    SUBROUTINE moz_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
      REAL(r8), INTENT(IN) :: T, H
      INTEGER, INTENT(IN)  :: Code
      INTEGER, INTENT(OUT) :: IERR
   
      IERR = Code
      PRINT * , &
       'Forced exit from Mozart due to the following error:' 
     
   SELECT CASE (Code)
    CASE (-1)    
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-3)    
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4)    
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5) 
      PRINT * , '--> Improper tolerance values'
    CASE (-6) 
      PRINT * , '--> No of steps exceeds maximum bound'
    CASE (-7) 
      PRINT * , '--> Step size too small: T + 10*H = T or H < Roundoff'
    CASE (-8)    
      PRINT * , '--> Matrix is repeatedly singular'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT
   
   WRITE(*,'(a,1p,g0,a,g0)') "Time = ", T, " and H = ", H
     
 END SUBROUTINE moz_ErrorMsg

  FUNCTION moz_NRErrorNorm ( this, deltaY, Ynew, Yerr ) result( Converged )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Input arguments
   class(mozartSolver)  :: this
   REAL(r8), INTENT(IN) :: deltaY(:), Ynew(:), Yerr(:)

   LOGICAL :: Msk(size(deltaY))
   LOGICAL :: Converged

   Msk(:) = Ynew(:) > this%AbsTol(:)
!  Converged = .not. ANY( ABS(deltaY(:)) > this%RelTol(:)*ABS(Ynew(:)) )
   Converged = .not. ANY( ABS(deltaY(:)) > this%RelTol(:)*ABS(Ynew(:)) .and. Msk(:) )

  END FUNCTION moz_NRErrorNorm

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  FUNCTION moz_TruncErrorNorm ( this, Y, Ynew, Yerr ) result( Error )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Input arguments
   class(MozartSolver)  :: this
   REAL(r8), INTENT(IN) :: Y(:), Ynew(:), Yerr(:)
   REAL(r8)             :: Error

! Local variables
   REAL(r8) :: Scale(this%N), Ymax(this%N)
   REAL(r8) :: Ytmp(this%N)

   Ymax(:)  = MAX( ABS(Y(:)),ABS(Ynew(:)) )
!  Scale(:) = this%AbsTol(:) + this%RelTol(:)*Ymax(:)
   Scale(:) = this%RelTol(:)*Ymax(:)
   WHERE( Ymax(:) >= this%AbsTol(:) )
     Ytmp(:) = Yerr(:)
   ELSEWHERE
     Ytmp(:) = ZERO
   ENDWHERE
   Error    = MAX( SQRT( sum( (Ytmp(:)/Scale(:))**2 )/real(this%N,kind=r8) ),ErrMin )

  END FUNCTION moz_TruncErrorNorm

  subroutine Finalize( this )

  type(MozartSolver) :: this

  if( allocated( this%absTol ) ) then
    deallocate( this%absTol )
  endif
  if( allocated( this%relTol ) ) then
    deallocate( this%relTol )
  endif

  end subroutine Finalize

END MODULE Mozart_Solver
