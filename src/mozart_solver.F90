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
      this%FacRej = 0.5_r8
   ELSEIF (RCNTRL(6) > ZERO) THEN
      this%FacRej = MIN( RCNTRL(6),0.5_r8 )
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

    write(*,*) ' '
    write(*,*) 'Mozart ODE solver package initialized'
    write(*,*) 'Autonomous = ',this%Autonomous
    write(*,*) 'Vectol     = ',this%VectorTol
    write(*,*) 'N,Max_NR{iter,rej} = ',this%N,this%iterMax,this%rejMax
    write(*,*) 'Hmin,Hmax,Hstart     = ',this%Hmin,this%Hmax,this%Hstart
    write(*,*) 'Fac{Rej,Acc} = ',this%FacRej,this%FacAcc
    write(*,*) 'RelTol       = ',RelTol(:)
    write(*,*) ' '

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
      INTEGER  :: Direction, j
      INTEGER  :: rejCnt
      INTEGER  :: Pivot(this%N)
      INTEGER  :: istat(20)
      REAL(r8) :: H, Hinv, Hnew
      REAL(r8) :: presentTime
      REAL(r8) :: truncError
      REAL(r8), target  :: residual(this%N)
      REAL(r8), pointer :: deltaY(:)
      REAL(r8) :: Ynew(this%N), Yerr(this%N)
      REAL(r8) :: Fcn0(this%N)
      REAL(r8) :: Jac0(this%N,this%N), Ghimj(this%N,this%N)
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

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   H = Direction*H

   Ierr   = 0
   rejCnt = 0
   rejectLastH = .TRUE.
   istat(:) = 0
   rstat(:) = ZERO

!~~~> Time loop begins below
TimeLoop: DO WHILE ( (Direction > 0).AND.((presentTime-Tend)+this%Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-presentTime)+this%Roundoff <= ZERO) )

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
       CALL JacTemplate( N, presentTime, Ynew, Jac0, theKinetics )
       ISTAT(Njac) = ISTAT(Njac) + 1
       CALL moz_PrepareMatrix( N, Hinv, Direction, Jac0, Ghimj, &
                               Pivot, Singular, istat(Nsng), istat(Ndec) )
       IF (Singular) THEN ! More than 5 consecutive failed decompositions
         Ierr = -8
         CALL moz_ErrorMsg(-8,presentTime,H,IERR)
         RETURN
       END IF
!~~~>   Compute the function at current time
       CALL FunTemplate( N, presentTime, Ynew, Fcn0, theKinetics )
       ISTAT(Nfun) = ISTAT(Nfun) + 1
       residual(1:N) = Fcn0(1:N) - (Ynew(1:N) - Y(1:N))*Hinv
!~~~>   Compute the iteration delta
       CALL moz_Solve( N, Ghimj, Pivot, deltaY, istat(Nsol) )
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
       Yerr(:) = MATMUL( Jac0,Fcn0 )
       truncError = moz_TruncErrorNorm ( this, Y, Ynew, Yerr )
       Hnew = SQRT( 2.0_r8/truncError )
       Hnew = MAX(this%Hmin,MIN(Hnew,this%Hmax))
!~~~>   Check step truncation error
acceptStep: &
       IF( truncError <= 2.0_r8/(H*H) ) THEN
         Y(1:N) = Ynew(1:N)
         presentTime = presentTime + Direction*H
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
         IF (ISTAT(Nacc) >= 1)  ISTAT(Nrej) = ISTAT(Nrej) + 1
         RejectLastH = .TRUE.
       ENDIF acceptStep
     ELSE nrHasConverged
!~~~> N-R failed to converge
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

SUBROUTINE FunTemplate( N, T, Y, Ydot, theKinetics )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the ODE function call.
!  Updates the rate coefficients (and possibly the fixed species) at each call
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~> Input variables
   INTEGER, intent(in) :: N
   REAL(r8), intent(in) :: T
   REAL(r8), intent(in) :: Y(N)
!~~~> Output variables
   REAL(r8), intent(out) :: Ydot(N)

   TYPE(kinetics_type) :: theKinetics

   Ydot(:) =  theKinetics%force( Y )

END SUBROUTINE FunTemplate

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE JacTemplate( N, T, Y, Jcb, theKinetics )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the ODE Jacobian call.
!  Updates the rate coefficients (and possibly the fixed species) at each call
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~> Input variables
   INTEGER, intent(in) :: N
   REAL(r8), intent(in) :: T                     ! time
   REAL(r8), intent(in) :: Y(N)
!~~~> Output variables
   REAL(r8), intent(inout) :: Jcb(N,N)
   TYPE(kinetics_type) :: theKinetics

    Jcb(:,:) = theKinetics%Jac( Y )

END SUBROUTINE JacTemplate

  SUBROUTINE moz_PrepareMatrix ( N, Hinv, Direction, Jac0, Ghimj, &
                                 Pivot, Singular, Nsng, Ndec )
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/H - Jac0
! --- --- --- --- --- --- --- --- --- --- --- --- ---

!~~~> Input arguments
   INTEGER, INTENT(IN)  ::  N
   INTEGER, INTENT(IN)  ::  Direction
   INTEGER, INTENT(OUT) ::  Pivot(N)
   REAL(r8), INTENT(IN) ::  Hinv   ! inverse step size
   INTEGER, INTENT(INOUT) ::  Nsng, Ndec
   REAL(r8), INTENT(IN) ::  Jac0(N,N)
!~~~> Output arguments
   REAL(r8), INTENT(OUT) :: Ghimj(N,N)
   LOGICAL, INTENT(OUT)  :: Singular

   INTEGER  :: i, ISING
   REAL(r8) :: ghinv

!~~~>    Construct Ghimj = 1/H - Jac0
   Ghimj(:,:) = -Jac0(:,:)
   ghinv = Direction*Hinv
   FORALL( I = 1:N ) 
     Ghimj(i,i) = Ghimj(i,i) + ghinv
   ENDFORALL
!~~~>    Compute LU decomposition
   CALL moz_Decomp( N, Ghimj, Pivot, ISING, Ndec )
   IF (ISING == 0) THEN
!~~~>    If successful done
     Singular = .FALSE.
   ELSE ! ISING .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
     Nsng = Nsng + 1
     Singular = .TRUE.
     PRINT*,'Warning: LU Decomposition returned ISING = ',ISING
   ENDIF

  END SUBROUTINE moz_PrepareMatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
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
   Scale(:) = this%AbsTol(:) + this%RelTol(:)*Ymax(:)
   WHERE( Ymax(:) >= this%AbsTol(:) )
     Ytmp(:) = Yerr(:)
   ELSEWHERE
     Ytmp(:) = 0.0_r8 
   ENDWHERE
   Error    = MAX( SQRT( sum( (Ytmp(:)/Scale(:))**2 )/real(this%N,kind=r8) ),ErrMin )

  END FUNCTION moz_TruncErrorNorm

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE moz_Decomp( N, A, Pivot, ISING, Ndec )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the LU decomposition
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~> Inout variables
   INTEGER, INTENT(IN) :: N
   INTEGER, INTENT(INOUT) :: Ndec
   REAL(r8), INTENT(INOUT) :: A(N,N)
!~~~> Output variables
   INTEGER, INTENT(OUT) :: ISING
   INTEGER, INTENT(OUT) :: Pivot(N)

   CALL DGEFA( A, N, Pivot, ISING )
   Ndec = Ndec + 1
   IF ( ISING /= 0 ) THEN
      PRINT*,"Error in DGEFA. ISING=",ISING
   END IF  

  END SUBROUTINE moz_Decomp


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE moz_Solve( N, A, Pivot, b, Nsol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the forward/backward substitution (using pre-computed LU decomposition)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~> Input variables
   INTEGER, INTENT(IN)  :: N
   INTEGER, INTENT(INOUT)  :: Nsol
   REAL(r8), INTENT(IN) :: A(N,N)
   INTEGER, INTENT(IN)  :: Pivot(N)
!~~~> InOut variables
   REAL(r8), INTENT(INOUT) :: b(N)

   CALL DGESL( A, N, Pivot, b, 0 )

   Nsol = Nsol + 1

  END SUBROUTINE moz_Solve

    SUBROUTINE DGEFA (A, N, IPVT, INFO)
!***BEGIN PROLOGUE  DGEFA
!***PURPOSE  Factor a matrix using Gaussian elimination.
!***CATEGORY  D2A1
!***TYPE      DOUBLE PRECISION (SGEFA-S, DGEFA-D, CGEFA-C)
!***KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
!             MATRIX FACTORIZATION
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DGEFA factors a double precision matrix by Gaussian elimination.
!
!     DGEFA is usually called by DGECO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .
!
!     On Entry
!
!        A       DOUBLE PRECISION(N, N)
!                the matrix to be factored.
!
!        N       INTEGER
!                the order of the matrix  A .
!
!     On Return
!
!        A       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.
!
!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.
!
!        INFO    INTEGER
!                = 0  normal value.
!                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
!                     condition for this subroutine, but it does
!                     indicate that DGESL or DGEDI will divide by zero
!                     if called.  Use  RCOND  in DGECO for a reliable
!                     indication of singularity.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DGEFA

!
!     DUMMY ARGUMENTS
!
      INTEGER, INTENT(IN)    ::  N
      INTEGER, INTENT(OUT)   ::  INFO
      INTEGER, INTENT(INOUT) ::  IPVT(N)

      REAL(r8), INTENT(INOUT) :: A(N,N)
!
!     LOCAL VARIABLES
!
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.D0
      DOUBLE PRECISION, PARAMETER :: ONE  = 1.D0

      INTEGER :: J, K, KP1, NM1
      INTEGER, POINTER :: L
      INTEGER, TARGET :: MAXNDX(1)
      DOUBLE PRECISION :: T
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!

      INFO = 0
      NM1 = N - 1
      IF( N >= 2 ) THEN
        L => MAXNDX(1)
        DO K = 1, NM1
          KP1 = K + 1
!
!        FIND L = PIVOT INDEX
!
          MAXNDX(:) = MAXLOC( ABS(A(K:N,K)) ) + K - 1
          IPVT(K) = L
!
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!
          IF (A(L,K) /= ZERO) THEN
!
!           INTERCHANGE IF NECESSARY
!
             IF (L /= K) THEN
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
             ENDIF
!
!           COMPUTE MULTIPLIERS
!
             T = -ONE/A(K,K)
             A(K+1:N,K) = T*A(K+1:N,K)
!
!           ROW ELIMINATION WITH COLUMN INDEXING
!
             DO J = KP1, N
               T = A(L,J)
               IF (L /= K) THEN
                 A(L,J) = A(K,J)
                 A(K,J) = T
               ENDIF   
               A(K+1:N,J) = A(K+1:N,J) + T*A(K+1:N,K)
             END DO
          ELSE
           INFO = K
          ENDIF
        END DO  
      ELSE
        IPVT(N) = N
        IF(A(N,N) == ZERO ) INFO = N
      ENDIF

      END SUBROUTINE DGEFA

      SUBROUTINE DGESL (A, N, IPVT, B, JOB)
!***BEGIN PROLOGUE  DGESL
!***PURPOSE  Solve the real system A*X=B or TRANS(A)*X=B using the
!            factors computed by DGECO or DGEFA.
!***CATEGORY  D2A1
!***TYPE      DOUBLE PRECISION (SGESL-S, DGESL-D, CGESL-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DGESL solves the double precision system
!     A * X = B  or  TRANS(A) * X = B
!     using the factors computed by DGECO or DGEFA.
!
!     On Entry
!
!        A       DOUBLE PRECISION(N, N)
!                the output from DGECO or DGEFA.
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        IPVT    INTEGER(N)
!                the pivot vector from DGECO or DGEFA.
!
!        B       DOUBLE PRECISION(N)
!                the right hand side vector.
!
!        JOB     INTEGER
!                = 0         to solve  A*X = B ,
!                = nonzero   to solve  TRANS(A)*X = B  where
!                            TRANS(A)  is the transpose.
!
!     On Return
!
!        B       the solution vector  X .
!
!     Error Condition
!
!        A division by zero will occur if the input factor contains a
!        zero on the diagonal.  Technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of LDA .  It will not occur if the subroutines are
!        called correctly and if DGECO has set RCOND .GT. 0.0
!        or DGEFA has set INFO .EQ. 0 .
!
!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
!           IF (RCOND is too small) GO TO ...
!           DO 10 J = 1, P
!              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
!        10 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DGESL

      INTEGER, INTENT(IN)   ::  N
      INTEGER, INTENT(IN)   ::  JOB
      INTEGER, INTENT(IN)   ::  IPVT(N)

      REAL(r8), INTENT(IN)    :: A(N,N)
      REAL(r8), INTENT(INOUT) :: B(N)

!
      INTEGER :: K, KM1, KP1, KB, L, NM1
      DOUBLE PRECISION :: T

      NM1 = N - 1
      IF (JOB == 0) THEN
!
!        JOB = 0 , SOLVE  A * X = B
!        FIRST SOLVE  L*Y = B
!
         IF (N >= 2) THEN
           DO K = 1, NM1
             L = IPVT(K)
             T = B(L)
             IF (L /= K) THEN
               B(L) = B(K)
               B(K) = T
            ENDIF
            B(K+1:N) = B(K+1:N) + T*A(K+1:N,K)
           END DO
         ENDIF
!
!        NOW SOLVE  U*X = Y
!
         DO KB = 1, N
           K = N + 1 - KB
           KM1 = K - 1
           B(K) = B(K)/A(K,K)
           T = -B(K)
           B(1:KM1) = B(1:KM1) + T*A(1:KM1,K)
         END DO
      ELSE
!
!        JOB = NONZERO, SOLVE  TRANS(A) * X = B
!        FIRST SOLVE  TRANS(U)*Y = B
!
         DO K = 1, N
           KM1 = K - 1
           T = DOT_PRODUCT( A(1:KM1,K),B(1:KM1) )
           B(K) = (B(K) - T)/A(K,K)
         END DO
!
!        NOW SOLVE TRANS(L)*X = Y
!
         IF (NM1 > 0) THEN
           DO KB = 1, NM1
             K = N - KB
             KP1 = K + 1
             B(K) = B(K) + DOT_PRODUCT( A(KP1:N,K),B(KP1:N) )
             L = IPVT(K)
             IF (L /= K) THEN
               T = B(L)
               B(L) = B(K)
               B(K) = T
             ENDIF
           END DO
         ENDIF
      ENDIF   

      END SUBROUTINE DGESL

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
