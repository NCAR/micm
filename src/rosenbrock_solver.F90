!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Rosenbrock - Implementation of several Rosenbrock methods:             !
!               * Ros2                                                    !
!               * Ros3                                                    !
!               * Ros4                                                    !
!               * Rodas3                                                  !
!               * Rodas4                                                  !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

MODULE Rosenbrock_Solver

  USE ODE_solver
  use kinetics_module, only : kinetics_type

  IMPLICIT NONE

  PUBLIC
  
!~~~>  Statistics on the work performed by the Rosenbrock method
  INTEGER, PARAMETER :: Nfun=1, Njac=2, Nstp=3, Nacc=4, &
                        Nrej=5, Ndec=6, Nsol=7, Nsng=8, &
                        Ntotstp= 9, &
                        Ntexit=1, Nhexit=2, Nhnew = 3

  INTEGER, PARAMETER  :: RS2=1, RS3=2, RS4=3, RD3=4, RD4=5, RG3=6
  REAL(r8), PARAMETER :: ZERO = 0.0_r8, HALF = .5_r8, ONE = 1.0_r8, TEN = 10._r8
  REAL(r8), PARAMETER :: ErrMin   = 1.e-10_r8
  REAL(r8), PARAMETER :: DeltaMin = 1.0E-5_r8

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Rosenbrock extension
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TYPE, extends(baseOdeSolver) :: RosenbrockSolver
    INTEGER  :: N
    INTEGER  :: ros_S, rosMethod
    INTEGER  :: UplimTol, Max_no_steps
    REAL(r8) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
    REAL(r8) :: Hmin, Hmax, Hstart
    REAL(r8) :: ros_ELO
    REAL(r8) :: ros_A(15), ros_C(15), ros_M(6), ros_E(6), &
                ros_Alpha(6), ros_Gamma(6)
    REAL(r8), allocatable :: AbsTol(:), RelTol(:)
    logical  :: Autonomous
    logical  :: VectorTol
    LOGICAL  :: ros_NewF(6)
    CHARACTER(LEN=12) :: ros_Name
    CONTAINS
      procedure :: Initialize => RosenbrockInit
      procedure :: Run        => RosenbrockRun
      final     :: Finalize
  END TYPE RosenbrockSolver

CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RosenbrockInit( this, Tstart, Tend, AbsTol, RelTol, RCNTRL, ICNTRL, IERR)

    class(RosenbrockSolver) :: this
    integer, intent(out) :: ierr
    real(r8), optional, intent(in) :: Tstart
    real(r8), optional, intent(in) :: Tend
    REAL(r8), optional, INTENT(IN) :: AbsTol(:),RelTol(:)
    INTEGER,  optional, INTENT(IN) :: ICNTRL(:)
    REAL(r8), optional, INTENT(IN) :: RCNTRL(:)

!   local variables
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
      this%UplimTol  = N
      this%AbsTol(:N) = AbsTol(:N)
      this%RelTol(:N) = RelTol(:N)
    ELSE
      this%VectorTol = .FALSE.
      this%UplimTol  = 1
      this%AbsTol(2:N) = AbsTol(1)
      this%RelTol(2:N) = RelTol(1)
    END IF
!~~~>   Initialize the particular Rosenbrock method selected
    SELECT CASE (ICNTRL(3))
      CASE (1)
       CALL Ros2( this )
      CASE (2)
       CALL Ros3( this )
      CASE (3)
       CALL Ros4( this )
      CASE (0,4)
       CALL Rodas3( this )
      CASE (5)
       CALL Rodas4( this )
      CASE (6)
       CALL Rang3( this )
      CASE DEFAULT
       IERR = -2
       WRITE(*,*) 'Unknown Rosenbrock method: ICNTRL(3) = ',ICNTRL(3) 
       CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
    END SELECT

!~~~>   The maximum number of steps admitted
   IF (ICNTRL(4) == 0) THEN
      this%Max_no_steps = 200000
   ELSEIF (ICNTRL(4) > 0) THEN
      this%Max_no_steps=ICNTRL(4)
   ELSE
      IERR = -1
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
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
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Upper bound on the step size: (positive value)
   IF( .not. present(Tend) .or. .not. present(Tstart) ) THEN
     Ierr = -10
     RETURN
   ENDIF
   IF (RCNTRL(2) == ZERO) THEN
      this%Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      this%Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      IERR = -3
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
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
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax
   IF (RCNTRL(4) == ZERO) THEN
      this%FacMin = 0.2_r8
   ELSEIF (RCNTRL(4) > ZERO) THEN
      this%FacMin = RCNTRL(4)
   ELSE
      IERR = -4
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      this%FacMax = 6.0_r8
   ELSEIF (RCNTRL(5) > ZERO) THEN
      this%FacMax = RCNTRL(5)
   ELSE
      IERR = -4
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      this%FacRej = 0.1_r8
   ELSEIF (RCNTRL(6) > ZERO) THEN
      this%FacRej = RCNTRL(6)
   ELSE
      IERR = -4
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
      this%FacSafe = 0.9_r8
   ELSEIF (RCNTRL(7) > ZERO) THEN
      this%FacSafe = RCNTRL(7)
   ELSE
      IERR = -4
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Check if tolerances are reasonable
    DO i = 1,N
      IF( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= TEN*this%Roundoff) &
                              .OR. (RelTol(i) >= ONE) ) THEN
        PRINT * , ' AbsTol(',i,') = ',this%AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',this%RelTol(i)
        IERR = -5
        CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO

    write(*,*) ' '
    write(*,*) 'Rosenbrock ODE solver package initialized'
    write(*,*) 'Autonomous = ',this%Autonomous
    write(*,*) 'Vectol     = ',this%VectorTol
    write(*,*) 'N,ros_S,Max_no_steps = ',this%N,this%ros_S,this%Max_no_steps
    write(*,*) 'Hmin,Hmax,Hstart     = ',this%Hmin,this%Hmax,this%Hstart
    write(*,*) 'Fac{Min,Max,Rej,Safe} = ',this%FacMin,this%FacMax,this%FacRej,this%FacSafe
    write(*,*) ' '

    end subroutine RosenbrockInit

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine RosenbrockRun( this, Y, Tstart, Tend, T, &
                              theKinetics, istatus, rstatus, Ierr )

      class(RosenbrockSolver) :: this
      integer, intent(out)    :: Ierr
      real(r8), intent(inout)           :: Y(:)
      integer, optional, intent(inout)  :: istatus(:)
      real(r8), optional, intent(inout) :: rstatus(:)
      real(r8), optional, intent(out)   :: T
      real(r8), optional, intent(in)    :: Tstart
      real(r8), optional, intent(in)    :: Tend
      TYPE(kinetics_type), optional     :: theKinetics

! ~~~~ Local variables
      INTEGER  :: N
      INTEGER  :: S_ndx
      INTEGER  :: Direction, ioffset, j, istage
      integer  :: istat(20)
      INTEGER  :: Pivot(this%N)
      REAL(r8) :: H, Hnew, HC, HG, Fac, Tau, Err
      REAL(r8) :: presentTime
      REAL(r8) :: Ynew(this%N)
      REAL(r8) :: Fcn0(this%N), Fcn(this%N), dFdT(this%N)
      REAL(r8) :: K(this%N,this%ros_S)
      REAL(r8) :: Jac0(this%N,this%N), Ghimj(this%N,this%N)
      REAL(r8) :: Yerr(this%N)
      real(r8) :: rstat(20)
      LOGICAL  :: RejectLastH, RejectMoreH, Singular

!~~~>  Initial preparations
   IF( .not. present(theKinetics) .or. .not. present(Tstart) .or. &
       .not. present(Tend) ) THEN
     Ierr = -10
     RETURN
   ENDIF

   N = this%n
   presentTime = Tstart
   RSTAT(Nhexit) = ZERO
   H = MIN( MAX(ABS(this%Hmin),ABS(this%Hstart)) , ABS(this%Hmax) )
   IF (ABS(H) <= TEN*this%Roundoff) H = DeltaMin

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.

   istat(:) = 0
   rstat(:) = ZERO

!~~~> Time loop begins below

TimeLoop: DO WHILE ( (Direction > 0).AND.((presentTime-Tend)+this%Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-presentTime)+this%Roundoff <= ZERO) )

   IF ( istat(Nstp) > this%Max_no_steps ) THEN  ! Too many steps
      Ierr = -6
      CALL ros_ErrorMsg(-6,presentTime,H,IERR)
      RETURN
   END IF
   IF ( ((presentTime+0.1_r8*H) == presentTime).OR.(H <= this%Roundoff) ) THEN  ! Step size too small
      Ierr = -7
      CALL ros_ErrorMsg(-7,presentTime,H,IERR)
      RETURN
   END IF

!~~~>  Limit H if necessary to avoid going beyond Tend
   H = MIN(H,ABS(Tend-presentTime))

!~~~>   Compute the function at current time
   CALL FunTemplate( N, presentTime, Y, Fcn0, theKinetics )
   istat(Nfun) = istat(Nfun) + 1

!~~~>  Compute the function derivative with respect to T
   IF (.NOT. this%Autonomous) THEN
      CALL ros_FunTimeDerivative( N, presentTime, this%Roundoff, Y, Fcn0, &
                                  dFdT, istat(Nfun), theKinetics )
   END IF

!~~~>   Compute the Jacobian at current time
   CALL JacTemplate( N, presentTime, Y, Jac0, theKinetics )
   istat(Njac) = istat(Njac) + 1

!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO

   CALL ros_PrepareMatrix( N, H, Direction, this%ros_Gamma(1), Jac0, &
                           Ghimj, Pivot, Singular, istat )
   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       Ierr = -8
       CALL ros_ErrorMsg(-8,presentTime,H,IERR)
       RETURN
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, this%ros_S
         S_ndx = (istage - 1)*(istage - 2)/2
      ! For the 1st istage the function has been computed previously
         IF ( istage == 1 ) THEN
            Fcn(1:N) = Fcn0(1:N)
      ! istage>1 and a new function evaluation is needed at the current istage
         ELSEIF ( this%ros_NewF(istage) ) THEN
           Ynew(1:N) = Y(1:N)
           DO j = 1, istage-1
             Ynew(1:N) = Ynew(1:N) + this%ros_A(S_ndx+j)*K(1:N,j)
           END DO
           Tau = presentTime + this%ros_Alpha(istage)*Direction*H
           CALL FunTemplate( N, Tau, Ynew, Fcn, theKinetics )
           istat(Nfun) = istat(Nfun) + 1
         ENDIF
         K(1:N,istage) = Fcn(1:N)
         DO j = 1, istage-1
           HC = this%ros_C(S_ndx+j)/(Direction*H)
           K(1:N,istage) = K(1:N,istage) + HC*K(1:N,j)
         END DO
         IF ((.NOT. this%Autonomous).AND.(this%ros_Gamma(istage) /= ZERO)) THEN
           HG = Direction*H*this%ros_Gamma(istage)
           K(1:N,istage) = K(1:N,istage) + HG*dFdt(1:N)
         END IF
         CALL ros_Solve( N, Ghimj, Pivot, K(1,istage), istat(Nsol) )
   END DO Stage


!~~~>  Compute the new solution
   Ynew(1:N) = Y(1:N)
   DO j=1,this%ros_S
     Ynew(1:N) = Ynew(1:N) + this%ros_M(j)*K(1:N,j)
   END DO

!~~~>  Compute the error estimation
   Yerr(1:N) = ZERO
   DO j=1,this%ros_S
     Yerr(1:N) = Yerr(1:N) + this%ros_E(j)*K(1:N,j)
   END DO
   Err = ros_ErrorNorm( this, Y, Ynew, Yerr )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(this%FacMax,MAX(this%FacMin,this%FacSafe/Err**(ONE/this%ros_ELO)))
   Hnew = H*Fac

!~~~>  Check the error magnitude and adjust step size
   istat(Nstp) = istat(Nstp) + 1
   istat(Ntotstp) = istat(Ntotstp) + 1
Accepted: &
   IF ( (Err <= ONE).OR.(H <= this%Hmin) ) THEN
      istat(Nacc) = istat(Nacc) + 1
      Y(1:N) = Ynew(1:N)
      presentTime = presentTime + Direction*H
      Hnew = MAX(this%Hmin,MIN(Hnew,this%Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H)
      END IF
      RSTAT(Nhexit) = H
      RSTAT(Nhnew)  = Hnew
      RSTAT(Ntexit) = presentTime
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      H = Hnew
      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
   ELSE Accepted  !~~~> Reject step
      IF (RejectMoreH) THEN
         Hnew = H*this%FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (istat(Nacc) >= 1)  istat(Nrej) = istat(Nrej) + 1
   ENDIF Accepted ! Err <= 1

   END DO UntilAccepted

   END DO TimeLoop

!~~~> Succesful exit
   IERR = 0  !~~~> The integration was successful
   IF( present(istatus) ) THEN
     istatus(:) = istat(:)
   ENDIF
   IF( present(rstatus) ) THEN
     rstatus(:) = rstat(:)
   ENDIF
   IF( present(T) ) THEN
     T = presentTime
   ENDIF

    end subroutine RosenbrockRun

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_FunTimeDerivative ( N, T, Roundoff, Y, Fcn0, &
                                     dFdT, Nfun, theKinetics )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the function by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~> Input arguments
   INTEGER, INTENT(IN)    :: N
   INTEGER, INTENT(INOUT) :: Nfun
   REAL(r8), INTENT(IN)   :: T, Roundoff, Y(N), Fcn0(N)
   TYPE(kinetics_type )   :: theKinetics

!~~~> Output arguments
   REAL(r8), INTENT(OUT) :: dFdT(N)
!~~~> Local variables
   REAL(r8) :: Delta
   REAL(r8) :: factor

   REAL(r8), PARAMETER :: DeltaMin = 1.0E-6_r8

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL FunTemplate( N, T+Delta, Y, dFdT, theKinetics )
   Nfun = Nfun + 1
   factor = ONE/Delta
   dFdT(1:N) = factor*(dFdT(1:N) - Fcn0(1:N))

  END SUBROUTINE ros_FunTimeDerivative

SUBROUTINE FunTemplate( N, T, Y, Ydot, theKinetics )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the ODE function call.
!  Updates the rate coefficients (and possibly the fixed species) at each call
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~> Input variables
   INTEGER, intent(in)  :: N                  ! number of constituents
   REAL(r8), intent(in) :: T                  ! current time
   REAL(r8), intent(in) :: Y(N)               ! constituent concentration
   TYPE(kinetics_type)  :: theKinetics
!~~~> Output variables
   REAL(r8), intent(out) :: Ydot(N)           ! dY(:)/dT

   Ydot(:) =  theKinetics%force( Y )

END SUBROUTINE FunTemplate

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE JacTemplate( N, T, Y, Jcb, theKinetics )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the ODE Jacobian call.
!  Updates the rate coefficients (and possibly the fixed species) at each call
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~> Input variables
   INTEGER, intent(in)  :: N                  ! number of constituents
   REAL(r8), intent(in) :: T                  ! current time
   REAL(r8), intent(in) :: Y(N)               ! constituent concentration
   TYPE(kinetics_type)  :: theKinetics
!~~~> Output variables
   REAL(r8), intent(inout) :: Jcb(N,N)        ! jacobian matrix

   Jcb(:,:) = theKinetics%Jac( Y )

END SUBROUTINE JacTemplate

  SUBROUTINE ros_PrepareMatrix ( N, H, Direction, gam, &
                                 Jac0, Ghimj, Pivot, Singular, istatus )
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*ham) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---

!~~~> Input arguments
   INTEGER, INTENT(IN)  ::  N
   INTEGER, INTENT(IN)  ::  Direction
   INTEGER, INTENT(OUT) ::  Pivot(N)
   INTEGER, INTENT(INOUT) ::  istatus(:)
   REAL(r8), INTENT(IN) ::  Jac0(N,N)
   REAL(r8), INTENT(IN) ::  gam
!~~~> Output arguments
   REAL(r8), INTENT(OUT) :: Ghimj(N,N)
   LOGICAL, INTENT(OUT)  :: Singular
!~~~> Inout arguments
   REAL(r8), INTENT(INOUT) :: H   ! step size is decreased when LU fails

   INTEGER  :: i, ISING, Nconsecutive
   REAL(r8) :: ghinv

   Nconsecutive = 0
   Singular = .TRUE.

   DO WHILE (Singular)
!~~~>    Construct Ghimj = 1/(H*gam) - Jac0
     Ghimj(:,:) = -Jac0(:,:)
     ghinv = ONE/(Direction*H*gam)
     FORALL( I = 1:N ) 
       Ghimj(i,i) = Ghimj(i,i) + ghinv
     ENDFORALL
!~~~>    Compute LU decomposition
     CALL ros_Decomp( N, Ghimj, Pivot, ISING, istatus(Ndec) )
     IF (ISING == 0) THEN
!~~~>    If successful done
        Singular = .FALSE.
     ELSE ! ISING .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
        ISTATUS(Nsng) = ISTATUS(Nsng) + 1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ISING = ',ISING
        IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decompositions
           H = H*HALF
        ELSE  ! More than 5 consecutive failed decompositions
           RETURN
        END IF  ! Nconsecutive
      END IF

   END DO

  END SUBROUTINE ros_PrepareMatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
    SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
      REAL(r8), INTENT(IN) :: T, H
      INTEGER, INTENT(IN)  :: Code
      INTEGER, INTENT(OUT) :: IERR
   
      IERR = Code
      PRINT * , &
       'Forced exit from Rosenbrock due to the following error:' 
     
   SELECT CASE (Code)
    CASE (-1)    
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)    
      PRINT * , '--> Selected Rosenbrock method not implemented'
    CASE (-3)    
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4)    
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5) 
      PRINT * , '--> Improper tolerance values'
    CASE (-6) 
      PRINT * , '--> No of steps exceeds maximum bound'
    CASE (-7) 
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
    CASE (-8)    
      PRINT * , '--> Matrix is repeatedly singular'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT
   
   PRINT *, "T=", T, "and H=", H
     
 END SUBROUTINE ros_ErrorMsg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  FUNCTION ros_ErrorNorm ( this, Y, Ynew, Yerr ) result( Error )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Input arguments
   class(RosenbrockSolver) :: this
   REAL(r8), INTENT(IN) :: Y(:), Ynew(:), Yerr(:)

   REAL(r8) :: Error

! Local variables
   REAL(r8) :: Scale(this%N), Ymax(this%N)

   Ymax(:)  = MAX( ABS(Y(:)),ABS(Ynew(:)) )
   Scale(:) = this%AbsTol(:) + this%RelTol(:)*Ymax(:)
   Error    = MAX( SQRT( sum( (Yerr(:)/Scale(:))**2 )/real(this%N,kind=r8) ),ErrMin )

  END FUNCTION ros_ErrorNorm

  SUBROUTINE ros_Decomp( N, A, Pivot, ISING, Ndec )
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

   CALL  DGEFA( A, N, Pivot, ISING )
   Ndec = Ndec + 1
   IF ( ISING /= 0 ) THEN
      PRINT*,"Error in DGEFA. ISING=",ISING
   END IF  

  END SUBROUTINE ros_Decomp


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Solve( N, A, Pivot, b, Nsol )
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

   CALL  DGESL( A, N, Pivot, b, 0 )

   Nsol = Nsol + 1

  END SUBROUTINE ros_Solve

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros2( this )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    class(RosenbrockSolver) :: this

    real(r8),parameter :: g = 1.0_r8 + 1.0_r8/SQRT(2.0_r8)

    this%rosMethod = RS2
!~~~> Name of the method
    this%ros_Name = 'ROS-2'
!~~~> Number of stages
    this%ros_S = 2

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    this%ros_A(1) = (1.0_r8)/g
    this%ros_C(1) = (-2.0_r8)/g
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    this%ros_NewF(1:2) = .TRUE.
!~~~> M_i = Coefficients for new step solution
    this%ros_M(1)= (3.0_r8)/(2.0_r8*g)
    this%ros_M(2)= (1.0_r8)/(2.0_r8*g)
! E_i = Coefficients for error estimator
    this%ros_E(1) = 1.0_r8/(2.0_r8*g)
    this%ros_E(2) = 1.0_r8/(2.0_r8*g)
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus one
    this%ros_ELO = 2.0_r8
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    this%ros_Alpha(1:2) = (/ 0.0_r8,1.0_r8 /)
!~~~> Gamma_i = \sum_j  gamma_{i,j}
    this%ros_Gamma(1:2) = (/ g,-g /)

 END SUBROUTINE Ros2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros3( this )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   class(RosenbrockSolver) :: this

   this%rosMethod = RS3
!~~~> Name of the method
   this%ros_Name = 'ROS-3'
!~~~> Number of stages
   this%ros_S = 3

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   this%ros_A(1:3)= (/ 1.0_r8, 1.0_r8, 0.0_r8 /)

   this%ros_C(1) = -0.10156171083877702091975600115545E+01_r8
   this%ros_C(2) =  0.40759956452537699824805835358067E+01_r8
   this%ros_C(3) =  0.92076794298330791242156818474003E+01_r8
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   this%ros_NewF(1:3) = (/ .TRUE.,.TRUE.,.FALSE. /)
!~~~> M_i = Coefficients for new step solution
   this%ros_M(1) =  0.1E+01_r8
   this%ros_M(2) =  0.61697947043828245592553615689730E+01_r8
   this%ros_M(3) = -0.42772256543218573326238373806514_r8
! E_i = Coefficients for error estimator
   this%ros_E(1) =  0.5_r8
   this%ros_E(2) = -0.29079558716805469821718236208017E+01_r8
   this%ros_E(3) =  0.22354069897811569627360909276199_r8
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   this%ros_ELO = 3.0_r8
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   this%ros_Alpha(1)= 0.0_r8
   this%ros_Alpha(2)= 0.43586652150845899941601945119356_r8
   this%ros_Alpha(3)= 0.43586652150845899941601945119356_r8
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   this%ros_Gamma(1)= 0.43586652150845899941601945119356_r8
   this%ros_Gamma(2)= 0.24291996454816804366592249683314_r8
   this%ros_Gamma(3)= 0.21851380027664058511513169485832E+01_r8

  END SUBROUTINE Ros3

  SUBROUTINE Ros4( this )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
!     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1990)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   class(RosenbrockSolver) :: this

   this%rosMethod = RS4
!~~~> Name of the method
   this%ros_Name = 'ROS-4'
!~~~> Number of stages
   this%ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   this%ros_A(1) = 0.2000000000000000E+01_r8
   this%ros_A(2) = 0.1867943637803922E+01_r8
   this%ros_A(3) = 0.2344449711399156_r8
   this%ros_A(4) = this%ros_A(2)
   this%ros_A(5) = this%ros_A(3)
   this%ros_A(6) = 0.0_r8

   this%ros_C(1) =-0.7137615036412310E+01_r8
   this%ros_C(2) = 0.2580708087951457E+01_r8
   this%ros_C(3) = 0.6515950076447975_r8
   this%ros_C(4) =-0.2137148994382534E+01_r8
   this%ros_C(5) =-0.3214669691237626_r8
   this%ros_C(6) =-0.6949742501781779_r8
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   this%ros_NewF(1:3)  = .TRUE.
   this%ros_NewF(4)  = .FALSE.
!~~~> M_i = Coefficients for new step solution
   this%ros_M(1) = 0.2255570073418735E+01_r8
   this%ros_M(2) = 0.2870493262186792_r8
   this%ros_M(3) = 0.4353179431840180_r8
   this%ros_M(4) = 0.1093502252409163E+01_r8
!~~~> E_i  = Coefficients for error estimator
   this%ros_E(1) =-0.2815431932141155_r8
   this%ros_E(2) =-0.7276199124938920E-01_r8
   this%ros_E(3) =-0.1082196201495311_r8
   this%ros_E(4) =-0.1093502252409163E+01_r8
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   this%ros_ELO  = 4.0_r8
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   this%ros_Alpha(1) = 0.0_r8
   this%ros_Alpha(2) = 0.1145640000000000E+01_r8
   this%ros_Alpha(3) = 0.6552168638155900_r8
   this%ros_Alpha(4) = this%ros_Alpha(3)
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   this%ros_Gamma(1) = 0.5728200000000000_r8
   this%ros_Gamma(2) =-0.1769193891319233E+01_r8
   this%ros_Gamma(3) = 0.7592633437920482_r8
   this%ros_Gamma(4) =-0.1049021087100450_r8

  END SUBROUTINE Ros4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas3( this )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   class(RosenbrockSolver) :: this

   this%rosMethod = RD3
!~~~> Name of the method
   this%ros_Name = 'RODAS-3'
!~~~> Number of stages
   this%ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   this%ros_A(1:6) = (/ 0.0_r8,2.0_r8,0.0_r8,2.0_r8,0.0_r8,1.0_r8 /)

   this%ros_C(1:5) = (/ 4.0_r8,1.0_r8,-1.0_r8,1.0_r8,-1.0_r8 /)
   this%ros_C(6)   = -(8.0_r8/3.0_r8)

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   this%ros_NewF(1:4) = .TRUE.
   this%ros_NewF(2)   = .FALSE.
!~~~> M_i = Coefficients for new step solution
   this%ros_M(1:4) = (/ 2.0_r8,0.0_r8,1.0_r8,1.0_r8 /)
!~~~> E_i  = Coefficients for error estimator
   this%ros_E(1:3) = 0.0_r8
   this%ros_E(4)   = 1.0_r8
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   this%ros_ELO  = 3.0_r8
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   this%ros_Alpha(1:2) = 0.0_r8
   this%ros_Alpha(3:4) = 1.0_r8
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   this%ros_Gamma(1:2) = (/ 0.5_r8,1.5_r8 /)
   this%ros_Gamma(3:4) = 0.0_r8

  END SUBROUTINE Rodas3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas4( this )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   class(RosenbrockSolver) :: this

    this%rosMethod = RD4
!~~~> Name of the method
    this%ros_Name = 'RODAS-4'
!~~~> Number of stages
    this%ros_S = 6

!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    this%ros_Alpha(1) = 0.000_r8
    this%ros_Alpha(2) = 0.386_r8
    this%ros_Alpha(3) = 0.210_r8
    this%ros_Alpha(4) = 0.630_r8
    this%ros_Alpha(5:6) = 1.000_r8

!~~~> Gamma_i = \sum_j  gamma_{i,j}
    this%ros_Gamma(1) = 0.2500000000000000_r8
    this%ros_Gamma(2) =-0.1043000000000000_r8
    this%ros_Gamma(3) = 0.1035000000000000_r8
    this%ros_Gamma(4) =-0.3620000000000023E-01_r8
    this%ros_Gamma(5:6) = 0.0_r8

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    this%ros_A(1) = 0.1544000000000000E+01_r8
    this%ros_A(2) = 0.9466785280815826_r8
    this%ros_A(3) = 0.2557011698983284_r8
    this%ros_A(4) = 0.3314825187068521E+01_r8
    this%ros_A(5) = 0.2896124015972201E+01_r8
    this%ros_A(6) = 0.9986419139977817_r8
    this%ros_A(7) = 0.1221224509226641E+01_r8
    this%ros_A(8) = 0.6019134481288629E+01_r8
    this%ros_A(9) = 0.1253708332932087E+02_r8
    this%ros_A(10) =-0.6878860361058950_r8
    this%ros_A(11) = this%ros_A(7)
    this%ros_A(12) = this%ros_A(8)
    this%ros_A(13) = this%ros_A(9)
    this%ros_A(14) = this%ros_A(10)
    this%ros_A(15) = 1.0_r8

    this%ros_C(1) =-0.5668800000000000E+01_r8
    this%ros_C(2) =-0.2430093356833875E+01_r8
    this%ros_C(3) =-0.2063599157091915_r8
    this%ros_C(4) =-0.1073529058151375_r8
    this%ros_C(5) =-0.9594562251023355E+01_r8
    this%ros_C(6) =-0.2047028614809616E+02_r8
    this%ros_C(7) = 0.7496443313967647E+01_r8
    this%ros_C(8) =-0.1024680431464352E+02_r8
    this%ros_C(9) =-0.3399990352819905E+02_r8
    this%ros_C(10) = 0.1170890893206160E+02_r8
    this%ros_C(11) = 0.8083246795921522E+01_r8
    this%ros_C(12) =-0.7981132988064893E+01_r8
    this%ros_C(13) =-0.3152159432874371E+02_r8
    this%ros_C(14) = 0.1631930543123136E+02_r8
    this%ros_C(15) =-0.6058818238834054E+01_r8

!~~~> M_i = Coefficients for new step solution
    this%ros_M(1) = this%ros_A(7)
    this%ros_M(2) = this%ros_A(8)
    this%ros_M(3) = this%ros_A(9)
    this%ros_M(4) = this%ros_A(10)
    this%ros_M(5:6) = 1.0_r8

!~~~> E_i  = Coefficients for error estimator
    this%ros_E(1:5) = 0.0_r8
    this%ros_E(6) = 1.0_r8

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    this%ros_NewF(1:6) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    this%ros_ELO = 4.0_r8

  END SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rang3( this )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! STIFFLY-STABLE W METHOD OF ORDER 3, WITH 4 STAGES
!
! J. RANG and L. ANGERMANN
! NEW ROSENBROCK W-METHODS OF ORDER 3
! FOR PARTIAL DIFFERENTIAL ALGEBRAIC
!        EQUATIONS OF INDEX 1                                             
! BIT Numerical Mathematics (2005) 45: 761-787
!  DOI: 10.1007/s10543-005-0035-y
! Table 4.1-4.2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   class(RosenbrockSolver) :: this

    this%rosMethod = RG3
!~~~> Name of the method
    this%ros_Name = 'RANG-3'
!~~~> Number of stages
    this%ros_S = 4

    this%ros_A(1) = 5.09052051067020e+00_r8
    this%ros_A(2) = 5.09052051067020e+00_r8
    this%ros_A(3) = 0.0_r8
    this%ros_A(4) = 4.97628111010787e+00_r8
    this%ros_A(5) = 2.77268164715849e-02_r8
    this%ros_A(6) = 2.29428036027904e-01_r8

    this%ros_C(1) = -1.16790812312283e+01_r8
    this%ros_C(2) = -1.64057326467367e+01_r8
    this%ros_C(3) = -2.77268164715850e-01_r8
    this%ros_C(4) = -8.38103960500476e+00_r8
    this%ros_C(5) = -8.48328409199343e-01_r8
    this%ros_C(6) =  2.87009860433106e-01_r8

    this%ros_M(1) =  5.22582761233094e+00_r8 
    this%ros_M(2) = -5.56971148154165e-01_r8 
    this%ros_M(3) =  3.57979469353645e-01_r8 
    this%ros_M(4) =  1.72337398521064e+00_r8 

    this%ros_E(1) = -5.16845212784040e+00_r8 
    this%ros_E(2) = -1.26351942603842e+00_r8 
    this%ros_E(3) = -1.11022302462516e-16_r8 
    this%ros_E(4) =  2.22044604925031e-16_r8 

    this%ros_Alpha(1) = 0.0_r8
    this%ros_Alpha(2) = 2.21878746765329e+00_r8
    this%ros_Alpha(3) = 2.21878746765329e+0_r8
    this%ros_Alpha(4) = 1.55392337535788e+00_r8

    this%ros_Gamma(1) =  4.35866521508459e-01_r8 
    this%ros_Gamma(2) = -1.78292094614483e+00_r8 
    this%ros_Gamma(3) = -2.46541900496934e+00_r8 
    this%ros_Gamma(4) = -8.05529997906370e-01_r8 


!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    this%ros_NewF(1:4) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    this%ros_ELO = 3.0_r8

  END SUBROUTINE Rang3

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

  type(RosenbrockSolver) :: this

  if( allocated( this%absTol ) ) then
    deallocate( this%absTol )
  endif
  if( allocated( this%relTol ) ) then
    deallocate( this%relTol )
  endif

  end subroutine Finalize

END MODULE Rosenbrock_Solver
