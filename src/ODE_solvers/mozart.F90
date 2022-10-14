! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \todo Consider rewriting micm_ODE_solver_mozart to follow musica style conventions

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Backward Euler using N-R iteration
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

MODULE micm_ODE_solver_mozart

  USE micm_kinetics,                   only : kinetics_t
  USE micm_ODE_solver,                 only : ODE_solver_t
  USE musica_constants,                only : r8 => musica_dk, musica_ik

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ODE_solver_mozart_t

!~~~>  Statistics on the work performed by the Mozart method
  INTEGER, PARAMETER :: Nfun=1, Njac=2, Nstp=3, Nacc=4, &
                        Nrej=5, Ndec=6, Nsol=7, Nsng=8, &
                        Ntotstp=9, Nred=10, &
                        Ntexit=1, Nhexit=2, Nhnew = 3

  REAL(r8), PARAMETER :: ZERO = 0.0_r8, HALF = .5_r8, ONE = 1.0_r8, TEN = 10._r8
  REAL(r8), PARAMETER :: FIVEPCNT = .05_r8
  REAL(r8), PARAMETER :: ErrMin   = 1.e-10_r8
  REAL(r8), PARAMETER :: DeltaMin = 1.0E-5_r8

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Mozart extension
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TYPE, extends(ODE_solver_t) :: ODE_solver_mozart_t
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
      procedure :: solve
      procedure :: preprocess_input
      final     :: Finalize
  END TYPE ODE_solver_mozart_t

  !> Mozart solver constructor
  interface ODE_solver_mozart_t
    module procedure :: constructor
  end interface ODE_solver_mozart_t

CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function constructor( config ) result( this )

    use musica_config,                 only : config_t

    !> New solver
    class(ODE_solver_mozart_t), pointer  :: this
    !> Solver configuration data
    type(config_t), intent(inout) :: config

    character(len=*), parameter :: my_name = 'Mozart ODE solver constructor'
    real(r8) :: Tstart, Tend, time_step__s, abs_tol
    integer(musica_ik) :: rel_tol
    integer :: i, N, ierr
    integer(musica_ik), pointer :: icntrl(:)
    real(r8), pointer :: rcntrl(:)

    allocate( this )

    call config%get( "chemistry time step", "s", time_step__s, my_name )
    call config%get( "number of variables", N, my_name )
    call config%get( "absolute tolerance", abs_tol, my_name )
    call config%get( "relative tolerance", rel_tol, my_name )

    icntrl => this%icntrl
    rcntrl => this%rcntrl

    Tstart = 0.0_r8
    Tend   = time_step__s

    icntrl(:) = 0
    rcntrl(:) = 0.0_r8

    icntrl(1) = 1                                 ! autonomous, F depends only on Y
    rcntrl(2) = time_step__s                      ! Hmax, max time step
    rcntrl(3) = .01_r8*time_step__s               ! Hstart, initial time step

    this%N = N

    allocate( this%AbsTol(N),this%RelTol(N) )

!~~~>  Autonomous or time dependent ODE. Default is time dependent.
    this%Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
    IF (ICNTRL(2) == 0) THEN
      this%VectorTol = .TRUE.
      this%AbsTol(:N) = abs_tol
      this%RelTol(:N) = rel_tol
    ELSE
      this%VectorTol = .FALSE.
      this%AbsTol(:N) = abs_tol  !> \todo support vector tolerances in input data
      this%RelTol(:N) = rel_tol
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
      IF( (this%AbsTol(i) <= ZERO) .OR. (this%RelTol(i) <= TEN*this%Roundoff) &
                                   .OR. (this%RelTol(i) >= ONE) ) THEN
        PRINT * , ' AbsTol(',i,') = ',this%AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',this%RelTol(i)
        IERR = -5
        CALL moz_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO

    !> \todo develop log process for MUSICA, and update MOZART constructor
    if (.false.) then
       write(*,*) ' '
       write(*,*) 'Mozart ODE solver package initialized'
       write(*,*) 'Autonomous = ',this%Autonomous
       write(*,*) 'Vectol     = ',this%VectorTol
       write(*,*) 'N,Max_NR{iter,rej} = ',this%N,this%iterMax,this%rejMax
       write(*,*) 'Hmin,Hmax,Hstart     = ',this%Hmin,this%Hmax,this%Hstart
       write(*,*) 'Fac{Rej,Acc} = ',this%FacRej,this%FacAcc
       write(*,*) 'RelTol       = ',this%RelTol(:)
       write(*,*) 'AbsTol       = ',this%AbsTol(:)
       write(*,*) ' '
    endif

    end function constructor

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine solve( this, Y, Tstart, Tend, T, theKinetics, Ierr )

      class(ODE_solver_mozart_t), intent(inout) :: this
      integer,                    intent(out)   :: Ierr
      real(r8),                   intent(inout) :: Y(:,:)
      real(r8), optional,         intent(out)   :: T
      real(r8), optional,         intent(in)    :: Tstart
      real(r8), optional,         intent(in)    :: Tend
      TYPE(kinetics_t), optional, intent(inout) :: theKinetics

! ~~~~ Local variables
      INTEGER  :: N
      INTEGER  :: nIter
      INTEGER  :: j
      INTEGER  :: rejCnt
      REAL(r8) :: H, Hinv, Hnew
      REAL(r8) :: presentTime
      REAL(r8) :: truncError
      REAL(r8), target  :: residual(this%N)
      REAL(r8), pointer :: deltaY(:)
      REAL(r8) :: Ynew(this%N), Yerr(this%N)
      REAL(r8) :: Fcn(this%N)
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
   this%rcntrl(Nhexit) = ZERO
   H = MIN( MAX(ABS(this%Hmin),ABS(this%Hstart)) , ABS(this%Hmax) )
   IF (ABS(H) <= TEN*this%Roundoff) H = DeltaMin

   Ierr   = 0
   rejCnt = 0
   rejectLastH = .FALSE.
   this%icntrl(:) = 0
   this%rcntrl(:) = ZERO

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
     Ynew(1:N) = Y(1,1:N)
!~~~>  Newton-Raphson iteration loop
NRloop: DO nIter = 1,this%iterMax
!~~~>   Compute the Jacobian
       this%icntrl(Njac) = this%icntrl(Njac) + 1
       CALL theKinetics%LinFactor( H, ONE, Ynew, Singular, this%icntrl )
       IF (Singular) THEN ! More than 5 consecutive failed decompositions
         Ierr = -8
         CALL moz_ErrorMsg(-8,presentTime,H,IERR)
         RETURN
       END IF
!~~~>   Compute the function at current time
       Fcn(:) = theKinetics%force( Ynew )
       this%icntrl(Nfun) = this%icntrl(Nfun) + 1
       residual(1:N) = Fcn(1:N) - (Ynew(1:N) - Y(1,1:N))*Hinv
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
       truncError = moz_TruncErrorNorm ( this, Y(1,:), Ynew, Yerr )
       Hnew = SQRT( 2.0_r8/truncError )
       Hnew = MAX(this%Hmin,MIN(Hnew,this%Hmax))
!~~~>   Check step truncation error
acceptStep: &
!      IF( truncError <= 2.0_r8/(H*H) ) THEN
       IF( .5_r8*H*H*truncError <= (ONE + 100._r8*this%Roundoff) ) THEN
         !write(*,*) 'mozartRun: step accepted'
         Y(1,1:N) = Ynew(1:N)
         presentTime = presentTime + H
         IF (RejectLastH) THEN  ! No step size increase after a rejected step
           Hnew = MIN(Hnew,H)
         END IF
         this%rcntrl(Nhexit) = H
         this%rcntrl(Nhnew)  = Hnew
         this%rcntrl(Ntexit) = presentTime
         rejCnt = 0
         RejectLastH    = .FALSE.
         this%icntrl(Nstp)    = this%icntrl(Nstp) + 1
         this%icntrl(Ntotstp) = this%icntrl(Ntotstp) + 1
         this%icntrl(Nacc)    = this%icntrl(Nacc) + 1
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
         IF (this%icntrl(Nacc) >= 1)  this%icntrl(Nrej) = this%icntrl(Nrej) + 1
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
         this%icntrl(Nred) = this%icntrl(Nred) + 1
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

!   stop 'Debugging'

    end subroutine solve

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine preprocess_input(this,config,output_path)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Preprocesses solver input data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      use musica_assert,               only : assert
      use musica_config,               only : config_t

      class(ODE_solver_mozart_t), intent(in) :: this
      type(config_t), intent(out) :: config
      character(len=*), intent(in) :: output_path

      character(len=*), parameter :: my_name = "Mozart solver preprocessor"

      !! \todo Update Mozart preprocessor once vector tolerances are supported
      call assert( 905188514, allocated( this%AbsTol ) )
      call assert( 959668300, allocated( this%RelTol ) )
      call assert( 619354492, size( this%AbsTol ) .ge. 1 )
      call assert( 561515933, size( this%RelTol ) .ge. 1 )

      call config%add( "chemistry time step", "s", this%rcntrl(2), my_name )

      call config%add( "type",               "Mozart",       my_name )
      call config%add( "absolute tolerance", this%AbsTol(1), my_name )
      call config%add( "relative tolerance", this%RelTol(1), my_name )

    end subroutine preprocess_input

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE moz_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      use musica_assert,               only : die_msg

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

   call die_msg( 252212991, "Mozart solver error" )
 END SUBROUTINE moz_ErrorMsg

  FUNCTION moz_NRErrorNorm ( this, deltaY, Ynew, Yerr ) result( Converged )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Input arguments
   class(ODE_solver_mozart_t)  :: this
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
   class(ODE_solver_mozart_t)  :: this
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

  type(ODE_solver_mozart_t) :: this

  if( allocated( this%absTol ) ) then
    deallocate( this%absTol )
  endif
  if( allocated( this%relTol ) ) then
    deallocate( this%relTol )
  endif

  end subroutine Finalize

END MODULE micm_ODE_solver_mozart
