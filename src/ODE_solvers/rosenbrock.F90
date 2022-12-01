! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! TODO consider rewriting micm_ODE_solver_rosenbrock to follow musica style conventions

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Rosenbrock - Implementation of several Rosenbrock methods:             !
!               * Ros2                                                    !
!               * Ros3                                                    !
!               * Ros4                                                    !
!               * Rodas3                                                  !
!               * Rodas4                                                  !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

MODULE micm_ODE_solver_rosenbrock

  use micm_ODE_solver,                 only : ODE_solver_t
  use micm_kinetics,                   only : kinetics_t
  use musica_constants,                only : r8=>musica_dk, musica_ik
  use constants,                       only : ncell=>kNumberOfGridCells, &
                                              VLEN, STREAM0

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ODE_solver_rosenbrock_t

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
  TYPE, extends(ODE_solver_t) :: ODE_solver_rosenbrock_t
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
      procedure :: solve
      procedure :: preprocess_input
      final     :: Finalize
  END TYPE ODE_solver_rosenbrock_t

  !> Rosenbrock solver constructor
  interface ODE_solver_rosenbrock_t
    module procedure :: constructor
  end interface ODE_solver_rosenbrock_t

CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function constructor( config ) result( this )

    use musica_config,                 only : config_t

    !> New solver
    class(ODE_solver_rosenbrock_t), pointer :: this
    !> Solver configuration data
    type(config_t), intent(inout) :: config

!   local variables
    character(len=*), parameter :: my_name = 'Rosenbrock ODE solver constructor'
    integer(musica_ik), pointer :: icntrl(:)
    real(r8), pointer :: rcntrl(:)
    real(r8) :: Tstart, Tend, time_step__s, abs_tol, rel_tol
    integer :: i, N, ierr

    allocate( this )

    call config%get( "chemistry time step", "s", time_step__s, my_name )
    call config%get( "number of variables", N, my_name )
    call config%get( "absolute tolerance", abs_tol, my_name )
    call config%get( "relative tolerance", rel_tol, my_name )

    Tstart = 0.0_r8
    Tend   = time_step__s

    icntrl => this%icntrl
    rcntrl => this%rcntrl

    icntrl(:) = 0
    rcntrl(:) = 0.0_r8

    icntrl(1) = 1                                 ! autonomous, F depends only on Y
    icntrl(3) = 2                                 ! ros3 solver
    rcntrl(2) = time_step__s                      ! Hmax
    rcntrl(3) = .01_r8*time_step__s               ! Hstart

    this%N = N

    allocate( this%AbsTol(N),this%RelTol(N) )

!~~~>  Autonomous or time dependent ODE. Default is time dependent.
    this%Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
    IF (ICNTRL(2) == 0) THEN
      this%VectorTol = .TRUE.
      this%UplimTol  = N
      this%AbsTol(:N) = abs_tol
      this%RelTol(:N) = rel_tol
    ELSE
      this%VectorTol = .FALSE.
      this%UplimTol  = 1
      this%AbsTol(:N) = abs_tol !> \todo support vector tolerances in input data
      this%RelTol(:N) = rel_tol
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
      IF( (this%AbsTol(i) <= ZERO) .OR. (this%RelTol(i) <= TEN*this%Roundoff) &
                                   .OR. (this%RelTol(i) >= ONE) ) THEN
        PRINT * , ' AbsTol(',i,') = ',this%AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',this%RelTol(i)
        IERR = -5
        CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO

    !> \todo develop logging strategy for MUSICA and update ROSENBROCK constructor
    if (.false.) then
       write(*,*) ' '
       write(*,*) 'Rosenbrock ODE solver package initialized'
       write(*,*) 'Autonomous = ',this%Autonomous
       write(*,*) 'Vectol     = ',this%VectorTol
       write(*,*) 'N,ros_S,Max_no_steps = ',this%N,this%ros_S,this%Max_no_steps
       write(*,*) 'Hmin,Hmax,Hstart     = ',this%Hmin,this%Hmax,this%Hstart
       write(*,*) 'Fac{Min,Max,Rej,Safe} = ',this%FacMin,this%FacMax,this%FacRej,this%FacSafe
       write(*,*) 'RelTol                = ',this%RelTol(:)
       write(*,*) 'AbsTol                = ',this%AbsTol(:)
       write(*,*) ' '
    endif

    !$acc enter data copyin(this,this%ros_A,this%ros_M,this%ros_E, &
    !$acc                   this%AbsTol,this%RelTol,this%ros_Gamma) &
    !$acc            async(STREAM0)

    end function constructor

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine solve( this, Y, Tstart, Tend, T, theKinetics, Ierr )

      class(ODE_solver_rosenbrock_t), intent(inout) :: this
      integer,                        intent(out)   :: Ierr
      real(r8),                       intent(inout) :: Y(ncell,this%N) ! (grid cell, species)
      real(r8), optional,             intent(out)   :: T
      real(r8), optional,             intent(in)    :: Tstart
      real(r8), optional,             intent(in)    :: Tend
      TYPE(kinetics_t), optional,     intent(inout) :: theKinetics

! ~~~~ Local variables
      INTEGER  :: N
      INTEGER  :: S_ndx
      INTEGER  :: i, j, m, istage
      integer  :: istat(20)
      REAL(r8) :: H, Hnew, HC, Fac, Tau, Err
      REAL(r8) :: presentTime
      REAL(r8) :: Ynew(ncell,this%N)
      REAL(r8) :: Fcn0(ncell,this%N), Fcn(ncell,this%N)
      REAL(r8) :: K(ncell,this%N,this%ros_S)
      REAL(r8) :: Yerr(ncell,this%N)
      real(r8) :: rstat(20)
      LOGICAL  :: RejectLastH, RejectMoreH, Singular

   !$acc enter data create(Ynew,Fcn0,Fcn,K,Yerr) async(STREAM0)
   !$acc update device(Y,theKinetics%rateConst,theKinetics%environment) async(STREAM0)

!~~~>  Initial preparations
   IF( .not. present(theKinetics) .or. .not. present(Tstart) .or. &
       .not. present(Tend) ) THEN
     Ierr = -10
     RETURN
   ENDIF

!~~~>  Initializations, check incoming time step
   N = this%n
   presentTime = Tstart
   this%rcntrl(Nhexit) = ZERO
   H = MIN( MAX(ABS(this%Hmin),ABS(this%Hstart)) , ABS(this%Hmax) )
   IF (ABS(H) <= TEN*this%Roundoff) H = DeltaMin

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.

   this%icntrl(:) = 0
   this%rcntrl(:) = ZERO

!~~~> Time loop
TimeLoop: DO WHILE ( (presentTime-Tend)+this%Roundoff <= ZERO )

   IF ( this%icntrl(Nstp) > this%Max_no_steps ) THEN  ! Too many steps
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
   call theKinetics%calc_force( Y, Fcn0 )
   this%icntrl(Nfun) = this%icntrl(Nfun) + 1

!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO

!~~~>  Form and factor the rosenbrock ode jacobian
   CALL theKinetics%LinFactor( H, this%ros_Gamma(1), Y, Singular, this%icntrl )
   this%icntrl(Njac) = this%icntrl(Njac) + 1
   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       Ierr = -8
       CALL ros_ErrorMsg(-8,presentTime,H,IERR)
       RETURN
   END IF

!~~~>   Compute the stages
Stage_loop: &
   DO istage = 1, this%ros_S
     IF ( istage /= 1 ) THEN
       S_ndx = (istage - 1)*(istage - 2)/2
       IF ( this%ros_NewF(istage) ) THEN
         !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
         !$acc loop gang vector collapse(2)
         do m = 1, N
            do i = 1, ncell
               Ynew(i,m) = Y(i,m)
               DO j = 1, istage-1
                  Ynew(i,m) = Ynew(i,m) + this%ros_A(S_ndx+j)*K(i,m,j)
               END DO
            end do
         end do
         !$acc end parallel
         Tau = presentTime + this%ros_Alpha(istage)*H
         call theKinetics%calc_force( Ynew, Fcn )
         this%icntrl(Nfun) = this%icntrl(Nfun) + 1
       ENDIF
       !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
       !$acc loop gang vector collapse(2)
       do m = 1, N
          do i = 1, ncell
             K(i,m,istage) = Fcn(i,m)
             DO j = 1, istage-1
                HC = this%ros_C(S_ndx+j)/H
                K(i,m,istage) = K(i,m,istage) + HC*K(i,m,j)
             END DO
          end do
       end do
       !$acc end parallel
     ELSE
       !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
       !$acc loop gang vector collapse(2)
       do m = 1, N
          do i = 1, ncell
             K(i,m,1) = Fcn0(i,m)
             Fcn(i,m) = Fcn0(i,m)
          end do
       end do
       !$acc end parallel
     ENDIF
     CALL theKinetics%LinSolve( K(1:ncell,1:N,istage) )
     this%icntrl(Nsol) = this%icntrl(Nsol) + 1
   END DO Stage_loop

!~~~>  Compute the new solution & the error estimation

   !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
   !$acc loop gang vector collapse(2)
   do m = 1, N
      do i = 1, ncell
         Ynew(i,m) = Y(i,m)
         Yerr(i,m) = ZERO
         DO j=1,this%ros_S
            Ynew(i,m) = Ynew(i,m) + this%ros_M(j)*K(i,m,j)
            Yerr(i,m) = Yerr(i,m) + this%ros_E(j)*K(i,m,j)
         END DO
      end do
   end do
   !$acc end parallel

   Err = ros_ErrorNorm( this, Y, Ynew, Yerr )

   !$acc wait (STREAM0)

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(this%FacMax,MAX(this%FacMin,this%FacSafe/Err**(ONE/this%ros_ELO)))
   Hnew = H*Fac

!~~~>  Check the error magnitude and adjust step size
   this%icntrl(Nstp) = this%icntrl(Nstp) + 1
   this%icntrl(Ntotstp) = this%icntrl(Ntotstp) + 1
Accepted: &
   IF ( (Err <= ONE).OR.(H <= this%Hmin) ) THEN
      this%icntrl(Nacc) = this%icntrl(Nacc) + 1
      !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
      !$acc loop gang vector collapse(2)
      do m = 1, N 
         do i = 1, ncell
            Y(i,m) = Ynew(i,m)
         end do
      end do
      !$acc end parallel
      presentTime = presentTime + H
      Hnew = MAX(this%Hmin,MIN(Hnew,this%Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H)
      END IF
      this%rcntrl(Nhexit) = H
      this%rcntrl(Nhnew)  = Hnew
      this%rcntrl(Ntexit) = presentTime
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
      IF (this%icntrl(Nacc) >= 1)  this%icntrl(Nrej) = this%icntrl(Nrej) + 1
   ENDIF Accepted ! Err <= 1

   END DO UntilAccepted

   END DO TimeLoop

   !$acc update self(Y) wait(STREAM0)
   !$acc exit data delete(Ynew,Fcn0,Fcn,K,Yerr) async(STREAM0)
 
!~~~> Succesful exit
   IERR = 0  !~~~> The integration was successful
   IF( present(T) ) THEN
     T = presentTime
   ENDIF

    end subroutine solve
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine preprocess_input(this,config,output_path)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Preprocesses solver input data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      use musica_assert,               only : assert
      use musica_config,               only : config_t

      class(ODE_solver_rosenbrock_t), intent(in) :: this
      type(config_t), intent(out) :: config
      character(len=*), intent(in) :: output_path

      character(len=*), parameter :: my_name = "Rosenbrock solver preprocessor"

      !! \todo Update Rosenbrock preprocessor once vector tolerances are supported
      call assert( 194331941, allocated( this%AbsTol ) )
      call assert( 924175036, allocated( this%RelTol ) )
      call assert( 418968631, size( this%AbsTol ) .ge. 1 )
      call assert( 248811727, size( this%RelTol ) .ge. 1 )

      call config%add( "chemistry time step", "s", this%rcntrl(2), my_name )

      call config%add( "type",               "Rosenbrock",   my_name )
      call config%add( "absolute tolerance", this%AbsTol(1), my_name )
      call config%add( "relative tolerance", this%RelTol(1), my_name )

    end subroutine preprocess_input
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      use musica_assert,               only : die_msg

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

   call die_msg( 752607384, "Rosenbrock solver failure" )

 END SUBROUTINE ros_ErrorMsg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  FUNCTION ros_ErrorNorm ( this, Y, Ynew, Yerr ) result( Error )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Input arguments
   class(ODE_solver_rosenbrock_t) :: this
   REAL(r8), INTENT(IN) :: Y(ncell,this%N), Ynew(ncell,this%N), &
                           Yerr(ncell,this%N)

   REAL(r8) :: Error

! Local variables
   REAL(r8) :: Scale, Ymax, sum_tmp
   integer  :: i, m

   Error = 0._r8

   !$acc parallel default(present) vector_length(VLEN) async(STREAM0)
   !$acc loop gang reduction(max:Error)
   do i = 1, ncell
      sum_tmp = 0._r8
      !$acc loop vector reduction(+:sum_tmp)
      do m = 1, this%N
         Ymax = MAX( ABS(Y(i,m)),ABS(Ynew(i,m)) )
         Scale  = this%AbsTol(m) + this%RelTol(m)*Ymax
         sum_tmp = sum_tmp + (Yerr(i,m)/Scale)**2
      end do
      Error     = MAX( SQRT( sum_tmp / real(this%N,kind=r8) ), ErrMin, Error )
   end do
   !$acc end parallel

  END FUNCTION ros_ErrorNorm

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros2( this )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    class(ODE_solver_rosenbrock_t) :: this

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
   class(ODE_solver_rosenbrock_t) :: this

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
   class(ODE_solver_rosenbrock_t) :: this

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
   class(ODE_solver_rosenbrock_t) :: this

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
   class(ODE_solver_rosenbrock_t) :: this

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
   class(ODE_solver_rosenbrock_t) :: this

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

  subroutine Finalize( this )

  type(ODE_solver_rosenbrock_t) :: this

  if( allocated( this%absTol ) ) then
    deallocate( this%absTol )
  endif
  if( allocated( this%relTol ) ) then
    deallocate( this%relTol )
  endif

  end subroutine Finalize

END MODULE micm_ODE_solver_rosenbrock
