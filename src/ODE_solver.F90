
   module ODE_solver

   implicit none

   integer, parameter :: r8 = selected_real_kind( 15 )

   TYPE, abstract :: baseOdeSolver
     CONTAINS
       procedure(OdeSolver_init), deferred :: Initialize
       procedure(OdeSolver_run), deferred  :: Run
   END TYPE baseOdeSolver

  abstract interface
    subroutine OdeSolver_init( this, Tstart, Tend, AbsTol, RelTol, RCNTRL, ICNTRL, IERR)

      use kinetics_module, only : kinetics_type

      import baseOdeSolver
      import r8

      class(baseOdeSolver) :: this
      integer, intent(out) :: Ierr
      real(r8), optional, intent(in) :: Tstart
      real(r8), optional, intent(in) :: Tend
      REAL(r8), optional, INTENT(IN) :: AbsTol(:), RelTol(:)
      INTEGER,  optional, INTENT(IN) :: ICNTRL(:)
      REAL(r8), optional, INTENT(IN) :: RCNTRL(:)

    end subroutine OdeSolver_init

    subroutine OdeSolver_run( this, Y, Tstart, Tend, T, &
                              theKinetics, istatus, rstatus, Ierr )

      use kinetics_module, only : kinetics_type

      import baseOdeSolver
      import r8

      class(baseOdeSolver)  :: this
      integer, intent(out)  :: Ierr
      integer, optional, intent(inout)  :: istatus(:)
      real(r8), optional, intent(inout) :: rstatus(:)
      real(r8),        intent(inout)  :: Y(:)
      real(r8), optional, intent(out) :: T
      real(r8), optional, intent(in)  :: Tstart
      real(r8), optional, intent(in)  :: Tend
      TYPE(kinetics_type), optional   :: theKinetics

    end subroutine OdeSolver_run
  end interface

   end module ODE_solver
