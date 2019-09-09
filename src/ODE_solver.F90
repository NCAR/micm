
   module ODE_solver

   use ccpp_kinds, only: r8 => kind_phys
   use ccpp_kinds, only:  kind_phys
   
   implicit none

   TYPE, abstract :: baseOdeSolver
     CONTAINS
       procedure(OdeSolver_init), deferred :: Initialize
       procedure(OdeSolver_run), deferred  :: Run
       logical :: print_log_message = .false.
   END TYPE baseOdeSolver

  abstract interface
    subroutine OdeSolver_init( this, Tstart, Tend, AbsTol, RelTol, RCNTRL, ICNTRL, IERR)

      use kinetics_module, only : kinetics_type

      import baseOdeSolver
      import r8
      import kind_phys

      class(baseOdeSolver) :: this
      integer, intent(out) :: Ierr
      real(kind_phys), optional, intent(in) :: Tstart
      real(kind_phys), optional, intent(in) :: Tend
      REAL(kind_phys), optional, INTENT(IN) :: AbsTol(:), RelTol(:)
      INTEGER,  optional, INTENT(IN) :: ICNTRL(:)
      REAL(kind_phys), optional, INTENT(IN) :: RCNTRL(:)

    end subroutine OdeSolver_init

    subroutine OdeSolver_run( this, Y, Tstart, Tend, T, &
                              theKinetics, istatus, rstatus, Ierr )

      use kinetics_module, only : kinetics_type

      import baseOdeSolver
      import r8
      import kind_phys

      class(baseOdeSolver)  :: this
      integer, intent(out)  :: Ierr
      integer, optional, intent(inout)  :: istatus(:)
      real(kind_phys), optional, intent(inout) :: rstatus(:)
      real(kind_phys),        intent(inout)  :: Y(:)
      real(kind_phys), optional, intent(out) :: T
      real(kind_phys), optional, intent(in)  :: Tstart
      real(kind_phys), optional, intent(in)  :: Tend
      TYPE(kinetics_type), optional   :: theKinetics

    end subroutine OdeSolver_run
  end interface

   end module ODE_solver
