!-----------------------------------------------------------------------------------------------
! chemistry driver; uses either the Rosenbrock or Mozart solver method
!-----------------------------------------------------------------------------------------------
module chemistry_driver

use kinetics_module,  only       : kinetics_type
use ccpp_kinds,       only       : kind_phys

use ODE_solver, only             : baseOdeSolver
use Rosenbrock_Solver, only      : RosenbrockSolver
use Mozart_Solver, only          : MozartSolver

implicit none

type(RosenbrockSolver), target :: aRosenbrockSolver
type(MozartSolver), target     :: aMozartSolver
class(baseOdeSolver), pointer  :: theSolver => null()
type(kinetics_type), pointer   :: theKinetics => null()

contains

!> \section arg_table_chemistry_driver_init Argument Table
!! \htmlinclude chemistry_driver_init.html
!!
subroutine chemistry_driver_init(nSpecies, nkRxt, njRxt, TimeStart, TimeEnd, dt, errmsg, errflg)

  implicit none
!-----------------------------------------------------------
!  these dimension parameters will be set by the cafe/configurator
!-----------------------------------------------------------
  real(kind_phys), intent(in)     :: TimeStart, TimeEnd
  real(kind_phys), intent(in)     :: dt
  integer, intent(in)             :: nSpecies   ! number prognostic constituents
  integer, intent(in)             :: nkRxt      ! number gas phase reactions
  integer, intent(in)             :: njRxt      ! number of photochemical reactions
  character(len=512),intent(out)  :: errmsg
  integer, intent(out)            :: errflg          ! error index from CPF

  real(kind_phys), parameter :: NOT_SET = -huge(1._kind_phys)

  integer            :: nTotRxt    ! total number of chemical reactions

  integer  :: icntrl(20)     ! integer control array for ODE solver
  real(kind_phys) :: Hstart
  real(kind_phys) :: rcntrl(20)     ! real control array for ODE solver
  real(kind_phys)              :: absTol, relTol
  real(kind_phys), allocatable :: abs_tol(:), rel_tol(:)
  character(len=80) :: model_name
  character(len=80) :: Solver_method = ' '

  namelist /options/ Solver_method
  namelist /errCntrl/ absTol, relTol
  namelist /timeCntrl/ Hstart
  
!-----------------------------------------------------------
!  get and set the Solver method
!-----------------------------------------------------------
  open(unit=10,file='../Solver_options')
  read(unit=10,nml=options)
  read(unit=10,nml=errCntrl)
  read(unit=10,nml=timeCntrl)
  close(unit=10)

  select case( Solver_method )
    case( 'mozart', 'Mozart', 'MOZART' )
      theSolver => aMozartSolver
    case( 'rosenbrock', 'Rosenbrock', 'ROSENBROCK' )
      theSolver => aRosenbrockSolver
    case default
      write(errmsg,*) 'chemistry_driver: solver method must be {Mozart,Rosenbrock}; not ',trim(Solver_method)
      errflg = -1
      return
  end select
  
  !--- initialize CCPP error handling variables
  errmsg = ''
  errflg = 0

  nTotRxt =  nkRxt + njRxt

!-----------------------------------------------------------
!  allocate error controls
!-----------------------------------------------------------
  allocate(abs_tol(nSpecies))
  allocate(rel_tol(nSpecies))

  icntrl(:) = 0
  rcntrl(:) = 0._kind_phys

  select type(theSolver)
    class is (RosenbrockSolver)
!-----------------------------------------------------------
!  set ode solver "control" variables for Rosenbrock solver
!-----------------------------------------------------------
      icntrl(1) = 1                                 ! autonomous, F depends only on Y
      icntrl(3) = 2                                 ! ros3 solver
      rcntrl(2) = dt                                ! Hmax
      rcntrl(3) = .01_kind_phys*dt                  ! Hstart
    class is (MozartSolver)
!-----------------------------------------------------------
!  set ode solver "control" variables for MOZART solver
!-----------------------------------------------------------
      icntrl(1) = 1                                 ! autonomous, F depends only on Y
      rcntrl(2) = dt                                ! Hmax, max time step
      rcntrl(3) = .01_kind_phys*dt                  ! Hstart, initial time step
  end select

  abs_tol(:) = absTol
  rel_tol(:) = relTol

  write(*,*) ' '
  write(*,*) 'icntrl settings'
  write(*,'(10i6)') icntrl(1:10)
  write(*,*) 'rcntrl settings'
  write(*,'(1p,10(1x,g0))') rcntrl(1:10)
  write(*,*) 'Absolute error tolerances'
  write(*,*) abs_tol
  write(*,*) 'Relative error tolerances'
  write(*,*) rel_tol
  write(*,*) ' '

  call theSolver%Initialize( Tstart=TimeStart, Tend=TimeEnd, AbsTol=Abs_tol, RelTol=Rel_tol, &
                             ICNTRL=icntrl, RCNTRL=rcntrl, Ierr=errflg )

!-----------------------------------------------------------
!  allocate & initialize the kinetics
!-----------------------------------------------------------
  allocate( theKinetics )
  call thekinetics%init( nTotRxt, nSpecies )

end subroutine chemistry_driver_init

!> \section arg_table_chemistry_driver_run Argument Table
!! \htmlinclude chemistry_driver_run.html
!!
subroutine chemistry_driver_run(vmr, TimeStart, TimeEnd, j_rateConst,  k_rateConst, number_density_air, errmsg, errflg)

  use kinetics_utilities, only: kinetics_init, kinetics_final

  implicit none
!-----------------------------------------------------------
!  these dimension parameters will be set by the cafe/configurator
!-----------------------------------------------------------
  real(kind=kind_phys), intent(inout)    :: vmr(:)                ! "working" concentration passed thru CPF
  real(kind_phys), intent(in)            :: TimeStart, TimeEnd
  real(kind=kind_phys), intent(in)       :: j_rateConst(:)        ! host model provides photolysis rates for now
  real(kind=kind_phys), intent(in)       :: k_rateConst(:)        ! host model provides photolysis rates for now
  real(kind=kind_phys), intent(in)       :: number_density_air    ! number density
  character(len=512), intent(out)        :: errmsg
  integer, intent(out)                   :: errflg                ! error index from CPF
  real(kind=kind_phys)                   :: number_density(size(vmr))     ! "working" number density of each molecule

  integer :: i, k, n

  !--- initialize CCPP error handling variables
  errmsg = ''
  errflg = 0

  call kinetics_init(vmr, number_density, number_density_air)

!-----------------------------------------------------------
!  update the kinetics
!-----------------------------------------------------------
  call theKinetics%rateConst_update( k_rateConst, j_rateConst, number_density_air )

  call theKinetics%rateConst_print()

!-----------------------------------------------------------
!  solve the current timestep's chemistry
!-----------------------------------------------------------
  call theSolver%Run( Tstart=TimeStart, Tend=TimeEnd, y=number_density, theKinetics=theKinetics, Ierr=errflg )

  call kinetics_final(vmr, number_density, number_density_air)

  if (errflg /= 0) then
    errmsg = 'chemistry_driver_run: ERROR in theSolver%Run'
  end if
  
end subroutine chemistry_driver_run

!> \section arg_table_chemistry_driver_finalize Argument Table
!! \htmlinclude chemistry_driver_finalize.html
!!
subroutine chemistry_driver_finalize(errmsg,errflg)
  character(len=512), intent(out) :: errmsg
  integer, intent(out)            :: errflg                ! error index from CPF
  errmsg = ''
  errflg = 0
  
  deallocate( theKinetics )
  nullify( theKinetics )
  
end subroutine chemistry_driver_finalize

end module chemistry_driver
