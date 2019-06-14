!-----------------------------------------------------------------------------------------------
! chemistry driver which uses MOZART solver method
!-----------------------------------------------------------------------------------------------
module chemistry_driver_moz

use kinetics_module,  only       : kinetics_type
use kinetics,         only       : kinetics_init, kinetics_run
! use ccpp_kinds,       only       : r8 => kind_phys
use ccpp_kinds,       only       : kind_phys

use const_props_mod,  only       : const_props_type
use Mozart_Solver,    only: MozartSolver

implicit none

type(MozartSolver) :: theSolver
type(kinetics_type), pointer :: theKinetics => null()

contains

!> \section arg_table_chemistry_driver_moz_init Argument Table
!! \htmlinclude chemistry_driver_moz_init.html
!!
subroutine chemistry_driver_moz_init(nSpecies, nkRxt, njRxt, TimeStart,TimeEnd, dt, errmsg, errflg)

  implicit none
!-----------------------------------------------------------
!  these dimension parameters will be set by the cafe/configurator
!-----------------------------------------------------------
  integer, intent(in)                    :: nSpecies   ! number prognostic constituents
  integer, intent(in)                    :: nkRxt      ! number gas phase reactions
  integer, intent(in)                    :: njRxt      ! number of photochemical reactions
  real(kind_phys), intent(in)            :: TimeStart, TimeEnd
  real(kind_phys), intent(in)            :: dt
  character(len=512),intent(out)         :: errmsg
  integer, intent(out)                   :: errflg          ! error index from CPF

  type(const_props_type), allocatable :: cnst_info(:)
  integer            :: nTotRxt    ! total number of chemical reactions

  integer  :: icntrl(20)     ! integer control array for ODE solver
  real(kind_phys) :: rcntrl(20)     ! real control array for ODE solver
  real(kind_phys), allocatable :: absTol(:), relTol(:)
  character(len=80) :: model_name
  
  write(0,*) ' Entered chemistry_driver_moz_init'
  !--- initialize CCPP error handling variables
  errmsg = ''
  errflg = 0

  nTotRxt =  nkRxt + njRxt

!-----------------------------------------------------------
!  initialize ode solver "control" variable defaults
!-----------------------------------------------------------
  allocate(absTol(nSpecies))
  allocate(relTol(nSpecies))

  absTol(:) = 1.e-8_kind_phys
  relTol(:) = 1.e-3_kind_phys
  icntrl(:) = 0
  rcntrl(:) = 0._kind_phys

!-----------------------------------------------------------
!  set ode solver "control" variables for MOZART solver
!-----------------------------------------------------------
  icntrl(1) = 1                                 ! autonomous, F depends only on Y
  rcntrl(2) = dt                                ! Hmax
  rcntrl(3) = .01_kind_phys*dt                         ! Hstart

  write(*,*) ' '
  write(*,*) 'icntrl settings'
  write(*,'(10i6)') icntrl(1:10)
  write(*,*) 'rcntrl settings'
  write(*,'(1p,10(1x,g0))') rcntrl(1:10)
  write(*,*) ' '

  call theSolver%Initialize( Tstart=TimeStart, Tend=TimeEnd, AbsTol=AbsTol, RelTol=RelTol, &
                                       ICNTRL=icntrl, RCNTRL=rcntrl, Ierr=errflg )

!-----------------------------------------------------------
!  initialize the kinetics
!-----------------------------------------------------------
  allocate( theKinetics )
  call kinetics_init(nTotRxt, theKinetics, errmsg, errflg)

end subroutine chemistry_driver_moz_init

!> \section arg_table_chemistry_driver_moz_run Argument Table
!! \htmlinclude chemistry_driver_moz_run.html
!!
subroutine chemistry_driver_moz_run(vmr, TimeStart, TimeEnd, j_rateConst,  k_rateConst, c_m, errmsg, errflg)


  implicit none
!-----------------------------------------------------------
!  these dimension parameters will be set by the cafe/configurator
!-----------------------------------------------------------

  real(kind=kind_phys), intent(inout)    :: vmr(:)                ! "working" concentration passed thru CPF
  real(kind_phys), intent(in)            :: TimeStart, TimeEnd
  real(kind=kind_phys), intent(in)       :: j_rateConst(:)        ! host model provides photolysis rates for now
  real(kind=kind_phys), intent(in)       :: k_rateConst(:)        ! host model provides photolysis rates for now
  real(kind=kind_phys), intent(in)       :: c_m                   ! total number density
  character(len=512), intent(out) :: errmsg
  integer, intent(out)            :: errflg                ! error index from CPF

  integer            :: i, k, n

  !--- initialize CCPP error handling variables
  errmsg = ''
  errflg = 0

!-----------------------------------------------------------
!  update the kinetics
!-----------------------------------------------------------
  call kinetics_run(theKinetics, k_rateConst, j_rateConst, c_m, errmsg, errflg)
  if (errflg /= 0) return

  call theKinetics%rateConst_print()

!-----------------------------------------------------------
!  solve the current timestep's chemistry
!-----------------------------------------------------------
  call theSolver%Run( Tstart=TimeStart, Tend=TimeEnd, y=vmr, theKinetics=theKinetics, Ierr=errflg  )

  if (errflg /= 0) then
     errmsg = 'ERROR: theSolver%Run'
     return
  end if
  
end subroutine chemistry_driver_moz_run

!> \section arg_table_chemistry_driver_moz_finalize Argument Table
!! \htmlinclude chemistry_driver_moz_finalize.html
!!
subroutine chemistry_driver_moz_finalize(errmsg,errflg)
  character(len=512), intent(out) :: errmsg
  integer, intent(out)            :: errflg                ! error index from CPF
  errmsg = ''
  errflg = 0
  
  deallocate( theKinetics )
  nullify( theKinetics )
  
end subroutine chemistry_driver_moz_finalize

end module chemistry_driver_moz

