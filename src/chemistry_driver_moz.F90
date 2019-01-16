!-----------------------------------------------------------------------------------------------
! chemistry driver which uses MOZART solver method
!-----------------------------------------------------------------------------------------------
module chemistry_driver_moz

use kinetics_module,  only       : kinetics_type
use kinetics,         only       : kinetics_init, kinetics_run
use machine,          only       : r8 => kind_phys
use const_props_mod,  only       : const_props_type
use prepare_chemistry_mod, only  : prepare_chemistry_init
use Mozart_Solver,    only: MozartSolver

implicit none

type(MozartSolver) :: theSolver
type(kinetics_type), pointer :: theKinetics => null()

contains

!> \section arg_table_chemistry_driver_moz_init Argument Table
!! | local_name | standard_name                                    | long_name                               | units   | rank | type            | kind      | intent | optional |
!! |------------|--------------------------------------------------|-----------------------------------------|---------|------|-----------------|-----------|--------|----------|
!! | TimeStart  | chem_step_start_time                             | Chem step start time                    | s       |    0 | real            | kind_phys | in     | F        |
!! | TimeEnd    | chem_step_end_time                               | Chem step end time                      | s       |    0 | real            | kind_phys | in     | F        |
!! | dt         | time_step_for_physics                            | time_step_for_physics                   | s       |    0 | real            | kind_phys | in     | F        |
!! | errmsg     | ccpp_error_message                               | CCPP error message                      | none    |    0 | character       | len=512   | out    | F        |
!! | errflg     | ccpp_error_flag                                  | CCPP error flag                         | flag    |    0 | integer         |           | out    | F        |
!!
subroutine chemistry_driver_moz_init(TimeStart,TimeEnd, dt, errmsg, errflg)

  implicit none
!-----------------------------------------------------------
!  these dimension parameters will be set by the cafe/configurator
!-----------------------------------------------------------
  real(r8), intent(in)            :: TimeStart, TimeEnd
  real(r8), intent(in)            :: dt
  character(len=512),intent(out)  :: errmsg
  integer, intent(out)            :: errflg          ! error index from CPF

  type(const_props_type), pointer :: cnst_info(:)
  integer            :: nSpecies   ! number prognostic constituents
  integer            :: nkRxt      ! number gas phase reactions
  integer            :: njRxt      ! number of photochemical reactions
  integer            :: nTotRxt    ! total number of chemical reactions

  integer  :: icntrl(20)     ! integer control array for ODE solver
  real(r8) :: rcntrl(20)     ! real control array for ODE solver
  real(r8), allocatable :: absTol(:), relTol(:)
  character(len=40) :: model_name
  
  write(0,*) ' Entered chemistry_driver_moz_init'
  !--- initialize CCPP error handling variables
  errmsg = ''
  errflg = 0

  !   This routine should be only called here when the main program no longer needs to allocate variables
  call prepare_chemistry_init(cnst_info, model_name, nSpecies, nkRxt, njRxt )

  nTotRxt =  nkRxt + njRxt

!-----------------------------------------------------------
!  initialize ode solver "control" variable defaults
!-----------------------------------------------------------
  allocate(absTol(nSpecies))
  allocate(relTol(nSpecies))

  absTol(:) = 1.e-9_r8
  relTol(:) = 1.e-4_r8
  icntrl(:) = 0
  rcntrl(:) = 0._r8

!-----------------------------------------------------------
!  set ode solver "control" variables for MOZART solver
!-----------------------------------------------------------
  icntrl(1) = 1                                 ! autonomous, F depends only on Y
  rcntrl(2) = dt                                ! Hmax
  rcntrl(3) = .01_r8*dt                         ! Hstart

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
!! | local_name | standard_name                                    | long_name                               | units   | rank | type         | kind      | intent | optional |
!! |------------|--------------------------------------------------|-----------------------------------------|---------|------|--------------|-----------|--------|----------|
!! | vmr        | concentration                                    | species concentration                   | mole/mole |    1 |real        | kind_phys | inout  | F        |
!! | TimeStart  | chem_step_start_time                             | Chem step start time                    | s       |    0 | real         | kind_phys | in     | F        |
!! | TimeEnd    | chem_step_end_time                               | Chem step end time                      | s       |    0 | real         | kind_phys | in     | F        |
!! | j_rateConst| photo_rate_constants                             | photochemical rates constants           | s-1     |    1 | real         | kind_phys | in     | F        |
!! | k_rateConst| gasphase_rate_constants                          | k rate constants                        | s-1     |    1 | real         | kind_phys | in     | F        |
!! | c_m        | total_number_density                             | total number density                    | molecules/cm3 | 0 | real      | kind_phys | in     | F        |
!! | errmsg     | ccpp_error_message                               | CCPP error message                      | none    |    0 | character    | len=512   | out    | F        |
!! | errflg     | ccpp_error_flag                                  | CCPP error flag                         | flag    |    0 | integer      |           | out    | F        |
!!
subroutine chemistry_driver_moz_run(vmr, TimeStart, TimeEnd, j_rateConst,  k_rateConst, c_m, errmsg, errflg)


  implicit none
!-----------------------------------------------------------
!  these dimension parameters will be set by the cafe/configurator
!-----------------------------------------------------------

  real(kind=r8), intent(inout)    :: vmr(:)                ! "working" concentration passed thru CPF
  real(r8), intent(in)            :: TimeStart, TimeEnd
  real(kind=r8), intent(in)       :: j_rateConst(:)        ! host model provides photolysis rates for now
  real(kind=r8), intent(in)       :: k_rateConst(:)        ! host model provides photolysis rates for now
  real(kind=r8)                   :: c_m                   ! total number density
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
!! | local_name | standard_name                                    | long_name                               | units   | rank | type         | kind      | intent | optional |
!! |------------|--------------------------------------------------|-----------------------------------------|---------|------|--------------|-----------|--------|----------|
!! | errmsg     | ccpp_error_message                               | CCPP error message                      | none    |    0 | character    | len=512   | out    | F        |
!! | errflg     | ccpp_error_flag                                  | CCPP error flag                         | flag    |    0 | integer      |           | out    | F        |
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

