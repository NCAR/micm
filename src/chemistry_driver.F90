module chemistry_driver

use solver_var_defs,  only: Solver_type
use kinetics_module,  only: kinetics_type
use kinetics,         only: kinetics_init, kinetics_run
use chem_solve,       only: chem_solve_init, chem_solve_run
use k_rateConst,      only: k_rateConst_init, k_rateConst_run
use machine,          only: r8 => kind_phys
use const_props_mod,  only: const_props_type
use prepare_chemistry_mod, only: prepare_chemistry_init


implicit none

contains


!> \section arg_table_chemistry_driver_init Argument Table
!! | local_name | standard_name                                    | long_name                               | units   | rank | type            | kind      | intent | optional |
!! |------------|--------------------------------------------------|-----------------------------------------|---------|------|-----------------|-----------|--------|----------|
!! | TimeStart  | chem_step_start_time                             | Chem step start time                    | s       |    0 | real            | kind_phys | in     | F        |
!! | TimeEnd    | chem_step_end_time                               | Chem step end time                      | s       |    0 | real            | kind_phys | in     | F        |
!! | nTotRxt    | Number_chemical_reactions                        |                                         | none    |    0 | integer         |           | out    | F        |
!! | theKinetics | kinetics_data                                   | chemistry kinetics                      | DDT     |    0 | kinetics_type   |           | out    | F        |
!! | ODE_obj    | ODE_ddt                                          | ODE derived data type                   | DDT     |    0 | Solver_type     |           | none   | F        |
!! | k_rateConst| gasphase_rate_constants                          | k rate constants                        | s-1     |    1 | real            | kind_phys | inout  | F        |
!! | AbsTol     | abs_trunc_error                                  | ODE absolute step truncation error      | none    |    1 | real            | kind_phys | in     | F        |
!! | RelTol     | rel_trunc_error                                  | ODE relative step truncation error      | none    |    1 | real            | kind_phys | in     | F        |
!! | icntrl     | ODE_icontrol                                     | ODE integer controls                    | flag    |    1 | integer         |           | in     | F        |
!! | rcntrl     | ODE_rcontrol                                     | ODE real controls                       | none    |    1 | real            | kind_phys | in     | F        |
!! | errmsg     | ccpp_error_message                               | CCPP error message                      | none    |    0 | character       | len=512   | out    | F        |
!! | cnst_info  | chemistry_constituent_info                       | chemistry_constituent_info              | DDT     |    1 | const_props_type|           | out    | F        |
!! | errflg     | ccpp_error_flag                                  | CCPP error flag                         | flag    |    0 | integer         |           | out    | F        |
!!
subroutine chemistry_driver_init(TimeStart,TimeEnd, nTotRxt,  theKinetics, ODE_obj, k_rateConst, AbsTol, RelTol, icntrl, rcntrl, cnst_info, errmsg, errflg)

  implicit none
!-----------------------------------------------------------
!  these dimension parameters will be set by the cafe/configurator
!-----------------------------------------------------------
  integer, intent(out) :: nTotRxt    ! total number of chemical reactions

  integer            :: nSpecies   ! number prognostic constituents
  integer            :: nkRxt      ! number gas phase reactions
  integer            :: njRxt      ! number of photochemical reactions
  integer            :: i, k, n
  integer            :: errflg          ! error index from CPF
  integer            :: ierr
  real(kind=r8),pointer      :: k_rateConst(:)  ! host model provides photolysis rates for now
  character(len=512) :: errmsg

  integer  :: icntrl(20)     ! integer control array for ODE solver
  real(r8) :: rcntrl(20)     ! real control array for ODE solver
  real(r8) :: TimeStart, TimeEnd, Time, dt
  real(r8), intent(in) :: AbsTol(:), RelTol(:)

! declare the types
  type(Solver_type),    pointer  :: ODE_obj
  type(kinetics_type),  pointer  :: theKinetics
  type(const_props_type), pointer :: cnst_info(:)

  character(len=16) :: cnst_name

  write(0,*) ' Entered chemistry_driver_init'

!   This routine should be called here when the main program no longer needs to allocate variables
!   call prepare_chemistry_init(cnst_info, nSpecies, nkRxt, njRxt, nTotRxt)

  call k_rateConst_init(k_rateConst, errflg, errmsg)
  write(0,*) 'inside init, k_rateConst(1)=',k_rateConst(1)
  call kinetics_init(nTotRxt, theKinetics, errmsg, errflg)
  call chem_solve_init(TimeStart, TimeEnd, AbsTol, RelTol, icntrl, rcntrl, ODE_obj, errmsg, errflg)

end subroutine chemistry_driver_init

!> \section arg_table_chemistry_driver_run Argument Table
!! | local_name | standard_name                                    | long_name                               | units   | rank | type         | kind      | intent | optional |
!! |------------|--------------------------------------------------|-----------------------------------------|---------|------|--------------|-----------|--------|----------|
!! | vmr        | concentration                                    | species concentration                   | mole/mole |    1 |real        | kind_phys | inout  | F        |
!! | TimeStart  | chem_step_start_time                             | Chem step start time                    | s       |    0 | real         | kind_phys | in     | F        |
!! | TimeEnd    | chem_step_end_time                               | Chem step end time                      | s       |    0 | real         | kind_phys | in     | F        |
!! | Time       | Simulation_time                                  | Present simulation time                 | s       |    0 | real         | kind_phys | in     | F        |
!! | theKinetics | kinetics_data                                   | chemistry kinetics                      | DDT     |    0 | kinetics_type|           | in     | F        |
!! | ODE_obj    | ODE_ddt                                          | ODE derived data type                   | DDT     |    0 | Solver_type  |           | none   | F        |
!! | j_rateConst| photo_rate_constants                             | photochemical rates constants           | s-1     |    1 | real         | kind_phys | in     | F        |
!! | k_rateConst| gasphase_rate_constants                          | k rate constants                        | s-1     |    1 | real         | kind_phys | inout  | F        |
!! | errmsg     | ccpp_error_message                               | CCPP error message                      | none    |    0 | character    | len=512   | out    | F        |
!! | errflg     | ccpp_error_flag                                  | CCPP error flag                         | flag    |    0 | integer      |           | out    | F        |
!!
subroutine chemistry_driver_run(vmr, TimeStart, TimeEnd, Time, theKinetics, ODE_obj, j_rateConst,  k_rateConst, errmsg, errflg)


  implicit none
!-----------------------------------------------------------
!  these dimension parameters will be set by the cafe/configurator
!-----------------------------------------------------------
  integer :: nTotRxt = 0    ! total number of chemical reactions

  ! Temporary hardwiring of environmental conditions
  
  
  integer            :: i, k, n
  integer            :: errflg          ! error index from CPF
  integer            :: ierr
  real(kind=r8)              :: j_rateConst(:)  ! host model provides photolysis rates for now
  real(kind=r8),pointer      :: k_rateConst(:)  ! host model provides photolysis rates for now
  real(kind=r8)              :: vmr(:)          ! "working" concentration passed thru CPF
  character(len=512) :: errmsg

  integer  :: icntrl(20)     ! integer control array for ODE solver
  real(r8) :: rcntrl(20)     ! real control array for ODE solver
  real(r8) :: TimeStart, TimeEnd, Time, dt
  real(r8), allocatable :: absTol(:), relTol(:)

! declare the types
  type(Solver_type),    pointer  :: ODE_obj
  type(kinetics_type),  pointer  :: theKinetics

!  write(0,*) ' at point run 1'
  call k_rateConst_run(k_rateConst, errflg, errmsg)
!  write(0,*) ' at point  run3'
!  write(0,*) 'inside run, k_rateConst(1)=',k_rateConst(1)
! write(0,*) 'inside run, j_rateConst(1)=',j_rateConst(1)
  call kinetics_run(theKinetics, k_rateConst, j_rateConst, errmsg, errflg)
  call theKinetics%rateConst_print()
!  write(0,*) ' at point run 4'
  call chem_solve_run(TimeStart, TimeEnd, Time, vmr, theKinetics, ODE_obj, errmsg, errflg)
!  write(0,*) ' at point run 5'

end subroutine chemistry_driver_run

subroutine chemistry_driver_finalize()
end subroutine chemistry_driver_finalize
end module chemistry_driver

