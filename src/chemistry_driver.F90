module chemistry_driver

use kinetics_module,  only       : kinetics_type
use kinetics,         only       : kinetics_init, kinetics_run
use chem_solve,       only       : chem_solve_run
use machine,          only       : r8 => kind_phys
use const_props_mod,  only       : const_props_type
use prepare_chemistry_mod, only  : prepare_chemistry_init

implicit none

contains


!> \section arg_table_chemistry_driver_init Argument Table
!! | local_name | standard_name                                    | long_name                               | units   | rank | type            | kind      | intent | optional |
!! |------------|--------------------------------------------------|-----------------------------------------|---------|------|-----------------|-----------|--------|----------|
!! | TimeStart  | chem_step_start_time                             | Chem step start time                    | s       |    0 | real            | kind_phys | in     | F        |
!! | TimeEnd    | chem_step_end_time                               | Chem step end time                      | s       |    0 | real            | kind_phys | in     | F        |
!! | dt         | time_step_for_physics                            | time_step_for_physics                   | s       |    0 | real            | kind_phys | in     | F        |
!! | theKinetics | kinetics_data                                   | chemistry kinetics                      | DDT     |    0 | kinetics_type   |           | out    | F        |
!! | k_rateConst| gasphase_rate_constants                          | k rate constants                        | s-1     |    1 | real            | kind_phys | inout  | F        |
!! | AbsTol     | abs_trunc_error                                  | ODE absolute step truncation error      | none    |    1 | real            | kind_phys | in     | F        |
!! | RelTol     | rel_trunc_error                                  | ODE relative step truncation error      | none    |    1 | real            | kind_phys | in     | F        |
!! | icntrl     | ODE_icontrol                                     | ODE integer controls                    | flag    |    1 | integer         |           | inout  | F        |
!! | rcntrl     | ODE_rcontrol                                     | ODE real controls                       | none    |    1 | real            | kind_phys | inout  | F        |
!! | cnst_info  | chemistry_constituent_info                       | chemistry_constituent_info              | DDT     |    1 | const_props_type|           | out    | F        |
!! | errmsg     | ccpp_error_message                               | CCPP error message                      | none    |    0 | character       | len=512   | out    | F        |
!! | errflg     | ccpp_error_flag                                  | CCPP error flag                         | flag    |    0 | integer         |           | out    | F        |
!!
subroutine chemistry_driver_init(TimeStart,TimeEnd, dt, theKinetics,  k_rateConst,  AbsTol, RelTol, icntrl, rcntrl, cnst_info, errmsg, errflg)

  implicit none
!-----------------------------------------------------------
!  these dimension parameters will be set by the cafe/configurator
!-----------------------------------------------------------
  real(r8), intent(in)           :: TimeStart, TimeEnd
  real(r8), intent(in)           :: dt
  type(kinetics_type),  pointer  :: theKinetics
  real(kind=r8),intent(inout)    :: k_rateConst(:)  ! host model provides photolysis rates for now
  real(r8), intent(inout)        :: AbsTol(:), RelTol(:)
  integer, intent(inout)         :: icntrl(:)     ! integer control array for ODE solver
  real(r8), intent(inout)        :: rcntrl(:)     ! real control array for ODE solver
  character(len=512),intent(out) :: errmsg
  integer, intent(out)           :: errflg          ! error index from CPF

  integer            :: nTotRxt    ! total number of chemical reactions
  integer            :: nSpecies   ! number prognostic constituents
  integer            :: nkRxt      ! number gas phase reactions
  integer            :: njRxt      ! number of photochemical reactions

! declare the types
  type(const_props_type), pointer :: cnst_info(:)

  character(len=16) :: cnst_name

  write(0,*) ' Entered chemistry_driver_init'

!   This routine should be only called here when the main program no longer needs to allocate variables
    call prepare_chemistry_init(cnst_info, nSpecies, nkRxt, njRxt)

    nTotRxt =  nkRxt + njRxt

!-----------------------------------------------------------
!  initialize ode solver "control" variable defaults
!-----------------------------------------------------------
  absTol(:) = 1.e-9_r8
  relTol(:) = 1.e-4_r8
  icntrl(:) = 0
  rcntrl(:) = 0._r8

!-----------------------------------------------------------
!  set ode solver "control" variables for Rosenbrock solver
!-----------------------------------------------------------
  icntrl(1) = 1                                 ! autonomous, F depends only on Y
  icntrl(3) = 2                                 ! ros3 solver
  rcntrl(2) = dt                                ! Hmax
  rcntrl(3) = .01_r8*dt                         ! Hstart

  write(*,*) ' '
  write(*,*) 'icntrl settings'
  write(*,'(10i6)') icntrl(1:10)
  write(*,*) 'rcntrl settings'
  write(*,'(1p,10(1x,g0))') rcntrl(1:10)
  write(*,*) ' '


!-----------------------------------------------------------
!  initialize the kinetics
!-----------------------------------------------------------
  call kinetics_init(nTotRxt, theKinetics, errmsg, errflg)

end subroutine chemistry_driver_init

!> \section arg_table_chemistry_driver_run Argument Table
!! | local_name | standard_name                                    | long_name                               | units   | rank | type         | kind      | intent | optional |
!! |------------|--------------------------------------------------|-----------------------------------------|---------|------|--------------|-----------|--------|----------|
!! | vmr        | concentration                                    | species concentration                   | mole/mole |    1 |real        | kind_phys | inout  | F        |
!! | TimeStart  | chem_step_start_time                             | Chem step start time                    | s       |    0 | real         | kind_phys | in     | F        |
!! | TimeEnd    | chem_step_end_time                               | Chem step end time                      | s       |    0 | real         | kind_phys | in     | F        |
!! | theKinetics | kinetics_data                                   | chemistry kinetics                      | DDT     |    0 | kinetics_type|           | in     | F        |
!! | j_rateConst| photo_rate_constants                             | photochemical rates constants           | s-1     |    1 | real         | kind_phys | in     | F        |
!! | k_rateConst| gasphase_rate_constants                          | k rate constants                        | s-1     |    1 | real         | kind_phys | in     | F        |
!! | AbsTol     | abs_trunc_error                                  | ODE absolute step truncation error      | none    |    1 | real         | kind_phys | in     | F        |
!! | RelTol     | rel_trunc_error                                  | ODE relative step truncation error      | none    |    1 | real         | kind_phys | in     | F        |
!! | icntrl     | ODE_icontrol                                     | ODE integer controls                    | flag    |    1 | integer      |           | in     | F        |
!! | rcntrl     | ODE_rcontrol                                     | ODE real controls                       | none    |    1 | real         | kind_phys | in     | F        |
!! | errmsg     | ccpp_error_message                               | CCPP error message                      | none    |    0 | character    | len=512   | out    | F        |
!! | errflg     | ccpp_error_flag                                  | CCPP error flag                         | flag    |    0 | integer      |           | out    | F        |
!!
subroutine chemistry_driver_run(vmr, TimeStart, TimeEnd, theKinetics,  j_rateConst,  k_rateConst, AbsTol, RelTol, &
                                icntrl, rcntrl, errmsg, errflg)


  implicit none
!-----------------------------------------------------------
!  these dimension parameters will be set by the cafe/configurator
!-----------------------------------------------------------

  real(kind=r8), intent(inout)    :: vmr(:)                ! "working" concentration passed thru CPF
  real(r8), intent(in)            :: TimeStart, TimeEnd
  type(kinetics_type),  pointer   :: theKinetics
  real(kind=r8), intent(in)       :: j_rateConst(:)        ! host model provides photolysis rates for now
  real(kind=r8), intent(in)       :: k_rateConst(:)        ! host model provides photolysis rates for now
  real(r8), intent(in)            :: absTol(:), relTol(:)
  integer,intent(in)              :: icntrl(:)             ! integer control array for ODE solver
  real(r8),intent(in)             :: rcntrl(:)             ! real control array for ODE solver
  character(len=512), intent(out) :: errmsg
  integer, intent(out)            :: errflg                ! error index from CPF

  integer            :: i, k, n

!-----------------------------------------------------------
!  update the kinetics
!-----------------------------------------------------------
  call kinetics_run(theKinetics, k_rateConst, j_rateConst, errmsg, errflg)
  if (errflg /= 0) return

  call theKinetics%rateConst_print()

!-----------------------------------------------------------
!  solve the current timestep's chemistry
!-----------------------------------------------------------
  call chem_solve_run(TimeStart, TimeEnd, vmr, theKinetics, icntrl, rcntrl, AbsTol, RelTol, errmsg, errflg)
  if (errflg /= 0) return

end subroutine chemistry_driver_run

subroutine chemistry_driver_finalize()
end subroutine chemistry_driver_finalize

end module chemistry_driver

