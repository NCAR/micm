
module chem_solve

  use solver_var_defs, only: Solver_type
  use kinetics_module, only: kinetics_type
  use machine,         only: rk => kind_phys

  
  implicit none

  private
  public :: chem_solve_init 
  public :: chem_solve_finalize
  public :: chem_solve_run

contains

!> \section arg_table_chem_solve_init Argument Table
!! | local_name | standard_name                                    | long_name                               | units   | rank | type      | kind      | intent | optional |
!! |------------|--------------------------------------------------|-----------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | icntrl     | ODE_icontrol                                     | ODE integer controls                    | flag      |    1 | integer     |           | in     | F        |
!! | rcntrl     | ODE_rcontrol                                     | ODE real controls                       | none      |    1 | real        | kind_phys | in     | F        |
!! | AbsTol     | abs_trunc_error                                  | ODE absolute step truncation error      | none      |    1 | real        | kind_phys | in     | F        |
!! | RelTol     | rel_trunc_error                                  | ODE relative step truncation error      | none      |    1 | real        | kind_phys | in     | F        |
!! | TimeStart  | chem_step_start_time                             | Chem step start time                    | s         |    0 | real        | kind_phys | in     | F        |
!! | TimeEnd    | chem_step_end_time                               | Chem step end time                      | s         |    0 | real        | kind_phys | in     | F        |
!! | errmsg     | ccpp_error_message                               | CCPP error message                      | none    |    0 | character | len=512   | out    | F        |
!! | errflg     | ccpp_error_flag                                  | CCPP error flag                         | flag    |    0 | integer   |           | out    | F        |
!!
  subroutine chem_solve_init( TimeStart, TimeEnd, AbsTol, RelTol, &
                              icntrl, rcntrl,  errmsg, errflg)

    implicit none

    !--- arguments
    integer,  intent(in)            :: icntrl(:)
    real(rk), intent(in)            :: rcntrl(:)
    real(rk), intent(in)            :: TimeStart
    real(rk), intent(in)            :: TimeEnd
    real(rk), intent(in)            :: AbsTol(:)
    real(rk), intent(in)            :: RelTol(:)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

  end subroutine chem_solve_init

!> \section arg_table_chem_solve_run Argument Table
!! | local_name | standard_name                                    | long_name                               | units   | rank | type      | kind      | intent | optional |
!! |------------|--------------------------------------------------|-----------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | TimeStart  | chem_step_start_time                             | Chem step start time                    | s       |    0 | real          | kind_phys | in     | F        |
!! | TimeEnd    | chem_step_end_time                               | Chem step end time                      | s       |    0 | real          | kind_phys | in     | F        |
!! | vmr        | concentration                                    | species concentration                   | mole/mole |    1 | real        | kind_phys | inout  | F        |
!! | theKinetics | kinetics_data                                   | chemistry kinetics                      | DDT     |    0 | kinetics_type |           | none   | F        |
!! | icntrl     | ODE_icontrol                                     | ODE integer controls                    | flag      |    1 | integer     |           | in     | F        |
!! | rcntrl     | ODE_rcontrol                                     | ODE real controls                       | none      |    1 | real        | kind_phys | in     | F        |
!! | AbsTol     | abs_trunc_error                                  | ODE absolute step truncation error      | none      |    1 | real        | kind_phys | in     | F        |
!! | RelTol     | rel_trunc_error                                  | ODE relative step truncation error      | none      |    1 | real        | kind_phys | in     | F        |
!! | errmsg     | ccpp_error_message                               | CCPP error message                      | none    |    0 | character     | len=512   | out    | F        |
!! | errflg     | ccpp_error_flag                                  | CCPP error flag                         | flag    |    0 | integer       |           | out    | F        |
!!
  subroutine chem_solve_run ( TimeStart, TimeEnd, vmr, theKinetics, icntrl, rcntrl, AbsTol, RelTol, errmsg, errflg)

  use Rosenbrock_Integrator, only: rosenbrock_run

    implicit none

    !--- arguments
    real(rk), intent(in)         :: TimeStart
    real(rk), intent(in)         :: TimeEnd
    real(rk), intent(inout)      :: vmr(:)
    type(kinetics_type), pointer :: theKinetics
    integer,  intent(in)            :: icntrl(:)
    real(rk), intent(in)            :: rcntrl(:)
    real(rk), intent(in)            :: AbsTol(:)
    real(rk), intent(in)            :: RelTol(:)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer :: nSpecies
    real(rk) :: rstatus(20)
    integer :: istatus(20)

    !--- local variables
    integer :: Ierr
 
    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    nSpecies = size(AbsTol)
    !--- initialize intent(out) variables
    ! initialize all intent(out) variables here

    !--- actual code
    ! add your code here

    call Rosenbrock_run(nSpecies,vmr,TimeStart,TimeEnd, &
           theKinetics, AbsTol,RelTol,&
           rcntrl,icntrl,RSTATUS,ISTATUS,IERR)

    ! in case of errors, set errflg to a value != 0,
    ! create a meaningfull error message and return

  end subroutine chem_solve_run

!> \section arg_table_chem_solve_finalize Argument Table
!! | local_name | standard_name                                    | long_name                               | units   | rank | type      | kind      | intent | optional |
!! |------------|--------------------------------------------------|-----------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | errmsg     | ccpp_error_message                               | CCPP error message                      | none    |    0 | character | len=512   | out    | F        |
!! | errflg     | ccpp_error_flag                                  | CCPP error flag                         | flag    |    0 | integer   |           | out    | F        |
!!
  subroutine chem_solve_finalize( errmsg, errflg )

    !--- arguments
    character(len=512), intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    write(*,*) 'chem_solve_finalize: CALLED'

  end subroutine chem_solve_finalize

end module chem_solve
