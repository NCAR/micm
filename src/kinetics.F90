
module kinetics

  use kinetics_module, only: kinetics_type
  
  implicit none

  private
  public :: kinetics_init 
  public :: kinetics_run
  public :: kinetics_finalize

  integer, parameter :: rk = selected_real_kind( 15 )
  
contains

!> \section arg_table_kinetics_init Argument Table
!! | local_name | standard_name                                    | long_name                               | units   | rank | type      | kind      | intent | optional |
!! |------------|--------------------------------------------------|-----------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | nkRxt      | no_gas_phase_reactions                           | ngas_phase_reactions                    | count   |    0 | integer   |           | none   | F        |
!! | theKinetics | kinetics_data                                   | chemistry kinetics                      | DDT     |    0 | kinetics_type |       | none   | F        |
!! | errmsg     | error_message                                    | CCPP error message                      | none    |    0 | character | len=512   | out    | F        |
!! | errflg     | error_flag                                       | CCPP error flag                         | flag    |    0 | integer   |           | out    | F        |
!!
  subroutine kinetics_init( nkRxt, theKinetics, errmsg, errflg )

    !--- arguments
    integer,            intent(in)  :: nkRxt
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg
    type(kinetics_type), pointer    :: theKinetics

    theKinetics%nkReact = nkRxt
    call theKinetics%k_rateConst_init( nkRxt )

    errmsg = ''
    errflg = 0

  end subroutine kinetics_init

!> \section arg_table_kinetics_run Argument Table
!! | local_name | standard_name                                    | long_name                               | units   | rank | type      | kind      | intent | optional |
!! |------------|--------------------------------------------------|-----------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | theKinetics | kinetics_data                                   | chemistry kinetics                      | DDT     |    0 | kinetics_type |       | none   | F        |
!! | errmsg     | error_message                                    | CCPP error message                      | none    |    0 | character | len=512   | out    | F        |
!! | errflg     | error_flag                                       | CCPP error flag                         | flag    |    0 | integer   |           | out    | F        |
!!
  subroutine kinetics_run( theKinetics, errmsg, errflg )

    !--- arguments
    character(len=512), intent(out)   :: errmsg
    integer,            intent(out)   :: errflg
    type(kinetics_type), pointer      :: theKinetics

    !--- local variables
    integer :: Ierr

    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    !--- set the gas phase rate constants
    call theKinetics%k_rateConst_run()

  end subroutine kinetics_run

!> \section arg_table_kinetics_finalize Argument Table
!! | local_name | standard_name                                    | long_name                               | units   | rank | type      | kind      | intent | optional |
!! |------------|--------------------------------------------------|-----------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | errmsg     | error_message                                    | CCPP error message                      | none    |    0 | character | len=512   | out    | F        |
!! | errflg     | error_flag                                       | CCPP error flag                         | flag    |    0 | integer   |           | out    | F        |
!!
  subroutine kinetics_finalize( errmsg, errflg )

    !--- arguments
    character(len=512), intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

  end subroutine kinetics_finalize

end module kinetics
