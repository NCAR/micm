module kinetics

  use kinetics_module, only: kinetics_type
  use ccpp_kinds, only: rk => kind_phys

  
  implicit none

  private
  public :: kinetics_init 
  public :: kinetics_run
  public :: kinetics_finalize


  
contains

!> \section arg_table_kinetics_init Argument Table
!! \htmlinclude kinetics_init.html
!!
  subroutine kinetics_init( nTotRxt, theKinetics, errmsg, errflg )

    !--- arguments
    integer,            intent(in)  :: nTotRxt
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg
    type(kinetics_type), pointer, intent(inout)    :: theKinetics

    call theKinetics%rateConst_init( nTotRxt )

    errmsg = ''
    errflg = 0

  end subroutine kinetics_init

!> \section arg_table_kinetics_run Argument Table
!! \htmlinclude kinetics_run.html
!!
  subroutine kinetics_run( theKinetics, k_rateConst, j_rateConst, c_m, errmsg, errflg )

    !--- arguments
    type(kinetics_type), pointer, intent(inout)      :: theKinetics
    real(rk),           intent(in)    :: k_rateConst(:)
    real(rk),           intent(in)    :: j_rateConst(:)
    real(rk),           intent(in)    :: c_m ! total number density
    character(len=512), intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    !--- local variables
    integer :: Ierr

    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    !--- set the gas phase rate constants
    call theKinetics%rateConst_update( k_rateConst, j_rateConst, c_m)

  end subroutine kinetics_run

!> \section arg_table_kinetics_finalize Argument Table
!! \htmlinclude kinetics_finalize.html
!!
  subroutine kinetics_finalize( errmsg, errflg )

    !--- arguments
    character(len=512), intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    errmsg = ''
    errflg = 0

  end subroutine kinetics_finalize

end module kinetics
