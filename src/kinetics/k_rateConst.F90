module k_rateConst
!--------------------------------------------------------
! k rate constants for 3component chemistry
!--------------------------------------------------------

  use ccpp_kinds, only: r8 => kind_phys

implicit none
private

public :: k_rateConst_init
public :: k_rateConst_run
public :: k_rateConst_finalize

! k_rateConst are computed at the beginning of the 
!   chemistry_box_solver time step.
!   They are not computed for internal steps of the
!   box-model time-step advancer
! k_rateConst will be thread-safe memory provided elsewhere.
! rate_constant_store will be an accessor to memory
! For now, it is allocated here. It is not thread safe

contains

!> \section arg_table_k_rateConst_init Argument Table
!! \htmlinclude k_rateConst_init.html
!!
  subroutine k_rateConst_init(errmsg, errflg)
      
    integer,            intent(out) :: errflg
    character(len=512), intent(out) :: errmsg

    errmsg=''
    errflg=0

    ! Nothing for the 3component chemistry to do at init time currently

  end  subroutine k_rateConst_init

  !---------------------------
  ! Compute k_rateConst, given M, P, T
  ! Execute once for the chemistry-time-step advance
  !---------------------------
!> \section arg_table_k_rateConst_run Argument Table
!! \htmlinclude k_rateConst_run.html
!!
  subroutine k_rateConst_run(k_rateConst, c_m, rh, c_h2o, temp, errmsg, errflg)
  
    real(r8),pointer, intent(out) :: k_rateConst(:)
    real(r8),           intent(in)  :: c_m
    real(r8),           intent(in)  :: rh 
    real(r8),           intent(in)  :: c_h2o     
    real(r8),           intent(in)  :: TEMP
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg


    ! retrieve the temperature used by tuv (the photolysis level)

    errmsg=''
    errflg=0
  
! These are probably set by the Chemistry Cafe

! Rate Constants

! N2_O1D
k_rateConst(1) =  2.150000e-11_r8 * exp(110.00_r8 / TEMP) 

! O1D_O2
k_rateConst(2) =  3.300000e-11_r8 * exp(55.00_r8 / TEMP) 

! O_O3
k_rateConst(3) =  8.000000e-12_r8 * exp(-2060.00_r8 / TEMP) 

! O_O2_M
k_rateConst(4) = usr_O_O2( temp )

  end subroutine k_rateConst_run
  
  subroutine k_rateConst_finalize
  end subroutine k_rateConst_finalize


! number of Functions: 1



REAL(KIND=r8) FUNCTION usr_O_O2( temp )
! for cesm-consistent reaction labels
! O+O2+M -> O3+M

    REAL(KIND=r8), INTENT(IN) :: temp

    usr_O_O2 = 6.00e-34_r8*(temp/300._r8)**(-2.4_r8)

END FUNCTION usr_O_O2

! Included shim as a number of WRF functions depend on this
!-------------------------------------------
! Troe equilibrium reactions (as in Stockwell et al, 1997)

    real(kind=r8) FUNCTION TROEE(A, B, k0_300K,n,kinf_300K,m,temp,cair)

    INTRINSIC LOG10

    real(kind=r8), INTENT(IN) :: temp      ! temperature [K]
    real(kind=r8), INTENT(IN) :: cair      ! air concentration [molecules/cm3]
    real(kind=r8),     INTENT(IN) :: k0_300K   ! low pressure limit at 300 K
    real(kind=r8),     INTENT(IN) :: n         ! exponent for low pressure limit
    real(kind=r8),     INTENT(IN) :: kinf_300K ! high pressure limit at 300 K
    real(kind=r8),     INTENT(IN) :: m         ! exponent for high pressure limit
    real(kind=r8),     INTENT(IN) :: A, B
    real(kind=r8)             :: zt_help, k0_T, kinf_T, k_ratio, troe


    zt_help = 300._r8/temp
    k0_T    = k0_300K   * zt_help**(n) * cair ! k_0   at current T
    kinf_T  = kinf_300K * zt_help**(m)        ! k_inf at current T
    k_ratio = k0_T/kinf_T
    troe   = k0_T/(1._r8+k_ratio)*0.6_r8**(1._r8/(1._r8+LOG10(k_ratio)**2))

    TROEE = A * EXP( - B / temp) * troe



  END FUNCTION TROEE
  
end module k_rateConst
