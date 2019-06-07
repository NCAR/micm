module photolysis_interstitial
! use ccpp_kinds, only: rk => kind_phys
 use ccpp_kinds, only: kind_phys

  implicit none

  integer :: level_number = 0
contains

!> \section arg_table_photolysis_interstitial_init Argument Table
!! \htmlinclude photolysis_interstitial_init.html
!!
  subroutine photolysis_interstitial_init(photo_lev, errmsg, errflg)
    integer,            intent(in)  :: photo_lev
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    level_number = photo_lev
  end subroutine photolysis_interstitial_init

!> \section arg_table_photolysis_interstitial_run Argument Table
!! \htmlinclude photolysis_interstitial_run.html
!!
  subroutine photolysis_interstitial_run(prates, j_rateConst, errmsg, errflg)
    real(kind_phys),           intent(in)  :: prates(:,:) ! /sec
    real(kind_phys),           intent(out) :: j_rateConst(:) ! /sec
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

!TUV Mapping rates

!    O2 ~ TUV_rate(O2 -> O + O)
    j_rateConst(1) = 1*prates(level_number,1)

!    O3 ~ TUV_rate(O3 -> O2 + O(1D))
    j_rateConst(2) = 1*prates(level_number,2)

!    O3 ~ TUV_rate(O3 -> O2 + O(3P))
    j_rateConst(3) = 1*prates(level_number,3)

  end subroutine photolysis_interstitial_run
  
!> \section arg_table_photolysis_interstitial_finalize Argument Table
!! \htmlinclude photolysis_interstitial_finalize.html
!!
  subroutine photolysis_interstitial_finalize( errmsg, errflg )

    !--- arguments
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine photolysis_interstitial_finalize
  
end module photolysis_interstitial
