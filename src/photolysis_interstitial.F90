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
    use rate_constants_utility, only: p_rate_mapping

    real(kind_phys),           intent(in)  :: prates(:,:) ! /sec
    real(kind_phys),           intent(out) :: j_rateConst(:) ! /sec
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    call p_rate_mapping(prates(level_number,:), j_rateConst)

  end subroutine photolysis_interstitial_run
end module photolysis_interstitial
