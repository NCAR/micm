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

    integer :: i

    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    j_rateConst(:) = prates(level_number,:)

    print *, 'j_rateConst'
    do i = 1, size(j_rateConst(:))
    print *, i, j_rateConst(i)
    end do

  end subroutine photolysis_interstitial_run
end module photolysis_interstitial
