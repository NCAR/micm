module photolysis_interstitial
  use machine, only: rk => kind_phys

  implicit none

  integer :: level_number = 0
contains

!> \section arg_table_photolysis_interstitial_init Argument Table
!! | local_name | standard_name               | long_name                             | units   | rank | type      | kind      | intent | optional |
!! |------------|-----------------------------|---------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | photo_lev  | level_number_for_photolysis | level number used to set j_rateConst  | count   |    0 | integer   |           | none   | F        |
!! | errmsg     | ccpp_error_message          | CCPP error message                    | none    |    0 | character | len=512   | out    | F        |
!! | errflg     | ccpp_error_flag             | CCPP error flag                       | flag    |    0 | integer   |           | out    | F        |
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
!! | local_name | standard_name         | long_name                      | units   | rank | type      | kind      | intent | optional |
!! |------------|-----------------------|--------------------------------|---------|------|-----------|-----------|--------|----------|
!! | prates     | photolysis_rates_col  | photolysis rates column        | s-1     |    2 | real      | kind_phys | in     | F        |
!! | j_rateConst| photo_rate_constants  | photochemical rates constants  | s-1     |    1 | real      | kind_phys | out    | F        |
!! | errmsg     | ccpp_error_message    | CCPP error message             | none    |    0 | character | len=512   | out    | F        |
!! | errflg     | ccpp_error_flag       | CCPP error flag                | flag    |    0 | integer   |           | out    | F        |
!!
  subroutine photolysis_interstitial_run(prates, j_rateConst, errmsg, errflg)
    real(rk),           intent(in)  :: prates(:,:) ! /sec
    real(rk),           intent(out) :: j_rateConst(:) ! /sec
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

#include "prates.inc"    

  end subroutine photolysis_interstitial_run
  
!> \section arg_table_photolysis_interstitial_finalize Argument Table
!! | local_name | standard_name         | long_name                      | units   | rank | type      | kind      | intent | optional |
!! |------------|-----------------------|--------------------------------|---------|------|-----------|-----------|--------|----------|
!! | errmsg     | ccpp_error_message    | CCPP error message             | none    |    0 | character | len=512   | out    | F        |
!! | errflg     | ccpp_error_flag       | CCPP error flag                | flag    |    0 | integer   |           | out    | F        |
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
