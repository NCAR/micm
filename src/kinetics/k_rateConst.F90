module k_rateConst
!--------------------------------------------------------
! k rate constants for 3component chemistry
!--------------------------------------------------------

use machine,         only: r8 => kind_phys

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
!! | local_name | standard_name                                    | long_name                               | units   | rank | type      | kind      | intent | optional |
!! |------------|--------------------------------------------------|-----------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | k_rateConst| gasphase_rate_constants                          | k rate constants                        | s-1     |    1 | real      | kind_phys | inout  | F        |
!! | errmsg     | ccpp_error_message                               | CCPP error message                      | none    |    0 | character | len=512   | out    | F        |
!! | errflg     | ccpp_error_flag                                  | CCPP error flag                         | flag    |    0 | integer   |           | out    | F        |
!!
  subroutine k_rateConst_init(k_rateConst, errflg, errmsg)
      
    real(r8), pointer, intent(inout) :: k_rateConst(:)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg=''
    errflg=0

    ! Nothing for the 3component chemistry to do at init time currently

  end  subroutine k_rateConst_init

  !---------------------------
  ! Compute k_rateConst, given M, P, T
  ! Execute once for the chemistry-time-step advance
  !---------------------------
!> \section arg_table_k_rateConst_run Argument Table
!! | local_name | standard_name                                    | long_name                               | units   | rank | type      | kind      | intent | optional |
!! |------------|--------------------------------------------------|-----------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | k_rateConst| gasphase_rate_constants                          | k rate constants                        | s-1     |    1 | real      | kind_phys | inout  | F        |
!! | temp_arr   | layer_temperature                                | mid-point layer temperature             | K       |    1 | real      | kind_phys | in     | F        |
!! | photo_lev  | level_number_for_photolysis                      | level number used to set j_rateConst    | count   |    0 | integer   |           | none   | F        | 
!! | errmsg     | ccpp_error_message                               | CCPP error message                      | none    |    0 | character | len=512   | out    | F        |
!! | errflg     | ccpp_error_flag                                  | CCPP error flag                         | flag    |    0 | integer   |           | out    | F        |
!!
  subroutine k_rateConst_run(k_rateConst, temp_arr, photo_lev, errflg, errmsg)
  
    real(r8),pointer, intent(inout) :: k_rateConst(:)
    real(r8),           intent(in)  :: temp_arr(:)
    integer,            intent(in)  :: photo_lev
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    real(r8)                        :: TEMP

    ! retrieve the temperature used by tuv (the photolysis level)
    TEMP = temp_arr(photo_lev)

    errmsg=''
    errflg=0
  
! These are probably set by the Chemistry Cafe
#include "k_rateConst.inc"

  end subroutine k_rateConst_run
  
  subroutine k_rateConst_finalize
  end subroutine k_rateConst_finalize
  
end module k_rateConst
