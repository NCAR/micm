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
!! | errmsg     | ccpp_error_message                               | CCPP error message                      | none    |    0 | character | len=512   | out    | F        |
!! | errflg     | ccpp_error_flag                                  | CCPP error flag                         | flag    |    0 | integer   |           | out    | F        |
!!
  subroutine k_rateConst_init(errflg, errmsg)
      
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
!! | c_m        | total_number_density                             | total number density                    | molecules/cm3 | 0    | real| kind_phys | in     | F        |
!! | rh         | relative humidity                                | relative humidity                       | percent |    0 | real      | kind_phys | in     | F        |
!! | c_h2o      | water_vapor                                      | water_vapor                             | mole/mole |  0 | real      | kind_phys | in     | F        |
!! | TEMP       | temperature                                      | mid-point layer temperature             | K       |    0 | real      | kind_phys | in     | F        |
!! | errmsg     | ccpp_error_message                               | CCPP error message                      | none    |    0 | character | len=512   | out    | F        |
!! | errflg     | ccpp_error_flag                                  | CCPP error flag                         | flag    |    0 | integer   |           | out    | F        |
!!
  subroutine k_rateConst_run(k_rateConst, c_m, rh, c_h2o, temp, errflg, errmsg)
  
    real(r8),pointer, intent(inout) :: k_rateConst(:)
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
#include "k_rateConst.inc"

  end subroutine k_rateConst_run
  
  subroutine k_rateConst_finalize
  end subroutine k_rateConst_finalize

#include "rate_functions.inc"
  
end module k_rateConst
