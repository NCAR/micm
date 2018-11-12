module prepare_chemistry_mod

use json_loader,            only: json_loader_read
use const_props_mod,        only: const_props_type

implicit none

public :: prepare_chemistry_init, prepare_chemistry_run, prepare_chemistry_finalize

contains

subroutine prepare_chemistry_init(cnst_info, nSpecies, nkRxt, njRxt, nTotRxt)

! This routine reads in the chemistry json file 

  integer, intent(out) :: nSpecies    ! number prognostic constituents
  integer, intent(out) :: nkRxt       ! number gas phase reactions
  integer, intent(out) :: njRxt       ! number of photochemical reactions
  integer, intent(out) :: nTotRxt     ! total number of chemical reactions

  character(len=*), parameter :: jsonfile = '../../MICM_chemistry/generated/terminator/terminator.json'
  type(const_props_type), pointer :: cnst_info(:)

  call json_loader_read( jsonfile, cnst_info, nSpecies, nkRxt, njRxt )

  nTotRxt = nkRxt + njRxt

end subroutine prepare_chemistry_init

subroutine prepare_chemistry_run 
! This routine is intentionally empty 
end subroutine prepare_chemistry_run 

subroutine prepare_chemistry_finalize 
! This routine is intentionally empty 
end subroutine prepare_chemistry_finalize 

end module prepare_chemistry_mod
