module prepare_chemistry

use json_loader,            only: json_loader_read
use const_props_mod,        only: const_props_type

implicit none

public :: prepare_chemistry_init
!public :: prepare_chemistry_run

contains

!> \section arg_table_prepare_chemistry_init Argument Table
!! \htmlinclude prepare_chemistry_init.html
!!
subroutine prepare_chemistry_init(cnst_info, model_name, nSpecies, nkRxt, njRxt, errmsg, errflg)

! This routine reads in the chemistry json file

  type(const_props_type), pointer, intent(out) :: cnst_info(:)
  character(len=80), intent(out) :: model_name
  integer,           intent(out) :: nSpecies    ! number prognostic constituents
  integer,           intent(out) :: nkRxt       ! number gas phase reactions
  integer,           intent(out) :: njRxt       ! number of photochemical reactions
  character(len=512),intent(out) :: errmsg
  integer,           intent(out) :: errflg      ! error index for CPF


  character(len=120) :: jsonfile

  ! Read chemistry namelist (eventually)
  ! Will include model_name - hardcoded for now

  errmsg = ''
  errflg = 0
  model_name = 'Chapman_v3_1547831703456'

  jsonfile = '../../MICM_chemistry/generated/'//trim(model_name)//'/molec_info.json'
  call json_loader_read( jsonfile, cnst_info, nSpecies, nkRxt, njRxt )

end subroutine prepare_chemistry_init

!> \section arg_table_prepare_chemistry_run Argument Table
!! \htmlinclude prepare_chemistry_run.html
!!
subroutine prepare_chemistry_run(errmsg, errflg)

! This routine does nothing

  character(len=512),intent(out) :: errmsg
  integer,           intent(out) :: errflg      ! error index for CPF


  errmsg = ''
  errflg = 0

end subroutine prepare_chemistry_run


end module prepare_chemistry
