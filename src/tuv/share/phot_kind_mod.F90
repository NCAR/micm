module phot_kind_mod

use ccpp_kind, only: kind_phot => kind_phys

  implicit none

  integer, parameter :: dp = selected_real_kind(14,300)


end module phot_kind_mod
