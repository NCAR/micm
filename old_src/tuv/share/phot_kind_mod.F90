module phot_kind_mod

  implicit none

  integer, parameter :: dp = selected_real_kind(14,300)

  integer,parameter :: kind_phot = selected_real_kind(12) ! 8 byte real

end module phot_kind_mod
