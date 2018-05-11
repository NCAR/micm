module external_fields

use precision, only : r8

implicit none


! These will eventually be provided by cpf

real(r8), parameter :: temperature = 273
real(r8), parameter :: pressure = 100000
real(r8), parameter :: mass_density = 1

real(r8) :: some_global_data

contains

  !---------------------------
  ! Set global data 
  !---------------------------
  integer function set_externals()
  
    some_global_data = 3
    set_externals = 1

  end function set_externals
  
  
end module external_fields
