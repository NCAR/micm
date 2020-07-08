module environmental_state_mod

  type environmental_state_type
    real :: number_density_air, temperature, pressure 
    real :: aerosol_surface_area_density(4), aerosol_diameter(4)
    real :: h2ovmr, o2_number_density
  end type environmental_state_type

end module environmental_state_mod
