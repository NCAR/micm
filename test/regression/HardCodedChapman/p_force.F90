module p_force_mod
  use iso_c_binding

  implicit none

contains

  subroutine p_force(rate_constants, number_densities, number_density_air, force) bind(c)
    use iso_c_binding,      only : c_double
    use kinetics_utilities, only : old_p_force => p_force

    real(kind=c_double), intent(in) :: rate_constants(:), number_densities(:)
    real(kind=c_double), intent(in), value :: number_density_air
    real(kind=c_double), pointer, intent(out) :: force(:)

    allocate(force(size(number_densities)))

    call old_p_force(rate_constants, number_densities, number_density_air, force)

  end subroutine p_force

end module p_force_mod