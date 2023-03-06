module dforce_dy_mod
  use iso_c_binding

  implicit none

contains

  subroutine dforce_dy(rate_constants, number_densities, number_density_air, result) bind(c)
    use iso_c_binding,              only : c_double
    use kinetics_utilities, only : p_dforce_dy=>dforce_dy

    real(kind=c_double), value :: number_density_air
    real(kind=c_double), intent(in) :: rate_constants(:), number_densities(:)
    real(kind=c_double), pointer, intent(out) :: result(:)

    allocate(result(23))

    call p_dforce_dy(result, rate_constants, number_densities, number_density_air)

  end subroutine dforce_dy

end module dforce_dy_mod