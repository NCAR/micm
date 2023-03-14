module dforce_dy_times_vector_mod
  use iso_c_binding

  implicit none

contains

  subroutine dforce_dy_times_vector(dforce_dy, vector, product) bind(c)
    use iso_c_binding,      only : c_double
    use kinetics_utilities, only : old_dforce_dy_times_vector => dforce_dy_times_vector

    real(kind=c_double), intent(in) :: dforce_dy(:), vector(:)
    real(kind=c_double), pointer, intent(out) :: product(:)

    allocate(product(size(dforce_dy)))

    call old_dforce_dy_times_vector(dforce_dy, vector, product)

  end subroutine dforce_dy_times_vector

end module dforce_dy_times_vector_mod