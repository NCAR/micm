module factored_alpha_minus_jac_mod
  use iso_c_binding

  implicit none

contains

  subroutine factored_alpha_minus_jac(dforce_dy, alpha, LU) bind(c)
    use iso_c_binding,      only : c_double
    use kinetics_utilities, only : old_factored_alpha_minus_jac => factored_alpha_minus_jac

    real(kind=c_double), intent(in) :: dforce_dy(:)
    real(kind=c_double), intent(in), value :: alpha
    real(kind=c_double), pointer, intent(out) :: LU(:)

    allocate(LU(size(dforce_dy)))

    call old_factored_alpha_minus_jac(LU, alpha, dforce_dy)

  end subroutine factored_alpha_minus_jac

end module factored_alpha_minus_jac_mod