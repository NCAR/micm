    ! Jacobian
    jac(:,:) = 0._r8

    ! Jacobian
    jac(1,1)  = jac(1,1) - rateConstants(1)
    jac(1,2)  = jac(1,2) + rateConstants(2) * vmr(3)
    jac(1,3)  = jac(1,3) + rateConstants(2) * vmr(2)

    jac(2,1)  = jac(2,1) + rateConstants(1)
    jac(2,2)  = jac(2,2) - rateConstants(2) * vmr(3)
    jac(2,3)  = jac(2,3) - rateConstants(2) * vmr(2)
    jac(2,2)  = jac(2,2) - 4._r8 * rateConstants(3) * vmr(2)

    jac(3,2)  = jac(3,2) + 4._r8 * rateConstants(3) * vmr(2)

