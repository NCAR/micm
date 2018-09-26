   ! set the forcing array
    force(1) = rates(2) - rates(1)
    force(2) = rates(1) - rates(2) - 2._r8 * rates(3)
    force(3) = 2._r8 * rates(3)
