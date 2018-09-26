! Rates
   ! Y0_a
   rates(1) = rateConstants(1) * vmr(1)
   ! Y1_Y2_M_b
   rates(2) = rateConstants(2) * vmr(2) * vmr(3)
   ! Y1_Y1_a
   rates(3) = rateConstants(3) * vmr(2) * vmr(2)
