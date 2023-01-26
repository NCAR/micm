program rate_constants

  use micm_environment,             only : environment_t
  use musica_constants,             only : musica_dk
  use micm_rate_constant_arrhenius, only : rate_constant_arrhenius_t
  implicit none

  type(rate_constant_arrhenius_t) :: rate
  type(environment_t)             :: env

  env%temperature = 1
  env%pressure = 1

  rate = rate_constant_arrhenius_t( A=1.0_musica_dk, B=1.0_musica_dk,         &
    C=1.0_musica_dk, D=1.0_musica_dk, E=1.0_musica_dk)

  print *, rate%calculate(env)

end program rate_constants