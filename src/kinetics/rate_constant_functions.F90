module rate_constant_functions

  use environmental_state_mod

implicit none

  ! sufficient for accuracy of rate constants
  ! no need to bring in other packages.
  real, parameter :: boltzman = 1.38065e-23 !J/K/molecule
  real, parameter :: boltzman_cgs = boltzman * 1.e7 ! erg (cgs)
  real, parameter :: avogadro = 6.02214076e23 ! molecules / mole
  real, parameter :: c_pi = 3.14159265358979
  real, parameter :: c_300 = 300.
  real, parameter :: c_1 = 1.


  type, abstract :: rate_param_type
  end type rate_param_type

  type, extends( rate_param_type ) :: arrhenius_rate_param_type
     real :: A, C, D, B, E
  end type arrhenius_rate_param_type

! Troe and Troe-like
  type, extends( rate_param_type ) :: Troe_rate_param_type
     real :: A_k0, B_k0, C_K0, A_Kinf, B_Kinf, C_Kinf, F_C
  end type Troe_rate_param_type

  type, extends( rate_param_type ) :: Troe_low_pressure_rate_param_type
     real :: A_k0, B_k0, C_k0
  end type Troe_low_pressure_rate_param_type

  type, extends( rate_param_type ) :: Troe_Reverse_rate_param_type
     real :: A_k0, B_k0, C_K0, A_Kinf, B_Kinf, C_Kinf, F_C, A_r, C_r
  end type Troe_Reverse_rate_param_type

  type, extends( rate_param_type ) :: Troe_chemical_activation_rate_param_type
     real :: A_k0, B_k0, C_K0, A_Kinf, B_Kinf, C_Kinf, F_C
  end type Troe_chemical_activation_rate_param_type

  type, extends( rate_param_type ) :: MCO3_rate_param_type
     real :: A_k0, B_k0, C_K0, A_Kinf, B_Kinf, C_Kinf, F_C
  end type MCO3_rate_param_type

! het rates
  type, extends( rate_param_type) :: Simple_aerosol_heterogeneous_rate_rate_param_type
     real :: gamma, Dg, m_g
  end type Simple_aerosol_heterogeneous_rate_rate_param_type

  type, extends( rate_param_type) :: HET_glyoxyl_rate_param_type
     real :: gamma
  end type HET_glyoxyl_rate_param_type

! specialized 2-body (with possible M) and their reverse rates
  type, extends( rate_param_type ) :: combined_CO_OH_rate_param_type
     real :: A, B
  end type combined_CO_OH_rate_param_type

! specialized rate constants
  type, extends( rate_param_type ) :: CH3COCH3_OH_rate_param_type
     real :: A, B, C
  end type CH3COCH3_OH_rate_param_type

  type, extends( rate_param_type) :: HO2_HO2_rate_param_type
     real :: A_k0, C_k0, A_kinf, C_kinf, F
  end type HO2_HO2_rate_param_type

  type, extends( rate_param_type) :: DMS_OH_rate_param_type
     real :: A, B, C, D
  end type DMS_OH_rate_param_type

  type, extends( rate_param_type) :: HNO3_OH_rate_param_type
     real :: A_k0, C_k0, A_k1, C_k1, A_k2, C_k2
  end type HNO3_OH_rate_param_type

  type, extends( rate_param_type) :: SO2_OH_rate_param_type
     real :: A, B, C, F
  end type SO2_OH_rate_param_type

  type, extends( rate_param_type) :: MCO3_NO2_rate_param_type
     real :: A
  end type MCO3_NO2_rate_param_type

  type, extends( rate_param_type) :: MPAN_M_rate_param_type
     real :: A_f, A_r, B_r
  end type MPAN_M_rate_param_type


contains

function arrhenius(params, state) result(rate_constant)
  type(arrhenius_rate_param_type), intent(in):: params
  type(environmental_state_type), intent(in):: state
  real :: rate_constant

  rate_constant = params%A &
      * exp( -params%C / state%temperature) &
      * ( state%temperature / params%D) ** params%B &
      * ( c_1 + params%E *  state%pressure)

end function arrhenius

function Troe(params, state) result(rate_constant)
  type(Troe_rate_param_type), intent(in) :: params
  type(environmental_state_type), intent(in):: state
  real :: rate_constant

  real :: k_0, k_inf, power_of_F, log10_term

  k_0 =  params%A_k0 &
      * exp( -params%C_k0 / state%temperature) &
      * ( state%temperature / c_300 ) ** params%B_k0

  k_inf =  params%A_Kinf &
      * exp( -params%C_Kinf / state%temperature) &
      * ( state%temperature / c_300 ) ** params%B_Kinf

  log10_term = log10(k_0 * state%number_density_air / k_inf)
  power_of_F =  c_1 + log10_term * log10_term
  
  rate_constant = k_0 / ( c_1 + k_0 * state%number_density_air / k_inf)  &
      *  params%F_C ** ( c_1/power_of_F)

end function Troe

function Troe_low_pressure(params, state) result(rate_constant)
  type(Troe_low_pressure_rate_param_type), intent(in) :: params
  type(environmental_state_type), intent(in):: state
  real :: rate_constant

  real :: k_0, k_inf, power_of_F, log10_term, reverse_factor

  rate_constant =  params%A_k0 &
      * exp( -params%C_k0 / state%temperature) &
      * ( state%temperature / c_300 ) ** params%B_k0

end function Troe_low_pressure

function Troe_reverse(params, state) result(rate_constant)
  type(Troe_reverse_rate_param_type), intent(in) :: params
  type(environmental_state_type), intent(in):: state
  real :: rate_constant

  real :: k_0, k_inf, power_of_F, log10_term, reverse_factor

  k_0 =  params%A_k0 &
      * exp( -params%C_k0 / state%temperature) &
      * ( state%temperature / c_300 ) ** params%B_k0

  k_inf =  params%A_Kinf &
      * exp( -params%C_Kinf / state%temperature) &
      * ( state%temperature / c_300 ) ** params%B_Kinf

  log10_term = log10(k_0 * state%number_density_air / k_inf)
  power_of_F =  c_1 + log10_term * log10_term

  reverse_factor = params%A_r*exp( -params%C_r / state%temperature)

  rate_constant = reverse_factor &
      *  k_0 / ( c_1 + k_0 * state%number_density_air / k_inf) &
      *  params%F_C ** ( c_1/power_of_F)

end function Troe_reverse


function Troe_chemical_activation(params, state) result(rate_constant)
  type(Troe_chemical_activation_rate_param_type), intent(in) :: params
  type(environmental_state_type), intent(in):: state
  real :: rate_constant

  real :: k_0, k_inf, power_of_F, log10_term

  k_0 =  params%A_k0 &
      * exp( -params%C_k0 / state%temperature) &
      * ( state%temperature / c_300 ) ** params%B_k0

  k_inf =  params%A_Kinf &
      * exp( -params%C_Kinf / state%temperature) &
      * ( state%temperature / c_300 ) ** params%B_Kinf

  log10_term = log10(k_0 * state%number_density_air / k_inf)
  power_of_F =  c_1 + log10_term * log10_term

  rate_constant = k_0 / ( c_1 + k_0 * state%number_density_air / k_inf) &
      *  params%F_C ** ( c_1/power_of_F)

end function Troe_chemical_activation

function simple_aerosol_heterogeneous_rate(params, state) result (rate_constant)
  type(simple_aerosol_heterogeneous_rate_rate_param_type), intent(in) :: params
  type(environmental_state_type), intent(in):: state
  real :: rate_constant

  real :: gas_velocity  ! cm/sec

  gas_velocity = sqrt(8. * boltzman_cgs * avogadro * state%temperature &
     / (c_pi * params%m_g) )

  rate_constant = sum( &
    state%aerosol_surface_area_density(:) &
    / ( 0.5*state%aerosol_diameter(:) / params%Dg  &
      + 4.0/(params%gamma * gas_velocity) &
      ) &
  )

end function simple_aerosol_heterogeneous_rate

function HET_glyoxyl(params, state) result(rate_constant)
  type(HET_glyoxyl_rate_param_type), intent(in) :: params
  type(environmental_state_type), intent(in):: state
  real :: rate_constant

  real :: gas_velocity ! cm/sec

  ! first order uptake, Fuchs and Sutugin, 1971,  dCg = 1/4 * gamma * ! A * |v_mol| * Cg * dt
  ! 19267.9 = sqrt(8 boltzman_cgs * avogadro / pi / HO2_mass)

  ! mysterious special constant
  gas_velocity = 19267.9*sqrt(state%temperature )

  rate_constant = params%gamma * gas_velocity * sum( state%aerosol_surface_area_density(:) ) / 4.

end function

function combined_CO_OH(params, state) result (rate_constant)
  type( combined_CO_OH_rate_param_type ), intent(in) :: params
  type(environmental_state_type), intent(in) :: state
  real :: rate_constant


  rate_constant = params%A * &
      ( c_1 + params%B * boltzman_cgs * avogadro * state%number_density_air * state%temperature)

end function combined_CO_OH

function CH3COCH3_OH(params, state) result(rate_constant)
  type(CH3COCH3_OH_rate_param_type), intent(in) :: params
  type(environmental_state_type), intent(in):: state
  real :: rate_constant

  rate_constant = params%A + &
      params%B * exp(-params%C /state%temperature)

end function CH3COCH3_OH


function HO2_HO2(params, state) result (rate_constant)
  type(HO2_HO2_rate_param_type), intent(in) :: params
  type(environmental_state_type), intent(in) :: state
  real :: rate_constant

  real :: k0, kinf, fc

  k0 = params%A_k0 * exp(-params%C_k0 / state%temperature)
  kinf = params%A_kinf * exp(-params%C_kinf / state%temperature)

  fc =  c_1 + params%F * state%number_density_air * state%h2ovmr * exp( -params%C_kinf / state%temperature)
  rate_constant = (k0 + kinf) * fc

end function HO2_HO2

function DMS_OH(params, state) result (rate_constant)
  type(DMS_OH_rate_param_type), intent(in) :: params
  type(environmental_state_type), intent(in) :: state
  real :: rate_constant

  real :: k0, k1

  k0 = c_1 + params%A*exp(-params%B/state%temperature)*state%o2_number_density
  k1 = params%C*exp(-params%D/state%temperature)*state%o2_number_density
  rate_constant = k1/k0
  

end function DMS_OH

function HNO3_OH(params, state) result (rate_constant)
  type(HNO3_OH_rate_param_type), intent(in) :: params
  type(environmental_state_type), intent(in) :: state
  real :: rate_constant

  real :: k0, k1, k2

  k0 = params%A_k0 * exp(-params%C_k0 / state%temperature)
  k1 = params%A_k1 * exp(-params%C_k1 / state%temperature)
  k2 = params%A_k2 * exp(-params%C_k2 / state%temperature)

  rate_constant = k2 + state%number_density_air*k0 / &
    ( c_1 + state%number_density_air*k0 / k2)

end function HNO3_OH


function SO2_OH(params, state) result (rate_constant)
  type(SO2_OH_rate_param_type), intent(in) :: params
  type(environmental_state_type), intent(in) :: state
  real :: rate_constant

  real  f_c, k0, power_of_F, log_term_in_power

  f_c = params%A * (300./state%temperature)**params%B
  k0 = f_c*state%number_density_air / (c_1 + f_c * state%number_density_air/params%C)

  log_term_in_power = log10(f_c*state%number_density_air)-log10(params%C)
  power_of_F = c_1 / ( c_1 +  log_term_in_power*log_term_in_power)

  rate_constant = k0*params%F**power_of_F

end function SO2_OH

function MCO3_NO2(params, state) result (rate_constant)
  type(MCO3_NO2_rate_param_type), intent(in) :: params
  type(environmental_state_type), intent(in) :: state
  real :: rate_constant

  rate_constant = params%A / state%number_density_air * 300. / state%temperature

end function MCO3_NO2


function MPAN_M(params, state) result (rate_constant)
  type(MPAN_M_rate_param_type), intent(in) :: params
  type(environmental_state_type), intent(in) :: state
  real :: rate_constant

  real :: k_f, k_r

  k_f = params%A_f / state%number_density_air * 300. / state%temperature
  k_r = params%A_r * exp(-params%B_r / state%temperature)
 
  rate_constant = k_r * k_f

end function MPAN_M


end module rate_constant_functions
