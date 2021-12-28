import numpy as np

# --------- Numeric proximiation of nitrous properties -----------------
p_crit = 72.51*100000      #critical pressure [Pa]
rho_crit = 452.0    #critical density [kg/m^3]
t_crit = 309.57     #critical temperature [Kelvin]
z_crit = 0.28       #critical compressibility factor
gamma = 1.3         #average over subctritical range

def _nox_vp(T_Kelvin):
    """Calculate vapour pressure of nitrous oxide in Pa
    based on experimental data
    Arguments:
        T_Kelvin - temperature in Kelvin
    Return:
        bob - estimated vapour pressure"""
    p = [1, 1.5, 2.5, 5]
    b = [-6.71893, 1.35966, -1.3779, -4.051]

    Tr = T_Kelvin/309.57
    rab = 1 - Tr
    shona = 0

    for i in range(4):
        shona += b[i] * np.power(rab, p[i])

    bob = 72.51 * np.exp(shona/Tr) * 100000 #bar to Pa
    return bob

nox_vp = np.vectorize(_nox_vp)


def _nox_l_rho(T_Kelvin):
    """Calculate saturated liquid density of nitrous oxide
    based on experimental data
    Arguments:
        T_Kelvin - temperature in Kelvin
    Return:
        bob - estimated density"""
    b = [1.72328, -0.8395, 0.5106, -0.10412]
    Tr = T_Kelvin / t_crit
    rab = 1.0 - Tr
    shona = 0.0
    for i in range(4):
        shona += b[i] * np.power(rab, ((i + 1) / 3.0))
    bob = rho_crit * np.exp(shona)
    return bob

nox_l_rho = np.vectorize(_nox_l_rho)


def _nox_v_rho(T_Kelvin):
    """Calculate saturated vapour density of nitrous oxide
    based on experimental data
    Arguments:
        T_Kelvin - temperature in Kelvin
    Return:
        bob - estimated density"""
    b = [-1.009, -6.28792, 7.50332, -7.90463, 0.629427]
    Tr = T_Kelvin / t_crit
    rab = (1.0 / Tr) - 1.0
    shona = 0.0
    for i in range(5):
        shona += b[i] * np.power(rab, ((i+1) / 3.0))
        bob = rho_crit * np.exp(shona)
    return bob


nox_v_rho = np.vectorize(_nox_v_rho)


def _nox_enth_v(T_Kelvin):
    """Calculate nitrous liquid enthalpy (latent heat) of vaporisation, J/kg
    Arguments:
        T_Kelvin - temperature in Kelvin
    Return:
        bob - enthalpy of vaporisation"""
    b_l = [-200.0, 116.043, -917.225, 794.779, -589.587]
    b_v = [-200.0, 440.055, -459.701, 434.081, -485.338]
    Tr = T_Kelvin / t_crit
    rab = 1.0 - Tr
    shona_l = b_l[0]
    shona_v = b_v[0]

    for i in range(5):
        shona_l += b_l[i] * np.power(rab, (i / 3.0))  # saturated liquid enthalpy
        shona_v += b_v[i] * np.power(rab, (i / 3.0))  # saturated vapour enthalpy

    bob = (shona_v - shona_l) * 1000.0  #net during change from liquid to vapour
    return bob


nox_enth_v = np.vectorize(_nox_enth_v)


def _nox_l_Cp(T_Kelvin):
    """Calculate Nitrous saturated liquid isobaric heat capacity, J/kg K
    Arguments:
        T_Kelvin - temperature in Kelvin
    Return:
        bob - enthalpy of vaporisation"""
    b = [2.49973, 0.023454, -3.80136, 13.0945, -14.518]
    Tr = T_Kelvin / t_crit
    rab = 1.0 - Tr
    shona = 1.0 + b[1] / rab
    for i in range(4):
        shona += b[(i+1)] * rab**i
    bob = b[0] * shona * 1000.0; # convert from KJ to J 
    return bob
nox_l_Cp = np.vectorize(_nox_l_Cp)

def _nox_temp(pressure):
    """Calculate saturated nitrous temperature based on vapour pressure"""
    p = [1, 1.5, 2.5, 5]
    b = [-6.71803, 1.35966, -1.3779, -4.051]

    step = -1

    T = (t_crit-0.1)-step
    while True:
        while True:
            T += step
            Tr = T/t_crit
            rab = 1-Tr
            shona = 0
            for i in range(4):
                shona += b[i] * np.power(rab,p[i])
            pp_guess = p_crit * np.exp(shona/Tr)
            if (pp_guess - pressure) * np.sign(step) >= 0:
                break
        step = step/(-2)
        if abs(pp_guess - pressure) <= 0.01:
            break
    return T
    
nox_temp = np.vectorize(_nox_temp)