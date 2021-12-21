#from types import DynamicClassAttribute
import warnings
import math
import matplotlib.pyplot as plt
import json
import numpy as np

plt.style.use('ggplot')

# unit conversion
unit = {
    #length
    "mm": 1000,
    "cm": 100,
    "m": 1,
    "in": 39.37,
    "ft": 3.281,
    #volume
    "mm^3": 1000000000,
    "cm^3": 1000000,
    "m^3": 1,
    "ml": 1000000,
    "l": 1000,
    "in^3": 61024,
    "ft^3": 35.315,
    "gal": 264,
    #area
    "mm^2": 1000000,
    "cm^2": 10000,
    "m^2": 1,
    "in^2": 1550,
    "ft^2": 10.764,
    #velocity
    "mm/s": 1000,
    "cm/s": 100,
    "m/s": 1,
    "in/s": 39.37,
    "ft/s": 3.281,
    #force
    "N": 1,
    "kgf": 1/9.81,
    "lbf": 1/4.448,
    #impulse
    "Ns": 1,
    "kgfs": 1/9.81,
    "lbfs": 1/4.448,
    #pressure 
    "Pa": 1,
    "hPa": 1/100,
    "kPa": 1/1000,
    "MPa": 1/1000000,
    "atm": 1/101325,
    "bar": 1/100000,
    "psi": 1/6895,
    #mass
    "g": 1000,
    "kg": 1,
    "lb": 2.205,
    "oz": 35.274,
    #density
    "g/cm^3": 1/1000,
    "kg/m^3": 1,
    "lb/in^3": 1/27680,
    #flow
    "g/s": 1000,
    "kg/s": 1,
    "lb/s": 2.205,
    #flux
    "kg/m^2*s": 1,
    "lb/in^2*s": 1/703,
    #regression coefficeint
    "m*m^2/kg": 1,
    "mm*m^2/kg": 1000,
    "in*in^2/lb": 27680,
    # time
    "s": 1,
    # -
    "-": 1
}

# data storage for plotting
plots = {
    'chamber_pressure': [],
    'thrust': [],
    'ox_mass': [],
    'fuel_mass': [],
    'ox_flow': [],
    'fuel_flow': [],
    'mass_flow': [],
    'ox_flux': [],
    'fuel_flux': [],
    'mass_flux': [],
    'regression_rate': [],
    'port_diameter': [],
    'exhoust_pressure': [],
    'cf': [],
    'pressure_thrust': [],
    'of': [],
    'isp': []
}

# --------- Numeric proximiation of nitrous properties -----------------
p_crit = 72.51      #critical pressure [Bar]
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
        shona += b[i] * pow(rab, ((i+1) / 3.0))
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
        shona_l += b_l[i] * pow(rab, (i / 3.0))  # saturated liquid enthalpy
        shona_v += b_v[i] * pow(rab, (i / 3.0))  # saturated vapour enthalpy

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

# functions
def Unit(value, type):
    return unit[units[type]]*value

#----------- INIT ---------------

# load model data
with open("data.json", 'r') as f:
    data = json.load(f)
grain = data['grain']
propellant = data['propellant']
nozzle = data['nozzle']
tank = data ['tank']
config = data['config']

# load user preferences
with open('preferences.json', 'r') as f:
    data = json.load(f)
units = data['units']
plot = data['plot']

# calculate initial values
throat_area = math.pi/4 * nozzle['throat']**2
exit_area = math.pi/4 * nozzle['exit']**2
expansion_ration = exit_area/throat_area
fuel_mass = propellant['density'] * grain['length'] * math.pi/4 * (grain['diameter']**2 - grain['port']**2)

# initialize variables
chamber_pressure = config['ambient_pressure']
regression_rate = 0
port_diameter = grain['port']
ox_mass = tank['mass']
chamber_pressure = 0
head_pressure = 0
fuel_flow = 0
ox_flow = 0
of = 0
k = propellant['k']
mod_heat_ratio = math.sqrt(k * math.pow(2/(k+1), (k+1)/(k-1)))

#temporary
z = 100*100000
tank_pressure = 50*100000
exhoust_pressure = config['ambient_pressure']

# ---------------- LOOP ------------------

while fuel_mass > 0:
    # time step'
    port_diameter += regression_rate*config['time_step']    # increase port diameter
    fuel_mass -= fuel_flow*config['time_step']  # burn fuel
    ox_mass -= ox_flow*config['time_step']  # burn oxidizer

    #basic calculations
    port_area = math.pi/4 * port_diameter**2    # diameter -> area
    pt_ratio = port_area/throat_area            # port to throat area ratio
    if pt_ratio < config['min_pt']:     # sanity check
        warnings.warn("port to throat area ratio is too large (%s), increase port diameter!" %round(pt_ratio, 2))
    
    # oxidizer calculations
    w = config['ox_flow_lag']
    ox_flow = w*ox_flow + (1-w)*math.sqrt((tank_pressure-head_pressure) / z)
    ox_flux = ox_flow/port_area
    if ox_flux > config['max_ox_flux']: # sanity check
        warnings.warn("Oxidizer flux is too large! oxidizer mass flux: %s" %ox_flux)
    
    # fuel calculations
    regression_rate = propellant['a'] * math.pow(ox_flux, propellant['n'])
    w = config['fuel_flow_lag']
    fuel_flow = w*fuel_flow + (1-w)* propellant['density'] * regression_rate * grain['length'] * math.pi * port_diameter
    fuel_flux = fuel_flow/port_area

    # propellant calculations
    if fuel_flow != 0: # sanity check
        of = ox_flow/fuel_flow          # calculate O/F
    mass_flow = fuel_flow+ox_flow       # total mass flow rate
    mass_flux = mass_flow/port_area     # total mass flux
    if mass_flux > config['max_flux']:  # sanity check
        warnings.warn("propellant mass flux is too large! (%s)" %mass_flux)

    # pressure calculations
    w = config['pressure_lag']
    chamber_pressure = (w*chamber_pressure) + (1-w)*(mass_flow*propellant['cstar']/throat_area)
    head_pressure = chamber_pressure * (1 + 0.5 * (1/pt_ratio)**2 * mod_heat_ratio**2)
    if head_pressure > tank_pressure:   # sanity check
        #head_pressure = tank_pressure   # quick correction so the program can go on
        warnings.warn("Backpressure! Tank pressure: %s %s | head pressure: %s %s" %(round(Unit(tank_pressure, 'pressure'),2), units['pressure'], round(Unit(head_pressure, 'pressure'),2), units['pressure']))
        break

    # performance calculations
    pressure_thrust = expansion_ration * (exhoust_pressure-config['ambient_pressure'])/chamber_pressure
    cf = config['cf_eff'] * mod_heat_ratio * math.sqrt(math.pow(2/(k+1),(k+1)/(k-1)) * (1- math.pow(exhoust_pressure/chamber_pressure, (k-1)/k)))
    thrust = mass_flow * propellant['cstar'] * cf
    isp = thrust/mass_flow/9.81

    # save for plotting
    plots['chamber_pressure'].append(chamber_pressure)
    plots['thrust'].append(thrust)
    plots['ox_mass'].append(ox_mass)
    plots['fuel_mass'].append(fuel_mass)
    plots['ox_flow'].append(ox_flow)
    plots['fuel_flow'].append(fuel_flow)
    plots['mass_flow'].append(mass_flow)
    plots['ox_flux'].append(ox_flux)
    plots['fuel_flux'].append(fuel_flux)
    plots['mass_flux'].append(mass_flux)
    plots['regression_rate'].append(regression_rate)
    plots['port_diameter'].append(port_diameter)
    plots['exhoust_pressure'].append(exhoust_pressure)
    plots['cf'].append(cf)
    plots['pressure_thrust'].append(pressure_thrust)
    plots['of'].append(of)
    plots['isp'].append(isp)


# unit change
for element in plot:
    for i in range(len(plots[element])):
        plots[element][i] *= unit[units[plot[element]['type']]]

ts = np.arange(0, len(plots['port_diameter'])*config['time_step'], config['time_step'])


#plot
for element in plot:
    if(plot[element]['show']):
        plt.plot(ts, plots[element], label = plot[element]['label'] + " - " + units[plot[element]['type']])
plt.legend()
plt.show()