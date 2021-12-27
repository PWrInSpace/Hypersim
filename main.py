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
    #temperature
    "c": -273.15,
    "k": 0,
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
    'tank_pressure': [],
    'tank_temp': [],
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

# functions
def Unit(value, type):
    if type == 'temperature':
        return unit[units[type]]+value
    else:
        return unit[units[type]]*value

def sgn(x):
    if x >= 0:
        return 1
    else:
        return -1

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
                shona += b[i] * math.pow(rab,p[i])
            pp_guess = p_crit * math.exp(shona/Tr)
            if (pp_guess - pressure) * sgn(step) >= 0:
                break
        step = step/(-2)
        if abs(pp_guess - pressure) <= 0.01:
            break
    return T
nox_temp = np.vectorize(_nox_temp)

#----------- INIT ---------------

# load model data
with open("data.json", 'r') as f:
    data = json.load(f)
grain = data['grain']
propellant = data['propellant']
nozzle = data['nozzle']
tank = data['tank']
injector = data['injector']
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
head_pressure = chamber_pressure
regression_rate = 0
port_diameter = grain['port']
fuel_flow = 0
ox_flow = 0
of = 0
k = propellant['k']
mod_heat_ratio = math.sqrt(k * math.pow(2/(k+1), (k+1)/(k-1)))
stuck = 0
FAULT = []

# initialize oxidizer tank
""" TODO:
        - choose initial temperature or pressure
        - choose initial ullage, or oxidizer mass"""
ox_flow = 0
tank_pressure = 50*100000
tank_temp = float(nox_temp(tank_pressure))

# reality check
if tank_temp > t_crit-0.1:
    tank_temp = t_crit-0.1
    FAULT.append('supercritical')

# calculate initial oxidizer volume/density/mass
liquid_density = nox_l_rho(tank_temp)
vapour_density = nox_v_rho(tank_temp)
vapour_volume = tank['ullage'] * tank['volume']
liquid_volume = tank['volume'] - vapour_volume
liquid_mass = liquid_volume*liquid_density
vapour_mass = vapour_volume*vapour_density
ox_mass = liquid_mass+vapour_mass

# initialize some values for later use
old_liquid_mass = liquid_mass
old_vapour_mass = vapour_mass
old_vapourised_mass = 0.001 

# calculate injector
injector_area = injector['number'] * math.pi/4*injector['orfice_diameter']**2
d_loss = injector['loss_coef']/injector_area**2


#temporary
z = 100*100000
exhoust_pressure = config['ambient_pressure']
pressure_drop = tank_pressure-chamber_pressure
# ---------------- LOOP ------------------



while fuel_mass > 0:
    #basic calculations
    port_area = math.pi/4 * port_diameter**2    # diameter -> area
    pt_ratio = port_area/throat_area            # port to throat area ratio
    if pt_ratio < config['min_pt']:     # sanity check
        warnings.warn("port to throat area ratio is too large (%s), increase port diameter!" %round(pt_ratio, 2))
        FAULT.append('pt_ratio') 

    # # oxidizer calculations
    # w = config['ox_flow_lag']
    # ox_flow = w*ox_flow + (1-w)*math.sqrt((pressure_drop) / z)
    # ox_flux = ox_flow/port_area
    # if ox_flux > config['max_ox_flux']: # sanity check
    #     warnings.warn("Oxidizer flux is too large! oxidizer mass flux: %s" %ox_flux)
    #     FAULT.append('max_ox_flux')

    
    lagged_bob = 0
    old_ox_flow = ox_flow       
    enth_v = float(nox_enth_v(tank_temp))         # entalphy of vapourisation
    heat_capacity = nox_l_Cp(tank_temp)    # liquid heat capacity
    heat_removed = old_vapourised_mass*enth_v
    temp_drop = -(heat_removed/(liquid_mass*heat_capacity))
    tank_temp += temp_drop

    # reality check
    if tank_temp < 183:
        warnings.warn("tank temperature too low! %s %s" %(round(tank_temp-273.15, 2), "C"))
        tank_temp = 183
        FAULT.append['low_ox_temp']
    if tank_temp > 309:
        warnings.warn("nitrous supercritical! %s %s" %(round(tank_temp-273.15, 2), "C"))
        tank_temp = 309
        FAULT.append['supercritical']
    liquid_density = nox_l_rho(tank_temp)
    vapour_density = nox_v_rho(tank_temp)
    tank_pressure = float(nox_vp(tank_temp))
    pressure_drop = tank_pressure - head_pressure

    if pressure_drop < 0.000001:    # reality check
        pressure_drop = 0.000001
        warnings.warn("Backpressure: run tank pressure %s %s | head pressure: %s %s" %(round(Unit(tank_pressure, 'pressure'),2), units['pressure'], round(Unit(head_pressure, 'pressure'),2), units['pressure']))
        FAULT.append('backpressure')

    if pressure_drop/chamber_pressure < 0.2:    # safety check
        warnings.warn("Pressure drop is too low, backpressure is likely to happen: run tank pressure %s %s | head pressure: %s %s" %(round(Unit(tank_pressure, 'pressure'),2), units['pressure'], round(Unit(head_pressure, 'pressure'),2), units['pressure']))
        FAULT.append('pressure_drop')
        #break

    ox_flow = math.sqrt(2*liquid_density*pressure_drop/d_loss)
    delta_ox = 0.5*config['time_step']*(3*ox_flow-old_ox_flow)  # 2nd order adams integration
    ox_mass -= delta_ox
    old_liquid_mass -= delta_ox
    bob = (1/liquid_density) - (1/vapour_density)
    liquid_mass = (tank['volume']-(ox_mass/vapour_density))/bob
    vapour_mass = ox_mass-liquid_mass
    bob = old_liquid_mass-liquid_mass
    tc = config['time_step']/0.15
    lagged_bob = tc*(bob-lagged_bob) + lagged_bob
    old_vapourised_mass = lagged_bob
    if liquid_mass > old_liquid_mass:
        warnings.warn("no more liquid nitrous")
        break
    old_liquid_mass = liquid_mass
    ox_flux = ox_flow/port_area

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
        FAULT.append('max_flux')

    # pressure calculations
    w = config['pressure_lag']
    old_chamber_pressure = chamber_pressure
    chamber_pressure = (w*chamber_pressure) + (1-w)*(mass_flow*propellant['cstar']/throat_area)
    head_pressure = chamber_pressure * (1 + 0.5 * (1/pt_ratio)**2 * mod_heat_ratio**2)
    
    # temporary fix, TODO calculate exhoust pressure
    if chamber_pressure < exhoust_pressure:
        chamber_pressure = exhoust_pressure
    
    


    # performance calculations
    pressure_thrust = exit_area * (exhoust_pressure-config['ambient_pressure'])
    cf = config['cf_eff'] * mod_heat_ratio * math.sqrt(math.pow(2/(k+1),(k+1)/(k-1)) * (1- math.pow(exhoust_pressure/chamber_pressure, (k-1)/k)))
    cf = pressure_thrust/throat_area/chamber_pressure
    thrust = mass_flow * propellant['cstar'] * cf
    isp = thrust/mass_flow/9.81

    pressure_change = (chamber_pressure-old_chamber_pressure)/config['time_step']

    if(abs(pressure_change) > 100000):
        stuck += 1
        if stuck > 1000:
            warnings.warn("Error, infinite loop!")
            FAULT.append('stuck')
            break
    else:
        stuck = 0
    #print(pressure_change)

    #if pressure_change < 100000:
    if True:
        # time step
        port_diameter += regression_rate*config['time_step']    # increase port diameter
        fuel_mass -= fuel_flow*config['time_step']  # burn fuel
        ox_mass -= ox_flow*config['time_step']  # burn oxidizer

        # save for plotting
        plots['chamber_pressure'].append(chamber_pressure)
        plots['tank_pressure'].append(tank_pressure)
        plots['tank_temp'].append(tank_temp)
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
        
# print errors/warnings
printed = [None]
exists = False
for i in FAULT:
    for i2 in printed:
        if i == i2:
            exists = True
            continue
    if exists == False:
        printed.append(i)
        print("FAULT: ", i)
    else:
        exists = False

# unit change
for element in plot:
    for i in range(len(plots[element])):
        plots[element][i] = Unit(plots[element][i], plot[element]['type'])

ts = np.arange(0, len(plots['port_diameter'])*config['time_step'], config['time_step'])

#plot
for element in plot:
    if(plot[element]['show']):
        plt.plot(ts, plots[element], label = plot[element]['label'] + " - " + units[plot[element]['type']])
plt.legend()
plt.show()