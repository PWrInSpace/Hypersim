import numpy as np
import math
import json
import warnings
import matplotlib.pyplot as plt
from nitrous_oxide import *

# load user preferences
with open('preferences.json', 'r') as f:
    data = json.load(f)
units = data['units']
plot = data['plot']

# load model data
with open("data.json", 'r') as f:
    data = json.load(f)
grain = data['grain']
propellant = data['propellant']
nozzle = data['nozzle']
tank = data['tank']
injector = data['injector']
config = data['config']

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

    """Tranforms given value to ISO units
    Aruments:
        - value - value of the measure
        - type - unit of the value
    Return:
        - value in ISO units"""

    if type == 'temperature':
        return unit[units[type]]+value
    else:
        return unit[units[type]]*value


def integrate_mass_flowrate(ox_flow, old_ox_flow):
    """Intregration of oxidizer mass flowrate using 2nd order Adam's integration formula
    Arguments:
        - ox_flow - oxidier flow
        - old_ox_flow - old oxodizer flow
    Return:
        - delta_ox - numerical approximation of the integral"""

    delta_ox = 0.5*config['time_step']*(3*ox_flow-old_ox_flow)
    return delta_ox


def linear_interpolation(x, x1, y1, x2, y2):
    """Extrapolate the y value for postition x on 
    line created by points x1, y1, x2, y2.
    Arguments:
        - x - point to extrapolate the value of
        - x1 - minimum x bound
        - y1 - value of x1
        - x2 - maximum x bound
        - y2 - value of x2
    Return:
        - y - the extrapolated value"""

    if x1 < x2 and (x <= x1 or x >= x2):
        if x <= x1:
            return y1
        else:
            return y2
    
    elif x1 > x2 and (x >= x1 or x <= x2):
        if x >= x1:
            return y1
        else:
            return y2
    
    else:
        a = (y2 - y1)/(x2 - x1)
        b = y1 - a * x1
        y = a * x + b
        return y


def compress_factor(p, p_crit, z_crit):
    """Calculate compressibility factor of subctritical vapour on 
    the saturation line.
    Arguments:
        - p - pressure
        - p_crit - critical pressure
        - z_crit - critical compressibility factor
    Return:
        - z - compressibility factor"""

    p = abs(p)
    z = linear_interpolation(p, 0, 1, p_crit, z_crit)

    return z

def tank_temp_reality_check(tank_temp, FAULT):
    """Check if the temperature is at reasonable value"""
    if tank_temp < 183:
        warnings.warn("tank temperature too low! %s %s" %(round(tank_temp-273.15, 2), "C"))
        tank_temp = 183
        FAULT.append['low_ox_temp']
    if tank_temp > 309:
        warnings.warn("nitrous supercritical! %s %s" %(round(tank_temp-273.15, 2), "C"))
        tank_temp = 309
        FAULT.append['supercritical']

