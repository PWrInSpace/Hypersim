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
