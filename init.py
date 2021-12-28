from numpy.core.numeric import array_equal
from funs import *


#----------- INIT ---------------

# calculate initial values

"""Calculate initial values.
    Arguments:
        - nozzle - nozzle data
        - grain - grain data
        - propellant - propellant data
    Return:
        - throat_area - area of the nozzle throat
        - exit_area - area of the nozzle exit
        - expansion_ration
        - fuel_mass - mass of the fuel"""

throat_area = np.pi/4 * nozzle['throat']**2
exit_area = np.pi/4 * nozzle['exit']**2
expansion_ration = exit_area/throat_area
fuel_mass = propellant['density'] * grain['length'] * np.pi/4 * (grain['diameter']**2 - grain['port']**2)



"""Initialize variables
    Arguments:
        - config - basic configuration data
        - chamber_pressure
        - grain - grain data
        - propellant - propellant data
    Return:
        - chamber_pressure - set to ambient_pressure for now
        - head_pressure - # ? pressure after the injector set to chamber_pressure
        - regression_rate
        - port_diameter
        - fuel_flow  - flow of the fuel # ? jednostki?
        - ox_flow - flow of the oxidizer
        - of - oxidizer to fuel ratio
        - k - constant associated with propellant (or gamma)
        - mod_heat_ratio - modified specific heat ratio
        - stuck - count of idle iterations
        - FAULT - list of warnings
    """

chamber_pressure = config['ambient_pressure']
head_pressure = chamber_pressure
regression_rate = 0
port_diameter = grain['port']
fuel_flow = 0
ox_flow = 0
of = 0
k = propellant['k']
mod_heat_ratio = np.sqrt(k * np.power(2/(k+1), (k+1)/(k-1)))
stuck = 0
FAULT = []

# initialize oxidizer tank
""" TODO:
        - choose initial temperature or pressure
        - choose initial ullage, or oxidizer mass"""
ox_flow = 0
tank_pressure = 50*100000
tank_temp = float(nox_temp(tank_pressure))

"""Reality check. Check if temperature in nitrous oxide
 is supercritical. """
if tank_temp > t_crit-0.1:
    tank_temp = t_crit-0.1
    FAULT.append('supercritical')

"""Calculate initial oxidizer volume/density/mass using 
numerical approximation of those values. 
    Return:
    - liquid_density
    - vapour_density
    - vapour volume
    - liquid_mass
    - vapour_mass
    - ox_mass - oxidizer mass"""
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
"""Calculate injector
    - injector_area - total area  of the injectior
    - d_loss  - loss coefficient divided by inector area to the power of 2"""
injector_area = injector['number'] * math.pi/4*injector['orfice_diameter']**2
d_loss = injector['loss_coef']/injector_area**2


#temporary
z = 100*100000
exhoust_pressure = config['ambient_pressure']
pressure_drop = tank_pressure-chamber_pressure
