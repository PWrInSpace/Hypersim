from funs import *

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
