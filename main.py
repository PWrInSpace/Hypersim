#from types import DynamicClassAttribute
from nitrous_oxide import *
from init import *

plt.style.use('ggplot')

"""
TODO:
    - dysza
    - faza gazowa
    - oznaczenia nasze -> Aspirespace
    - wszystko co po '# @' usuwamy przed publikacją
"""


# ---------------- LOOP ------------------

while fuel_mass > 0 and liquid_mass > 0:
    #basic calculations
    port_area = math.pi/4 * port_diameter**2    # diameter -> area
    pt_ratio = port_area/throat_area            # port to throat area ratio
    if pt_ratio < config['min_pt']:     # sanity check
        warnings.warn("port to throat area ratio is too large (%s), increase port diameter!" %round(pt_ratio, 2))
        FAULT.append('pt_ratio') 

    
    lagged_bob = 0
    # @ old_ox_flow = Omdot_tank_outflow
    # @ ox_flow = mdot_tank_outflow
    old_ox_flow = ox_flow 
    # @ enth_v -> Enth_of_vap  
    enth_v = float(nox_enth_v(tank_temp))         # entalphy of vapourisation
    # @ heat_capacity -> Spec_heat_cap
    heat_capacity = nox_l_Cp(tank_temp)    # liquid heat capacity
    # @ heat_removed -> deltaQ
    # @ old_vapourised_mass -> vapourised_mass_old
    heat_removed = old_vapourised_mass*enth_v
    # @ temp_drop -> deltaTemp
    # @ liquid mass -> hybrid.tank_liquid_mass
    temp_drop = -(heat_removed/(liquid_mass*heat_capacity))
    # @ tank_temp hybrid.tank_fluid_temperature_K
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
    # @ liquid_density -> hybrid.tank_liquid_density
    liquid_density = nox_l_rho(tank_temp)
    # @ vapour_density -> hybrid.tank_vapour_density
    vapour_density = nox_v_rho(tank_temp)
    # @ tank pressure -> hybrid.tank_pressure_bar
    tank_pressure = float(nox_vp(tank_temp))
    # @ pressure_drop -> pressure_drop
    pressure_drop = tank_pressure - head_pressure

    if pressure_drop < 0.000001:    # reality check
        pressure_drop = 0.000001
        warnings.warn("Backpressure: run tank pressure %s %s | head pressure: %s %s" %(round(Unit(tank_pressure, 'pressure'),2), units['pressure'], round(Unit(head_pressure, 'pressure'),2), units['pressure']))
        FAULT.append('backpressure')

    if pressure_drop/chamber_pressure < 0.2:    # safety check
        warnings.warn("Pressure drop is too low, backpressure is likely to happen: run tank pressure %s %s | head pressure: %s %s" %(round(Unit(tank_pressure, 'pressure'),2), units['pressure'], round(Unit(head_pressure, 'pressure'),2), units['pressure']))
        FAULT.append('pressure_drop')
        #break

    # @ d_loss -> hybrid.injector_loss_coefficent
    ox_flow = math.sqrt(2*liquid_density*pressure_drop/d_loss)
    # @ delta_ox -> delta_outflow_mass
    # @ config['time_step'] -> delta_time
    # @ ox_flow -> mdot_tank_outflow
    # @ old_ox_flow -> Omdot_tank_outflaw
    delta_ox = 0.5*config['time_step']*(3*ox_flow-old_ox_flow)  # 2nd order adams integration
    # @ ox_mass -> hybrid.tank_propellant_cintents_mass
    ox_mass -= delta_ox
    # @ old_liquid_mass -> old_liquid_nox_mass
    old_liquid_mass -= delta_ox
    bob = (1/liquid_density) - (1/vapour_density)
    # @ liquid_mass -> hybrid.tank_liquid_mass
    # @ tank['volume'] -> hybrid.tank_volume
    liquid_mass = (tank['volume']-(ox_mass/vapour_density))/bob
    # @ vapour_mass -> hybrid.tank_vapour_mass
    vapour_mass = ox_mass-liquid_mass
    bob = old_liquid_mass-liquid_mass
    tc = config['time_step']/0.15
    lagged_bob = tc*(bob-lagged_bob) + lagged_bob
    # @ old_vapourised_mass -> vapourised_mass_old
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

    if pressure_change < 100000:
    #if True:
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