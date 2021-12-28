#from types import DynamicClassAttribute

from init import *

plt.style.use('ggplot')

"""
TODO:
    - dysza
    - faza gazowa
    - oznaczenia nasze -> Aspirespace
    - wszystko co po '# @' usuwamy przed publikacją
"""

## Comparison of our and Aspirespace's variables
    # @ old_ox_flow = Omdot_tank_outflow
    # @ ox_flow = mdot_tank_outflow
    # @ enth_v -> Enth_of_vap 
    # @ heat_capacity -> Spec_heat_cap
    # @ heat_removed -> deltaQ
    # @ old_vapourised_mass -> vapourised_mass_old
    # @ temp_drop -> deltaTemp
    # @ liquid mass -> hybrid.tank_liquid_mass
    # @ tank_temp hybrid.tank_fluid_temperature_K
    # @ liquid_density -> hybrid.tank_liquid_density
    # @ vapour_density -> hybrid.tank_vapour_density
    # @ tank pressure -> hybrid.tank_pressure_bar
    # @ d_loss -> hybrid.injector_loss_coefficent
    # @ delta_ox -> delta_outflow_mass
    # @ config['time_step'] -> delta_time
    # @ ox_flow -> mdot_tank_outflow
    # @ old_ox_flow -> Omdot_tank_outflaw
    # @ ox_mass -> hybrid.tank_propellant_cintents_mass
    # @ old_liquid_mass -> old_liquid_nox_mass
    # @ liquid_mass -> hybrid.tank_liquid_mass
    # @ tank['volume'] -> hybrid.tank_volume
    # @ vapour_mass -> hybrid.tank_vapour_mass
    # @ old_vapourised_mass -> vapourised_mass_old
    # @ p - P_bar_abs


## FFunctions - necessary to be in the main - inherit from main loop

def injector_flow(tank_pressure, head_pressure, density=liquid_density):
    """Calculate mass flow rate out of the injector and perform a reality check and safety check.
    Arguments:
        - tank_pressure - pressure inside the tank
        - head_pressure - pressure inside comubstion chamber immidietly after injector
    Return:
        - ox_flow - liquid flowrate through the injector based 
        on the pressure drop and d_loss between the run-tank and combustion chamber
        """

    pressure_drop = tank_pressure - head_pressure

    if pressure_drop < 0.000001:    # reality check
        pressure_drop = 0.000001
        warnings.warn("Backpressure: run tank pressure %s %s | head pressure: %s %s" %(round(Unit(tank_pressure, 'pressure'),2), units['pressure'], round(Unit(head_pressure, 'pressure'),2), units['pressure']))
        FAULT.append('backpressure')

    if pressure_drop/chamber_pressure < 0.2:    # safety check
        warnings.warn("Pressure drop is too low, backpressure is likely to happen: run tank pressure %s %s | head pressure: %s %s" %(round(Unit(tank_pressure, 'pressure'),2), units['pressure'], round(Unit(head_pressure, 'pressure'),2), units['pressure']))
        FAULT.append('pressure_drop')
        #break

    ox_flow = math.sqrt(2*density*pressure_drop/d_loss)

    return ox_flow


# ---------------- LOOP ------------------

while fuel_mass > 0 and liquid_mass > 0:
    #basic calculations

    """Basic calculations and check if pt_ratio is not too large
        - port area - area of the grain's port
        - pt_ratio - port throat area ratio
        """
    port_area = math.pi/4 * port_diameter**2    # diameter -> area
    pt_ratio = port_area/throat_area            # port to throat area ratio
    if pt_ratio < config['min_pt']:     # sanity check
        warnings.warn("port to throat area ratio is too large (%s), increase port diameter!" %round(pt_ratio, 2))
        FAULT.append('pt_ratio') 

    
    lagged_tmp = 0
    """Init:
            - old_ox_flow - old oxidizer flow as current oxidizer flow
            - enth_v - current enthalpy of nitrous vapour
            - heat_capacity - current heat capacity of liquid nitrous
            - heat_reamoved - heat removed in nitrous vaporisation
            - temp_drop - drop in temperature caused by removed heat
            - tank_temp - update after temp_drop"""


    old_ox_flow = ox_flow 
     
    enth_v = float(nox_enth_v(tank_temp))         # entalphy of vapourisation
    heat_capacity = nox_l_Cp(tank_temp)    # liquid heat capacity
    heat_removed = old_vapourised_mass*enth_v
    temp_drop = -(heat_removed/(liquid_mass*heat_capacity))
    tank_temp += temp_drop

    # reality check
    """Check if the temperature is at reasonable value"""
    if tank_temp < 183:
        warnings.warn("tank temperature too low! %s %s" %(round(tank_temp-273.15, 2), "C"))
        tank_temp = 183
        FAULT.append['low_ox_temp']
    if tank_temp > 309:
        warnings.warn("nitrous supercritical! %s %s" %(round(tank_temp-273.15, 2), "C"))
        tank_temp = 309
        FAULT.append['supercritical']

    """Update nitrous properties after temperature changed"""

    liquid_density = nox_l_rho(tank_temp)
    vapour_density = nox_v_rho(tank_temp)
    tank_pressure = float(nox_vp(tank_temp))


    """Calculate
        - delta_ox  - integral of the difference between old and new ox_flow calculated using 2nd 
        order adams integration.
        - ox_mass - update by extracting delta_ox
        - liquid mass - current mass of the liquid nitrous
        - vapour_mass - update vapour mass"""


    ox_flow = injector_flow(tank_pressure, head_pressure)
    delta_ox = integrate_mass_flowrate(ox_flow, old_ox_flow)  # 2nd order adams integration
    ox_mass -= delta_ox
    old_liquid_mass -= delta_ox
    tmp = (1/liquid_density) - (1/vapour_density)
    liquid_mass = (tank['volume']-(ox_mass/vapour_density))/tmp
    vapour_mass = ox_mass-liquid_mass  
    tmp = old_liquid_mass-liquid_mass
    tc = config['time_step']/0.15
    lagged_tmp = tc*(tmp-lagged_tmp) + lagged_tmp # first order time lag to help with numerical stability
    old_vapourised_mass = lagged_tmp

    """Check if there is liquid nitrous left"""
    if liquid_mass > old_liquid_mass:
        warnings.warn("no more liquid nitrous")
        break


    old_liquid_mass = liquid_mass
    ox_flux = ox_flow/port_area

    # fuel calculations
    """Fuel calculations:
        - regression_rate
        - fuel_flow
        - fluel_flux"""
    regression_rate = propellant['a'] * np.power(ox_flux, propellant['n'])
    fuel_flow = propellant['density'] * regression_rate * grain['length'] * np.pi * port_diameter
    fuel_flux = fuel_flow/port_area

    # propellant calculations
    """Propellant calculations:
        - of - oxodizer to fuel rate
        - mass_flow - total mass flow
        - mass_flux - total mass_flux"""

    if fuel_flow != 0: # sanity check
        of = ox_flow/fuel_flow          # calculate O/F
    mass_flow = fuel_flow+ox_flow       # total mass flow rate
    mass_flux = mass_flow/port_area     # total mass flux

    """Cehck if mass_flux does not exceed maximum value"""
    if mass_flux > config['max_flux']:  # sanity check
        warnings.warn("propellant mass flux is too large! (%s)" %mass_flux)
        FAULT.append('max_flux')

    # pressure calculations
    """Pressure calculations:
        - old_chamber_pressure - set to chamber pressure
        - chaber_pressure - new pressure calculated from mass_flow,
         a propellant constant and throat area
        - head_pressure - # ?"""
    old_chamber_pressure = chamber_pressure
    chamber_pressure = mass_flow*propellant['cstar']/throat_area
    head_pressure = chamber_pressure * (1 + 0.5 * (1/pt_ratio)**2 * mod_heat_ratio**2)
    
    # temporary fix, TODO calculate exhoust pressure
    if chamber_pressure < exhoust_pressure:
        chamber_pressure = exhoust_pressure
    
    


    # performance calculations
    """Calculate performance:
        - pressure_thrust
        - cf - efficiency 
        - thrust
        - isp - specific impulse"""
    pressure_thrust = exit_area * (exhoust_pressure-config['ambient_pressure'])
    cf = config['cf_eff'] * mod_heat_ratio * np.sqrt(np.power(2/(k+1),(k+1)/(k-1)) * (1- np.power(exhoust_pressure/chamber_pressure, (k-1)/k)))

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
        for key in plots.keys():
            plots[key].append(globals()[key])
        
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

#unit change
for element in plot:
    for i in range(len(plots[element])):
        plots[element][i] = Unit(plots[element][i], plot[element]['type'])

ts = np.arange(0, len(plots['port_diameter'])*config['time_step'], config['time_step'])


print(len(ts))
#plot
for element in plot:
    if(plot[element]['show']):
        plt.plot(ts, plots[element], label = plot[element]['label'] + " - " + units[plot[element]['type']])
plt.legend()
plt.show()


## Vapour faze

# initial conditions
init_vapour_temp = tank_temp
init_vapour_mass = vapour_mass
init_vapour_pressure = tank_pressure
init_vapour_density = vapour_density

init_Z = compress_factor(init_vapour_pressure, p_crit, z_crit)
old_ox_flow = 0 #reset
first = False

while vapour_mass > 0:

    ox_flow = injector_flow(tank_pressure, head_pressure, denisty=vapour_density)
    delta_ox = integrate_mass_flowrate(ox_flow, old_ox_flow)
    ox_mass -= delta_ox
    vapour_mass -= delta_ox

    #Guessing of Z

    current_z_guess = compress_factor(tank_pressure, p_crit, z_crit) #initial guess
    step = 1/0.9 #initial step size
    old_aim = 2
    aim = 0
    current_z = current_z_guess*2
    go = True

    while go:
        tmp = k - 1
        vapour_temp = init_vapour_temp * np.power(vapour_mass * current_z_guess\
            /(init_vapour_mass * init_Z), tmp)
        tmp = k/(k -1)
        tank_pressure = init_vapour_pressure * np.power(vapour_temp/init_vapour_temp, tmp)
        current_z = compress_factor(tank_pressure, p_crit, z_crit)

        old_aim = aim

        if current_z_guess < current_z: #Z guessed is to little
            current_z_guess *= step 
            aim = 1
        else:
            current_z_guess /= step #Z guessed is to large
            aim = -1

        if aim == -old_aim: #if the target is overshoot reduce step so not to create inifite loop
            step = np.sqrt(step)
    
        if current_z_guess/current_z > (1 + 10**(-6)) or  current_z_guess/current_z < 1/(1 + 10**(-6)):
            pass
        else:
            go = False
    # Nor sure if this part should be in while but is not used so i leave it after it

    tmp = 1/(k - 1)
    vapour_denisty = init_vapour_density * np.power(vapour_temp/init_vapour_temp, tmp)





    


    

    


