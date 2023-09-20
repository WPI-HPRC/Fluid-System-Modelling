import numpy as np
#import matplotlib as plt

def valve_incompressible(mdot,density,cv): 
    dP = (mdot/(cv*density))**2
    return dP

def orifice_incompressible(mdot,density,cd,A):
    dP = density*(mdot/(density*cd*A))**2
    return dP

p_injector = 300 *6894.76 # 300 psi in Pa
mdot_injector = 3 * 0.453592 # 3 lbm/s in kg/s
burntime = 5 #seconds
density_water = 1000 # kg/m^3
water_fill_volume = mdot_injector*burntime/density_water
tank_diameter = 6 * 0.0254 # 6in to meters


print(water_fill_volume)

component = np.array(["ball_1","check_1"])
cv = np.array([2,4.7])
chocked = np.zeros(component.size)
dP = np.zeros(component.size)
P = np.zeros(component.size*2-1)
P[0] = p_injector
dP[0] = valve_incompressible(mdot_injector,density_water,cv[0])
dP[1] = valve_incompressible(mdot_injector,density_water,cv[1])

test = orifice_incompressible(mdot_injector,density_water,0.7,1e-4)
diameter = 39.3701*(4*1e-4/3.14)**0.5
print(diameter)