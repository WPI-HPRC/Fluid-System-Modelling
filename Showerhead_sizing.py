import numpy as np

#import matplotlib as plt

def valve_incompressible(mdot,density,cv): 
    dP = (density/density_water)*(mdot/(cv*density))**2
    return dP

def orifice_incompressible(mdot,density,cd,A):
    dP = density*(mdot/(density*cd*A))**2
    return dP

def orifice_compressible_chocked(mdot,density,gamma,cd,A):
    P_up = (1/gamma*density)*((mdot/(cd*A))**2)*(2/(gamma+1))**((gamma+1)/(gamma-1))
    return P_up

def orifice_compressible_unchoked():
    return

def dia2area(diameter):
    area = (1/4)*np.pi*diameter**2
    return area

def area2dia(area):
    dia = (4*area/np.pi)**0.5
    return dia


# desire properties of injector
p_injector = 300 *6894.76 # 300 psi in Pa
mdot_injector = 3 * 0.453592 # 3 lbm/s in kg/s
burntime = 5 #seconds



# mechanical properties of water
density_water = 1000 # kg/m^3
viscosity_water = 1e-6 #dynamic viscosity of water in [Pa*s] @20C

# mechanical properties of nitrogen
density_nitrogen = 1.19 # kg/m^3
gamma_nitrogen = 1.4


vdot_injector = mdot_injector/density_water #volumetric flow rate


# tank properties
water_fill_volume = mdot_injector*burntime/density_water
tank_diameter = 6 * 0.0254 # 6in to meters
fill_height = water_fill_volume/dia2area(tank_diameter)
hydrostatic=density_water*10*fill_height


component = np.array(["orifice","ball_1","check_1",'tube','water'])
component_type = np.array([0,1,1,2,3]) # 0 = orifice, 1 = valve, 2 = tube, 3 = water
cv = np.array([0,2,4.7])
orifice_diameter = 0.0254 * np.array([0.5,0,0]) # diameter in inches converted to meters
cd = 0.6
orifice_area = dia2area(orifice_diameter)

tube_length = 0.0254 * np.array([0,0,0,10,0]) # tube length converted from in to m
tube_diameter = 0.0254 * np.array([0,0,0,1,0]) # tube diameter converted from in to m
tube_velocity = mdot_injector/(dia2area(tube_diameter)*density_water)
flow_regime = np.zeros(tube_length.size)
roughness = 0.0000025 # in meters

reynolds = np.zeros(tube_diameter.size)
for x in range(tube_diameter.size):
    reynolds[x] = density_water*tube_velocity[x]*tube_length[x]/viscosity_water
    if reynolds[x] <= 2000: # laminar flow
        flow_regime[x] = 0
    if reynolds[x] >= 2000 and flow_regime[x]<= 4000: # transition
        flow_regime[x] = 1
    if reynolds[x] >= 4000: # turbulent flow
        flow_regime[x] = 2



chocked = np.zeros(component.size)
dP = np.zeros(component.size)
P = np.zeros(component.size+1)
P[0] = p_injector


test = orifice_incompressible(mdot_injector,density_water,0.7,1e-4)
diameter = 39.3701*(4*1e-4/np.pi)**0.5



for x in range(component_type.size):
    print(x)
    print(component[x])
    if component_type[x] == 0: # 0 = orifice
        dP[x] = orifice_incompressible(mdot_injector,density_water,cd,orifice_area[x])
    if component_type[x] == 1: # 1 = valve
        dP[x] = valve_incompressible(mdot_injector,density_water,cv[x])
    if component_type[x] == 2: # 2 = tube
        if flow_regime[x] == 0: # laminar flow
            friction_factor = 64/reynolds[x]
        if flow_regime[x] == 1: # transition
            friction_factor = 99999999999 # dont feel like coding yet. 
        if flow_regime[x] == 2: # turbulent
            friction_factor = (1.14 + 2* np.log10(tube_diameter[x]/roughness))**-2
        headloss = friction_factor * tube_length*(tube_length[x]**2)/(2*tube_diameter[x]*9.81)
        dP[x] = headloss[x] * 9.81 * density_water
    if component_type[x] == 3: # 3 = water
        dP[x] = -1 * hydrostatic # negative dP bc head pressure increase
    print(0.000145038*dP[x])
    P[x+1] = P[x] + dP[x]
       
     
     
     

   
P_imperial = P*0.000145038
print(dP)
print(P_imperial)

hi = 2 
diameter_nitrogen = 0.1*.00254
area_nitrogen = dia2area(diameter_nitrogen)
orificechockpressureup = orifice_compressible_chocked(mdot_injector,density_nitrogen,gamma_nitrogen,0.7,area_nitrogen)
criticalpressure = orificechockpressureup*(2/(gamma_nitrogen+1))**(gamma_nitrogen/(gamma_nitrogen-1))
criticalpressureratio = criticalpressure/orificechockpressureup
print(orificechockpressureup*0.000145038)
# print(criticalpressure)
# print(criticalpressureratio)

