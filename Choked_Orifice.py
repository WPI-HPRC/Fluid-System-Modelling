import numpy as np

import matplotlib as plt

def mach(mdot,area,molarmass,density,gamma,temp):
    M = (mdot/(area*density))*(molarmass*(1+((gamma-1)/2))/(gamma*Rbar*temp))**0.5
    return M

def orificearea(mdot,molarmass,density,gamma,temp):
    area = (mdot/density)*(molarmass*(1+(gamma-1)/2)/(gamma*Rbar*temp))**0.5
    return area
    
def density_idealgas(pressure,molarmass,temp):
    density = pressure*molarmass/(Rbar*temp)
    return density   
    
def density_throat(density_up,gamma):
    density = density_up*(1+((gamma-1)/2))**(-1/(gamma-1))
    return density

def chokedcondition(area_up,area_orifice,mdot,gamma,molar_mass,density_up,temp):
    mach_up = mach(mdot,area_up,molar_mass,density_up,gamma,temp) 
    bruh = (area_up/area_orifice)-(((gamma+1)/2)**(-0.5*(gamma+1)/(gamma-1)))*(mach_up**(-1))*(1+((gamma-1)/2)*mach_up**2)**((gamma+1)/2*(gamma-1))
    return bruh

def epsilon(diameter_tube,diameter_orifice,pressure_1,pressure_2,gamma):
    beta = diameter_orifice/diameter_tube
    e = 1 - (0.351 + 0.256*beta**4 + 0.98*beta**8)*(1-(pressure_2/pressure_1)*1/gamma)
    return e

def compressible_massflow(cd,area_orifice,density_1,diameter_tube,diameter_orifice,pressure_1,pressure_2,gamma):
    e = epsilon(diameter_tube,diameter_orifice,pressure_1,pressure_2,gamma)
    mdot = cd * e * area_orifice * (2*density_1*(pressure_1-pressure_2))**0.5
    return mdot

def func_compressible_ploss(mdot,cd,area_orifice,density_1,diameter_tube,diameter_orifice,pressure_1,pressure_2,gamma):
    e = epsilon(diameter_tube,diameter_orifice,pressure_1,pressure_2,gamma)
    right = cd * e * area_orifice * (2*density_1*(pressure_1-pressure_2))**0.5
    left = mdot
    func = left - right
    return func
    
def dia2area(diameter):
    area = (1/4)*np.pi*diameter**2
    return area

def area2dia(area):
    dia = (4*area/np.pi)**0.5
    return dia

#constants
Rbar = 8.314 #J/mol*K

# flow properties
mdot_water = 3 * 0.453592 # 3 lbm/s in kg/s
density_water = 1000 # kg/m^3
vdot_water = mdot_water/density_water

# mechanical properties of nitrogen
density_nitrogen = 1.19 # kg/m^3
gamma_nitrogen = 1.4
temp_nitrogen = 20 +273 # room temp in celvin
molarmass_nitrogen = 28.02/1000 # g/mol to kg/mol

# press tank properties
pressure_tank = 350 * 6894.76 # Psi to Pa
density_press = density_idealgas(pressure_tank,molarmass_nitrogen,temp_nitrogen)

# ullage properties
pressure_ullage = 300 * 6894.76 # 300 psi to Pa,
density_ullage = density_idealgas(pressure_ullage,molarmass_nitrogen,temp_nitrogen)
vdot_ullage = vdot_water
mdot_ullage = vdot_ullage*density_ullage

# orifice properties
cd = 0.6

pressure_press_guess = 500 * 6894.76 #500 psi to Pa
density_press_guess = density_idealgas(pressure_press_guess,molarmass_nitrogen,temp_nitrogen)
density_orifice = density_throat(density_press_guess,gamma_nitrogen)

min_area_orifice = orificearea(mdot_ullage,molarmass_nitrogen,density_orifice,gamma_nitrogen,temp_nitrogen)
min_diameter_orifice = dia2area(min_area_orifice)

area_orifice = (mdot_ullage/cd)/(gamma_nitrogen*density_press_guess*pressure_press_guess*(2/(gamma_nitrogen+1))**((gamma_nitrogen+1)/(gamma_nitrogen-1)))
diameter_orifice = dia2area(area_orifice)


a2 = dia2area(0.1*0.0254)
m = mach(mdot_ullage,a2,molarmass_nitrogen,density_throat(pressure_press_guess,gamma_nitrogen),gamma_nitrogen,temp_nitrogen)
print(m)



# upstreampressures = np.linspace(150,700, num = 55)
# mindiameter = 39.3701*dia2area(orificearea(mdot_ullage,molarmass_nitrogen,density_orifice,gamma_nitrogen,temp_nitrogen))
# diameter = 39.3701*dia2area((mdot_ullage/cd)/(gamma_nitrogen*density_idealgas(upstreampressures,molarmass_nitrogen,temp_nitrogen)*upstreampressures*(2/(gamma_nitrogen+1))**((gamma_nitrogen+1)/(gamma_nitrogen-1))))
# print(mindiameter)
# plt.plot(upstreampressures,diameter, 'o')
# fig, ax = plt.subplots()

# ax.plot(upstreampressures, diameter, linewidth=2.0)


residual = 1

test_tube_diameter = 1 * 0.0254 # 1in to m
test_orifice_diameter = 0.1 * 0.0254 
test_orifice_area = dia2area(test_orifice_diameter)
test_pressure_1 = 500 * 6894.76 #500 psi to Pa
test_pressure_2 = 100 # guess
test_density_1 = density_idealgas(test_pressure_1,molarmass_nitrogen,temp_nitrogen)

while abs(residual) > 0.000000001:
    func = func_compressible_ploss(mdot_ullage,cd,test_orifice_area,test_density_1,test_tube_diameter,test_orifice_diameter,test_pressure_1,test_pressure_2,gamma_nitrogen)
    change = .0000001
    func_change = func_compressible_ploss(mdot_ullage,cd,test_orifice_area,test_density_1,test_tube_diameter,test_orifice_diameter,test_pressure_1,test_pressure_2 + change,gamma_nitrogen)
    
    test_pressure_2 = test_pressure_2 -(change*func)/(func_change-func)
    residual = func
    print(func)

print('hi')
print(test_pressure_1*0.000145038)
print(test_pressure_2*0.000145038)

# while abs(residual) > 0.000001:
#     func = chokedcondition(tube_area,orifice_area,mdot,gamma_nitrogen,molarmass_nitrogen,density_up,temp_nitrogen)
#     change = 0.0000000000001
#     func_change = chokedcondition(tube_area,orifice_area + change,mdot,gamma_nitrogen,molarmass_nitrogen,density_up,temp_nitrogen)
    
#     orifice_area = orifice_area -(change*func)/(func_change-func)
#     residual = func
#     print(func)

#print(area2dia(orifice_area)*39.3701)

hi = 2
