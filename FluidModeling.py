# FluidModeling.py
# Ryan Hunter
# HPRC-Propulsion Team (Experimental)

#need python to run this code and also need to install numpy (extrenal library fro arrays and stuff)


import numpy as arr
import math

#Assumptions && Extra Info that I am too lazy to get rid of/check claims of

# Goal -- Final PSI should be 550
# Cd = .6
# beta ~ 0.0   -->  C = Cd                      (beta = actual diameter / theoretical dimeter)
# Q* = Mdot / rho                               (QSTAR = mass flow / density)
# We finding delta P                            (same as p2 - p1)
# Qm = rho * Qv
# Can ignore G bc we using water (1)            (included as variable since this may change much later)
# A is cross-sectional Area
# Mdot = 3 lbm/s for now                        (target mass flow rate for rocket)

## Formulas

# deltaP = [ (QSTAR * sqrt(G) ) /Cv ]^2        # sqrt(G) is 1 since we are modeling water (if it were something else then it would matter)
# (FINAL_PSI - initial_PSI) = [ (QSTAR * sqrt(G) ) /Cv ]^2
# initial_PSI = FINAL_PSI - [ (QSTAR * sqrt(G) ) /Cv ]^2

# C = CD / math.sqrt( 1 - ( Beta**4 ) )

# Qv = C * A * math.sqrt( ( 2 (Delta_P--------------) ) / DENSITY)
# Qm = Qv * DENSITY


## really need to convert to SI units for equations (ratios are above)












## ------------- Most Important Information You Will Ever See -------------
# Arrays for values to make equations work (all the arrays must be the same length) (this is what changes pressures)

# These are the only variables you will have to change (unless you use different stuff going through the pipes)

valveCValues = arr.array([0.6, 1.15, 1.38, 0.6])                 # how pressure is changing

typeHole = arr.array(["Orifice", "Valve", "Valve", "Valve"])         # what is changing the prssure

holeArea = arr.array([7.0, 9.0, 10.0, 0.005])                      # size of hole

resultsOfPSI = arr.zeros(4)                                 # what PSI was before typeHole

CompressOrIncompress = arr.array(["Incompressible", "Incompressible", "Incompressible", "Compressible"])












##  ------------- Constants (As defined by water and nitrogen) -------------

# Ratios
PSI_TO_PA = 6894.76
LBS_TO_KGS = 0.453592
IN_TO_M = 0.0254


## Units we want stuff in terms of
# meter
# Kilograms
# kelvin
# seconds


# Important Numbers that change depending on the system/requirements        (Thanks Terence for the simple code unit conversion solution)
DENSITY = 997.0 #kg/m^3  # currently se to the density of water (greek letter of rho)
N2_DENSITY = 1.4
BURNTIME = 5.0 #s
TANK_DIAMETER = 6.0 * IN_TO_M #in (need to multiply by constant to make it work in equations)
FINAL_PSI = 550.0 * PSI_TO_PA #psi
MDOT = 2.5 * LBS_TO_KGS #Mdot stands for mass flow #lbs/sec (multiplied by LBS_TO_KGS so that the equations work) 
MDOT_N2 = 0.1 * LBS_TO_KGS
QSTAR = (MDOT / DENSITY) #gallons per second   # flow rate is equal to mass flow/density
G_VAR = 1 # Since we are modeling water, it is 1 (Specific gravity of fluid)
CD = 0.6        #orifice coefficient to show how good it is at flwoing (this value is only for the very final orifice since the other orifices will be different)
GAMMA = 1.4
How_TANK_FILLS_WATER = (MDOT * BURNTIME) / DENSITY
NEWTON_RAPHSON_DELTA_X = 0.001
Beta = 0                                                    # important if system did not have a beta which is 0 (we going too slow)












##  ------------- Incompressible Flow Modeling -------------

def calculatePsi_incompressible(final_psi, holeName, count):
    cVar = valveCValues[count]
    hArea = holeArea[count]
    if str(holeName) == "Orifice":
        #print(final_psi)
        resultsOfPSI[count] = ( (final_psi + ( ( QSTAR * (G_VAR**0.5) ) / cVar )**2.0) / PSI_TO_PA )
        count += 1
        #testSTR = str( resultsOfPSI[count - 1] )
        #print("" + testSTR + " top")
        return count
    elif str(holeName) == "Valve":
        #print(final_psi)
        resultsOfPSI[count] = (final_psi + ( ( MDOT**2 ) / ( 2.0 * ( cVar**2 ) * ( hArea**2 ) * DENSITY) ))
        count += 1
        #testSTR = str( resultsOfPSI[count - 1] )
        #print("" + testSTR + " bot")
        return count
    else:
        count = -2
        return count












##  ------------- Compressible Flow Modeling -------------

def calculatePsi_compressible(final_psi, holeName, count):     #would love to have a do-while loop here but it seems python does not have this
    
    if str(holeName) == "Orifice":
        #only place where fluids can choke (check if choked flow here somewhere)
        print("you should not be here as this is not done")
    elif str(holeName) == "Valve":
        #We are using the equation "Compressible Flow (Unchoked)" for this part of the simulation
        #upstream flow is p2 for eq so I am finding p2 (since I am working backwards)
        #Using the Newton Raphson Method to obtain a very close-to-actual value for p2
        resultsOfPSI[count] = funcForNewtRaphMethod(final_psi, count)
    else:
        count = -2
        return count
    return (count + 1)



def calcUsingFlowEquation(p1, count, p2):           #this method does the math for what value that would be computed by plugging pressure and other values into equation with new x
    cVar = valveCValues[count]
    hArea = holeArea[count]                         #Since I am working backwards, I want to solve for p1 instead of p2
    pressureApproximation = (cVar * hArea) * ( 2 * N2_DENSITY * p1 * ( GAMMA / (GAMMA - 1) ) * ( (p2 / p1)**(2/GAMMA) - (p2 / p1)**( (GAMMA + 1) / GAMMA) ) )**0.5 - MDOT_N2
    return pressureApproximation




def funcForNewtRaphMethod(p2, count):    #this method actually uses the Newton-Raphson Method to compute the value of pressure
    notAccurateEnough = True
    pressureOne_Old = p2
    loopCount = 0
    while notAccurateEnough:
        pressureOne_New = pressureOne_Old - ( (NEWTON_RAPHSON_DELTA_X * calcUsingFlowEquation(pressureOne_Old, count, p2)) / (calcUsingFlowEquation((pressureOne_Old + NEWTON_RAPHSON_DELTA_X), count, p2) - calcUsingFlowEquation(pressureOne_Old, count, p2) ) )
        if ( int(pressureOne_New * 1000000) == int(pressureOne_Old * 1000000) or loopCount == 20):          #if the Newton-Raphson Method made a negligible change to the value of X then you can stop (accurate enough to actual value)
            notAccurateEnough = False
            return pressureOne_New
        pressureOne_Old = pressureOne_New
        loopCount = loopCount + 1












##  ------------- Math For Everything (Delta PSI Control Room) -------------

# This is where all the output is displayed/Correct Algorithm is Used

counter = -1
#resultsOfPSI[0] = 300
endPsi = 300

for x in typeHole:
    if counter == -2:
        raise Exception("You need to either insert an orifice or valve for simulation")
    elif(counter == -1):
        counter = 0
        if (CompressOrIncompress[counter] == "Incompressible"):
            counter = 0
            counter = calculatePsi_incompressible(FINAL_PSI, x, counter)
        else:
            raise Exception("You are dumb or special")
    else:
        #print("Got here")
        if (CompressOrIncompress[counter] == "Incompressible"):
            counter = calculatePsi_incompressible(resultsOfPSI[counter-1], x, counter)
        else:
            counter = calculatePsi_compressible(resultsOfPSI[counter-1], x, counter)


counter2 = 0
print("\n")
print("Position\tPressure Value\t\tMedium Type\n")
print("Start:\t\t550")
for lcv in typeHole:
    print(counter2)
    print(lcv + " : \t" + str(resultsOfPSI[counter2]) + "\t" + CompressOrIncompress[counter2])
    counter2 += 1

print("\n")
