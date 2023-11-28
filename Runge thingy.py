#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#defining variables:
#Ndot equation:
Sg = 0.9737 #Specific gravity of N2
deltaP = 1 #Ask for pressure values
P1 = 1 #Ask for before oriface pressure
P3 = 1 #Ask for final pressure for the oriface
Cv = 1.2 #go over equation and get values needed (guess for now)


#E dot equation
dotW = 0 #no pump
h = 45.989 #enthalpy calculated from href
dotQ = 0 #0 as the temperature is the same at room temp
hRef = 191.609 #from janaf
T = 293.15 #Kelvin or 20* Celcius
e = 1 #ask for the values of mass + energy we are using
mDot = 1 #ask for our mass, mass flow should be constant as we have one inlet/outlet for the oriface
c = 10 #ask for velocity we are using/calculated
Cp = 29.124 #specific heat of constant pressure for N2
Tref = 298.15 #Kelvin
v = 0.872 #m^3/kg specific volume


#Other N dot equations
mBar = 0.02802 #kg/mol molar mass of N2
Cd = 0.60 # coefficent of discharge (we are using sharp edge)
rho = 808.4 #kg/m^3
Rbar = 8.3144598 #J/mol*K, universal gas constant
gamma = 1.40 #ratio of specific heats
Na = 6.023*(10**23) #Avogadro's constant
beta = 1 #What are our dimensions for the oriface + pipe?
P2 = 1 #ask for more info
A = 1 #ask if cross sectional or other




#Ndot equation Unchoked:
def functionNdotUnchoked(Na,mBar,Cd,beta,A,rho,P1,gamma,P3):
    val1 = ((P3/P1)**(2/gamma))-((P3/P1)**((gamma+1)/(gamma)))
    val2 = 2*rho*(gamma/(gamma-1))
    val3 = (Na/mBar)*(Cd/(1-beta*beta)**(1/2))*A
    return val3*((val2*val1)**(1/2)) 

#Energy Transport equation:
def functionEdot(dotQ,dotW,mDot,h,c):
    val4 = (dotQ+dotW)+((mDot*(h+((c**2)/2)))-(mDot*(h+((c**2)/2))))
    return val4                        

#Valve Mass Transport equation:
def functionNdot(Na,mBar,rho,Cv,P1,P3,deltaP,Sg,T):
    val5 = (964)*((Na*rho)/mBar)*Cv
    val6 = (((P1+P3)*deltaP)/(Sg*T))**(1/2)
    return val5*val6


# Python program to implement Runge Kutta method
# A sample differential equation "dy / dx = (x - y)/2"
def dydx(x, y):
	return ((x - y)/2)

# Finds value of y for a given x using step size h
# and initial value y0 at x0.

def rungeKutta(x0, y0, x, h):
	# Count number of iterations using step size or
	# step height h
    n = (int)((x - x0)/h) 
	# Iterate for number of iterations
    y = y0
    arr = [functionEdot(dotQ,dotW,mDot,h,c), functionNdot(Na,mBar,rho,Cv,P1,P3,deltaP,Sg,T), functionEdot(dotQ,dotW,mDot,h,c), functionNdot(Na,mBar,rho,Cv,P1,P3,deltaP,Sg,T), functionEdot(dotQ,dotW,mDot,h,c), functionNdot(Na,mBar,rho,Cv,P1,P3,deltaP,Sg,T)]
    for x in arr:
    	for i in range(1, n + 1):
    		"Apply Runge Kutta Formulas to find next value of y"
    		k1 = h * dydx(x0, y)
    		k2 = h * dydx(x0 + 0.5 * h, y + 0.5 * k1)
    		k3 = h * dydx(x0 + 0.5 * h, y + 0.5 * k2)
    		k4 = h * dydx(x0 + h, y + k3)
    
    		# Update next value of y
    		y = y + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
    
    		# Update next value of x
    		x0 = x0 + h
    return y

# Driver method
x0 = 0
y = 1
x = 2
h = 0.2
print ('The value of y at x is:', rungeKutta(x0, y, x, h))

# This code is contributed by Prateek Bhindwar
