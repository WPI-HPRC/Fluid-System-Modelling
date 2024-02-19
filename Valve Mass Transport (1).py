import numpy as np
import matplotlib.pyplot as plt

def f(X, t, const):

    N = X[0]
    E = X[1]
    

    Sg = const[0] #specific Gravity
    Cv = const[1] #Flow Coefficient
    Href = const[2] #reference Enthalpy
    Tref = const[3] #reference temperature
    Mbar = const[4] #molar mass of fluid
    Cd = const[5] #Discharge Coefficient
    Rbar = const[6] #Universal Gas Constant
    Gamma = const[7] #ratio of specific heats
    Na = const[8] #Avogadro's number
    Beta = const[9] #orifice diameter/pipe diameter
    A1 = const[10] #area in
    A2 = const[11] #area out
    P1 = const[12] #Pressure at point 1
    P3 = const[13] #Pressure at point 3
    V = const[14] #volume 
    min = const[15] #mass of fluid in 
    mout = const[16] #mass of fluid out 
    dt = const[17]

 #specifc energy
    ein = E / min
    eout = E / mout 
#specifc heat at constant pressure 
    Cp= ((Gamma*Rbar)/(Gamma-1))
#mass
    m = (N * Mbar)/Na
#Specific Volume 
    v = V/m
 #pressure 
    P1 = (1/(1-((N*Rbar*v)/(Na*V*Mbar*Cp))))*((((N*Rbar)/(Na*V*Mbar*Cp))*(ein - Href) )-((N*Rbar*Tref)/(Na*V*Mbar*Cp)))
    P3 = (1/(1-((N*Rbar*v)/(Na*V*Mbar*Cp))))*((((N*Rbar)/(Na*V*Mbar*Cp))*(eout - Href) )-((N*Rbar*Tref)/(Na*V*Mbar*Cp)))
#Enthalpy
    hin = ein + (P1*v)
    hout = eout + (P3*v)
#average specific energy
    eavg = (ein + eout)/2
# Average Density   
    rhoin = 1/v
    rhoout = 1/v
#Temperature 
    Tin = (hin - Href)/Cp
    Tout = (hout - Href)/Cp
    Tavg = (Tin + Tout)/2
# Mass flux
    mdotin = 964*(rhoin*(np.sqrt((P1 ** 2) - (P3 ** 2)) / (Sg * Tin)))
    mdotout = 964*(rhoout*(np.sqrt((P1 ** 2) - (P3 ** 2)) / (Sg * Tout)))
#fluid velocity
    cin = mdotin / (rhoin * A1)
    cout = mdotout / (rhoout * A2)
#Change in pressure
    DeltaP = P3 - P1
    Wdot = 0 #work on fluid
#test variable
    test = (N*Rbar*v) #(1/(1-((N*Rbar*v)/(Na*V*Mbar*Cp))))
    Qdot = 0 #heat flux
#Change in number of molecules (Valve mass Transport)
    NDot = 964*(((Na*rhoin)/Mbar)*Cv)*((((P1 + P3)*DeltaP)/(Sg*Tavg))**0.5)
#Energy Transport
    Edot = Qdot + Wdot + ((mdotin*(hin+((cin**2)/2))-(mdotout*(hout+((cout**2)/2)))))


#mass
    
#Orifice mass transport
#Unchoked mass transport
    NdotU = (Na/Mbar)*(Cd/np.sqrt(1-Beta**2))*A2*((2*rhoin*P1*(Gamma/(Gamma-1))*((P3/P1)**(2/Gamma)-(P3/P1)**((Gamma+1)/Gamma)))**0.5)  
#choked Mass transport
    NdotC = (Na/Mbar)*(Cd/np.sqrt(1-Beta**2))*A2*((Gamma*rhoin*P1*(2/(1+Gamma))**((Gamma+1)/(Gamma-1)))**0.5)
#Choked Condition
    choked_cond = ((A2*rhoin)/m) * (((Gamma*Rbar*Tavg)/(Mbar*(1+((Gamma-1)/2))))*((2/(Gamma+1))*(1+(Gamma-1)/2)*((mdotin/(A2*rhoin))**2)*((Mbar*(1+(Gamma -1)/2))/(Gamma*Rbar*Tavg)))**0.5)
#check on whether orifice mass transport is 
    if choked_cond < (A1/A2):
        NdotO = NdotC
    else:
        NdotO = NdotU

# temp function

    return [NDot, dt]

def RK4(func ,X ,const ,t0):

    XCurr = np.copy(X)
    XRec = np.zeros((X.size,1))
    XRec[: ,0] = np.copy(XCurr)

    tRec = np.zeros(1)
    tRec[0] = t0
    t = t0
    while(t < 10):
        [k1 ,t1] = func(XCurr ,t ,const);
        [k2 ,t2] = func(XCurr +t1 *(k1/2) , t +(t1/2) ,const);
        [k3 ,t3] = func(XCurr +t2 *(k2/2) , t +(t2/2) ,const);
        [k4 ,t4] = func(XCurr +t3 *k3 , t +t3 ,const);

        XCurr = XCurr + (( 1 /6 ) *(t1 *k1 + 2 *t2 *k2 + 2 *t3 *k3 +t4 *k4))
        t = t + (( 1 /6 ) * (t1 + 2 *t2 + 2 *t3 +t4))

        XRec = np.append(XRec ,np.reshape(XCurr, (XCurr.shape[0], 1)) ,axis=1)
        tRec = np.append(tRec, t)

    return [XRec, tRec]
X = np.array([3,4])
const = [1,50,1511.8,20,18,6.9,.287,1.05, 6.22*10**23,.5,4,5,101.3,600,15,20,17,0.1]
[XRec, tRec] = RK4(f,X,const,0)



fig,ax1 = plt.subplots()
ax1.plot(tRec,XRec[0])
plt.show()