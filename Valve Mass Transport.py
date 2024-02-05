import numpy as np
import matplotlib.pyplot as plt

def f(X, dt, const):

    Ndotin = X[0]
    Ndotout = X[2]
    Edotin = X[1]
    Edotout = X[3]

    Sgin = const[0]
    Sgout = const[15]
    Cv = const[1]
    Href = const[2]
    Tref = const[3]
    Mbarin = const[4]
    Mbarout = const[19]
    Cd = const[5]
    Rbar = const[6]
    Gamma = const[7]
    Nain = const [8]
    Naout = const[23]
    Beta = const[9]
    Ain = const[10]
    Aout = const[25]
    P1in = const[11]
    P1out = const[26]
    P3in = const[12]
    P3out = const[27]
    Vin = const[13]
    Vout = const[28]
    min = const[14]
    mout = const[29]
    Tin = const[30]
    Tout = const[31]

    rhoin = (Ndotin*Mbarin)/(Nain*Vin)
    rhoout = (Ndotout*Mbarout)/(Naout*Vout)
    mdotin = 964(rhoin(np.sqrt((P1in ^ 2) - (P3in ^ 2)) / (Sgin * Tin)))
    mdotout = 964(rhoout(np.sqrt((P1out ^ 2) - (P3out ^ 2)) / (Sgout * Tout)))
    ein = Edotin /min
    eout = Edotout / mout
    hin = ein + (P1in * Vin)
    hin = eout + (P1out * Vout)
    cin = mdotin / (rhoin * Ain)
    cout = mdotout / (rhoout * Aout)

    DeltaPin = P1in - P3in
    DeltaPout= P1out - P3out
    Wdot = #work on fluid

    Qdot = #heat flux
    Cpin = ((Gamma*Rbar)/(Gamma-1))
    Cpout =
    v = ((Vin + Vout)/2/m)
    NDot = 964*(((Na*p1)/Mbar)*Cv)*np.sqrt(((P1 - P3)*deltaP)/(Sg*T))
    Edot = Qdot +Wdot +((mdot(h+((c^2)/2))-()
    P = (Ndot*Rbar*Tref)/(Na*V*Mbar)
    m = (Ndot*Mbar)/(Na)



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
        t = t + (( 1 /6 ) *(t1 + 2 *t2 + 2 *t3 +t4))

        XRec = np.append(XRec ,np.reshape(XCurr, (XCurr.shape[0], 1)) ,axis=1)
        tRec = np.append(tRec, t)

    return [XRec, tRec]
X = np.array([1])
const = [1,1,1,1,1,1,1,1,1,1,1,1]
[XRec, tRec] = RK4(f,X,const,0)



fig,ax1 = plt.subplots()
ax1.plot(tRec,XRec[0])
plt.show()

