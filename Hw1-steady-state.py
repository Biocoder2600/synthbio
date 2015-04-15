# -*- coding: utf-8 -*-
"""
Created on Thu Apr 09 08:30:18 2015

@author: Raandrew
"""
import tellurium as te
import math

# simulates the model given by antModelStr and the given constants
# returns the steady state value of the protein
def simulate(antModelStr, constants):
    rr = te.loada(antModelStr.format(**constants))
    rr.simulate(0, 800000, 1000)
    rr.plot()
    rr.getSteadyStateValues()
    return rr.P
    
# design parameters for pPro-ref device
kinittxn = 0.12
kpol = 25.0
ORF = 678.0
kinittrans = 2.0
ktrans = 21.0
kfold = 0.0
kobsminus = 0.0 / 60
kobsplus = 0.0 / 60
RNAthalf = 1.2 * 60
Proteinthalf = 72.0 * 3600
d = 275.0

# constants for the pPro-ref model
k1 = kinittxn
k2 = math.log(2) / (0.5 * (d / kpol))
k3 = math.log(2) / (0.5 * (ORF / kpol))
k4 = 2.0
k5 = math.log(2) / (0.5 * (ORF / ktrans))
k6 = kfold
k7 = kobsminus
k8 = kobsplus
k9 = math.log(2) / RNAthalf
k10 = 0.1745  * k9
k11 = math.log(2) / Proteinthalf

# dictionary of constants for pPro-ref model
constants_ref = {'k1': str(k1), 'k2': str(k2), 'k3': str(k3), 'k4': str(k4), 
'k5': str(k5), 'k6': str(k6), 'k7': str(k7), 'k8': str(k8), 'k9': str(k9),
'k10': str(k10), 'k11': str(k11)}

# design parameters for rRed13 static device
kinittxn = 0.12
kpol = 25.0
ORF = 678.0
kinittrans = 2.0
ktrans = 21.0
kfold = 0.052
kobsminus = 3.64 / 60
kobsplus = 0.0 / 60
RNAthalf = 1.2 * 60
Proteinthalf = 72.0 * 3600
d = 275.0

# constants for rRed13 model
k1 = kinittxn
k2 = math.log(2) / (0.5 * (d / kpol))
k3 = math.log(2) / (0.5 * (ORF / kpol))
k4 = 2.0
k5 = math.log(2) / (0.5 * (ORF / ktrans)) # 21.0
k6 = kfold
k7 = kobsminus
k8 = kobsplus
k9 = math.log(2) / RNAthalf
k10 = 0.1745 * k9
k11 = math.log(2) / Proteinthalf

# dictionary of constants for rRed13 model
constants_rRed13 = {'k1': str(k1), 'k2': str(k2), 'k3': str(k3), 'k4': str(k4), 
'k5': str(k5), 'k6': str(k6), 'k7': str(k7), 'k8': str(k8), 'k9': str(k9),
'k10': str(k10), 'k11': str(k11)}

# string that defines the model
antModel = """
    
    -> I; k1
    I -> Uppp; k2 * I
    Uppp -> Rppp; k3 * Uppp
    UOH -> ROH; k3 * UOH
    Upppf ->Rpppf; k3 * Upppf
    Uppp -> Uppp + T; k4 * Uppp
    Upppf -> Upppf + T; k4 * Upppf
    UOH -> UOH + T; k4 * UOH
    Rppp -> Rppp + T; k4 * Rppp
    Rpppf -> Rpppf + T; k4 * Rpppf
    ROH -> ROH + T; k4 * ROH
    T -> P; k5 * T
    
    Uppp -> Upppf; k6 * Uppp
    Upppf -> UOH; k7 * Upppf
    Upppf -> UOH; k8 * Upppf
    Rppp -> Rpppf; k6 * Rppp
    Rpppf -> ROH; k7 * Rpppf
    Rpppf -> ROH; k8 * Rpppf
    
    Uppp -> ; k9 * Uppp
    Upppf -> ; k9 * Upppf
    Rppp -> ; k9 * Rppp
    Rpppf -> ; k9 * Rpppf
    UOH -> ; k10 * UOH
    ROH -> ; k10 * ROH
    P -> ; k11 * P

    I = 0;
    Uppp = 0;
    Upppf = 0;
    Rppp = 0;
    Rpppf = 0;
    UOH = 0;
    ROH = 0;
    P = 0;
    
    k1 = {k1}; k2 = {k2}; k3 = {k3}; k4 = {k4}; k5 = {k5};
    k6 = {k6}; k7 = {k7}; k8 = {k8}; k9 = {k9}; k10 = {k10};
    k11 = {k11}"""

# run simulations for the pPro-ref device and rRed13 device
P_ref = simulate(antModel, constants_ref)
print "The above plot does not go to steady state because Tellurium cannot plot any more points."
print "At steady state, P_ref =", P_ref
P_rRed13 = simulate(antModel, constants_rRed13)
print "The above plot does not go to steady state because Tellurium cannot plot any more points."
print "At steady state, P_rRed13 =", P_rRed13
print "gamma_rel =", (P_rRed13 / P_ref)


