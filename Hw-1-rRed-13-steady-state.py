# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 21:46:42 2015

@author: mmurthy2
"""

import tellurium as te
import math

def simulate(antModelStr, constants):
    rr = te.loada(antModelStr.format(**constants))
    rr.simulate(0, 800, 100)
    rr.plot()
    print "P =", rr.P


kinittxn = 0.12
kpol = 25
ORF = 678
kinittrans = 2
ktrans = 21
kfold = 0.052
kobsminus = 3.64 * 60
kobsplus = 0 * 60
RNAthalf = 1.2 * 60
Proteinthalf = 72 * (60 * 60)
d = 380

k1 = kinittxn
k2 = math.log(2) / (0.5 * (d / kpol))
k3 = math.log(2) / (0.5 * (d / kpol))
k4 = (math.log(2) / (0.5 * 1/kinittxn))
k5 = math.log(2) / (0.5 * (d / ktrans))
k6 = kfold
k7 = kobsminus
k8 = kobsplus
k9 = math.log(2) / RNAthalf
k10 = 0.1745 * k9
k11 = math.log(2) / Proteinthalf

constants = {'k1': str(k1), 'k2': str(k2), 'k3': str(k3), 'k4': str(k4), 
'k5': str(k5), 'k6': str(k6), 'k7': str(k7), 'k8': str(k8), 'k9': str(k9),
'k10': str(k10), 'k11': str(k11)}

antModel = """
    -> I; k1
    I -> Uppp; k2 * I
    Uppp -> Rppp; k3 * Uppp
    UOH -> ROH; k3 * UOH
    Upppf -> Rpppf; k3 * Upppf
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

simulate(antModel, constants)
