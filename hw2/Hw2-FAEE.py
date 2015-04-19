# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 18:29:33 2015

@author: mmurthy2
"""

import tellurium as te
import math
import random
import csv

def simulate(antModel, constants):
    rr = te.loada(antModel.format(**constants))
    rr.getSteadyStateValues()
    return rr.P

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
    
# reference protein value
P_ref = 9322401.71306
    
ORF = 678
d = 275
RNAthalf = 1.2 * 60
Proteinthalf = 72.0 * 3600

result = []

for i in range(0, 100):
    kinittxn = random.uniform(0.0001, 1)
    kpol = random.uniform(25, 230)
    ktrans = random.uniform(21, 63)
    kfold = random.uniform(0.0001, 3)
    kobsminus = random.uniform(0.001 / 60, 10 / 60)
    kobsplus = random.uniform(0.001 / 60, 10 / 60)
    EC = random.uniform(0.000000001, 1)
    L = random.uniform(0.0000001, 1)
    
    if (kobsplus > kobsminus):
        # constants for the pPro-ref model
        k1 = kinittxn
        k2 = math.log(2) / (0.5 * (d / kpol))
        k3 = math.log(2) / (0.5 * (ORF / kpol))
        k4 = 2.0
        k5 = math.log(2) / (0.5 * (ORF / ktrans))
        k6 = kfold
        k7 = 0
        k8 = (kobsplus - kobsminus) * (1 / (1 + (EC/L))) + kobsminus
        k9 = math.log(2) / RNAthalf
        k10 = 0.1745  * k9
        k11 = math.log(2) / Proteinthalf
        
        # dictionary of constants for pPro-ref model
        constants = {'k1': str(k1), 'k2': str(k2), 'k3': str(k3), 'k4': str(k4), 
        'k5': str(k5), 'k6': str(k6), 'k7': str(k7), 'k8': str(k8), 'k9': str(k9),
        'k10': str(k10), 'k11': str(k11)}

        P_aRED = simulate(antModel, constants)
        gamma_rel = (P_aRED / P_ref)
        print "At steady state, P_aRED = ", P_aRED
        print "gamma_rel = ", gamma_rel
        
        if (gamma_rel > 3.00):
            result.append((kinittxn, kpol, ktrans, kfold, kobsminus, kobsplus, EC, L, gamma_rel))

with open('result.csv', 'w') as out:
    csv_out = csv.writer(out)
    csv_out.writerow(('kinittxn', 'kpol', 'ktrans', 'kfold', 'kobsminus', 'kobsplus', 'EC', 'L', 'gamma_rel'))
    for row in result:
        csv_out.writerow(row)

    
        
