# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 09:09:32 2015

@author: mmurthy2
"""

import tellurium as te
import math
import random
import csv
import numpy as np

def simulate(antModel, constants):
    rr = te.loada(antModel.format(**constants))
    rr.getSteadyStateValues()
    return rr.P

def reactionModel(coefficients, result):
    kinittxn = coefficients[0]
    kpol = coefficients[1]
    ktrans = coefficients[2]
    kfold = coefficients[3]
    kobsminus = coefficients[4]
    kobsplus = coefficients[5]
    RNAthalf = coefficients[6]
    Proteinthalf = coefficients[7]
    EC = coefficients[8]
    
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
        
    # same as reference model
    ORF = 678
    d = 275
    L = (0.0000000000000001, 1)
   
    if (kobsplus > kobsminus):
        # constants for the pPro-ref model
        k1 = kinittxn
        k2 = math.log(2) / (0.5 * (d / kpol))
        k3 = math.log(2) / (0.5 * (ORF / kpol))
        k4 = 2.0
        k5 = math.log(2) / (0.5 * (ORF / ktrans))
        k6 = kfold
        k7 = 0
        k8 = (kobsplus - kobsminus) * (1 / (1 + (EC/L[0]))) + kobsminus
        k9 = math.log(2) / RNAthalf
        k10 = 0.1745  * k9
        k11 = math.log(2) / Proteinthalf
        
        # dictionary of constants for ref model
        constants_ref = {'k1': str(k1), 'k2': str(k2), 'k3': str(k3), 'k4': str(k4), 
        'k5': str(k5), 'k6': str(k6), 'k7': str(k7), 'k8': str(k8), 'k9': str(k9),
        'k10': str(k10), 'k11': str(k11)}

        # compute outputs for sampled parameter sets when FAEE production is zero.
        P_ref = simulate(antModel, constants_ref)
        # print "At steady state, P_ref = ", P_ref
        
        k8 = (kobsplus - kobsminus) * (1 / (1 + (EC/L[1]))) + kobsminus
        
        constants_aRED = {'k1': str(k1), 'k2': str(k2), 'k3': str(k3), 'k4': str(k4), 
        'k5': str(k5), 'k6': str(k6), 'k7': str(k7), 'k8': str(k8), 'k9': str(k9),
        'k10': str(k10), 'k11': str(k11)}
        
        # compute outputs for the same sampled parameter sets when FAEE production is maximal.
        P_aRED = simulate(antModel, constants_aRED)
        # print "At steady state, P_aRED = ", P_aRED
        
        gamma_rel = (P_aRED / P_ref)
        
        result.append((kinittxn, kpol, ktrans, kfold, kobsminus, kobsplus, RNAthalf, Proteinthalf, EC, gamma_rel))
    
    
# use the optimal coefficients determined from R
baseoptcoef = [0.161934555, 213.54424, 51.17770, 1.954204876, 8.480704e-08, 
               1.572963e-05, 57976.99421,  127708.1747, 0.3464432996]

print 'Base opt coeffs', baseoptcoef

for i in range(0, 9):
    uncert = np.zeros(9)
    # get uncertainty range
    for j in range(0, len(baseoptcoef)):
        if (i == j):
            uncert[j] = baseoptcoef[j] * 0.50
        else:
            uncert[j] = baseoptcoef[j] * 0.10
    
    print 'Uncertainty', uncert
    
    result = []
    for k in range(0, 50):
        params = np.zeros(9)
        for l in range(0, len(baseoptcoef)):
            params[l] = random.uniform(baseoptcoef[l] - uncert[l], baseoptcoef[l] + uncert[l])
        
        reactionModel(params, result)

    np_result = np.array(result)
    
    print 'Result size', np.shape(np_result)
    
    # find the top and bottom 2.5% gamma rel, remove them
    top = np.percentile(np_result[:, 9], 97.5)
    bottom = np.percentile(np_result[:, 9], 2.5)
    
    np_result = np_result[np_result[:, 9] < top]
    np_result = np_result[np_result[:, 9] > bottom]
    
    print 'top', top
    print 'bottom', bottom
    print 'Result', np_result
    
    filename = 'results/result' + str(i) + '.csv'
    
    print 'filename', filename
    
    result = list(np_result)
    
    with open(filename, 'w') as out:
        csv_out = csv.writer(out)
        csv_out.writerow(('kinittxn', 'kpol', 'ktrans', 'kfold', 'kobsminus', 'kobsplus', 'RNAthalf', 'Proteinthalf', 'EC', 'gamma_rel'))
        for row in result:
            print row
            csv_out.writerow(tuple(row))
