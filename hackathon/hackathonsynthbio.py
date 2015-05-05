# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 11:49:56 2015

@author: mmurthy2
"""

import tellurium as te
import math
import random
import csv

model = """
    -> T1; k1 / (1 + (inter * k2))
    -> T2; k1 / (1 + (inter * k2))
    -> T3; k2 / (1 + (inter * k1))
    T1 -> P1 + T1; k3 * T1
    T2 -> P2 + T2; k4 * T2
    T3 -> P3 + T3; k5 * T3
    T1 -> ; k6 * T1
    T2 -> ; k6 * T2
    T3 -> ; k6 * T3
    P1 -> ; k7 * P1
    P2 -> ; k7 * P2
    P3 -> ; k7 * P3

    $S1 + P1 -> S2 + P1; k9 * P1 * (S1 / (Km9 + S1))
    S2 + P2 -> S3 + P2; k10 * P2 * (S2 / (Km10 + S2))
    S3 + P3 -> S4 + P3; k11 * P3 * (S3 / (Km11 + S3))
    
    S2 -> ; k12 * S2
    S3 -> ; k12 * S3
    S4 -> ; k12 * S4
    
    S1 = {S1}
    
    k1 = {k1}; k2 = {k2}; k3 = {k3}; k4 = {k4}; k5 = {k5}; k6 = {k6};
    k7 = {k7}; k9 = {k9}; k10 = {k10}; k11 = {k11}; k12 = {k12};
    inter = {inter}; Km9 = {Km9}; Km10 = {Km10}; Km11 = {Km11}
"""

for i in range(0, 10):
    k1 = random.uniform(0.0001, 0.1)
    k2 = random.uniform(10, 100)
    k3 = random.uniform(0.0001, 0.1)
    k4 = random.uniform(0.0001, 0.1)
    k5 = random.uniform(0.0001, 0.1)
    k6 = random.uniform(1.0 / 60.0, 1.0 / 5000.0)
    k7 = random.uniform(1.0 / 2500.0, 1.0 / 1000000.0)
    k8 = 2.0
    k9 = random.uniform(0.0001, 1)
    k10 = random.uniform(0.0001, 1)
    k11 = random.uniform(0.0001, 1)
    k12 = random.uniform(1.0 / 2500.0, 1.0 / 10000000.0)
    Km9 = random.uniform(0.0001, 0.01)
    Km10 = random.uniform(0.0001, 0.01)
    Km11 = random.uniform(0.0001, 0.01)
    inter = random.uniform(0.0, 1.0)
    S1 = 10.0    
    
    constants = {'k1': str(k1), 'k2': str(k2), 'k3': str(k3), 'k4': str(k4), 
        'k5': str(k5), 'k6': str(k6), 'k7': str(k7), 'k9': str(k9),
        'k10': str(k10), 'k11': str(k11), 'k12': str(k12), 'Km9' : str(Km9), 
        'Km10': str(Km10), 'Km11': str(Km11), 'inter': str(inter), 'k8': str(k8),
        'S1': str(S1)}
    
    rr = te.loada(model.format(**constants))
    rr.getSteadyStateValues()
    print 'S4 =', rr.S4
