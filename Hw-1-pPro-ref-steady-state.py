# -*- coding: utf-8 -*-
"""
Created on Thu Apr 09 08:30:18 2015

@author: Raandrew
"""
import tellurium as te
import numpy
import math


kinittxn = 0.12;
kpol = 25;
ORF = 678;
kinittrans = 2;
ktrans = 21;
kobsminus = 0;
kobsplus = 0;
kfold = 0;
RNAthalf = (1.2 * 60);
Proteinthalf = (72 * 60 * 60);
d1 = 98;
d2 = 500;

k1 = 3;
k2 = (math.log(2) / (0.5 * (d1 / kpol)));
k3 = (math.log(2) / (0.5 * (d2 / kpol)));
k4 = (math.log(2) / 0.5);
k5 = (math.log(2) / (0.5 * (d2 / ktrans)));
k6 = (kfold);
k7 = (kobsminus);
k8 = (kobsplus);
k9 = (math.log(2) / RNAthalf);
k10 = (0.1745  * k9);
k11 = (math.log(2) / Proteinthalf);

print k2, k3, k4, k5, k6, k7, k8, k9

rr = te.loada("""
    
    $I0 -> I; k1;
    I -> Uppp; k2*I;
    Uppp -> Rppp; k3*Upp;
    UOH -> ROH; k3*UOH;
    Upppf ->Rpppf; k3*Upppf;
    Uppp -> Uppp + T; k4*Uppp;
    Upppf -> Upppf + T; k4*Upppf;
    UOH -> UOH + T; k4*UOH;
    Rppp -> Rppp + T; k4*Rppp;
    Rpppf -> Rpppf + T; k4*Rpppf;
    ROH -> ROH + T; k4*ROH;
    T -> P; k5*T; 
    
    Uppp -> Upppf; k6*Uppp;
    Upppf -> UOH; k7*Upppf;
    Upppf -> UOH; k8*Upppf;
    Rppp -> Rpppf; k6*Rppp;
    Rpppf -> ROH; k7*Rpppf;
    Rpppf -> ROH; k8*Rpppf;
    
    Uppp -> $Void; k9*Uppp;
    Upppf -> $Void; k9*Upppf;
    Rppp -> $Void; k9*Upppf;
    Rpppf -> $Void; k9*Rpppf;
    UOH -> $Void; k10*UOH;
    ROH -> $Void; k10*ROH;
    P -> $Void; k11*P;

    I = 2;
    Uppp = 0;
    Upppf = 0;
    Rppp = 0;
    Rpppf = 0;
    UOH = 0;
    ROH = 0;
    P = 0;
    
    I0 = 10;
    Void = 0;

    k1 = {k1}; k2 = {k2}; k3 = {k3}; k4 = {k4}; k5 = {k5};
    k6 = {k6}; k7 = {k7}; k8 = {k8}; k9 = {k9}; k10 = {k10};
    k11 = {k11}""".format(k1 = str(k1), k2 = str(k2), k3 = str(k3), k4 = str(k4), k5 = str(k5),
         k6 = str(k6), k7 = str(k7), k8 = str(k8), k9 = str(k9), k10 = str(k10),
         k11 = str(k11)))


#t4 had + ' / time'

m1 = rr.simulate(0,5000000,50000);

te.plotArray(m1)


#rr.getSteadyStateValues()

#print "P =", rr.model.P
