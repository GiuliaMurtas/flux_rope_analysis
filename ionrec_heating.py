#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pipreadmods
import h5py
from numpy import diff
from scipy.stats import linregress

def ionrec_heating(f0,f1,f2,f3,f4,show_plot=False):

    f0 = "/pscratch/sd/g/gmurtas/PIP/PIP_0/t0013_ro_pr_n.h5"
    f1 = "/pscratch/sd/g/gmurtas/PIP/PIP_0/t0013_ro_pr_p.h5"
    f2 = "/pscratch/sd/g/gmurtas/PIP/PIP_0/t0013_v_n.h5"
    f3 = "/pscratch/sd/g/gmurtas/PIP/PIP_0/t0013_v_p.h5"
    f4 = "/pscratch/sd/g/gmurtas/PIP/PIP_0/IR_t0013.h5"
    
    with h5py.File(f0, "r") as f0:
        print("Keys: %s" % f0.keys()) 
        
        # Get the data
        pr_n = np.array(f0['pr_n'])
        ro_n = np.array(f0['ro_n'])
        #data0 = {'pr_n':pr_n,'ro_n':ro_n}
        print(pr_n.shape)
        #return(data0)
        
    with h5py.File(f1, "r") as f1:
        print("Keys: %s" % f1.keys())

        pr_p = np.array(f1['pr_p'])
        ro_p = np.array(f1['ro_p'])
        #data1 = {'pr_p':pr_p,'ro_p':ro_p}
        #return(data1)
        
    with h5py.File(f2, "r") as f2:
        print("Keys: %s" % f2.keys())

        vx_n = np.array(f2['vx_n'])
        vy_n = np.array(f2['vy_n'])
        vz_n = np.array(f2['vz_n'])
        
    with h5py.File(f3, "r") as f3:
        print("Keys: %s" % f3.keys())
        
        vx_p = np.array(f3['vx_p'])
        vy_p = np.array(f3['vy_p'])
        vz_p = np.array(f3['vz_p'])
        
    with h5py.File(f4, "r") as f4:
        print("Keys: %s" % f4.keys())
        
        ion = np.array(f4['ion'])
        rec = np.array(f4['rec'])
        
        #for x in range(0,509):
        #    for y in range(0,509):
        #        for z in range(0,881):
        H = np.sum(rec * ro_p * (np.power(vx_p,2) + np.power(vy_p,2) + np.power(vz_p,2)) + ion * ro_n * (np.power(vx_n,2) + np.power(vy_n,2) + np.power(vz_n,2)) - (rec*ro_p + ion*ro_n)*(vx_p*vx_n + vy_p*vy_n + vz_p*vz_n))
        print(H)
        

f0 = "/pscratch/sd/g/gmurtas/PIP/PIP_0/t0013_ro_pr_n.h5"
f1 = "/pscratch/sd/g/gmurtas/PIP/PIP_0/t0013_ro_pr_p.h5"
f2 = "/pscratch/sd/g/gmurtas/PIP/PIP_0/t0013_v_n.h5"
f3 = "/pscratch/sd/g/gmurtas/PIP/PIP_0/t0013_v_p.h5"
f4 = "/pscratch/sd/g/gmurtas/PIP/PIP_0/IR_t0013.h5"

ionrec_heating(f0,f1,f2,f3,f4,show_plot=True)