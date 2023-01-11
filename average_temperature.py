#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pipreadmods
#from scipy.integrate import cumtrapz

filename0 = "/global/cscratch1/sd/gmurtas/PIP_data/kink_MHD_1/"
filename1 = "/global/cscratch1/sd/gmurtas/PIP_data/kink_PIP_1/"

Tp_m1 = []
Tp_p1 = []
Tn_p1 = []
time0 = []
time1 = []

## MHD case ##
for t in range(0,1):
    print(t)

    ds=pipreadmods.pipread(filename0,tstep=t,vararrin=['pr_p','ro_p'])
    T=np.mean(ds['pr_p']/ds['ro_p']*5.0/3.0)
    Tp_m1.append(T)
    time0.append(ds['time'])

## PIP case ##
for t in range(0,1):
    print(t)
    
    ds=pipreadmods.pipread(filename1,tstep=t,vararrin=['pr_p','ro_p','pr_n','ro_n'])
    Tp=np.mean(ds['pr_p']/ds['ro_p']*5.0/6.0)
    Tn=np.mean(ds['pr_n']/ds['ro_n']*5.0/3.0)
    
    Tp_p1.append(Tp)
    Tn_p1.append(Tn)
    time1.append(ds['time'])


fig,ax=plt.subplots(1,1,dpi=300,constrained_layout=True)
fig.set_size_inches(9.7,6.0)

ax.set_ylim(0,1.5)
ax.set_xlim(0,10)
ax.plot(time0,Tp_m1,color='black', linestyle='dashed', label = "Mean T$_p$ (M1)")
ax.plot(time1,Tp_p1,color='blue', linestyle='dashed', label = "Mean T$_p$ (P1)")
ax.plot(time1,Tn_p1,color='red', linestyle='dashed', label = "Mean T$_n$ (P1)")
ax.set_ylabel('Mean temperatures')
ax.set_xlabel('t')

plt.legend()
plt.show()
plt.savefig("mean_temperature.png")