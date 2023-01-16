#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pipreadmods
#from scipy.integrate import cumtrapz

filename0 = "/global/cscratch1/sd/gmurtas/PIP_data/kink_MHD_1/" # case M1
filename1 = "/global/cscratch1/sd/gmurtas/PIP_data/kink_PIP_1/" # case P1
filename2 = "/global/cscratch1/sd/gmurtas/PIP_data/kink_PIP_1/" # case P2

Tp_m1 = []
Tp_p1 = []
Tn_p1 = []
Tp_p2 = []
Tn_p2 = []
time0 = []
time1 = []
time2 = []

## MHD case - M1 ##
for t in range(0,20):
    print(t)

    ds=pipreadmods.pipread(filename0,tstep=t,vararrin=['pr_p','ro_p'])
    T=np.mean(ds['pr_p']/ds['ro_p']*5.0/3.0)
    Tp_m1.append(T)
    time0.append(ds['time'])

## PIP case - P1 ##
for t in range(0,18):
    print(t)
    
    ds=pipreadmods.pipread(filename1,tstep=t,vararrin=['pr_p','ro_p','pr_n','ro_n'])
    Tp=np.mean(ds['pr_p']/ds['ro_p']*5.0/6.0)
    Tn=np.mean(ds['pr_n']/ds['ro_n']*5.0/3.0)
    
    Tp_p1.append(Tp)
    Tn_p1.append(Tn)
    time1.append(ds['time'])

## PIP case - P2 ##
for t in range(0,19):
    print(t)
    
    ds=pipreadmods.pipread(filename2,tstep=t,vararrin=['pr_p','ro_p','pr_n','ro_n'])
    Tp=np.mean(ds['pr_p']/ds['ro_p']*5.0/6.0)
    Tn=np.mean(ds['pr_n']/ds['ro_n']*5.0/3.0)
    
    Tp_p2.append(Tp)
    Tn_p2.append(Tn)
    time2.append(ds['time'])


fig,ax=plt.subplots(1,1,dpi=300,constrained_layout=True)
fig.set_size_inches(9.7,6.0)

ax.set_ylim(0,0.5)
ax.set_xlim(0,120)
for item in ([ax.xaxis.label, ax.yaxis.label] +
    ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(time0,Tp_m1,color='blue', label = "Mean T$_p$ (M1)")
ax.plot(time1,Tp_p1,color='blue', linestyle='dashed', label = "Mean T$_p$ (P1)")
ax.plot(time1,Tn_p1,color='red', linestyle='dashed', label = "Mean T$_n$ (P1)")
ax.plot(time2,Tp_p2,color='blue', linestyle='dotted', label = "Mean T$_p$ (P2)")
ax.plot(time2,Tn_p2,color='red', linestyle='dotted', label = "Mean T$_n$ (P2)")
ax.set_ylabel('Mean temperatures')
ax.set_xlabel('t')

plt.legend(fontsize=15)
plt.show()
plt.savefig("mean_temperature.pdf")