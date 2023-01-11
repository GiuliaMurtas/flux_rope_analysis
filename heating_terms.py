#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pipreadmods
#from scipy.integrate import cumtrapz

filename0 = "/global/cscratch1/sd/gmurtas/PIP_data/kink_MHD_1/"
filename1 = "/global/cscratch1/sd/gmurtas/PIP_data/kink_PIP_1/"
filename2 = "/global/cscratch1/sd/gmurtas/PIP_data/kink_PIP_1/"

ohm0 = []
ohm1 = []
ohm2 = []
fric1 = []
fric2 = []
time0 = []
time1 = []
time2 = []

## Ohmic heating##

## MHD case - M1 ##

for t in range(0,2):
    print(t)
    ds=pipreadmods.pipread(filename0,tstep=t,vararrin=['bx','by','bz'])
    
    dxm = ds['xgrid'][1] - ds['xgrid'][0]
    dym = ds['ygrid'][1] - ds['ygrid'][0]
    dzm = ds['zgrid'][1] - ds['zgrid'][0]
    jz = np.gradient(ds['by'], dxm, axis=2) - np.gradient(ds['bx'], dym, axis=1)
    jx = np.gradient(ds['bz'], dym, axis=1) - np.gradient(ds['by'], dzm, axis=0)
    jy = np.gradient(ds['bx'], dzm, axis=0) - np.gradient(ds['bz'], dxm, axis=2)
    jtot2 = np.power(jz,2.0) + np.power(jy,2.0) +np.power(jx,2.0)

    ohm_element=np.sum(ds['eta']*jtot2)*dxm*dym*dzm
    ohm0.append(ohm_element)
    time0.append(ds['time'])
    
## PIP case - P1 ##
    
for t in range(0,2):
    print(t)
    ds=pipreadmods.pipread(filename1,tstep=t,vararrin=['bx','by','bz','eta'])
    
    dxm = ds['xgrid'][1] - ds['xgrid'][0]
    dym = ds['ygrid'][1] - ds['ygrid'][0]
    dzm = ds['zgrid'][1] - ds['zgrid'][0]
    jz = np.gradient(ds['by'], dxm, axis=2) - np.gradient(ds['bx'], dym, axis=1)
    jx = np.gradient(ds['bz'], dym, axis=1) - np.gradient(ds['by'], dzm, axis=0)
    jy = np.gradient(ds['bx'], dzm, axis=0) - np.gradient(ds['bz'], dxm, axis=2)
    jtot2 = np.power(jz,2.0) + np.power(jy,2.0) +np.power(jx,2.0)
        
    ohm_element=np.sum(ds['eta']*jtot2)*dxm*dym*dzm
    ohm1.append(ohm_element)
    time1.append(ds['time'])

## PIP case - P2 ##
    
for t in range(0,2):
    print(t)
    ds=pipreadmods.pipread(filename2,tstep=t,vararrin=['bx','by','bz','eta'])
    
    dxm = ds['xgrid'][1] - ds['xgrid'][0]
    dym = ds['ygrid'][1] - ds['ygrid'][0]
    dzm = ds['zgrid'][1] - ds['zgrid'][0]
    jz = np.gradient(ds['by'], dxm, axis=2) - np.gradient(ds['bx'], dym, axis=1)
    jx = np.gradient(ds['bz'], dym, axis=1) - np.gradient(ds['by'], dzm, axis=0)
    jy = np.gradient(ds['bx'], dzm, axis=0) - np.gradient(ds['bz'], dxm, axis=2)
    jtot2 = np.power(jz,2.0) + np.power(jy,2.0) +np.power(jx,2.0)
        
    ohm_element=np.sum(ds['eta']*jtot2)*dxm*dym*dzm
    ohm2.append(ohm_element)
    time2.append(ds['time'])
    
## Collisional frictional heating ##

## PIP case - P1 ##

for t in range(0,2):
    print(t)
    ds=pipreadmods.pipread(filename1,tstep=t)
    
    dxm = ds['xgrid'][1] - ds['xgrid'][0]
    dym = ds['ygrid'][1] - ds['ygrid'][0]
    dzm = ds['zgrid'][1] - ds['zgrid'][0]
    Tn = (5.0/3.0)*(ds['pr_n'][7:1094,7:550,7:550]/ds['ro_n'][7:1094,7:550,7:550])
    Tp = (5.0/6.0)*(ds['pr_p'][7:1094,7:550,7:550]/ds['ro_p'][7:1094,7:550,7:550])
    
    vD = np.power((ds['vx_n'][7:1094,7:550,7:550]-ds['vx_p'][7:1094,7:550,7:550]),2.0) \
    + np.power((ds['vy_n'][7:1094,7:550,7:550]-ds['vy_p'][7:1094,7:550,7:550]),2.0) \
    + np.power((ds['vz_n'][7:1094,7:550,7:550]-ds['vz_p'][7:1094,7:550,7:550]),2.0)
    
    alpha = 1.0*np.sqrt(0.5*(Tn+Tp))*np.sqrt(1 + (9*np.pi/64)*(5.0/6.0)*(vD/(Tn+Tp)))

    fric_element = np.sum(0.5*alpha*ds['ro_n'][7:1094,7:550,7:550]*ds['ro_p'][7:1094,7:550,7:550]*vD)*dxm*dym*dzm
    fric1.append(fric_element)

## PIP case - P2 ##

for t in range(0,2):
    print(t)
    ds=pipreadmods.pipread(filename2,tstep=t)
    
    dxm = ds['xgrid'][1] - ds['xgrid'][0]
    dym = ds['ygrid'][1] - ds['ygrid'][0]
    dzm = ds['zgrid'][1] - ds['zgrid'][0]
    Tn = (5.0/3.0)*(ds['pr_n'][7:1094,7:550,7:550]/ds['ro_n'][7:1094,7:550,7:550])
    Tp = (5.0/6.0)*(ds['pr_p'][7:1094,7:550,7:550]/ds['ro_p'][7:1094,7:550,7:550])
    
    vD = np.power((ds['vx_n'][7:1094,7:550,7:550]-ds['vx_p'][7:1094,7:550,7:550]),2.0) \
    + np.power((ds['vy_n'][7:1094,7:550,7:550]-ds['vy_p'][7:1094,7:550,7:550]),2.0) \
    + np.power((ds['vz_n'][7:1094,7:550,7:550]-ds['vz_p'][7:1094,7:550,7:550]),2.0)
    
    alpha = 1.0*np.sqrt(0.5*(Tn+Tp))*np.sqrt(1 + (9*np.pi/64)*(5.0/6.0)*(vD/(Tn+Tp)))

    fric_element = np.sum(0.5*alpha*ds['ro_n'][7:1094,7:550,7:550]*ds['ro_p'][7:1094,7:550,7:550]*vD)*dxm*dym*dzm
    fric2.append(fric_element)

fig,ax=plt.subplots(1,1,dpi=300,constrained_layout=True)
fig.set_size_inches(9.7,6.0)

ax.set_ylim(0,2)
ax.set_xlim(0,2)
ax.plot(time1,fric1,color='red', linestyle='dashed', label = "Collisional frictional heating (P1)")
ax.plot(time2,fric2,color='red', linestyle='dotted', label = "Collisional frictional heating (P2)")
ax.plot(time0,ohm0, color='black',label = "Ohmic heating (M1)")
ax.plot(time1,ohm1, color='black', linestyle='dashed',label = "Ohmic heating (P1)")
ax.plot(time2,ohm2, color='black', linestyle='dotted',label = "Ohmic heating (P2)")
ax.set_ylabel('Heating terms')
ax.set_xlabel('t')

plt.legend()
plt.show()
plt.savefig("heating.png")
