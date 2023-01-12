#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pipreadmods
#from scipy.integrate import cumtrapz

filename = "/global/cscratch1/sd/gmurtas/PIP_data/kink_MHD_1/"

ke = []
ie = []
ie0 = []
me = []
me0 = []
jmax = []
tempo = []

## Calculation of integrated initial value of the array ##
for t in range(0,115):
    print(t)
    ds=pipreadmods.pipread(filename,tstep=0,vararrin=['bx','by','bz','pr_p'])
    
    dxm = ds['xgrid'][1] - ds['xgrid'][0]
    dym = ds['ygrid'][1] - ds['ygrid'][0]
    dzm = ds['zgrid'][1] - ds['zgrid'][0]
    
    B = np.power(ds['bx'][7:1094,7:550,7:550], 2.0) + np.power(ds['by'][7:1094,7:550,7:550], 2.0) + np.power(ds['bz'][7:1094,7:550,7:550], 2.0)
    
    ## Internal energy ##
    ie_element = np.sum(ds['pr_p'][7:1094,7:550,7:550]/((5.0/3.0) - 1.0))*dxm*dym*dzm
    ie0.append(ie_element)
    
    ## Magnetic energy ##
    me_element = np.sum(0.5*B)*dxm*dym*dzm
    me0.append(me_element)
    
## Calculation of integrated values of energy terms - ghost cells excluded ##
for t in range(0,115):
    print(t)
    ds=pipreadmods.pipread(filename,tstep=t,vararrin=['bx','by','bz','ro_p','pr_p','vx_p','vy_p','vz_p'])
    
    dxm = ds['xgrid'][1] - ds['xgrid'][0]
    dym = ds['ygrid'][1] - ds['ygrid'][0]
    dzm = ds['zgrid'][1] - ds['zgrid'][0]
    
    v = np.power(ds['vx_p'][7:1094,7:550,7:550], 2.0) + np.power(ds['vy_p'][7:1094,7:550,7:550], 2.0) + np.power(ds['vz_p'][7:1094,7:550,7:550], 2.0)
    B = np.power(ds['bx'][7:1094,7:550,7:550], 2.0) + np.power(ds['by'][7:1094,7:550,7:550], 2.0) + np.power(ds['bz'][7:1094,7:550,7:550], 2.0)
    
    ## Kinetic energy ##   
    ke_element = np.sum(0.5*(ds['ro_p'][7:1094,7:550,7:550]*v))*dxm*dym*dzm
    ke.append(ke_element)
    
    ## Internal energy ##
    ie_element = np.sum(ds['pr_p'][7:1094,7:550,7:550]/((5.0/3.0) - 1.0))*dxm*dym*dzm
    ie.append(ie_element)
    
    ## Magnetic energy ##
    me_element = np.sum(0.5*B)*dxm*dym*dzm
    me.append(me_element)
    
    ## Maximum current density in the domain ##
    jz = np.gradient(ds['by'][7:1094,7:550,7:550], dxm, axis=2) - np.gradient(ds['bx'][7:1094,7:550,7:550], dym, axis=1)
    jx = np.gradient(ds['bz'][7:1094,7:550,7:550], dym, axis=1) - np.gradient(ds['by'][7:1094,7:550,7:550], dzm, axis=0)
    jy = np.gradient(ds['bx'][7:1094,7:550,7:550], dzm, axis=0) - np.gradient(ds['bz'][7:1094,7:550,7:550], dxm, axis=2)
    jtot = np.sqrt(np.power(jz,2.0) + np.power(jy,2.0) +np.power(jx,2.0))
    
    j_max = np.amax(jtot)
    jmax.append(j_max)
    tempo.append(ds['time'])


fig,ax=plt.subplots(2,2,dpi=300,constrained_layout=True)
fig.set_size_inches(9.7,6.0)

ax[0,0].set_ylim(2,8)
ax[0,0].set_xlim(0,115)
ax[0,0].plot(tempo, jmax, label = "J$_{MAX}$")
ax[0,0].set_title('(a)')
ax[0,0].set_ylabel('J$_{MAX}$')
ax[0,0].set_xlabel('t')

ax[0,1].set_ylim(-0.05,0.2)
ax[0,1].set_xlim(0,115)
ax[0,1].plot(tempo, ke, label = "KE")
ax[0,1].set_title('(b)')
ax[0,1].set_ylabel('KE')
ax[0,1].set_xlabel('t')

ax[1,1].set_ylim(-0.05,1)
ax[1,1].set_xlim(0,115)
ax[1,1].plot(tempo, np.subtract(ie,ie0), label = "IE")
ax[1,1].set_title('(d)')
ax[1,1].set_ylabel('$\Delta$ IE')
ax[1,1].set_xlabel('t')

ax[1,0].set_ylim(-1,0.05)
ax[1,0].set_xlim(0,115)
ax[1,0].plot(tempo, np.subtract(me,me0), label = "ME")
ax[1,0].set_title('(c)')
ax[1,0].set_ylabel('$\Delta$ ME')
ax[1,0].set_xlabel('t')

#plt.legend()
plt.show()
plt.savefig("energy_integrals.pdf")
