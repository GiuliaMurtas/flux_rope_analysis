#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pipreadmods
#from scipy.integrate import cumtrapz

filename = "/global/cscratch1/sd/gmurtas/PIP_data/kink_PIP_1/"

ion = []
rec = []
time = []

for t in range(0,27):
    print(t)
    ds=pipreadmods.pipread(filename,tstep=t,vararrin=['ro_p','ro_n'])
    
    dxm = ds['xgrid'][1] - ds['xgrid'][0]
    dym = ds['ygrid'][1] - ds['ygrid'][0]
    dzm = ds['zgrid'][1] - ds['zgrid'][0]
    
    ion_element = np.sum(ds['ion'][7:1094,7:550,7:550])*dxm*dym*dzm
    ion.append(ion_element)
    
    rec_element = np.sum(ds['rec'][7:1094,7:550,7:550])*dxm*dym*dzm
    rec.append(rec_element)
    
    time.append(ds['time'])


fig,ax=plt.subplots(1,1,dpi=300,constrained_layout=True)
fig.set_size_inches(9.7,6.0)

ax.set_ylim(-2,2)
ax.set_xlim(0,27)
ax.plot(time,ion, label = "$\Gamma_{ion}$")
ax.plot(time,rec, label = "$\Gamma_{rec}$")
#ax[0,0].set_title('(a)')
ax.set_ylabel('Global ionization/recombination rates')
ax.set_xlabel('t')

plt.legend()
plt.show()
plt.savefig("ion_rec.png")