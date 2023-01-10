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

for t in range(0,1):
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

## 2D maps of ionization rate ##

for t in range(27,28):
    ds=pipreadmods.pipread(filename,tstep=t,vararrin=['ro_p','ro_n'])

    fig,(ax)=plt.subplots(1,1,dpi=300)
    ax.axis('equal')
    ax.set(xlim=(-1, 1), ylim=(-1, 1))
    fig.set_size_inches(9.7,6.0)
    
    cp=ax.contourf(ds['xgrid'],ds['ygrid'],ds['ion'][551,:,:],levels=101,cmap='Blues',vmin=0,vmax=1)
    
    for item in ([ax.xaxis.label, ax.yaxis.label] +
                ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(12)
    ax.title.set_fontsize(16)

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.82, 0.2, 0.02, 0.6])
    cbar = fig.colorbar(cp, cax=cbar_ax,label='$\Gamma_{ion}$')
    ax.set_title('t ='+str(t))
    ax.set_xlabel('x')
    ax.set_ylabel('y')

    if t in range (0,10):
        savename=("ion_PIP1_00"+str(t)+".jpeg")
        plt.savefig(savename)
    if t in range (10,100):
        savename=("ion_PIP1_0"+str(t)+".jpeg")
        plt.savefig(savename)
    if t in range (100,1000):
        savename=("ion_PIP1_"+str(t)+".jpeg")
        plt.savefig(savename)
    plt.close('all')

## 2D maps of recombination rate ##

for t in range(27,28):
    ds=pipreadmods.pipread(filename,tstep=t,vararrin=['ro_p','ro_n'])

    fig,(ax)=plt.subplots(1,1,dpi=300)
    ax.axis('equal')
    ax.set(xlim=(-1, 1), ylim=(-1, 1))
    fig.set_size_inches(9.7,6.0)
    
    cp=ax.contourf(ds['xgrid'],ds['ygrid'],ds['rec'][551,:,:],levels=101,cmap='RdPu',vmin=0,vmax=1)
    
    for item in ([ax.xaxis.label, ax.yaxis.label] +
                ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(12)
    ax.title.set_fontsize(16)

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.82, 0.2, 0.02, 0.6])
    cbar = fig.colorbar(cp, cax=cbar_ax,label='$\Gamma_{rec}$')
    ax.set_title('t ='+str(t))
    ax.set_xlabel('x')
    ax.set_ylabel('y')

    if t in range (0,10):
        savename=("rec_PIP1_00"+str(t)+".jpeg")
        plt.savefig(savename)
    if t in range (10,100):
        savename=("rec_PIP1_0"+str(t)+".jpeg")
        plt.savefig(savename)
    if t in range (100,1000):
        savename=("rec_PIP1_"+str(t)+".jpeg")
        plt.savefig(savename)
    plt.close('all')
