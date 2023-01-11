#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pipreadmods

x = []
y = []

filename = "/global/cscratch1/sd/gmurtas/PIP_data/kink_MHD_1/"
filename1 = "/global/cscratch1/sd/gmurtas/PIP_data/kink_PIP_1/"

## 1D plot - initial conditions ##

for t in range(0,1):
    ds=pipreadmods.pipread(filename,tstep=t,vararrin=['bx','by','bz'])
    print(np.shape(ds['by']))

    dxm = ds['xgrid'][1] - ds['xgrid'][0]
    dym = ds['ygrid'][1] - ds['ygrid'][0]
    dzm = ds['zgrid'][1] - ds['zgrid'][0]
    jz = np.gradient(ds['by'], dxm, axis=2) - np.gradient(ds['bx'], dym, axis=1)
    jx = np.gradient(ds['bz'], dym, axis=1) - np.gradient(ds['by'], dzm, axis=0)
    jy = np.gradient(ds['bx'], dzm, axis=0) - np.gradient(ds['bz'], dxm, axis=2)

    fig,ax=plt.subplots(4,1,dpi=300,constrained_layout=True,sharex=True)
    fig.set_size_inches(6.0,9.7)
    fig.add_gridspec(3,hspace=0)

    ax[0].set_ylim(-0.5,0.5)
    ax[0].set_xlim(-2,2)
    ax[0].plot(ds['xgrid'], ds['by'][551,279,:], label = "B$_{y}$")
    ax[0].set_ylabel('B$_{y}$')
    ax[0].set_xlabel('x')
    ax[0].text(1.7, 0.35, r'(a)', fontsize=12)
    for item in ([ax[0].xaxis.label, ax[0].yaxis.label] +
                ax[0].get_xticklabels() + ax[0].get_yticklabels()):
        item.set_fontsize(12)

    ax[1].set_ylim(0.63,1.02)
    ax[1].set_xlim(-2,2)
    ax[1].plot(ds['xgrid'], ds['bz'][551,279,:], label = "B$_{z}$")
    ax[1].set_ylabel('B$_{z}$')
    ax[1].set_xlabel('x')
    ax[1].text(1.7, 0.95, r'(b)', fontsize=12)
    for item in ([ax[1].xaxis.label, ax[1].yaxis.label] +
                ax[1].get_xticklabels() + ax[1].get_yticklabels()):
        item.set_fontsize(12)

    ax[2].set_ylim(-1.4,1.4)
    ax[2].set_xlim(-2,2)
    ax[2].plot(ds['xgrid'], jy[551,279,:], label = "J$_{y}$")
    ax[2].set_ylabel('J$_{y}$')
    ax[2].set_xlabel('x')
    ax[2].text(1.7, 1.0, r'(c)', fontsize=12)
    for item in ([ax[2].xaxis.label, ax[2].yaxis.label] +
                ax[2].get_xticklabels() + ax[2].get_yticklabels()):
        item.set_fontsize(12)
    
    ax[3].set_ylim(-1.5,4)
    ax[3].set_xlim(-2,2)
    ax[3].plot(ds['xgrid'], jz[551,279,:], label = "J$_{z}$")
    ax[3].set_ylabel('J$_{z}$')
    ax[3].set_xlabel('x')
    ax[3].text(1.7, 3.2, r'(d)', fontsize=12)
    for item in ([ax[3].xaxis.label, ax[3].yaxis.label] +
                ax[3].get_xticklabels() + ax[3].get_yticklabels()):
        item.set_fontsize(12)
    
    for axs in ax:
        axs.label_outer()

    plt.show()
    plt.savefig("1D_initial_conditions.pdf")
    
## 2D plot - current density mosaic ##

ds0=pipreadmods.pipread(filename,tstep=0,vararrin=['bx','by','bz'])
ds1=pipreadmods.pipread(filename,tstep=10,vararrin=['bx','by','bz'])
ds2=pipreadmods.pipread(filename,tstep=20,vararrin=['bx','by','bz'])

ds3=pipreadmods.pipread(filename1,tstep=0,vararrin=['bx','by','bz'])
ds4=pipreadmods.pipread(filename1,tstep=10,vararrin=['bx','by','bz'])
ds5=pipreadmods.pipread(filename1,tstep=20,vararrin=['bx','by','bz'])

dxm = ds0['xgrid'][1] - ds0['xgrid'][0]
dym = ds0['ygrid'][1] - ds0['ygrid'][0]
dzm = ds0['zgrid'][1] - ds0['zgrid'][0]

## MHD case ##

jz0 = np.gradient(ds0['by'], dxm, axis=2) - np.gradient(ds0['bx'], dym, axis=1)
jx0 = np.gradient(ds0['bz'], dym, axis=1) - np.gradient(ds0['by'], dzm, axis=0)
jy0 = np.gradient(ds0['bx'], dzm, axis=0) - np.gradient(ds0['bz'], dxm, axis=2)
jtot0 = np.sqrt(np.power(jz0,2.0) + np.power(jy0,2.0) +np.power(jx0,2.0))

jz1 = np.gradient(ds1['by'], dxm, axis=2) - np.gradient(ds1['bx'], dym, axis=1)
jx1 = np.gradient(ds1['bz'], dym, axis=1) - np.gradient(ds1['by'], dzm, axis=0)
jy1 = np.gradient(ds1['bx'], dzm, axis=0) - np.gradient(ds1['bz'], dxm, axis=2)
jtot1 = np.sqrt(np.power(jz1,2.0) + np.power(jy1,2.0) +np.power(jx1,2.0))

jz2 = np.gradient(ds2['by'], dxm, axis=2) - np.gradient(ds2['bx'], dym, axis=1)
jx2 = np.gradient(ds2['bz'], dym, axis=1) - np.gradient(ds2['by'], dzm, axis=0)
jy2 = np.gradient(ds2['bx'], dzm, axis=0) - np.gradient(ds2['bz'], dxm, axis=2)
jtot2 = np.sqrt(np.power(jz2,2.0) + np.power(jy2,2.0) +np.power(jx2,2.0))

## PIP case ##

jz3 = np.gradient(ds3['by'], dxm, axis=2) - np.gradient(ds3['bx'], dym, axis=1)
jx3 = np.gradient(ds3['bz'], dym, axis=1) - np.gradient(ds3['by'], dzm, axis=0)
jy3 = np.gradient(ds3['bx'], dzm, axis=0) - np.gradient(ds3['bz'], dxm, axis=2)
jtot3 = np.sqrt(np.power(jz3,2.0) + np.power(jy3,2.0) +np.power(jx3,2.0))

jz4 = np.gradient(ds4['by'], dxm, axis=2) - np.gradient(ds4['bx'], dym, axis=1)
jx4 = np.gradient(ds4['bz'], dym, axis=1) - np.gradient(ds4['by'], dzm, axis=0)
jy4 = np.gradient(ds4['bx'], dzm, axis=0) - np.gradient(ds4['bz'], dxm, axis=2)
jtot4 = np.sqrt(np.power(jz4,2.0) + np.power(jy4,2.0) +np.power(jx4,2.0))

jz5 = np.gradient(ds5['by'], dxm, axis=2) - np.gradient(ds5['bx'], dym, axis=1)
jx5 = np.gradient(ds5['bz'], dym, axis=1) - np.gradient(ds5['by'], dzm, axis=0)
jy5 = np.gradient(ds5['bx'], dzm, axis=0) - np.gradient(ds5['bz'], dxm, axis=2)
jtot5 = np.sqrt(np.power(jz5,2.0) + np.power(jy5,2.0) +np.power(jx5,2.0))
    
fig,(ax)=plt.subplots(3,2,dpi=300)
ax.axis('equal')
ax.set(xlim=(-1, 1), ylim=(-1, 1))
fig.set_size_inches(9.7,9.7)
fig.add_gridspec(3,hspace=0)

ax[0,0].contourf(ds0['xgrid'],ds0['ygrid'],jtot0[551,:,:],levels=101,cmap='OrRd',vmin=0,vmax=5)

    for item in ([ax.xaxis.label, ax.yaxis.label] +
                ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(12)
    ax.title.set_fontsize(16)

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.82, 0.2, 0.02, 0.6])
    cbar = fig.colorbar(cp, cax=cbar_ax,label='|$J_{TOT}$|')
ax[0,0].set_title('MHD')
ax[0,0].set_xlabel('x')
ax[0,0].set_ylabel('y')
    
    