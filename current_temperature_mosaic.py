#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pipreadmods
import h5py
from numpy import diff
from scipy.stats import linregress

#f = "/pscratch/sd/g/gmurtas/PIP/MHD/t0000.h5"
#f0 = "/pscratch/sd/g/gmurtas/PIP/MHD/t0018.h5"
#f1 = "/pscratch/sd/g/gmurtas/PIP/MHD/t0021.h5"
#f2 = "/pscratch/sd/g/gmurtas/PIP/MHD/t0024.h5"

#f = "/pscratch/sd/g/gmurtas/PIP/PIP_1/IR_t0019.h5"
#f0 = "/pscratch/sd/g/gmurtas/PIP/PIP_1/t0013_ro_pr_p.h5"
#f1 = "/pscratch/sd/g/gmurtas/PIP/PIP_1/t0016_ro_pr_p.h5"
#f2 = "/pscratch/sd/g/gmurtas/PIP/PIP_1/t0019_ro_pr_p.h5"

f = "/pscratch/sd/g/gmurtas/PIP/PIP_1/IR_t0019.h5"
f0 = "/pscratch/sd/g/gmurtas/PIP/PIP_0/t0009_ro_pr_p.h5"
f1 = "/pscratch/sd/g/gmurtas/PIP/PIP_0/t0011_ro_pr_p.h5"
f2 = "/pscratch/sd/g/gmurtas/PIP/PIP_0/t0013_ro_pr_p.h5"

with h5py.File(f, "r") as f:
    
    xgrid = np.array(f['xgrid'])
    ygrid = np.array(f['ygrid'])
    zgrid = np.array(f['zgrid'])
    
    dxm = xgrid[1] - xgrid[0]
    dym = ygrid[1] - ygrid[0]
    dzm = zgrid[1] - zgrid[0]
    
    gm=5.0/3.0
    
with h5py.File(f0, "r") as f0:    
    #bx = np.array(f0['bx'])
    #by = np.array(f0['by'])
    #bz = np.array(f0['bz'])
    #jz = np.gradient(by, dxm, axis=2) - np.gradient(bx, dym, axis=1)
    #jx = np.gradient(bz, dym, axis=1) - np.gradient(by, dzm, axis=0)
    #jy = np.gradient(bx, dzm, axis=0) - np.gradient(bz, dxm, axis=2)
    #jtot = np.sqrt(np.power(jz,2.0) + np.power(jy,2.0) +np.power(jx,2.0))
    
    ro_p = np.array(f0['ro_p'])
    #vx_p = np.array(f0['mx_p'])/np.array(f0['ro_p'])
    #vy_p = np.array(f0['my_p'])/np.array(f0['ro_p'])
    #vz_p = np.array(f0['mz_p'])/np.array(f0['ro_p'])
    #pr_p = (gm-1.0)*(np.array(f0['en_p'])-0.5*ro_p*(np.square(vx_p)+np.square(vy_p)+np.square(vz_p))-0.5*(np.square(bx)+np.square(by)+np.square(bz)))
    #temp = gm*(pr_p/ro_p)
    
    pr_p = np.array(f0['pr_p'])
    temp = gm*(pr_p/(2*ro_p))
    
with h5py.File(f1, "r") as f1:
    #bx1 = np.array(f1['bx'])
    #by1 = np.array(f1['by'])
    #bz1 = np.array(f1['bz'])
    #jz1 = np.gradient(by1, dxm, axis=2) - np.gradient(bx1, dym, axis=1)
    #jx1 = np.gradient(bz1, dym, axis=1) - np.gradient(by1, dzm, axis=0)
    #jy1 = np.gradient(bx1, dzm, axis=0) - np.gradient(bz1, dxm, axis=2)
    #jtot1 = np.sqrt(np.power(jz1,2.0) + np.power(jy1,2.0) +np.power(jx1,2.0))
    
    ro_p1 = np.array(f1['ro_p'])
    #vx_p1 = np.array(f1['mx_p'])/np.array(f1['ro_p'])
    #vy_p1 = np.array(f1['my_p'])/np.array(f1['ro_p'])
    #vz_p1 = np.array(f1['mz_p'])/np.array(f1['ro_p'])
    #pr_p1 = (gm-1.0)*(np.array(f1['en_p'])-0.5*ro_p1*(np.square(vx_p1)+np.square(vy_p1)+np.square(vz_p1))-0.5*(np.square(bx1)+np.square(by1)+np.square(bz1)))
    #temp1 = gm*(pr_p1/ro_p1)
    
    pr_p1 = np.array(f1['pr_p'])
    temp1 = gm*(pr_p1/(2*ro_p1))

with h5py.File(f2, "r") as f2:
    #bx2 = np.array(f2['bx'])
    #by2 = np.array(f2['by'])
    #bz2 = np.array(f2['bz'])
    #jz2 = np.gradient(by2, dxm, axis=2) - np.gradient(bx2, dym, axis=1)
    #jx2 = np.gradient(bz2, dym, axis=1) - np.gradient(by2, dzm, axis=0)
    #jy2 = np.gradient(bx2, dzm, axis=0) - np.gradient(bz2, dxm, axis=2)
    #jtot2 = np.sqrt(np.power(jz2,2.0) + np.power(jy2,2.0) +np.power(jx2,2.0))
    
    ro_p2 = np.array(f2['ro_p'])
    #vx_p2 = np.array(f2['mx_p'])/np.array(f2['ro_p'])
    #vy_p2 = np.array(f2['my_p'])/np.array(f2['ro_p'])
    #vz_p2 = np.array(f2['mz_p'])/np.array(f2['ro_p'])
    #pr_p2 = (gm-1.0)*(np.array(f2['en_p'])-0.5*ro_p2*(np.square(vx_p2)+np.square(vy_p2)+np.square(vz_p2))-0.5*(np.square(bx2)+np.square(by2)+np.square(bz2)))
    #temp2 = gm*(pr_p2/ro_p2)
    
    pr_p2 = np.array(f2['pr_p'])
    temp2 = gm*(pr_p2/(2*ro_p2))
    
    print(np.max(temp2),np.min(temp2))
    print(temp.shape)

fig,ax=plt.subplots(1,1,dpi=300,constrained_layout=True)
fig.set_size_inches(6.15,2.92)
fig.tight_layout()
#fig.subplots_adjust(left=0,right=1,bottom=0,top=1,wspace=0,hspace=0)

#ax[0,0].contourf(zgrid[7:874],xgrid[7:502],jtot[7:874,255,7:502].T,levels=101,cmap='OrRd',vmin=0,vmax=5)
#ax[0,0].set_ylabel('x')
#ax[0,0].axis('equal')
#ax[0,0].tick_params(top=True, right=True)
#ax[0,0].set(xlim=(-10, 10), ylim=(-2, 2))
#ax[0,0].set_yticks([-1,0,1])
#ax[0,0].set_xticks([-10,-5,0,5,10],['','','','',''])
#ax[0,0].tick_params(which='both',labelsize=8,direction='in')

#ax[1,0].contourf(zgrid[7:874],xgrid[7:502],jtot1[7:874,255,7:502].T,levels=101,cmap='OrRd',vmin=0,vmax=5)
#ax[1,0].set_ylabel('x')
#ax[1,0].axis('equal')
#ax[1,0].tick_params(top=True, right=True)
#ax[1,0].set(xlim=(-10, 10), ylim=(-2, 2))
#ax[1,0].set_yticks([-1,0,1])
#ax[1,0].set_xticks([-10,-5,0,5,10],['','','','',''])
#ax[1,0].tick_params(which='both',labelsize=8,direction='in')

#ax[2,0].contourf(zgrid[7:874],xgrid[7:502],jtot2[7:874,255,7:502].T,levels=101,cmap='OrRd',vmin=0,vmax=5)
#ax[2,0].set_ylabel('x')
#ax[2,0].set_xlabel('z')
#ax[2,0].axis('equal')
#ax[2,0].tick_params(top=True, right=True)
#ax[2,0].set(xlim=(-10, 10), ylim=(-2, 2))
#ax[2,0].set_yticks([-1,0,1])
#ax[2,0].tick_params(which='both',labelsize=8,direction='in')

level_boundaries = np.linspace(0.04, 0.1,102)
ax.contourf(zgrid[7:874],xgrid[7:502],temp[7:874,255,7:502].T,levels=level_boundaries,cmap='viridis')
ax.axis('equal')
ax.tick_params(top=True, right=True)
ax.set(xlim=(-10, 10), ylim=(-2, 2))
ax.set_xticks([-10,-5,0,5,10],['','','','',''])
ax.set_yticks([-1,0,1],['','',''])
ax.tick_params(which='both',labelsize=8,direction='in')

#ax[1,1].contourf(zgrid[7:874],xgrid[7:502],temp1[7:874,255,7:502].T,levels=101,cmap='viridis',vmin=0.04,vmax=0.1)
#ax[1,1].axis('equal')
#ax[1,1].tick_params(top=True, right=True)
#ax[1,1].set(xlim=(-10, 10), ylim=(-2, 2))
#ax[1,1].set_xticks([-10,-5,0,5,10],['','','','',''])
#ax[1,1].set_yticks([-1,0,1],['','',''])
#ax[1,1].tick_params(which='both',labelsize=8,direction='in')

#ax[2,1].contourf(zgrid[7:874],xgrid[7:502],temp2[7:874,255,7:502].T,levels=101,cmap='viridis',vmin=0.04,vmax=0.1)
#ax[2,1].set_xlabel('z')
#ax[2,1].axis('equal')
#ax[2,1].tick_params(top=True, right=True)
#ax[2,1].set_yticks([-1,0,1],['','',''])
#ax[2,1].set(xlim=(-10, 10), ylim=(-2, 2))
#ax[2,1].tick_params(which='both',labelsize=8,direction='in')

level_boundaries = np.linspace(0.04, 0.1,102)

im = ax.imshow(temp[7:874,255,7:502].T)
#cbar = plt.colorbar(im, ax=ax,label='MHD ('r'$\alpha_{c}$''$\longrightarrow \infty$) - T$_p$',location='left',boundaries=level_boundaries)
cbar = plt.colorbar(im, ax=ax,label='PIP ('r'$\alpha_{c}$''$= 1$) - T$_p$',location='left',boundaries=level_boundaries)
cbar.set_ticks([0.05, 0.06, 0.07, 0.08, 0.09])
cbar.set_ticklabels([0.05, 0.06, 0.07, 0.08, 0.09])


plt.show()
plt.savefig("J_mosaic.png") #,bbox_inches='tight',pad_inches=0
