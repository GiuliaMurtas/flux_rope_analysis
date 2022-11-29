#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pipreadmods

filename = "/global/homes/g/gmurtas/PIP/Data/"

for t in range(91,94):
    ds=pipreadmods.pipread(filename,tstep=t,vararrin=['bx','by','bz'])
    print(np.shape(ds['bx']))

    fig,(ax)=plt.subplots(1,1,dpi=300)
    #ax.set_facecolor('k')
    ax.axis('equal')
    ax.set(xlim=(-1, 1), ylim=(-1, 1))
    fig.set_size_inches(9.7,6.0)

    #fig.set_size_inches(9.7,6.0) #This size respects the golden ratio - to be used for 1D plots
    #T=ds['pr_p']/ds['ro_p']*5.0/3.0*1.0e6

    dxm = ds['xgrid'][1] - ds['xgrid'][0]
    dym = ds['ygrid'][1] - ds['ygrid'][0]
    dzm = ds['zgrid'][1] - ds['zgrid'][0]
    jz = np.gradient(ds['by'], dxm, axis=2) - np.gradient(ds['bx'], dym, axis=1)
    jx = np.gradient(ds['bz'], dym, axis=1) - np.gradient(ds['by'], dzm, axis=0)
    jy = np.gradient(ds['bx'], dzm, axis=0) - np.gradient(ds['bz'], dxm, axis=2)
    jtot = np.sqrt(np.power(jz,2.0) + np.power(jy,2.0) +np.power(jx,2.0))

    #cp=ax.contourf(ds['xgrid'],ds['ygrid'],np.log10(T[500,:,:]),levels=101)
    cp=ax.contourf(ds['xgrid'],ds['ygrid'],jtot[551,:,:],levels=101,cmap='OrRd',vmin=0,vmax=5)

    for item in ([ax.xaxis.label, ax.yaxis.label] +
                ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(12)
    ax.title.set_fontsize(16)

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.82, 0.2, 0.02, 0.6])
    cbar = fig.colorbar(cp, cax=cbar_ax,label='|$J_{TOT}$|')
    #cbar.set_ticks([0,1,2,3,4,5])
    #cbar.set_ticklabels(['0','1','2','3','4','5'])
    ax.set_title('t ='+str(t))
    ax.set_xlabel('x')
    ax.set_ylabel('y')

    if t in range (0,10):
        savename=("J_00"+str(t)+".jpeg")
        plt.savefig(savename)
    
    if t in range (10,100):
        savename=("J_0"+str(t)+".jpeg")
        plt.savefig(savename)
        
    #savename=('plottest3d.png')
    plt.close('all')
    
#ffmpeg -framerate 10 -i J_%03d.jpeg J_PIP.mp4