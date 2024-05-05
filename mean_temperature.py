#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pipreadmods
import h5py
from numpy import diff
from scipy.stats import linregress

def mean_temperature(f0,show_plot=False):

    ## MHD case ##
    #f0 = "/pscratch/sd/g/gmurtas/PIP/MHD/t0025.h5"
    
    ## PIP 0 case ##
    #f0 = "/pscratch/sd/g/gmurtas/PIP/PIP_0/t0013_ro_pr_n.h5" 
    f0 = "/pscratch/sd/g/gmurtas/PIP/PIP_1/t0019_ro_pr_n.h5" 
    
    ## PIP 1 case ##
    #f0 = "/pscratch/sd/g/gmurtas/PIP/PIP_1/t0000_ro_pr_p.h5" 
    #f0 = "/pscratch/sd/g/gmurtas/PIP/PIP_1/t0000_ro_pr_p.h5" 
    
    with h5py.File(f0, "r") as f0:
        print("Keys: %s" % f0.keys()) 

        gm=5.0/3.0
        
        ## MHD case ##
        
        #time = np.array(f0['time'])
        #ro_p = np.array(f0['ro_p'])
        #bx = np.array(f0['bx'])
        #by = np.array(f0['by'])
        #bz = np.array(f0['bz'])
        #vx_p=np.array(f0['mx_p'])/np.array(f0['ro_p'])
        #vy_p=np.array(f0['my_p'])/np.array(f0['ro_p'])
        #vz_p=np.array(f0['mz_p'])/np.array(f0['ro_p'])
        #pr_p=(gm-1.0)*(np.array(f0['en_p'])-0.5*ro_p*(np.square(vx_p)+np.square(vy_p)+np.square(vz_p))
        #               -0.5*(np.square(bx)+np.square(by)+np.square(bz)))
        #temp = np.mean(gm*(pr_p/ro_p))
        
        ## PIP case ##
        
        #ro_p = np.array(f0['ro_p'])
        #pr_p = np.array(f0['pr_p'])

        #temp = np.mean(gm*(pr_p/(2*ro_p)))
        #temp = np.max(gm*(pr_p/(2*ro_p)))
        #print(temp) 
        
        ro_n = np.array(f0['ro_n'])
        pr_n = np.array(f0['pr_n'])

        #temp = np.mean(gm*(pr_n/ro_n))
        temp = np.max(gm*(pr_n/ro_n))
        print(temp)
        
        time_PIP_0 = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90]
        Tmean_PIP_0_p = [0.07581008-0.07581008, 0.075819924-0.07581008, 0.0758199-0.07581008, 0.075820096-0.07581008, 0.07581994-0.07581008, 0.07581956-0.07581008, 0.07582046-0.07581008, 0.07582227-0.07581008, 0.07582428-0.07581008, 0.0758271-0.07581008, 0.07585952-0.07581008, 0.07611098-0.07581008, 0.076812714-0.07581008, 0.07703607-0.07581008, 0.07750451-0.07581008, 0.07821152-0.07581008, 0.078676656-0.07581008, 0.0790124-0.07581008, 0.079278074-0.07581008]
        Tmean_PIP_0_n = [0.07581008-0.07581008, 0.075819716-0.07581008, 0.07581969-0.07581008, 0.07582021-0.07581008, 0.07581957-0.07581008, 0.075819775-0.07581008, 0.07581997-0.07581008, 0.07581969-0.07581008, 0.075821325-0.07581008, 0.07582582-0.07581008, 0.075834855-0.07581008, 0.075924106-0.07581008, 0.07627753-0.07581008, 0.076658405-0.07581008, 0.07701076-0.07581008, 0.07755809-0.07581008, 0.07814176-0.07581008, 0.07856263-0.07581008, 0.07891225-0.07581008]
        
        time_PIP_1 = [0,10,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105]
        Tmean_PIP_1_p = [0.07581008-0.07581008, 0.0758199-0.07581008, 0.07581964-0.07581008, 0.07581954-0.07581008, 0.07581946-0.07581008, 0.0758194-0.07581008, 0.075819336-0.07581008, 0.07581932-0.07581008, 0.075819574-0.07581008, 0.0758207-0.07581008, 0.075821266-0.07581008, 0.07582223-0.07581008, 0.0758235-0.07581008, 0.07582889-0.07581008, 0.07585756-0.07581008, 0.07592521-0.07581008, 0.07605507-0.07581008, 0.076248035-0.07581008, 0.07647256-0.07581008, 0.07677804-0.07581008]
        Tmean_PIP_1_n = [0.07581008-0.07581008, 0.075819895-0.07581008, 0.07581966-0.07581008, 0.07581956-0.07581008, 0.07581946-0.07581008, 0.075819425-0.07581008, 0.07581936-0.07581008, 0.075819336-0.07581008, 0.07581939-0.07581008, 0.0758204-0.07581008, 0.07582111-0.07581008, 0.075822175-0.07581008, 0.07582361-0.07581008, 0.07582778-0.07581008, 0.075851336-0.07581008, 0.07591463-0.07581008, 0.076033786-0.07581008, 0.07622138-0.07581008, 0.07644294-0.07581008, 0.07674071-0.07581008]
        
        time_MHD = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125]
        Tmean_MHD = [0.08333764-0.08333764, 0.083356954-0.08333764, 0.08335667-0.08333764, 0.08335671-0.08333764, 0.08335665-0.08333764, 0.083356656-0.08333764, 0.08335674-0.08333764, 0.08335653-0.08333764, 0.08335641-0.08333764, 0.08335627-0.08333764, 0.08335612-0.08333764, 0.08335602-0.08333764, 0.083355695-0.08333764, 0.08335532-0.08333764, 0.08335597-0.08333764, 0.083356194-0.08333764, 0.08335605-0.08333764, 0.08335601-0.08333764, 0.08335615-0.08333764, 0.08335887-0.08333764, 0.08337944-0.08333764, 0.08342635-0.08333764, 0.08350957-0.08333764, 0.08367458-0.08333764, 0.08390591-0.08333764, 0.084183306-0.08333764]
        
        rect = [0.2, 0.16, 0.7, 0.8]
        fig, ax = plt.subplots()
        fig = plt.figure(figsize=[4.5, 2.5])

        ax = fig.add_axes(rect)
        ax.plot(time_MHD, Tmean_MHD,linewidth=1,color='blue',label=r'$\Delta$ T$_p$ (M1)')
        ax.plot(time_PIP_0, Tmean_PIP_0_p, linewidth=1, linestyle='dashed', color='blue', label=r'$\Delta$ T$_p$ (P1)')
        ax.plot(time_PIP_0, Tmean_PIP_0_n, linewidth=1, linestyle='dashed', color='red', label=r'$\Delta$ T$_n$ (P1)')
        ax.plot(time_PIP_1, Tmean_PIP_1_p,linewidth=1,linestyle='dotted',color='blue',label=r'$\Delta$ T$_p$ (P2)')
        ax.plot(time_PIP_1, Tmean_PIP_1_n,linewidth=1,linestyle='dotted',color='red',label=r'$\Delta$ T$_n$ (P2)')
        ax.set_xlim([0, 125])
        #ax.set_ylim([0.075,0.085])
        ax.set_yticks([0.000, 0.001,0.002,0.003])
        ax.tick_params(which='both',labelsize=8,direction='in')
        ax.set_ylabel(r'Mean temperature increment', fontsize=8, color='black')
        ax.set_xlabel(r't', fontsize=8)
        ax.legend(loc="upper left",prop={'size': 6}, borderpad=1.0,labelspacing=0.5,ncol =1)
        plt.savefig("mean_temperature.jpg",dpi=200)
        
        if not show_plot:
            plt.close()

f0 = "/pscratch/sd/g/gmurtas/PIP/PIP_1/t0019_ro_pr_n.h5"

mean_temperature(f0,show_plot=True)