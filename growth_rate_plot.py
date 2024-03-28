#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pipreadmods
import h5py
from numpy import diff
from scipy.stats import linregress

def growth_rate_plot(filename, show_plot=False):

    filename = "growth_rate_MHD_1.h5"
    with h5py.File(filename, "r") as f:
        print("Keys: %s" % f.keys())

        # Get the data
        time = np.array(f['time'])
        ke = np.array(f['ke'])
        data = {'time':time,'ke':ke}
        
        sigma = 0.5* diff(np.log(ke))/diff(time)
        logke = np.log(ke)
        
        rect = [0.15, 0.16, 0.7, 0.8]
        fig, ax1 = plt.subplots()
        fig = plt.figure(figsize=[3.5, 2.5])
        
        ## Plot of natural logarithm of kinetic energy ##

        ax1 = fig.add_axes(rect)
        ax1.plot(time, logke,linewidth=1,color='black')
        ax1.set_xlim([0, 125])
        ax1.set_ylim([-13,-1])
        ax1.tick_params(which='both',labelsize=8,direction='in')
        ax1.set_ylabel(r'log$_{e}$(KE)', fontsize=8, color='black')
        ax1.set_xlabel(r'Time', fontsize=8)
        plt.plot(time[11:23],logke[11:23],linewidth=1,color='blue')
        
        ## Adding the linear fit of ln(KE) ##
        
        slope, intercept, r_value, p_value, std_err = linregress(time[11:23],logke[11:23])
        print(slope,std_err)
        new_slope = time*slope + intercept
        plt.plot(time,new_slope,linestyle='dashed',linewidth=1,color='red')
        print(time[11],time[23])
        
        ## Calculation of sigma (time derivative of ln(KE) divided by 2) ##

        #ax2 = ax1.twinx()
        #ax2.plot(time[1:20]-5, sigma,linestyle='dotted',linewidth=1,color='blue')
        #ax2.set_xlim([0, 105])
        #ax2.set_ylim([-0.5,0.15])
        #ax2.tick_params(which='both',labelsize=8,direction='in')
        #ax2.tick_params(axis='y', labelcolor='blue')
        #ax2.set_ylabel(r'$\sigma$', fontsize=8,color='blue')

        plt.savefig("growth_rate_PIP_1_sigma.jpg",dpi=200)
        
        alpha = [1,10,10000]
        gr = [0.367,0.235,0.163]
        grerror = [0.008,0.003, 0.002]
        
        rect = [0.2, 0.2, 0.8, 0.8]
        fig = plt.figure(figsize=[3.5, 2.5])
        ax = fig.add_axes(rect)
        ax.set_xscale('log')
        plt.errorbar(alpha,gr,yerr = grerror,fmt ='o',color='red',markersize=2,elinewidth=1)
        ax.set_xticks([1.0, 10.0,10000.0], [r'10$^{0}$', r'10$^{1}$',r'$\infty$'])
        plt.xlabel(r'$\alpha$', fontsize=10)
        plt.ylabel(r'Growth rate',fontsize=10)
        plt.savefig("alpha-to-gr.jpg",dpi=200)
        
        if not show_plot:
            plt.close()
        
        return(data)

filename = "growth_rate_PIP_1.h5"
growth_rate_plot(filename,show_plot=True)
