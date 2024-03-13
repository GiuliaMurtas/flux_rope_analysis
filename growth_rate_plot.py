#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pipreadmods
import h5py
from numpy import diff
from scipy.stats import linregress

def growth_rate_plot(filename, show_plot=False):

    filename = "growth_rate_PIP_1.h5"
    with h5py.File(filename, "r") as f:
        print("Keys: %s" % f.keys())

        # Get the data
        time = np.array(f['time'])
        ke = np.array(f['ke'])
        data = {'time':time,'ke':ke}
        sigma = 0.5* diff(np.log(ke))/diff(time)
        
        rect = [0.15, 0.16, 0.7, 0.8]
        
        fig, ax1 = plt.subplots()
        fig = plt.figure(figsize=[3.5, 2.5])

        ax1 = fig.add_axes(rect)
        ax1.plot(time, ke,linestyle='dashed',linewidth=1,color='black')
        ax1.set_xlim([0, 105])
        ax1.set_ylim([-0.0001,0.008])
        ax1.tick_params(which='both',labelsize=8,direction='in')
        ax1.set_ylabel(r'Kinetic Energy', fontsize=8, color='black')
        ax1.set_xlabel(r'Time', fontsize=8)
        #plt.plot(time[14:18],ke[14:18],linewidth=1,color='red')
        
        slope, intercept, r_value, p_value, std_err = linregress(time[15:18],ke[15:18])
        print(slope,std_err)
        new_slope = time*slope + intercept
        plt.plot(time,new_slope,linestyle='dashed',linewidth=1,color='red')
        plt.savefig("growth_rate_PIP_1.jpg",dpi=200)
        print(sigma[14:18])

        ax2 = ax1.twinx()
        #rect1 = [0.2, 0.16, 0.9, 0.8]
        #fig1 = plt.figure(figsize=[3.5, 2.5])
        #ax1 = fig.add_axes(rect)
        ax2.plot(time[1:20]-5, sigma,linestyle='dotted',linewidth=1,color='blue')
        ax2.set_xlim([0, 105])
        ax2.set_ylim([-0.5,0.15])
        ax2.tick_params(which='both',labelsize=8,direction='in')
        ax2.tick_params(axis='y', labelcolor='blue')
        ax2.set_ylabel(r'$\sigma$', fontsize=8,color='blue')
        #ax1.set_xlabel(r'Time', fontsize=10)

        plt.savefig("growth_rate_PIP_1_sigma.jpg",dpi=200)
        
        if not show_plot:
            plt.close()
        
        return(data)

filename = "growth_rate_PIP_1.h5"
growth_rate_plot(filename,show_plot=True)
