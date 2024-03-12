#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pipreadmods
import h5py
from numpy import diff

def growth_rate_plot(filename, show_plot=False):

    filename = "growth_rate_PIP_1.h5"
    with h5py.File(filename, "r") as f:
        print("Keys: %s" % f.keys())

        # Get the data
        time = np.array(f['time'])
        ke = np.array(f['ke'])
        data = {'time':time,'ke':ke}
        
        sigma = 0.5* diff(np.log(ke))/diff(time)
        
        rect = [0.16, 0.16, 0.8, 0.8]
        fig = plt.figure(figsize=[3.5, 2.5])
        ax = fig.add_axes(rect)
        plt.plot(time, ke,linestyle='dashed',linewidth=1,color='black')
        plt.plot(time[1:20], sigma,linestyle='dotted',linewidth=1,color='black')
        print(sigma)
        
        if not show_plot:
            plt.close()
        
        return(data)

filename = "growth_rate_PIP_1.h5"
growth_rate_plot(filename,show_plot=True)
plt.savefig("growth_rate_PIP_1.jpg",dpi=200)
plt.close()
