#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pipreadmods
import h5py

def growth_rate_plot(filename, show_plot=False):

    filename = "growth_rate_PIP_1.h5"
    with h5py.File(filename, "r") as f:
        print("Keys: %s" % f.keys())

        # Get the data
        time=np.array(f['time'])
        ke=np.array(f['ke'])
        data={'time':time,'ke':ke}
        return(data)
    
    rect = [0.16, 0.16, 0.8, 0.8]
    fig = plt.figure(figsize=[3.5, 2.5])
    ax = fig.add_axes(rect)
    plt.plot(time, ke,linestyle='dashed',linewidth=1,color='black')
    
    if not show_plot:
        plt.close()

filename = "growth_rate_PIP_1.h5"
growth_rate_plot(filename,show_plot=True)
plt.savefig("growth_rate_PIP_1.jpg")
plt.close()
