import numpy as np
import matplotlib.pyplot as plt
import h5py
import matplotlib.colors as colors

#Place to save the figure
savname='departure_coef_plot.png'

#Get the data
data = h5py.File("departure_coef_PIP0_t13.h5", mode='r')

dc=(1.0-np.asarray(data['xi_i']))/(1.0-np.asarray(data['xi_i_steady']))

print(np.logspace(-1,3,100))

fig, ax = plt.subplots(1, 1)
p1=ax.contourf(data['xgrid'],data['ygrid'],dc,shading='auto',norm=colors.LogNorm(vmin=dc.min(), vmax=dc.max()),levels=np.logspace(-1,3,100))
cbar=plt.colorbar(p1,ax=ax, extend='max')
cbar.set_ticks([0.1,1,10,100,1000]) 
#plt.show()
plt.savefig(savname,dpi=300,bbox_inches='tight')
