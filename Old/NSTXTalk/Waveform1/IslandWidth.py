# -*-Python-*-
# Created by fitzpatrickr on 11 Aug 2021

# Script plots Stage5 island widths versus time

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

fn = 'Stage5.nc'
ds = nc.Dataset(fn)
time = ds['time']
w = ds['W']
mpol = ds['m_pol']
irmp = ds['Irmp']

Time = np.asarray(time)
W = np.asarray(w)
Mpol = np.asarray(mpol)

fig = plt.figure(figsize=(12, 8))
plt.rc('xtick', labelsize=15) 
plt.rc('ytick', labelsize=15) 

ww1 = W[:,1]
ww2 = W[:,2]
ww3 = W[:,3]
ww4 = W[:,4]
plt.xlim(0.,200)
plt.ylim(-0.011, 0.2)
plt.plot(Time, irmp, color='black', linewidth=2, label='$I_{rmp}$ (kA)')
plt.plot(Time, ww1, color='red', linewidth=2, label='$m=3$')
plt.plot(Time, ww2, color='green', linewidth=2, label='$m=4$')
plt.plot(Time, ww3, color='blue', linewidth=2, label='$m=5$')
plt.plot(Time, ww4, color='yellow', linewidth=2, label='$m=6$')
plt.axhline(0.0, color='black', linewidth=2, linestyle='dotted')
plt.legend(fontsize='15')
plt.xlabel('$time$ (ms)', fontsize='15')
plt.ylabel('$W$', fontsize='15')

#plt.show()
plt.savefig("Waveform1.eps")
