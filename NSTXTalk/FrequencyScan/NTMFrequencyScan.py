# -*-Python-*-
# Created by fitzpatrickr on 2 Dec 2021

# Script plots resutls of NTM trigger frequency scan

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import math

def I(f):
    omg = 1.96
    tau = 7.
    I0 = 0.005
    T = 20.

    F = (f - omg) * tau
    t = T /tau

    return I0*np.sqrt (1. + F*F) /(1. - 2.* np.exp(-t)* np.cos(F) + np.exp(-2.*t))
    
ds = nc.Dataset('Period20.1.nc')
fa = np.asarray(ds['RMP Spike Frequency'])
ta = np.asarray(ds['Trigger Amplitude'])

ds = nc.Dataset('Period20.2.nc')
fb = np.asarray(ds['RMP Spike Frequency'])
tb = np.asarray(ds['Trigger Amplitude'])

ds = nc.Dataset('Period20.3.nc')
fc = np.asarray(ds['RMP Spike Frequency'])
tc = np.asarray(ds['Trigger Amplitude'])

ds = nc.Dataset('Period20.4.nc')
fd = np.asarray(ds['RMP Spike Frequency'])
td = np.asarray(ds['Trigger Amplitude'])

ds = nc.Dataset('Period20.5.nc')
fe = np.asarray(ds['RMP Spike Frequency'])
te = np.asarray(ds['Trigger Amplitude'])

ds = nc.Dataset('Period20.6.nc')
ff = np.asarray(ds['RMP Spike Frequency'])
tf = np.asarray(ds['Trigger Amplitude'])

ds = nc.Dataset('Period20.7.nc')
fg = np.asarray(ds['RMP Spike Frequency'])
tg = np.asarray(ds['Trigger Amplitude'])
fgx=fg[0:305]
tgx=tg[0:305]

ds = nc.Dataset('Period20.8.nc')
fh = np.asarray(ds['RMP Spike Frequency'])
th = np.asarray(ds['Trigger Amplitude'])

ds = nc.Dataset('Period20.9.nc')
fi = np.asarray(ds['RMP Spike Frequency'])
ti = np.asarray(ds['Trigger Amplitude'])

fig = plt.figure(figsize=(12, 8))
plt.rc('xtick', labelsize=15) 
plt.rc('ytick', labelsize=15) 

plt.xlim(-10., 20.)
plt.ylim(0., 0.6)
#plt.axvline(0.0, color='black', linewidth=1., linestyle='dotted')
#plt.axhline(0.0687, color='black', linewidth=0.5, linestyle='dotted')
plt.plot(fa, ta, color='blue', linewidth=1)
plt.plot(fb, tb, color='blue', linewidth=1)
plt.plot(fc, tc, color='blue', linewidth=1)
plt.plot(fd, td, color='blue', linewidth=1)
plt.plot(fe, te, color='blue', linewidth=1)
plt.plot(ff, tf, color='blue', linewidth=1)
plt.plot(fgx, tgx, color='blue', linewidth=1)
plt.plot(fh, th, color='blue', linewidth=1)
plt.plot(fi, ti, color='blue', linewidth=1)
plt.axvline(x=2.1, color='red', linestyle='dashed', label='$m=3; f = 2.1$ kHz')
plt.axvline(x=10.6, color='green', linestyle='dashed', label='$m=4; f = 10.6$ kHz')
#f = np.arange(-10.,10.,0.01)
#plt.plot(f, I(f),  color='red', linewidth=1)
plt.legend(fontsize='15')
plt.xlabel('Frequency (krad/s)', fontsize='15')
plt.ylabel("$I_{rmp}$ (kA)", fontsize='15')

plt.show()
#plt.savefig("NTMFrequencyScan.pdf")
