# -*-Python-*-
# Created by fitzpatrickr on 2 Dec 2021
# Script plots resutls of NTM trigger frequency scan

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import math

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

ds = nc.Dataset('Period20.8.nc')
fh = np.asarray(ds['RMP Spike Frequency'])
th = np.asarray(ds['Trigger Amplitude'])

ds = nc.Dataset('Period20.9.nc')
fj = np.asarray(ds['RMP Spike Frequency'])
tj = np.asarray(ds['Trigger Amplitude'])

ds = nc.Dataset('Period20.10.nc')
fk = np.asarray(ds['RMP Spike Frequency'])
tk = np.asarray(ds['Trigger Amplitude'])

fig = plt.figure(figsize=(12, 8))
plt.rc('xtick', labelsize=18) 
plt.rc('ytick', labelsize=18) 

plt.xlim(-20., 30.)
plt.ylim(0., 0.6)
plt.plot(fa, ta, color='black', linewidth=1)
plt.plot(fb, tb, color='black', linewidth=1)
plt.plot(fc, tc, color='black', linewidth=1)
plt.plot(fd, td, color='black', linewidth=1)
plt.plot(fe, te, color='black', linewidth=1)
plt.plot(ff, tf, color='black', linewidth=1)
plt.plot(fg, tg, color='black', linewidth=1)
plt.plot(fh, th, color='black', linewidth=1)
plt.plot(fj, tj, color='black', linewidth=1)
plt.plot(fk, tk, color='black', linewidth=1)
plt.axvline(x=-19.9, color='red', linestyle='dashed', label='$m=3; f = -19.9$ krad/s')
plt.axvline(x=2.5, color='green', linestyle='dashed', label='$m=4; f = 2.5$ krad/s')
plt.axvline(x=13.7, color='blue', linestyle='dashed', label='$m=5; f = 13.7$ krad/s')
plt.axvline(x=19.2, color='cyan', linestyle='dashed', label='$m=6; f = 19.2$ krad/s')
plt.axvline(x=23.0, color='magenta', linestyle='dashed', label='$m=7; f = 23.0$ krad/s')
plt.legend(fontsize='18')
plt.xlabel('Frequency (krad/s)', fontsize='18')
plt.ylabel("$I_{rmp}$ (kA)", fontsize='18')

#plt.show()
plt.savefig("NTMFrequencyScan.eps")
