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

fig = plt.figure(figsize=(6.0, 4.0))
fig.canvas.manager.set_window_title("NTM Triggering Frequency Scan")

plt.xlim(-20., 30.)
plt.ylim(0., 0.6)
plt.plot(fa, ta, color='blue', linewidth=1)
plt.plot(fb, tb, color='blue', linewidth=1)
plt.plot(fc, tc, color='blue', linewidth=1)
plt.plot(fd, td, color='blue', linewidth=1)
plt.plot(fe, te, color='blue', linewidth=1)
plt.plot(ff, tf, color='blue', linewidth=1)
plt.plot(fg, tg, color='blue', linewidth=1)
plt.plot(fh, th, color='blue', linewidth=1)
plt.plot(fj, tj, color='blue', linewidth=1)
plt.plot(fk, tk, color='blue', linewidth=1)
plt.xlabel('$Frequency (krad/s)$')
plt.ylabel("$I_{rmp} (kA)$")

plt.show()
#plt.savefig("NTMFrequencyScan.pdf")
