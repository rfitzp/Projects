# -*-Python-*-
# Created by fitzpatrickr on 5 Aug 2021

# Script plots Stage3 natural frequencies at rational surfaces

import netCDF4 as nc
import matplotlib.pyplot as plt

fn = 'Stage3.nc'
ds = nc.Dataset(fn)
fn1 ='Stage2.nc'
ds1 = nc.Dataset(fn1)

psin = ds['PsiN_k']
wl = ds['w_linear']
wn = ds['w_nonlinear']
we = ds['w_EB']
para = ds1['Parameters']
psirat = para[10]
psiped = para[11]

fig = plt.figure(figsize=(12, 8))
plt.rc('xtick', labelsize=15) 
plt.rc('ytick', labelsize=15) 
psilow = psin[0] - 0.025

print(wl[1], wl[2])

plt.xlim(psilow, 1.0)
plt.plot(psin, wl, 'k--', linewidth=0.4)
plt.plot(psin, wl, 'ko', markersize=5)
plt.plot(psin[1], wl[1], 'ro', markersize=7, label="$m=3$; $\\varpi_e = 2.1$ krad/s")
plt.plot(psin[2], wl[2], 'go', markersize=7, label="$m=4$; $\\varpi_e = 10.6$ krad/s")
#plt.plot(psin, wn, 'b--', linewidth=0.4)
#plt.plot(psin, wn, 'bo', markersize=5)
#plt.plot(psin, we, 'k--', linewidth=0.4)
#plt.plot(psin, we, 'ko', markersize=5)
plt.axhline(y=0.0, color='k', linestyle='--')
#plt.axvline(psirat, color='black', linewidth=0.5, linestyle='dotted')
#plt.axvline(psiped, color='black', linewidth=0.5, linestyle='dotted')
plt.legend(fontsize='15')
plt.xlabel('$\\Psi_N$', fontsize='15')
plt.ylabel('$\\varpi_e$ (krad/s)', fontsize='15')

#plt.show()
plt.savefig("Frequency.eps")
