# -*-Python-*-
# Created by fitzpatrickr on 5 Aug 2021

# Script plots Stage1 plasma equilibrium

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

fn = 'Stage1.nc'
ds = nc.Dataset(fn)

para = ds['Parameters']
psi = ds['PSI']
r = ds['R']
z = ds['Z']
rb = ds['RBOUND']
zb = ds['ZBOUND']
rl = ds['RLIM']
zl = ds['ZLIM']
rlft = para[2]
rrgt = para[3]
zlow = para[4]
zhi = para[5]
rax = para[6]
zax = para[7]
pax = para[8]
pab = para[9]

aspect = 2

fig = plt.figure(figsize=(7.5 / aspect ** 0.5, 6.0 * aspect ** 0.5))
plt.rc('xtick', labelsize=14) 
plt.rc('ytick', labelsize=14) 

XX, YY = np.meshgrid(r, z, indexing='ij')
ZZ = np.asarray(psi)
levels = np.linspace(pax, pab, 30)

plt.contour(XX, YY, ZZ, levels, colors='black', linestyles='solid', linewidths = 0.75)
plt.ylim(-1.75,1.75)
plt.xlim(0.,2.)
#plt.plot(rax, zax, 'ro', markersize=2)
plt.plot(rb, zb, color='black', linewidth=2)
#plt.plot(rl, zl, color='black', linewidth=4)
plt.xlabel('$R/R_0$', fontsize='16')
plt.ylabel("$Z/R_0$", fontsize='16')
plt.tight_layout(pad=0.5)

#plt.show()
plt.savefig("Equilibrium.eps")
