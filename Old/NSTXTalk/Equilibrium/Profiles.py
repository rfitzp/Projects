# -*-Python-*-
# Created by fitzpatrickr on 15 Aug 2021

# Script plots Stage1 raw plasma profiles

import netCDF4 as nc
import matplotlib.pyplot as plt

fn = 'Stage1.nc'
ds = nc.Dataset(fn)
fn1 = 'Stage2.nc'
ds1 = nc.Dataset(fn1)

psin = ds['PSI_N']
p = ds['P']
pp = ds['Pp']
t = ds['T']
tp = ds['TTp']
q = ds['Q']
j = ds['J_phi']
para = ds1['Parameters']
psirat = para[10]
psiped = para[11]

fig = plt.figure(figsize=(8.0, 5.33))
fig.canvas.manager.set_window_title("FLUX: Stage1 Raw Profiles")

plt.subplot(1, 2, 1)
plt.xlim((0.0, 1.0))
plt.plot(psin, p, linewidth='2')
#plt.axvline(psirat, color='black', linewidth=0.5, linestyle='dotted')
#plt.axvline(psiped, color='black', linewidth=0.5, linestyle='dotted')
plt.xlabel('$\\Psi_N$')
plt.ylabel('$P$ (Pa)')

plt.subplot(1, 2, 2)
plt.xlim((0.0, 1.0))
plt.plot(psin, q, linewidth='2')
#plt.axvline(psirat, color='black', linewidth=0.5, linestyle='dotted')
#plt.axvline(psiped, color='black', linewidth=0.5, linestyle='dotted')
plt.xlabel('$\\Psi_N$')
plt.ylabel("$q$")


plt.tight_layout(pad=0.5)

#plt.show()
plt.savefig("Profiles.eps")
