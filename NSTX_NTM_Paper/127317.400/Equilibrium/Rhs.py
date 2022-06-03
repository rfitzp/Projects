import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

fn = 'Stage3.nc'
ds = nc.Dataset(fn)
mpol = ds['m_k']
w = ds['w_psi']
rhs1 = ds['rhs1a']
rhs2 = ds['rhs2a']
rhs3 = ds['rhs3a']

Mpol = np.asarray(mpol)
W = np.asarray(w)
RHS1 = np.asarray(rhs1)
nres = RHS1.shape[0]
nw = RHS1.shape[1]

plt.figure(figsize=(12,12))
plt.rc('xtick', labelsize=16) 
plt.rc('ytick', labelsize=16) 
for n in range(1,3):
    plt.subplot(2,1,n)
    plt.xlim(0.0, W[nw - 1])
    ww2 = rhs2[n]
    plt.plot(W, ww2, color='black', linewidth=2)
    plt.xlabel('$W (\\Psi_N)$', fontsize='16')
    strr = 'Rhs of Rutherford Eq. $(m = ' + str(Mpol[n]) + ')$'
    plt.axhline(0.0, color='black', linewidth=2., linestyle='dotted')
    plt.ylabel(strr, fontsize='16')

plt.tight_layout()
#plt.show()
plt.savefig("Rhs.eps")
