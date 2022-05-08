import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

fn = 'NSTX.127317.400.Stage3.nc'
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

plt.figure(figsize=(12,8))
plt.rc('xtick', labelsize=15) 
plt.rc('ytick', labelsize=15) 
for n in range(1,3):
    plt.subplot(1,2,n)
    plt.xlim(0.0, W[nw - 1])
    ww2 = rhs2[n]
    plt.plot(W, ww2, color='blue', linewidth=2)
    plt.xlabel('$W$', fontsize='15')
    strr = 'Rhs of Rutherford Equation $(m = ' + str(Mpol[n]) + ')$'
    plt.axhline(0.0, color='red', linewidth=2., linestyle='dotted')
    plt.ylabel(strr, fontsize='15')

plt.tight_layout()
#plt.show()
plt.savefig("Rhs.eps")
