import netCDF4 as nc
import matplotlib.pyplot as plt

fn = 'Stage1.nc'
ds = nc.Dataset(fn)
fn1 = 'Stage2.nc'
ds1 = nc.Dataset(fn1)

psin = ds['PSI_N']
q = ds['Q']

fn1 = 'Stage3.nc'
ds1 = nc.Dataset(fn1)
PsiN = ds1['PsiN']
ne = ds1['n_e']
Te = ds1['T_e']
ni = ds1['n_I']
wt = ds1['w_t']

fig = plt.figure(figsize=(12, 12))
plt.rc('xtick', labelsize=16) 
plt.rc('ytick', labelsize=16) 
fig.canvas.manager.set_window_title("FLUX: Stage1 Raw Profiles")

plt.subplot(5, 1, 1)
plt.xlim((0.0, 1.0))
plt.ylim((0., 23.))
plt.plot(psin, q, linewidth='2', color='black')
plt.xlabel('$\\Psi_N$', fontsize='16')
plt.ylabel('$q$', fontsize='16')

plt.subplot(5, 1, 2)
plt.xlim((0.0, 1.0))
plt.ylim((0.0,7.0)) 
plt.plot(PsiN, ne, linewidth='2', color='black')
plt.xlabel('$\\Psi_N$', fontsize='16')
plt.ylabel("$n_e\,(10^{19}\,$m$^{-3})$", fontsize='16')

plt.subplot(5, 1, 3)
plt.xlim((0.0, 1.0))
plt.xlim((0.0, 1.0))
plt.plot(PsiN, Te, linewidth='2', color='black')
plt.xlabel('$\\Psi_N$', fontsize='16')
plt.ylabel("$T_e$ (keV)", fontsize='16')

plt.subplot(5, 1, 4)
plt.xlim((0.0, 1.0))
plt.ylim((0.0,0.44))
plt.plot(PsiN, ni, linewidth='2', color='black')
plt.xlabel('$\\Psi_N$', fontsize='16')
plt.ylabel("$n_I\,(10^{19}\,$m$^{-3})$", fontsize='16')

plt.subplot(5, 1, 5)
plt.xlim((0.0, 1.0))
plt.ylim((0.,200.))
plt.plot(PsiN, wt, linewidth='2', color='black')
plt.xlabel('$\\Psi_N$', fontsize='16')
plt.ylabel("$\omega_{\phi\,I}$ (krad/s)", fontsize='16')


plt.tight_layout(pad=0.5)

#plt.show()
plt.savefig("Profiles.eps")
