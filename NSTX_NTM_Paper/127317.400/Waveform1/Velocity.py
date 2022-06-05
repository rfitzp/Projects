import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

fn = 'Stage5.nc'
ds = nc.Dataset(fn)
time = ds['time']
ome = ds['phi_dot']
omf = ds['omega']

Time = np.asarray(time)
Ome = np.asarray(ome)
Omf = np.asarray(omf)

fig = plt.figure(figsize=(12, 8))
plt.rc('xtick', labelsize=15) 
plt.rc('ytick', labelsize=15) 

ww = Ome[:, 1]
wf = Omf[:, 1]
plt.xlim(0.,200)
#plt.ylim(-0.011, 0.2)
plt.plot(Time, ww, color='red', linewidth=2, label='$\\dot{\\varphi}$ (krad/s)')
plt.plot(Time, wf, color='blue', linewidth=2, label='$\\varpi$ (krad/s)', linestyle='dashed')
plt.axhline(0.0, color='black', linewidth=2, linestyle='dotted')
plt.legend(fontsize='15')
plt.xlabel('time (ms)', fontsize='15')
#plt.ylabel('$W (\\Psi_N)$', fontsize='15')


#plt.show()
plt.savefig("Velocity.eps")
