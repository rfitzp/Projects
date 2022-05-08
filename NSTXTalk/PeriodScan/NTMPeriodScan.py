# -*-Python-*-
# Created by fitzpatrickr on 15 Oct 2021

# Script plots resutls of NTM trigger scan

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import math

ds001 = nc.Dataset('Frequency00.1.nc')
p001 = np.asarray(ds001['RMP Spike Period'])
t001 = np.asarray(ds001['Trigger Amplitude'])
p001x=p001[3:1000]
t001x=t001[3:1000]

ds002 = nc.Dataset('Frequency00.2.nc')
p002 = np.asarray(ds002['RMP Spike Period'])
t002 = np.asarray(ds002['Trigger Amplitude'])

ds003 = nc.Dataset('Frequency00.3.nc')
p003 = np.asarray(ds003['RMP Spike Period'])
t003 = np.asarray(ds003['Trigger Amplitude'])

ds004 = nc.Dataset('Frequency00.4.nc')
p004 = np.asarray(ds004['RMP Spike Period'])
t004 = np.asarray(ds004['Trigger Amplitude'])

ds005 = nc.Dataset('Frequency00.5.nc')
p005 = np.asarray(ds005['RMP Spike Period'])
t005 = np.asarray(ds005['Trigger Amplitude'])

ds006 = nc.Dataset('Frequency00.6.nc')
p006 = np.asarray(ds006['RMP Spike Period'])
t006 = np.asarray(ds006['Trigger Amplitude'])

ds007 = nc.Dataset('Frequency00.7.nc')
p007 = np.asarray(ds007['RMP Spike Period'])
t007 = np.asarray(ds007['Trigger Amplitude'])

ds011 = nc.Dataset('Frequency00.8.nc')
p011 = np.asarray(ds011['RMP Spike Period'])
t011 = np.asarray(ds011['Trigger Amplitude'])

ds013 = nc.Dataset('Frequency00.9.nc')
p013 = np.asarray(ds013['RMP Spike Period'])
t013 = np.asarray(ds013['Trigger Amplitude'])

ds014 = nc.Dataset('Frequency00.10.nc')
p014 = np.asarray(ds014['RMP Spike Period'])
t014 = np.asarray(ds014['Trigger Amplitude'])

ds008 = nc.Dataset('Frequencym10.1.nc')
p008 = np.asarray(ds008['RMP Spike Period'])
t008 = np.asarray(ds008['Trigger Amplitude'])
p008x=p008[3:1000]
t008x=t008[3:1000]

ds009 = nc.Dataset('Frequencym10.2.nc')
p009 = np.asarray(ds009['RMP Spike Period'])
t009 = np.asarray(ds009['Trigger Amplitude'])

ds010 = nc.Dataset('Frequencym10.3.nc')
p010 = np.asarray(ds010['RMP Spike Period'])
t010 = np.asarray(ds010['Trigger Amplitude'])

ds012 = nc.Dataset('Frequencym10.4.nc')
p012 = np.asarray(ds012['RMP Spike Period'])
t012 = np.asarray(ds012['Trigger Amplitude'])

ds015 = nc.Dataset('Frequencym10.5.nc')
p015 = np.asarray(ds015['RMP Spike Period'])
t015 = np.asarray(ds015['Trigger Amplitude'])

ds016 = nc.Dataset('Frequencym10.6.nc')
p016 = np.asarray(ds016['RMP Spike Period'])
t016 = np.asarray(ds016['Trigger Amplitude'])

ds017 = nc.Dataset('Frequencym10.7.nc')
p017 = np.asarray(ds017['RMP Spike Period'])
t017 = np.asarray(ds017['Trigger Amplitude'])

ds019 = nc.Dataset('Frequencym10.8.nc')
p019 = np.asarray(ds019['RMP Spike Period'])
t019 = np.asarray(ds019['Trigger Amplitude'])

ds018 = nc.Dataset('Frequencym20.1.nc')
p018 = np.asarray(ds018['RMP Spike Period'])
t018 = np.asarray(ds018['Trigger Amplitude'])
p018x= p018[1:10000]
t018x= t018[1:10000]

ds020 = nc.Dataset('Frequencym20.2.nc')
p020 = np.asarray(ds020['RMP Spike Period'])
t020 = np.asarray(ds020['Trigger Amplitude'])

ds022 = nc.Dataset('Frequencym20.3.nc')
p022 = np.asarray(ds022['RMP Spike Period'])
t022 = np.asarray(ds022['Trigger Amplitude'])

ds025 = nc.Dataset('Frequencym20.4.nc')
p025 = np.asarray(ds025['RMP Spike Period'])
t025 = np.asarray(ds025['Trigger Amplitude'])

ds027 = nc.Dataset('Frequencym20.5.nc')
p027 = np.asarray(ds027['RMP Spike Period'])
t027 = np.asarray(ds027['Trigger Amplitude'])

ds021 = nc.Dataset('Frequencyp10.1.nc')
p021 = np.asarray(ds021['RMP Spike Period'])
t021 = np.asarray(ds021['Trigger Amplitude'])
p021x= p021[3:10000]
t021x= t021[3:10000]

ds023 = nc.Dataset('Frequencyp10.2.nc')
p023 = np.asarray(ds023['RMP Spike Period'])
t023 = np.asarray(ds023['Trigger Amplitude'])

ds024 = nc.Dataset('Frequencyp10.3.nc')
p024 = np.asarray(ds024['RMP Spike Period'])
t024 = np.asarray(ds024['Trigger Amplitude'])

ds026 = nc.Dataset('Frequencyp10.4.nc')
p026 = np.asarray(ds026['RMP Spike Period'])
t026 = np.asarray(ds026['Trigger Amplitude'])

ds028 = nc.Dataset('Frequencyp10.5.nc')
p028 = np.asarray(ds028['RMP Spike Period'])
t028 = np.asarray(ds028['Trigger Amplitude'])

ds029 = nc.Dataset('Frequencyp20.1.nc')
p029 = np.asarray(ds029['RMP Spike Period'])
t029 = np.asarray(ds029['Trigger Amplitude'])
p029x= p029[3:10000]
t029x= t029[3:10000]

ds030 = nc.Dataset('Frequencyp20.2.nc')
p030 = np.asarray(ds030['RMP Spike Period'])
t030 = np.asarray(ds030['Trigger Amplitude'])

ds031 = nc.Dataset('Frequencyp20.3.nc')
p031 = np.asarray(ds031['RMP Spike Period'])
t031 = np.asarray(ds031['Trigger Amplitude'])

ds032 = nc.Dataset('Frequencyp20.4.nc')
p032 = np.asarray(ds032['RMP Spike Period'])
t032 = np.asarray(ds032['Trigger Amplitude'])

ds033 = nc.Dataset('Frequencyp20.5.nc')
p033 = np.asarray(ds033['RMP Spike Period'])
t033 = np.asarray(ds033['Trigger Amplitude'])

ds034 = nc.Dataset('Frequencyp20.0.nc')
p034 = np.asarray(ds034['RMP Spike Period'])
t034 = np.asarray(ds034['Trigger Amplitude'])
p034x= p034[4:10000]
t034x= t034[4:10000]

fig = plt.figure(figsize=(12, 8))
plt.rc('xtick', labelsize=15) 
plt.rc('ytick', labelsize=15) 
#fig.canvas.manager.set_window_title("NTM Triggering Period Scan")

plt.xlim(0., 20.)
plt.ylim(0., 0.5)
plt.axhline(0.0, color='white', linewidth=0.5, linestyle='dotted')
plt.axvline(0.0, color='white', linewidth=0.5, linestyle='dotted')



plt.plot(p008x, t008x, color='red', linewidth=2, label="-1 kHz")
plt.plot(p009,  t009,  color='red', linewidth=2)
plt.plot(p010,  t010,  color='red', linewidth=2)
plt.plot(p012,  t012,  color='red', linewidth=2)
plt.plot(p015,  t015,  color='red', linewidth=2)
plt.plot(p016,  t016,  color='red', linewidth=2)
plt.plot(p017,  t017,  color='red', linewidth=2)
plt.plot(p019,  t019,  color='red', linewidth=2)

plt.plot(p001x, t001x, color='blue', linewidth=2, label="0 kHz")
plt.plot(p002,  t002,  color='blue', linewidth=2)
plt.plot(p003,  t003,  color='blue', linewidth=2)
plt.plot(p004,  t004,  color='blue', linewidth=2)
plt.plot(p005,  t005,  color='blue', linewidth=2)
plt.plot(p006,  t006,  color='blue', linewidth=2)
plt.plot(p007,  t007,  color='blue', linewidth=2)
plt.plot(p011,  t011,  color='blue', linewidth=2)
plt.plot(p013,  t013,  color='blue', linewidth=2)
plt.plot(p014,  t014,  color='blue', linewidth=2)

plt.plot(p021x, t021x, color='black', linewidth=2, label="+1 kHz")
plt.plot(p023,  t023,  color='black', linewidth=2)
plt.plot(p024,  t024,  color='black', linewidth=2)
plt.plot(p026,  t026,  color='black', linewidth=2)
plt.plot(p028,  t028,  color='black', linewidth=2)

plt.plot(p034x, t034x, color='cyan', linewidth=2, label = "+2 kHz")
plt.plot(p029x, t029x, color='cyan', linewidth=2)
plt.plot(p030,  t030,  color='cyan', linewidth=2)
plt.plot(p031,  t031,  color='cyan', linewidth=2)
plt.plot(p032,  t032,  color='cyan', linewidth=2)
plt.plot(p033,  t033,  color='cyan', linewidth=2)

plt.plot(p018x, t018x, color='green', linewidth=2, label="-2 kHz")
plt.plot(p020,  t020,  color='green', linewidth=2)
plt.plot(p022,  t022,  color='green', linewidth=2)
plt.plot(p025,  t025,  color='green', linewidth=2)
plt.plot(p027,  t027,  color='green', linewidth=2)

plt.xlabel('$Period (ms)$', fontsize="20")
plt.ylabel("$I_{rmp} (kA)$", fontsize="20")
plt.legend(fontsize="20")

plt.show()
#plt.savefig("PeriodScan.pdf")
