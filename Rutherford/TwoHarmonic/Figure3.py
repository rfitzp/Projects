import warnings
warnings.filterwarnings('ignore')

import matplotlib.pyplot as plt
   
with open('Data.txt' ,"r") as file:
    lines = file.readlines()
    file.close()

pc = []
mu = []
eps2 = []
for x in lines:
    data = x.split()
    pc.append(float(data[2]))
    mu.append(float(data[3]))
    eps2.append(float(data[4]))

with open('SingleHarmonic.txt' ,"r") as file:
    lines = file.readlines()
    file.close()

ppc = []
kk1 = []
mmu = []
for x in lines:
    data = x.split()
    ppc.append(float(data[0]))
    kk1.append(float(data[1]))
    mmu.append(0.5/float(data[1]))
        
plt.figure(figsize=(10,7))
plt.rc('xtick', labelsize=16) 
plt.rc('ytick', labelsize=16) 
plt.xlim(0., 2.)
plt.ylim(-0.1, 0.4)
plt.plot(pc, eps2, color='red',  linewidth='1.5', label="$\epsilon_2$")
plt.plot(pc, mu, color='blue',  linewidth='1.5', label="$\mu_{crit}$ (2-harmonic)")
plt.plot(ppc, mmu, color='green',  linewidth='1.5', linestyle="dashed", label="$\mu_{crit}$ (1-harmonic)")
plt.axhline(y=0.25, linestyle='dotted', linewidth='1.5', color='black')
plt.axhline(y=0., linestyle='dotted', linewidth='1.5', color='black')
plt.xlabel("$p_c$", fontsize='20')
#plt.ylabel("$\epsilon_2$", fontsize='20')
plt.legend(fontsize='20')

#plt.show()
plt.savefig("Figure7.eps")
