import warnings
warnings.filterwarnings('ignore')

import matplotlib.pyplot as plt

with open('pc02.txt' ,"r") as file:
    lines = file.readlines()
    file.close()

mu = []
eps2 = []
lam = []
lam1 = []
for x in lines:
    data = x.split()
    mu.append(float(data[0]))
    eps2.append(float(data[1]))
    lam.append(float(data[2]))
    lam1.append((0.5-float(data[0])*0.63623)/0.8227)

with open('SingleHarmonic.txt' ,"r") as file:
    lines = file.readlines()
    file.close()

plt.figure(figsize=(10,7))
plt.rc('xtick', labelsize=16) 
plt.rc('ytick', labelsize=16) 
plt.xlim(0., 0.3)
#plt.ylim(-0.1, 0.4)
plt.plot(mu, eps2, color='red',  linewidth='1.5', label="$\epsilon_2$")
plt.plot(mu, lam, color='blue',  linewidth='1.5', label="$\lambda$ (2-harmonic)")
plt.plot(mu, lam1, color='green',  linewidth='1.5', linestyle="dashed", label="$\lambda$ (1-harmonic)")
plt.axhline(y=0.25, linestyle='dotted', linewidth='1.5', color='black')
plt.axhline(y=0., linestyle='dotted', linewidth='1.5', color='black')
plt.xlabel("$\mu$", fontsize='20')
#plt.ylabel("$\epsilon_2$", fontsize='20')
plt.legend(fontsize='20')

#plt.show()
plt.savefig("Figure7.eps")
