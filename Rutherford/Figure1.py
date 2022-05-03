import warnings
warnings.filterwarnings('ignore')

import matplotlib.pyplot as plt

with open('TwoHarmonic.m100.txt' ,"r") as file:
    lines = file.readlines()
    file.close()

Delta_1 = []
epsilon_2  = []
lam = []
lamd  = []
for x in lines:
    data = x.split()
    Delta_1.append(float(data[0]))
    epsilon_2.append(float(data[1]))
    lam.append(float(data[2]))
    lamd.append(float(data[3]))

with open('TwoHarmonic.m200.txt' ,"r") as file:
    lines = file.readlines()
    file.close()

Delta_1a = []
epsilon_2a  = []
lama = []
lamda  = []
for x in lines:
    data = x.split()
    Delta_1a.append(float(data[0]))
    epsilon_2a.append(float(data[1]))
    lama.append(float(data[2]))
    lamda.append(float(data[3]))

plt.figure(figsize=(10,7))
plt.rc('xtick', labelsize=16) 
plt.rc('ytick', labelsize=16) 
plt.xlim(0., 10.)
plt.ylim(0., 0.10)
plt.plot(Delta_1, epsilon_2, color='red',  linewidth='1.5', label="$\Delta_2'=-1.")
plt.plot(Delta_1a, epsilon_2a, color='green',  linewidth='1.5', label="$\Delta_2'=-1.")
plt.xlabel("$\Delta'_1$", fontsize='20')
plt.ylabel("$\epsilon_2$", fontsize='20')
#plt.legend(fontsize='20')

plt.show()    
