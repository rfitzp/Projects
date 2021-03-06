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
    lamd.append(float(data[3])-1.2159)

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
    lamda.append(float(data[3])-1.2159)

with open('TwoHarmonic.000.txt' ,"r") as file:
    lines = file.readlines()
    file.close()

Delta_1b = []
epsilon_2b  = []
lamb = []
lamdb  = []
for x in lines:
    data = x.split()
    Delta_1b.append(float(data[0]))
    epsilon_2b.append(float(data[1]))
    lamb.append(float(data[2]))
    lamdb.append(float(data[3])-1.2159)   

with open('TwoHarmonic.m300.txt' ,"r") as file:
    lines = file.readlines()
    file.close()

Delta_1c = []
epsilon_2c  = []
lamc = []
lamdc  = []
for x in lines:
    data = x.split()
    Delta_1c.append(float(data[0]))
    epsilon_2c.append(float(data[1]))
    lamc.append(float(data[2]))
    lamdc.append(float(data[3])-1.2159)   
    
plt.figure(figsize=(10,7))
plt.rc('xtick', labelsize=16) 
plt.rc('ytick', labelsize=16) 
plt.xlim(0., 20.)
plt.ylim(0., 0.1)
plt.plot(Delta_1, lamd, color='red',  linewidth='1.5', label="$\Delta_2'=-1.")
plt.plot(Delta_1a, lamda, color='green',  linewidth='1.5', label="$\Delta_2'=-1.")
plt.plot(Delta_1c, lamdc, color='blue',  linewidth='1.5', label="$\Delta_2'=-1.")
plt.plot(Delta_1b, lamdb, color='black',  linewidth='1.5', label="$\Delta_2'=-1.")
plt.xlabel("$\Delta'_1$", fontsize='20')
plt.ylabel("$\lambda/\Delta_1'-I_1^{-1}$", fontsize='20')
#plt.legend(fontsize='20')

plt.show()    
