import warnings
warnings.filterwarnings('ignore')

import matplotlib.pyplot as plt

with open('Intrinsic.000.txt' ,"r") as file:
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

with open('Intrinsic.m050.txt' ,"r") as file:
    lines = file.readlines()
    file.close()

Delta_1e = []
epsilon_2e  = []
lame = []
lamde  = []
for x in lines:
    data = x.split()
    Delta_1e.append(float(data[0]))
    epsilon_2e.append(float(data[1]))
    lame.append(float(data[2]))
    lamde.append(float(data[3])-1.2159) 

with open('Intrinsic.m100.txt' ,"r") as file:
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

with open('Intrinsic.m200.txt' ,"r") as file:
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
    
with open('Intrinsic.m400.txt' ,"r") as file:
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
plt.plot(Delta_1b, lamdb, color='black',  linewidth='1.5')
plt.plot(Delta_1e, lamde, color='red',  linewidth='1.5')
plt.plot(Delta_1, lamd, color='green',  linewidth='1.5')
plt.plot(Delta_1a, lamda, color='blue',  linewidth='1.5')
plt.plot(Delta_1c, lamdc, color='magenta',  linewidth='1.5')
plt.xlabel("$\Delta'_1$", fontsize='20')
plt.ylabel("$(\lambda/\Delta_1')_{2-harmonic}-(\lambda/\Delta_1')_{1-harmonic}$", fontsize='20')
#plt.legend(fontsize='20')

#plt.show()    
plt.savefig("Figure6.eps")
