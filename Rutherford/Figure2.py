import warnings
warnings.filterwarnings('ignore')

import matplotlib.pyplot as plt

with open('TwoHarmonic_1.txt' ,"r") as file:
    lines = file.readlines()
    file.close()

Delta_1 = []
epsilon_2  = []
for x in lines:
    data = x.split()
    Delta_1.append(float(data[5]))
    epsilon_2.append(float(data[2]))

with open('TwoHarmonic_2.txt' ,"r") as file:
    lines = file.readlines()
    file.close()

Delta_1a = []
epsilon_2a  = []
for x in lines:
    data = x.split()
    Delta_1a.append(float(data[5]))
    epsilon_2a.append(float(data[2]))

with open('TwoHarmonic_3.txt' ,"r") as file:
    lines = file.readlines()
    file.close()

Delta_1b = []
epsilon_2b  = []
for x in lines:
    data = x.split()
    Delta_1b.append(float(data[5]))
    epsilon_2b.append(float(data[2]))        

with open('TwoHarmonic_4.txt' ,"r") as file:
    lines = file.readlines()
    file.close()

Delta_1c = []
epsilon_2c  = []
for x in lines:
    data = x.split()
    Delta_1c.append(float(data[5]))
    epsilon_2c.append(float(data[2]))        

with open('TwoHarmonic_5.txt' ,"r") as file:
    lines = file.readlines()
    file.close()

Delta_1d = []
epsilon_2d  = []
for x in lines:
    data = x.split()
    Delta_1d.append(float(data[5]))
    epsilon_2d.append(float(data[2]))

Delta_1d.append(17.793)
epsilon_2d.append(0.15731)
    
plt.figure(figsize=(10,7))
plt.rc('xtick', labelsize=16) 
plt.rc('ytick', labelsize=16) 
plt.xlim(0., 20.)
plt.ylim(0., 0.32)
plt.plot(Delta_1, epsilon_2, color='black',  linewidth='1.5')
plt.plot(Delta_1a, epsilon_2a, color='black',  linewidth='1.5')
plt.plot(Delta_1b, epsilon_2b, color='black',  linewidth='1.5')
plt.plot(Delta_1c, epsilon_2c, color='black',  linewidth='1.5')
plt.plot(Delta_1d, epsilon_2d, color='black',  linewidth='1.5')
plt.axhline(y=0.15769, linewidth='1.5', linestyle='dashed', color='black')
plt.legend(fontsize='20')
plt.xlabel("$\lambda$", fontsize='20')
plt.ylabel("$\epsilon_2$", fontsize='20')

plt.show()    
