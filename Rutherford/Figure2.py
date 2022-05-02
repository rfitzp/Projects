import warnings
warnings.filterwarnings('ignore')

import matplotlib.pyplot as plt

with open('TwoHarmonic_2.txt' ,"r") as file:
    lines = file.readlines()
    file.close()

Delta_1 = []
epsilon_2  = []
for x in lines:
    data = x.split()
    Delta_1.append(float(data[1]))
    epsilon_2.append(float(data[2]))

plt.figure(figsize=(10,7))
plt.rc('xtick', labelsize=16) 
plt.rc('ytick', labelsize=16) 
plt.xlim(0., 14.)
plt.ylim(0., 0.16)
plt.plot(Delta_1, epsilon_2, color='black',  linewidth='1.5')
plt.axhline(y=0.15769, linewidth='1.5', linestyle='dashed', color='black')
plt.xlabel("$\Delta'_1$", fontsize='20')
plt.ylabel("$\epsilon_2$", fontsize='20')

plt.show()    
