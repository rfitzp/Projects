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
    epsilon_2.append(float(data[5]))

plt.figure(figsize=(10,7))
plt.rc('xtick', labelsize=16) 
plt.rc('ytick', labelsize=16) 
plt.xlim(0., 14.)
plt.ylim(0., 20.)
plt.plot(Delta_1, epsilon_2, color='black',  linewidth='1.5')
plt.xlabel("$\Delta'_1$", fontsize='20')
plt.ylabel("$\lambda$", fontsize='20')

plt.show()    
