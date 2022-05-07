# Script to plot magnetic flux surfaces in single-harmonic approximation

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import math
import matplotlib.pyplot as plt

val = input("epsilon_2 ? ")
eps2 = float(val)

def Psi(x, y):
    xi = y*math.pi
    return x*x/2. + math.cos(xi) + eps2*math.cos(2.*xi)

xmax = 4.
nxpoint = 1001
nypoint = 1001
X = np.linspace(-xmax, xmax, nxpoint)
Y = np.linspace(0., 6., nypoint)
Z = np.zeros((len(X), len(Y)))
for i in range(len(X)):
    for j in range(len(Y)):
        Z[i,j] = Psi(X[i], Y[j])
            
plt.figure(figsize=(12,8))
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.contour(Y, X, Z, levels=40, linewidths=1.5, colors='black')
if eps2 > -0.25:
    plt.contour(Y, X, Z, levels=[1.+eps2], linewidths=3, colors='black')
if eps2 > 0.25:
    plt.contour(Y, X, Z, levels=[-1.+eps2], linewidths=3, colors='black', linestyle='solid')
    plt.contour(Y, X, Z, levels=[0.25*(-1.+eps2)+0.75*(-eps2-1./8./eps2)], linewidths=1.5, colors='black', linestyle='solid')
if eps2 < -0.25:
    plt.contour(Y, X, Z, levels=[-eps2-1./8./eps2], linewidths=3, colors='black', linestyle='solid')
    plt.contour(Y, X, Z, levels=[0.75*(1.+eps2)+0.25*(-eps2-1./8./eps2)], linewidths=1.5, colors='black', linestyle='solid')
plt.ylabel("$X$", fontsize='20')
plt.xlabel(r"$\xi / \pi$", fontsize='20')

#plt.show()
plt.savefig("Flux.eps")        
