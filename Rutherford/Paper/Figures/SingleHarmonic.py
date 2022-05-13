# ######################################################################
# Program to implement single-harmonic Rutherford island width evolution
# ######################################################################

import warnings
warnings.filterwarnings('ignore')

import scipy.integrate as integrate
import scipy.special as special
import numpy as np
import numpy.linalg as lin
import math
import matplotlib.pyplot as plt

print ("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nSingle-Harmonic Rutherford Island Calculation\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

# ######################
# Calculation parameters
# ######################
pc_start = 1.e-3   # p_c scan start value
pc_end   = 4.      # p_c scan end value
Npc      = 501     # Number of points Delta_1 scan

eta      = 1.e-7   # k integrand regularization parameter
pmax     = 20.     # Maximum p value 

# ################
# Define functions
# ################

def Cos_m_Integrand_Outer(y, p, m):
    siny = math.sin(y)
    fun = p*p - siny*siny
    if fun < 0.:
        return 0.
    else:
        return math.cos(2.*m *math.acos(siny)) /math.sqrt(fun) /math.pi

def Cos_m_Integrand_Inner(y, p, m):
    sinp = p * math.sin(y)
    fun = 1. - sinp*sinp
    if fun < 0.:
        return 0
    else:
        return math.cos(2.*m *math.acos(sinp)) /math.sqrt(fun) /math.pi

def Cos_m(p, m):
    if p >= 1.:
        return integrate.quad(Cos_m_Integrand_Outer, 0., math.pi/2., args=(p, m))[0]
    else :
        return integrate.quad(Cos_m_Integrand_Inner, 0., math.pi/2., args=(p, m))[0]

def I_Integrand(p, m, mp):
    return 8. * p * Cos_m(p, m) * Cos_m(p, mp) /Cos_m(p, 0)

def I(m, mp):
    return integrate.quad(I_Integrand, 0., 1.-eta, args=(m, mp))[0] \
    + integrate.quad(I_Integrand, 1.+eta, pmax, args=(m, mp))[0]

def L_Integrand(p, m):
    return -64. * p * Cos_m(p, m)

def L(m, pc):
    if pc > 1.:
        return integrate.quad(L_Integrand, 0., 1.-eta, args=(m))[0] \
        + integrate.quad(L_Integrand, 1.+eta, pc, args=(m))[0]
    else:
        return integrate.quad(L_Integrand, 0., pc, args=(m))[0]

# #################
# Calculate I_{1,1}
# #################                        
I_11 = I(1, 1)
print ("I_{1,1} = %11.4e" %I_11)

# #######
# Scan pc
# #######
pcc = []
ll1 = []
for pc in np.linspace(pc_start, pc_end, Npc):
    L1 =  L(1, pc)
    print ("p_c = %11.4e  L_1 = %11.4e" %(pc, L1))
    pcc.append(pc)
    ll1.append(L1)

# ###########
# Output data
# ###########
with open("SingleHarmonic.txt", "w") as file:
    for n in range(len(pcc)):
        file.write("%11.4e %11.4e\n" % (pcc[n], ll1[n]))

# ##########    
# Graph data
# ##########
plt.figure(figsize=(10,7))
plt.rc('xtick', labelsize=16) 
plt.rc('ytick', labelsize=16) 
plt.xlim(0.,pc_end)
plt.ylim(0., 8.)
plt.plot(pcc, ll1, color='black',  linewidth='1.5')
plt.xlabel("$p_c$", fontsize='20')
plt.ylabel("$K_1$", fontsize='20')

#plt.show()
plt.savefig("SingleHarmonic.eps")    
    
