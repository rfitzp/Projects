# ###################################################################
# Program to implement two-harmonic Rutherford island width evolution
# ###################################################################

import warnings
warnings.filterwarnings('ignore')

import scipy.integrate as integrate
import scipy.special as special
import numpy as np
import numpy.linalg as lin
import math
import matplotlib.pyplot as plt

# ###########
# Get Delta_2
# ###########
print ("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nTwo-Harmonic Rutherford Island Calculation\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
val1, val2, val3 = input("Delta_2 epsilon_2 verbose ?? ").split()
Delta_2 = float(val1)
Eps2_old = float(val2)
verb = int(val3)

# ######################
# Calculation parameters
# ######################
d1_start = 1.e-3   # Delta_1 scan start value
d1_mid   = 5.      # Delta_1 scan middle value
d1_end   = 20.     # Delta_1 scan end value
Nd1      = 51      # Number of points Delta_1 scan

eps      = 1.e-4   # phi integrand regularization parameter
eta      = 1.e-5   # k integrand regularization parameter
kmax     = 20.     # Maximum k value 
tol      = 1.e-10  # Root finding accuracy
step     = 1.e-3   # Initial step in root finding

# ################
# Define functions
# ################

def Cosxi02(k, eps2):
    if k >= 1.:
        return 1.
    elif k >= 0.:
        if eps2 == 0.:
            return k
        else:
            return math.sqrt((4.*eps2-1. + math.sqrt((4.*eps2-1.)*(4.*eps2-1.) + 16.*eps2*k*k)) /8./eps2)
    else:
        return math.sqrt((4.*eps2-1. + math.sqrt((4.*eps2-1.)*(4.*eps2-1.) - 16.*eps2*k*k)) /8./eps2)

def Cosxi12(k, eps2):
    if k >= 0.:
        return 0.
    else:
        return math.sqrt((4.*eps2-1. - math.sqrt((4.*eps2-1.)*(4.*eps2-1.) - 16.*eps2*k*k)) /8./eps2)

def Get_kmin(eps2):
    if eps2 <= 0.25:
        return 0.
    else:
        return - (4.*eps2-1.) /4./math.sqrt(eps2)

def Cos_m_Integrand_1(y, k, m, eps2):
    siny = math.sin(y)
    fun = k*k - (1.-4.*eps2)*siny*siny - 4.*eps2*siny*siny*siny*siny
    if fun < 0.:
        return 0.
    else:
        return math.cos(2.*m*(math.pi/2.-y)) /math.sqrt(fun) /math.pi

def Cos_m_Integrand_2(y, k, m, eps2):
    kk = Cosxi02(k, eps2)
    siny = math.sin(y)
    sink = kk*siny
    fun = k*k - (1.-4.*eps2)*sink*sink - 4.*eps2*sink*sink*sink*sink
    if fun < 0. or sink > 1.:
        return 0.
    else:
        return math.cos(2.*m*math.acos(sink)) *kk*math.cos(y) /math.sqrt(fun) \
            /math.sqrt(1.-sink*sink) /math.pi

def Cos_m_Integrand_3(y, k, m, eps2):
    kk = Cosxi02(k, eps2)
    siny = math.sin(y)
    sink = kk*siny
    fun = - k*k - (1.-4.*eps2)*sink*sink - 4.*eps2*sink*sink*sink*sink
    if fun < 0. or sink > 1.:
        return 0.
    else:
        return math.cos(2.*m*math.acos(sink)) *kk*math.cos(y) /math.sqrt(fun) \
            /math.sqrt(1.-sink*sink) /math.pi    

def Cos_m(k, m, eps2):
    if k >= 1.:
        return integrate.quad(Cos_m_Integrand_1, 0., math.pi/2., args=(k, m, eps2))[0]
    elif k >= 0.:
        return integrate.quad(Cos_m_Integrand_2, eps, math.pi/2., args=(k, m, eps2))[0]
    else:
        y0 = math.asin(Cosxi12(k, eps2) /Cosxi02(k, eps2))
        return integrate.quad(Cos_m_Integrand_3, y0+eps, math.pi/2., args=(k, m, eps2))[0]

def I_Element_Integrand(k, m, mp, eps2):
    return 8. * abs(k) * Cos_m(k, m, eps2) * Cos_m(k, mp, eps2) /Cos_m(k, 0, eps2)

def I_Element(m, mp, eps2):
    kmin = Get_kmin(eps2)
    if kmin < 0.:
        return integrate.quad(I_Element_Integrand, kmin, -eta, args=(m, mp, eps2))[0] \
            + integrate.quad(I_Element_Integrand, eps, 1.-eta, args=(m, mp, eps2))[0] \
            + integrate.quad(I_Element_Integrand, 1.+eta, kmax, args=(m, mp, eps2))[0]
    else:
        return integrate.quad(I_Element_Integrand, eta, 1.-eta, args=(m, mp, eps2))[0] \
            + integrate.quad(I_Element_Integrand, 1.+eta, kmax, args=(m, mp, eps2))[0]

def I_Matrix(eps2):
    for m in range(2):
        for mp in range(m+1):
            if mp == m:
                I[m,mp] = I_Element(M[m], M[mp], eps2)
            else:
                I[m,mp] = I_Element(M[m], M[mp], eps2)
                I[mp,m] = I[m,mp]
             
def lambda_plus(eps2):
    I_Matrix(eps2)
    
    I11 = I[0,0]
    I12 = I[0,1]
    I22 = I[1,1]
    D1 = Delta[0]
    D2 = Delta[1]
    fac1 = I11*D2 + I22*D1
    fac2 = I11*D2 - I22*D1
    fac3 = 4.*D1*D2*I12*I12

    return (fac1 + math.sqrt(fac2*fac2 + fac3)) /2./(I11*I22 - I12*I12)

def f_plus(eps2):
    Lambda_plus = lambda_plus(eps2)

    return I[0,0] + I[0,1] * eps2 - Delta[0] /Lambda_plus, Lambda_plus

def FindRoot(x1):
    x2 = x1 + step
    f1, l1 = f_plus(x1)
    if verb:
        print ("x = %11.4e f = %11.4e" % (x1, f1))
    f2, l2 = f_plus(x2)
    if verb:
        print ("x = %11.4e f = %11.4e" % (x2, f2))

    feval = 2
    while True and feval < 15:
        x = (f1*x2 - f2*x1) /(f1 - f2)
        f, l = f_plus(x)
        feval += 1   
        if verb:
            print ("x = %11.4e f = %11.4e" % (x, f))

        if abs(f) < tol:
            break
        
        if f*f1 < 0.:
            x2 = x
            f2 = f
        else:
            x1 = x
            f1 = f

    return x, f, l, feval       
                   
# #################
# Construct m array
# #################
I = np.zeros((2,2))
M = np.zeros(2)
for m in range(2):
    M[m] = (m+1)*1.
   
# ############
# Scan Delta_1
# ############
ddd = []
ep2 = []
lpl = []
lpd = []
for Delta_1 in np.linspace(d1_start, d1_mid, Nd1):

    # ######################
    # Construct Delta arrays
    # ######################
    Delta = np.zeros(2)
    Delta[0] = Delta_1
    Delta[1] = Delta_2

    # #########    
    # Find root
    # #########
    Eps2, F_plus, Lambda_plus, feval = FindRoot(Eps2_old)
    Eps2_old = Eps2

    ddd.append(Delta[0])
    ep2.append(Eps2)
    lpl.append(Lambda_plus)
    lpd.append(Lambda_plus/Delta[0])
    print ("Delta_1 = %11.4e  F_+ = %11.4e  feval = %2d  epsilon_2 = %11.4e  lambda = %11.4e  lambda/Delta_1 = %11.4e" \
           % (Delta[0], F_plus, feval, Eps2, Lambda_plus, Lambda_plus/Delta[0]))

for Delta_1 in np.linspace(d1_mid, d1_end, Nd1):

    # ######################
    # Construct Delta arrays
    # ######################
    Delta = np.zeros(2)
    Delta[0] = Delta_1
    Delta[1] = Delta_2

    # #########    
    # Find root
    # #########
    Eps2, F_plus, Lambda_plus, feval = FindRoot(Eps2_old)
    Eps2_old = Eps2

    ddd.append(Delta[0])
    ep2.append(Eps2)
    lpl.append(Lambda_plus)
    lpd.append(Lambda_plus/Delta[0])
    print ("Delta_1 = %11.4e  F_+ = %11.4e  feval = %2d  epsilon_2 = %11.4e  lambda = %11.4e  lambda/Delta_1 = %11.4e" \
           % (Delta[0], F_plus, feval, Eps2, Lambda_plus, Lambda_plus/Delta[0]))
    
# ###########
# Output data
# ###########
with open("TwoHarmonic.txt", "w") as file:
    for n in range(len(ddd)):
        file.write("%11.4e %11.4e %11.4e %11.4e\n" % (ddd[n], ep2[n], lpl[n], lpd[n]))
    
# ##########    
# Graph data
# ##########
plt.figure(figsize=(10,7))
plt.rc('xtick', labelsize=16) 
plt.rc('ytick', labelsize=16) 
plt.xlim(0.,10.)
plt.ylim(0., 1.3)
plt.plot(ddd, ep2, color='blue',  linewidth='1.5', label='$\epsilon_2$')
plt.plot(ddd, lpd, color='red',  linewidth='1.5', label="$\lambda/\Delta_1$'")
plt.legend(fontsize='20')
plt.xlabel("$\Delta_1'$", fontsize='20')

plt.show()
#plt.savefig("TwoHarmonic.eps")    
