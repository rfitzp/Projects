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

# ######################
# Calculation parameters
# ######################
eps = 1.e-4       # Integrand regularization parameter
kmax = 20.        # Maximum k value in k integrals
tol = 1.e-10      # Root finding accuracy
step = 1.e-4      # Initial step in root finding

k1_start = 0.20    # Lowest wavenumber scan start value
k1_end   = 0.15   # Lowest wavenumber scan end value
Nk1      = 101    # Number of points in lowest wavenumber scan

# ################
# Define functions
# ################

def Delta_Element(k):
    return 2.*(1./k - k)

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
    return math.cos(2.*m*(math.pi/2.-y)) /math.sqrt(k*k-(1.-4.*eps2)*siny*siny-4.*eps2*siny*siny*siny*siny) /math.pi

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
        return integrate.quad(I_Element_Integrand, kmin, -eps, args=(m, mp, eps2))[0] \
            + integrate.quad(I_Element_Integrand, eps, 1.-eps, args=(m, mp, eps2))[0] \
            + integrate.quad(I_Element_Integrand, 1.+eps, kmax, args=(m, mp, eps2))[0]
    else:
        return integrate.quad(I_Element_Integrand, eps, 1.-eps, args=(m, mp, eps2))[0] \
            + integrate.quad(I_Element_Integrand, 1.+eps, kmax, args=(m, mp, eps2))[0]

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
    fac3 = 4.*I11*I22*I12*I12

    return (fac1 + math.sqrt(fac2*fac2 + fac3)) /2./I11/I22

def f_plus(eps2):
    Lambda_plus = lambda_plus(eps2)
    return Lambda_plus * I[0,0] - Delta[0] + I[0,1]*eps2, Lambda_plus

def I12(x):
    return I_Element(1, 2, x)

def FindRoot0(x1):
    x2 = x1 + step
    f1 = I12(x1)
    f2 = I12(x2)

    while True:
        x = (f1*x2 - f2*x1) /(f1 - f2)
        f = I12(x)
        #print ("x = %11.4e f = %11.4e" % (x, f))

        if abs(f) < tol:
            break
        
        if f*f1 < 0.:
            x2 = x
            f2 = f
        else:
            x1 = x
            f1 = f

    return x, f   

def FindRoot(x1):
    x2 = x1 + step
    f1, l1 = f_plus(x1)
    f2, l2 = f_plus(x2)

    while True:
        x = (f1*x2 - f2*x1) /(f1 - f2)
        f, l = f_plus(x)
        #print ("x = %11.4e f = %11.4e" % (x, f))

        if abs(f) < tol:
            break
        
        if f*f1 < 0.:
            x2 = x
            f2 = f
        else:
            x1 = x
            f1 = f

    return x, f, l       
                   
print ("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nTwo-Harmonic Rutherford Island Calculation\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

# #################
# Construct m array
# #################
M = np.zeros(2)
for m in range(2):
    M[m] = (m+1)*1.

# #######################
# Calculate critical Eps2
# #######################
I = np.zeros((2,2))
Eps2_crit = 0.157
Eps2_crit, I12_crit = FindRoot0(Eps2_crit)
I11_crit = I_Element(M[0], M[0], Eps2_crit)
I22_crit = I_Element(M[1], M[1], Eps2_crit)

print ("Eps2_crit = %11.4e  I_11_crit = %11.4e  I_12_crit = %11.4e  I_22_crit = %11.5e" % (Eps2_crit, I11_crit, I12_crit, I22_crit))

# #######
# Scan k1
# #######
kk1 = []
ddd = []
ep2 = []
fpl = []
i12 = []
lpl = []
lc1 = []
lc2 = []
Eps2_old = 9.3e-2
for k1 in np.linspace(k1_start, k1_end, Nk1):

    # ############################
    # Construct k and Delta arrays
    # ############################
    K = np.zeros(2)
    Delta = np.zeros(2)
    for m in range(2):
        K[m] = (m+1)*k1
        Delta[m] = Delta_Element(K[m])

    # #########    
    # Find root
    # #########
    Eps2, F_plus, Lambda_plus = FindRoot(Eps2_old)
    Eps2_old = Eps2

    kk1.append(k1)
    ddd.append(Delta[0])
    ep2.append(Eps2)
    fpl.append(F_plus)
    i12.append(-I[0,1])
    lpl.append(Lambda_plus)
    lc1.append(Delta[0]/I11_crit)
    lc2.append(Delta[1]/I22_crit)
    print ("k_1 = %11.4e Del_1 = %11.4e eps_2 = %11.4e F_+ = %11.4e I_12 = %11.4e lam = %11.4e lam_1 = %11.4e lam_2 = %11.4e" \
           % (k1, Delta[0], Eps2, F_plus, I[0,1], Lambda_plus, Delta[0]/I11_crit, Delta[1]/I22_crit))

# ###########
# Output data
# ###########
with open("TwoHarmonic_2.txt", "w") as file:
    for n in range(len(kk1)):
        file.write("%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n" % (kk1[n], ddd[n], ep2[n], fpl[n], i12[n], lpl[n], lc1[n], lc2[n]))
    
# ##########    
# Graph data
# ##########
plt.figure(figsize=(10,7))
plt.rc('xtick', labelsize=16) 
plt.rc('ytick', labelsize=16) 
plt.xlim(k1_end,k1_start)
plt.ylim(0., 0.5)
lpl1 = np.asarray(lpl)/50.
lc11 = np.asarray(lc1)/50.
lc21 = np.asarray(lc2)/50.
plt.plot(kk1, ep2, color='blue',  linewidth='1.5', label='$\epsilon_2$')
plt.plot(kk1, lpl1, color='red',  linewidth='1.5', label='$\lambda/50$')
plt.plot(kk1, lc11, color='yellow',  linewidth='1.5', label='$\lambda_1$/50')
plt.plot(kk1, lc21, color='magenta',  linewidth='1.5', label='$\lambda_2$/50')
plt.plot(kk1, i12, color='green',  linewidth='1.5', label='-$A_{12}$')
plt.axhline(y=Eps2_crit, linewidth='1.5', linestyle='dashed', color='blue')
plt.legend(fontsize='20')
plt.xlabel('$k_1$', fontsize='20')

plt.show()
#plt.savefig("TwoHarmonic.eps")    
