# ##################################################################
# Script to implement two-harmonic Rutherford island width evolution
# ##################################################################

import warnings
warnings.filterwarnings('ignore')

import scipy.integrate as integrate
import scipy.special as special
import numpy as np
import numpy.linalg as lin
import math
import matplotlib.pyplot as plt

# #####
# p_min
# #####
def Get_pmin(eps2):

    if eps2 <= 0.25:
        return 0.
    else:
        return - (4.*eps2-1.) /4./math.sqrt(eps2)

# #####
# p_sep
# #####
def Get_psep(eps2):

    if eps2 >= -0.25:
        return 1.
    else:
        return (1. - 4.*eps2) /4./math.sqrt(abs(eps2))    

# ###########
# cos(xi_0/2)
# ###########
def Cosxi02(p, eps2):

    if p >= 0.:
        arg = (4.*eps2-1. + math.sqrt((4.*eps2-1.)*(4.*eps2-1.) + 16.*eps2*p*p)) /8./eps2
    else:
        arg = (4.*eps2-1. + math.sqrt((4.*eps2-1.)*(4.*eps2-1.) - 16.*eps2*p*p)) /8./eps2

    if arg >= 1.:
        return 1.
    elif arg <= 0.:
        return 0.
    else:
        return math.sqrt(arg)

# ###########
# cos(xi_1/2)
# ###########
def Cosxi12(p, eps2):

    if p >= 0.:
        arg = (4.*eps2-1. - math.sqrt((4.*eps2-1.)*(4.*eps2-1.) + 16.*eps2*p*p)) /8./eps2
    else:
        arg = (4.*eps2-1. - math.sqrt((4.*eps2-1.)*(4.*eps2-1.) - 16.*eps2*p*p)) /8./eps2

    if arg >= 1.:
        return 1.
    elif arg <= 0.:
        return 0.
    else:
        return math.sqrt(arg)

# ################################
# Type-1 integrand for <cos(m xi)>
# ################################
def Cos_m_Integrand_1(y, p, m, eps2):

    siny = math.sin(y)
    fun = p*p - (1.-4.*eps2)*siny*siny - 4.*eps2*siny*siny*siny*siny

    if fun <= 0.:
        return 0.
    else:
        return math.cos(2.*m * math.acos(siny)) /math.sqrt(fun) /math.pi

# ################################
# Type-2 integrand for <cos(m xi)>
# ################################
def Cos_m_Integrand_2(y, p, m, eps2):

    pp = Cosxi02(p, eps2)
    siny = math.sin(y)
    sinp = pp*siny

    fun = p*p - (1.-4.*eps2)*sinp*sinp - 4.*eps2*sinp*sinp*sinp*sinp
    if fun <= 0. or sinp >= 1.:
        return 0.
    else:
        return math.cos(2.*m * math.acos(sinp)) * pp*math.cos(y) /math.sqrt(fun) \
            /math.sqrt(1.-sinp*sinp) /math.pi

# ################################
# Type-2 integrand for <cos(m xi)>
# ################################
def Cos_m_Integrand_3(y, p, m, eps2):
    
    pp = Cosxi02(p, eps2)
    siny = math.sin(y)
    sinp = pp*siny
    
    fun = - p*p - (1.-4.*eps2)*sinp*sinp - 4.*eps2*sinp*sinp*sinp*sinp
    if fun < 0. or sinp > 1.:
        return 0.
    else:
        return math.cos(2.*m * math.acos(sinp)) * pp*math.cos(y) /math.sqrt(fun) \
            /math.sqrt(1.-sinp*sinp) /math.pi    

# ###########
# <cos(m xi)>
# ###########
def Cos_m(p, m, eps2):

    psep = Get_psep(p)
    
    if p >= psep:
        return integrate.quad(Cos_m_Integrand_1, 0., math.pi/2., args=(p, m, eps2))[0]
    elif eps2 > 0.25:
        if p >= 0.:
            return integrate.quad(Cos_m_Integrand_2, 0., math.pi/2., args=(p, m, eps2))[0]
        else:
            y0 = math.asin(Cosxi12(p, eps2) /Cosxi02(p, eps2))
            return integrate.quad(Cos_m_Integrand_3, y0, math.pi/2., args=(p, m, eps2))[0]
    elif eps2 < -0.25:
        if p > 1.:
            y1 = math.asin(Cosxi12(p, eps2))
            return integrate.quad(Cos_m_Integrand_2, 0., math.pi/2., args=(p, m, eps2))[0] \
                + integrate.quad(Cos_m_Integrand_1, y1, math.pi/2., args=(p, m, eps2))[0]
        else:
            return integrate.quad(Cos_m_Integrand_2, 0., math.pi/2., args=(p, m, eps2))[0]
    else:
        return integrate.quad(Cos_m_Integrand_2, 0., math.pi/2., args=(p, m, eps2))[0]

# ######################
# Integrand for I_{m,m'}
# #####################
def I_mmp_Integrand(p, m, mp, eps2):

    if m == mp:
        Cosm = Cos_m(p, m, eps2)
        return 8. * abs(p) * Cosm * Cosm /Cos_m(p, 0, eps2)
    else:
        return 8. * abs(p) * Cos_m(p, m, eps2) * Cos_m(p, mp, eps2) /Cos_m(p, 0, eps2)

# ########
# I_{m,m'}
# ########
def I_mmp(m, mp, eps2, eta, pmax):

    pmin = Get_pmin(eps2)
    psep = Get_psep(eps2)

    if pmin < 0.:
        return integrate.quad(I_mmp_Integrand, pmin+eta, -eta, args=(m, mp, eps2))[0] \
            + integrate.quad(I_mmp_Integrand, eta, psep-eta, args=(m, mp, eps2))[0] \
            + integrate.quad(I_mmp_Integrand, psep+eta, pmax, args=(m, mp, eps2))[0]
    elif psep > 1.:
        return integrate.quad(I_mmp_Integrand, pmin+eta, 1.-eta, args=(m, mp, eps2))[0] \
            + integrate.quad(I_mmp_Integrand, 1.+eta, psep-eta, args=(m, mp, eps2))[0] \
            + integrate.quad(I_mmp_Integrand, psep+eta, pmax, args=(m, mp, eps2))[0]
    else:
        return integrate.quad(I_mmp_Integrand, pmin+eta, psep-eta, args=(m, mp, eps2))[0] \
            + integrate.quad(I_mmp_Integrand, psep+eta, pmax, args=(m, mp, eps2))[0]

# #################
# Integrand for K_m
# #################
def K_m_Integrand(p, m, eps2):

    return -64. * abs(p) * Cos_m(p, m, eps2)

# ###
# K_m
# ###
def K_m(eps2, m, pc, eta):

    pmin = Get_pmin(eps2)
    psep = Get_psep(eps2)

    if pmin < 0.:
        if pc < 0.:
            return integrate.quad(K_m_Integrand, pmin+eta, pc, args=(m, eps2))[0]
        elif pc < psep:
            return integrate.quad(K_m_Integrand, pmin+eta, -eta, args=(m, eps2))[0] \
                + integrate.quad(K_m_Integrand, eta, pc, args=(m, eps2))[0]
        else:
            return integrate.quad(K_m_Integrand, pmin+eta, -eta, args=(m, eps2))[0] \
                + integrate.quad(K_m_Integrand, eta, psep-eta, args=(m, eps2))[0] \
                + integrate.quad(K_m_Integrand, psep+eta, pc, args=(m, eps2))[0]
    elif psep > 1.:
        if pc < 1.:
            return integrate.quad(K_m_Integrand, pmin+eta, pc, args=(m, eps2))[0]
        elif pc < psep:
            return integrate.quad(K_m_Integrand, pmin+eta, 1.-eta, args=(m, eps2))[0] \
                + integrate.quad(K_m_Integrand, 1.+eta, pc, args=(m, eps2))[0]
        else:
            return integrate.quad(K_m_Integrand, pmin+eta, 1.-eta, args=(m, eps2))[0] \
                + integrate.quad(K_m_Integrand, 1.+eta, psep-eta, args=(m, eps2))[0] \
                + integrate.quad(K_m_Integrand, psep+eta, pc, args=(m, eps2))[0]
    else:
        if pc < psep:
            return integrate.quad(K_m_Integrand, pmin+eta, pc, args=(m, eps2))[0]
        else:
           return integrate.quad(K_m_Integrand, pmin+eta, psep-eta, args=(m, eps2))[0] \
                + integrate.quad(K_m_Integrand, psep+eta, pc, args=(m, eps2))[0] 

# ########
# lambda_+
# ########
def lambda_plus(eps2, Delta1, Delta2, mu, pc, eta, pmax):

    I11 = I_mmp(1, 1, eps2, eta, pmax)
    I12 = I_mmp(1, 2, eps2, eta, pmax)
    I22 = I_mmp(2, 2, eps2, eta, pmax)
    D1 = Delta1
    D2 = Delta2

    if mu > 0.:
        K1 = K_m(eps2, 1, pc, eta)
        K2 = K_m(eps2, 2, pc, eta)

        A = I11*I22 - I12*I12
        B =  - I11*D2 - I22*D1 + mu*I22*K1 - mu*I12*K2
        C = D2*(D1 - mu*K1)
    else:
        K1 = 0.
        K2 = 0.
        
        A = I11*I22 - I12*I12
        B = - I11*D2 - I22*D1
        C = D2*D1

    return (-B + math.sqrt(B*B - 4.*A*C)) /2./A, I11, I12, I22, K1, K2

# ######
# F_plus
# ######
def F_plus(eps2, Delta1, Delta2, mu, pc, eta, pmax):

    Lambda_plus, I11, I12, I22, K1, K2 = lambda_plus(eps2, Delta1, Delta2, mu, pc, eta, pmax)

    return  I11 + I12 * eps2 + (mu*K1 - Delta1) /Lambda_plus, Lambda_plus, I11, I12, I22, K1, K2

# #####################################
# Function to find epsilon_2 and lambda
# #####################################
def FindRoot(x1, step, tol, itermax, Delta1, Delta2, mu, pc, eta, pmax, verb):

    x2 = x1 + step
    f1, l1, i11, i12, i22, k1, k2 = F_plus(x1, Delta1, Delta2, mu, pc, eta, pmax)
    if verb:
        print ("x = %11.4e f = %11.4e" % (x1, f1))
    f2, l2, i11, i12, i22, k1, k2 = F_plus(x2, Delta1, Delta2, mu, pc, eta, pmax)
    if verb:
        print ("x = %11.4e f = %11.4e" % (x2, f2))

    feval = 2
    while True and feval < itermax:
        x = (f1*x2 - f2*x1) /(f1 - f2)
        f, l, i11, i12, i22, k1, k2 = F_plus(x, Delta1, Delta2, mu, pc, eta, pmax)
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

    return x, f, l, feval, i11, i12, i22, k1, k2       

# ######################
# Calculation parameters
# ######################
eps2_start = -0.3  # epsilon_2 scan start value
eps2_end   = +0.3  # epsilon_2 scan end value
Neps2      = 500   # Number of points in epssilon_2 scan

eta      = 1.e-6   # p integrand regularization parameter
pmax     = 20.     # Maximum p value 
tol      = 1.e-10  # Root finding accuracy parameter
step     = 1.e-3   # Initial step in root finding
itermax  = 15      # Maximum allowed number of root finding iterations

# ####################
# Get input parameters
# ####################
print ("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nTwo-Harmonic Rutherford Island Calculation\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
val1 = input("pc ? ")
pc = float(val1)

# #######
# Scan epsilon_2
# #######
ep2 = []
kk1 = []
kk2 = []
gg  = []
for eps2 in np.linspace(eps2_start, eps2_end, Neps2):

    K1  = K_m(eps2, 1, pc, eta)
    K2  = K_m(eps2, 2, pc, eta)
    G   = - eps2 * K1 /K2

    ep2.append(eps2)
    kk1.append(K1)
    kk2.append(K2)
    gg.append(G)
    print ("epsilon_2 = %11.4e K1 = %11.4e K2 = %11.4e G = %11.4e" \
           % (eps2, K1, K2, G))

# ###########
# Output data
# ###########
with open("Marginal.txt", "w") as file:
    for n in range(len(ep2)):
        file.write("%11.4e %11.4e %11.4e %11.4e\n" % (ep2[n], kk1[n], kk2[n], gg[n]))
    
# ##########    
# Graph data
# ##########
plt.figure(figsize=(10,7))
plt.rc('xtick', labelsize=16) 
plt.rc('ytick', labelsize=16) 
plt.xlim(eps2_start,eps2_end)
#plt.ylim(0., 1.3)
plt.plot(ep2, gg, color='blue',  linewidth='1.5')
plt.axhline(y=0., linewidth='1.5', linestyle='dotted', color='red')
plt.axvline(x=0., linewidth='1.5', linestyle='dotted', color='red')
plt.axvline(x=0.25, linewidth='1.5', linestyle='dotted', color='red')
plt.axvline(x=-0.25, linewidth='1.5', linestyle='dotted', color='red')
plt.xlabel("$\epsilon_2$", fontsize='20')
plt.ylabel("$G$", fontsize='20')

plt.show()
#plt.savefig("TwoHarmonic.eps")    
