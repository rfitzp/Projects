# Program to implement multi-harmonic Rutherford island width evolution 

import warnings
warnings.filterwarnings('ignore')

import scipy.integrate as integrate
import scipy.special as special
import numpy as np
import math
import matplotlib.pyplot as plt

kmax = 20.

def Cos_m_Integrand(x, k, m):
    if k < 1.:
        return math.cos(2.*m*math.acos(k*math.sin(x))) /math.sqrt(1.-k*k*math.sin(x)*math.sin(x)) /math.pi
    else:
        return math.cos(2.*m*math.acos(  math.sin(x))) /math.sqrt(k*k-math.sin(x)*math.sin(x))   /math.pi

def Cos_m(k, m):
    return integrate.quad(Cos_m_Integrand, 0., math.pi/2., args=(k, m))[0]

def X0(k):
    if k < 1.:
        return special.ellipk(k*k) /math.pi
    else:
        return special.ellipk(1./k/k) /k /math.pi

def M_Element_Integrand(k, m, mp):
    return 8. * Cosn(k, m) * Cos(k, mp) * k /X0(k)

def M_Element(m, mp):
    return integrate.quad(M_Element_Integrand, 0., kmax, args=(m, mp))[0]

