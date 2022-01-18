import numpy as np
import scipy.signal
import scipy
import math
from scipy.integrate import quad_vec
import matplotlib.pyplot as plt
import sys
import scipy.constants as cst
wl = np.arange(200,1000,0.1)
data1 = np.loadtxt("K1NIST.txt")
a_PF = np.loadtxt("PF.txt")
mask1 = (data1[:,0]>200) & (data1[:,0]<1000)
data = data1[mask1]
#T = 1
lam = data[:,0]
aki = data[:,1]
f_ij = data[:,2]
ei = data[:,3]
ej = data[:,4]
g_i = data[:,5]
gj = data[:,6]
gauss_width = .0001
lorentzian_width = .0001
f_inst = 1000
l_abs = 1e+7
n_0 = 10e+16
#T = T_eV / 11600
T = 0.5
def Z(T):
	lnZ = sum([a_PF[i]*np.log(T*11604)**i for i in range(0,6)])
	return np.exp(lnZ)
for l in lam:
    ind = lam == l
    if np.size(lam[ind]) > 1:
        for dw in range(0, np.size(lam[ind])):
         def a(gauss_width, lorentzian_width):
             return(np.sqrt(np.log(2))*((lorentzian_width)/(gauss_width)))
         def x(wl, l, gauss_width):
             const = 2*np.sqrt(np.log(2))
             return np.array(const*(wl - l)/gauss_width)
         def k0(f_ij, n_0, g_i, ei, T, wl, gauss_width): 
             return np.array((2*((cst.e)**2)*np.sqrt(np.pi*np.log(2)*np.exp((-ei[ind][dw])/(100*T)))*n_0*g_i[ind][dw]*f_ij[ind][dw]*wl**2)/(cst.m_e*Z(T)*gauss_width*cst.c**2))
         def kappa(f_ij, n_0, g_i, ei, T, wl, gauss_width, lorentzian_width, l):
             k01 = k0(f_ij, n_0, g_i, ei, T, wl, gauss_width)
             a1 = a(gauss_width, lorentzian_width)
             x1 = x(wl, l, gauss_width)
             func = lambda t: (k01*a1/cst.pi)*((np.exp(-t**2))/(((t-x1)**2)+a1**2))
             kappa,_ = np.array(quad_vec(func, -np.inf, np.inf))
             return kappa
         def I(wl, l_abs, kappa, T, ej, f_inst, f_ij, n_0, g_i, ei, gauss_width, lorentzian_width, l):
             return np.array(((f_inst*8*cst.pi*cst.h*(cst.c)**2)/l**5)*(np.exp((ei[ind][dw]-ej[ind][dw])/(100*T)))*(1-np.exp(-kappa(f_ij, n_0, g_i, ei, T, wl, gauss_width, lorentzian_width, l)*l_abs)))
         print("abcd ",  l)
    else:
         def a(gauss_width, lorentzian_width):
             return(np.sqrt(np.log(2))*((lorentzian_width)/(gauss_width)))
         def x(wl, l, gauss_width):
             const = 2*np.sqrt(np.log(2))
             return np.array(const*(wl - l)/gauss_width)
         def k0(f_ij, n_0, g_i, ei, T, wl, gauss_width): 
             return np.array((2*((cst.e)**2)*np.sqrt(np.pi*np.log(2)*np.exp((-ei[ind])/(100*T)))*n_0*g_i[ind]*f_ij[ind]*wl**2)/(cst.m_e*Z(T)*gauss_width*cst.c**2))
         def kappa(f_ij, n_0, g_i, ei, T, wl, gauss_width, lorentzian_width, l):
             k01 = k0(f_ij, n_0, g_i, ei, T, wl, gauss_width)
             a1 = a(gauss_width, lorentzian_width)
             x1 = x(wl, l, gauss_width)
             func = lambda t: (k01*a1/cst.pi)*((np.exp(-t**2))/(((t-x1)**2)+a1**2))
             kappa,_ = np.array(quad_vec(func, -np.inf, np.inf))
             return kappa
         def I(wl, l_abs, kappa, T, ej, f_inst, f_ij, n_0, g_i, ei, gauss_width, lorentzian_width, l):
             return np.array(((f_inst*8*cst.pi*cst.h*(cst.c)**2)/l**5)*(np.exp((ei[ind]-ej[ind])/(100*T)))*(1-np.exp(-kappa(f_ij, n_0, g_i, ei, T, wl, gauss_width, lorentzian_width, l)*l_abs)))
         #print(1-np.exp(-kappa(f_ij, n_0, g_i, ei, T, wl, gauss_width, lorentzian_width, l)*l_abs))
         #print(I(wl, l_abs, kappa, T, ej, f_inst, f_ij, n_0, g_i, ei, gauss_width, lorentzian_width, l))
         plt.plot(wl, I(wl, l_abs, kappa, T, ej, f_inst, f_ij, n_0, g_i, ei, gauss_width, lorentzian_width, l))
###############################################
#generation
###############################################
    
    
###############################################
#fitting
###############################################





















