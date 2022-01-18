import numpy as np
import scipy.signal
import scipy
import os
from scipy.integrate import quad
import matplotlib.pyplot as plt
import sys
obs = np.loadtxt("expt_data.txt")
x = obs[:,0]
data1 = np.loadtxt("K1DATA.txt")
mask1 = (data1[:,0]>566) & (data1[:,0]<870)
data = data1[mask1]
T = 1
a0 = 0.56534
a1 = 0.0245
a2 = 0.0157
Stark_w = .6240e-1
l = data[:,0]
aki = data[:,1]
ei = data[:,2]
ek = data[:,3]
gi = data[:,4]
gk = data[:,5]
I =((aki*gk)/(l))*np.exp(-(ek)/(T*100))
peak,_ = scipy.signal.find_peaks(I, prominence=5000)
for lam in l[peak]:
    Ei = ei[l == lam]/8065.73
    Ek = ek[l == lam]/8065.73
    mask = (x>lam-1.0) & (x<lam+1.0)
    c_o = []
    b_l = []
    area = 0
    dip,_ = scipy.signal.find_peaks(-obs[:,1][mask], prominence=500)
    dip_wl = obs[:,0][mask][dip]
    dip_int = obs[:,1][mask][dip]
    import numpy as np
    for ind in range(np.size(obs[mask][:,1])):
        base_l = (((dip_int[-1] - dip_int[0])/(dip_wl[-1] - dip_wl[0]))*(obs[mask][ind,0] - dip_wl[0])) + dip_int[0]
        correct_peak = obs[mask][ind,0],obs[mask][ind,1] - base_l
        c_o.append(correct_peak)
        correct_out = np.stack(c_o, axis = 0)
        correct_out = np.array(correct_out)
        b_l.append(base_l)
    plt.plot(obs[mask][:,0],obs[mask][:,1], label = 'obs')
    plt.plot(correct_out[:,0],correct_out[:,1], label = 'corrected')
    #plt.plot(obs[mask][:,0],b_l, '--',  label = 'background')
    plt.axvline(lam, label = 'Standerd Peak Position')
    maskl = correct_out[:,0] <= lam
    maskr = correct_out[:,0] >= lam
    dipl,_ = scipy.signal.find_peaks(-correct_out[:,1][maskl])
    dipr,_ = scipy.signal.find_peaks(-correct_out[:,1][maskr])
    plt.legend(loc = 'upper left')
    lim_left = dipl[-1]
    lim_right = dipr[0]
    mask_int = (correct_out[:,0][maskl][lim_left] <= correct_out[:,0])&(correct_out[:,0] <= correct_out[:,0][maskr][lim_right])
    x_int = correct_out[:,0][mask_int]
    y_int = correct_out[:,1][mask_int]
    for x1 in range(np.size(correct_out[:,0][mask_int]) - 1):
        area1 = ((x_int[x1+1]-x_int[x1])*y_int[x1+1]) - ((1/2)*(x_int[x1+1]-x_int[x1])*(y_int[x1+1]-y_int[x1]))
        area += area1
    plt.fill_between(correct_out[:,0][mask_int],correct_out[:,1][mask_int], color = 'gray')
    plt.axvline(correct_out[:,0][maskl][lim_left], color = 'r')
    plt.axvline(correct_out[:,0][maskr][lim_right], color = 'r')
    plt.title("Peak at $\lambda = $"+str(lam))
    plt.ylabel(r'$I \rightarrow$')
    plt.xlabel(r'$\lambda \longrightarrow$')
    plt.grid()
    plt.show()
    while True:
        sel = input("Remove peak at "+str(lam)+ ":           [Input y]           ")
        if sel == 'y':
            print("This Line is removed")
            break
        else:
            param_bounds=([-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,0],[np.inf,np.inf,np.inf,np.inf,np.inf,1])
            p = 0.5
            popt_pv = []
            import numpy as np
            [amp_g, amp_l, cen, width_g, width_l, p] = [np.max(y_int),np.max(y_int), lam, x_int[-1]-x_int[0] , x_int[-1]-x_int[0], .5]
            def pv(x_int, amp_g, amp_l, cen, width_g, width_l, p ):
                return p*((amp_g/np.sqrt(2*np.pi*width_g**2))*np.exp(-((x_int-cen)**2)/(2*width_g**2)))+(1-p)*((amp_l/np.pi)*((width_l/2)**2/((x_int-cen)**2+(width_l/2)**2)))
            popt_pv,_ = scipy.optimize.curve_fit(pv, x_int, y_int, p0 = [np.max(y_int),np.max(y_int), lam, x_int[-1]-x_int[0] , x_int[-1]-x_int[0], p] ,bounds = param_bounds, maxfev=500000)
            x_int_smooth = np.arange(x_int[0],x_int[-1],0.0001)
            PV = pv(x_int_smooth, *popt_pv)
            C = 4.700e+21
            FWHM_g=popt_pv[3]*np.sqrt(4*np.log(2))
            FWHM_l=popt_pv[4]
            N_e = (FWHM_l*(10**16))/(2*Stark_w)
            #N_e = (C)/((popt_pv[2]/10)**2*(10**(-8))*(Ek-Ei))
            def T_e(FWHM_l):
                return np.exp((-a1+np.sqrt(np.abs((a1**2)-(4*a2*(a0-np.log(FWHM_l))))))/(2*a2))
            plt.plot(x_int, y_int)
            plt.plot(x_int_smooth, PV)
            plt.show()
            while True:
                sel1 = input("do you want to change following parameters    :    [Input y] \n amp_g = "+str(popt_pv[0])+" \n amp_l = "+str(popt_pv[1])+"\n cen = "+str(popt_pv[2])+" \n width_g = "+str(popt_pv[3])+" \n width_l = "+str(popt_pv[4])+" \n p = "+str(popt_pv[5])+"\n \n               Input response: ")
                if sel1 == 'y':
                    k = input("press the index of parameter to be changed(1 for amp_g, 2 for amp_l and so on... ...):")
                    if k == '1':
                        param = input("input the new value of parameter:  ")
                        opt_param = popt_pv
                        opt_param[0] = param
                        PV = pv(x_int, *opt_param)
                        plt.plot(x_int, y_int)
                        plt.plot(x_int, PV)
                        plt.show()
                        continue
                    elif k == '2':
                        param = input("input the new value of parameter:  ")
                        opt_param = popt_pv
                        opt_param[1] = param
                        PV = pv(x_int, *opt_param)
                        plt.plot(x_int, y_int)
                        plt.plot(x_int, PV)
                        plt.show()
                        continue
                    elif k == '3':
                        param = input("input the new value of parameter:  ")
                        opt_param = popt_pv
                        opt_param[2] = param
                        PV = pv(x_int, *opt_param)
                        plt.plot(x_int, y_int)
                        plt.plot(x_int, PV)
                        plt.show()
                        continue
                    elif k == '4':
                        param = input("input the new value of parameter:  ")
                        opt_param = popt_pv
                        opt_param[3] = param
                        PV = pv(x_int, *opt_param)
                        plt.plot(x_int, y_int)
                        plt.plot(x_int, PV)
                        plt.show()
                        continue
                    elif k == '5':
                        param = input("input the new value of parameter:  ")
                        opt_param = popt_pv
                        opt_param[4] = param
                        PV = pv(x_int, *opt_param)
                        plt.plot(x_int, y_int)
                        plt.plot(x_int, PV)
                        plt.show()
                        continue
                    elif k == '6':
                        param = input("input the new value of parameter:  ")
                        opt_param = popt_pv
                        opt_param[5] = param
                        PV = pv(x_int, *opt_param)
                        plt.plot(x_int, y_int)
                        plt.plot(x_int, PV)
                        plt.show()
                        continue
                    elif k == 'n':
                        break
#p0 = [np.max(y_int),np.max(y_int), lam, x_int[-1]-x_int[0]] , x_int[-1]-x_int[0], p] , 


