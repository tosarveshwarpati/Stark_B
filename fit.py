import numpy as np
import scipy.signal
import scipy
from scipy.integrate import quad
import matplotlib.pyplot as plt
import sys
obs = np.loadtxt("expt_data.txt")
x = obs[:,0]
data1 = np.loadtxt("K1DATA.txt")
mask1 = (data1[:,0]>566) & (data1[:,0]<870)
data = data1[mask1]
T = 1
l = data[:,0]
aki = data[:,1]
ei = data[:,2]
ek = data[:,3]
gi = data[:,4]
gk = data[:,5]
I =((aki*gk)/(l))*np.exp(-(ek)/(T*100))
#def _1gaussian(x,amp,cen,sigma):
 #   return amp*(np.exp((-1.0/2.0)*(((x-cen)/sigma)**2)))
b_d = []
b_d_g = []
b_d_l = []
b_l = []
bp_x = []
bp_y = []
area = 0
peak,_ = scipy.signal.find_peaks(I, prominence=1000)
for lam in l[peak]:
    popt_g = []
    pcov_g = []
    mask = (x>lam-.5) & (x<lam+.5)
    c_o = []
    dip,_ = scipy.signal.find_peaks(-obs[:,1][mask], prominence=100)
    dip_wl = obs[:,0][mask][dip]
    dip_int = obs[:,3][mask][dip]
    for ind in range(np.size(obs[mask][:,1])):
        base_l = (((dip_int[-1] - dip_int[0])/(dip_wl[-1] - dip_wl[0]))*(obs[mask][ind,0] - dip_wl[0])) + dip_int[0]
        correct_peak = obs[mask][ind,0],obs[mask][ind,1] - base_l
        c_o.append(correct_peak)
        correct_out = np.stack(c_o, axis = 0)
        correct_out = np.array(correct_out)
        b_l.append(base_l)
    plt.plot(obs[mask][:,0],obs[mask][:,1], label = 'obs')
    plt.plot(correct_out[:,0],correct_out[:,1], label = 'corrected')
    plt.plot(obs[mask][:,0],b_l, '--',  label = 'background')
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
    def gaussian(x_int,amp,cen,sigma):
        return amp*np.exp(-(x_int-cen)**2/(2*sigma**2))
    def lorentzian(x_int, a1, cen1, width1):
        return (a1*(width1**2)/(((x_int-cen1)**2)+width1**2))
    #plt.plot(lam,((baki*bgk)/(lam))*np.exp(-(bek)/(T*100)), "x")
    absolute_difference_function = lambda list_value : abs(list_value - lam)
    bl = min(correct_out[:,0], key=absolute_difference_function)
    bi = correct_out[:,1][correct_out[:,0]==bl]
    baki = aki[l==lam]
    bei = ei[l==lam]
    bek = ek[l==lam]
    bgi = gi[l==lam]
    bgk = gk[l==lam]
    b_l =[]
    bp_data = bek, np.log((area*lam)/(baki*bgk)) #area of corrected peak
    #bp_data = bek, np.log((bi*lam)/(baki*bgk))  #intensity of corrected peak
    amp = np.max(y_int)
    cen = lam
    sigma = x_int[-1]-x_int[0]
    while True:
        sel = input("Remove peak at "+str(lam)+ ":           [Input y]           ")
        if sel == 'y':
            print("This Line is removed")
            break
        else:
            popt_g, pcov_g = scipy.optimize.curve_fit(gaussian, x_int, y_int, p0 = [np.max(y_int), lam, x_int[-1]-x_int[0]], maxfev=500000)
            popt_lo, pcov_lo = scipy.optimize.curve_fit(lorentzian, x_int, y_int, p0 = [np.max(y_int), lam, x_int[-1]-x_int[0]], maxfev=500000)
            lor_err =  np.sqrt(np.diag(pcov_lo))
            gau_err =  np.sqrt(np.diag(pcov_g))
            err_g = abs(gau_err[0]/popt_g[0])+abs(2*(gau_err[1]/popt_g[1]))+abs(2*(gau_err[2]/popt_g[2]))
            err_lor = abs(lor_err[0]/popt_lo[0])+abs(lor_err[1]/popt_lo[1])+abs(4*(lor_err[2]/popt_lo[2]))
            lorentzian1 = lorentzian(correct_out[:,0], *popt_lo)
            gaussian1 = gaussian(x_int, *popt_g)
            plt.plot(correct_out[:,0],correct_out[:,1], label = 'corrected')
            plt.axvline(lam, label = 'Standerd Peak Position')
            gaussian2 = gaussian(correct_out[:,0], *popt_g)
            plt.plot(correct_out[:,0], gaussian2, label = 'Gaussian Fitting')
            plt.plot(correct_out[:,0], lorentzian1, label = 'Lorentzian Fitting')
            plt.legend(loc = 'upper left')
            plt.grid()
            plt.ylabel(r'$I \rightarrow$')
            plt.xlabel(r'$\lambda \longrightarrow$')
            plt.show()
            IG = quad(gaussian, correct_out[:,0][0],correct_out[:,0][-1], args = (popt_g[0], popt_g[1], popt_g[2]))
            IL = quad(lorentzian, correct_out[:,0][0],correct_out[:,0][-1], args = (popt_lo[0], popt_lo[1], popt_lo[2]))
            bp_data_L = bek, np.log((IL[0]*lam)/(baki*bgk))
            bp_data_G = bek, np.log((IG[0]*lam)/(baki*bgk))
            b_d_g.append(bp_data_G)
            b_d_l.append(bp_data_L)
            b_d.append(bp_data)
            boltzmann_g = np.stack(b_d_g, axis = 0)
            boltzmann_lo = np.stack(b_d_l, axis = 0)
            boltzmann = np.stack(b_d, axis = 0)
            with open('area'+str(lam)+'.txt', 'w') as f1:
                f1.write(str(area))
            with open('correct_out'+str(lam+1)+'.txt', 'w') as f2:
                f2.write(str(correct_out))
            break
bpx = np.array(boltzmann[:,0]).ravel()
bpy = np.array(boltzmann[:,1]).ravel()
order = np.argsort(bpx)
BPX = bpx[order]
BPY = bpy[order]
def sline(BPX, m, c):
	return m*BPX + c
popt_l, pcov_l = scipy.optimize.curve_fit(sline, BPX, BPY, p0=[BPY[0]/BPX[-1],BPY[0]])
lin1 = sline(BPX, *popt_l)
Te = -(1/popt_l[0])
#plt.plot(BPX,lin1, label = 'BP Linear Fit with $T_e$ ='+str(Te)+' eV')
#plt.plot(BPX,BPY,"x", label = 'BP points')
#plt.title("Boltzamnn Plot with $T_e$ ="+str(Te)+" eV")
#######################################################Gaussian
bpx_g = np.array(boltzmann_g[:,0]).ravel()
bpy_g = np.array(boltzmann_g[:,1]).ravel()
order_g = np.argsort(bpx_g)
BPX_g = bpx_g[order_g]
BPY_g = bpy_g[order_g]
def sline_g(BPX_g, m_g, c_g):
	return m_g*BPX_g + c_g
popt_l_g, pcov_l_g = scipy.optimize.curve_fit(sline_g, BPX_g, BPY_g, p0=[BPY_g[0]/BPX_g[-1],BPY_g[0]])
lin1_g = sline_g(BPX_g, *popt_l_g)
Te_g = -(1/popt_l_g[0])
plt.plot(BPX_g,lin1_g, label = 'BP Linear Fit(Gaussian Lineshape) with $T_e$ ='+str(Te_g)+' eV')
plt.plot(BPX_g,BPY_g,"bo", label = 'BP points(Gaussian)')
########################################################Lorentzian
bpx_lo = np.array(boltzmann_lo[:,0]).ravel()
bpy_lo = np.array(boltzmann_lo[:,1]).ravel()
order_lo = np.argsort(bpx_lo)
BPX_lo = bpx_lo[order_lo]
BPY_lo = bpy_lo[order_lo]
def sline_lo(BPX_lo, m_lo, c_lo):
	return m_lo*BPX_lo + c_lo 
popt_l_lo, pcov_l_lo = scipy.optimize.curve_fit(sline_lo, BPX_lo, BPY_lo, p0=[BPY_lo[0]/BPX_lo[-1],BPY_lo[0]])
lin1_lo = sline_lo(BPX_lo, *popt_l_lo)
Te_lo = -(1/popt_l_lo[0])
plt.plot(BPX_lo,lin1_lo, label = 'BP Linear Fit(Lorentzianian Lineshape) with $T_e$ ='+str(Te_lo)+' eV')
plt.plot(BPX_lo,BPY_lo,"ro", label = 'BP points(Lorentzian)')
########################################################
plt.title("Boltzamnn Plot")
plt.legend(loc = 'upper right')
plt.ylabel(r'$ln \frac{I\lambda}{A_{ki}g_k} \uparrow$',rotation=0)
plt.xlabel(r'$E_{k} \longrightarrow$', color='red')
plt.show()
BP_data = ek, np.log((I*l)/(aki*gk))












