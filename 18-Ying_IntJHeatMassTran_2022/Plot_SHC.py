# -*- coding: utf-8 -*-
"""
Created on Mon May 16 15:09:47 2022

@author: benward
"""
from thermo.gpumd.data import load_shc, load_compute
from pylab import *

def hnemd_shc_Nc(Nc, path, Fe):
    shc = load_shc(Nc=[500]*5, num_omega=[500]*5, directory=path)
    for keys in shc:
        shc[keys]["K"] = (shc[keys]['Ki']+shc[keys]['Ko'])
    shc_K = np.zeros((len(shc["run0"]["K"]),5))
    shc_J = np.zeros((len(shc["run0"]["jwi"]),5))
    i = 0
    for keys in shc:
        shc_K[:,i] = shc[keys]["K"]
        shc_J[:,i] = shc[keys]['jwi'] + shc[keys]['jwo']
        i += 1
    shc_K_ave = np.average(shc_K, axis=1)
    shc_J_ave = np.average(shc_J, axis=1)
    ### This affects the resolution 
    Nc = Nc
    tc = shc['run0']['t'][500 - Nc: Nc - 500] 
    K = shc_K_ave[500 - Nc: Nc - 500]
    if Nc == 500:
        tc = shc['run0']['t']
        K = shc_K_ave
    # plot(tc, K)
    ### Hann window
    K = K*(np.cos(np.pi*np.arange(-Nc + 1, Nc) / Nc) + 1) * 0.5
    # plot(tc, K)
    ### the Fourier transform
    nu = np.arange(0.1, 16.1, 0.1)
    q = np.zeros(len(nu))
    for n in range(len(nu)):
        q[n] = 2 * 0.01 * sum(K * np.cos(2 * np.pi * nu[n] * np.arange(-Nc + 1, Nc) * 0.01))
    q[q < 0] = 0
    with open(path + "xyz.in", "r") as fin:
        fin.readline()
        line = fin.readline().split()
    lx = float(line[-3])
    ly = float(line[-2])
    lz = float(line[-1])
    V = lx*ly*lz 
    T = 300
    Fe = Fe
    convert = 1602.17662
    shc = q*convert/(Fe*T*V)
    # plot(nu, shc)
    return shc

def nemd_shc_Nc(Nc, path, l_transport, axis):
    ### Obtain the original K and correlation time
    shc = load_shc(Nc=[500]*5, num_omega=[500]*5, directory=path)
    for keys in shc:
        shc[keys]["K"] = (shc[keys]['Ki']+shc[keys]['Ko'])
    shc_K = np.zeros((len(shc["run0"]["K"]),5))
    i = 0
    for keys in shc:
        shc_K[:,i] = shc[keys]["K"]
        i += 1
    shc_K_ave = np.average(shc_K, axis=1)
    ### This affects the resolution 
    Nc = Nc
    tc = shc['run0']['t'][500 - Nc: Nc - 500] 
    K = shc_K_ave[500 - Nc: Nc - 500]
    if Nc == 500:
        tc = shc['run0']['t']
        K = shc_K_ave
    ### Hann window
    K = K*(np.cos(np.pi*np.arange(-Nc + 1, Nc) / Nc) + 1) * 0.5
    ### the Fourier transform
    nu = np.arange(0.1, 16.1, 0.1)
    q = np.zeros(len(nu))
    for n in range(len(nu)):
        q[n] = 2 * 0.01 * sum(K * np.cos(2 * np.pi * nu[n] * np.arange(-Nc + 1, Nc) * 0.01))
    with open(path + "xyz.in", "r") as fin:
        fin.readline()
        line = fin.readline().split()
    lx = float(line[-3])
    ly = float(line[-2])
    lz = float(line[-1])
    if axis == "x":
        A = ly*lz  
    elif axis == "y":
        A = lx*lz   
    # print(lz)
    compute = load_compute(['T'], directory = path)
    T = compute['T']
    ndata = int(T.shape[0]/5)
    # print(ndata)
    temp_ave = mean(T[int(ndata/2)+1:, 1:], axis=0)
    deltaT = temp_ave[0] - temp_ave[-1]
    # print(deltaT)
    shc = 1.6e4*q/A/deltaT/l_transport
    Ein = compute['Ein'][:ndata]
    Eout= compute['Eout'][:ndata]
    dt = 0.001
    Ns = 1000
    t = dt*np.arange(1,ndata+1)
    Q1 = (Ein[int(ndata/2)] - Ein[-1])/(ndata/2)/dt/Ns
    Q2 = (Eout[-1] - Eout[int(ndata/2)])/(ndata/2)/dt/Ns
    # plot(t, Ein)
    # plot(t, Eout)
    Q = mean([Q1, Q2])
    G_thermo = 1.6e4*Q/A/deltaT
    error = format(trapz(shc, nu) - G_thermo, ".3f")
    print("Error between direct method and spectral method: %s GW/m2/K"%error)
    return shc

def calc_lambda(shc_hnemd, shc_nemd):
    lambda_i = shc_hnemd/shc_nemd #unit in nm
    length = logspace(0.1, 6, 100) #consider length from 10 nm to 1 mm
    k_L = np.zeros_like(length)
    nu = np.arange(0.1, 16.1, 0.1)
    for i, el in enumerate(length):
        k_L[i] = np.trapz(shc_hnemd/(1+lambda_i/el), nu)
    output = dict()
    output["L"] = length
    output["lambda"] = lambda_i
    output["K"] = k_L
    return output

aw = 1.5
fs = 18
font = {'size'   : fs}
matplotlib.rc('font', **font)
matplotlib.rc('axes' , linewidth=aw)

def set_fig_properties(ax_list):
    tl = 6
    tw = 1.5
    tlm = 3
    
    for ax in ax_list:
        ax.tick_params(which='major', length=tl, width=tw)
        ax.tick_params(which='minor', length=tlm, width=tw)
        ax.tick_params(which='both', axis='both', direction='out', right=False, top=False)
        



figure(figsize=(14, 12))


black_arm = hnemd_shc_Nc(100, "HNEMD_diffusive/Black/300K/x_1e-4/", 1e-4)
black_zig = hnemd_shc_Nc(100, "HNEMD_diffusive/Black/300K/y_1e-4/", 1e-4)
blue_arm =  hnemd_shc_Nc(100, "HNEMD_diffusive/Blue/300K/x_1e-5/", 1e-5) 
blue_zig =  hnemd_shc_Nc(100, "HNEMD_diffusive/Blue/300K/y_1e-5/", 1e-5) 
violet_x = hnemd_shc_Nc(100, "HNEMD_diffusive/Violet/300K/x_1e-4/", 1e-4)
violet_y = hnemd_shc_Nc(100, "HNEMD_diffusive/Violet/300K/y_1e-4/", 1e-4)
blue = (blue_arm + blue_zig) / 2
violet = (violet_x + violet_y) / 2
black_arm[78: 105] = 0
black_zig[78: 105] = 0
blue[76:123] = 0

subplot(2,2,1)
set_fig_properties([gca()])
nu = np.arange(0.1, 16.1, 0.1)
plot(nu, black_arm, c = "black", ls = "-", lw = 2, label = "Black-P armchair")
plot(nu, black_zig, c = "black", ls = "--", lw = 2, label = "Black-P zigzag")
plot(nu, blue, c = "blue", ls = "-.", lw = 2, label = "Blue-P")
plot(nu, violet, c = "purple", ls = ":", lw = 3, label = "Violet-P")
xlim([0, 16])
gca().set_xticks(linspace(0,16,5))
ylim([1e-2, 1e2])
yscale("log")
ylabel(r'$\kappa$($\omega$) (W m$^{-1}$ K$^{-1}$ THz$^{-1}$)')
xlabel(r'$\omega$/2$\pi$ (THz)')
legend(frameon = False, fontsize = 18, labelspacing = 0.1)
text(1, 50, "(a)", fontsize = 18)


black_arm1 = nemd_shc_Nc(100, "NEMD_ballstic/Black/x/", 13.576, "x")   
black_zig1 = nemd_shc_Nc(100, "NEMD_ballstic/Black/y/", 13.181, "y")   
blue_arm1 = nemd_shc_Nc(100, "NEMD_ballstic/Blue/x/", 13.051, "x") 
blue_zig1 = nemd_shc_Nc(100, "NEMD_ballstic/Blue/y/", 11.300, "y") 
violet_x1 = nemd_shc_Nc(100, "NEMD_ballstic/Violet/x/", 18.250, "x")  
violet_y1 = nemd_shc_Nc(100, "NEMD_ballstic/Violet/y/", 18.364, "y") 
blue1 = (blue_arm1 + blue_zig1) / 2
violet1 = (violet_x1 + violet_y1) / 2
black_arm_nemd = np.array([i for i in black_arm1])
black_arm_nemd[78: 105] = 0
black_zig_nemd = np.array([i for i in black_zig1])
black_zig_nemd[78: 105] = 0
blue_nemd = np.array([i for i in blue1])
blue_nemd[76:123] = 0


subplot(2, 2, 2)
set_fig_properties([gca()])
plot(nu, black_arm_nemd, c = "black", ls = "-", lw = 2, label = "Black armchair")
plot(nu, black_zig_nemd, c = "black", ls = "--", lw = 2, label = "Black zigzag")
plot(nu, blue_nemd, c = "blue", ls = "-.", lw = 2, label = "Blue")
plot(nu, violet1, c = "purple", ls = ":", lw = 3, label = "Violet")
xlim([0, 16])
gca().set_xticks(linspace(0,16,5))
ylim([1e-4, 1])
yscale("log")
ylabel(r'G($\omega$) (GW m$^-2$ K$^{-1}$ THz$^{-1}$)')
xlabel(r'$\omega$/2$\pi$ (THz)')
text(1, 0.5, "(b)", fontsize = 18)


black_arm2 = calc_lambda(black_arm, black_arm1)
black_zig2 = calc_lambda(black_zig, black_zig1)
blue_arm2 = calc_lambda(blue_arm, blue_arm1)
blue_zig2 = calc_lambda(blue_zig, blue_zig1)
violet_x2 = calc_lambda(violet_x, violet_x1)
violet_y2 = calc_lambda(violet_y, violet_y1)

black_arm_mfp = black_arm2['lambda']
black_zig_mfp = black_zig2['lambda']
black_arm_mfp[78: 105] = 0
black_zig_mfp[78: 105] = 0
blue_mfp = (blue_arm2['lambda'] + blue_zig2['lambda'])/2
blue_mfp[76:123] = 0
violet_mfp = (violet_x2['lambda'] + violet_y2['lambda'])/2
subplot(2,2,3) 
set_fig_properties([gca()])
plot(nu, black_arm_mfp, c = "black", ls = "-", lw = 2, label = "Black armchair")
plot(nu, black_zig_mfp, c = "black", ls = "--", lw = 2, label = "Black zigzag")
plot(nu, blue_mfp, c = "blue", ls = "-.", lw = 2, label = "Blue")
plot(nu, violet_mfp, c = "purple", ls = ":", lw = 3, label = "Violet")
xlim([0, 16])
gca().set_xticks(linspace(0, 16, 5))
ylim([1, 10000])
# gca().set_yticks(linspace(0, 10000, 6))
yscale("log")
ylabel(r'$\lambda$($\omega$) (nm)')
xlabel(r'$\omega$/2$\pi$ (THz)')
text(1, 5000, "(c)", fontsize = 18)

subplot(2,2,4) 
set_fig_properties([gca()])
length = logspace(0.1, 6, 100)/1000
blue_kl = (blue_arm2['K'] + blue_zig2['K'])/2
violet_kl = (violet_x2['K'] + violet_y2['K'])/2
semilogx(length, black_arm2['K'], c = "black", ls = "-", lw = 2, label = "Black armchair")
semilogx(length[42], black_arm2['K'][42], c = "black", lw = 0, marker = "o", markersize = 10)
semilogx(length, black_zig2['K'], c = "black", ls = "--", lw = 2, label = "Black zigzag")
semilogx(length[53], black_zig2['K'][53], c = "black", lw = 0, marker = "s", markersize = 10)
semilogx(length, blue_kl, c = "blue", ls = "-.", lw = 2, label = "Blue")
semilogx(length[67], blue_kl[67], c = "blue", lw = 0, marker = "p", markersize = 10)
semilogx(length, violet_kl, c = "purple", ls = ":", lw = 3, label = "Violet")
semilogx(length[40], violet_kl[40], c = "purple", lw = 0, marker = "^", markersize = 10)
xlim([1e-3, 1e3])
ylim([0.1, 500])
yscale("log")
ylabel(r'$\kappa$ (W m$^{-1}$ K$^{-1}$)')
xlabel(r'L ($\mu$m)')
text(0.003, 250, "(d)", fontsize = 18)

subplots_adjust(wspace=0.35, hspace=0.25)
savefig("./fig6-shc.pdf", bbox_inches='tight')



