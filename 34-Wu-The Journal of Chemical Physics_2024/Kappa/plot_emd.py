# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 11:45:03 2021

@author: benward (hityingph@163.com)
"""
from scipy import integrate
from pylab import *
from glob import glob

##set figure properties
aw = 1.5
fs = 16
lw = 2
ms = 0
font = {'size'   : fs}
matplotlib.rc('font', **font)
matplotlib.rc('axes' , lw=aw)

def set_fig_properties(ax_list):
    tl = 6
    tw = 1.5
    tlm = 3
    
    for ax in ax_list:
        ax.tick_params(which='major', length=tl, width=tw)
        ax.tick_params(which='minor', length=tlm, width=tw)
        ax.tick_params(which='both', axis='both', direction='out', right=False, top=False)

def calc_kappa(hacfile):
    ### Obtain the hac data
    with open(hacfile, "r") as fin:
        raw_data = fin.readlines()                      # read the original data
    
    raw_data = raw_data[3:]                             # skip the first three head lines
    T_start = int(raw_data[0].split()[0])
    Nc = int(raw_data[0].split()[1])
    
    # skip the initial zero values, we can not calculate thermal conductivity when the correlation time equals zero.
    raw_data = raw_data[Nc+1:]                           
    nall  = int(len(raw_data)/Nc)
    
    hac_matrix = np.zeros((Nc, nall))                   # obtain the hac matrix
    for i in range(nall):
        # print("\n")
        for j in range(Nc):
            line = raw_data[i*(Nc+1)+j+1]
            line = line.split()
            hac_xx = float(line[-3])
            hac_yy = float(line[-2])
            hac_zz = float(line[-1])
            hac_xyz    = 1/3*(hac_xx + hac_yy + hac_zz)
            # print(hac_ave)
            hac_matrix[j, i] = hac_xyz                  # each column denotes each independent simulations
            
    hac_ave = np.mean(hac_matrix, axis=1)

       
    ### Calculate the thermal conductivity
    # some parameters in MD simulations
    dt = 0.001                                          # tim step, unit in ps
    T  = 300                                            # temperature, unit in K
    kB = 8.617e-5                                       # Boltzmann’s constant, unit in eV/K
    V  = 54.74**3                                       # volume, unit in \AA**3
    Ns = 10                                             # sampling range
    scale = 1600 / (kB*T*T*V)*(Ns*dt)                   # we use the real unit in lammps, be careful for the unit conversion 
    
    t = np.arange(1,Nc+1)*Ns*dt                         # the correlation time
    rtc = scale*integrate.cumtrapz(hac_matrix, axis =0, initial=0) # obtain the running thermal conductivity for all simulations
    rtc_ave = np.mean(rtc, axis=1)
    return t, rtc

def output_kappa(path):
    fs = glob(f"{path}/hac_*.txt")
    N = len(fs)
    print(f"We have performed {N} independent runs!")
    rtcs = np.zeros((100000, N))
    for i in range(N):
        t, rtc = calc_kappa(fs[i])
        rtcs[:, i] = rtc.squeeze()
    rtc_ave = np.mean(rtcs, axis=1)
    rtc_error =  np.std(rtcs, axis=1)/np.sqrt(N)
    return N, rtcs, rtc_ave, rtc_error
    
nhc_N, nhc_rtcs, nhc_rtc_ave, nhc_rtc_error = output_kappa("NHC")
lan40_N, lan40_rtcs, lan40_rtc_ave, lan40_rtc_error = output_kappa("LAN_40")
lan100_N, lan100_rtcs, lan100_rtc_ave, lan100_rtc_error = output_kappa("LAN_100")
lan250_N, lan250_rtcs, lan250_rtc_ave, lan250_rtc_error = output_kappa("LAN_250")
lan350_N, lan350_rtcs, lan350_rtc_ave, lan350_rtc_error = output_kappa("LAN_350")

t = np.arange(1, 100000+1) * 10 / 1000
figure(figsize=(10, 20))
subplot(5, 1, 1)
set_fig_properties([gca()])
for i in range(nhc_N):
    plot(t, nhc_rtcs[:,i], lw = 0.5, alpha=0.5, c = "grey")
plot(t, nhc_rtc_ave, lw=lw, c = "C0", label="NHC")
plot(t, nhc_rtc_ave + nhc_rtc_error, lw=1, ls = "--", c = "black")
plot(t, nhc_rtc_ave - nhc_rtc_error, lw=1, ls = "--", c = "black")
xlim([0, 1000])
gca().set_xticks([])
ylim([0, 100])
gca().set_yticks(linspace(10, 90, 3))
text(700, 10, f"{np.mean(nhc_rtc_ave[60000]):.1f} ({np.mean(nhc_rtc_error[60000]):.1f}) W/m/K")
legend(loc="upper left")

subplot(5, 1, 2)
set_fig_properties([gca()])
for i in range(lan350_N):
    plot(t, lan350_rtcs[:,i], lw = 0.5, alpha=0.5, c = "grey")
plot(t, lan350_rtc_ave, lw=lw, c = "C2", label="NHC + LAN (350 ps)")
plot(t, lan350_rtc_ave + lan350_rtc_error, lw=1, ls = "--", c = "black")
plot(t, lan350_rtc_ave - lan350_rtc_error, lw=1, ls = "--", c = "black")
xlim([0, 1000])
gca().set_xticks([])
ylim([0, 100])
gca().set_yticks(linspace(10, 90, 3))
text(700, 10, f"{np.mean(lan350_rtc_ave[50000]):.1f} ({np.mean(lan350_rtc_error[50000]):.1f}) W/m/K")
legend(loc="upper left")

subplot(5, 1, 3)
set_fig_properties([gca()])
for i in range(lan250_N):
    plot(t, lan250_rtcs[:,i], lw = 0.5, alpha=0.5, c = "grey")
plot(t, lan250_rtc_ave, lw=lw, c = "C1", label="NHC + LAN (250 ps)")
plot(t, lan250_rtc_ave + lan250_rtc_error, lw=1, ls = "--", c = "black")
plot(t, lan250_rtc_ave - lan250_rtc_error, lw=1, ls = "--", c = "black")
xlim([0, 1000])
gca().set_xticks([])
ylim([0, 100])
gca().set_yticks(linspace(10, 90, 3))
text(700, 10, f"{np.mean(lan250_rtc_ave[40000]):.1f} ({np.mean(lan250_rtc_error[40000]):.1f}) W/m/K")
ylabel('Thermal conductivity (W/m/K))')
legend(loc="upper left")

subplot(5, 1, 4)
set_fig_properties([gca()])
for i in range(lan100_N):
    plot(t, lan100_rtcs[:,i], lw = 0.5, alpha=0.5, c = "grey")
plot(t, lan100_rtc_ave, lw=lw, c = "C2", label="NHC + LAN (100 ps)")
plot(t, lan100_rtc_ave + lan100_rtc_error, lw=1, ls = "--", c = "black")
plot(t, lan100_rtc_ave - lan100_rtc_error, lw=1, ls = "--", c = "black")
xlim([0, 1000])
gca().set_xticks([])
ylim([0, 100])
gca().set_yticks(linspace(10, 90, 3))
text(700, 10, f"{np.mean(lan100_rtc_ave[30000]):.1f} ({np.mean(lan100_rtc_error[30000]):.1f}) W/m/K")
legend(loc="upper left")

subplot(5, 1, 5)
set_fig_properties([gca()])
for i in range(lan40_N):
    plot(t, lan40_rtcs[:,i], lw = 0.5, alpha=0.5, c = "grey")
plot(t, lan40_rtc_ave, lw=lw, c = "C1", label="NHC + LAN (40 ps)")
plot(t, lan40_rtc_ave + lan40_rtc_error, lw=1, ls = "--", c = "black")
plot(t, lan40_rtc_ave - lan40_rtc_error, lw=1, ls = "--", c = "black")
xlim([0, 1000])
gca().set_xticks(linspace(0, 1000, 6))
ylim([0, 100])
gca().set_yticks(linspace(10, 90, 3))
text(700, 10, f"{np.mean(lan40_rtc_ave[20000]):.1f} ({np.mean(lan40_rtc_error[20000]):.1f}) W/m/K")
xlabel('Correlation Time (ps)')
legend(loc="upper left")


subplots_adjust(wspace=0.3, hspace=0)
savefig("Si_DP_EMD.png", bbox_inches='tight')

