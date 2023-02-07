# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 23:17:35 2021

@author: benward
"""

from pylab import *
from thermo.gpumd.data import load_kappa, load_shc, load_compute
from thermo.gpumd.calc import running_ave, hnemd_spectral_kappa

def kappa_std(kappa):
    std = []
    for i in range(len(kappa)):
        std.append(np.std(kappa[i])/sqrt(10))
    return std

time_step = 1.0  #fs
output_inter = 1000

def output_kappa(path, direction, output_num):
    kappa_raw = np.loadtxt(path + "kappa.out")
    run_num = int(len(kappa_raw)/output_num)
    kappa = dict()
    for i in range(run_num):
        data = kappa_raw[output_num*i:output_num*(i+1)]
        out = dict()
        labels = ['kxi', 'kxo', 'kyi', 'kyo', 'kz']
        for j, key in enumerate(labels):
            out[key] = data[:,j]
        kappa["run%s"%i] = out
    
    kappa_in = np.zeros((output_num, run_num))
    kappa_out = np.zeros((output_num, run_num))
    kappa_tol = np.zeros((output_num, run_num))
    t = np.arange(1,kappa["run0"]['kxi'].shape[0]+1)*output_inter*time_step/1000  # ps
    i = 0
    if direction == "x":
        for keys in kappa:
            kappa[keys]['kxi_ra'] = running_ave(kappa[keys]['kxi'],t)
            kappa[keys]['kxo_ra'] = running_ave(kappa[keys]['kxo'],t)
            kappa[keys]['kx_ra'] = running_ave(kappa[keys]['kxi'],t) + running_ave(kappa[keys]['kxo'],t)
            kappa_in[:,i] = kappa[keys]['kxi_ra']
            kappa_out[:,i] = kappa[keys]['kxo_ra']
            kappa_tol[:,i] = kappa[keys]['kx_ra']
            i += 1
        
        kappa_in_ave = np.average(kappa_in, axis = 1)
        kappa_out_ave = np.average(kappa_out, axis = 1)
        kappa_tol_ave = np.average(kappa_tol, axis = 1)
        
        kappa_in_std = kappa_std(kappa_in)
        kappa_out_std = kappa_std(kappa_out)
        kappa_tol_std = kappa_std(kappa_tol)
        
        kappa_ave = dict()
        kappa_ave["in"] = kappa_in_ave
        kappa_ave["in_std"] = kappa_in_std
        kappa_ave["out"] = kappa_out_ave
        kappa_ave["out_std"] = kappa_out_std
        kappa_ave["tol"] = kappa_tol_ave
        kappa_ave["tol_std"] = kappa_tol_std
    elif direction == "y":
        for keys in kappa:
            kappa[keys]['kyi_ra'] = running_ave(kappa[keys]['kyi'],t)
            kappa[keys]['kyo_ra'] = running_ave(kappa[keys]['kyo'],t)
            kappa[keys]['ky_ra'] = running_ave(kappa[keys]['kyi'],t) + running_ave(kappa[keys]['kyo'],t)
            kappa_in[:,i] = kappa[keys]['kyi_ra']
            kappa_out[:,i] = kappa[keys]['kyo_ra']
            kappa_tol[:,i] = kappa[keys]['ky_ra']
            i += 1
        
        kappa_in_ave = np.average(kappa_in, axis = 1)
        kappa_out_ave = np.average(kappa_out, axis = 1)
        kappa_tol_ave = np.average(kappa_tol, axis = 1)
        
        kappa_in_std = kappa_std(kappa_in)
        kappa_out_std = kappa_std(kappa_out)
        kappa_tol_std = kappa_std(kappa_tol)
        
        kappa_ave = dict()
        kappa_ave["in"] = kappa_in_ave
        kappa_ave["in_std"] = kappa_in_std
        kappa_ave["out"] = kappa_out_ave
        kappa_ave["out_std"] = kappa_out_std
        kappa_ave["tol"] = kappa_tol_ave
        kappa_ave["tol_std"] = kappa_tol_std
    t = t * 0.001
    print(path + "\n")
    print(str(run_num) + "\n")
    print("kin = " + format(kappa_in_ave[-1], ".3f") + " ± " 
          + format(kappa_in_std[-1], ".3f") + "\n")
    print("kout = " + format(kappa_out_ave[-1], ".3f") + " ± " 
          + format(kappa_out_std[-1], ".3f") + "\n")
    print("ktol = " + format(kappa_tol_ave[-1], ".3f") + " ± " 
          + format(kappa_tol_std[-1], ".3f") + "\n")
    
    return t, kappa, kappa_ave


   
aw = 1.5
fs = 20
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
c_in = (0/255, 112/255, 177/255)
c_out = (218/255, 66/255, 72/255)
c_tol = (18/255, 156/255, 71/255)   
     
figure(figsize=(12, 12))

t, black_arm, black_arm_ave = output_kappa("HNEMD_diffusive/Black/300K/x_1e-4/", "x", 2000)
t, black_zig, black_zig_ave = output_kappa("HNEMD_diffusive/Black/300k/y_1e-4/", "y", 2000)
subplot(3, 2, 1)
set_fig_properties([gca()])
plot(t, black_arm_ave["in"], linewidth=2.0, color=c_in, label = "in")
plot(t, black_arm_ave["out"], linewidth=2.0, color=c_out, label = "out")
plot(t, black_arm_ave["tol"],  linewidth=2.0, color=c_tol, label = "total")
plot(t, black_arm_ave["in"]+black_arm_ave["in_std"], linewidth=1.0, linestyle = "--", color=c_in, alpha=0.7)
plot(t, black_arm_ave["in"]-black_arm_ave["in_std"], linewidth=1.0, linestyle = "--", color=c_in, alpha=0.7)
plot(t, black_arm_ave["out"]+black_arm_ave["out_std"], linewidth=1.0, linestyle = "--", color=c_out, alpha=0.7)
plot(t, black_arm_ave["out"]-black_arm_ave["out_std"], linewidth=1.0, linestyle = "--", color=c_out, alpha=0.7)
plot(t, black_arm_ave["tol"]+black_arm_ave["tol_std"], linewidth=1.0, linestyle = "--", color=c_tol, alpha=0.7)
plot(t, black_arm_ave["tol"]-black_arm_ave["tol_std"], linewidth=1.0, linestyle = "--", color=c_tol, alpha=0.7)
# xlabel('Time (ps)')
ylabel(r'$\kappa$(Wm$^{-1}$K$^{-1}$)')       
xlim([0, 2])
gca().set_xticks(linspace(0, 2, 5))
ylim([-2, 20])
gca().set_yticks(linspace(0, 20, 5))
text(0.45, 14, "Black-P, amrchair, total", color=c_tol, fontsize = 18)
text(0.45, 9, "Black-P, amrchair, in", color=c_in, fontsize = 18)
text(0.45, 1, "Black-P, amrchair, out", color=c_out, fontsize = 18)
text(0.05, 16, "(a)", fontsize = 20)

subplot(3, 2, 2)
set_fig_properties([gca()])
plot(t, black_zig_ave["in"], linewidth=2.0, color=c_in, label = "in")
plot(t, black_zig_ave["out"], linewidth=2.0, color=c_out, label = "out")
plot(t, black_zig_ave["tol"],  linewidth=2.0, color=c_tol, label = "total")
plot(t, black_zig_ave["in"]+black_zig_ave["in_std"], linewidth=1.0, linestyle = "--", color=c_in, alpha=0.7)
plot(t, black_zig_ave["in"]-black_zig_ave["in_std"], linewidth=1.0, linestyle = "--", color=c_in, alpha=0.7)
plot(t, black_zig_ave["out"]+black_zig_ave["out_std"], linewidth=1.0, linestyle = "--", color=c_out, alpha=0.7)
plot(t, black_zig_ave["out"]-black_zig_ave["out_std"], linewidth=1.0, linestyle = "--", color=c_out, alpha=0.7)
plot(t, black_zig_ave["tol"]+black_zig_ave["tol_std"], linewidth=1.0, linestyle = "--", color=c_tol, alpha=0.7)
plot(t, black_zig_ave["tol"]-black_zig_ave["tol_std"], linewidth=1.0, linestyle = "--", color=c_tol, alpha=0.7)
# xlabel('Time (ps)')
# ylabel(r'$\kappa$(Wm$^{-1}$K$^{-1}$)')       
xlim([0, 2])
gca().set_xticks(linspace(0, 2, 5))
ylim([-10, 120])
gca().set_yticks(linspace(0, 120, 5))
text(0.65, 85, "Black-P, zigzag, total", color=c_tol, fontsize = 18)
text(0.65, 40, "Black-P, zigzag, in", color=c_in, fontsize = 18)
text(0.65, 25, "Black-P, zigzag, out", color=c_out, fontsize = 18)
text(0.05, 96, "(b)", fontsize = 20)


t1, blue_arm, blue_arm_ave = output_kappa("HNEMD_diffusive/Blue/300K/x_1e-5/", "x", 5000)
t2, blue_zig, blue_zig_ave = output_kappa("HNEMD_diffusive/Blue/300k/y_1e-5/", "y", 5000)
blue = []
for keys in blue_arm:
    blue.append(blue_arm[keys]["kx_ra"][-1])
for keys in blue_zig:
    blue.append(blue_zig[keys]["ky_ra"][-1])
subplot(3, 2, 3)
set_fig_properties([gca()])
plot(t1, blue_arm_ave["in"], linewidth=2.0, color=c_in, label = "in")
plot(t1, blue_arm_ave["out"], linewidth=2.0, color=c_out, label = "out")
plot(t1, blue_arm_ave["tol"],  linewidth=2.0, color=c_tol, label = "total")
plot(t1, blue_arm_ave["in"]+blue_arm_ave["in_std"], linewidth=1.0, linestyle = "--", color=c_in, alpha=0.7)
plot(t1, blue_arm_ave["in"]-blue_arm_ave["in_std"], linewidth=1.0, linestyle = "--", color=c_in, alpha=0.7)
plot(t1, blue_arm_ave["out"]+blue_arm_ave["out_std"], linewidth=1.0, linestyle = "--", color=c_out, alpha=0.7)
plot(t1, blue_arm_ave["out"]-blue_arm_ave["out_std"], linewidth=1.0, linestyle = "--", color=c_out, alpha=0.7)
plot(t1, blue_arm_ave["tol"]+blue_arm_ave["tol_std"], linewidth=1.0, linestyle = "--", color=c_tol, alpha=0.7)
plot(t1, blue_arm_ave["tol"]-blue_arm_ave["tol_std"], linewidth=1.0, linestyle = "--", color=c_tol, alpha=0.7)
# xlabel('Time (ps)')
ylabel(r'$\kappa$(Wm$^{-1}$K$^{-1}$)')       
xlim([0, 5])
gca().set_xticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
ylim([-2, 200])
gca().set_yticks(linspace(0,200,5))
text(1.25, 140, "Blue-P, armchair, total", color=c_tol, fontsize = 18)
text(1.25, 90, "Blue-P, armchair, in", color=c_in, fontsize = 18)
text(1.25, 50, "Blue-P, armchair, out", color=c_out, fontsize = 18)
text(0.125, 160, "(c)", fontsize=20)



subplot(3, 2, 4)
set_fig_properties([gca()])
plot(t2, blue_zig_ave["in"], linewidth=2.0, color=c_in, label = "in")
plot(t2, blue_zig_ave["out"], linewidth=2.0, color=c_out, label = "out")
plot(t2, blue_zig_ave["tol"],  linewidth=2.0, color=c_tol, label = "total")
plot(t2, blue_zig_ave["in"]+blue_zig_ave["in_std"], linewidth=1.0, linestyle = "--", color=c_in, alpha=0.7)
plot(t2, blue_zig_ave["in"]-blue_zig_ave["in_std"], linewidth=1.0, linestyle = "--", color=c_in, alpha=0.7)
plot(t2, blue_zig_ave["out"]+blue_zig_ave["out_std"], linewidth=1.0, linestyle = "--", color=c_out, alpha=0.7)
plot(t2, blue_zig_ave["out"]-blue_zig_ave["out_std"], linewidth=1.0, linestyle = "--", color=c_out, alpha=0.7)
plot(t2, blue_zig_ave["tol"]+blue_zig_ave["tol_std"], linewidth=1.0, linestyle = "--", color=c_tol, alpha=0.7)
plot(t2, blue_zig_ave["tol"]-blue_zig_ave["tol_std"], linewidth=1.0, linestyle = "--", color=c_tol, alpha=0.7)
# xlabel('Time (ps)')
# ylabel(r'$\kappa$(Wm$^{-1}$K$^{-1}$)')        
xlim([0, 5])
gca().set_xticks(linspace(0, 5, 6))
ylim([-2, 200])
gca().set_yticks(linspace(0,200,5))
text(1.25, 145, "Blue-P, zigzag, total", color=c_tol, fontsize = 18)
text(1.25, 95, "Blue-P, zigzag, in", color=c_in, fontsize = 18)
text(1.25, 55, "Blue-P, zigzag, out", color=c_out, fontsize = 18)
text(0.125, 160, "(d)", fontsize = 20)


t, violet_arm, violet_arm_ave = output_kappa("HNEMD_diffusive/Violet/300K/x_1e-4/", "x", 2000)
t, violet_zig, violet_zig_ave = output_kappa("HNEMD_diffusive/Violet/300k/y_1e-4/", "y", 2000)
subplot(3, 2, 5)
set_fig_properties([gca()])
plot(t, violet_arm_ave["in"], linewidth=2.0, color=c_in, label = "in")
plot(t, violet_arm_ave["out"], linewidth=2.0, color=c_out, label = "out")
plot(t, violet_arm_ave["tol"],  linewidth=2.0, color=c_tol, label = "total")
plot(t, violet_arm_ave["in"]+violet_arm_ave["in_std"], linewidth=1.0, linestyle = "--", color=c_in, alpha=0.7)
plot(t, violet_arm_ave["in"]-violet_arm_ave["in_std"], linewidth=1.0, linestyle = "--", color=c_in, alpha=0.7)
plot(t, violet_arm_ave["out"]+violet_arm_ave["out_std"], linewidth=1.0, linestyle = "--", color=c_out, alpha=0.7)
plot(t, violet_arm_ave["out"]-violet_arm_ave["out_std"], linewidth=1.0, linestyle = "--", color=c_out, alpha=0.7)
plot(t, violet_arm_ave["tol"]+violet_arm_ave["tol_std"], linewidth=1.0, linestyle = "--", color=c_tol, alpha=0.7)
plot(t, violet_arm_ave["tol"]-violet_arm_ave["tol_std"], linewidth=1.0, linestyle = "--", color=c_tol, alpha=0.7)
xlabel('Time (ns)')
ylabel(r'$\kappa$(Wm$^{-1}$K$^{-1}$)')       
xlim([0, 2])
gca().set_xticks(linspace(0, 2, 5))
ylim([-0.1, 3])
gca().set_yticks(linspace(0, 3, 4))
text(0.55, 2.5, "Violet-P, x, total", color=c_tol, fontsize = 18)
text(0.55, 1.3, "Violet-P, x, in", color=c_in, fontsize = 18)
text(0.55, 0.5, "Violet-P, x, out", color=c_out, fontsize = 18)
text(0.05, 2.4, "(e)", fontsize = 20)

subplot(3,2,6)
set_fig_properties([gca()])
plot(t, violet_zig_ave["in"], linewidth=2.0, color=c_in, label = "in")
plot(t, violet_zig_ave["out"], linewidth=2.0, color=c_out, label = "out")
plot(t, violet_zig_ave["tol"],  linewidth=2.0, color=c_tol, label = "total")
plot(t, violet_zig_ave["in"]+violet_zig_ave["in_std"], linewidth=1.0, linestyle = "--", color=c_in, alpha=0.7)
plot(t, violet_zig_ave["in"]-violet_zig_ave["in_std"], linewidth=1.0, linestyle = "--", color=c_in, alpha=0.7)
plot(t, violet_zig_ave["out"]+violet_zig_ave["out_std"], linewidth=1.0, linestyle = "--", color=c_out, alpha=0.7)
plot(t, violet_zig_ave["out"]-violet_zig_ave["out_std"], linewidth=1.0, linestyle = "--", color=c_out, alpha=0.7)
plot(t, violet_zig_ave["tol"]+violet_zig_ave["tol_std"], linewidth=1.0, linestyle = "--", color=c_tol, alpha=0.7)
plot(t, violet_zig_ave["tol"]-violet_zig_ave["tol_std"], linewidth=1.0, linestyle = "--", color=c_tol, alpha=0.7)
xlabel('Time (ns)')
# ylabel(r'$\kappa$(Wm$^{-1}$K$^{-1}$)')       
xlim([0, 2])
gca().set_xticks(linspace(0, 2, 5))
ylim([-0.1, 3])
gca().set_yticks(linspace(0, 3, 4))
text(0.65, 2.5, "Violet-P, y, total", color=c_tol, fontsize = 18)
text(0.65, 1.3, "Violet-P, y, in", color=c_in, fontsize = 18)
text(0.65, 0.5, "Violet-P, y, out", color=c_out, fontsize = 18)
text(0.05, 2.4, "(f)", fontsize = 20)
subplots_adjust(wspace=0.3, hspace=0.3)
savefig("./fig5-hnemd.pdf", bbox_inches='tight')
