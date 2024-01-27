# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 16:15 2022

@author: hityingph
"""

from thermo.gpumd.data import load_thermo
from pylab import *
from scipy import optimize

aw = 1.5
fs = 16
font = {'size'   : fs}
matplotlib.rc('font', **font)
matplotlib.rc('axes' , linewidth=aw)

def f_1(x, A, B):
    return A * x + B

def avg_10(data):
    length = int(len(data)/50)*50
    output = []
    for i in range(int(len(data)/50)):
        output.append(np.average(data[i * 50 : i * 50 + 50]))
    output = np.array(output)
    return output

def tension_calc(path):
    thermo = load_thermo(path)
    output = dict()
    for i in range(int(len(thermo["T"]) / 5000)):
        temp = dict()
        for keys in thermo.keys():
            temp[keys] = avg_10(thermo[keys][5000 * i : 5000 * (i + 1)])
            output["run%s"%i] = temp
    
    sigma_x_max = []
    sigma_y_max = []
    epsilon_x_max = []
    epsilon_y_max = []
    modulus_x = []
    modulus_y = []
    for keys in output.keys():
        output[keys]["sigma_x"] = output[keys]["Px"] * -1
        output[keys]["sigma_y"] = output[keys]["Py"] * -1
        output[keys]["epsilon_x"] = output[keys]["Lx"]/output[keys]["Lx"][0] - 1
        output[keys]["epsilon_y"] = output[keys]["Ly"]/output[keys]["Ly"][0] - 1
        output[keys]["poisson_xy"] = -output[keys]["epsilon_y"][1:] / output[keys]["epsilon_x"][1:]
        output[keys]["poisson_yx"] = -output[keys]["epsilon_x"][1:] / output[keys]["epsilon_y"][1:]
        sigma_x_max.append(max(output[keys]["sigma_x"]))
        sigma_y_max.append(max(output[keys]["sigma_y"]))
        epsilon_x_max.append(output[keys]["epsilon_x"][np.argmax(output[keys]["sigma_x"])])
        epsilon_y_max.append(output[keys]["epsilon_y"][np.argmax(output[keys]["sigma_y"])])
        A_x, B_x = optimize.curve_fit(f_1, output[keys]["epsilon_x"][:20], output[keys]["sigma_x"][:20])[0]
        A_y, B_y = optimize.curve_fit(f_1, output[keys]["epsilon_y"][:20], output[keys]["sigma_y"][:20])[0]
        modulus_x.append(A_x)
        modulus_y.append(A_y)
    
    output["modulus_x"] = np.array(modulus_x)
    output["modulus_y"] = np.array(modulus_y)
    output["sigma_x_max"] = np.array(sigma_x_max)
    output["sigma_y_max"] = np.array(sigma_y_max)
    output["epsilon_x_max"] = np.array(epsilon_x_max)
    output["epsilon_y_max"] = np.array(epsilon_y_max)
    output["sigma_x_max_ave"] = np.average(np.array(sigma_x_max))
    output["sigma_y_max_ave"] = np.average(np.array(sigma_y_max))
    output["epsilon_x_max_ave"] = np.average(np.array(epsilon_x_max))
    output["epsilon_y_max_ave"] = np.average(np.array(epsilon_y_max))
    output["modulus_x_ave"] = np.average(np.array(modulus_x))
    output["modulus_y_ave"] = np.average(np.array(modulus_y))
    output["sigma_x_max_bar"] = np.std(np.array(sigma_x_max)) / np.sqrt(len(np.array(sigma_x_max)))
    output["sigma_y_max_bar"] = np.std(np.array(sigma_y_max)) / np.sqrt(len(np.array(sigma_y_max)))
    output["epsilon_x_max_bar"] = np.std(np.array(epsilon_x_max)) / np.sqrt(len(np.array(epsilon_x_max)))
    output["epsilon_y_max_bar"] = np.std(np.array(epsilon_y_max)) / np.sqrt(len(np.array(epsilon_y_max)))
    output["modulus_x_bar"] = np.std(np.array(modulus_x)) / np.sqrt(len(np.array(modulus_x)))
    output["modulus_y_bar"] = np.std(np.array(modulus_y)) / np.sqrt(len(np.array(modulus_y)))
    return output

lw = 1.5
def set_fig_properties(ax_list):
    tl = 6
    tw = 1.5
    tlm = 3
    
    for ax in ax_list:
        ax.tick_params(which='major', length=tl, width=tw)
        ax.tick_params(which='minor', length=tlm, width=tw)
        ax.tick_params(which='both', axis='both', direction='out', right=False, top=False)

x_100K = tension_calc("100K/x")  
y_100K = tension_calc("100K/y") 
x_200K = tension_calc("200K/x")  
y_200K = tension_calc("200K/y") 
x_300K = tension_calc("300K/x") 
y_300K = tension_calc("300K/y") 
x_400K = tension_calc("400K/x") 
y_400K = tension_calc("400K/y") 
x_500K = tension_calc("500K/x") 
y_500K = tension_calc("500K/y")
x_600K = tension_calc("600K/x") 
y_600K = tension_calc("600K/y")
x_700K = tension_calc("700K/x") 
y_700K = tension_calc("700K/y")
x_800K = tension_calc("800K/x") 
y_800K = tension_calc("800K/y")
clist = plt.get_cmap("rainbow")(np.linspace(0, 1, 8))
 
figure(figsize=(15, 10))
subplot(2, 1, 1)
set_fig_properties([gca()])
for keys in x_800K.keys():
    if keys == "run0":
        plot(x_100K[keys]["epsilon_x"] * 100, x_100K[keys]["sigma_x"], ls = '-', color = clist[0], lw = lw, label=r"100 K, $x$")
        plot(x_200K[keys]["epsilon_x"] * 100, x_200K[keys]["sigma_x"], ls = '-', color = clist[1], lw = lw, label=r"200 K, $x$")
        plot(x_300K[keys]["epsilon_x"] * 100, x_300K[keys]["sigma_x"], ls = '-', color = clist[2], lw = lw, label=r"300 K, $x$")
        plot(x_400K[keys]["epsilon_x"] * 100, x_400K[keys]["sigma_x"], ls = '-', color = clist[3], lw = lw, label=r"400 K, $x$")
        plot(x_500K[keys]["epsilon_x"] * 100, x_500K[keys]["sigma_x"], ls = '-', color = clist[4], lw = lw, label=r"500 K, $x$")
        plot(x_600K[keys]["epsilon_x"] * 100, x_600K[keys]["sigma_x"], ls = '-', color = clist[5], lw = lw, label=r"600 K, $x$")
        plot(x_700K[keys]["epsilon_x"] * 100, x_700K[keys]["sigma_x"], ls = '-', color = clist[6], lw = lw, label=r"700 K, $x$")
        plot(x_800K[keys]["epsilon_x"] * 100, x_800K[keys]["sigma_x"], ls = '-', color = clist[7], lw = lw, label=r"800 K, $x$")

    elif "run" in keys:
        plot(x_100K[keys]["epsilon_x"] * 100, x_100K[keys]["sigma_x"], ls = '-', color = clist[0], lw = lw)
        plot(x_200K[keys]["epsilon_x"] * 100, x_200K[keys]["sigma_x"], ls = '-', color = clist[1], lw = lw)
        plot(x_300K[keys]["epsilon_x"] * 100, x_300K[keys]["sigma_x"], ls = '-', color = clist[2], lw = lw)
        plot(x_400K[keys]["epsilon_x"] * 100, x_400K[keys]["sigma_x"], ls = '-', color = clist[3], lw = lw)
        plot(x_500K[keys]["epsilon_x"] * 100, x_500K[keys]["sigma_x"], ls = '-', color = clist[4], lw = lw)
        plot(x_600K[keys]["epsilon_x"] * 100, x_600K[keys]["sigma_x"], ls = '-', color = clist[5], lw = lw)
        plot(x_700K[keys]["epsilon_x"] * 100, x_700K[keys]["sigma_x"], ls = '-', color = clist[6], lw = lw)
        plot(x_800K[keys]["epsilon_x"] * 100, x_800K[keys]["sigma_x"], ls = '-', color = clist[7], lw = lw)

for keys in y_800K.keys():
    if keys == "run0":
        plot(y_100K[keys]["epsilon_y"] * 100, y_100K[keys]["sigma_y"], ls = '--', color = clist[0], lw = lw, label=r"100 K, $y$")
        plot(y_200K[keys]["epsilon_y"] * 100, y_200K[keys]["sigma_y"], ls = '--', color = clist[1], lw = lw, label=r"200 K, $y$")        
        plot(y_300K[keys]["epsilon_y"] * 100, y_300K[keys]["sigma_y"], ls = '--', color = clist[2], lw = lw, label=r"300 K, $y$")
        plot(y_400K[keys]["epsilon_y"] * 100, y_400K[keys]["sigma_y"], ls = '--', color = clist[3], lw = lw, label=r"400 K, $y$")
        plot(y_500K[keys]["epsilon_y"] * 100, y_500K[keys]["sigma_y"], ls = '--', color = clist[4], lw = lw, label=r"500 K, $y$")        
        plot(y_600K[keys]["epsilon_y"] * 100, y_600K[keys]["sigma_y"], ls = '--', color = clist[5], lw = lw, label=r"600 K, $y$")
        plot(y_700K[keys]["epsilon_y"] * 100, y_700K[keys]["sigma_y"], ls = '--', color = clist[6], lw = lw, label=r"700 K, $y$")
        plot(y_800K[keys]["epsilon_y"] * 100, y_800K[keys]["sigma_y"], ls = '--', color = clist[7], lw = lw, label=r"800 K, $y$")        



    elif "run" in keys:
        plot(y_100K[keys]["epsilon_y"] * 100, y_100K[keys]["sigma_y"], ls = '--', color = clist[0], lw = lw)
        plot(y_200K[keys]["epsilon_y"] * 100, y_200K[keys]["sigma_y"], ls = '--', color = clist[1], lw = lw)        
        plot(y_300K[keys]["epsilon_y"] * 100, y_300K[keys]["sigma_y"], ls = '--', color = clist[2], lw = lw)
        plot(y_400K[keys]["epsilon_y"] * 100, y_400K[keys]["sigma_y"], ls = '--', color = clist[3], lw = lw)
        plot(y_500K[keys]["epsilon_y"] * 100, y_500K[keys]["sigma_y"], ls = '--', color = clist[4], lw = lw)        
        plot(y_600K[keys]["epsilon_y"] * 100, y_600K[keys]["sigma_y"], ls = '--', color = clist[5], lw = lw)
        plot(y_700K[keys]["epsilon_y"] * 100, y_700K[keys]["sigma_y"], ls = '--', color = clist[6], lw = lw)
        plot(y_800K[keys]["epsilon_y"] * 100, y_800K[keys]["sigma_y"], ls = '--', color = clist[7], lw = lw)        

xlabel(r"Strain (%)")
ylabel(r"Stress (GPa)")   
xlim([0, 6])
gca().set_xticks(linspace(0, 6, 7))
ylim([0, 18])
gca().set_yticks(linspace(0, 18, 7))
legend(loc = "upper left", ncol = 4, fontsize = 16, columnspacing =1, frameon = True)
title("(a)", fontsize=16)

subplot(2, 3, 4)
set_fig_properties([gca()])
temps = np.arange(100, 900, 100)
modulus_output = np.zeros((8, 5))
for i in range(len(temps)):
    modulus_output[i, 0] = temps[i]
    modulus_output[i, 1] = tension_calc(str(temps[i]) + "K/x")["modulus_x_ave"]
    modulus_output[i, 2] = tension_calc(str(temps[i]) + "K/x")["modulus_x_bar"]
    modulus_output[i, 3] = tension_calc(str(temps[i]) + "K/y")["modulus_y_ave"]
    modulus_output[i, 4] = tension_calc(str(temps[i]) + "K/y")["modulus_y_bar"]    

errorbar(modulus_output[:, 0], modulus_output[:, 1], modulus_output[:, 2], 
          fmt = "bo",
          markersize = 7.5,
          label=r"$x$")

errorbar(modulus_output[:, 0], modulus_output[:, 3], modulus_output[:, 4], 
          fmt = "rs",
          markersize = 7.5,
          label=r"$y$")
# hlines(162.8, 0, 1000, colors = "blue", ls = "-", label = r"DFT, 0K, $x$")
# hlines(207, 0, 1000, colors = "red", ls = "--", label = r"DFT, 0K, $y$")

xlabel(r"Temperature (K)")
ylabel(r"Tensile modulus (GPa)")    
xlim([50, 850])
gca().set_xticks(linspace(200, 800, 4))
ylim([120, 240])
gca().set_yticks(linspace(120, 240, 7))
legend(loc = "lower left")
title("(b)", fontsize=16)

subplot(2, 3, 5)
set_fig_properties([gca()])
temps = np.arange(100, 900, 100)
sigma_output = np.zeros((8, 5))
for i in range(len(temps)):
    sigma_output[i, 0] = temps[i]
    sigma_output[i, 1] = tension_calc(str(temps[i]) + "K/x")["sigma_x_max_ave"]
    sigma_output[i, 2] = tension_calc(str(temps[i]) + "K/x")["sigma_x_max_bar"]
    sigma_output[i, 3] = tension_calc(str(temps[i]) + "K/y")["sigma_y_max_ave"]
    sigma_output[i, 4] = tension_calc(str(temps[i]) + "K/y")["sigma_y_max_bar"]    

errorbar(sigma_output[:, 0], sigma_output[:, 1], sigma_output[:, 2], 
          fmt = "bo",
          markersize = 7.5,
          label=r"$x$")

errorbar(sigma_output[:, 0], sigma_output[:, 3], sigma_output[:, 4], 
          fmt = "rs",
          markersize = 7.5,
          label=r"$y$")

xlabel(r"Temperature (K)")
ylabel(r"Tensile strength (GPa)")     
xlim([50, 850])
gca().set_xticks(linspace(200, 800, 4))
ylim([0, 13])
gca().set_yticks(linspace(0, 12, 7))
legend(loc = "lower left")
title("(c)", fontsize=16)

subplot(2, 3, 6)
set_fig_properties([gca()])
temps = np.arange(100, 900, 100)
epsilon_output = np.zeros((8, 5))
for i in range(len(temps)):
    epsilon_output[i, 0] = temps[i]
    epsilon_output[i, 1] = tension_calc(str(temps[i]) + "K/x")["epsilon_x_max_ave"]
    epsilon_output[i, 2] = tension_calc(str(temps[i]) + "K/x")["epsilon_x_max_bar"]
    epsilon_output[i, 3] = tension_calc(str(temps[i]) + "K/y")["epsilon_y_max_ave"]
    epsilon_output[i, 4] = tension_calc(str(temps[i]) + "K/y")["epsilon_y_max_bar"]    

errorbar(epsilon_output[:, 0], epsilon_output[:, 1] * 100, epsilon_output[:, 2] * 100, 
          fmt = "bo",
          markersize = 7.5,
          label=r"$x$")
errorbar(epsilon_output[:, 0], epsilon_output[:, 3] * 100, epsilon_output[:, 4] * 100, 
          fmt = "rs",
          markersize = 7.5,
          label=r"$y$")

xlabel(r"Temperature (K)")
ylabel(r"Fracture strain (%)")    
xlim([50, 850])
gca().set_xticks(linspace(200, 800, 4))
ylim([0, 6])
gca().set_yticks(linspace(0, 6, 7))
legend(loc = "lower left")
title("(d)", fontsize=16)

subplots_adjust(wspace = 0.35, hspace = 0.35)
savefig("Tension_temp.pdf", bbox_inches='tight')
    
    
    
    
    
    
    
    
    