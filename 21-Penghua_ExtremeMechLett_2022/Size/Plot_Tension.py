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
        A_x, B_x = optimize.curve_fit(f_1, output[keys]["epsilon_x"][:10], output[keys]["sigma_x"][:10])[0]
        A_y, B_y = optimize.curve_fit(f_1, output[keys]["epsilon_y"][:10], output[keys]["sigma_y"][:10])[0]
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

x_5nm = tension_calc("5nm/x")  
y_5nm = tension_calc("5nm/y") 
x_11nm = tension_calc("11nm/x")  
y_11nm = tension_calc("11nm/y") 
x_22nm = tension_calc("22nm/x") 
y_22nm = tension_calc("22nm/y") 
x_33nm = tension_calc("33nm/x") 
y_33nm = tension_calc("33nm/y")
x_44nm = tension_calc("44nm/x") 
y_44nm = tension_calc("44nm/y")
clist = plt.get_cmap("rainbow")(np.linspace(0, 1, 5))

figure(figsize=(15, 10))
subplot(2, 1, 1)
set_fig_properties([gca()])
for keys in x_5nm.keys():
    if keys == "run0":
        plot(x_5nm[keys]["epsilon_x"] * 100, x_5nm[keys]["sigma_x"], ls = '-', color = clist[0], lw = lw, label=r"5 nm, $x$")
        plot(x_11nm[keys]["epsilon_x"] * 100, x_11nm[keys]["sigma_x"], ls = '-', color = clist[1], lw = lw, label=r"11 nm, $x$")
        plot(x_22nm[keys]["epsilon_x"] * 100, x_22nm[keys]["sigma_x"], ls = '-', color = clist[2], lw = lw, label=r"22 nm, $x$")
        plot(x_33nm[keys]["epsilon_x"] * 100, x_33nm[keys]["sigma_x"], ls = '-', color = clist[3], lw = lw, label=r"33 nm, $x$")
        plot(x_44nm[keys]["epsilon_x"] * 100, x_44nm[keys]["sigma_x"], ls = '-', color = clist[4], lw = lw, label=r"44 nm, $x$")
    elif "run" in keys:
        plot(x_5nm[keys]["epsilon_x"] * 100, x_5nm[keys]["sigma_x"], ls = '-', color = clist[0], lw = lw)
        plot(x_11nm[keys]["epsilon_x"] * 100, x_11nm[keys]["sigma_x"], ls = '-', color = clist[1], lw = lw)
        plot(x_22nm[keys]["epsilon_x"] * 100, x_22nm[keys]["sigma_x"], ls = '-', color = clist[2], lw = lw)
        plot(x_33nm[keys]["epsilon_x"] * 100, x_33nm[keys]["sigma_x"], ls = '-', color = clist[3], lw = lw)
        plot(x_44nm[keys]["epsilon_x"] * 100, x_44nm[keys]["sigma_x"], ls = '-', color = clist[4], lw = lw)   

for keys in y_5nm.keys():
    if keys == "run0":
        plot(y_5nm[keys]["epsilon_y"] * 100, y_5nm[keys]["sigma_y"], ls = '--', color = clist[0], lw = lw, label=r"5 nm, $y$")
        plot(y_11nm[keys]["epsilon_y"] * 100, y_11nm[keys]["sigma_y"], ls = '--', color = clist[1], lw = lw, label=r"11 nm, $y$")
        plot(y_22nm[keys]["epsilon_y"] * 100, y_22nm[keys]["sigma_y"], ls = '--', color = clist[2], lw = lw, label=r"22 nm, $y$")
        plot(y_33nm[keys]["epsilon_y"] * 100, y_33nm[keys]["sigma_y"], ls = '--', color = clist[3], lw = lw, label=r"33 nm, $y$")
        plot(y_44nm[keys]["epsilon_y"] * 100, y_44nm[keys]["sigma_y"], ls = '--', color = clist[4], lw = lw, label=r"44 nm, $y$")
    elif "run" in keys:
        plot(y_5nm[keys]["epsilon_y"] * 100, y_5nm[keys]["sigma_y"], ls = '--', color = clist[0], lw = lw)
        plot(y_11nm[keys]["epsilon_y"] * 100, y_11nm[keys]["sigma_y"], ls = '--', color = clist[1], lw = lw)
        plot(y_22nm[keys]["epsilon_y"] * 100, y_22nm[keys]["sigma_y"], ls = '--', color = clist[2], lw = lw)
        plot(y_33nm[keys]["epsilon_y"] * 100, y_33nm[keys]["sigma_y"], ls = '--', color = clist[3], lw = lw)
        plot(y_44nm[keys]["epsilon_y"] * 100, y_44nm[keys]["sigma_y"], ls = '--', color = clist[4], lw = lw)  

# xlabel(r"$\varepsilon$ (%)")
# ylabel(r"$\sigma$ (GPa)")    
xlabel(r"Strain (%)")
ylabel(r"Stress (GPa)")    
xlim([0, 5])
gca().set_xticks(linspace(0, 5, 6))
ylim([0, 10])
gca().set_yticks(linspace(0, 10, 6))
legend(loc = "upper left", ncol = 2, labelspacing = 0.2, frameon = True)
title("(a)", fontsize=18)

subplot(2, 3, 4)
set_fig_properties([gca()])
sizes = [5, 11, 22, 33, 44]
modulus_output = np.zeros((5, 5))
for i in range(len(sizes)):
    modulus_output[i, 0] = sizes[i]
    modulus_output[i, 1] = tension_calc(str(sizes[i]) + "nm/x")["modulus_x_ave"]
    modulus_output[i, 2] = tension_calc(str(sizes[i]) + "nm/x")["modulus_x_bar"]
    modulus_output[i, 3] = tension_calc(str(sizes[i]) + "nm/y")["modulus_y_ave"]
    modulus_output[i, 4] = tension_calc(str(sizes[i]) + "nm/y")["modulus_y_bar"]    

errorbar(modulus_output[:, 0], modulus_output[:, 1], modulus_output[:, 2], 
          fmt = "bo",
          markersize = 7.5,
          label=r"$x$")

errorbar(modulus_output[:, 0], modulus_output[:, 3], modulus_output[:, 4], 
          fmt = "rs",
          markersize = 7.5,
          label=r"$y$")

# xlabel(r"$a$ (nm)")
# ylabel(r'$E$ (GPa)')    

xlabel(r"Membrane length (nm)")
ylabel(r"Tensile modulus (GPa)") 
  
xlim([0, 50])
gca().set_xticks(linspace(0, 50, 6))
ylim([140, 240])
gca().set_yticks(linspace(140, 240, 6))
legend(loc = "center left")
title("(b)", fontsize=16)

subplot(2, 3, 5)
set_fig_properties([gca()])
sizes = [5, 11, 22, 33, 44]
sigma_max_output = np.zeros((5, 5))
for i in range(len(sizes)):
    sigma_max_output[i, 0] = sizes[i]
    sigma_max_output[i, 1] = tension_calc(str(sizes[i]) + "nm/x")["sigma_x_max_ave"]
    sigma_max_output[i, 2] = tension_calc(str(sizes[i]) + "nm/x")["sigma_x_max_bar"]
    sigma_max_output[i, 3] = tension_calc(str(sizes[i]) + "nm/y")["sigma_y_max_ave"]
    sigma_max_output[i, 4] = tension_calc(str(sizes[i]) + "nm/y")["sigma_y_max_bar"]    

errorbar(sigma_max_output[:, 0], sigma_max_output[:, 1], sigma_max_output[:, 2], 
          fmt = "bo",
          markersize = 7.5,
          label=r"$x$")

errorbar(sigma_max_output[:, 0], sigma_max_output[:, 3], sigma_max_output[:, 4], 
          fmt = "rs",
          markersize = 7.5,
          label=r"$y$")

xlabel(r"Membrane length (nm)")
ylabel(r"Tensile strength (GPa)")     
xlim([0, 50])
gca().set_xticks(linspace(0, 50, 6))
ylim([3, 10])
gca().set_yticks(linspace(3, 10, 8)) 
legend(loc = "center left")
title("(c)", fontsize=16)

subplot(2, 3, 6)
set_fig_properties([gca()])
sizes = [5, 11, 22, 33, 44]
epsilon_max_output = np.zeros((5, 5))
for i in range(len(sizes)):
    epsilon_max_output[i, 0] = sizes[i]
    epsilon_max_output[i, 1] = tension_calc(str(sizes[i]) + "nm/x")["epsilon_x_max_ave"]
    epsilon_max_output[i, 2] = tension_calc(str(sizes[i]) + "nm/x")["epsilon_x_max_bar"]
    epsilon_max_output[i, 3] = tension_calc(str(sizes[i]) + "nm/y")["epsilon_y_max_ave"]
    epsilon_max_output[i, 4] = tension_calc(str(sizes[i]) + "nm/y")["epsilon_y_max_bar"]    

errorbar(epsilon_max_output[:, 0], epsilon_max_output[:, 1] * 100, epsilon_max_output[:, 2] * 100, 
          fmt = "bo",
          markersize = 7.5,
          label=r"$x$")

errorbar(epsilon_max_output[:, 0], epsilon_max_output[:, 3] * 100, epsilon_max_output[:, 4] * 100, 
          fmt = "rs",
          markersize = 7.5,
          label=r"$y$")

xlabel(r"Membrane length (nm)")
ylabel(r"Fracture strain (%)")    
xlim([0, 50])
gca().set_xticks(linspace(0, 50, 6))
ylim([2, 4.5])
gca().set_yticks(linspace(2, 4.5, 6))
legend(loc = "center left")
title("(d)", fontsize=16)

subplots_adjust(wspace = 0.35, hspace = 0.35)
savefig("Tension_size.pdf", bbox_inches='tight')
  


# modulus_x = np.array([x_5nm["modulus_x_ave"], x_11nm["modulus_x_ave"], x_22nm["modulus_x_ave"], x_33nm["modulus_x_ave"], x_44nm["modulus_x_ave"]])
# modulus_y = np.array([y_5nm["modulus_y_ave"], y_11nm["modulus_y_ave"], y_22nm["modulus_y_ave"], y_33nm["modulus_y_ave"], y_44nm["modulus_y_ave"]])
        
        
    
    
    
    
    
    
    
    
    
    
    
    