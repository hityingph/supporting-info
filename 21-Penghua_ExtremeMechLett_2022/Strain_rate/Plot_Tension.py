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

def f_2(x, A, B):
    return A * x ** B

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

def tension_quasi(path):
    thermo = load_thermo(path)
    stress_strain = np.zeros((int(len(thermo["T"])/200), 4))
    for i in range(int(len(thermo["T"])/200)):
        stress_strain[i, 0] = np.mean(thermo["Lx"][200 * i + 100 : 200 * i + 200])
        stress_strain[i, 1] = -np.mean(thermo["Px"][200 * i + 100 : 200 * i + 200])
        stress_strain[i, 2] = np.mean(thermo["Ly"][200 * i + 100 : 200 * i + 200])
        stress_strain[i, 3] = -np.mean(thermo["Py"][200 * i + 100 : 200 * i + 200])
    stress_strain[:, 0] = stress_strain[:, 0] / stress_strain[0, 0] - 1
    stress_strain[:, 2] = stress_strain[:, 2] / stress_strain[0, 2] - 1
    stress_strain[:, 1] = stress_strain[:, 1] - stress_strain[0, 1]
    stress_strain[:, 3] = stress_strain[:, 3] - stress_strain[0, 3]
    
    modulus_x, B_x = optimize.curve_fit(f_1, stress_strain[:5, 0], stress_strain[:5, 1])[0]
    modulus_y, B_y = optimize.curve_fit(f_1, stress_strain[:5, 2], stress_strain[:5, 3])[0]
    
    output = dict()
    output["epsilon_x"] = stress_strain[:, 0]
    output["sigma_x"] = stress_strain[:, 1]
    output["epsilon_y"] = stress_strain[:, 2]
    output["sigma_y"] = stress_strain[:, 3]
    output["modulus_x"] = modulus_x
    output["modulus_y"] = modulus_y
    output["sigma_x_max"] = max(output["sigma_x"])
    output["sigma_y_max"] = max(output["sigma_y"])
    output["epsilon_x_max"] = output["epsilon_x"][np.argmax(output["sigma_x"])]
    output["epsilon_y_max"] = output["epsilon_y"][np.argmax(output["sigma_y"])]
    
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

x_1e7 = tension_calc("1e7/x")  
y_1e7 = tension_calc("1e7/y") 
x_5e7 = tension_calc("5e7/x")  
y_5e7 = tension_calc("5e7/y") 
x_1e8 = tension_calc("1e8/x")  
y_1e8 = tension_calc("1e8/y")  
x_5e8 = tension_calc("5e8/x")  
y_5e8 = tension_calc("5e8/y")  
x_1e9 = tension_calc("1e9/x")  
y_1e9 = tension_calc("1e9/y")  
clist = plt.get_cmap("rainbow")(np.linspace(0, 1, 6))
 
figure(figsize=(15, 10))
subplot(2, 1, 1)
set_fig_properties([gca()])
for keys in x_1e7.keys():
    if keys == "run0":
        plot(x_1e7[keys]["epsilon_x"] * 100, x_1e7[keys]["sigma_x"], ls = '-', color = clist[1], lw = lw, label=r"$1\times 10^{7}$/s, $x$")
        plot(x_5e7[keys]["epsilon_x"] * 100, x_5e7[keys]["sigma_x"], ls = '-', color = clist[2], lw = lw, label=r"$5\times 10^{7}$/s, $x$")
        plot(x_1e8[keys]["epsilon_x"] * 100, x_1e8[keys]["sigma_x"], ls = '-', color = clist[3], lw = lw, label=r"$1\times 10^{8}$/s, $x$")
        plot(x_5e8[keys]["epsilon_x"] * 100, x_5e8[keys]["sigma_x"], ls = '-', color = clist[4], lw = lw, label=r"$5\times 10^{8}$/s, $x$")
        plot(x_1e9[keys]["epsilon_x"] * 100, x_1e9[keys]["sigma_x"], ls = '-', color = clist[5], lw = lw, label=r"$1\times 10^{9}$/s, $x$")

    elif "run" in keys:
        plot(x_1e7[keys]["epsilon_x"] * 100, x_1e7[keys]["sigma_x"], ls = '-', color = clist[1], lw = lw)
        plot(x_5e7[keys]["epsilon_x"] * 100, x_5e7[keys]["sigma_x"], ls = '-', color = clist[2], lw = lw)
        plot(x_1e8[keys]["epsilon_x"] * 100, x_1e8[keys]["sigma_x"], ls = '-', color = clist[3], lw = lw)
        plot(x_5e8[keys]["epsilon_x"] * 100, x_5e8[keys]["sigma_x"], ls = '-', color = clist[4], lw = lw)
        plot(x_1e9[keys]["epsilon_x"] * 100, x_1e9[keys]["sigma_x"], ls = '-', color = clist[5], lw = lw)

for keys in y_1e7.keys():
    if keys == "run0":
        plot(y_1e7[keys]["epsilon_y"] * 100, y_1e7[keys]["sigma_y"], ls = '--', color = clist[1], lw = lw, label=r"$1\times 10^{7}$/s, $y$")
        plot(y_5e7[keys]["epsilon_y"] * 100, y_5e7[keys]["sigma_y"], ls = '--', color = clist[2], lw = lw, label=r"$5\times 10^{7}$/s, $y$")
        plot(y_1e8[keys]["epsilon_y"] * 100, y_1e8[keys]["sigma_y"], ls = '--', color = clist[3], lw = lw, label=r"$1\times 10^{8}$/s, $y$")   
        plot(y_5e8[keys]["epsilon_y"] * 100, y_5e8[keys]["sigma_y"], ls = '--', color = clist[4], lw = lw, label=r"$5\times 10^{8}$/s, $y$")   
        plot(y_1e9[keys]["epsilon_y"] * 100, y_1e9[keys]["sigma_y"], ls = '--', color = clist[5], lw = lw, label=r"$1\times 10^{9}$/s, $y$")  
   



    elif "run" in keys:
        plot(y_1e7[keys]["epsilon_y"] * 100, y_1e7[keys]["sigma_y"], ls = '--', color = clist[1], lw = lw)
        plot(y_5e7[keys]["epsilon_y"] * 100, y_5e7[keys]["sigma_y"], ls = '--', color = clist[2], lw = lw)
        plot(y_1e8[keys]["epsilon_y"] * 100, y_1e8[keys]["sigma_y"], ls = '--', color = clist[3], lw = lw)   
        plot(y_5e8[keys]["epsilon_y"] * 100, y_5e8[keys]["sigma_y"], ls = '--', color = clist[4], lw = lw)   
        plot(y_1e9[keys]["epsilon_y"] * 100, y_1e9[keys]["sigma_y"], ls = '--', color = clist[5], lw = lw)  
    

# xlabel(r"$\varepsilon$ (%)")
# ylabel(r"$\sigma$ (GPa)")       
xlabel(r"Strain (%)")
ylabel(r"Stress (GPa)")     
xlim([0, 5])
gca().set_xticks(linspace(0, 5, 6))
ylim([0, 15])
gca().set_yticks(linspace(0, 15, 6))
legend(loc = "upper left", 
       ncol = 2, 
       fontsize = 16, 
       columnspacing =1, 
       frameon = True)
title("(a)", fontsize=16)

subplot(2, 3, 4)
set_fig_properties([gca()])
rates = ["1e7", "5e7", "1e8", "5e8", "1e9"]
modulus_output = np.zeros((5, 5))
for i in range(len(rates)):
    modulus_output[i, 0] = float(rates[i])/1e7
    modulus_output[i, 1] = tension_calc(rates[i] + "/x")["modulus_x_ave"]
    modulus_output[i, 2] = tension_calc(rates[i] + "/x")["modulus_x_bar"]
    modulus_output[i, 3] = tension_calc(rates[i] + "/y")["modulus_y_ave"]
    modulus_output[i, 4] = tension_calc(rates[i] + "/y")["modulus_y_bar"]    

errorbar(modulus_output[:, 0], modulus_output[:, 1], modulus_output[:, 2], 
          fmt = "bo",
          markersize = 7.5,
          label=r"$x$", clip_on=False)

errorbar(modulus_output[:, 0], modulus_output[:, 3], modulus_output[:, 4], 
          fmt = "rs",
          markersize = 7.5,
          label=r"$y$", clip_on=False)

# hlines(quasi_x["modulus_x"], 1, 1e2, colors = "blue", ls = "-", label = r"Quasi-static, $x$")
# hlines(quasi_y["modulus_y"], 1, 1e2, colors = "red", ls = "--", label = r"Quasi-static, $y$")
xlabel(r"Strain rate ($\times$10$^7$/s)")
ylabel(r'Tensile modulus (GPa)')      
xscale("log")
xlim([1, 1e2])
# gca().set_xticks([1, 5, 10, 50, 100])
ylim([140, 240])
gca().set_yticks(linspace(140, 240, 6))
legend(loc = "center left")
title("(b)", fontsize=16)

subplot(2, 3, 5)
set_fig_properties([gca()])
rates = ["1e7",  "5e7", "1e8", "5e8", "1e9"]
sigma_output = np.zeros((5, 5))
for i in range(len(rates)):
    sigma_output[i, 0] = float(rates[i])/1e7
    sigma_output[i, 1] = tension_calc(rates[i] + "/x")["sigma_x_max_ave"]
    sigma_output[i, 2] = tension_calc(rates[i] + "/x")["sigma_x_max_bar"]
    sigma_output[i, 3] = tension_calc(rates[i] + "/y")["sigma_y_max_ave"]
    sigma_output[i, 4] = tension_calc(rates[i] + "/y")["sigma_y_max_bar"]    

errorbar(sigma_output[:, 0], sigma_output[:, 1], sigma_output[:, 2], 
          fmt = "bo",
          markersize = 7.5,
          label=r"$x$", clip_on=False)

errorbar(sigma_output[:, 0], sigma_output[:, 3], sigma_output[:, 4], 
          fmt = "rs",
          markersize = 7.5,
          label=r"$y$", clip_on=False)

x = np.arange(1, 1e2, 0.1)
C_x, m_x = optimize.curve_fit(f_2, sigma_output[:, 0], sigma_output[:, 1])[0]
C_y, m_y = optimize.curve_fit(f_2, sigma_output[:, 0], sigma_output[:, 3])[0]
y_x = C_x * x ** m_x
y_y = C_y * x ** m_y
plot(x, y_x, ls = "-", color = "blue")
plot(x, y_y, ls = "-", color = "red")

text(10, 3, r"$\sigma_{f} = 3.58\dot{\epsilon}^{0.036}$", color = "blue")
text(2, 9, r"$\sigma_{f} = 7.09\dot{\epsilon}^{0.054}$", color = "red")
# hlines(8.9, 1, 1e2, colors = "blue", ls = "-", label = r"DFT, $x$")
# hlines(14.1, 1, 1e2, colors = "red", ls = "--", label = r"DFT, $y$")

xlabel(r"Strain rate ($\times$10$^7$/s)")
ylabel(r'Tensile strength (GPa)')      
xscale("log")
xlim([1, 1e2])
ylim([2, 10])
gca().set_yticks(linspace(2, 10, 5))
# legend(loc = "lower left")
title("(c)", fontsize=16)

subplot(2, 3, 6)
set_fig_properties([gca()])
rates = ["1e7",  "5e7", "1e8", "5e8", "1e9"]
epsilon_output = np.zeros((5, 5))
for i in range(len(rates)):
    epsilon_output[i, 0] = float(rates[i])/1e7
    epsilon_output[i, 1] = tension_calc(rates[i] + "/x")["epsilon_x_max_ave"]
    epsilon_output[i, 2] = tension_calc(rates[i] + "/x")["epsilon_x_max_bar"]
    epsilon_output[i, 3] = tension_calc(rates[i] + "/y")["epsilon_y_max_ave"]
    epsilon_output[i, 4] = tension_calc(rates[i] + "/y")["epsilon_y_max_bar"]    

errorbar(epsilon_output[:, 0], epsilon_output[:, 1] * 100, epsilon_output[:, 2] * 100, 
          fmt = "bo",
          markersize = 7.5,
          label=r"$x$", clip_on=False)

errorbar(epsilon_output[:, 0], epsilon_output[:, 3] * 100, epsilon_output[:, 4] * 100, 
          fmt = "rs",
          markersize = 7.5,
          label=r"$y$", clip_on=False)

# hlines(6.5, 1, 1e2, colors = "blue", ls = "-", label = r"DFT, $x$")
# hlines(8, 1, 1e2, colors = "red", ls = "--", label = r"DFT, $y$")

xlabel(r"Strain rate ($\times$10$^7$/s)")
ylabel(r'Fracture strain (%)')      
xscale("log")
xlim([1, 1e2])
ylim([2, 5])
gca().set_yticks(linspace(2, 5, 4))
# legend(loc = "lower left")
title("(d)", fontsize=16)

subplots_adjust(wspace = 0.35, hspace = 0.35)
savefig("Tension_Strain_rate.pdf", bbox_inches='tight')
  



        
        
    
    
    
    
    
    
    
    
    
    
    
    