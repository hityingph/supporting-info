# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 18:22:39 2022

@author: benward
"""

import numpy as np
import os
from ase.io import read, write
from pylab import *
import pandas as pd

output = np.zeros((251, 9))
for i in range(251):
    print(i)
    atoms = read("MD/%s.dump"%i, format="lammps-dump-text")
    bottom_xyz = atoms.get_positions()[:24000]
    up_xyz = atoms.get_positions()[24000:]
    bottom_x = bottom_xyz[:,0]
    bottom_y = bottom_xyz[:,1]
    up_x = up_xyz[:,0]
    up_y = up_xyz[:,1]
    
    up_lx = np.average(np.sort(up_x)[-120:] - np.sort(up_x)[:120])  
    up_ly = np.average(np.sort(up_y)[-100:] - np.sort(up_y)[:100])  
    bottom_lx = np.average(np.sort(bottom_x)[-120:] - np.sort(bottom_x)[:120])  
    bottom_ly = np.average(np.sort(bottom_y)[-100:] - np.sort(bottom_y)[:100])
    
    total_force = atoms.get_forces()
    total_force_x = -sum(total_force[:,0]) #unit in eV/A
    total_force_y = -sum(total_force[:,1]) #unit in eV/A
    total_stress_x = total_force_x / 3.35 / 2 / bottom_ly * 160.2 #unit in GPa
    total_stress_y = total_force_y / 3.35 / 2 / bottom_lx * 160.2 #unit in GPa
    
    bottom_force = atoms.get_forces()[:24000]
    bottom_force_x = -sum(bottom_force[:,0]) #unit in eV/A
    bottom_force_y = -sum(bottom_force[:,1]) #unit in eV/A
    bottom_stress_x = bottom_force_x / 3.35 / bottom_ly * 160.2 #unit in GPa
    bottom_stress_y = bottom_force_y / 3.35 / bottom_lx * 160.2 #unit in GPa
    
    
    output[i] = np.array([i, 
                          bottom_lx, bottom_ly, 
                          up_lx, up_ly, 
                          bottom_stress_x, bottom_stress_y,
                          total_stress_x, total_stress_y])

aw = 1.5
fs = 16
lw = 1.5
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
        
figure(figsize=(12, 10))
subplot(2,2,1)
set_fig_properties([gca()])  
strain = np.arange(251)*2/10000
up_strain_x = output[:, 3] / output[0, 3] - 1
up_strain_y = output[:, 4] / output[0, 4] - 1

plot(strain*100, up_strain_y*100,
      color='k',
      linestyle='-',
      linewidth=1.0, 
      marker='o', 
      markersize=8, 
      markeredgecolor=None, 
      markerfacecolor= "C2", label = r"Armchair")

plot(strain*100, up_strain_x*100, 
      color='k',
      linestyle='-',
      linewidth=1.0, 
      marker='o', 
      markersize=8, 
      markeredgecolor=None, 
      markerfacecolor= "C3", label = "Zigzag")

xlim([0, 5])
gca().set_xticks(linspace(0, 5, 6))
ylim([0, 2])
gca().set_yticks(linspace(0, 2, 5))
xlabel(r'Strain of bottom layer (%)')
ylabel(r'Strain of up layer (%)')
legend(loc = 'upper right')
title('(a)')

subplot(2,2,2)
set_fig_properties([gca()])  
strain = np.arange(251)*2/10000
bottom_stress_x = output[:, 5] 
bottom_stress_y = output[:, 6] 

plot(strain*100, bottom_stress_y,
      color='k',
      linestyle='-',
      linewidth=1.0, 
      marker='o', 
      markersize=8, 
      markeredgecolor=None, 
      markerfacecolor= "C2", label = r"Armchair")

plot(strain*100, bottom_stress_x, 
      color='k',
      linestyle='-',
      linewidth=1.0, 
      marker='o', 
      markersize=8, 
      markeredgecolor=None, 
      markerfacecolor= "C3", label = "Zigzag")

xlim([0, 5])
gca().set_xticks(linspace(0, 5, 6))
ylim([-1, 2])
gca().set_yticks(linspace(-1, 2, 4))
hlines(0, 0, 5, linestyles='--', linewidth = 1.5, color = 'grey')
xlabel(r'Strain of bottom layer (%)')
ylabel(r'Stress of bottom layer (GPa)')
legend(loc = 'upper right')
title('(b)')

subplot(2,2,2)
set_fig_properties([gca()])  
strain = np.arange(251)*2/10000
Total_stress_x = output[:, 7] 
Total_stress_y = output[:, 8] 

plot(strain*100, Total_stress_y,
      color='k',
      linestyle='-',
      linewidth=1.0, 
      marker='o', 
      markersize=8, 
      markeredgecolor=None, 
      markerfacecolor= "C2", label = r"Armchair")

plot(strain*100, Total_stress_x, 
      color='k',
      linestyle='-',
      linewidth=1.0, 
      marker='o', 
      markersize=8, 
      markeredgecolor=None, 
      markerfacecolor= "C3", label = "Zigzag")

xlim([0, 5])
gca().set_xticks(linspace(0, 5, 6))
ylim([-1, 10])
gca().set_yticks(linspace(0, 10, 6))
hlines(0, 0, 5, linestyles='--', linewidth = 1.5, color = 'grey')
xlabel(r'Strain of bottom layer (%)')
ylabel(r'Total stress (GPa)')
legend(loc = 'upper right')
title('(b)')
plt.subplots_adjust(wspace = 0.35, hspace=0.3)
savefig("Bilayer_Dislocation.pdf", bbox_inches='tight')

data = {'strain': strain,
        'up_strain_armchair':up_strain_y,
        'up_strain_zigzag':up_strain_y,
        'bottom_stress_armchair':bottom_stress_y,
        'bottom_stress_zigzag':bottom_stress_x,
        'total_stress_armchair':Total_stress_y,
        'total_stress_zigzag':Total_stress_x}
frame = pd.DataFrame(data)
frame.to_csv("BLG_Dislocation.csv")
