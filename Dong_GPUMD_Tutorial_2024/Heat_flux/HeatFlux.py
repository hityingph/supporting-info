# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 14:13:33 2022

@author: benward
"""
from pylab import *
from thermo.gpumd.data import load_shc,load_compute

aw = 1.5
fs = 16
lw = 4
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
        ax.tick_params(axis='both', direction='out', right=False, top=False)

def lammps_heatflux(path):
    dt = 0.001  # ps 
    Ns = 1000  # Sample interval
    thermo = loadtxt(path + "/compute_Energy_Temp.out")
    jp = loadtxt(path + "/compute_HeatFlux.out")
    BLOCK_LENGTH = 426
    Ein = thermo[500:, 1]
    Eout = thermo[500:, 2]
    Etol = (Eout - Ein) / 2 / 1000 #in units of KeV
    Etol = Etol - Etol[0]
    t = dt * np.arange(1, len(Etol) + 1) * Ns / 1000 #unit in ns
    jpy = jp[500:, 2] - jp[500:, 5]
    jpy = jpy / BLOCK_LENGTH / 10 * 1000 #in units of eV/ns
    accum_jpy = cumsum(jpy) * 0.001 / 1000 #in units of KeV
    return t, accum_jpy, Etol

def gpumd_heatflux(path):
    dt = 0.001  # ps 
    Ns = 1000  # Sample interval
    BLOCK_LENGTH = 426
    compute = load_compute(['T', 'jp'], filename=path + "/compute.out")
    Ein = compute['Ein']
    Eout = compute['Eout']
    jp = compute['jp']
    Ein = compute['Ein'][500:]
    Eout = compute['Eout'][500:]
    Etol = (Eout - Ein) / 2 / 1000 #in units of KeV
    Etol = Etol - Etol[0]
    t = dt * np.arange(1, len(Etol) + 1) * Ns / 1000 #unit in ns
    jp = compute['jp'][500:]
    jpy = jp[:, 15:25] # in units of eV^{3/2} amu^{-1/2}
    jpy = jpy / BLOCK_LENGTH * 98207.9 #in units of eV/ns
    jpy_ave = np.mean(jpy, axis = 1)
    accum_jpy = cumsum(jpy_ave) * 0.001 / 1000 #in units of KeV
    return t, accum_jpy, Etol

t_nep, jp_nep, etol_nep = gpumd_heatflux("GPUMD_NEP")
t_dp, jp_dp, etol_dp = lammps_heatflux("LAMMPS_DP")
t_mtp, jp_mtp, etol_mtp = lammps_heatflux("LAMMPS_MTP")
t_nepcpu, jp_nepcpu, etol_nepcpu = lammps_heatflux("LAMMPS_NEP")


figure(figsize=(12, 10))
subplot(2, 2, 1)
set_fig_properties([gca()])
plot(t_nep, jp_nep, lw = lw, ls = "-", label = "NEP, from atoms")
plot(t_nep, etol_nep, lw = lw, ls = "--", label = "NEP, from thermostats")
legend(frameon=False)
xlim([0, 1.5])
gca().set_xticks(linspace(0, 1.5, 4))
ylim([0, 2])
gca().set_yticks(linspace(0, 2, 5))
xlabel('Time (ns)')
ylabel('Accumulated heat (keV)')
title("(a) GPUMD-NEP")

subplot(2, 2, 2)
set_fig_properties([gca()])
plot(t_nepcpu, -1 * jp_nepcpu, lw = lw, ls = "-", label = "NEP, from atoms")
plot(t_nepcpu, etol_nepcpu, lw = lw, ls = "--", label = "NEP, from thermostats")
legend(frameon=False)
xlim([0, 1.5])
gca().set_xticks(linspace(0, 1.5, 4))
ylim([0, 2])
gca().set_yticks(linspace(0, 2, 5))
xlabel('Time (ns)')
ylabel('Accumulated heat (keV)')
title("(b) LAMMPS-NEP")

subplot(2, 2, 3)
set_fig_properties([gca()])
plot(t_dp, jp_dp, lw = lw, ls = "-", label = "DP, from atoms")
plot(t_dp, etol_dp, lw = lw, ls = "--", label = "DP, from thermostats")
legend(frameon=False)
xlim([0, 1.5])
gca().set_xticks(linspace(0, 1.5, 4))
ylim([0, 2])
gca().set_yticks(linspace(0, 2, 5))
xlabel('Time (ns)')
ylabel('Accumulated heat (keV)')
title("(c) LAMMPS-DP")

subplot(2, 2, 4)
set_fig_properties([gca()])
plot(t_mtp, -1 * jp_mtp, lw = lw, ls = "-", label = "MTP, from atoms")
plot(t_mtp, etol_mtp, lw = lw, ls = "--", label = "MTP, from thermostats")
legend(frameon=False)
xlim([0, 1.5])
gca().set_xticks(linspace(0, 1.5, 4))
ylim([0, 2])
gca().set_yticks(linspace(0, 2, 5))
xlabel('Time (ns)')
ylabel('Accumulated heat (keV)')
title("(d) LAMMPS-MTP")

subplots_adjust(wspace = 0.3, hspace = 0.35)
savefig("MLPs_HeatFluex.pdf", bbox_inches='tight')
