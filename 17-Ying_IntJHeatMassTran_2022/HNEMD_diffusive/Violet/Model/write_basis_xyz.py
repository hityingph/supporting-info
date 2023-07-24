# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 20:28:33 2021

@author: benward
"""
from ase.io import read ,write
from thermo.gpumd.preproc import add_basis, repeat
from thermo.gpumd.io import create_basis, create_kpoints, ase_atoms_to_gpumd
from thermo.gpumd.data import load_omega2
from pylab import *


BP_uc = read("POSCAR")
l = BP_uc.cell.lengths()
BP_uc.pbc = [True, True, False]
# BP_uc.euler_rotate(theta=90,psi=90.0)
# BP_uc.set_cell([l[2], l[0], l[1]])
# BP_uc.center()
# write("BP_cell.vasp",BP_uc)
target_size = 250 #corresponding 25 nm
rx = int(target_size/l[0])+1
ry = int(target_size/l[1])+1
BP = BP_uc.repeat([rx,ry,1])
BP.wrap()
write("BP_cubic.vasp",BP)

ase_atoms_to_gpumd(BP, M=200, cutoff=9)