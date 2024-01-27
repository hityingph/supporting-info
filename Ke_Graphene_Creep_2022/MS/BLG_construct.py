# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 15:21:32 2022

@author: benward
"""

import numpy as np
import math
from ase.build import graphene_nanoribbon
from pylab import *
from scipy import interpolate
from ase.io import read, write


gr = graphene_nanoribbon(1, 1, type='armchair', sheet=True, vacuum=3.35/2, C_C=1.42)
gr.euler_rotate(theta=90)
l = gr.cell.lengths()
gr.cell = gr.cell.new((l[0], l[2], l[1]))
gr.center()
write("gr.vasp", gr)
l = gr.cell.lengths()
lx = l[0]
ly = l[1]

# bottom_layer = gr.get_positions()
# up_layer = gr.get_positions()
# for i in range(len(up_layer)):
#     y = up_layer[i,1] + 1.42
#     z = up_layer[i,2] + 3.35
#     if y > ly:
#         up_layer[i,1] = y - ly
#     else:
#         up_layer[i,1] = y
#     up_layer[i,2] = z

# max_atom_id = len(up_layer) + len(bottom_layer)
# max_atom_type = 2

# fout = open("BilayerGraphene_ABstacking.data", "w")
# fout.write("# LAMMPS data of bilayerGraphene with AB stacking mode" + 
#             " written by Penghua Ying (hityingph@163.com) \n")
# fout.write(str(max_atom_id) + "\t atoms \n")
# fout.write(str(max_atom_type) + "\t atom types \n\n")
# fout.write(format(0, ".4f") + "\t" + format(lx, ".4f") + "\t xlo xhi \n")
# fout.write(format(0, ".4f") + "\t" + format(ly, ".4f") + "\t ylo yhi \n")
# fout.write(format(0, ".4f") + "\t" + format(20, ".4f") + "\t zlo zhi \n\n\n")
# fout.write("Atoms \n\n")
# lines = ""
# atom_id = 0
# for i in range(len(bottom_layer)):
#     atom_id += 1
#     lines += " ".join([str(atom_id), 
#                      str(1), 
#                      format(bottom_layer[i,0], ".4f"), 
#                      format(bottom_layer[i,1], ".4f"), 
#                      format(bottom_layer[i,2], ".4f"), 
#                      "\n"])
# for i in range(len(up_layer)):
#     atom_id += 1
#     lines += " ".join([str(atom_id), 
#                      str(2), 
#                      format(up_layer[i,0], ".4f"), 
#                      format(up_layer[i,1], ".4f"), 
#                      format(up_layer[i,2], ".4f"), 
#                      "\n"])
# fout.write(lines)
# fout.close()

