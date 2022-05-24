# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 10:11:08 2022

@author: benward
"""
import numpy as np
from ase.io import read, write
from pylab import *

def calculate_overlap(coord1, coord2, r12, r21):
    """
    Calculate the overlap area for center atom and atom on other layer
    coord1: atom on up layer
    coord2: atom on down layer
    r12: radius for center atom on up layer
    r21: radius for center atom on down layer
    """
    # print(coord1)
    # print(coord2)
    x1, y1 = coord1[0], coord1[1]
    x2, y2 = coord2[0], coord2[1] 
    l = np.sqrt((x1-x2)**2 + (y1-y2)**2)
    if l >= r12 + r21:                              #tangent or apart between two circles
        overlap_area = 0
    elif l <= max(r12,r21) - min(r12,r21):          #big circle include small circle 
        overlap_area = math.pi*min(r12,r21)**2
    else:                                           #intersect between two circles                                 
        alpha1 = math.acos((r12**2+l**2-r21**2) / (2*r12*l)) #The Law of Cosines
        alpha2 = math.acos((r21**2+l**2-r12**2) / (2*r21*l))
        area1 = alpha1*r12**2 - r12**2*math.sin(2*alpha1)/2
        area2 = alpha2*r21**2 - r21**2*math.sin(2*alpha2)/2
        overlap_area = area1 + area2 
    return overlap_area

def LRI_calc(up_coord):    
    overlap_area = 0
    for i in range(len(bottom_positions)):
        overlap_area += calculate_overlap(up_coord, bottom_positions[i], 0.71, 0.71)
    return overlap_area

atom_num = 24000
list = [122]
for m in list:
    atoms = read("%s.dump"%m, format="lammps-dump-text")
    up_positions = atoms.get_positions()[24000:][:,:2]
    bottom_positions = atoms.get_positions()[:24000][:,:2]
    LRI_output = np.zeros(atom_num)
    for i in range(atom_num):
        up_i = up_positions[i]
        print(i)
        S_i = LRI_calc(up_i)
        if i == 0:
            up_neighbor = up_positions[i+1:i+300]
        elif i < 300:
            up_neighbor = np.r_[up_positions[:i], up_positions[i+1:i+300]] 
        else:
            up_neighbor = np.r_[up_positions[i-300:i], up_positions[i+1:i+300]]       
        pair_distance = np.sqrt(np.sum((up_neighbor - up_i)**2, axis=1))
        index = np.where(np.array(pair_distance) < 1.6)
        index = index[0]
        up_jkl = np.zeros((len(index), 2))
        S_jkl = np.zeros(len(index))
        for j in range(len(index)):
            coord = up_neighbor[index[j]]
            lri = LRI_calc(coord)
            up_jkl[j] = coord
            S_jkl[j] = lri
        LRI = np.average(((S_i + S_jkl) - (0 + 1.58367)) / ((1.58367 + 1.58367) - (0 + 1.58367)))
        print(LRI)
        LRI_output[i] = LRI
    
    
    up_positions = atoms.get_positions()[24000:]
    with open("%s_LRI.dump"%m, "w") as fid:
        fid.write('ITEM: TIMESTEP\n')
        fid.write('{} \n'.format(60))
        fid.write('ITEM: NUMBER OF ATOMS\n')
        fid.write('{}\n'.format(24000))
        fid.write('ITEM: BOX BOUNDS pp pp pp\n')
    
    
        fid.write('{0:.10f}  {1:20.10f}\n'.format(-100, 445.951))
        fid.write('{0:.10f}  {1:20.10f}\n'.format(-100, 455.6))
        fid.write('{0:.10f}  {1:20.10f}\n'.format(0, 20))
    
        fid.write('ITEM: ATOMS id type x y z lri\n')
        for i in range(atom_num):
            fid.write('{0}   {1:.0f} {2:20.10f} {3:20.10f} {4:20.10f} {5:0.10f}\n'.format(i + 1, 1, up_positions[i, 0],
                                                                                up_positions[i, 1],
                                                                                up_positions[i, 2],
                                                                                LRI_output[i]))
    fid.close()    


    
    
    