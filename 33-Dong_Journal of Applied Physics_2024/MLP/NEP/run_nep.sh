#!/bin/bash
#SBATCH -J gr_nep
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -o slurm_%J.out
#SBATCH -e slurm_%J.err

export PATH=//share/apps/GPUMD-3.4/src:$PATH
stdbuf -i0 -o0 -e0 nep
