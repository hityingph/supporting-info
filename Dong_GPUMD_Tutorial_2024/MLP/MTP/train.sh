#!/bin/bash
#SBATCH -J mtp_gr
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o outfile_%J.mtp
#SBATCH -e errfile_%J.mtp

source /public3/soft/modules/module.sh
module load mpi/oneAPI/2022.1
MLP_EXE=~/MTP/mlip-2-master/bin/mlp
TMP_DIR=./out
mkdir -p $TMP_DIR

mpirun -n 64 $MLP_EXE train 18.mtp train.cfg --trained-pot-name=$TMP_DIR/pot.mtp 
