#!/bin/bash
#不同应变下的串行计算
#SBATCH -J Serial-computation
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64
source /public3/soft/modules/module.sh
module load mpi/intel/17.0.7-thc-public3
export PATH=/public3/home/sc52255/vasp.6.1.0/vasp.6.1.0/bin:$PATH

#脚本文件路径
filepath=`pwd`
# 计算的应变个数
n=50
for ((i=0; i<$n; i++))
do
	mkdir strain-$i%
	cd strain-$i%
	mkdir 1.optimize
	cd 1.optimize
##############################################################################	
# INCAR设置
##############################################################################
	# 如果是最初始结构
	if [[ $i == 0 ]];then
		echo 'SYSTEM = Everything
LWAVE = F               		!Output wave function: off
ICHARG = 2             			!Read CHGCAR (2:off) (11:on): off
LCHARG = F             			!Output charge density: off 
PREC = Accurate         		!Precision: Normal; Accurate 
LREAL = F
NPAR = 8                		!Applicable to Beijing Supercomputing

GGA = PE                		!With PBE potential
METAGGA = SCAN          		!SCAN
LASPH = T
LMIXTAU = T
ADDGRID = T

LUSE_VDW = .TRUE.       		!vdW: rVV10
BPARAM = 6.3            		!default but can be overwritten by this tag
CPARAM = 0.0093         		!default but can be overwritten by this tag

ISIF = 3
IOPTCELL = 1 0 0 0 1 0 0 0 1
#VCA = 0.667 0.333 1 1 1 1		!Atomic impurity
ISYM = 0                		!Symmetry: off

ENCUT = 500 eV
ISPIN = 1               		!Spin (1:off) (2:on): off
#MAGMOM = 52*0 24*4 96*0 4 0 0 4 0 4 4 0
NELM = 60
EDIFF = 1E-05          			!Electronic convergence

ISMEAR = 0             			!(0:semiconductor) (-1:metal)
#SIGMA = 0.05
ALGO = F                		!Fast

IBRION = 2              		!Ion relaxation (CG): on
NSW = 300              			!Ion relaxation: on  
#POTIM = 0.1            		!Step size
EDIFFG = -0.01         			!Ion convergence

#EMIN = -10              		!DOS calculation
#EMAX = 5
#NEDOS = 3001
#LORBIT = 11             		!Projected to px, py and pz orbitals
'>INCAR
	# 如果不是初始结构
	else
		echo 'SYSTEM = Everything
LWAVE = F               		!Output wave function: off
ICHARG = 2             			!Read CHGCAR (2:off) (11:on): off
LCHARG = F             			!Output charge density: off 
PREC = Accurate         		!Precision: Normal; Accurate 
LREAL = F
NPAR = 8                		!Applicable to Beijing Supercomputing

GGA = PE                		!With PBE potential
METAGGA = SCAN          		!SCAN
LASPH = T
LMIXTAU = T
ADDGRID = T

LUSE_VDW = .TRUE.       		!vdW: rVV10
BPARAM = 6.3            		!default but can be overwritten by this tag
CPARAM = 0.0093         		!default but can be overwritten by this tag

ISIF = 3
IOPTCELL = 0 0 0 0 1 0 0 0 1
#VCA = 0.667 0.333 1 1 1 1		!Atomic impurity
ISYM = 0                		!Symmetry: off

ENCUT = 500 eV
ISPIN = 1               		!Spin (1:off) (2:on): off
#MAGMOM = 52*0 24*4 96*0 4 0 0 4 0 4 4 0
NELM = 60
EDIFF = 1E-05          			!Electronic convergence

ISMEAR = 0             			!(0:semiconductor) (-1:metal)
#SIGMA = 0.05
ALGO = F                		!Fast

IBRION = 2              		!Ion relaxation (CG): on
NSW = 300              			!Ion relaxation: on  
#POTIM = 0.1            		!Step size
EDIFFG = -0.01         			!Ion convergence

#EMIN = -10              		!DOS calculation
#EMAX = 5
#NEDOS = 3001
#LORBIT = 11             		!Projected to px, py and pz orbitals
'>INCAR
	fi
##############################################################################
# POSCAR设置
##############################################################################
	# 如果是最初始结构
	if [[ $i == 0 ]];then
		echo 'bulk-VP
1
9.21 0.0 0.0
0.0 9.128 0.0
-2.9621350956 0.0 21.6916851507
P 
84 
Direct
0.28456 0.63814 0.57338
0.21544 0.63814 0.92662
0.71544 0.36186 0.42662
0.78456 0.36186 0.07338
0.03108 0.89041 0.57221
0.46892 0.89041 0.92779
0.96892 0.10959 0.42779
0.53108 0.10959 0.07221
0.19227 0.28052 0.61625
0.30773 0.28052 0.88375
0.80773 0.71948 0.38375
0.69227 0.71948 0.11625
0.69208 0.7775 0.61646
0.80792 0.7775 0.88354
0.30792 0.2225 0.38354
0.19208 0.2225 0.11646
0.52618 0.14747 0.53369
0.97382 0.14747 0.96631
0.47382 0.85253 0.46631
0.02618 0.85253 0.03369
0.79388 0.90516 0.54389
0.70612 0.90516 0.95611
0.20612 0.09484 0.45611
0.29388 0.09484 0.04389
0.29619 0.40883 0.54459
0.20381 0.40883 0.95541
0.70381 0.59117 0.45541
0.79619 0.59117 0.04459
0.53294 0.37302 0.56862
0.96706 0.37302 0.93138
0.46706 0.62698 0.43138
0.03294 0.62698 0.06862
0.37906 0.17906 0.67427
0.12094 0.17906 0.82573
0.62094 0.82094 0.32573
0.87906 0.82094 0.17427
0.60742 0.94853 0.67255
0.89258 0.94853 0.82745
0.39258 0.05147 0.32745
0.10742 0.05147 0.17255
0.1298 0.44736 0.68018
0.3702 0.44736 0.81982
0.8702 0.55264 0.31982
0.6298 0.55264 0.18018
0.42066 0.01547 0.60449
0.07934 0.01547 0.89551
0.57934 0.98453 0.39551
0.92066 0.98453 0.10449
0.04542 0.65811 0.54454
0.45458 0.65811 0.95546
0.95458 0.34189 0.45546
0.54542 0.34189 0.04454
0.02663 0.8891 0.67228
0.47337 0.8891 0.82772
0.97337 0.1109 0.32772
0.52663 0.1109 0.17228
0.94246 0.53228 0.61668
0.55754 0.53228 0.88332
0.05754 0.46772 0.38332
0.44246 0.46772 0.11668
0.76273 0.13663 0.56655
0.73727 0.13663 0.93345
0.23727 0.86337 0.43345
0.26273 0.86337 0.06655
0.25465 0.82242 0.70066
0.24535 0.82242 0.79934
0.74535 0.17758 0.29934
0.75465 0.17758 0.20066
0.3124 0.60138 0.67345
0.1876 0.60138 0.82655
0.6876 0.39862 0.32655
0.8124 0.39862 0.17345
0.8803 0.70036 0.68115
0.6197 0.70036 0.81885
0.1197 0.29964 0.31885
0.3803 0.29964 0.18115
0.7492 0.14559 0.66767
0.7508 0.14559 0.83233
0.2508 0.85441 0.33233
0.2492 0.85441 0.16767
0.57248 0.3252 0.66938
0.92752 0.3252 0.83062
0.42752 0.6748 0.33062
0.07248 0.6748 0.16938   
'>POSCAR
	# 如果不是初始结构
	else
		# 将上一步CONTCAR复制为这一步的POSCAR
		j=`echo "$i-1" | bc`
		cp $filepath/strain-$j%/1.optimize/CONTCAR temporary-1
		# 读取上一步CONTCAR中的晶格常数a和b沿x轴的分量
		ax=`awk 'NR==3 {print $1}' temporary-1`
		bx=`awk 'NR==4 {print $1}' temporary-1`
		# 读取第1步的CONTCAR中的晶格常数a和b沿x轴的分量
		a_ax=`awk 'NR==3 {print $1}' $filepath/strain-0%/1.optimize/CONTCAR`
		a_bx=`awk 'NR==4 {print $1}' $filepath/strain-0%/1.optimize/CONTCAR`
		# 计算工程应变
		delta_ax=`echo "scale=8; $a_ax/100" | bc`
		delta_bx=`echo "scale=8; $a_bx/100" | bc`
		# 施加x轴方向的应变
		new_ax=`echo "scale=8; $ax+$delta_ax" | bc`
		new_bx=`echo "scale=8; $bx+$delta_bx" | bc`
		awk 'NR==3 {$1="'$new_ax'"} {print}' temporary-1 > temporary-2
		awk 'NR==4 {$1="'$new_bx'"} {print}' temporary-2 > POSCAR
		rm temporary-1 temporary-2
	fi
##############################################################################
# KPOINTS，POTCAR设置
##############################################################################
	# 使用vaspkit生成KPOINTS和POTCAR，此时K点密度为0.04	
	echo -e "102\n2\n0.04\n" |vaspkit

##############################################################################	
# 提交vasp作业
##############################################################################
	srun vasp_std
	echo "计算strain-$i%已经结束！。。。。。。3"
	echo "计算strain-$i%已经结束！。。。。。。2"
	echo "计算strain-$i%已经结束！。。。。。。1"
	echo "计算strain-$i%已经结束！"
	# 退出当前文件夹
	cd ../..
done
