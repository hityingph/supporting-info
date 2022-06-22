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
IOPTCELL = 1 0 0 0 1 0 0 0 0
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
IOPTCELL = 0 0 0 0 1 0 0 0 0
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
		echo 'Double-VP
1
9.21 0.0 0.0
0.0 9.128 0.0
-2.9621350956 0.0 41.6917
P 
84 
Direct
0.273238451 0.63814 0.5381785447
0.1496185364 0.63814 0.7219650265
0.7267614344 0.36186 0.4618210991
0.850381349 0.36186 0.2780346174
0.0199389654 0.89041 0.5375698079
0.402918022 0.89041 0.7225737633
0.9800609201 0.10959 0.4624298359
0.5970818634 0.10959 0.2774258806
0.1743342189 0.28052 0.560483285
0.2485227685 0.28052 0.6996602861
0.8256656665 0.71948 0.4395163588
0.751477117 0.71948 0.3003393577
0.6741118189 0.7775 0.5605925455
0.7487451685 0.7775 0.6995510257
0.3258880665 0.2225 0.4394070984
0.2512547169 0.2225 0.3004486181
0.5209820543 0.14747 0.5175283197
0.9018749332 0.14747 0.7426152514
0.4790178312 0.85253 0.4824713241
0.0981249523 0.85253 0.2573843924
0.7871083391 0.90516 0.5228352559
0.6357486483 0.90516 0.7373083153
0.2128915463 0.09484 0.477164388
0.3642512372 0.09484 0.2626913285
0.2893103391 0.40883 0.5231994574
0.1335466484 0.40883 0.7369441138
0.7106895464 0.59117 0.4768001865
0.8664532371 0.59117 0.26305553
0.5223528514 0.37302 0.5357019745
0.900504136 0.37302 0.7244415967
0.477647034 0.62698 0.4642976693
0.0994957494 0.62698 0.2755580472
0.352172557 0.17906 0.5906703863
0.0706844304 0.17906 0.6694731848
0.6478273284 0.82094 0.4093292575
0.9293154551 0.82094 0.330526459
0.5807979286 0.94853 0.5897754912
0.8420590588 0.94853 0.6703680799
0.4192019568 0.05147 0.4102241526
0.1579408266 0.05147 0.3296315639
0.102000728 0.44736 0.5937452876
0.3208562594 0.44736 0.6663982836
0.8979991575 0.55264 0.4062543563
0.679143626 0.55264 0.3336013602
0.4045386199 0.01547 0.5543646999
0.0183183675 0.01547 0.7057788713
0.5954612656 0.98453 0.445634944
0.9816815179 0.98453 0.2942207725
0.0385480534 0.65811 0.523173443
0.3843089341 0.65811 0.7369701282
0.9614518321 0.34189 0.4768262009
0.6156909514 0.34189 0.2630295156
4.9585783713e-05 0.8891 0.5896350135
0.4228074016 0.8891 0.6705085576
0.9999502997 0.1109 0.4103646303
0.5771924838 0.1109 0.3294910862
0.924457876 0.53228 0.5607070088
0.4983991114 0.53228 0.6994365624
0.0755420094 0.46772 0.439292635
0.5016007741 0.46772 0.3005630815
0.752462223 0.13663 0.5346249786
0.6703947644 0.13663 0.7255185925
0.2475376625 0.86337 0.4653746652
0.329605121 0.86337 0.2744810513
0.2236909549 0.82242 0.6044007828
0.1991660325 0.82242 0.6557427883
0.7763089306 0.17758 0.395598861
0.8008338529 0.17758 0.3442568555
0.2856390714 0.60138 0.5902437503
0.137217916 0.60138 0.6698998209
0.714360814 0.39862 0.4097558935
0.8627819694 0.39862 0.330099823
0.8523510708 0.70036 0.5942499668
0.5705059167 0.70036 0.6658936044
0.1476488147 0.29964 0.405749677
0.4294939688 0.29964 0.3341060395
0.7233308433 0.14559 0.5872364865
0.6995261441 0.14559 0.6729070847
0.2766690421 0.85441 0.4127631573
0.3004737413 0.85441 0.3270925592
0.5463470146 0.3252 0.5881261787
0.8765099728 0.3252 0.6720173924
0.4536528709 0.6748 0.4118734651
0.1234899126 0.6748 0.3279822514   
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
