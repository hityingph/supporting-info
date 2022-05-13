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
n=55
for ((i=50; i<$n; i++))
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
		echo 'monolayer-VP
   1.00000000000000     
     8.8911380349672111    0.0000000000000000    0.0000000000000000
     0.0000000000000000    8.9970467601671125    0.0000000000000000
     0.0000000000000000    0.0000000000000000   30.1621093750000000
   P 
    42
Direct
  0.6617947065489598  0.9770449962078170  0.5351187680652832
  0.7873283346545401  0.3480799736103480  0.5452049872515875
  0.9648624601894595  0.7248976285427647  0.4093295536146359
  0.9648631043370401  0.2751019335512322  0.5906706646430460
  0.2171814121431756  0.9769653935618682  0.4080590838932578
  0.2171819804123750  0.0230346015035136  0.5919406852280837
  0.0439728509082492  0.9006954296278353  0.4529915230545707
  0.0439730999313108  0.0993046653977471  0.5470083473655898
  0.3077316066021879  0.1563066801440111  0.4482305221377872
  0.3077321752839072  0.8436928458024071  0.5517692304828135
  0.4783118674426748  0.2140361132963536  0.3987631728327609
  0.4783125993648573  0.7859641230166530  0.6012366569435729
  0.5343619029101007  0.3838992686339031  0.4477603129636820
  0.5343625670806671  0.6161000125303978  0.5522399913015505
  0.7173078262753904  0.4771838794040495  0.4095482441090392
  0.7173084858694315  0.5228158185105085  0.5904521910026799
  0.6054558304279891  0.8293949102785654  0.3715167388627660
  0.6054551505802211  0.1706042638868950  0.6284827664556274
  0.8427238311156988  0.8479140081604551  0.3541719717760414
  0.8427228417422570  0.1520852915490814  0.6458277915510799
  0.8569269958540713  0.0826868135612048  0.3724213673879125
  0.8569276839401748  0.9173125844205805  0.6275783513467412
  0.0905789182034756  0.0993207698512592  0.3530166677213579
  0.0905795473229328  0.9006800034478266  0.6469831121865456
  0.1192771864359387  0.3338256691888113  0.3676929770832722
  0.1192790118953985  0.6661749534243908  0.6323074284336916
  0.3507440845443257  0.3391119046735592  0.3445904043058525
  0.3507462482800714  0.6608894861434602  0.6554097522184571
  0.3556346174363351  0.5683004273954767  0.3673827413688367
  0.3556353709169248  0.4317003147018064  0.6326177363814989
  0.5936418011730190  0.5972824223631280  0.3533205191401946
  0.5936423655209957  0.4027162054704500  0.6466796771466129
  0.5950541830087500  0.8009892446716437  0.4437993650096093
  0.5950553943249688  0.1990115418158974  0.5562002735995845
  0.8877740640165537  0.0894939648607967  0.4446692465342424
  0.8877749654930486  0.9105064831407124  0.5553304256113104
  0.1543885801225172  0.3517768874731196  0.4408652333558505
  0.1543888561566519  0.6482224945659582  0.5591351670370981
  0.3350105215728100  0.5325298500876795  0.4405443475727400
  0.3350111148322221  0.4674695819030552  0.5594560144359730
  0.6617938950679619  0.0229558084186144  0.4648809022374822
  0.7873269480603780  0.6519207422042020  0.4547952773496790    
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
