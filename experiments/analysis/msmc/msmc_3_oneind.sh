#!/bin/sh

#SBATCH --job-name=msmc3
#SBATCH --account=latesgenomics
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem=124G

BASEDIR=/project/latesgenomics/jrick/latesWGS_2018/msmc

### Variables:
#Had to put loop in script because it wouldn't carry the values over
#for i in `cat amoena_IND`
#	do echo $i
#	IND=$i

#	for s in `cat SCAFFOLDS.txt`
#		do echo $s
#		ls input/msmc_input.${IND}.${s}.txt >> input/SCAFS_INPUT_${IND}
#		done
#	done

P_PAR=1*3+3*2+30*1+3*2+1*3
NR_IND=1
POP_OR_IND=$1
DATE=`date +%m%d%y`
RUN_NAME=msmc_latesWGS_2.bamlist_${POP_OR_IND}_${DATE}_30seg
METHOD="samtools + mask"

find input/msmc_input.latesWGS_2.bamlist.${POP_OR_IND}.ASB*.txt -size 0 -delete
#ls input/msmc_input.${POP_OR_IND}.*.txt | grep -v 'LG9' > input/SCAFS_INPUT_${POP_OR_IND}_noLG9
ls input/msmc_input.latesWGS_2.bamlist.${POP_OR_IND}.ASB*.txt > input/SCAFS_INPUT_${POP_OR_IND}

module load msmc2

### Report settings/parameters:
date
echo "Script: msmc_3_onepop.sh"
echo "Run name: $RUN_NAME"
echo "SNP calling method: $METHOD"
echo "Period setting: $P_PAR"
echo "Nr of individuals (1 or 2+): $NR_IND"
echo "Population or individuals ID: $POP_OR_IND"
echo "Individual: "
echo "Scaffolds: SCAFS_INPUT_${POP_OR_IND}"
echo "Iterations: 100"


if [ $NR_IND == 1 ]
	then
	echo "Running MSMC for one individual"
	MSMC_INPUT=`cat input/SCAFS_INPUT_${POP_OR_IND}`
	MSMC_OUTPUT=${BASEDIR}/output/msmc_output.${RUN_NAME}

	if [ -f input/SCAFS_INPUT_${POP_OR_IND} ]
		then
			echo "MSMC_INPUTS: SCAFS_INPUT_${POP_OR_IND}"
			echo "MSMC_OUTPUT: $MSMC_OUTPUT"
		else
			echo "MSMC_INPUT does not exist! Exiting now"
			exit 1
	fi

	msmc2 -t 16 -p $P_PAR -i 100 -o $MSMC_OUTPUT -I 0,1 $MSMC_INPUT
	
	#mv $MSMC_OUTPUT*loop.txt $BASEDIR/output/$METHOD/log_and_loop/
	#mv $MSMC_OUTPUT*log $BASEDIR/output/$METHOD/log_and_loop/
else
	echo "Running MSMC for $NR_IND individuals"
	MSMC_INPUT=`cat input/SCAFS_INPUT_${POP_OR_IND}_noEmpty`
	MSMC_OUTPUT=output/msmc_output.${POP_OR_IND}.${RUN_NAME}
#	n=$(expr ${NR_IND} - 1)
#	INDEX=$(for num in `seq 1 ${n}`; do echo -n '${num},'; done; echo ${NR_IND})
	msmc2 -t 16 -p $P_PAR -i 100 -o ${MSMC_OUTPUT} $MSMC_INPUT	
fi


echo "DONE_WITH_SCRIPT"
date
