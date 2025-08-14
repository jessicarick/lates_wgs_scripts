#!/bin/sh

#SBATCH --job-name=msmc3
#SBATCH --account=latesgenomics
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
# #SBATCH --output=outs/stdout_msmc3_LMIC
# #SBATCH --error=outs/stderr_msmc3_LMIC
#SBATCH --time=0-12:00:00
#SBATCH --mem=124G
#SBATCH --mail-type=ALL

BASEDIR=/project/latesgenomics/jrick/latesWGS_2018/msmc

POP=$1

### Variables:
#Had to put loop in script because it wouldn't carry the values over
#for i in `cat lmic.ind`
#	do echo $i
#	IND=$i

	for s in `cat CHROM_only`
		do echo $s
		ls ${BASEDIR}/input/msmc_input.latesWGS_2.bamlist.${POP}.${s}.txt >> SCAFS_INPUT_${POP}
		done
#	done

P_PAR=1*3+1*2+30*1+3*2+1*3
NR_IND=1
POP_OR_IND=$POP
RUN_NAME=msmc_050823_${POP}
METHOD="samtools + whatshap + mask"

#find input/pseudohap/msmc_input.${IND}.*.txt -size 0 -delete
#ls input/pseudohap/msmc_input.${IND}.*.txt > SCAFS_INPUT_${IND}_pseudohap

module load msmc2

### Report settings/parameters:
date
echo "Script: msmc_3_runMSMC_onepop.sh"
echo "Run name: $RUN_NAME"
echo "SNP calling method: $METHOD"
echo "Period setting: $P_PAR"
echo "Nr of individuals (1 or 2+): $NR_IND"
echo "Population or individuals ID: $POP_OR_IND"
#echo "Individual: $IND"
echo "Scaffolds: `cat SCAFS_INPUT_${POP}`"
echo "Iterations: 100"


if [ $NR_IND == 1 ]
	then
	echo "Running MSMC for one individual"
	MSMC_INPUT=`cat SCAFS_INPUT_${POP}`
	MSMC_OUTPUT=${BASEDIR}/output/msmc_output.${RUN_NAME}

	#if [ -f $MSMC_INPUT ]
	#	then
			echo "MSMC_INPUTS: $MSMC_INPUT"
			echo "MSMC_OUTPUT: $MSMC_OUTPUT"
	#	else
	#		echo "MSMC_INPUT does not exist! Exiting now"
	#		exit 1
	#fi

	msmc2 -t 16 -p $P_PAR -i 100 -o $MSMC_OUTPUT -I 0,1 $MSMC_INPUT
	
	#mv $MSMC_OUTPUT*loop.txt $BASEDIR/output/$METHOD/log_and_loop/
	#mv $MSMC_OUTPUT*log $BASEDIR/output/$METHOD/log_and_loop/
else
	echo "Running MSMC for $NR_IND individuals"
	MSMC_INPUT=`cat SCAFS_INPUT_${POP}`
	MSMC_OUTPUT=${BASEDIR}/output/msmc_output.${POP_OR_IND}.${RUN_NAME}

	msmc2 -t 16 -p $P_PAR -i 100 -o ${MSMC_OUTPUT} $MSMC_INPUT	
fi

echo "DONE_WITH_SCRIPT"
date
