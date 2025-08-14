#!/bin/sh

#SBATCH --job-name=msmc3_cyanea
#SBATCH --account=passerinagenome
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --output=outs/stdout_msmc3_INBU_first4
#SBATCH --error=outs/stderr_msmc3_INBU_first4
#SBATCH --time=1-0:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem=0

BASEDIR=/project/passerinagenome/jrick/msmc

### Variables:
#Had to put loop in script because it wouldn't carry the values over
#for i in `cat INDS`; \
#do echo $i;\
#IND1=CEW16_022
#IND2=CEW17_120

#for s in `cat SCAFFOLDS`; \
#do echo $s;\
#`ls input/beagle/msmc_input.${IND}.${s}.gatk.beagle.phased.mask.txt` >> SCAFS_INPUT_${IND}; \
#done
#done

POP=INBU

P_PAR=1*2+25*1+1*2+1*3
NR_IND=4
POP_OR_IND=$POP
RUN_NAME=msmcc_${POP}_firstfour_071418
METHOD="samtools + mask"

find input/msmc_input.${POP}.*.txt -size 0 -delete
ls input/msmc_input.${POP}.*.txt > SCAFS_INPUT_${POP}

module load msmc2

### Report settings/parameters:
date
echo "Script: msmc_3_runMSMC_onepop.sh"
echo "Run name: $RUN_NAME"
echo "SNP calling method: $METHOD"
echo "Period setting: $P_PAR"
echo "Nr of individuals (1 or 2+): $NR_IND"
echo "Population or individuals ID: $POP_OR_IND"
echo "Individual: $IND"
echo "Scaffolds: SCAFS_INPUT_$POP"
echo "Iterations: 100"


if [ $NR_IND == 1 ]
	then
	echo "Running MSMC for one individual"
	MSMC_INPUT=`cat SCAFS_INPUT_${IND}_pseudohap`
	MSMC_OUTPUT=${BASEDIR}/output/msmc_output.${POP_OR_IND}.samtools.beagle.pseudohap.${RUN_NAME}

	#if [ -f $MSMC_INPUT ]
	#	then
			echo "MSMC_INPUTS: $MSMC_INPUT"
			echo "MSMC_OUTPUT: $MSMC_OUTPUT"
	#	else
	#		echo "MSMC_INPUT does not exist! Exiting now"
	#		exit 1
	#fi

	msmc2 -t 8 -p $P_PAR -i 100 -o $MSMC_OUTPUT -I 0,1 $MSMC_INPUT
	
	#mv $MSMC_OUTPUT*loop.txt $BASEDIR/output/$METHOD/log_and_loop/
	#mv $MSMC_OUTPUT*log $BASEDIR/output/$METHOD/log_and_loop/
else
	echo "Running MSMC for $NR_IND individuals"
#	MSMC_INPUT=`ls input/${POP}*.txt`
	MSMC_OUTPUT=${BASEDIR}/output/msmc_output.${POP_OR_IND}.${RUN_NAME}

	msmc2 -t 16 -p $P_PAR -i 100 -o ${MSMC_OUTPUT} -I 1,2,3,4,5,6,7,8 `cat SCAFS_INPUT_$POP`	
fi

echo "DONE_WITH_SCRIPT"
date
