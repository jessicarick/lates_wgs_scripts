#!/bin/sh

#SBATCH --job-name=msmc3
#SBATCH --account=latesgenomics
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=2-00:00:00
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

POP_FILE=$1
P_PAR=1*3+3*2+20*1+3*2+1*3
#NR_IND=`cat $POP_FILE | wc -l`
NR_IND=2
POP_OR_IND=`basename $POP_FILE .ind`
DATE=`date +%m%d%y`
RUN_NAME=msmc_latesWGS_2.bamlist_${POP_OR_IND}_${DATE}_20seg
METHOD="samtools + mask"

find input/msmc_input.latesWGS_2.bamlist.${POP_OR_IND}.*.txt -size 0 -delete
ls input/msmc_input.latesWGS_2.bamlist.${POP_OR_IND}.*.txt > input/SCAFS_INPUT_${POP_OR_IND}

module load gcc
module load msmc2
module load miniconda3

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


if [ "$NR_IND" -eq 4 ]; then
	echo "Running MSMC for $NR_IND individual"
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

	msmc2 -t 16 -p $P_PAR -i 100 -o ${MSMC_OUTPUT}.within1 -I 0,1,2,3,4,5 $MSMC_INPUT
	msmc2 -t 16 -p $P_PAR -i 100 -o ${MSMC_OUTPUT}.within2 -I 6,7 $MSMC_INPUT
	msmc2 -t 16 -p $P_PAR -i 100 -o ${MSMC_OUTPUT}.across -I 0-6,1-6,2-6,3-6,4-6,5-6,0-7,1-7,2-7,3-7,4-7,5-7 $MSMC_INPUT
	
	python3 /project/passerinagenome/scripts/msmc-tools-master/combineCrossCoal.py ${MSMC_OUTPUT}.across.final.txt ${MSMC_OUTPUT}.within1.final.txt ${MSMC_OUTPUT}.within2.final.txt > ${MSMC_OUTPUT}.crosscoal.final.txt

	#mv $MSMC_OUTPUT*loop.txt $BASEDIR/output/$METHOD/log_and_loop/
	#mv $MSMC_OUTPUT*log $BASEDIR/output/$METHOD/log_and_loop/
elif [ "$NR_IND" -eq 6 ]; then
	echo "Running MSMC for $NR_IND individuals"
	MSMC_INPUT=`cat input/SCAFS_INPUT_${POP_OR_IND}`
	MSMC_OUTPUT=${BASEDIR}/output/msmc_output.${RUN_NAME}
#	n=$(expr ${NR_IND} - 1)
#	INDEX=$(for num in `seq 1 ${n}`; do echo -n '${num},'; done; echo ${NR_IND})
        
	msmc2 -t 16 -p $P_PAR -i 100 -o ${MSMC_OUTPUT}.within1 -I 0,1,2,3,4,5 $MSMC_INPUT
        msmc2 -t 16 -p $P_PAR -i 100 -o ${MSMC_OUTPUT}.within2 -I 6,7,8,9,10,11 $MSMC_INPUT
        msmc2 -t 16 -p $P_PAR -i 100 -o ${MSMC_OUTPUT}.across -I 0-6,1-6,2-6,3-6,4-6,5-6,0-7,1-7,2-7,3-7,4-7,5-7,0-8,1-8,2-8,3-8,4-8,5-8,0-9,1-9,2-9,3-9,4-9,5-9,0-10,1-10,2-10,3-10,4-10,5-10,0-11,1-11,2-11,3-11,4-11,5-11 $MSMC_INPUT

	python3 /project/passerinagenome/scripts/msmc-tools-master/combineCrossCoal.py ${MSMC_OUTPUT}.across.final.txt ${MSMC_OUTPUT}.within1.final.txt ${MSMC_OUTPUT}.within2.final.txt > ${MSMC_OUTPUT}.crosscoal.final.txt

elif [ "$NR_IND" -eq 2 ]; then
	echo "Running MSMC for $NR_IND individuals"
        MSMC_INPUT=`cat input/SCAFS_INPUT_${POP_OR_IND}`
        MSMC_OUTPUT=${BASEDIR}/output/msmc_output.${RUN_NAME}
#       n=$(expr ${NR_IND} - 1)
#       INDEX=$(for num in `seq 1 ${n}`; do echo -n '${num},'; done; echo ${NR_IND})

       msmc2 -t 16 -p $P_PAR -i 100 -o ${MSMC_OUTPUT}.within1 -I 0,1 $MSMC_INPUT
       msmc2 -t 16 -p $P_PAR -i 100 -o ${MSMC_OUTPUT}.within2 -I 2,3 $MSMC_INPUT
       msmc2 -t 16 -p $P_PAR -i 100 -o ${MSMC_OUTPUT}.across -I 0-2,0-3,1-2,1-3  $MSMC_INPUT

       python3 /home/jrick/bin/msmc-tools/combineCrossCoal.py ${MSMC_OUTPUT}.across.final.txt ${MSMC_OUTPUT}.within1.final.txt ${MSMC_OUTPUT}.within2.final.txt > ${MSMC_OUTPUT}.crosscoal.final.txt

else
	echo "unknown number of individuals; exiting now"
	exit 1
fi


echo "DONE_WITH_SCRIPT"
date
