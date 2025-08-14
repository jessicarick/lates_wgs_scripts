#!/bin/bash

module load gcc miniconda3 msmc2

date

POP_OR_IND=$1
P_PAR=1*3+3*2+30*1+3*2+1*3
NR_IND=1
#DATE=`date +%m%d%y`
DATE=021722
RUN_NAME=msmc_${POP_OR_IND}_${DATE}_30seg
METHOD="samtools + mask"
OUTDIR=/project/latesgenomics/jrick/latesWGS_2018/msmc
nchr=24

#input for the bootstrapping
BS_INPUT=$(for s in `cat CHROM_only`; do find ${OUTDIR}/input/ -maxdepth 1 -name "msmc_input.${POP_OR_IND}.${s}*.txt"; done)

#output from the bootstrapping 
BS_OUTPUT=${OUTDIR}/bootstrap/${POP_OR_IND}.${RUN_NAME}.bootstrap

echo "generating bootstraps for ${POP_OR_IND}.${RUN_NAME}"
#/home/jrick/bin/msmc-tools/multihetsep_bootstrap.py -n 50 -s 200000 --chunks_per_chromosome 10 --nr_chromosomes $nchr $BS_OUTPUT $BS_INPUT

cd ${OUTDIR}/bootstrap
for i in `seq 50`; do
	echo "working with bootstrap $i"
	ls -d ${POP_OR_IND}.${RUN_NAME}.bootstrap_$i/ > ${POP_OR_IND}.${RUN_NAME}.${i}.bs_file_list.txt
	##### run msmc ####
	MSMC_BS=$(for x in `cat ${POP_OR_IND}.${RUN_NAME}.${i}.bs_file_list.txt`; do find $x -maxdepth 2 -name "bootstrap_multihetsep*.txt"; done)

	#MSMC_BS=`ls ${OUTDIR}/bootstrap/${FOLDER}/*.txt`
	MSMC_OUTPUT=${OUTDIR}/bootstrap/msmc_output.${RUN_NAME}.bs$i

	echo "running msmc2 on bootstraps for ${POP_OR_IND}.${RUN_NAME}"
	msmc2 -t 16 -p $P_PAR -o $MSMC_OUTPUT -I 0,1 $MSMC_BS
done

echo "done with msmc bootstraps for ${POP_OR_IND}.${RUN_NAME}"
date
