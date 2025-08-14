#!/bin/bash -l

### Variables:
IND=`cat $1`
POP=$2

#IFS='.' read -ra SUBSTRINGS <<< "$VCF"
#IFS='/' read -ra SUBSTRINGS2 <<< "${SUBSTRINGS[0]}"

#TOT=`echo ${#SUBSTRINGS2[@]}`
#IND_NR=`expr $TOT - 1`
#IND=${SUBSTRINGS2[$IND_NR]}

#SCAFFOLD=${SUBSTRINGS[1]}
#PHASING=${SUBSTRINGS[3]}
PHASING=unphased
METHOD=samtools
#METHOD=${SUBSTRINGS[4]}

module load samtools
module load miniconda3
PATH=$PATH:/project/passerinagenome/scripts/msmc-tools-master

BASEDIR=/project/latesgenomics/jrick/latesWGS_2018/msmc


for s in `cat CHROM_only`
	do SCAFFOLD=$s

	MSMC_INPUT=${BASEDIR}/input/msmc_input.${POP}.${SCAFFOLD}.txt

	printf "\n \n \n \n"
	date
	echo "Script: msmc_2_generateInput_multiInd"
	echo "Individuals: ${IND}"
	echo "Population: $POP"
	echo "Scaffold: ${SCAFFOLD}"
	echo "Method: ${METHOD}"
	echo "MSMC input file: ${MSMC_INPUT}"

#if [ -f $VCF1 ]
#	then
#		echo "VCF1 exists, moving on!"
#	else
#		echo "VCF1 does not exist, stopping now"
#		exit 1
#fi

#if [ -f $VCF2 ]
#	then
#		echo "VCF2 exists, moving on !"
#	else
#		echo "VCF2 does not exist, stopping now"
#		exit 1
#fi

	for ind in $IND
		do INDMASK=`ls ${BASEDIR}/mask/ind_mask.lcal.${ind}.${SCAFFOLD}.samtools.bed.gz`
		echo "--mask=$INDMASK " >> mask/${POP}.mask_file.lcal.$SCAFFOLD
		INDVCF=`ls ${BASEDIR}/vcf/phased.${ind}.${SCAFFOLD}.lcal.samtools.vcf.gz`
		echo $INDVCF >> vcf/${POP}.vcf_file.${SCAFFOLD}.lcal
	done

### Generate MSMC input files:
	if [ $METHOD == samtools ]
   		then
   		#MASK_INDIV=`ls ${BASEDIR}/mask/ind_mask.${IND}.${SCAFFOLD}.bed.gz` # store indiv.mask file path
   		MASK_GENOME=${BASEDIR}/mask/lcalcarifer_${SCAFFOLD}.mask.35.50.bed.gz

   		echo "MAPPABILITY MASK: ${MASK_GENOME}"
   		echo "Creating MSMC input file WITH individual mask (samtools)"
   		#${MSMCTOOLS}/generate_multihetsep.py --negative_mask=$MASK_REPEATS --mask=$MASK_INDIV $VCF > $MSMC_INPUT # with repeat mask
   		generate_multihetsep.py `cat mask/${POP}.mask_file.lcal.${SCAFFOLD}` --mask=$MASK_GENOME `cat vcf/${POP}.vcf_file.${SCAFFOLD}.lcal` > ${MSMC_INPUT} # without repeat mask

	elif [ $METHOD == gatk ]
		then
		echo "Creating MSMC input file WITHOUT individual mask (gatk)"
		MASK_GENOME=`ls /project/WagnerLab/jrick/msmc_Sept2017/mask/lcalcarifer_v3_ASB_${SCAFFOLD}.mask.35.50.bed.gz`
		#msmc-tools/generate_multihetsep.py --negative_mask=$MASK_REPEATS $VCF > $MSMC_INPUT # with repeat mask
		${MSMCTOOLS}/generate_multihetsep.py --mask=$MASK_GENOME $VCF2 > $MSMC_INPUT # without repeat mask
	fi

#echo "Editing scaffold names in input files.." # Remove .1 from scaffold names, msmc program won't accept this
#sed s/\\.1//g $MSMC_INPUT > $MSMC_INPUT

done

echo "Done with script."
date
exit 0

####
