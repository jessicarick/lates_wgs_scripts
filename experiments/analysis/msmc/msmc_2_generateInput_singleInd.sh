#!/bin/bash -l

### Variables:
#VCF=$1
PHASING=$1
IND=$2
#SCAFFOLD=$4

#IFS='.' read -ra SUBSTRINGS <<< "$VCF"
#IFS='/' read -ra SUBSTRINGS2 <<< "${SUBSTRINGS[0]}"

#TOT=`echo ${#SUBSTRINGS2[@]}`
#IND_NR=`expr $TOT - 1`
#IND=${SUBSTRINGS2[$IND_NR]}

#SCAFFOLD=${SUBSTRINGS[1]}
#PHASING=${SUBSTRINGS[3]}
#PHASING=unphased
METHOD=samtools
#METHOD=${SUBSTRINGS[4]}

BASEDIR=/project/latesgenomics/jrick/latesWGS_2018/msmc

module load gcc
module load samtools
module load python/3.6.3
PATH=$PATH:/project/passerinagenome/scripts/msmc-tools-master

module list

for s in `cat CHROM_only`
	do echo "working on scaffold $s"
	SCAFFOLD=$s
	VCF=`ls /project/latesgenomics/jrick/latesWGS_2018/msmc/vcf/${IND}.${SCAFFOLD}.lcal.samtools.vcf.gz`

	#MASK_REPEATS=repeats.bed.gz # Needs to be gzipped
	MSMC_INPUT=${BASEDIR}/input/msmc_input.${IND}.${SCAFFOLD}.txt
	MSMCTOOLS=/project/passerinagenome/scripts/msmc-tools-master

	printf "\n \n \n \n"
	date
	echo "Script: msmc_2_generateInput_singleInd"
	echo "Individual: ${IND}"
	echo "Scaffold: ${SCAFFOLD}"
	echo "Phasing: ${PHASING}"
	echo "Method: ${METHOD}"
	echo "MSMC input file: ${MSMC_INPUT}"
	echo "VCF: ${VCF}"

	if [ -f $VCF ]
		then
			echo "VCF exists, moving on!"

### Generate MSMC input files:
		if [ $METHOD == samtools ]
   			then
   				MASK_INDIV=${BASEDIR}/mask/ind_mask.lcal.${IND}.${SCAFFOLD}.samtools.bed.gz # store indiv.mask file path
		   		MASK_GENOME=${BASEDIR}/mask/lcalcarifer_${SCAFFOLD}.mask.35.50.bed.gz
   				echo "MASK: ${MASK_INDIV}"
				echo "MAPPABILITY MASK: ${MASK_GENOME}"
				echo "Creating MSMC input file WITH individual mask (samtools)"
				#${MSMCTOOLS}/generate_multihetsep.py --negative_mask=$MASK_REPEATS --mask=$MASK_INDIV $VCF > $MSMC_INPUT # with repeat mask
   				generate_multihetsep.py --mask=${MASK_INDIV} --mask=${MASK_GENOME} $VCF > $MSMC_INPUT # without repeat mask

		elif [ $METHOD == gatk ]
			then
				echo "Creating MSMC input file WITHOUT individual mask (gatk)"
				#msmc-tools/generate_multihetsep.py --negative_mask=$MASK_REPEATS $VCF > $MSMC_INPUT # with repeat mask
				${MSMCTOOLS}/generate_multihetsep.py $VCF > $MSMC_INPUT # without repeat mask
		fi

	else
		echo "VCF does not exist, moving on to next scaffold"
	fi

#echo "Editing scaffold names in input files.." # Remove .1 from scaffold names, msmc program won't accept this
#sed s/\\.1//g $MSMC_INPUT > $MSMC_INPUT
	echo "Done with $SCAFFOLD; moving on to next scaffold"
done;

echo "Done with script."
date
exit 0

####
