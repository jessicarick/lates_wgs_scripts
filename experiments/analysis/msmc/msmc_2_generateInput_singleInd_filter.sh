#!/bin/bash -l

IND=$1
PHASING=$2
METHOD=samtools

for s in `cat CHROM_LM`;
do echo $s ;
SCAFFOLD=$s ;
VCF=`ls vcf/lmariae_MEGA/${IND}.gatk.beagle.${SCAFFOLD}.vcf.gz`;

#IFS='.' read -ra SUBSTRINGS <<< "$VCF"
#IFS='/' read -ra SUBSTRINGS2 <<< "${SUBSTRINGS[0]}"

#TOT=`echo ${#SUBSTRINGS2[@]}`
#IND_NR=`expr $TOT - 1`
#IND=${SUBSTRINGS2[$IND_NR]}

#SCAFFOLD=${SUBSTRINGS[1]}
#PHASING=${SUBSTRINGS[3]}
#PHASING=unphased
#METHOD=${SUBSTRINGS[4]}

#MASK_REPEATS=repeats.bed.gz # Needs to be gzipped

MSMC_INPUT=/project/WagnerLab/jrick/msmc_Sept2017/input/beagle/lmariae_MEGA/msmc_input.${IND}.${SCAFFOLD}.${METHOD}.beagle.phased.txt
MSMCTOOLS=/project/WagnerLab/jrick/msmc-tools-master

module load samtools
module load python/3.4.3
PATH=$PATH:/project/WagnerLab/jrick/msmc-tools-master

printf "\n \n \n \n"
date
echo "Script: msmc_2_generateInput_singleInd_filter"
echo "Individual: ${IND}"
echo "Scaffold: ${SCAFFOLD}"
echo "Phasing: ${PHASING}"
echo "Method: ${METHOD} + beagle"
echo "MSMC input file: ${MSMC_INPUT}"
echo "VCF: ${VCF}"

if [ -f $VCF ]
	then
		echo "VCF exists, moving on!"
	else
		echo "VCF does not exist, stopping now"
		exit 1
fi

#gunzip -c $VCF | sed 's/^ *//g' | gzip -c > $VCF

### Generate MSMC input files:
if [ $METHOD == samtools ]
   then
   MASK_INDIV=`ls /project/WagnerLab/jrick/msmc_Sept2017/mask/ind_mask.${IND}.${SCAFFOLD}.bed.gz` # store indiv.mask file path
   MASK_GENOME=`ls /project/WagnerLab/jrick/msmc_Sept2017/mask/lcalcarifer_v3_ASB_${SCAFFOLD}.mask.35.50.bed.gz`
   echo "MASK: ${MASK_INDIV}"
   echo "MAPPABILITY MASK: ${MASK_GENOME}"
   echo "Creating MSMC input file WITH individual mask (samtools)"
   #${MSMCTOOLS}/generate_multihetsep.py --negative_mask=$MASK_REPEATS --mask=$MASK_INDIV $VCF > $MSMC_INPUT # with repeat mask
   generate_multihetsep.py --mask=$MASK_INDIV --mask=$MASK_GENOME $VCF > $MSMC_INPUT # without repeat mask

elif [ $METHOD == gatk ]
	then
   	#MASK_GENOME=`ls /project/WagnerLab/jrick/msmc_Sept2017/mask/lcalcarifer_v3_ASB_${SCAFFOLD}.mask.35.50.bed.gz`
	echo "Creating MSMC input file WITHOUT mask (gatk)"
	#msmc-tools/generate_multihetsep.py --negative_mask=$MASK_REPEATS $VCF > $MSMC_INPUT # with repeat mask
	generate_multihetsep.py $VCF > $MSMC_INPUT # without repeat mask
fi

#echo "Editing scaffold names in input files.." # Remove .1 from scaffold names, msmc program won't accept this
#sed s/\\.1//g $MSMC_INPUT > $MSMC_INPUT

echo "Done with scaffold ${s}, moving to next scaffold";

done;

echo "Done with script."
date
exit 0

####
