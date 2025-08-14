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
module load samtools/1.14
module load python/3.10.6
module load parallel
PATH=$PATH:/project/passerinagenome/scripts/msmc-tools-master

module list

export IND
export BASEDIR

parallel -j 8 --delay 2 --env BASEDIR --env IND "echo 'working on scaffold {}' && generate_multihetsep.py --mask=${BASEDIR}/mask/ind_mask.latesWGS_2.bamlist$IND.{}.minDP5.whatshap.samtools.bed.gz --mask=${BASEDIR}/mask/lcalcarifer_{}.mask.35.50.bed.gz /gscratch/jrick/lates_wgs_vcf/latesWGS_2.bamlist.${IND}.{}.minDP5.whatshap.samtools.vcf.gz > ${BASEDIR}/input/msmc_input.latesWGS_2.bamlist.${IND}.{}.txt" ::: `cat CHROM_only`

echo "Done with script."
date
exit 0

####
