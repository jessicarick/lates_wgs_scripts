#!/bin/bash -l

### Variables:
#VCF=$1
PHASING=$1
IND_LIST=$2
SPP=`basename $IND_LIST .ind`
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
module load samtools/1.6
module load python/3.6.3
module load parallel
PATH=$PATH:/project/passerinagenome/scripts/msmc-tools-master

module list

mask_call=$(for ind in `cat $IND_LIST`; do echo "--mask=${BASEDIR}/mask/ind_mask.latesWGS_2.bamlist${ind}.{}.minDP5.${PHASING}.samtools.bed.gz "; done)
ind_vcfs=$(for ind in `cat $IND_LIST`; do echo "/gscratch/jrick/msmc/ind_mask.latesWGS_2.bamlist.${ind}.{}.minDP5.${PHASING}.samtools.vcf.gz "; done)

echo $mask_call
echo $ind_vcfs

export mask_call
export ind_vcfs
export BASEDIR
export SPP

parallel -j 8 --delay 2 --env BASEDIR --env SPP --env mask_call --env ind_vcfs "echo 'working on scaffold {}' && 
	 generate_multihetsep.py \
		`echo $mask_call` \
		--mask=${BASEDIR}/mask/lcalcarifer_{}.mask.35.50.bed.gz \
		`echo $ind_vcfs` > ${BASEDIR}/input/msmc_input.latesWGS_2.bamlist.${SPP}.{}.txt" ::: `cat CHROM_only`

echo "Done with script."
date
exit 0

####
