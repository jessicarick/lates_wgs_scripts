#!/bin/sh

#scriptdir=/project/ltcichlidgenomics/cichlid_demography/msmc2/msmc2_scripts
#source ${scriptdir}/msmc_params.sh

#input is either the bam file name or file of bame names (if you are phasing >1 inds)
INPUT=$1
input_type=$2
PHASING=whatshap
SCAFFOLD_FILE=CHROM_only
OUTDIR=/project/latesgenomics/jrick/latesWGS_2018/msmc
method=samtools

#BASE=`basename $BAMFILE .sorted.bam`
#PHASING=$2 #already defined in msmc_params.sh

module load gcc
module load miniconda3
module load samtools
module load bcftools

source activate new_env #needs to have whatshap installed via conda


#pull out the scaffold to work on
scaf=$(awk -v chr=${SLURM_ARRAY_TASK_ID} 'NR == chr {print $1}' $SCAFFOLD_FILE)

#figure out the number of samples being input
if [ $input_type == "file" ]
then
  number_samples=$(cat $INPUT | wc -w)
else
  number_samples=$(echo $INPUT | wc -w)
fi

if [ $number_samples -eq 1 ]
then
  BAMFILE=${BAMDIR}/${INPUT}
  BASE=$(basename $BAMFILE .sorted.bam)
else
  BASE=$(basename $bamlist)  
fi

echo "working with ${OUTDIR}/vcf/${BASE}.${scaf}.${method}.vcf"

#how many samples in the vcf?
number_samples_vcf=$(bcftools query -l ${OUTDIR}/vcf/${BASE}.${scaf}.${method}.vcf | wc -l)
if [ "$number_samples" -ne "$number_samples_vcf" ]
then
  echo "The bam input does not have the number of samples as the vcf. Exiting script."
  exit 1
fi

echo "starting phasing"
date

if [ "$PHASING" == "whatshap" ]
then
  
  if [ $number_samples -eq 1 ]
  then
    echo "phasing one sample (${BASE}) for the following chromosome: $scaf"
    if [ -f $OUTDIR/vcf/${BASE}.${scaf}.${prefix}.minDP5.${PHASING}.${method}.vcf.gz ]; then
      echo "phased VCF already exists; moving onto next scaffold"
    else
      echo "phased VCF does not exist; phasing VCF for scaffold ${scaf}"
      sed -i 's/^ //g' ${OUTDIR}/vcf/${BASE}.${scaf}.${method}.vcf
      
      whatshap phase --reference $GENOME --ignore-read-groups -o ${OUTDIR}/vcf/${BASE}.${scaf}.${prefix}.minDP5.${PHASING}.${method}.vcf.gz ${OUTDIR}/vcf/${BASE}.${scaf}.${method}.vcf $BAMFILE
      whatshap stats --tsv=$OUTDIR/stats/${BASE}.${scaf}.${prefix}.minDP5.${PHASING}.stats.tsv $OUTDIR/vcf/${BASE}.${scaf}.${prefix}.minDP5.${PHASING}.${method}.vcf.gz
    fi
  else
  #FOR MULTIPLE SAMPLES
    echo "phasing $number_samples samples for the following chromosome: $scaf"
    bam_list=$(cat $INPUT | tr '\n' ' ')
    
    if [ -f $OUTDIR/vcf/${BASE}.${scaf}.${prefix}.minDP5.${PHASING}.${method}.vcf.gz ]; then
      echo "phased VCF already exists; moving onto next scaffold"
    else
      echo "phased VCF does not exist; phasing VCF for scaffold ${scaf}"
      sed -i 's/^ //g' ${OUTDIR}/vcf/${BASE}.${scaf}.${method}.vcf
    
      whatshap phase --reference $GENOME -o ${OUTDIR}/vcf/${BASE}.${scaf}.${prefix}.minDP5.${PHASING}.${method}.vcf.gz ${OUTDIR}/vcf/${BASE}.${scaf}.${method}.vcf $bam_list
      whatshap stats --tsv=$OUTDIR/stats/${BASE}.${scaf}.${prefix}.minDP5.${PHASING}.stats.tsv $OUTDIR/vcf/${BASE}.${scaf}.${prefix}.minDP5.${PHASING}.${method}.vcf.gz
    fi
  fi
else
  echo "Phasing parameter is not whatshap; exiting now"
  exit 1
fi

echo "finished phasing"
date
