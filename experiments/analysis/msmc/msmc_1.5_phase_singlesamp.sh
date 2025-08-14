#!/bin/sh

#####################
### SCRIPT SET-UP ###
##################### 
#load variables and modules in msmc_params.sh
#scriptdir=/project/ltcichlidgenomics/cichlid_demography/msmc2/msmc2_scripts
#source ${scriptdir}/msmc_params.sh

#activate conda environment
module load gcc
module load miniconda3
module load bcftools/1.14
#source activate new_env #needs to have whatshap installed via conda
source activate /project/latesgenomics/jrick/latesWGS_2018/msmc/conda_env/msmc_env

#assign input variables to informative names
bam_list=$1
scaf=$2
BASE=$3
OUTDIR=/project/latesgenomics/jrick/latesWGS_2018/msmc
GENOME=/project/latesgenomics/reference_genomes/Lates_calcarifer/lcalcarifer_genome_v3_chromosomal.fa
TEMPDIR=/gscratch/jrick/lates_wgs_vcf/new_042723 # directory for files to be written -- I use gscratch because of space issues
PHASING="whatshap"
method=samtools

##########################
### INITIAL PROCESSING ###
##########################
#pull out the bam file to work on and identify sample name from bam file
bam_samp=$(awk -v samp=${SLURM_ARRAY_TASK_ID} 'NR == samp {print $1}' $bam_list)
sample_name=$(basename ${bam_samp} .sorted.bam | sed 's/\.trimmomatic//g' | sed 's/aln_//g') ## WARNING: SENSITIVE TO FILE NAMING CONVENTION -- make sure this matches your bamfile names

## check if phased output exists and is not empty
if [ -s /gscratch/jrick/lates_wgs_vcf/${BASE}.${sample_name}.${scaf}.minDP5.${PHASING}.${method}.vcf.gz ]; then
	echo "file ${BASE}.${sample_name}.${scaf}.minDP5.${PHASING}.${method}.vcf.gz already exists; not creating a new one"
	exit 1
else

#subset vcf down to the focal sample
bcftools view -s $sample_name ${OUTDIR}/vcf/${BASE}.${scaf}.${method}.vcf > ${TEMPDIR}/TEMP_SUBSET_${sample_name}_${scaf}.vcf


####################
### PHASE SAMPLE ###
####################
echo "starting phasing ${scaf} for ${sample_name}"
echo "bam fed to WhatsHap: $bam_samp"
echo "vcf fed to WhatsHap as the following sample $(bcftools query -l ${TEMPDIR}/TEMP_SUBSET_${sample_name}_${scaf}.vcf)" 
date

whatshap phase --reference $GENOME -o ${TEMPDIR}/${BASE}.${sample_name}.${scaf}.minDP5.${PHASING}.${method}.vcf.gz ${TEMPDIR}/TEMP_SUBSET_${sample_name}_${scaf}.vcf $bam_samp
tabix ${TEMPDIR}/${BASE}.${sample_name}.${scaf}.minDP5.${PHASING}.${method}.vcf.gz
whatshap stats --tsv=${TEMPDIR}/stats/${BASE}.${sample_name}.${scaf}.minDP5.${PHASING}.stats.tsv ${TEMPDIR}/${BASE}.${sample_name}.${scaf}.minDP5.${PHASING}.${method}.vcf.gz

#remove the temporary subset vcf after phasing
#cp ${TEMPDIR}/${BASE}.${sample_name}.${scaf}.minDP5.${PHASING}.${method}.vcf.gz $OUTDIR/vcf/
#rm -f ${TEMPDIR}/TEMP_SUBSET_${sample_name}_${scaf}.vcf

#whatshap phase --reference $GENOME -o ${OUTDIR}/vcf/${BASE}.${sample_name}.${scaf}.${prefix}.minDP5.${PHASING}.${method}.vcf.gz ${OUTDIR}/vcf/${BASE}.${scaf}.${method}.vcf $bam_samp
#whatshap stats --tsv=$OUTDIR/stats/${BASE}.${sample_name}.${scaf}.${prefix}.minDP5.${PHASING}.stats.tsv ${OUTDIR}/vcf/${BASE}.${sample_name}.${scaf}.${prefix}.minDP5.${PHASING}.${method}.vcf.gz

echo "finished phasing"

fi
date
