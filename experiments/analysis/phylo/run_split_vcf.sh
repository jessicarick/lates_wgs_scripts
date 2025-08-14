#!/bin/sh

#SBATCH --account=latesgenomics
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --array=23

date
module load gcc vcftools
module list

vcf=$1
nbp=$2 # number of bp for windows

readarray -t chr_list < chrom
chr=${chr_list[$SLURM_ARRAY_TASK_ID]}
#vcf=`ls ${vcf_dir}/*.recode.vcf | grep $chr | grep -v beagle`

vcftools --gzvcf $vcf --chr $chr --recode --out ${chr}.tmp

echo "starting split_vcf_coord.sh script for $chr; vcf is $vcf"
./scripts/split_vcf_coord.sh ${chr}.tmp.recode.vcf $chr $nbp

date
