#!/bin/sh

#SBATCH --account=latesgenomics
#SBATCH --time=6-00:00:00
#SBATCH --job-name=allsites
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --array=0-23

module load gcc
module load samtools
module load bcftools

date

#chr=$1

readarray -t chr < chrom

chrom=${chr[$SLURM_ARRAY_TASK_ID]}

echo "creating allsites VCF for $chrom"

bcftools mpileup -f /project/latesgenomics/reference_genomes/Lates_calcarifer/lcalcarifer_genome_v3_chromosomal.fa -b bamfiles.txt -a FORMAT/DP,FORMAT/AD -r $chrom -Ou --max-depth 1000 | bcftools call -m -Oz -f GQ,GP -o latesWGS_july2020_allsites.nocontam.${chrom}.vcf.gz

date
