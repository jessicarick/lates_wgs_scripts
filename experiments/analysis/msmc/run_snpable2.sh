#!/bin/sh

#SBATCH --job-name snpable
#SBATCH --account=latesgenomics
#SBATCH -o stdout_snpable2_cal_%A
#SBATCH -e stderr_snpable2_cal_%A
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jrick@uwyo.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=2-00:00:00
#SBATCH --mem=0

date

echo "loading modules"

module load gcc
module load bwa
module load samtools
PATH=$PATH:/project/passerinagenome/jrick/msmc/snpable/scripts

echo "modules loaded, read to go!"

ref='/project/latesgenomics/reference_genomes/Lates_calcarifer/lcalcarifer_genome_v3_chromosomal.fa'
k=35
out=lcalcarifer
basedir=/project/latesgenomics/jrick/latesWGS_2018/msmc

mkdir ${basedir}/snpable
cd ${basedir}/snpable

echo "Starting extraction of overlapping ${k}-mer subsequences"

splitfa $ref $k | split -l 20000000 
cat x* >> ${out}_split.$k
rm -f x*

# if it can't find splitfa, try adding seqbility to the path using 'PATH=$PATH:/project/WagnerLab/jrick/msmc_Sept2017/snpable/scripts'

echo "Aligning ${k}-mer reads to the genome with BWA, then converting to sam file"

# the genome needs to be indexed prior to this step-- if it has not already been indexed, run:
#echo "indexing $ref"
#bwa index $ref

echo "aligning reads to genome with BWA and converting to sam"
bwa aln -t 8 -R 1000000 -O 3 -E 3 ${ref} ${out}_split.${k} > ${out}_split.${k}.sai
bwa samse -f ${out}_split.${k}.sam $ref ${out}_split.${k}.sai ${out}_split.${k} 

echo "reads aligned, starting to generate rawMask"
gen_raw_mask.pl ${out}_split.${k}.sam > ${out}_rawMask.${k}.fa

echo "raw mask created as ${out}_rawMask.35.fa, now generating final mask with stringency 50%"
gen_mask -l ${k} -r 0.5 ${out}_rawMask.${k}.fa > ${out}_mask.${k}.50.fa

echo "all done! final mask saved as ${out}_mask.${k}.50.fa"
echo "done!"
date
