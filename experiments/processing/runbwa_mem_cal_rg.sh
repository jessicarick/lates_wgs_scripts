#!/bin/sh

#SBATCH --account=latesgenomics
#SBATCH --time=5-00:00:00
#SBATCH --ntasks-per-node=16
#SBATCH --nodes=1
#SBATCH --mem=120G

date
module load gcc bwa samtools

fastq1=$1
fastq2=$2
info=$3
genome=/project/latesgenomics/reference_genomes/Lates_calcarifer/lcalcarifer_genome_v3_chromosomal.fa

id=$(basename $fastq1 .trimmomatic_1P.fastq)

echo "Mapping reads for $id using info file $info. fastq files are $fastq1 and $fastq2"

# pulling info from info file-- sensitive to order of info in file
ind=$(cat $info | grep $id | head -n 1 | cut -f 10 -d',')
lib=$(cat $info | grep $id | head -n 1 | cut -f 3 -d',')
bcode=$(cat $info | grep $id | head -n 1 | cut -f 1 -d',')
echo "$ind $lib $bcode"
rg=$(echo "@RG\tID:1\tLB:$lib\tPL:illumina\tSM:$ind\tPU:$bcode")
echo "$rg"

bwa mem -t 16 -P -R $rg $genome $fastq1 $fastq2 > /gscratch/jrick/lates_wgs_bamfiles/aln_${id}.trimmomatic.sam

echo "Converting sam to bam for $id"
samtools view -b -S -@ 16 -o /gscratch/jrick/lates_wgs_bamfiles/aln_${id}.trimmomatic.bam /gscratch/jrick/lates_wgs_bamfiles/aln_${id}.trimmomatic.sam
rm -f /gscratch/jrick/lates_wgs_bamfiles/aln_${id}.sam

echo "Sorting and indexing bam files for $id"
samtools sort -@ 16 /gscratch/jrick/lates_wgs_bamfiles/aln_${id}.trimmomatic.bam -o /gscratch/jrick/lates_wgs_bamfiles/aln_${id}.trimmomatic.sorted.bam
#samtools index -@ 16 /project/latesgenomics/jrick/latesWGS_2018/march2020/bamfiles/new_041322/aln_${id}.trimmomatic.sorted.bam

cp /gscratch/jrick/lates_wgs_bamfiles/aln_${id}.trimmomatic.sorted.bam* /project/latesgenomics/jrick/latesWGS_2018/march2020/bamfiles/new_060822/

echo "Done creating sorted bam for $id"
date
