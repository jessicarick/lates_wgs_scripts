#!/bin/sh

IND=$1
GENOME=/project/latesgenomics/reference_genomes/Lates_calcarifer/lcalcarifer_genome_v3_chromosomal.fa
OUTDIR=/project/latesgenomics/jrick/admixture/july2020/msmc/msmc
BAMFILE=/project/latesgenomics/jrick/latesWGS_2018/march2020/bamfiles/aln_${IND}.sorted.bam

module load gcc
module load miniconda3

source activate new_env # needs to have whatshap installed via conda

echo "working with individual $IND"

for s in `cat CHROM_only`
	do echo "working with scaffold $s"
	if [ -f $OUTDIR/vcf/${IND}.${s}.lcal.minDP10.whatshap.samtools.vcf.gz ]; then
		echo "phased VCF already exists; moving onto next scaffold"
	else
		echo "phased VCF does not exist; phasing VCF for scaffold $s"
		sed -i 's/^ //g' $OUTDIR/vcf/${IND}.${s}.lcal.minDP10.samtools.vcf

		whatshap phase --reference $GENOME --ignore-read-groups -o ${OUTDIR}/vcf/${IND}.${s}.lcal.minDP10.whatshap.samtools.vcf.gz ${OUTDIR}/vcf/${IND}.${s}.lcal.minDP10.samtools.vcf $BAMFILE
		whatshap stats --tsv=$OUTDIR/stats/${IND}.${s}.lcal.minDP10.whatshap.stats.tsv $OUTDIR/vcf/${IND}.${s}.lcal.minDP10.whatshap.samtools.vcf.gz 
	fi
done

echo "finished with individual $IND"
date
