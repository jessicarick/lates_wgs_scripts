#!/bin/sh

vcf=$1
chr=$2
nbp=$3

module load vcftools

echo "working on chrom $chr and vcf $vcf"

if [ ${vcf: -7} == ".vcf.gz" ]; then
		maxbp=`zcat $vcf | tail -n 1 | cut -f 2`
	else
		maxbp=`tail -n 1 $vcf | cut -f 2`
fi

start_bp=1
end_bp=$nbp
n_part=1

while [ "$start_bp" -lt "$maxbp" ]; do 
    echo "creating partition $n"
		vcftools --vcf $vcf --chr $chr --from-bp $start_bp --to-bp $end_bp --recode --out partition_vcfs/latesWGS.${chr}_${n_part}
		start_bp=$((start_bp + nbp))
		end_bp=$((end_bp + nbp))
		n_part=$((n_part + 1))
done
