#!/bin/sh

module load gcc
module load python/2.7.18

basedir=/project/latesgenomics/jrick/latesWGS_2018/msmc
msmc_tools_path=/project/ltcichlidgenomics/cichlid_demography/msmc2/msmc-tools
init_vcf_path=${basedir}/vcf
mask_vcf_path=${basedir}/mask_vcf
mask_bed_path=${basedir}/mask

samp_file=${basedir}/latesWGS_2.bamlist
scaf_file=${basedir}/CHROM_only

SAMP_BAM=$(awk -v samp=${SLURM_ARRAY_TASK_ID} 'NR == samp {print $1}' $samp_file)
SAMP=$(basename $SAMP_BAM .trimmomatic.sorted.bam | sed 's/aln_//g')
echo -e "=========================\n=== STARTING $SAMP\n========================="
   
cat $scaf_file | while read SCAF
do
  echo "$SCAF"
#  if [ ! -f /gscratch/jrick/msmc/ind_mask.latesWGS_2.bamlist.$SAMP.$SCAF.minDP5.whatshap.samtools.vcf.gz ]; then
  zcat /gscratch/jrick/lates_wgs_vcf/latesWGS_2.bamlist.$SAMP.$SCAF.minDP5.whatshap.samtools.vcf.gz | python $msmc_tools_path/vcfAllSiteParser.py $SCAF $mask_bed_path/ind_mask.latesWGS_2.bamlist$SAMP.$SCAF.minDP5.whatshap.samtools.bed.gz | gzip -c > /gscratch/jrick/msmc/ind_mask.latesWGS_2.bamlist.$SAMP.$SCAF.minDP5.whatshap.samtools.vcf.gz
  #zcat $init_vcf_path/bamlist.$SAMP.$SCAF.pnye.minDP5.whatshap.samtools.vcf.gz | python $msmc_tools_path/vcfAllSiteParser.py $SCAF $mask_bed_path/ind_mask.$SAMP.$SCAF.pnye.minDP5.whatshap.samtools.bed.gz > $mask_vcf_path/ind_mask.$SAMP.$SCAF.pnye.minDP5.whatshap.samtools.vcf
#  else
#  echo "mask file already exists for $SAMP $SCAF"
#  fi
done



echo "done with $SAMP"

#cat $samp_file | while read SAMP_BAM 
#do
#   SAMP=$(basename $SAMP_BAM .rg.bam)
#   echo -e "=========================\n=== STARTING $SAMP\n========================="
#   
#   cat $scaf_file | while read SCAF
#   do
#     echo "$SCAF"
#     zcat $init_vcf_path/bamlist.$SAMP.$SCAF.pnye.minDP5.whatshap.samtools.vcf.gz | python $msmc_tools_path/vcfAllSiteParser.py $SCAF $mask_bed_path/ind_mask.$SAMP.$SCAF.pnye.minDP5.whatshap.samtools.bed.gz | gzip -c > $mask_vcf_path/ind_mask.$SAMP.$SCAF.pnye.minDP5.whatshap.samtools.vcf.gz
#     #zcat $init_vcf_path/bamlist.$SAMP.$SCAF.pnye.minDP5.whatshap.samtools.vcf.gz | python $msmc_tools_path/vcfAllSiteParser.py $SCAF $mask_bed_path/ind_mask.$SAMP.$SCAF.pnye.minDP5.whatshap.samtools.bed.gz > $mask_vcf_path/ind_mask.$SAMP.$SCAF.pnye.minDP5.whatshap.samtools.vcf
#   done   
#done






#msmc_tools_path=/project/ltcichlidgenomics/cichlid_demography/msmc2/msmc-tools
#init_vcf_path=/project/ltcichlidgenomics/cichlid_demography/msmc2/output/vcf
#output_path=/gscratch/alewansk/prac_msmc_mask

#zcat $init_vcf_path/bamlist.PBM2015_566.chr20.pnye.minDP5.whatshap.samtools.vcf.gz | python $msmc_tools_path/vcfAllSiteParser.py chr20 $output_path/chr3Mask_output_bed.gz > $output_path/chr3_output.vcf

