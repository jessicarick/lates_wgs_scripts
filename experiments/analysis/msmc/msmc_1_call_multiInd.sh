#!/bin/sh

# all parameters come from the msmc_param control file
# make edits there before using this script!
#scriptdir=$(dirname "$0") # note: doesn't work with slurm wrapper!
#source ${scriptdir}/msmc_params.sh

## VARIABLES:
#IND=$1
#BAMFILE=${BAMDIR}/aln_${IND}.sorted.bam
BAMLIST=$1
OUTDIR=/project/latesgenomics/jrick/admixture/july2020/msmc/msmc
METHOD=samtools
GENOME=/project/latesgenomics/reference_genomes/Lates_calcarifer/lcalcarifer_genome_v3_chromosomal.fa

module load gcc samtools bcftools

printf "\n \n \n \n"
date
echo "Current script: msmc_1_call_multiInd.sh"
#echo "Individual: $IND"
echo "Bamlist: $BAMLIST"

#if [ -f "${BAMFILE}.bai" ]
#        then
#                echo "Bamfile index already exists, moving on!"
#        else
#                echo "Bamfile index does not exist, creating index"
#                samtools index $BAMFILE ${BAMFILE}.bai
#fi

readarray -t scafs < CHROM_only
s=${scafs[$SLURM_ARRAY_TASK_ID]}

#for s in `cat CHROM_only`; \
#        do echo "Handling scaffold $s"; \
echo "Handling scaffold $s"
#        ### Calculate mean coverage (to be used as input for bamCaller.py):
#        MEANCOV=`samtools depth -r $s $BAMFILE | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` # calculate mean coverage
#        echo ${IND}.${s} $MEANCOV >> ${OUTDIR}/coverage_samtoolsDepth_${IND}.txt # save mean coverage in separate file
#        echo "Mean coverage for this individual, scaffold ${s}: $MEANCOV"

#        ### Generate a single-sample VCF and a mask-file:
#        MASK_IND=${OUTDIR}/mask/ind_mask.${IND}.${s}.${METHOD}.bed.gz # Individual mask file to be created
        VCF=${OUTDIR}/vcf/${BAMLIST}.${s}.${METHOD}.vcf # VCF file to be created

        #If genome isn't indexed, add:
        #samtools faidx $GENOME
        if [ "$METHOD" == "samtools" ]; then
                echo "starting samtools alignment"
		# filters for quality > 20, no depth filter currently
                bcftools mpileup -Ou -q 20 -Q 20 -A -r ${s} --threads 16 -f $GENOME -a INFO/AD --bam-list $BAMLIST | bcftools call -c --threads 16 -V indels | bcftools view -Ov > ${VCF}

        fi
        echo "done with scaffold ${s}"
#done

### Report:
echo "Filtered VCF created for scaffold $s for ${BAMLIST}."
echo "Done with script."
date

