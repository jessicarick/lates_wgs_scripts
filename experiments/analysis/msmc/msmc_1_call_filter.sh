#!/bin/sh -l

#SBATCH --account=latesgenomics 
#SBATCH --job-name=msmc1_aln
#SBATCH --mail-type=ALL 
# #SBATCH --output=outs/stdout_%a_msmc1 
# #SBATCH --error=outs/stderr_%a_msmc1
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --time=1-00:00:00 
# #SBATCH --array=0-13

#readarray -t ind_arr < INDS.txt

## VARIABLES:
PHASING=unphased # PHASING=unphased
#SCAFFOLD=$2 # SCAFFOLD=NT_168080.1 #.1 added .1
#IND=${ind_arr[$SLURM_ARRAY_TASK_ID]} # IND=Ceja262
IND=$1

GENOME=/project/latesgenomics/reference_genomes/Lates_calcarifer/lcalcarifer_genome_v3_chromosomal.fa # reference genome fasta
BAMFILE=/project/latesgenomics/jrick/latesWGS_2018/march2020/bamfiles/aln_${IND}.sorted.bam # focal bamfile
MSMCTOOLS=/project/passerinagenome/scripts/msmc-tools-master # folder with msmc-tools binaries

module load gcc
module load samtools/1.6
module load bcftools
module load gatk
module load python/3.6.3
PATH=$PATH:/project/passerinagenome/scripts/msmc-tools-master/

printf "\n \n \n \n"
date
echo "Current script: msmc_1_call_filter.sh"
echo "Individual: $IND"
#echo "Chromosome: ${SCAFFOLD}"
echo "Bamfile: $BAMFILE"
echo "Loaded modules:" 
module list

#if [ -f $GENOME ]
#	then
#		echo "Genome does not need to be unzipped, moving on!"
#	else
#		echo "Genome needs to be unzipped, unzipping now"
#		gunzip -c ${GENOME}.gz > ${GENOME}
#fi

if [ -f "${BAMFILE}.bai" ]
	then
		echo "Bamfile index already exists, moving on!"
	else
		echo "Bamfile index does not exist, creating index"
		samtools index $BAMFILE ${BAMFILE}.bai
fi

export GENOME
export BAMFILE
export MSMCTOOLS

parallel -j 16 --dry-run --delay 2 --env GENOME --env BAMFILE --env MSMCTOOLS 'MEANCOV=$(samtools depth -r $s $BAMFILE | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'); echo "starting samtools alignment"; samtools mpileup -u -r {} -q 20 -Q 20 -C 50 $GENOME $BAMFILE | bcftools call -c -V indels | python3 $MSMCTOOLS/bamCaller_min10.py `echo "$MEANCOV"` /project/latesgenomics/jrick/latesWGS_2018/msmc/mask/ind_mask.lcal.minDP10.${IND}.${s}.samtools.bed.gz > /project/latesgenomics/jrick/latesWGS_2018/msmc/vcf/${IND}.${s}.lcal.minDP10.samtools.vcf' ::: `cat CHROM_only`

### Report:
echo "Filtered VCF and mask created for all scaffolds for $IND ."
echo "Done with script."
date

