module load gcc
module load msmc2
module load samtools
module load bwa
module load python/3.6.3

MSMCTOOLS=/project/passerinagenome/scripts/msmc-tools-master # folder with msmc-tools binaries
PATH=$PATH:$MSMCTOOLS

OUTDIR=/project/latesgenomics/jrick/admixture/july2020/msmc # main directory for output files
SCAFFOLDS=SCAFFOLDS.txt

# for msmc_1_call.sh
GENOME=/project/latesgenomics/reference_genomes/Lates_calcarifer/lcalcarifer_genome_v3_chromosomal.fa # reference genome fasta
prefix=lcal # prefix of genome masks
BAMDIR=/project/latesgenomics/jrick/latesWGS_2018/may2020/bamfiles # directory with bamfiles
k=35

# for msmc_3_generateInput.sh
NR_IND=1 # number of individuals in analysis
POP_OR_IND=POPNAME # name of individual or population being analyzed for script 3
DATE=`date +%m%d%y`
RUN_NAME=msmc_${POP_OR_IND}_${DATE}
P_PAR=1*2+25*1+1*2+1*3 

nchr=`wc -l $SCAFFOLDS` # number of chromosomes in reference genome
sex_chr=LG9 # name of sex chromosome to omit in analyses

METHOD="samtools"
PHASING="whatshap"
