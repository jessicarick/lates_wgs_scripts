#!/bin/sh

#scriptdir=/project/ltcichlidgenomics/cichlid_demography/msmc2/msmc2_scripts
#source ${scriptdir}/msmc_params.sh

bam_list=$1 #BAM list
BASE=$(basename $bam_list)
PHASING="whatshap"
SCAFFOLD_FILE=CHROM_only

#pull out the scaffold to work on
scaf=$(awk -v chr=${SLURM_ARRAY_TASK_ID} 'NR == chr {print $1}' $SCAFFOLD_FILE)

bam_count=$(wc -l < $bam_list)
slurm_throttle=$((bam_count))


echo "starting phasing"
date

if [ "$PHASING" == "whatshap" ]
then
sbatch --account=latesgenomics \
--job-name=msmc2_1.5 \
--mail-type=END,FAIL \
--mail-user=jrick@uwyo.edu \
--output=outs/std_msmc2_1.5_singlesamp_%A_%a.out_err \
--error=outs/std_msmc2_1.5_singlesamp_%A_%a.out_err \
--mem=50G \
--nodes=1 \
--ntasks-per-node=1 \
--time=0-02:00:00 \
--array=1-$bam_count%3 \
msmc_1.5_phase_singlesamp.sh $bam_list $scaf $BASE

echo "finished launching phasing for all samples for the following chromosome: $scaf"
date
else

echo "phasing is not whatshap --> exiting script"
fi
