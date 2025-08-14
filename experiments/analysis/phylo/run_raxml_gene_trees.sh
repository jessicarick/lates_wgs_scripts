#!/bin/sh

#SBATCH --account=latesgenomics
#SBATCH --time=0-08:00:00
#SBATCH --ntasks-per-node=16
#SBATCH --nodes=1
#SBATCH --mail-type=END
#SBATCH --array=1-23

date

module load gcc openmpi raxml miniconda3 parallel
module list

readarray -t chrom_list < chrom

chr=${chrom_list[$SLURM_ARRAY_TASK_ID]}

cd partition_trees/
ls ../partition_phylips/*${chr}_*.noInv.phy | sed 's/\.noInv\.phy//g' | parallel -j 8 "echo 'working with {}' && raxmlHPC-PTHREADS -s {}.noInv.phy -m ASC_GTRCAT -n {/} -T 2 -V -x 12345 -m ASC_GTRCAT --asc-corr=lewis -p 12345 -# 100 -f a"

echo "done"
date
