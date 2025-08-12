#!/bin/sh

### Job Name
#SBATCH --job-name fsc

### Declare an account for the job to run under - You might need to change this to Carling Lab (check designation with arccquota on Moran)
#SBATCH --account=latesgenomics

### The directive below directs that the standard output and error streams are to be merged, \
### intermixed, as standard output. 
### mailing options
#SBATCH --mail-type=END

### Specify Resources
### 1 nodes, 8 processors (cores) each node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64G
#SBATCH --array=1-300%20

### Set max walltime (days-hours:minutes:seconds)
#SBATCH --time=0-3:00:00

date
module load gcc
module load htslib
module load parallel

###############################
### SETTING UP THE ANALYSES ###
###############################
PREFIX=$1 # prefix defined on command line
MODEL=$2

file_path=/project/latesgenomics/jrick/admixture/july2020/fastsimcoal/six_taxon_models_noContam_jun2022
fsc26=/home/jrick/bin/fsc26 
number_of_runs=300
number_of_sims=200000

i=$SLURM_ARRAY_TASK_ID

echo "running $MODEL , run number $i"

#########################
### RUNNING  ANALYSES ###
#########################

#make directory within the specified directory based on the PREFIX
#mkdir $file_path/$MODEL

#move .tpl, .est, and .obs files into new file ${PREFIX}
#cp $file_path/*_MSFS.obs $file_path/${MODEL}/${MODEL}_MSFS.obs

#enter newly-created directory (this is where modelling results will be located)
cd $file_path/$MODEL

#slightly modified for loop header to allow for a variable-defined number of iterations (i.e. number of runs):
#https://www.cyberciti.biz/faq/unix-linux-iterate-over-a-variable-range-of-numbers-in-bash/

#export file_path
#export PREFIX
#export number_of_sims

#parallel --env file_path --env PREFIX --env number_of_sims "mkdir run{}/ && cp $file_path/${PREFIX}/${PREFIX}.tpl $file_path/${PREFIX}/${PREFIX}.est $file_path/${PREFIX}/${PREFIX}_jointMAFpop1_0.obs run{}/ && cd run{}/ && fsc26 -t ${PREFIX}.tpl -e ${PREFIX}.est -m -0 -C 10 -n $number_of_sims -L 50 -s0 -M -q -c 16" ::: `seq $number_of_runs`

#if [ -d "run1/" ]; then
#	begin=`ls -d run* | sort -V | tail -n 1 | sed 's/run//g'`
#	echo "model previously run; starting again with run $begin"
#else
#	begin=1
#fi

#for i in $(eval echo "{$start..$end}") # make parallel later?
#do
    mkdir run$i
    cp $file_path/${MODEL}/${MODEL}.tpl $file_path/${MODEL}/${MODEL}.est $file_path/${MODEL}/${MODEL}_MSFS.obs run$i"/"
    cd run$i
    if [ -f ${MODEL}_MSFS.obs ]
	then 
    	echo "using multi SFS for running fastsimcoal"
	fsc26 -t ${MODEL}.tpl -e ${MODEL}.est -m -C 10 -n $number_of_sims -L 50 -s0 -M -q -c 8 -B 8 --multiSFS
    elif [ -f ${MODEL}_jointMAFpop1_0.obs ]
	then
	echo "using joint SFS for running fastsimcoal"
	fsc26 -t ${MODEL}.tpl -e ${MODEL}.est -m -C 4 -n $number_of_sims -L 50 -s0 -M -q -c 8 -B 8
	
    else
	echo "could not find file named ${MODEL}_MSFS.obs or ${MODEL}_jointMAFpop1_0.obs! exiting now."
	exit 1
    fi
    cd ..


#identify run with best likelihood
nrun=`ls -d run*/ | wc -l`
echo "number of completed runs = $nrun"
if [ "$nrun" == "$number_of_runs" ]; then 
	echo "number of completed runs equals number of desired runs; selecting best now"
	/project/latesgenomics/jrick/admixture/july2020/fastsimcoal/scripts/fsc_selectbestrun.sh
fi

#if [ -f bestrun/ ]; then
#	rm -rf run*/
#fi

date
