#!/bin/sh

### Job Name
#SBATCH --job-name single_bs

### Declare an account for the job to run under - You might need to change this to Carling Lab (check designation with arccquota on Moran)
#SBATCH --account=latesgenomics

### Specify Resources
### 1 nodes, 1 processors (cores) each node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8

### Set max walltime (days-hours:minutes:seconds)
#SBATCH --time=5-00:00:00
#SBATCH --array=1-100

module load gcc
module load htslib

bs_num=$SLURM_ARRAY_TASK_ID
file_path=$1
storage_folder=$2
number_of_runs=$3
PREFIX=$4
number_of_sims=$5

remove_run_dirs="YES" # change if you don't want individual run directories deleted at the end

#enter $bs bootstrap folder
cd $file_path/$storage_folder/bs$bs_num
  
#complete $number_of_runs runs for the bootstrap
restart=FALSE
if [ -d "run1/" ]; then
	restart=TRUE
fi

if [ $restart == TRUE ]; then
	last=`ls -d run* | sed 's/run//g' | sort -n | tail -n 1`
	if [ -f "bs${last}/*/*_maxL.txt" ]; then
		start=$(($last + 1))
		echo "restarting with run $start"
	else
		start=$last
		echo "restarting with run $start"
	fi
else
	start=1
fi

for i in `seq $start $number_of_runs`
do
  mkdir run$i #create folder for run i
  #copy and rename tpl and est files into the run folder
  cp $file_path/$storage_folder/${PREFIX}.tpl run$i/${PREFIX}.bs.${bs_num}.tpl
  cp $file_path/$storage_folder/${PREFIX}.est run$i/${PREFIX}.bs.${bs_num}.est
  cp $file_path/$storage_folder/bs${bs_num}/${PREFIX}.bs.${bs_num}_MSFS.obs run$i"/"

  cd run$i #enter the run folder
        
  #run faststimcoal
  fsc26 -t ${PREFIX}.bs.${bs_num}.tpl -e ${PREFIX}.bs.${bs_num}.est -m -C 10 -n $number_of_sims -L 50 -s0 -M -q -c 8 -B 8 --multiSFS --initvalues $file_path/$storage_folder/bestrun/${PREFIX}.pv
  cd .. #leave the run folder (go up 1 level in the hierarchy)
        
done
 
/home/jrick/bin/fastsimcoal_scripts/scripts/bootstrap_selectbestrun.sh

echo 'Number of runs: ' > run_confirmation.txt
find . -type d -name 'run*' | wc -l >> run_confirmation.txt
n_completed=`find . -type d -name 'run*' | wc -l`

#Delete all the run directories after the best run directory is identified and extracted (these are never used again)
#helpful source: https://unix.stackexchange.com/questions/23576/how-do-i-recursively-delete-directories-with-wildcard
if [ $remove_run_dirs = 'YES' & -d "bestrun/" & "$n_completed" -eq "$number_of_runs" ]
then
     find . -type d -name 'run*' -exec rm -r {} +
else
     echo "did not remove run directories"
fi

echo "done with all bootstraps for bs $bs_num"
