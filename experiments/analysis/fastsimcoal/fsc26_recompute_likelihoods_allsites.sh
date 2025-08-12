#!/bin/sh

### Job Name
#SBATCH --job-name fsc_recompute

### Declare an account for the job to run under - You might need to change this to Carling Lab (check designation with arccquota on Moran)
#SBATCH --account=latesgenomics

### The directive below directs that the standard output and error streams are to be merged, \
### intermixed, as standard output. 
### mailing options
#SBATCH --mail-type=END

### Specify Resources
### 1 nodes, 1 processors (cores) each node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
# #SBATCH --mem=64G
#SBATCH --array=0-30

### Set max walltime (days-hours:minutes:seconds)
#SBATCH --time=0:30:00

module load gcc
module load htslib

#spp=$1

###############################
### SETTING UP THE ANALYSES ###
###############################

#PUT THE NAME(S) OF THE MODEL(S) THAT YOU WANT SUMMARIZED IN THE ARRAY (without quotes, separated by spaces):
#currently implemented models: mod_list=(geneflow  early_geneflow recent_geneflow no_geneflow)
#mod_list=(no_geneflow_threepop earlyadmix_threepop earlygeneflow_threepop lateadmix_threepop lategeneflow_threepop singleadmix_pop1_threepop singleadmix_pop2_threepop singlegeneflow_pop1_threepop singlegeneflow_pop2_threepop alladmix_threepop)
#mod_list=(`cat models_allsites.txt | tr '\n' ' '`)
#mod_list=(`cat mod_uni`)
mod_list=(fbranch_model_new_new)

#path to directory where the model directories are located
file_path=/project/latesgenomics/jrick/admixture/july2020/fastsimcoal/six_taxon_models_noContam_jun2022

#location to put directory that will store the likelihood files
storage_path=/project/latesgenomics/jrick/admixture/july2020/fastsimcoal/six_taxon_models_noContam_jun2022

#number of rounds of simulations
#100 is a good number of iterations
number_of_iters=200


#########################
### RUNNING  ANALYSES ###
#########################

#make directory to store the likelihood files
mkdir $storage_path/mods_recompute_lhoods_latesWGS_2_082922

#loop through each model
#for mod in ${mod_list[*]}
mod=${mod_list[$SLURM_ARRAY_TASK_ID]}
echo $mod
#do
	cd $file_path/$mod/bestrun
	
	# create temporary obs file with name _maxL_MSFS.obs
	cp ${mod}_MSFS.obs ${mod}_maxL_MSFS.obs
	
	# Run fastsimcoal $number_of_iters times to get the likelihood of the observed SFS under the best parameter values with 1 million simulated SFS.
	for iter in $(eval echo "{1..$number_of_iters}")
	do
		fsc26 -i ${mod}_maxL.par -n1000000 -m -q -C 10 -L 100 -B 4 -c 4 --multiSFS
		# Fastsimcoal will generate a new folder called ${model}_maxL and write files in there
		
		# collect the lhood values (Note that >> appends to the file, whereas > would overwrite it)
		sed -n '2,3p' ${mod}_maxL/${mod}_maxL.lhoods >> ${mod}.lhoods
		
		# delete the folder with results
		rm -r ${mod}_maxL/
	done
	
	#move likelihood file to the final storage directory
	mv ${mod}.lhoods $storage_path/mods_recompute_lhoods_latesWGS_2_082922/
#done
