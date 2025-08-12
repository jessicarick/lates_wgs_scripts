#!/bin/sh

### Job Name
#SBATCH --job-name fsc_AIC

### Declare an account for the job to run under
#SBATCH --account=latesgenomics

### The directive below directs that the standard output and error streams are to be merged, \
### intermixed, as standard output. 
### mailing options
#SBATCH --mail-type=END
#SBATCH --mail-user=jrick@uwyo.edu

### Specify Resources
### 1 nodes, 1 processors (cores) each node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
# #SBATCH --mem=0

### Set max walltime (days-hours:minutes:seconds)
#SBATCH --time=0-01:00:00

#load libraries
module load gcc
module load r

#spp=$1

###############################
### SETTING UP THE ANALYSES ###
###############################

#PUT THE NAME(S) OF THE MODEL(S) THAT YOU WANT SUMMARIZED IN THE ARRAY (without quotes, separated by spaces):
#currently implemented models: mod_list=(geneflow  early_geneflow recent_geneflow no_geneflow)
#mod_list=(mod1 mod2 mod3)
#mod_list=(no_geneflow_threepop earlyadmix_threepop earlygeneflow_threepop lateadmix_threepop lategeneflow_threepop singleadmix_pop1_threepop singleadmix_pop2_threepop singlegeneflow_pop1_threepop singlegeneflow_pop2_threepop alladmix_threepop)
mod_list=(`cat models_allsites.txt | tr '\n' ' '`)
#mod_list=(fbranch_model_new_new)

#path to directory where the model directories are located
file_path=/project/latesgenomics/jrick/admixture/july2020/fastsimcoal/six_taxon_models_noContam_jun2022


#########################
### RUNNING  ANALYSES ###
#########################

#add header to file that stores AIC values
echo -e  model"\t"deltaL"\t"AIC"\t"AIC25"\t"AIC50"\t"AIC200 > $file_path/sixtaxa_noContam_allsites_latesWGS_2.AIC

#loop through each model
for mod in ${mod_list[*]}
do
	cd $file_path/$mod/bestrun/
	/home/jrick/bin/fastsimcoal_scripts/scripts/calculateAIC.sh $mod
	
	echo -e $mod"\t"`tail -n 1 $mod.AIC` >> $file_path/sixtaxa_noContam_allsites_latesWGS_2.AIC
done

echo "done calculating AIC"
