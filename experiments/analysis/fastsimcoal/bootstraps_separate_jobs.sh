#!/bin/sh

### Job Name
#SBATCH --job-name fsc_boot

### Declare an account for the job to run under - You might need to change this to Carling Lab (check designation with arccquota on Moran)
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
#SBATCH --mem=120G

### Set max walltime (days-hours:minutes:seconds)
#SBATCH --time=2-01:00:00

module load gcc htslib miniconda3 python/3.6.3 r

source activate dadi-env

pop_file=$1
PREFIX=$2
vcf=$3

##################################
### 1. SETTING UP THE ANALYSES ###
##################################

### 1a. PATHS ###
#the pop spec file with the full file path
pop_spec=/project/latesgenomics/jrick/admixture/july2020/fastsimcoal/six_taxon_models_noContam_jun2022/$pop_file

#path to directory where the model directories are located
file_path=/project/latesgenomics/jrick/admixture/july2020/fastsimcoal/six_taxon_models_noContam_jun2022/${PREFIX}/bestrun

#path to the tpl and est files
est_tpl_path=$file_path

## 1b. BOOTSTRAP DATASET SPECIFICATION ###
#the number of sites per block for creating the bootstrap files
nsites=`grep -v '^#' $vcf | wc -l`
#nsites=`cat /project/latesgenomics/jrick/admixture/july2020/fastsimcoal/six_taxon_models_noContam_jun2022/tmp/*allSites | wc -l`
block_size=100000
nblocks=`echo $(( nsites / block_size ))`

#should the last site file be removed before creating the bootstrap vcf files? (options: "yes", "no")
#NOTE: this should only be yes if the last block has very few sites compared to the other files (the "leftover" sites)
remove_last_file="yes"


### 1c. ANALYSES SPECIFICATIONS ###
#number of runs per bootstrap dataset
number_of_runs=100

number_of_sims=100000

#number of bootstrap datasets
bootstraps=100


### 1d. NAMING DETAILS ###
#what do you want to name the bootstrap results folder to be named?
storage_folder="bootstrap_confint"

proj=`head -n 2 $file_path/*MSFS.obs | tail -n 1 | cut -f 2 | tr ' ' ','`

echo "projection is $proj"

#if false; then # DEBUG

############################
### 2. RUNNING  ANALYSES ###
############################

#make file to store analyses
mkdir ${file_path}/${storage_folder}/

#enter into store folder
cd $file_path/$storage_folder

if false; then
### 2a. CREATE BOOTSTRAP DATASETS ###
# Get all lines with genomic data
zgrep -v "^#" $vcf > $PREFIX.allSites

# Get the header
zgrep "^#" $vcf > header

# get 100 files with total_number_sites/block_size sites each
split -l $block_size $PREFIX.allSites $PREFIX.sites.

if [ $remove_last_file = "yes" ]
then
  ls $PREFIX.sites.* | sort | tail -1 | xargs rm
fi

#fi # debug

# create SFS for each block
# loop could be parallelized to make it faster
for file in `ls ${PREFIX}.sites* | grep -v ".obs"`
	do echo $file
	cat header $file > tmp
	/home/jrick/bin/easySFS/easySFS.py -i tmp -p $pop_spec -a -f --proj=$proj --prefix $file
	cp output/fastsimcoal2/*_MSFS.obs .
	rm -f $file
done

fi

# Generate $bootstraps files each with randomly concatenated blocks and compute the SFS for each:
for i in `seq 1 100`
do
  # Make a new folder for each bootstrapping iteration:
  mkdir bs$i
  cd bs$i

  # Randomly add nblocks blocks
  # keep only lines > 100 characters (these are the actual SFS; this omits header info)
  if [ $remove_last_file = "yes" ]
	then
	nblocks=`echo $(( nblocks - 1 ))`
  fi

  cat `shuf -e -r -n $nblocks ${file_path}/${storage_folder}/*_MSFS.obs` | grep '.\{100\}' > ${PREFIX}.bs.${i}.combinedMSFS.obs
  Rscript /project/latesgenomics/jrick/admixture/july2020/fastsimcoal/six_taxon_models_noContam/combine_MSFS.R ${PREFIX}.bs.${i}.combinedMSFS.obs tmp.bs.${i}_MSFS.obs
  head -n 2 ${file_path}/${storage_folder}/*.aa_MSFS.obs > ${PREFIX}.bs.${i}_MSFS.obs
  cat tmp.bs.${i}_MSFS.obs >> ${PREFIX}.bs.${i}_MSFS.obs

  # Copy the observed SFS file into this folder renaming it to match the .tpl prefix
  #cp ../${PREFIX}_jointMAFpop1_0.obs ${PREFIX}.bs.${i}_jointMAFpop1_0.obs

  cd ..
done

#fi #DEBUG

### 2b. RUN FASTSIMCOAL2 ON EACH BOOTSTRAP DATASET ###

#copy the .est and .tpl files into the main bootstrap folder
cp $est_tpl_path/${PREFIX}.est $est_tpl_path/${PREFIX}.tpl $file_path/$storage_folder/

#make sure that you are in the main bootstrap folder
cd $file_path/$storage_folder

## cleanup!

#rm -f ${PREFIX}.sites.*
#rm -f ${PREFIX}.allSites

#submit separate jobs for each bootstrap
if false; then
for bs in `seq 1 $bootstraps`
do
  sbatch /project/latesgenomics/jrick/admixture/july2020/fastsimcoal/six_taxon_models/run_single_bootstrap.sh $bs $file_path $storage_folder $number_of_runs $PREFIX $number_of_sims
  sleep 1
done

echo "done submitting bootstrap jobs"

fi
echo "not submitting bootstrap jobs, but they're ready!"

#loop through each bootstrap dataset
#for bs in $(eval echo "{1..$bootstraps}")
#do
	
  #enter $bs bootstrap folder
#  cd bs$bs
	
  #complete $number_of_runs runs for the bootstrap
#  for i in $(eval echo "{1..$number_of_runs}")
#  do
#    mkdir run$i #create folder for run i
    	
    #copy and rename tpl and est files into the run folder
#    cp $file_path/$storage_folder/${PREFIX}.tpl run$i/${PREFIX}.bs.$bs.tpl
#    cp $file_path/$storage_folder/${PREFIX}.est run$i/${PREFIX}.bs.$bs.est
#    cp $file_path/$storage_folder/bs$bs/${PREFIX}.bs.${bs}_jointMAFpop1_0.obs run$i"/"

#    cd run$i #enter the run folder
    	
    #run faststimcoal
#    /project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/fastsimcoal/fsc26_linux64/fsc26 -t ${PREFIX}.bs.$bs.tpl -e ${PREFIX}.bs.$bs.est -m -0 -C 10 -n $number_of_sims -L 50 -s0 -M -q -c 16
#    cd .. #leave the run folder (go up 1 level in the hierarchy)
    	
#  done
	
  #identify run with best likelihood (after all the runs are complete)
#  /project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/scripts/fsc_selectbestrun.sh
#  /project/ltcichlidgenomics/alewansk/Petrochromis_FULLDATASET/demographic_mod/scripts/bootstrap_selectbestrun.sh
#  cd ..
#done


### 2c. COMPILE THE PARAMETER VALUES OF THE BEST RUN FOR EACH BOOTSTRAP INTO A SINGLE FILE ###

#make sure that you are in the bootstrap folder
#cd $file_path/$storage_folder

#for bs in $(eval echo "{1..$bootstraps}")
#do 

  #use the .bestlhoods file from bs1 as the start of the bs_bestrun_param.txt file (in order to include the parameter names in the file)
#  if [ $bs == 1 ]
#  then
#    cp bs$bs/bestrun/recent_geneflow.bs.$bs.bestlhoods bs_bestrun_param.txt
   
  #add the parameter values of the rest of the bs bestruns to bs_bestrun_param.txt
#  else
#    tail -1 bs$bs/bestrun/recent_geneflow.bs.$bs.bestlhoods >> bs_bestrun_param.txt
#  fi
#done
