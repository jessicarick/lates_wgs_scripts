#!/bin/sh

date
## workflow for TWISST using simon martin's genomics_general and twisst scripts
module load gcc miniconda3 jdk
module load openmpi raxml
  
source activate biopy-env
  
input=$1 # vcf file
output=`basename $input | sed 's/\.vcf/\.beagle/g'`

# impute using 10000bp windows with 1000bp overlap
java -Xmx12g -jar beagle.r1399.jar gt=$input out=${output} impute=true nthreads=16 window=10000 overlap=1000 gprobs=TRUE

# before this step, genotypes need to be phased! it won't work if gentoypes are 0/0 instead of 0|0
python genomics_general/VCF_processing/parseVCF.py -i ${output}.vcf.gz --skipIndel | gzip > ${output}.phased.geno.gz

# create 100bp sliding window trees
python genomics_general/raxml_sliding_windows_py3.py -T 16 -g ${output}.phased.geno.gz --prefix ${output}.raxml.w100 -w 100 --windType sites --model GTRCAT --raxml raxmlHPC-PTHREADS

echo "all done!"
date
