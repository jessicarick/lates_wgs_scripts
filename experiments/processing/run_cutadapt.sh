#!/bin/sh

#SBATCH --account=latesgenomics
#SBATCH --job-name="run_cutadapt"
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load gcc
module load py-cutadapt

fastq1=$1
fastq2=$2

date

prefix=`basename $fastq1 | sed 's/\(\w\)_s[0-9].*/\1/'`
echo "working with $prefix, which has $# arguments"

if [[ $# -eq 4 ]]; then 
	adapt="-a $3 -A $4"
elif [[ $# -eq 5 ]]; then
	adapt="-a $3 -A $4 -A $5"
elif [[ $# -eq 6 ]]; then
	adapt="-a $3 -a $4 -A $5 -A $6"
else
	echo "I do not understand your adapter commands. Please try again."
	exit 0
fi	

cutadapt $adapt -o ${prefix}_1.clean.fastq -p ${prefix}_2.clean.fastq $fastq1 $fastq2

echo "done with $prefix"
date
