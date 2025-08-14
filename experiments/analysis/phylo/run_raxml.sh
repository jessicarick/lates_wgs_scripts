#!/bin/bash

date
echo "loading modules"

module load gcc
module load raxml

echo "modules loaded, ready to go!"

##### Usage: sbatch run_raxml.sh file.phy outprefix

PHYLIP=$1
OUT=$2

echo "This is your shell script reporting for duty! Time to create a tree!!"
echo "Running RAxML on ${PHYLIP}"

SEED=`date +%N | sed -e 's/000$//' -e 's/^0//'` #sets seed based on clock

raxmlHPC-PTHREADS-AVX -s $PHYLIP -f a -p $SEED -n $OUT -m ASC_GTRGAMMA --asc-corr=lewis -T 32 -x $SEED -N 100

echo "Done RAxMLing!" | mail -s "done with $PHYLIP" jrick@uwyo.edu
date
