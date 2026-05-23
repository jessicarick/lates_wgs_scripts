#!/bin/sh
# script documenting calls for running Dsuite
# dsuite located at /project/latesgenomics/software/Dsuite/Build/Dsuite
alias Dsuite='/project/latesgenomics/software/Dsuite/Build/Dsuite'

Dsuite Dtrios -k 100 \
  -t ../tree_wgs.nwk \
  -o latesWGS_2_052326 \
  -c --ABBAclustering \
  ../../latesWGS_2.bamlist.all_miss0.5_maf0.01_minDP30.recode.vcf sets.txt

Dsuite Dinvestigate \
  -w 100,20 \
  ../../latesWGS_2.bamlist.all_miss0.5_maf0.01_minDP30.recode.vcf \
  sets.txt test_trios.txt

Dsuite Fbranch -Z \
  ../tree_wgs.nwk \
  latesWGS_2_052326_tree.txt
