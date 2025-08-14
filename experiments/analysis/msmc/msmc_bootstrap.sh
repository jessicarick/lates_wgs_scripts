#!/bin/bash -l

#for x in `cat IND`; \
#do echo $i; \
#IND=`echo ${i[@]}`; \
#for s in `cat SCAFFOLDS`; \
#do echo $s;\
#SCAFS=`echo ${s[@]}`; \
#done; done

FOLDER=$1
POP_OR_IND=$2

BASEDIR=/project/passerinagenome/jrick/msmc
RUN_NAME=${POP_OR_IND}_bootstrap
P_PAR=1*2+25*1+1*2+1*3

BS_LIST=${BASEDIR}/bootstrap/${POP_OR_IND}_bs_file_list.txt

module load msmc2

##### run msmc alone #### 
MSMC_BS=`ls ${BASEDIR}/bootstrap/${FOLDER}/*.txt`
MSMC_OUTPUT=${BASEDIR}/bootstrap/msmc_output.${RUN_NAME}.${FOLDER}

msmc2 -t 16 -p $P_PAR -o $MSMC_OUTPUT -I 0,1 $MSMC_BS

