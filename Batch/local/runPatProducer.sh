#!/bin/bash
# Run PAT-production jobs with environment setup and output redirection.
# This file is part of https://github.com/hh-italian-group/h-tautau.

if [ $# -ne 9 ] ; then
    echo "Usage: job_name working_path file_list_path output_path global_tag is_mc include_sim n_events prefix"
    exit
fi

NAME=$1
WORKING_PATH=$2
FILE_LIST_PATH=$3
OUTPUT_PATH=$4
GLOBAL_TAG=$5
IS_MC=$6
INCLUDE_SIM=$7
N_EVENTS=$8
PREFIX=$9

KEEP_PAT=False
RUN_TREE=True

cd $WORKING_PATH

echo "$NAME $( date )" >> $OUTPUT_PATH/job_start.log

source cmsenv.sh
eval $( scramv1 runtime -sh )

echo "$NAME $( date )" >> $OUTPUT_PATH/job_cmsRun_start.log

if [ $PREFIX = "none" ] ; then
    cmsRun PatProduction/python/patTuple.py globalTag=$GLOBAL_TAG isMC=$IS_MC includeSim=$INCLUDE_SIM \
           fileList=$FILE_LIST_PATH/${NAME}.txt maxEvents=$N_EVENTS outputFile=$OUTPUT_PATH/${NAME}_Pat.root \
           keepPat=$KEEP_PAT runTree=$RUN_TREE treeOutput=$OUTPUT_PATH/${NAME}_Tree.root \
           > $OUTPUT_PATH/${NAME}_detail.log 2> $OUTPUT_PATH/${NAME}.log
    RESULT=$?
else
    cmsRun PatProduction/python/patTuple.py globalTag=$GLOBAL_TAG isMC=$IS_MC includeSim=$INCLUDE_SIM \
           fileList=$FILE_LIST_PATH/${NAME}.txt maxEvents=$N_EVENTS outputFile=$OUTPUT_PATH/${NAME}_Pat.root \
           keepPat=$KEEP_PAT runTree=$RUN_TREE treeOutput=$OUTPUT_PATH/${NAME}_Tree.root \
           fileNamePrefix=$PREFIX > $OUTPUT_PATH/${NAME}_detail.log 2> $OUTPUT_PATH/${NAME}.log
    RESULT=$?
fi


echo "$RESULT $NAME $( date )" >> $OUTPUT_PATH/job_result.log

exit $RESULT
