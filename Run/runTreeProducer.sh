#!/bin/bash
# Run TREE production jobs with environment setup and output redirection.
# This file is part of https://github.com/hh-italian-group/h-tautau.

if [ $# -ne 8 ] ; then
    echo "Usage: job_name working_path file_list_path output_path global_tag include_sim n_events prefix"
    exit
fi

NAME=$1
WORKING_PATH=$2
FILE_LIST_PATH=$3
OUTPUT_PATH=$4
GLOBAL_TAG=$5
INCLUDE_SIM=$6
N_EVENTS=$7
PREFIX=$8

cd $WORKING_PATH

echo "$NAME $( date )" >> $OUTPUT_PATH/job_start.log

source cmsenv.sh
eval $( scramv1 runtime -sh )

echo "$NAME $( date )" >> $OUTPUT_PATH/job_cmsRun_start.log

if [ $PREFIX = "none" ] ; then
    nohup cmsRun TreeProduction/python/treeProducer.py globalTag=$GLOBAL_TAG includeSim=$INCLUDE_SIM \
           fileList=$FILE_LIST_PATH/${NAME}.txt maxEvents=$N_EVENTS outputFile=$OUTPUT_PATH/${NAME}_Tree.root \
           > $OUTPUT_PATH/${NAME}_detail.log 2> $OUTPUT_PATH/${NAME}.log
    RESULT=$?

else
    nohup cmsRun TreeProduction/python/treeProducer.py globalTag=$GLOBAL_TAG includeSim=$INCLUDE_SIM \
           fileList=$FILE_LIST_PATH/${NAME}.txt maxEvents=$N_EVENTS outputFile=$OUTPUT_PATH/${NAME}_Tree.root \
           fileNamePrefix=$PREFIX > $OUTPUT_PATH/${NAME}_detail.log 2> $OUTPUT_PATH/${NAME}.log
    RESULT=$?
fi


echo "$RESULT $NAME $( date )" >> $OUTPUT_PATH/job_result.log

exit $RESULT
