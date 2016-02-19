#!/bin/bash
# Run analysis jobs with environment setup and output redirection.
# This file is part of https://github.com/hh-italian-group/h-tautau.

if [ $# -lt 5 ] ; then
    echo "Usage: job_name working_path output_path exe_name set_cmsenv [other_args_for_analyzer]"
    exit
fi

NAME=$1
WORKING_PATH=$2
OUTPUT_PATH=$3
EXE_NAME=$4
SET_CMSENV=$5

cd $WORKING_PATH

echo "$NAME $( date )" >> $OUTPUT_PATH/job_start.log

if [ $SET_CMSENV = "yes" ] ; then
    source cmsenv.sh
    eval $( scramv1 runtime -sh )
    echo "$NAME $( date )" >> $OUTPUT_PATH/job_cmsRun_start.log
fi

$EXE_NAME "${@:6}" > $OUTPUT_PATH/${NAME}_detail.log 2> $OUTPUT_PATH/${NAME}.log
RESULT=$?

echo "$RESULT $NAME $( date )" >> $OUTPUT_PATH/job_result.log

exit $RESULT
