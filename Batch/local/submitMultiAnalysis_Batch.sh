#!/bin/bash
# Submit analysis jobs for a given dataset on batch for all analysis channels.
# This file is part of https://github.com/hh-italian-group/h-tautau.

if [ $# -ne 10 ] ; then
    echo "Usage: queue storage max_n_parallel_jobs file_list_path output_path_ref output_path_1 output_path_2 output_path_3 config_name exe_name"
    exit
fi

QUEUE=$1
STORAGE=$2
MAX_N_PARALLEL_JOBS=$3
FILE_LIST_PATH=$4
OUTPUT_PATH_REF=$5
OUTPUT_PATH_1=$6
OUTPUT_PATH_2=$7
OUTPUT_PATH_3=$8
CONFIG_NAME=$9
EXE_NAME=${10}

WORKING_PATH=$CMSSW_BASE/src/HHbbTauTau
RUN_SCRIPT_PATH=$WORKING_PATH/RunTools/runAnalysis.sh

MAX_N_EVENTS=0

if [ $STORAGE = "Pisa" ] ; then
    PREFIX="/gpfs/ddn/srm/cms"
elif [ $STORAGE = "Bari" ] ; then
    PREFIX="/lustre/cms"
elif [ $STORAGE = "Local" ] ; then
    PREFIX=$CMS_STORE
else
    echo "ERROR: unknown storage"
    exit
fi

if [ ! -d "$WORKING_PATH" ] ; then
        echo "ERROR: working path '$WORKING_PATH' does not exist."
        exit
fi

if [ ! -d "$WORKING_PATH/$FILE_LIST_PATH" ] ; then
        echo "ERROR: file list path '$WORKING_PATH/$FILE_LIST_PATH' does not exist."
        exit
fi

if [ ! -d "$OUTPUT_PATH_REF" ] ; then
    echo "ERROR: output ref path '$OUTPUT_PATH_REF' does not exist."
    exit
fi
OUTPUT_PATH_REF=$( cd "$OUTPUT_PATH_REF" ; pwd )

if [ ! -d "$OUTPUT_PATH_1" ] ; then
    echo "ERROR: output path 1 '$OUTPUT_PATH_1' does not exist."
        exit
fi
OUTPUT_PATH_1=$( cd "$OUTPUT_PATH_1" ; pwd )

if [ ! -d "$OUTPUT_PATH_2" ] ; then
    echo "ERROR: output path 2 '$OUTPUT_PATH_2' does not exist."
    exit
fi
OUTPUT_PATH_2=$( cd "$OUTPUT_PATH_2" ; pwd )

if [ ! -d "$OUTPUT_PATH_3" ] ; then
    echo "ERROR: output path 3 '$OUTPUT_PATH_3' does not exist."
    exit
fi
OUTPUT_PATH_3=$( cd "$OUTPUT_PATH_3" ; pwd )

if [ ! -f "$RUN_SCRIPT_PATH" ] ; then
        echo "ERROR: script '$RUN_SCRIPT_PATH' does not exist."
        exit
fi

JOBS=$( find $WORKING_PATH/$FILE_LIST_PATH -maxdepth 1 -name "*.txt" -print0 | xargs -0 -n 1 basename | sed "s/\.txt//" )

if [ "x$JOBS" = "x" ] ; then
        echo "ERROR: directory '$FILE_LIST_PATH' does not contains any job description."
        exit
fi

N_JOBS=$( echo "$JOBS" | wc -l )
echo "Following jobs will be submited:" $JOBS
echo "Total number of jobs to submit: $N_JOBS"

read -p "Compile these jobs and then submit (yes/no)? " -r REPLAY
if [ "$REPLAY" != "y" -a "$REPLAY" != "yes" -a "$REPLAY" != "Y" ] ; then
    echo "No jobs have been compiled or submitted."
    exit
fi

echo "Using pre-compiled executable $EXE_NAME."

i=0
n=0

source $WORKING_PATH/RunTools/batch.sh

if [ $STORAGE = "Local" ] ; then
    SET_CMS_ENV="dont_set_cmsenv"
else
    SET_CMS_ENV="yes"
fi

for NAME in $JOBS ; do
    submit_batch $QUEUE $STORAGE $NAME $OUTPUT_PATH_REF $RUN_SCRIPT_PATH $NAME $WORKING_PATH $OUTPUT_PATH_REF \
                 $EXE_NAME $SET_CMS_ENV $FILE_LIST_PATH/${NAME}.txt \
                 $OUTPUT_PATH_1/${NAME}.root $OUTPUT_PATH_2/${NAME}.root $OUTPUT_PATH_3/${NAME}.root \
                 $CONFIG_NAME $PREFIX @$MAX_N_EVENTS
    if [ $MAX_N_PARALLEL_JOBS -ne 0 ] ; then
        i=$(($i + 1))
        n=$(($n + 1))
        echo "job $n started"
        if [[ $i == $MAX_N_PARALLEL_JOBS ]] ; then
                wait
                i=0
        fi
    fi
done

if [ $MAX_N_PARALLEL_JOBS -eq 0 ] ; then
    echo "$N_JOBS have been submitted."
else
    wait
    echo "$N_JOBS finished"
fi
