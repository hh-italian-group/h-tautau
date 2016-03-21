#!/bin/bash
# Resubmit failed analysis jobs on batch.
# This file is part of https://github.com/hh-italian-group/h-tautau.

if [ $# -ne 7 ] ; then
    echo "Usage: analyzer_cfg dataset_cfg output_path queue storage n_parallel_jobs multi_ana_name"
    exit
fi

ANA_CFG_FILE=$1
if [ ! -f "$ANA_CFG_FILE" ] ; then
    echo "ERROR: config file '$ANA_CFG_FILE' with analyzer list does not exists."
    exit
fi

CFG_FILE=$2
if [ ! -f "$CFG_FILE" ] ; then
    echo "ERROR: config file '$CFG_FILE' with dataset list does not exists."
    exit
fi

OUTPUT_PATH=$3
if [ ! -d "$OUTPUT_PATH" ] ; then
    echo "ERROR: output path '$OUTPUT_PATH' does not exists."
    exit
fi
OUTPUT_PATH=$( cd "$OUTPUT_PATH" ; pwd )

QUEUE=$4
STORAGE=$5
MAX_N_PARALLEL_JOBS=$6
#MULTI_ANA_NAME="H_BaselineSync"
MULTI_ANA_NAME=$7

REFERENCE_INPUT_PATH="Analysis/dataset"
REFERENCE_CFG_PATH="Analysis/config"
DATASET_ARRAY=( $( cat $CFG_FILE | awk '{ print $1 }' ) )
CFG_ARRAY=( $( cat $CFG_FILE | awk '{ print $2 }' ) )
N_DATASET=${#DATASET_ARRAY[@]}


ANALYZER_PATH=( $( cat $ANA_CFG_FILE ) )
WORKING_PATH=$CMSSW_BASE/src/HHbbTauTau
MAKE_PATH=$WORKING_PATH/RunTools/make_withFactory.sh

EXE_NAME=$OUTPUT_PATH/$MULTI_ANA_NAME
if [ ! -f $EXE_NAME ] ; then
    echo "Compiling executable file $EXE_NAME..."
    $MAKE_PATH $OUTPUT_PATH $MULTI_ANA_NAME $MULTI_ANA_NAME
    RESULT=$?
    if [ $RESULT -ne 0 ] ; then
        echo "Compilation of $EXE_NAME failed with an error code $RESULT. No jobs will be submited."
        exit
    fi
    echo "Executable file $EXE_NAME is compiled."
fi

for (( i=0; i<$N_DATASET; i++ )) ; do
    DATASET=${DATASET_ARRAY[i]}
    echo "Processing dataset $DATASET"

    FILE_LIST_PATH=$REFERENCE_INPUT_PATH/$DATASET
    if [ ! -d $FILE_LIST_PATH ] ; then
        echo "ERROR: file list path '$FILE_LIST_PATH' does not exists."
        exit
    fi
    JOBS=$( find $FILE_LIST_PATH -maxdepth 1 -name "*.txt" -printf "%f\n" | sed "s/\.txt//" | sort )
    if [ "x$JOBS" = "x" ] ; then
            echo "ERROR: directory '$FILE_LIST_PATH' does not contains any job description."
            exit
    fi
    N_JOBS=$( echo "$JOBS" | wc -l )
    echo "Total number of jobs: $N_JOBS"

    FILE_JOB_RESULT="$OUTPUT_PATH/job_result.log"
    SUCCESSFULL_JOBS=""
    N_SUCCESSFULL_JOBS=0
    if [ -f $FILE_JOB_RESULT ] ; then
        SUCCESSFULL_JOBS=$( cat $FILE_JOB_RESULT | sed -n "s/\(^0 \)\($DATASET[^ ]*\)\(.*\)/\2/p" | sort )
        N_SUCCESSFULL_JOBS=$( echo "$SUCCESSFULL_JOBS" | sed '/^\s*$/d' | wc -l )
        echo "Number of successfull jobs: $N_SUCCESSFULL_JOBS"
    else
        echo "job_results.log not found. Considering that there are no finished jobs yet."
    fi

    if [ $STORAGE = "Pisa" ] ; then
        JOB_STAT_RESULT=$( bjobs -w | grep $DATASET )
        JOBS_IN_QUEUE=$( echo "$JOB_STAT_RESULT" | awk '{print $7}' )
        N_RUNNING_JOBS=$( echo "$JOB_STAT_RESULT" | grep " RUN " | wc -l )
        N_PENDING_JOBS=$( echo "$JOB_STAT_RESULT" | grep " PEND " | wc -l )
    elif [ $STORAGE = "Bari" ] ; then
        JOB_STAT_RESULT=$( qstat -u $(whoami) | grep $DATASET )
        JOBS_IN_QUEUE=$( echo "$JOB_STAT_RESULT" | awk '{print $4}' )
        N_RUNNING_JOBS=$( echo "$JOB_STAT_RESULT" | grep " R " | wc -l )
        N_PENDING_JOBS=$( echo "$JOB_STAT_RESULT" | grep " Q " | wc -l )
    else
        echo "ERROR: unknown storage"
        exit
    fi

    echo "Number of running jobs: $N_RUNNING_JOBS"
    echo "Number of pending jobs: $N_PENDING_JOBS"

    JOBS_TO_RESUBMIT=$( echo -e "$JOBS"\\n"$SUCCESSFULL_JOBS"\\n"$JOBS_IN_QUEUE" | sort | sed '/^\s*$/d' | uniq -u )
    N_JOBS_TO_RESUBMIT=$( echo "$JOBS_TO_RESUBMIT" | sed '/^\s*$/d' | wc -l )
    if [ $N_JOBS_TO_RESUBMIT -eq 0 ] ; then
        echo "There are no jobs to resubmit."
        continue
    fi
    echo "Jobs to resubmit: "$JOBS_TO_RESUBMIT
    read -p "Number of jobs to resubmit: $N_JOBS_TO_RESUBMIT. Continue? (yes/no) " -r REPLAY
    if [ "$REPLAY" != "y" -a "$REPLAY" != "yes" -a "$REPLAY" != "Y" ] ; then
        echo "No jobs have been resubmited."
        continue
    fi

    NEW_FILE_LIST_PATH="Analysis/dataset/retry/${DATASET}"
    n=1
    while [ -d $NEW_FILE_LIST_PATH ] ; do
        n=$(($n + 1))
        NEW_FILE_LIST_PATH="Analysis/dataset/retry/${DATASET}_${n}"
    done

    echo "Resubmit number $n"
    mkdir -p $NEW_FILE_LIST_PATH
    echo "$JOBS_TO_RESUBMIT" | xargs -n 1 printf "$FILE_LIST_PATH/%b.txt $NEW_FILE_LIST_PATH/\n" | xargs -n 2 cp

    ANA_OUTPUT_PATH_1=$OUTPUT_PATH/${ANALYZER_PATH[0]}/$DATASET
    mkdir -p $ANA_OUTPUT_PATH_1

    ANA_OUTPUT_PATH_2=$OUTPUT_PATH/${ANALYZER_PATH[1]}/$DATASET
    mkdir -p $ANA_OUTPUT_PATH_2

    ANA_OUTPUT_PATH_3=$OUTPUT_PATH/${ANALYZER_PATH[2]}/$DATASET
    mkdir -p $ANA_OUTPUT_PATH_3

    yes | ./RunTools/submitMultiAnalysis_Batch.sh $QUEUE $STORAGE $MAX_N_PARALLEL_JOBS $NEW_FILE_LIST_PATH \
        $OUTPUT_PATH $ANA_OUTPUT_PATH_1 $ANA_OUTPUT_PATH_2 $ANA_OUTPUT_PATH_3 \
        $REFERENCE_CFG_PATH/${CFG_ARRAY[i]} $EXE_NAME
done

echo "All datasets has been processed."
