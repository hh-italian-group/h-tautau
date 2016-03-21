#!/bin/bash
# Submit all available analysis datasets on batch.
# This file is part of https://github.com/hh-italian-group/h-tautau.

if [ $# -ne 8 ] ; then
    echo "Usage: analyzer_cfg dataset_cfg output_path queue storage n_parallel_jobs use_multi_ana multi_ana_name"
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
if [ -d "$OUTPUT_PATH" ] ; then
    echo "ERROR: output path '$OUTPUT_PATH' already exists."
    exit
fi
mkdir -p $OUTPUT_PATH
OUTPUT_PATH=$( cd "$OUTPUT_PATH" ; pwd )

QUEUE=$4
STORAGE=$5
MAX_N_PARALLEL_JOBS=$6
USE_MULTI_ANA=$7
MULTI_ANA_NAME=$8
#MULTI_ANA_NAME="H_BaselineSync"

WORKING_PATH=$CMSSW_BASE/src/HHbbTauTau
MAKE_PATH=$WORKING_PATH/RunTools/make_withFactory.sh

REFERENCE_INPUT_PATH="Analysis/dataset"
REFERENCE_CFG_PATH="Analysis/config"
DATASET_ARRAY=( $( cat $CFG_FILE | awk '{ print $1 }' ) )
CFG_ARRAY=( $( cat $CFG_FILE | awk '{ print $2 }' ) )

N_DATASET=${#DATASET_ARRAY[@]}

ANALYZER_PATH=( $( cat $ANA_CFG_FILE ) )

if [ $USE_MULTI_ANA = "yes" ] ; then
    EXE_NAME=$OUTPUT_PATH/$MULTI_ANA_NAME

    echo "Compiling executable file $EXE_NAME..."
    $MAKE_PATH $OUTPUT_PATH $MULTI_ANA_NAME $MULTI_ANA_NAME
    RESULT=$?
    if [ $RESULT -ne 0 ] ; then
        echo "Compilation of $EXE_NAME failed with an error code $RESULT. No jobs will be submited."
        exit
    fi
    echo "Executable file $EXE_NAME is compiled."


    for (( i=0; i<$N_DATASET; i++ )) ; do
        INPUT_PATH=$REFERENCE_INPUT_PATH/${DATASET_ARRAY[i]}

        if [ ! -d "$INPUT_PATH" ] ; then
            echo "ERROR: dataset '${DATASET_ARRAY[i]}' does not exist."
            exit
        fi

        ANA_OUTPUT_PATH_1=$OUTPUT_PATH/${ANALYZER_PATH[0]}/${DATASET_ARRAY[i]}
        mkdir -p $ANA_OUTPUT_PATH_1

        ANA_OUTPUT_PATH_2=$OUTPUT_PATH/${ANALYZER_PATH[1]}/${DATASET_ARRAY[i]}
        mkdir -p $ANA_OUTPUT_PATH_2

        ANA_OUTPUT_PATH_3=$OUTPUT_PATH/${ANALYZER_PATH[2]}/${DATASET_ARRAY[i]}
        mkdir -p $ANA_OUTPUT_PATH_3

        yes | ./RunTools/submitMultiAnalysis_Batch.sh $QUEUE $STORAGE $MAX_N_PARALLEL_JOBS $INPUT_PATH \
            $OUTPUT_PATH $ANA_OUTPUT_PATH_1 $ANA_OUTPUT_PATH_2 $ANA_OUTPUT_PATH_3 \
            $REFERENCE_CFG_PATH/${CFG_ARRAY[i]} $EXE_NAME
    done

else

    for ANA_FOLDER in $ANALYZER_PATH ; do
        mkdir -p $OUTPUT_PATH/$ANA_FOLDER
        EXE_NAME=$OUTPUT_PATH/$ANA_FOLDER/$ANA_FOLDER

        echo "Compiling executable file $EXE_NAME..."
        $MAKE_PATH $OUTPUT_PATH/$ANA_FOLDER $ANA_FOLDER $ANA_FOLDER
        RESULT=$?
        if [ $RESULT -ne 0 ] ; then
            echo "Compilation of $EXE_NAME failed with an error code $RESULT. No jobs will be submited."
            exit
        fi
        echo "Executable file $EXE_NAME is compiled."

        for (( i=0; i<$N_DATASET; i++ )) ; do
            INPUT_PATH=$REFERENCE_INPUT_PATH/${DATASET_ARRAY[i]}

            if [ ! -d "$INPUT_PATH" ] ; then
                echo "ERROR: dataset '${DATASET_ARRAY[i]}' does not exist."
                exit
            fi

            ANA_OUTPUT_PATH=$OUTPUT_PATH/$ANA_FOLDER/${DATASET_ARRAY[i]}
            mkdir -p $ANA_OUTPUT_PATH

            echo $ANA_FOLDER ${DATASET_ARRAY[i]} ${CFG_ARRAY[i]}

            yes | ./RunTools/submitAnalysis_Batch.sh $QUEUE $STORAGE $MAX_N_PARALLEL_JOBS $INPUT_PATH $ANA_OUTPUT_PATH \
                $ANA_FOLDER $REFERENCE_CFG_PATH/${CFG_ARRAY[i]} $EXE_NAME

        done

    done

fi

