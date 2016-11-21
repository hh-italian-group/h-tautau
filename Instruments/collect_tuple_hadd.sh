#!/bin/bash
# Collect hadd outputs outputs.
# This file is part of https://github.com/hh-italian-group/h-tautau.

if [ $# -ne 2 ] ; then
    echo "Usage: input_path output_path" >&2
    exit 1
fi

INPUT_PATH="$1"
OUTPUT_PATH="$2"

if [ ! -d "$INPUT_PATH" ] ; then
    echo "ERROR: can't find input path '$INPUT_PATH'."  >&2
    exit 1
fi

if [ ! -d "$OUTPUT_PATH" ] ; then
    echo "ERROR: can't find output path '$OUTPUT_PATH'."  >&2
    exit 1
fi

LOG_NAME="run_job.log"
cd "$INPUT_PATH"
for JOB in $(ls) ; do
    JOB_NAME=$( echo $JOB | sed -E 's/^[0-9]+_[0-9]+_//' )
	JOB_OUT_NAME=$( echo $JOB_NAME | sed -E 's/^crab_//' )
    if [ ! -f "$JOB/$LOG_NAME" ] ; then
        echo "ERROR: can't find log file for job '$JOB'." >&2
        exit 1
    fi
    LAST_LOG_ENTRY="$( tail -n 1 "$JOB/$LOG_NAME" )"

    if [[ $LAST_LOG_ENTRY =~ ^Job\ successfully\ ended\ at.* ]] ; then
        if [ -f "$OUTPUT_PATH/$JOB_OUT_NAME.root" ] ; then
            #echo "$JOB_NAME has already been transfered."
            continue
        fi
        if [ ! -f "$JOB/$JOB_NAME.root" ] ; then
            echo "ERROR: can't find root file for successfully finished job '$JOB'." >&2
            exit 1
        fi
        mv "$JOB/$JOB_NAME.root" "$OUTPUT_PATH/$JOB_OUT_NAME.root"
        RESULT=$?
        if [ $RESULT -ne 0 ] ; then
            echo "ERROR: can't transfer root file '$JOB/$JOB_NAME.root' into output directory '$OUTPUT_PATH'." >&2
            exit 1
        fi
        #rm -r "$JOB"
        echo "$JOB_NAME transfered into '$OUTPUT_PATH'."
    elif [[ $LAST_LOG_ENTRY =~ ^Job\ failed\ at.* ]] ; then
        echo "$JOB failed."
    elif [[ $LAST_LOG_ENTRY =~ ^Job\ submission\ at.* ]] ; then
        echo "$JOB not yet started."
    else
        echo "$JOB not yet ended."
    fi

done
