#!/bin/bash
# Submit hadd on CRAB tuple outputs.
# This file is part of https://github.com/hh-italian-group/h-tautau.

if [ $# -lt 4 ] ; then
    echo "Usage: queue job_list_file tuple_output_path crab_output_path [crab_output_path] ..." >&2
    exit 1
fi

QUEUE="$1"
JOB_LIST_FILE="$2"
TUPLE_OUTPUT_PATH="$3"
CRAB_OUTPUT_PATHS="${@:4}"

if [ ! -f "$JOB_LIST_FILE" ] ; then
    echo "ERROR: can't find job list file '$JOB_LIST_FILE'." >&2
    exit 1
fi

mkdir -p "$TUPLE_OUTPUT_PATH"
if [ ! -d "$TUPLE_OUTPUT_PATH" ] ; then
    echo "ERROR: can't create tuple output path '$TUPLE_OUTPUT_PATH'."  >&2
    exit 1
fi

for JOB in $(cat "$JOB_LIST_FILE") ; do
    JOB_ID=$( echo $JOB | sed -E 's/([0-9]+_[0-9]+).*/\1/' )
    JOB_NAME=$( echo $JOB | sed -E 's/[0-9]+_[0-9]+:[^_]*_//' )
    JOB_PATHS=( $( find $CRAB_OUTPUT_PATHS -maxdepth 3 -mindepth 3 -type d -path '*'"/$JOB_NAME/$JOB_ID") )
    if [ ${#JOB_PATHS[@]} -eq 0 ] ; then
        echo "ERROR: can't find a path with the crab ouput for job '$JOB'." >&2
        exit 1
    fi
    if [ ${#JOB_PATHS[@]} -gt 1 ] ; then
        echo "ERROR: found more than one path for job '$JOB':" >&2
        echo "${JOB_PATHS[@]}" >&2
        exit 1
    fi
    JOB_PATH="$(cd "${JOB_PATHS[0]}" && pwd)"

    LOCAL_JOB_NAME="${JOB_ID}_${JOB_NAME}"
    ./AnalysisTools/Run/submit_job.sh "$QUEUE" "$LOCAL_JOB_NAME" "$TUPLE_OUTPUT_PATH/$LOCAL_JOB_NAME" \
                                      hadd -f9 -n 11 "${JOB_NAME}.root" "$JOB_PATH/*/*.root"
done
