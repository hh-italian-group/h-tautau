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
MERGE_EXE=$CMSSW_BASE/build/h-tautau/TupleMerger

if [ ! -f $MERGE_EXE ] ; then
	echo "ERROR: TupleMerger executable not found."
	exit 1
fi

if [ ! -f "$JOB_LIST_FILE" ] ; then
    echo "ERROR: can't find job list file '$JOB_LIST_FILE'." >&2
    exit 1
fi

mkdir -p "$TUPLE_OUTPUT_PATH"
if [ ! -d "$TUPLE_OUTPUT_PATH" ] ; then
    echo "ERROR: can't create tuple output path '$TUPLE_OUTPUT_PATH'."  >&2
    exit 1
fi

cat "$JOB_LIST_FILE" | while read -r JOB ; do
    JOB_ID=$( echo $JOB | sed -E 's/([0-9]+_[0-9]+).*/\1/' )
    JOB_NAME=$( echo $JOB | sed -E 's/[0-9]+_[0-9]+:[^_]*_([^ ]*).*/\1/' )
    JOB_PATHS=( $( find $CRAB_OUTPUT_PATHS -maxdepth 3 -mindepth 3 -type d -path '*'"/$JOB_NAME/$JOB_ID") )
    FAILED_IDS=$( echo $JOB | awk '{print $2}' )
    if [ "x$FAILED_IDS" = "x" ] ; then
        EXCLUDE_LIST=""
    else
        EXCLUDE_FILES=$( echo $FAILED_IDS | sed -E 's/([0-9]+)/eventTuple_\1\.root/g' )
        EXCLUDE_LIST="--exclude-list $EXCLUDE_FILES"
    fi
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
                                      $MERGE_EXE --output "${JOB_NAME}.root" --input-dir "$JOB_PATH" $EXCLUDE_LIST \
                                      --exclude-dir-list failed
done
