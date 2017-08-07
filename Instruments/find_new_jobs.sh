#!/bin/bash
# List jobs that do not yet been processed by tuple hadd.
# This file is part of https://github.com/hh-italian-group/h-tautau.

if [ $# -lt 2 ] ; then
    echo "Usage: job_list_file tuple_path [tuple_path] ..." >&2
    exit 1
fi

JOB_LIST_FILE="$1"
TUPLE_PATHS="${@:2}"

if [ ! -f "$JOB_LIST_FILE" ] ; then
    echo "ERROR: can't find job list file '$JOB_LIST_FILE'." >&2
    exit 1
fi

cat "$JOB_LIST_FILE" | while read -r JOB ; do
    JOB_ID=$( echo $JOB | sed -E 's/([0-9]+_[0-9]+).*/\1/' )
    JOB_NAME=$( echo $JOB | sed -E 's/[0-9]+_[0-9]+:[^_]*_([^ ]*).*/\1/' )
    JOB_OUT_NAME=$( echo $JOB_NAME | sed -E 's/^crab_//' )

    OUTPUT_FOUND=0
    for OUTPUT in $TUPLE_PATHS ; do
        if [ -f "$OUTPUT/$JOB_OUT_NAME.root" ] ; then
            OUTPUT_FOUND=1
            break
        fi
    done
    if [ $OUTPUT_FOUND -eq 0 ] ; then
        echo "$JOB"
    fi
done
