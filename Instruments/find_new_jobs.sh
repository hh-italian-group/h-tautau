#!/bin/bash
# List jobs that do not yet been processed by tuple hadd.
# This file is part of https://github.com/hh-italian-group/h-tautau.

if [ $# -ne 2 ] ; then
    echo "Usage: job_list_file tuple_path" >&2
    exit 1
fi

JOB_LIST_FILE="$1"
TUPLE_PATH="$2"

if [ ! -f "$JOB_LIST_FILE" ] ; then
    echo "ERROR: can't find job list file '$JOB_LIST_FILE'." >&2
    exit 1
fi

if [ ! -d "$TUPLE_PATH" ] ; then
    echo "ERROR: can't find tuple path '$TUPLE_PATH'."  >&2
    exit 1
fi

for JOB in $(cat "$JOB_LIST_FILE") ; do
    JOB_ID=$( echo $JOB | sed -E 's/([0-9]+_[0-9]+).*/\1/' )
    JOB_NAME=$( echo $JOB | sed -E 's/[0-9]+_[0-9]+:[^_]*_//' )
    JOB_OUT_NAME=$( echo $JOB_NAME | sed -E 's/^crab_//' )

    if [ ! -f "$TUPLE_PATH/$JOB_OUT_NAME.root" ] ; then
        echo $JOB
    fi
done
