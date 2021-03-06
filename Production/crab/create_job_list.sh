#!/bin/bash
# Create list of CRAB jobs.
# This file is part of https://github.com/hh-italian-group/h-tautau.

if [ $# -ne 1 ] ; then
    echo "Usage: work_area"
    exit 1
fi

WORK_AREA="$1"

if [ ! -d "$WORK_AREA" ] ; then
    echo "ERROR: can't find crab work area '$WORK_AREA'." >&2
    exit 1
fi

for JOB in $(ls "$WORK_AREA") ; do
    JOB_ID_FILE="$WORK_AREA/$JOB/job_id.txt"
    if [ -f "$JOB_ID_FILE" ] ; then
        JOB_ID=$(cat "$JOB_ID_FILE")
    else
    	CRAB_LOG="$WORK_AREA/$JOB/crab.log"
    	if [ ! -f "$CRAB_LOG" ] ; then
    		echo "ERROR: can't find log for job '$JOB'." >&2
    		exit 1
    	fi
    	JOB_IDS=( $(grep -E ".*Task name:" "$CRAB_LOG" | sed -E 's/.*Task name\:[^0-9]*(.*)/\1/' ) )
    	if [ ${#JOB_IDS[@]} -lt 1 ] ; then
    		echo "ERROR: can't find job id inside '$CRAB_LOG'." >&2
    		exit 1
    	fi
        JOB_ID=${JOB_IDS[0]}
        echo $JOB_ID > "$JOB_ID_FILE"
    fi
	echo $JOB_ID
done
