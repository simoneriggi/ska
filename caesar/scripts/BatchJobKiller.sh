#!/bin/sh

STARTJOBID=$1
ENDJOBID=$2
STEP=1
echo "job start id $STARTJOBID"
echo "job end id $ENDJOBID"

for ((jobid=$STARTJOBID; jobid<=$ENDJOBID; jobid=$jobid+$STEP))
  do
    echo "$current job id $jobid"
    qdel $jobid

  done

echo "END"
