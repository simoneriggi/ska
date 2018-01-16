#!/bin/bash

NARGS="$#"
if [ "$NARGS" -lt 2 ]; then
	echo "ERROR: Invalid number of arguments...see script usage!"
	echo ""
	echo "**************************"
  echo "***     USAGE          ***"
	echo "**************************"
 	echo "$0 [FITS FILE] [REGION_FILE1] [REGION_FILE2] ..."
	echo ""
	echo "=========================="
  exit 1
fi


FILENAME="$1"
shift
REGION_FILES="$@"
echo "INFO: FILENAME: $FILENAME"
echo "INFO: REGION_FILES: $REGION_FILES"


## Drawing FITS
ds9 $FILENAME &

for item in "$REGION_FILES"
do
	echo "Loading regions present in file $item ..."
	ds9 -regions load $item &
done
