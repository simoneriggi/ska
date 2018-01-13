#!/bin/sh

NARGS="$#"
echo "INFO: NARGS= $NARGS"

if [ "$NARGS" -lt 2 ]; then
	echo "ERROR: Invalid number of arguments...see script usage!"
  echo ""
	echo "**************************"
  echo "***     USAGE          ***"
	echo "**************************"
 	echo "$0 [ARGS]"
	echo ""
	echo "=========================="
	echo "==    ARGUMENT LIST     =="
	echo "=========================="
	echo "*** MANDATORY ARGS ***"
	echo "--fileext=[FILE_EXT] - File extension (e.g. txt, root, fits) to be placed in list"
	echo ""
	echo ""
	echo "*** OPTIONAL ARGS ***"	
	echo "--rootdir=[ROOT_DIR] - Directory where to start searching file to be placed in list [default=pwd]"
	echo "--fileprefix=[FILE_PREFIX] - File prefix filter [default: none]"
	echo "--stripext - Strip file extension before writing to list? [default=no]"
	echo "--strippath - Strip file path before writing to list? [default=no]"
	echo "--recursive - Search recursively down from ROOT_DIR [default=no]"
	echo "--output=[OUTPUTFILE] - Output file name with file list [default=filelist.txt]"
	echo "=========================="
	exit 1
fi

#######################################
##         PARSE ARGS
#######################################
FILE_EXT=""
FILE_EXT_GIVEN=false
FILE_PREFIX=""
STRIP_EXT=false
STRIP_PATH=false
RECURSIVE=false
ROOT_DIR="$PWD"
OUTPUTFILE="filelist.txt"

for item in $*
do
	case $item in 
		## MANDATORY ##	
		--fileext=*)
    	FILE_EXT=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			if [ "$FILE_EXT" != "" ]; then
				FILE_EXT_GIVEN=true
			fi
    ;;
		## OPTIONAL ##	
		--rootdir=*)
    	ROOT_DIR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;	
		--fileprefix=*)
    	FILE_PREFIX=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--output=*)
    	OUTPUTFILE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		
		--stripext*)
    	STRIP_EXT=true
    ;;
		--strippath*)
    	STRIP_PATH=true
    ;;
		--recursive*)
    	RECURSIVE=true
    ;;

    *)
    # Unknown option
    echo "ERROR: Unknown option ($item)...exit!"
    exit 1
    ;;
	esac
done

CURRENTDIR=$PWD

echo "$CURRENTDIR"

if [ "$RECURSIVE" = true ]; then
	# Recursive search
	for file in $ROOT_DIR/$FILE_PREFIX*.$FILE_EXTENSION $ROOT_DIR/**/$FILE_PREFIX*.$FILE_EXTENSION ; do
		
		# Strip path?
		if [ "$STRIP_PATH" = true ]; then
			file_base=$(basename $file)
			file=$file_base
		fi
		
		# Strip extension?
		if [ "$STRIP_EXT" = true ]; then
			file_noExtension="${file//.$FILE_EXTENSION/}"
			file=$file_noExtension
		fi

  	echo $file
  	echo $file >> $CURRENTDIR/$OUTPUTFILE
	done

else
	# Normal search
	for file in $ROOT_DIR/$FILE_PREFIX*.$FILE_EXTENSION ; do
 		# Strip path?
		if [ "$STRIP_PATH" = true ]; then
			file_base=$(basename $file)
			file=$file_base
		fi
		
		# Strip extension?
		if [ "$STRIP_EXT" = true ]; then
			file_noExtension="${file//.$FILE_EXTENSION/}"
			file=$file_noExtension
		fi

  	echo $file
  	echo $file >> $CURRENTDIR/$OUTPUTFILE
	done
fi


