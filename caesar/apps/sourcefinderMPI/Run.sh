#!/bin/sh

echo "INFO: NARGS= $#"

CURRENTDIR="$PWD"
## Check if CAESAR environment variable exists and is properly set
##if [ "$CAESAR_DIR" = "" ]; then
##  echo "ERROR: Missing environment variable CAESAR_DIR --> Please set it to Caesar installation path!"
##  exit 1
##fi

EXE="$CAESAR_DIR/bin/FindSourceMPI"

NARGS="$#"

if [ "$NARGS" -lt 2 ]; then
	echo "ERROR: Invalid number of arguments...see script usage!"
  echo ""
  echo "**** USAGE ****"
 	echo "$0 --nproc=[NPROCESSOR] --config=[CONFIGFILE]"
  echo "****************"
  exit 1
fi


## Parse args
CONFIGFILE=""
NPROC=1

for i in $*
do
	case $i in 	
    --config=*)
    	CONFIGFILE=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--nproc=*)
      NPROC=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
    *)
    # unknown option
    echo "ERROR: Unknown option...exit!"
    exit 1
    ;;
	esac
done

echo ""
echo "*****  PARSED ARGUMENTS ****"
echo "CONFIGFILE= $CONFIGFILE"
echo "NPROC= $NPROC"
echo "****************************"


CMD="mpirun -np $NPROC $EXE --config=$CONFIGFILE"
echo "INFO: Running cmd: $CMD"

eval $CMD

