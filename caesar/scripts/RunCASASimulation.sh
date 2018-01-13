#!/bin/bash

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
	echo "--nruns=[NRUNS] - Number of simulation runs"
	echo "--envfile=[ENV_FILE] - File (.sh) with list of environment variables to be loaded by each processing node"
	echo "--mapsize=[MAP_SIZE] - Map size (in pixels)"
	echo "--pixsize=[PIX_SIZE] - Pixel size (in arcsec)"
	echo ""
	echo ""
	echo "*** OPTIONAL ARGS ***"	
	echo "--startid=[START_ID] - Run start id (default: 1)"
	echo "--sourcegenmargin=[SOURCE_GEN_MARGIN_SIZE] - Left/right margin in skymodel map for source generation (default: 0)"
	echo "--bmaj=[BMAJ] - Bmaj of sky model map (default: 9.8)"
	echo "--submit - Submit the script to the batch system using queue specified"
	echo "--containerrun - Run inside Caesar container"
	echo "--containerimg=[CONTAINER_IMG] - Singularity container image file (.simg) with CAESAR installed software"
	echo "--queue=[BATCH_QUEUE] - Name of queue in batch system" 
	echo "=========================="
	exit 1
fi



#######################################
##         PARSE ARGS
#######################################
NRUNS=1
NRUNS_GIVEN=false
START_ID=1
SUBMIT=false
BATCH_QUEUE=""
ENV_FILE=""
MAP_SIZE=""
MAP_SIZE_GIVEN=false
PIX_SIZE=""
PIX_SIZE_GIVEN=false
SOURCE_GEN_MARGIN_SIZE=0
GEN_SOURCES=true
GEN_EXT_SOURCES=false
CONTAINER_IMG=""
RUN_IN_CONTAINER=false

for item in $*
do
	case $item in 
		## MANDATORY ##	
		--nruns=*)
    	NRUNS=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			if [ "$NRUNS" != "" ]; then
				NRUNS_GIVEN=true
			fi
    ;;
		--mapsize=*)
    	MAP_SIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			if [ "$MAP_SIZE" != "" ]; then
				MAP_SIZE_GIVEN=true
			fi
    ;;
		--pixsize=*)
    	PIX_SIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			if [ "$PIX_SIZE" != "" ]; then
				PIX_SIZE_GIVEN=true
			fi
    ;;
		--envfile=*)
    	ENV_FILE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		

		## OPTIONAL ##	
		--startid=*)
    	START_ID=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;	
		--sourcegenmargin=*)
    	SOURCE_GEN_MARGIN_SIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--sources*)
    	GEN_SOURCES=true
    ;;
		--no-sources*)
    	GEN_SOURCES=false
    ;;
		--extsources*)
    	GEN_EXT_SOURCES=true
    ;;
		--no-extsources*)
    	GEN_EXT_SOURCES=false
    ;;

		--queue=*)
    	BATCH_QUEUE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--submit*)
    	SUBMIT=true
    ;;
		--containerimg=*)
    	CONTAINER_IMG=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--containerrun*)
    	RUN_IN_CONTAINER=true
    ;;
		

    *)
    # Unknown option
    echo "ERROR: Unknown option ($item)...exit!"
    exit 1
    ;;
	esac
done


echo ""
echo "*****  PARSED ARGUMENTS ****"
echo "SUBMIT? $SUBMIT, QUEUE=$BATCH_QUEUE"
echo "NRUNS: $NRUNS (START_ID=$START_ID)"
echo "ENV_FILE: $ENV_FILE"
echo "RUN_IN_CONTAINER? $RUN_IN_CONTAINER, CONTAINER_IMG=$CONTAINER_IMG"
echo "MAP_SIZE: $MAP_SIZE (pixels), PIX_SIZE=$PIX_SIZE (arcsec)"
echo "SOURCE_GEN_MARGIN_SIZE: $SOURCE_GEN_MARGIN_SIZE (pixels)"
echo "GEN_SOURCES? $GEN_EXT_SOURCES, GEN_EXT_SOURCES? $GEN_EXT_SOURCES"
echo "****************************"
echo ""

## Check arguments parsed
if [ "$NRUNS_GIVEN" = false ] ; then
  echo "ERROR: Missing or empty NRUNS args (hint: you should specify how many simulation to be performed)!"
  exit 1
fi

if [ "$MAP_SIZE_GIVEN" = false ] ; then
  echo "ERROR: Missing or empty MAP_SIZE args (hint: you should specify how big is your skymodel & sim map in pixels)!"
  exit 1
fi

if [ "$PIX_SIZE_GIVEN" = false ] ; then
  echo "ERROR: Missing or empty PIX_SIZE args (hint: you should specify how big is your skymodel & sim pixel size in arcsec)!"
  exit 1
fi

if [ "$BATCH_QUEUE" = "" ] && [ "$SUBMIT" = true ]; then
  echo "ERROR: Empty BATCH_QUEUE argument (hint: you must specify a queue if submit option is activated)!"
  exit 1
fi

if [ "$ENV_FILE" = "" ]; then
  echo "ERROR: Empty ENV_FILE arg!"
  exit 1
fi

if [ "$CONTAINER_IMG" = "" ] && [ "$RUN_IN_CONTAINER" = true ]; then
  echo "ERROR: Empty CONTAINER_IMG argument (hint: you must specify a container image if run in container option is activated)!"
  exit 1
fi

## Define flags
GEN_SOURCE_FLAG="--compactsources"
GEN_EXT_SOURCE_FLAG="--extsources"
if [ "$GEN_SOURCES" = false ] ; then
	GEN_SOURCE_FLAG="--no-compactsources"
fi
if [ "$GEN_EXT_SOURCES" = false ] ; then
	GEN_EXT_SOURCE_FLAG="--no-extsources"
fi

#######################################
##     DEFINE & LOAD ENV VARS
#######################################
export BASEDIR="$PWD"
export OUTPUT_DATADIR="$PWD"
export DATADIR=""

## Load env file
echo "INFO: Loading environment variables defined in file $ENV_FILE ..."
source "$ENV_FILE"

###################################################
###     GENERATE SUBMISSION SCRIPT
###################################################
RUN_ID=$START_ID
echo "Generate submission scripts (start run numbering from RUN_ID=$RUN_ID)"

for ((index=1; index<=$NRUNS; index=$index+1))
	do

	## Create job top directory
	JOB_DIR="$BASEDIR/RUN$RUN_ID"
	CASA_SIM_DIR="$JOB_DIR/sim"
	echo "INFO: Creating job top directory $JOB_DIR ..."
	mkdir -p "$JOB_DIR"

	echo "INFO: Creating sim directory $CASA_SIM_DIR ..."
	mkdir -p "$CASA_SIM_DIR"

	## Define output files
  #simoutfile='SimInfo-MuonWritingAtCenterScenario_Mu_n'"$Nev"'-RUN'"$index"
  #echo $simoutfile

	

  echo ""

	## Generate script
	shfile="Sim-RUN$RUN_ID.sh"
	echo "*** Creating sh file $shfile ***"
	(
		echo "#!/bin/bash"
		echo "#PBS -o $BASEDIR"
    echo "#PBS -o $BASEDIR"
    echo '#PBS -r n'
    echo '#PBS -S /bin/sh'
    echo "#PBS -N SimJob$RUN_ID"
    echo '#PBS -p 1'

    echo " "
    echo " "

    echo 'echo "*************************************************"'
    echo 'echo "****         PREPARE JOB                     ****"'
    echo 'echo "*************************************************"'
		echo "export JOBDIR=$JOB_DIR" 
		echo 'echo "INFO: Entering job directory $JOBDIR ..."'
		echo 'cd $JOBDIR'
    echo 'echo "INFO: Source the software environment ..."'
    echo "source $ENV_FILE"
		echo 'echo ""'
	
		echo " "
    echo " "

    echo 'echo "*************************************************"'
    echo 'echo "****         RUN SKYMODEL SIMULATION         ****"'
    echo 'echo "*************************************************"'
		
		if [ "$RUN_IN_CONTAINER" = true ] ; then
			echo "EXE=singularity run --app skymodelgenerator $CONTAINER_IMG"
		else
			echo "EXE=$CAESAR_SCRIPTS_DIR/map_simulator.py"
		fi

		echo 'EXE_ARGS="'"--nx=$MAP_SIZE --ny=$MAP_SIZE --pixsize=$PIX_SIZE --marginx=$SOURCE_GEN_MARGIN_SIZE --marginy=$SOURCE_GEN_MARGIN_SIZE $GEN_SOURCE_FLAG $GEN_EXT_SOURCE_FLAG"'"'
#--bmaj=9.8 --bmin=5.8 --bpa=-3 --crpix1=250 --crpix2=250 --bkg --bkg_level=10e-6 --bkg_rms=100e-6 --compactsources --zmin=9999 --zmax=10000 --source_density=1000 --no-extsources --zmin_model=1 --outputfile=sim_map.fits --outputfile_model=skymodel.fits --outputfile_sources=sources.root --outputfile_ds9region=dsregions.reg --outputfile_casaregion=casa_mask.dat 
		echo 'echo "Running command $EXE $EXE_ARGS"'
		echo '$EXE $EXE_ARGS'

		echo 'echo ""'

	 	echo " "
    echo " "
    
		echo 'echo "*************************************************"'
    echo 'echo "****         RUN CASA SIMULATION             ****"'
    echo 'echo "*************************************************"'
		echo 'echo ""'
    echo 'cd $JOBDIR'

		echo 'EXE="$CASAPATH/bin/casa --nologger --log2term --nogui -c $CAESAR_SCRIPTS_DIR/simulate_observation.py"'
		echo 'EXE_ARGS="'"--vis=$VIS --skymodel=$SKYMODEL"'"'
		#$CASAPATH/bin/casa --nologger --log2term --nogui -c $CAESAR_SCRIPTS_DIR/simulate_observation.py --vis=vis.ms --skymodel=skymodel.fits
		echo 'echo "Running command $EXE $EXE_ARGS"'
		echo '$EXE $EXE_ARGS'
    
    echo 'echo ""'
    
    echo 'echo "*** END RUN ***"'
				
	) > $shfile
	chmod +x $shfile

	####mv $shfile $CURRENTJOBDIR
	(( RUN_ID= $RUN_ID + 1 ))
	
	
	# Submits the job to batch system
	if [ "$SUBMIT" = true ] ; then
		echo "INFO: Submitting script $shfile to QUEUE $BATCH_QUEUE ..."
		qsub -q $BATCH_QUEUE $BASEDIR/$shfile
	fi

	
done

echo "*** END SUBMISSION ***"

