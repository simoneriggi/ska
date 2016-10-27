// ***********************************************************************
// * License and Disclaimer                                              *
// *                                                                     *
// * Copyright 2016 Simone Riggi																			   *
// *																																	   *
// * This file is part of Caesar. 																		   *
// * Caesar is free software: you can redistribute it and/or modify it   *
// * under the terms of the GNU General Public License as published by   *
// * the Free Software Foundation, either * version 3 of the License,    *
// * or (at your option) any later version.                              *
// * Caesar is distributed in the hope that it will be useful, but 			 *
// * WITHOUT ANY WARRANTY; without even the implied warranty of          * 
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                *
// * See the GNU General Public License for more details. You should     * 
// * have received a copy of the GNU General Public License along with   * 
// * Caesar. If not, see http://www.gnu.org/licenses/.                   *
// ***********************************************************************
/**
* @file SourceFinderMPI.cc
* @class SourceFinderMPI
* @brief Source finderMPI class
*
* Class to perform source finding in MPI
* @author S. Riggi
* @date 20/01/2015
*/

#include <SourceFinderMPI.h>
#include <Serializer.h>

//Caesar headers
#include <BlobFinder.h>
#include <Img.h>
#include <Source.h>
#include <Contour.h>
#include <ConfigParser.h>
#include <BkgData.h>
#include <CodeUtils.h>
#include <Logger.h>
#include <MathUtils.h>
#include <SLIC.h>
#include <SLICUtils.h>
#include <SLICSegmenter.h>

#include <TObject.h>
#include <TMatrixD.h>
#include <TFile.h>
#include <TCanvas.h>

//MPI header
#include <mpi.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <numeric>
#include <time.h>
#include <ctime>
#include <queue>

using namespace std;


namespace Caesar {

#define MASTER_ID 0

SourceFinderMPI::SourceFinderMPI() {

	//## Init MPI comm size & rank	
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

}//close costructor


SourceFinderMPI::~SourceFinderMPI(){
	
	//Delete task data
	for(unsigned int i=0;i<m_taskDataPerWorkers.size();i++) {
		for(unsigned int j=0;j<m_taskDataPerWorkers[i].size();j++) {
			if(m_taskDataPerWorkers[i][j]){
				delete m_taskDataPerWorkers[i][j];
				m_taskDataPerWorkers[i][j]= 0;
			}
		}
		m_taskDataPerWorkers[i].clear();
	}
	m_taskDataPerWorkers.clear();

	//Clear file
	if(myid==0){
		//if(m_SourceTree) m_SourceTree->Delete();
		if(m_OutputFile) m_OutputFile->Close();
		if(m_DS9CatalogFilePtr) fclose(m_DS9CatalogFilePtr);
		//if(m_Application) m_Application->Delete();
	}

}//close destructor


void SourceFinderMPI::InitOptions(){

	//Input file options
	m_InputFileName= "";
	m_InputImgName= "";
	m_InputFileExtension= "";
	m_InputImg= 0;

	//Output options
	m_OutputFile= 0;	
	m_OutputFileName= "";
	m_Application= 0;
	m_IsInteractiveRun= false;
	m_SaveToFile= true;
	m_SaveConfig= true;
	m_SaveSources= true;	
	m_SourceTree= 0;
	m_SaveDS9Region= true;
	m_DS9CatalogFileName= "";
	m_DS9CatalogFilePtr= 0;
	m_DS9RegionFormat= 1;
		
	//Source 
	m_Source= 0;
	m_SourceCollection.clear();
	m_CompactSources.clear();
	m_ExtendedSources.clear();

	//Read options
	m_ReadTile= false;
	m_TileMinX= 0;
	m_TileMaxX= 0;
	m_TileMinY= 0;
	m_TileMaxY= 0;
	
	//Bkg data
	m_BkgData= 0;
	m_SignificanceMap= 0;

	//Residual img
	m_ResidualImg= 0;

	//Saliency img
	m_SaliencyImg= 0;

}//close InitOptions()

int SourceFinderMPI::Init(){

	//## Init options
	InitOptions();

	//## Configure from parser
	if(Configure()<0){
		ERROR_LOG("Failed to configure options from parser!");
		return -1;
	}

	//## Create TApplication if interactive run is selected
	if(myid==0) {
		if(!m_Application && m_IsInteractiveRun){
			m_Application= new TApplication("Application", 0, 0);
		}	
	}

	//## Create output file
	if(myid==0) {
		if(m_SaveToFile){
			m_OutputFile= new TFile(m_OutputFileName.c_str(),"RECREATE");	
			m_OutputFile->cd();
	
			//Init source tree
			if(m_SaveSources){
				m_Source= 0;
				if(!m_SourceTree) m_SourceTree= new TTree("SourceInfo","SourceInfo");
				m_SourceTree->Branch("Source",&m_Source);
				m_SourceCollection.clear();
			}
	
			//Init DS9 catalog
			if(!m_DS9CatalogFilePtr) m_DS9CatalogFilePtr= fopen(m_DS9CatalogFileName.c_str(),"w");

		}//close if saveToFile
	}//close if nproc=0

	//## Init and fill task data (done by all processors)
	for(int i=0;i<nproc;i++){
		m_taskDataPerWorkers.push_back( std::vector<TaskData*>() );
	}
	if(PrepareWorkerTasks()<0){
		ERROR_LOG("Preparation of tasks per worker failed!");
		return -1;
	}
	
	return 0;

}//close Init()

int SourceFinderMPI::Configure(){

	//Get image read options
	if(GET_OPTION_VALUE(inputFile,m_InputFileName)<0){
		ERROR_LOG("Failed to get inputFile option!");
		return -1;
	}
	GET_OPTION_VALUE(inputImage,m_InputImgName);
	GET_OPTION_VALUE(readTileImage,m_ReadTile);
	if(m_ReadTile){
		GET_OPTION_VALUE(tileMinX,m_TileMinX);
		GET_OPTION_VALUE(tileMaxX,m_TileMaxX);
		GET_OPTION_VALUE(tileMinY,m_TileMinY);
		GET_OPTION_VALUE(tileMaxY,m_TileMaxY);
	}

	//Get distributed image options
	GET_OPTION_VALUE(tileSizeX,m_TileSizeX);
	GET_OPTION_VALUE(tileSizeY,m_TileSizeY);
	GET_OPTION_VALUE(useTileOverlap,m_UseTileOverlap);
	GET_OPTION_VALUE(tileStepSizeX,m_TileStepSizeX);
	GET_OPTION_VALUE(tileStepSizeY,m_TileStepSizeY);
	
	if(m_TileSizeX<=0 || m_TileSizeY<=0) {
		ERROR_LOG("Invalid tileSizeX/tileSizeY options!");
		return -1;	
	}
	if(m_TileStepSizeX<=0 || m_TileStepSizeY<=0 || m_TileStepSizeX>1 || m_TileStepSizeY>1){
		ERROR_LOG("Invalid tileStepSizeX/tileStepSizeY options!");
		return -1;
	}
	if(!m_UseTileOverlap){
		m_TileStepSizeX= 1;
		m_TileStepSizeY= 1;	
	}

	//Get output file options
	GET_OPTION_VALUE(outputFile,m_OutputFileName);
	GET_OPTION_VALUE(saveToFile,m_SaveToFile);
	GET_OPTION_VALUE(saveConfig,m_SaveConfig);
	GET_OPTION_VALUE(saveDS9Region,m_SaveDS9Region);
	GET_OPTION_VALUE(ds9RegionFile,m_DS9CatalogFileName);
	GET_OPTION_VALUE(DS9RegionFormat,m_DS9RegionFormat);
	GET_OPTION_VALUE(saveSources,m_SaveSources);
	GET_OPTION_VALUE(isInteractiveRun,m_IsInteractiveRun);
	GET_OPTION_VALUE(saveResidualMap,m_SaveResidualMap);
	
	//Get bkg options
	GET_OPTION_VALUE(useLocalBkg,m_UseLocalBkg);
	GET_OPTION_VALUE(use2ndPassInLocalBkg,m_Use2ndPassInLocalBkg);
	GET_OPTION_VALUE(skipOutliersInLocalBkg,m_SkipOutliersInLocalBkg);
	GET_OPTION_VALUE(localBkgMethod,m_LocalBkgMethod);
	GET_OPTION_VALUE(bkgEstimator,m_BkgEstimator);
	GET_OPTION_VALUE(boxSizeX,m_BoxSizeX);
	GET_OPTION_VALUE(boxSizeY,m_BoxSizeY);
	GET_OPTION_VALUE(gridSizeX,m_GridSizeX);
	GET_OPTION_VALUE(gridSizeY,m_GridSizeY);
	GET_OPTION_VALUE(useBeamInfoInBkg,m_UseBeamInfoInBkg);
	
	//Get source search options
	GET_OPTION_VALUE(searchCompactSources,m_SearchCompactSources);
	GET_OPTION_VALUE(minNPix,m_NMinPix);
	GET_OPTION_VALUE(seedBrightThr,m_SeedBrightThr);	
	GET_OPTION_VALUE(seedThr,m_SeedThr);
	GET_OPTION_VALUE(mergeThr,m_MergeThr);

	GET_OPTION_VALUE(mergeBelowSeed,m_MergeBelowSeed);
	GET_OPTION_VALUE(searchNegativeExcess,m_SearchNegativeExcess);
	GET_OPTION_VALUE(searchNestedSources,m_SearchNestedSources);
	GET_OPTION_VALUE(nestedBlobThrFactor,m_NestedBlobThrFactor);
	GET_OPTION_VALUE(applySourceSelection,m_ApplySourceSelection);
	
	//Get source selection options
	GET_OPTION_VALUE(applySourceSelection,m_ApplySourceSelection);
	GET_OPTION_VALUE(sourceMinBoundingBox,m_SourceMinBoundingBox);
	GET_OPTION_VALUE(psCircRatioThr,m_psCircRatioThr);
	GET_OPTION_VALUE(psElongThr,m_psElongThr);
	GET_OPTION_VALUE(psEllipseAreaRatioMinThr,m_psEllipseAreaRatioMinThr);
	GET_OPTION_VALUE(psEllipseAreaRatioMaxThr,m_psEllipseAreaRatioMaxThr);
	GET_OPTION_VALUE(psMaxNPix,m_psMaxNPix);
	
	//Get source residual options
	GET_OPTION_VALUE(dilateNestedSources,m_DilateNestedSources);
	GET_OPTION_VALUE(dilateKernelSize,m_DilateKernelSize);
	GET_OPTION_VALUE(dilatedSourceType,m_DilatedSourceType);
	GET_OPTION_VALUE(dilateSourceModel,m_DilateSourceModel);
	GET_OPTION_VALUE(dilateRandomize,m_DilateRandomize);
	
	//Get smoothing options
	GET_OPTION_VALUE(usePreSmoothing,m_UsePreSmoothing);
	GET_OPTION_VALUE(smoothFilter,m_SmoothFilter);
	GET_OPTION_VALUE(gausFilterKernSize,m_GausFilterKernSize);
	GET_OPTION_VALUE(gausFilterSigma,m_GausFilterSigma);
	GET_OPTION_VALUE(guidedFilterRadius,m_GuidedFilterRadius);
	GET_OPTION_VALUE(guidedFilterColorEps,m_GuidedFilterColorEps);
	
	//Get saliency options
	GET_OPTION_VALUE(saliencyThrFactor,m_SaliencyThrFactor);
	GET_OPTION_VALUE(saliencyBkgThrFactor,m_SaliencyBkgThrFactor);
	GET_OPTION_VALUE(saliencyImgThrFactor,m_SaliencyImgThrFactor);
	GET_OPTION_VALUE(saliencyResoMin,m_SaliencyResoMin);
	GET_OPTION_VALUE(saliencyResoMax,m_SaliencyResoMax);
	GET_OPTION_VALUE(saliencyResoStep,m_SaliencyResoStep);
	GET_OPTION_VALUE(saliencyUseRobustPars,m_SaliencyUseRobustPars);
	GET_OPTION_VALUE(saliencyUseBkgMap,m_SaliencyUseBkgMap);
	GET_OPTION_VALUE(saliencyUseNoiseMap,m_SaliencyUseNoiseMap);
	GET_OPTION_VALUE(saliencyUseCurvInDiss,m_SaliencyUseCurvInDiss);
	GET_OPTION_VALUE(saliencyNNFactor,m_SaliencyNNFactor);
	GET_OPTION_VALUE(saliencySpatialRegFactor,m_SaliencySpatialRegFactor);
	GET_OPTION_VALUE(saliencyMultiResoCombThrFactor,m_SaliencyMultiResoCombThrFactor);

	//Get extended source options
	GET_OPTION_VALUE(searchExtendedSources,m_SearchExtendedSources);
	GET_OPTION_VALUE(extendedSearchMethod,m_ExtendedSearchMethod);
	GET_OPTION_VALUE(wtScaleExtended,m_wtScaleExtended);
			
	//Get superpixel options
	GET_OPTION_VALUE(spSize,m_spSize);
	GET_OPTION_VALUE(spBeta,m_spBeta);
	GET_OPTION_VALUE(spMinArea,m_spMinArea);
	GET_OPTION_VALUE(spUseLogContrast,m_spUseLogContrast);

	//Chan-Vese options
	GET_OPTION_VALUE(cvTimeStepPar,m_cvTimeStepPar);
	GET_OPTION_VALUE(cvWindowSizePar,m_cvWindowSizePar);
	GET_OPTION_VALUE(cvLambda1Par,m_cvLambda1Par);
	GET_OPTION_VALUE(cvLambda2Par,m_cvLambda2Par);
	GET_OPTION_VALUE(cvMuPar,m_cvMuPar);
	GET_OPTION_VALUE(cvNuPar,m_cvNuPar);
	GET_OPTION_VALUE(cvPPar,m_cvPPar);
		
	return 0;

}//close Configure()


int SourceFinderMPI::PrepareWorkerTasks(){

	//## Generate a uuid for this job
	std::string jobId= CodeUtils::GenerateUUID();
	DEBUG_LOG("Generated jobId: "<<jobId);
	
	//## Get input image size
	FileInfo info;
	bool match_extension= false;
	if(!SysUtils::CheckFile(m_InputFileName,info,match_extension,"")){
		ERROR_LOG("Invalid input file name specified (filename="<<m_InputFileName<<"), invalid file path?!");
		return -1;
	}
	
	long int Nx= -1;
	long int Ny= -1;

	if(info.extension==".root"){// ROOT format
		if(m_ReadTile){
			Nx= m_TileMaxX-m_TileMinX+1;
			Ny= m_TileMaxY-m_TileMinY+1;
		}
		else{//READ FULL MAP
			TFile* inputFile= new TFile(m_InputFileName.c_str(),"READ");	
			if(!inputFile || inputFile->IsZombie()){
				ERROR_LOG("Failed to open input file image "<<m_InputFileName<<" and get image size!");
				return -1;
			}
			Img* inputImg= (Img*)inputFile->Get(m_InputImgName.c_str());
			if(!inputImg) {
				ERROR_LOG("Failed to open input file image "<<m_InputFileName<<" and get image size!");
				return -1;	
			}
			Nx= inputImg->GetNbinsX();
			Ny= inputImg->GetNbinsY();
		}
	}//close if
	else if(info.extension==".fits"){//FITS
		if(m_ReadTile){
			Nx= m_TileMaxX-m_TileMinX+1;
			Ny= m_TileMaxY-m_TileMinY+1;
		}
		else{//READ FULL MAP
			if(SysUtils::GetFITSImageSize(m_InputFileName,Nx,Ny)<0){
				ERROR_LOG("Failed to open input file image "<<m_InputFileName<<" and get image size!");
				return -1;
			}
			INFO_LOG("Input image opened with success: size="<<Nx<<"x"<<Ny);
		}
	}//close else if		
	else{
		ERROR_LOG("Invalid/unsupported file extension ("<<info.extension<<") detected!");
		return -1;
	}

	INFO_LOG("Image size: "<<Nx<<"x"<<Ny);

	//## Compute the 2D grid 
	std::vector<long int> ix_min;
	std::vector<long int> ix_max;
	std::vector<long int> iy_min;
	std::vector<long int> iy_max;

	long int tileOverlapX= 0;
	long int tileOverlapY= 0;
	if(m_UseTileOverlap){
		tileOverlapX= m_TileStepSizeX;
		tileOverlapY= m_TileStepSizeY;
	}

	INFO_LOG("Computing tile partition: tileOverlap("<<tileOverlapX<<","<<tileOverlapY<<")");
	if(MathUtils::Compute2DGrid(ix_min,ix_max,iy_min,iy_max,Nx,Ny,m_TileSizeX,m_TileSizeY,m_TileStepSizeX,m_TileStepSizeY)<0){
		WARN_LOG("Failed to compute a 2D image partition!");
		return -1;
	}
	int nExpectedTasks= ix_min.size()*iy_min.size();
	INFO_LOG("#"<<nExpectedTasks<<" expected number of tasks ("<<ix_min.size()<<"x"<<iy_min.size()<<")");

	//## Compute worker tasks (check max number of tasks per worker)
	INFO_LOG("Computing worker task list...");
	TaskData* aTaskData= 0;
	long int workerCounter= 0;

	for(unsigned int j=0;j<iy_min.size();j++){
		for(unsigned int i=0;i<ix_min.size();i++){
			//Assign worker
			DEBUG_LOG("Assign task ("<<i<<","<<j<<") to worker no. "<<workerCounter<<"...");
				
			aTaskData= new TaskData;
			aTaskData->filename= m_InputFileName;
			aTaskData->jobId= jobId;
			aTaskData->workerId= workerCounter;
			aTaskData->ix_min= ix_min[i] + m_TileMinX;
			aTaskData->ix_max= ix_max[i] + m_TileMinX;
			aTaskData->iy_min= iy_min[j] + m_TileMinY;
			aTaskData->iy_max= iy_max[j] + m_TileMinY;
			m_taskDataPerWorkers[workerCounter].push_back(aTaskData);

			if(workerCounter>=nproc-1) workerCounter= 0;
			else workerCounter++;
		}//end loop x
	}//end loop y
	
	//Fill neighbor task list
	for(unsigned int i=0;i<m_taskDataPerWorkers.size();i++){

		if(m_taskDataPerWorkers[i].size()==0) continue;//no tasks present

		for(unsigned int j=0;j<m_taskDataPerWorkers[i].size()-1;j++){
			long int ix_min= m_taskDataPerWorkers[i][j]->ix_min;
			long int ix_max= m_taskDataPerWorkers[i][j]->ix_max;
			long int iy_min= m_taskDataPerWorkers[i][j]->iy_min;
			long int iy_max= m_taskDataPerWorkers[i][j]->iy_max;

			for(unsigned int k=j+1;k<m_taskDataPerWorkers[i].size();k++){
				if(j==k) continue;
				long int next_ix_min= m_taskDataPerWorkers[i][k]->ix_min;
				long int next_ix_max= m_taskDataPerWorkers[i][k]->ix_max;
				long int next_iy_min= m_taskDataPerWorkers[i][k]->iy_min;
				long int next_iy_max= m_taskDataPerWorkers[i][k]->iy_max;
				
				//bool isNeighborInX= ( ix_max>next_ix_min || ix_min<next_ix_max || ix_max==next_ix_min-1 || ix_min==next_ix_max+1 ) && (next_iy_min==iy_min && next_iy_max==iy_max);
				//bool isNeighborInY= ( iy_max>next_iy_min || iy_min<next_iy_max || iy_max==next_iy_min-1 || iy_min==next_iy_max+1 ) && (next_ix_min==ix_min && next_ix_max==ix_max);
				bool isNeighborInX= ( ix_max>next_ix_min || ix_max==next_ix_min-1 || ix_min==next_ix_max+1 ) && (next_iy_min==iy_min && next_iy_max==iy_max);
				bool isNeighborInY= ( iy_max>next_iy_min || iy_max==next_iy_min-1 || iy_min==next_iy_max+1 ) && (next_ix_min==ix_min && next_ix_max==ix_max);

				if(isNeighborInX || isNeighborInY) {
					(m_taskDataPerWorkers[i][j]->neighborTaskId).push_back(k);
					(m_taskDataPerWorkers[i][k]->neighborTaskId).push_back(j);
					(m_taskDataPerWorkers[i][j]->neighborWorkerId).push_back(i);
					(m_taskDataPerWorkers[i][k]->neighborWorkerId).push_back(i);
				}

			}//end loop next task
		}//end loop tasks
	}//end loop workers

	//Print
	for(unsigned int i=0;i<m_taskDataPerWorkers.size();i++){
		if(m_taskDataPerWorkers[i].size()==0) continue;//no tasks present

		INFO_LOG("== Worker no. "<<i<<"==");
		for(unsigned int j=0;j<m_taskDataPerWorkers[i].size()-1;j++){
			long int ix_min= m_taskDataPerWorkers[i][j]->ix_min;
			long int ix_max= m_taskDataPerWorkers[i][j]->ix_max;
			long int iy_min= m_taskDataPerWorkers[i][j]->iy_min;
			long int iy_max= m_taskDataPerWorkers[i][j]->iy_max;
			INFO_LOG("Task no. "<<j<<"["<<ix_min<<","<<ix_max<<"] ["<<iy_min<<","<<iy_max<<"] neighbors{");	
			for(unsigned int k=0;k<m_taskDataPerWorkers[i][j]->neighborTaskId.size();k++){
				long int neighborWorkerId= m_taskDataPerWorkers[i][j]->neighborWorkerId[k];
				long int neighborTaskId= m_taskDataPerWorkers[i][j]->neighborTaskId[k];
				long int next_ix_min= m_taskDataPerWorkers[neighborWorkerId][neighborTaskId]->ix_min;
				long int next_ix_max= m_taskDataPerWorkers[neighborWorkerId][neighborTaskId]->ix_max;
				long int next_iy_min= m_taskDataPerWorkers[neighborWorkerId][neighborTaskId]->iy_min;
				long int next_iy_max= m_taskDataPerWorkers[neighborWorkerId][neighborTaskId]->iy_max;

				INFO_LOG("("<<neighborWorkerId<<","<<neighborTaskId<<") ["<<next_ix_min<<","<<next_ix_max<<"] ["<<next_iy_min<<","<<next_iy_max<<"]");
			}
			INFO_LOG("}");	
		}
	}


	bool hasTooManyTasks= false;
	int maxNTasksPerWorker_default= 100;
	for(unsigned int i=0;i<m_taskDataPerWorkers.size();i++){
		long int nTasksPerWorker= (long int)m_taskDataPerWorkers[i].size();
		if(nTasksPerWorker>maxNTasksPerWorker_default){
			hasTooManyTasks= true;
			break;
		}
	}

	if(hasTooManyTasks){
		WARN_LOG("Too many tasks per worker (thr="<<maxNTasksPerWorker_default<<")");
		return -1;
	}

	return 0;

}//close PrepareWorkerTasks()


int SourceFinderMPI::Run(){

	double globalTimerStart= MPI_Wtime();

	//## Init options (done by all processors)
	double initTimerStart= MPI_Wtime();
	INFO_LOG("Initializing source finder...");
	if(Init()<0){
		ERROR_LOG("Initialization failed!");
		return -1;
	}
	double initTimerEnd= MPI_Wtime();
	double initCPUTime= initTimerEnd-initTimerStart;

	//## Start loop on tasks per worker
	INFO_LOG("[PROC "<<myid<<"] - Start processing of #"<<m_taskDataPerWorkers[myid].size()<<" tasks...");
	
	for(unsigned int j=0;j<m_taskDataPerWorkers[myid].size();j++){
		long int ix_min= m_taskDataPerWorkers[myid][j]->ix_min;
		long int ix_max= m_taskDataPerWorkers[myid][j]->ix_max;
		long int iy_min= m_taskDataPerWorkers[myid][j]->iy_min;
		long int iy_max= m_taskDataPerWorkers[myid][j]->iy_max;

		//## Read input image
		INFO_LOG("Reading input image ["<<ix_min<<","<<ix_max<<"] ["<<iy_min<<","<<iy_max<<"]...");
		Img* taskImg= ReadImage(ix_min,ix_max,iy_min,iy_max);
		if(!taskImg){
			ERROR_LOG("Reading of input image failed, skip to next task...");
			continue;
		}

		
		//## Find compact sources
		INFO_LOG("Searching compact sources...");
		if(m_SearchCompactSources && FindCompactSources(m_taskDataPerWorkers[myid][j],taskImg)<0){
			ERROR_LOG("Compact source search failed!");
			return -1;
		}
		

		/*
		//## Find extended sources
		if(m_SearchExtendedSources){
			INFO_LOG("Searching extended sources...");

			// Find residual map
			if(FindResidualMap()<0){
				ERROR_LOG("Residual map computation failed!");
				return -1;
			}

			//Find extended sources
			if(FindExtendedSources()<0){
				ERROR_LOG("Extended source search failed!");
				return -1;
			}
		}//close if search extended sources
		*/

		if(!taskImg){
			delete taskImg;
			taskImg= 0;
		}
	}//end loop tasks per worker
	

	//## Update sources computed in each worker
	//## The updated list of task data is available in master processor
	//## (check if it is better to replace with MPI_Gather and put it available in all workers)
	if(UpdateSourceFromWorkers()<0){
		ERROR_LOG("Updating sources from workers failed!");
		return -1;
	}
	
	
	//## Deblend sources
	//if(fDeblendSources) DeblendSources(fInputImg);
	
	//## Draw & Store results (done by master processor)
	if(myid==0) {
		//## Draw final sources
		//if(m_IsInteractiveRun) DrawSources(m_InputImg,m_SourceCollection);
	
		//## Save to file
		if(m_SaveToFile) Save();	
		//if(m_Application && m_IsInteractiveRun) m_Application->Run();
	}//close if

	return 0;

}//close Run()


int SourceFinderMPI::UpdateSourceFromWorkers(){

	//## Put a barrier and collect all sources from workers in the master processor
	MPI_Barrier(MPI_COMM_WORLD);

	int MSG_TAG= 1;
	if (myid == 0) {//Receive data from the workers

		for (int i=1; i<nproc; i++) {
			//Check if this processor has tasks assigned, otherwise skip!
			if(m_taskDataPerWorkers[i].size()==0){
				INFO_LOG("No tasks assigned to worker no. "<<i<<", nothing to be collected, skip to next worker...");
				continue;
			}

			MPI_Status status;
    
			//## Probe for an incoming message from process zero
			INFO_LOG("Proc #"<<myid<<": Probing for message from proc "<<i);
    	//if(MPI_Probe(MPI_ANY_SOURCE, MSG_TAG, MPI_COMM_WORLD, &status)==MPI_SUCCESS){
			if(MPI_Probe(i, MSG_TAG, MPI_COMM_WORLD, &status)==MPI_SUCCESS){
				INFO_LOG("Proc #"<<myid<<": a message has been found with the probe, with tag " << status.MPI_TAG << ", source " << status.MPI_SOURCE);

    		//## When probe returns, the status object has the size and other
    		//## attributes of the incoming message. Get the message size
				INFO_LOG("Proc #"<<myid<<": Getting size of message... ");
    	
				int rcvMsgSize= 0;
    		MPI_Get_count(&status, MPI_CHAR, &rcvMsgSize);

				//## Allocate a buffer to hold the incoming numbers
				INFO_LOG("Proc #"<<myid<<": Allocating a message of size "<<rcvMsgSize);
				if(rcvMsgSize<=0){
					ERROR_LOG("Proc #"<<myid<<": rcvMsg size is negative!");
					continue;
				}
				char* recvBuffer= (char*)malloc(rcvMsgSize);
    		
    		//## Now receive the message with the allocated buffer
    		//MPI_Recv(recvBuffer, rcvMsgSize, MPI_CHAR, MPI_ANY_SOURCE, MSG_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(recvBuffer, rcvMsgSize, MPI_CHAR, i, MSG_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    		INFO_LOG("Proc #"<<myid<<": Received a message (msg="<<recvBuffer<<", size="<<rcvMsgSize<<") from proc "<<i);

				//## Update task data with received worker data	
				bool isTaskCollectionPreAllocated= true;
				if(Serializer::CharArrayToTaskDataCollection(m_taskDataPerWorkers[i],recvBuffer,rcvMsgSize,isTaskCollectionPreAllocated)<0 ){
					ERROR_LOG("Proc #"<<myid<<": Failed to decode recv message into task data list!");
    			if(recvBuffer) free(recvBuffer);
				}

				//## Free received buffer
    		if(recvBuffer) free(recvBuffer);
			}//close if
			else{
				ERROR_LOG("Proc #"<<myid<<": Message probing failed!");
				return -1;
				//continue;
			}
		}//end loop workers
	}//close if

	else {//Send data to master
		//## First encode taskData in protobuf
		long int msg_size= 0;
		char* msg= Serializer::TaskDataCollectionToCharArray(msg_size,m_taskDataPerWorkers[myid]);
		if(!msg){
			ERROR_LOG("Proc #"<<myid<<": Failed to encode task data to protobuf!");
			return -1;
		}

		//## Send buffer to master processor	
		INFO_LOG("Proc #"<<myid<<": Sending msg: "<<msg<<" (size="<<msg_size<<")");
		MPI_Send((void*)(msg),msg_size, MPI_CHAR, MASTER_ID, MSG_TAG, MPI_COMM_WORLD);

		//## Free buffer
		free(msg);
	}//close else

	MPI_Barrier(MPI_COMM_WORLD);

	//## Update sources in list
	if (myid == 0) {
		for(unsigned int i=0;i<m_taskDataPerWorkers.size();i++){
			for(unsigned int j=0;j<m_taskDataPerWorkers[i].size();j++){
				m_CompactSources.insert(m_CompactSources.end(),(m_taskDataPerWorkers[i][j]->sources).begin(),(m_taskDataPerWorkers[i][j]->sources).end());
				m_ExtendedSources.insert(m_ExtendedSources.end(),(m_taskDataPerWorkers[i][j]->ext_sources).begin(),(m_taskDataPerWorkers[i][j]->ext_sources).end());
			}//end loop tasks
		}//end loop workers
		m_SourceCollection.insert(m_SourceCollection.end(),m_CompactSources.begin(),m_CompactSources.end());
		m_SourceCollection.insert(m_SourceCollection.end(),m_ExtendedSources.begin(),m_ExtendedSources.end());
	}//close if

	return 0;

}//close UpdateSourceFromWorkers()


int SourceFinderMPI::MergeSourcesAtEdge(){
	
	//To be performed by the master processor
	if (myid != MASTER_ID) return 0;

	//Loop over all tasks and merge sources detected across tasks
	for(unsigned int i=0;i<m_taskDataPerWorkers.size();i++){
		for(unsigned int j=0;j<m_taskDataPerWorkers[i].size();j++){
				
			//Loop over detected sources in this task
			for(unsigned int k=0;k<(m_taskDataPerWorkers[i][j]->sources).size();k++){
				

			}//end loop sources			
			

		}//end loop tasks per workers
	}//end loop workers	
	

	

	return 0;

}//close MergeSourcesAtEdge()


int SourceFinderMPI::FindSources(std::vector<Source*>& sources,Img* inputImg,double seedThr,double mergeThr){

	if(!inputImg) return -1;

	//## Compute stats and bkg
	INFO_LOG("Computing image stats/bkg...");
	BkgData* bkgData= ComputeStatsAndBkg(inputImg);	
	if(!bkgData){
		ERROR_LOG("Failed to compute stats/bkg info!");
		return -1;
	}

	//## Compute significance map
	Img* significanceMap= inputImg->GetSignificanceMap(bkgData,m_UseLocalBkg);
	if(!significanceMap){
		ERROR_LOG("Failed to compute significance map!");
		return -1;
	}

	//## Find sources
	INFO_LOG("Finding sources...");	
	int status= inputImg->FindCompactSource(sources,significanceMap,bkgData,seedThr,mergeThr,m_NMinPix,m_SearchNegativeExcess,m_MergeBelowSeed,m_SearchNestedSources,m_NestedBlobThrFactor);

	if(status<0) {
		ERROR_LOG("Source finding failed!");
		delete bkgData;
		bkgData= 0;
		significanceMap->Delete();
		return -1;
	}
	int nSources= (int)sources.size();
	INFO_LOG(nSources<<" sources detected in input image...");	
	
	//## Clear data
	delete bkgData;
	bkgData= 0;
	significanceMap->Delete();

	return 0;

}//close FindSources()


int SourceFinderMPI::FindCompactSources(TaskData* taskData, Img* img){

	//## Check img
	if(!img || !taskData){
		ERROR_LOG("Null ptr to input img or task data!");
		return -1;
	}
	
	//## Compute stats and bkg
	INFO_LOG("Computing image stats/bkg...");	
	BkgData* bkgData= ComputeStatsAndBkg(img);	
	if(!bkgData){
		ERROR_LOG("Failed to compute stats/bkg info!");
		return -1;
	}

	//## Compute significance map
	Img* significanceMap= img->GetSignificanceMap(bkgData,m_UseLocalBkg);
	if(!significanceMap){
		ERROR_LOG("Failed to compute significance map!");
		return -1;
	}

	//## Find sources
	INFO_LOG("Finding compact sources...");	
	int status= img->FindCompactSource(taskData->sources,significanceMap,bkgData,m_SeedThr,m_MergeThr,m_NMinPix,m_SearchNegativeExcess,m_MergeBelowSeed,m_SearchNestedSources,m_NestedBlobThrFactor);
	if(status<0) {
		ERROR_LOG("Compact source finding failed!");
		if(significanceMap) significanceMap->Delete();
		return -1;
	}

	//## Retrieve found sources 
	int nSources= (int)(taskData->sources.size());
	INFO_LOG("#"<<nSources<<" bright sources detected in input image...");
	if(nSources<=0) return 0;

	//## Apply source selection?
	int nSelSources= nSources;

	if(m_ApplySourceSelection){
		if(SelectSources(taskData->sources)<0){
		//if(SelectSources(taskData)<0){
			ERROR_LOG("Failed to select sources!");
			return -1;
		}
		nSelSources= (int)(taskData->sources.size());
	}//close if source selection
		
	INFO_LOG("#"<<nSelSources<<" compact sources added to the list...");


	//## Find if sources are at tile edge
	double xmin, xmax, ymin, ymax;
	std::vector<Source*> sources_not_at_edges;

	for(unsigned int i=0;i<(taskData->sources).size();i++){
		//Get source coordinate range
		(taskData->sources)[i]->GetSourceRange(xmin,xmax,ymin,ymax);
		
		//Check if source is inside neighbour tile
		bool isAtEdge= false;
		for(unsigned int j=0;j<(taskData->neighborWorkerId).size();j++){	
			long int neighborTaskId= (taskData->neighborTaskId)[j];
			long int neighborWorkerId= (taskData->neighborWorkerId)[j];
			
			long int ix_min= (m_taskDataPerWorkers[neighborWorkerId][neighborTaskId])->ix_min;
			long int ix_max= (m_taskDataPerWorkers[neighborWorkerId][neighborTaskId])->ix_max;
			long int iy_min= (m_taskDataPerWorkers[neighborWorkerId][neighborTaskId])->iy_min;
			long int iy_max= (m_taskDataPerWorkers[neighborWorkerId][neighborTaskId])->iy_max;

			bool isAtEdgeX= ( (xmin<=ix_max && xmin>=ix_min) || (xmax<=ix_max && xmax>=ix_min) );
			bool isAtEdgeY= ( (ymin<=iy_max && ymin>=iy_min) || (ymax<=iy_max && ymax>=iy_min) );

			if( isAtEdgeX || isAtEdgeY ){
				isAtEdge= true;
				break;
			}
		}//end loop neighbors

		//Set edge flag in source
		if(isAtEdge) {
			(taskData->sources)[i]->SetEdgeFlag(true);
			(taskData->sources_edge).push_back( (taskData->sources)[i] );
		}
		else {
			(taskData->sources)[i]->SetEdgeFlag(false);
			sources_not_at_edges.push_back( (taskData->sources)[i] );
		}
	}//end loop selected sources	

	INFO_LOG("#"<<(taskData->sources_edge).size()<<"/"<<nSelSources<<" compact sources are found at tile edges...");
	
	//Clear initial vector (DO NOT CLEAR MEMORY!) and fill with selection (then reset selection)
	(taskData->sources).clear();
	(taskData->sources).insert((taskData->sources).end(),sources_not_at_edges.begin(),sources_not_at_edges.end());
	sources_not_at_edges.clear();

	//## Add detected sources to the list	
	//m_SourceCollection.insert(m_SourceCollection.end(),sources.begin(),sources.end());
	//m_CompactSources.insert(m_CompactSources.end(),sources.begin(),sources.end());

	
	//## Clear-up
	if(bkgData){
		delete bkgData;
		bkgData= 0;
	}
	if(significanceMap) significanceMap->Delete();

	return 0;

}//close FindCompactSources()

int SourceFinderMPI::FindResidualMap(){

	//Compute residual map
	Img* residualImg= 0;
	if(m_UseResidualInExtendedSearch && m_CompactSources.size()>0){
		residualImg= m_InputImg->GetSourceResidual(m_CompactSources,m_DilateKernelSize,m_DilateSourceModel,m_DilatedSourceType,m_DilateNestedSources,m_BkgData,m_UseLocalBkg,m_DilateRandomize,m_SeedBrightThr);
	}//close if
	else{
		residualImg= m_InputImg->GetCloned("residualImg",true,true);
	}

	if(!residualImg){
		cerr<<"FindResidualMap(): ERROR: Failed to compute residual map!"<<endl;
		return -1;
	}
	m_ResidualImg= residualImg;
	
	//Compute bkg & noise map for residual img
	m_ResidualBkgData= ComputeStatsAndBkg(m_ResidualImg);
	if(!m_ResidualBkgData){
		cerr<<"FindResidualMap(): ERROR: Failed to compute bkg data for residual map!"<<endl;
		return -1;
	}


	return 0;

}//close FindResidualMap()

int SourceFinderMPI::FindExtendedSources(){

	//## Set input map for extended source search
	Img* inputImg= 0;
	Img* smoothedImg= 0;	
	if(m_UsePreSmoothing){//Apply a smoothing stage?
		if(m_SmoothFilter==eGausFilter){
			smoothedImg= m_ResidualImg->GetSmoothedImage(m_GausFilterKernSize,m_GausFilterKernSize,m_GausFilterSigma,m_GausFilterSigma);
		}
		else if(m_SmoothFilter==eGuidedFilter){
			smoothedImg= m_ResidualImg->GetGuidedFilterImage(m_GuidedFilterRadius,m_GuidedFilterColorEps);
		}
		else{
			cerr<<"SourceFinder::FindResidualMap(): ERROR: Invalid smoothing algo selected!"<<endl;
			return -1;
		}

		if(!smoothedImg){
			cerr<<"SourceFinder::FindResidualMap(): ERROR: Source residual image smoothing failed!"<<endl;
			return -1;
		}
		smoothedImg->SetNameTitle("smoothedImg","smoothedImg");
		inputImg= smoothedImg;
	}//close if use smoothing
	else{
		inputImg= m_ResidualImg;
	}

	//## Run the segmentation	
	int status= 0;
	if(m_ExtendedSearchMethod==eHClust){
		status= FindExtendedSources_HClust(inputImg);
	}
	else if(m_ExtendedSearchMethod==eChanVese){
		status= FindExtendedSources_ChanVese(inputImg);
	}
	else if(m_ExtendedSearchMethod==eWT){
		status= FindExtendedSources_WT(inputImg);
	}
	else{
		cerr<<"FindSaliency(): ERROR: Invalid extended source method selected (method="<<m_ExtendedSearchMethod<<")!"<<endl;
		return -1;
	}
	
	if(status<0){
		cerr<<"FindSaliency(): ERROR: Failed to run the segmentation algorithm!"<<endl;
		return -1;
	}

	return 0;

}//close FindExtendedSources()

int SourceFinderMPI::FindExtendedSources_HClust(Img*){

	//## Compute saliency
	m_SaliencyImg= m_ResidualImg->GetSaliencyMap(
		m_SaliencyResoMin,m_SaliencyResoMax,m_SaliencyResoStep,
		m_spBeta,m_spMinArea,m_SaliencyNNFactor,
		m_SaliencySpatialRegFactor,m_SaliencyUseRobustPars,m_SaliencyUseCurvInDiss,m_SaliencyMultiResoCombThrFactor,
		m_SaliencyUseBkgMap,m_SaliencyUseNoiseMap,m_ResidualBkgData,
		m_SaliencyThrFactor,m_SaliencyImgThrFactor
	);
	if(!m_SaliencyImg){
		cerr<<"FindSaliency(): ERROR: Failed to compute saliency map!"<<endl;
		return -1;
	}

	//## Threshold saliency map and get signal and bkg markers
	int nbins= 100;
	double signalThr= m_SaliencyImg->FindOptimalGlobalThreshold(m_SaliencyThrFactor,nbins,true);
	double bkgThr= m_SaliencyImg->FindMedianThreshold(m_SaliencyBkgThrFactor);

	double fgValue= 1;
	Img* signalMarkerImg= m_SaliencyImg->GetBinarizedImage(signalThr,fgValue,false);
	Img* bkgMarkerImg= m_SaliencyImg->GetBinarizedImage(bkgThr,fgValue,true);
	
	//## Compute the Superpixel partition
	//SLICData* slicData= SLIC::SPGenerator(this,int regionSize,double regParam, int minRegionSize, bool useLogScaleMapping, Img* edgeImg);

	
	//## Tag the superpixel partition
	//...	

	//## Run the segmentation
	//...
	//...

	//## Clear-up

	return 0;

}//close FindExtendedSources_HClust()

int SourceFinderMPI::FindExtendedSources_ChanVese(Img* inputImg){

	if(!inputImg){
		ERROR_LOG("Null ptr to input image given!");
		return -1;
	}

	//## Perform segmentation
	std::vector<Source*> sources;
	int status= inputImg->FindExtendedSource_CV(sources,m_ResidualBkgData,m_NMinPix,m_SearchNegativeExcess,m_cvTimeStepPar,m_cvWindowSizePar,m_cvLambda1Par,m_cvLambda2Par,m_cvMuPar,m_cvNuPar,m_cvPPar);
	if(status<0){
		ERROR_LOG("ChanVese Segmentation failed!");
		return -1;
	}

	//## Add sources to extended sources
	m_ExtendedSources.insert(m_ExtendedSources.end(),sources.begin(),sources.end());		
	m_SourceCollection.insert(m_SourceCollection.end(),sources.begin(),sources.end());		

	return 0;

}//close FindExtendedSources_ChanVese()

int SourceFinderMPI::FindExtendedSources_WT(Img* inputImg){

	if(!inputImg){
		//cerr<<"SourceFinder::FindExtendedSources_WT(): ERROR: Null ptr to input image given!"<<endl;
		ERROR_LOG("Null ptr to input image given!");
		return -1;
	}
	
	//## Find extended sources in the W3, W5 scales of the residual image where ONLY POINT-LIKE SOURCES are removed
	//cout<<"SourceFinder::FindExtendedSources_WT(): INFO: Find extended sources in the residual image WT-"<<m_wtScaleExtended<<"  scale ..."<<endl;
	INFO_LOG("Find extended sources in the residual image WT-"<<m_wtScaleExtended<<"  scale ...");
	std::vector<Img*> wt_extended= inputImg->GetWaveletDecomposition(m_wtScaleExtended);
	

	std::vector<Source*> sources;
	int status= FindSources(sources,wt_extended[m_wtScaleExtended],m_SeedThr,m_MergeThr);
	if(status<0){
		ERROR_LOG("Extended source finding failed!");
		return -1;
	}
	int nSources= (int)sources.size();		

	INFO_LOG("#"<<nSources<<" found...");

	//## Tag sources as extended
	for(unsigned int i=0;i<sources.size();i++){
		sources[i]->SetType(Source::eExtended);
	}

	//## Add sources to extended sources
	m_ExtendedSources.insert(m_ExtendedSources.end(),sources.begin(),sources.end());		
	m_SourceCollection.insert(m_SourceCollection.end(),sources.begin(),sources.end());		
		
	//## Clear-up
	for(unsigned int i=0;i<wt_extended.size();i++){
		if(wt_extended[i]) wt_extended[i]->Delete();
	}

	return 0;

}//close FindExtendedSources_WT()

int SourceFinderMPI::SelectSources(std::vector<Source*>& sources){
//int SourceFinderMPI::SelectSources(TaskData* taskData){

	/*
	//## Check task data
	if(!taskData){
		ERROR_LOG("Given input task data is null ptr!");
		return -1;
	}

	//## Apply source selection?
	int nSources= (int)(taskData->sources.size());
	if(nSources<=0) return 0;
	
	int nSelSources= 0;
	std::vector<Source*> sources_sel;

	for(int i=0;i<nSources;i++){	
		std::string sourceName= (taskData->sources)[i]->Name;
		int sourceId= (taskData->sources)[i]->Id;
		long int NPix= (taskData->sources)[i]->NPix;
		double X0= (taskData->sources)[i]->X0;
		double Y0= (taskData->sources)[i]->Y0;

		//Is bad source (i.e. line-like blob, etc...)?
		if(!IsGoodSource((taskData->sources)[i])) {
			DEBUG_LOG("Source no. "<<i<<" (name="<<sourceName<<",id="<<sourceId<<", n="<<NPix<<"("<<X0<<","<<Y0<<")) tagged as bad source, skipped!");
			(taskData->sources)[i]->SetGoodSourceFlag(false);
			continue;
		}
			
		//Is point-like source?
		if( IsPointLikeSource((taskData->sources)[i]) ){
			DEBUG_LOG("Source no. "<<i<<" (name="<<sourceName<<",id="<<sourceId<<", n="<<NPix<<"("<<X0<<","<<Y0<<")) tagged as a point-like source ...");
			(taskData->sources)[i]->SetType(Source::ePointLike);
		}

		//Tag nested sources
		std::vector<Source*> nestedSources= (taskData->sources)[i]->GetNestedSources();
		for(unsigned int j=0;j<nestedSources.size();j++){
			std::string nestedSourceName= nestedSources[j]->Name;
			int nestedSourceId= nestedSources[j]->Id;
			long int nestedNPix= nestedSources[j]->NPix;
			double nestedX0= nestedSources[j]->X0;
			double nestedY0= nestedSources[j]->Y0;

			if(!IsGoodSource(nestedSources[j])) {
				DEBUG_LOG("Source no. "<<i<<": nested source no. "<<j<<" (name="<<nestedSourceName<<",id="<<nestedSourceId<<", n="<<nestedNPix<<"("<<nestedX0<<","<<nestedY0<<")) tagged as bad source, skipped!");
				nestedSources[j]->SetGoodSourceFlag(false);
			}
			if( IsPointLikeSource(nestedSources[j]) ){
				DEBUG_LOG("Source no. "<<i<<": nested source no. "<<j<<" (name="<<nestedSourceName<<",id="<<nestedSourceId<<", n="<<nestedNPix<<"("<<nestedX0<<","<<nestedY0<<")) tagged as a point-like source ...");
				nestedSources[j]->SetType(Source::ePointLike);
			}
		}//end loop nested sources
					
		//Add source to the list	
		sources_sel.push_back( (taskData->sources)[i] );
		nSelSources++;
	}//end loop sources

	INFO_LOG("Added "<<nSelSources<<" bright sources to the selected list...");
	

	//## Clear initial vector (DO NOT CLEAR MEMORY!) and fill with selection (then reset selection)
	(taskData->sources).clear();

	//## Find if sources are at tile edge
	double xmin, xmax, ymin, ymax;
	for(unsigned int i=0;i<sources_sel.size();i++){	
		//Get source coordinate range
		sources_sel[i]->GetSourceRange(xmin,xmax,ymin,ymax);
		
		//Check if source is inside neighbour tile
		bool isAtEdge= false;
		for(unsigned int j=0;j<(taskData->neighborWorkerId).size();j++){	
			long int neighborTaskId= (taskData->neighborTaskId)[j];
			long int neighborWorkerId= (taskData->neighborWorkerId)[j];
			
			long int ix_min= (m_taskDataPerWorkers[neighborWorkerId][neighborTaskId])->ix_min;
			long int ix_max= (m_taskDataPerWorkers[neighborWorkerId][neighborTaskId])->ix_max;
			long int iy_min= (m_taskDataPerWorkers[neighborWorkerId][neighborTaskId])->iy_min;
			long int iy_max= (m_taskDataPerWorkers[neighborWorkerId][neighborTaskId])->iy_max;

			bool isAtEdgeX= ( (xmin<=ix_max && xmin>=ix_min) || (xmax<=ix_max && xmax>=ix_min) );
			bool isAtEdgeY= ( (ymin<=iy_max && ymin>=iy_min) || (ymax<=iy_max && ymax>=iy_min) );

			if( isAtEdgeX || isAtEdgeY ){
				isAtEdge= true;
				break;
			}
		}//end loop neighbors

		//Set edge flag in source
		if(isAtEdge) {
			sources_sel[i]->SetEdgeFlag(true);
			(taskData->sources_edge).push_back(sources_sel[i]);
		}
		else {
			sources_sel[i]->SetEdgeFlag(false);
			(taskData->sources).push_back(sources_sel[i]);
		}
	}//end loop selected sources
	
	//Reset source selection (DO NOT CLEAR MEMORY)
	sources_sel.clear();
	*/

	
	//## Apply source selection?
	int nSources= (int)sources.size();
	if(nSources<=0) return 0;
	
	int nSelSources= 0;
	std::vector<Source*> sources_sel;

	for(int i=0;i<nSources;i++){	
		std::string sourceName= sources[i]->Name;
		int sourceId= sources[i]->Id;
		long int NPix= sources[i]->NPix;
		double X0= sources[i]->X0;
		double Y0= sources[i]->Y0;

		//Is bad source (i.e. line-like blob, etc...)?
		if(!IsGoodSource(sources[i])) {
			DEBUG_LOG("Source no. "<<i<<" (name="<<sourceName<<",id="<<sourceId<<", n="<<NPix<<"("<<X0<<","<<Y0<<")) tagged as bad source, skipped!");
			sources[i]->SetGoodSourceFlag(false);
			continue;
		}
			
		//Is point-like source?
		if( IsPointLikeSource(sources[i]) ){
			DEBUG_LOG("Source no. "<<i<<" (name="<<sourceName<<",id="<<sourceId<<", n="<<NPix<<"("<<X0<<","<<Y0<<")) tagged as a point-like source ...");
			sources[i]->SetType(Source::ePointLike);
		}

		//Tag nested sources
		std::vector<Source*> nestedSources= sources[i]->GetNestedSources();
		for(unsigned int j=0;j<nestedSources.size();j++){
			std::string nestedSourceName= nestedSources[j]->Name;
			int nestedSourceId= nestedSources[j]->Id;
			long int nestedNPix= nestedSources[j]->NPix;
			double nestedX0= nestedSources[j]->X0;
			double nestedY0= nestedSources[j]->Y0;

			if(!IsGoodSource(nestedSources[j])) {
				DEBUG_LOG("Source no. "<<i<<": nested source no. "<<j<<" (name="<<nestedSourceName<<",id="<<nestedSourceId<<", n="<<nestedNPix<<"("<<nestedX0<<","<<nestedY0<<")) tagged as bad source, skipped!");
				nestedSources[j]->SetGoodSourceFlag(false);
			}
			if( IsPointLikeSource(nestedSources[j]) ){
				DEBUG_LOG("Source no. "<<i<<": nested source no. "<<j<<" (name="<<nestedSourceName<<",id="<<nestedSourceId<<", n="<<nestedNPix<<"("<<nestedX0<<","<<nestedY0<<")) tagged as a point-like source ...");
				nestedSources[j]->SetType(Source::ePointLike);
			}
		}//end loop nested sources
			
		//Add source to the list	
		sources_sel.push_back(sources[i]);
		nSelSources++;
	}//end loop sources

	INFO_LOG("Added "<<nSelSources<<" bright sources to the selected list...");
	
	//Clear initial vector (DO NOT CLEAR MEMORY!) and fill with selection (then reset selection)
	sources.clear();
	sources.insert(sources.end(),sources_sel.begin(),sources_sel.end());
	sources_sel.clear();

	return 0;

}//close SelectSources()

bool SourceFinderMPI::IsGoodSource(Source* aSource){
	
	if(!aSource) return false;

	//## Check for pixels 	
	if(aSource->NPix<=0 || (aSource->GetPixels()).size()<=0) return false;

	//## Check for line-like source
	if( (aSource->GetContours()).size()<=0) {
		WARN_LOG("No contour stored for this source, cannot perform check!");
		return true;
	}

	double BoundingBoxMin= ((aSource->GetContours())[0])->BoundingBoxMin;
	if(BoundingBoxMin<m_SourceMinBoundingBox) {
		DEBUG_LOG("BoundingBox cut not passed (BoundingBoxMin="<<BoundingBoxMin<<"<"<<m_SourceMinBoundingBox<<")");
		return false;
	}

	//## Add other check here ...
	//...
	//...

	return true;

}//close IsGoodSource()

bool SourceFinderMPI::IsPointLikeSource(Source* aSource){

	if(!aSource) return false;
	if(!aSource->HasParameters()) {
		WARN_LOG("No parameters are available for this source (did you compute them?)...point-like check cannot be performed!");
		return true;
	}

	std::string sourceName= aSource->Name;
	int sourceId= aSource->Id;

	//Loop over contours and check if all of them have circular features
	bool isPointLike= true;
	std::vector<Contour*> contours= aSource->GetContours();

	for(unsigned int i=0;i<contours.size();i++){
		Contour* thisContour= contours[i];

		/*
		//Test circularity ratio: 1= circle
		if(thisContour->CircularityRatio<m_psCircRatioThr) {
			DEBUG_LOG("Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass CircularityRatio cut (CR="<<thisContour->CircularityRatio<<"<"<<psCircRatioThr<<")");
			isPointLike= false;
			break;
		}
		*/

		//Test elongation (how symmetrical is the shape): 0=circle,square
		if(thisContour->Elongation>m_psElongThr) {
			DEBUG_LOG("Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass Elongation cut (ELONG="<<thisContour->CircularityRatio<<">"<<m_psElongThr<<")");
			isPointLike= false;
			break;	
		}

		//Test ellipse fit
		if(thisContour->EllipseAreaRatio<m_psEllipseAreaRatioMinThr || thisContour->EllipseAreaRatio>m_psEllipseAreaRatioMaxThr) {
			DEBUG_LOG("Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass EllipseAreaRatio cut (EAR="<<thisContour->EllipseAreaRatio<<" outside range ["<<m_psEllipseAreaRatioMinThr<<","<<m_psEllipseAreaRatioMaxThr<<"])");
			isPointLike= false;
			break;	
		}

	}//end contour loop
	
	//Check number of pixels
	DEBUG_LOG("Source (name="<<sourceName<<","<<"id="<<sourceId<<") (NPix="<<aSource->NPix<<">"<<m_psMaxNPix<<")");
	if(aSource->NPix>m_psMaxNPix){
		DEBUG_LOG("Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass nMaxPix cut (NPix="<<aSource->NPix<<">"<<m_psMaxNPix<<")");
		isPointLike= false;
	}

	if(!isPointLike) return false;

	return true;

}//close IsPointLikeSource()



Img* SourceFinderMPI::ReadImage(long int ix_min,long int ix_max,long int iy_min,long int iy_max){

	//## Check file
	FileInfo info;
	bool match_extension= false;
	if(!SysUtils::CheckFile(m_InputFileName,info,match_extension,"")){
		ERROR_LOG("Invalid input file name specified (filename="<<m_InputFileName<<"), invalid file path?!");
		return 0;
	}
	m_InputFileExtension= info.extension;


	Img* img= 0;

	//=== ROOT reading ===
	if(m_InputFileExtension==".root"){// Read image from ROOT file
		TFile* inputFile = new TFile(m_InputFileName.c_str(),"READ");
		if(!inputFile || inputFile->IsZombie()){
			ERROR_LOG("Cannot open input file "<<m_InputFileName<<"!");
			return 0;
		}
		
		Img* fullImg= (Img*)inputFile->Get(m_InputImgName.c_str());
		if(!fullImg){
			ERROR_LOG("Cannot get image from input file "<<m_InputFileName<<"!");
			return 0;
		}
		img= fullImg->GetTile(ix_min,ix_max,iy_min,iy_max);	
		if(!img){
			ERROR_LOG("Null ptr to read image!");
			fullImg->Delete();
			return 0;
		}
		fullImg->Delete();
	}//close if

	//=== FITS reading ===
	else if(m_InputFileExtension==".fits"){// Read image from FITS file
		img= new Img;
		int status= img->ReadFITS(m_InputFileName,ix_min,ix_max,iy_min,iy_max);
		if(status<0){
			ERROR_LOG("Failed to read image from input file "<<m_InputFileName<<"!");
			if(img) img->Delete();
			return 0;
		}
	}//close else if

	//== Invalid extension ==
	else{
		ERROR_LOG("Invalid file extension detected (ext="<<m_InputFileExtension<<")!");
		return 0;
	}
	img->SetNameTitle("img","img");
	
	return img;

}//close ReadImage()


BkgData* SourceFinderMPI::ComputeStatsAndBkg(Img* img){

	//## Check input img
	if(!img){
		ERROR_LOG("Null ptr to input image given!");
		return 0;
	}

	//## Compute stats
	INFO_LOG("Computing image stats...");
	bool computeRobustStats= true;
	bool skipNegativePix= false;
	bool forceRecomputing= false;
	if(!img->ComputeStats(computeRobustStats,skipNegativePix,forceRecomputing)<0){
		ERROR_LOG("Stats computing failed!");
		return 0;
	}
	//img->PrintStats();
	img->LogStats("INFO");

	//## Set local bkg grid/box
	//## If MetaData & beam info are available, interpret grid&box options as multiple of beam
	//## If no info is available (or use of beam info is off) interpret grid&box options as fractions wrt image size
	double boxSizeX= m_BoxSizeX;
	double boxSizeY= m_BoxSizeY;
	int nPixelsInBeam= 0;
	if(m_UseBeamInfoInBkg && img->HasMetaData()){
		nPixelsInBeam= img->GetMetaData()->GetBeamSizeInPixel();	
	}
	
	if(m_UseBeamInfoInBkg && nPixelsInBeam>0){
		INFO_LOG("Setting bkg boxes as ("<<m_BoxSizeX<<","<<m_BoxSizeY<<") x beam (beam="<<nPixelsInBeam<<" pixels) ...");
		boxSizeX= nPixelsInBeam*m_BoxSizeX;
		boxSizeY= nPixelsInBeam*m_BoxSizeY;
	}
	else{
		WARN_LOG("Beam information is not available or its usage has been turned off, using image fractions...");
		double Nx= img->GetNbinsX();
		double Ny= img->GetNbinsY();
		boxSizeX= m_BoxSizeX*Nx;
		boxSizeY= m_BoxSizeY*Ny;
	}

	double gridSizeX= m_GridSizeX*boxSizeX;
	double gridSizeY= m_GridSizeY*boxSizeY;

	//## Compute Bkg
	BkgData* bkgData= img->ComputeBkg(m_BkgEstimator,m_UseLocalBkg,boxSizeX,boxSizeY,gridSizeX,gridSizeY,m_Use2ndPassInLocalBkg,m_SkipOutliersInLocalBkg,m_SeedThr,m_MergeThr,m_NMinPix);

	if(!bkgData) {
		ERROR_LOG("Bkg computing failed!");
		return 0;
	}
		
	return bkgData;

}//close ComputeStatsAndBkg()


int SourceFinderMPI::DrawSources(Img* image,std::vector<Source*>& sources){

	//cout<<"SourceFinder::DrawSources(): INFO: Drawing sources..."<<endl;
	INFO_LOG("Drawing sources...");
	if(!image) return -1;

	bool useCurrentCanvas= false;
	bool drawFull= false;
	image->Plot(sources,useCurrentCanvas,drawFull,eRAINBOW,true);
	
	return 0;

}//close DrawSources()


int SourceFinderMPI::Save(){

	INFO_LOG("Storing results to file & catalog...");

	//Save DS9 regions?
	if(m_SaveDS9Region && m_DS9CatalogFilePtr){
		DEBUG_LOG("Saving DS9 region header...");
	
		fprintf(m_DS9CatalogFilePtr,"global color=red font=\"helvetica 12 normal\" edit=1 move=1 delete=1 include=1\n");
		fprintf(m_DS9CatalogFilePtr,"image\n");

		DEBUG_LOG("Saving "<<m_SourceCollection.size()<<" sources to file...");

		for(unsigned int k=0;k<m_SourceCollection.size();k++){
			DEBUG_LOG("Dumping DS9 region info for source no. "<<k<<" ...");
			std::string regionInfo= "";
			if(m_DS9RegionFormat==1) regionInfo= m_SourceCollection[k]->GetDS9Region(true);
			else if(m_DS9RegionFormat==2) regionInfo= m_SourceCollection[k]->GetDS9EllipseRegion(true);
			else continue;

			fprintf(m_DS9CatalogFilePtr,"%s\n",regionInfo.c_str());
	  
			DEBUG_LOG("Set source ptr to source "<<k<<"...");
			m_Source= m_SourceCollection[k];
			DEBUG_LOG("Saving source no. "<<k<<" to tree...");
			m_SourceTree->Fill();
		}//end loop sources
		
		DEBUG_LOG("Closing DS9 file region...");
		fclose(m_DS9CatalogFilePtr);
		m_DS9CatalogFilePtr= 0;
	}//close if SaveDS9Region()


	//Check ROOT output file
	if(!m_OutputFile) {
		WARN_LOG("Null ptr to output file, nothing will be saved in ROOT file!");
		return -1;
	}
	m_OutputFile->cd();

	//Save source tree?
	if(m_SaveSources){
		DEBUG_LOG("Writing tree to file...");
		m_SourceTree->Write();
	}
	
	//Save config?
	if(m_SaveConfig){
		TTree* ConfigTree= ConfigParser::Instance().GetConfigTree("ConfigInfo");
		if(ConfigTree) ConfigTree->Write();
	}	
	
	//Save residual map
	if(m_SaveResidualMap && m_ResidualImg){
		m_ResidualImg->Write();
	}

	/*
	//Save image to file?
	if(m_SaveImageToFile){
		if(fSaveImageType==eInputImage) fInputImg->Write();
		else if(fSaveImageType==eResidualImage){
			Img* residualImg= fInputImg->GetSourceResidual(fUseLocalBackground,fSourceDilateKernelSize,fDilateNestedSources,fDilatedSourceType,fDilateSourceModel,fRandomizeInDilate,fRandSigmaInDilate);
			if(residualImg) residualImg->Write();
		}
		else if(fSaveImageType==eSegmentedImage && fSearchExtendedSources){
			if(fSegmentedImage) fSegmentedImage->Write();
			if(fFinalSegmentedImage) fFinalSegmentedImage->Write();
			if(fFinalSignalSegmentedImage) fFinalSignalSegmentedImage->Write();
			if(fSaliencyMap) fSaliencyMap->Write();
			if(fSumSaliencyMap) fSumSaliencyMap->Write();
			if(fInitSaliencyMap) fInitSaliencyMap->Write();
			if(fSegmentationContourGraph) fSegmentationContourGraph->Write();
		}
		else{
			cerr<<"SourceFinder::Save(): WARN: Invalid save image flag specified!"<<endl;
		}

		if(fSPMergingInfo) fSPMergingInfo->Write();
		if(fConfigInfo) fConfigInfo->Write();
	}//close if save image to file
	*/

	//cout<<"SourceFinder::Save(): INFO: Closing output file..."<<endl;	
	DEBUG_LOG("Closing output file...");
	m_OutputFile->Close();

	//cout<<"SourceFinder::Save(): INFO: End save to file"<<endl;
	INFO_LOG("End save to file");

	return 0;

}//close Save()

}//close namespace
