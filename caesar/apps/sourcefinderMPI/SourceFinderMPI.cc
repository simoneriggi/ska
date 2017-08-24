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
#include <Image.h>
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
#include <Consts.h>

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

	//Free comm & groups	
	MPI_Group_free(&m_WorkerGroup);
	MPI_Comm_free(&m_WorkerComm);
	MPI_Group_free(&m_WorldGroup);

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

	m_TaskInfoTree= 0;
		
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

			//Init task info tree
			if(!m_TaskInfoTree) m_TaskInfoTree= new TTree("TaskInfo","TaskInfo");
			m_TaskInfoTree->Branch("xmin",&m_xmin);
			m_TaskInfoTree->Branch("xmax",&m_xmax);
			m_TaskInfoTree->Branch("ymin",&m_ymin);
			m_TaskInfoTree->Branch("ymax",&m_ymax);
	

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
	else {
		m_TileMinX= 0;
		m_TileMaxX= 0;
		m_TileMinY= 0;
		m_TileMaxY= 0;
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
	GET_OPTION_VALUE(saliencyDissExpFalloffPar,m_SaliencyDissExpFalloffPar);
	GET_OPTION_VALUE(saliencySpatialDistRegPar,m_SaliencySpatialDistRegPar);

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
	DEBUG_LOG("[PROC "<<myid<<"] - Generated jobId: "<<jobId);
	
	//## Get input image size
	FileInfo info;
	bool match_extension= false;
	if(!SysUtils::CheckFile(m_InputFileName,info,match_extension,"")){
		ERROR_LOG("[PROC "<<myid<<"] - Invalid input file name specified (filename="<<m_InputFileName<<"), invalid file path?!");
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
				ERROR_LOG("[PROC "<<myid<<"] - Failed to open input file image "<<m_InputFileName<<" and get image size!");
				return -1;
			}
			Image* inputImg= (Image*)inputFile->Get(m_InputImgName.c_str());
			if(!inputImg) {
				ERROR_LOG("[PROC "<<myid<<"] - Failed to open input file image "<<m_InputFileName<<" and get image size!");
				return -1;	
			}
			Nx= inputImg->GetNx();
			Ny= inputImg->GetNy();
		}
	}//close if
	else if(info.extension==".fits"){//FITS
		if(m_ReadTile){
			Nx= m_TileMaxX-m_TileMinX+1;
			Ny= m_TileMaxY-m_TileMinY+1;
		}
		else{//READ FULL MAP
			if(SysUtils::GetFITSImageSize(m_InputFileName,Nx,Ny)<0){
				ERROR_LOG("[PROC "<<myid<<"] - Failed to open input file image "<<m_InputFileName<<" and get image size!");
				return -1;
			}
			INFO_LOG("[PROC "<<myid<<"] - Input image opened with success: size="<<Nx<<"x"<<Ny);
		}
	}//close else if		
	else{
		ERROR_LOG("[PROC "<<myid<<"] - Invalid/unsupported file extension ("<<info.extension<<") detected!");
		return -1;
	}

	INFO_LOG("[PROC "<<myid<<"] - Image size: "<<Nx<<"x"<<Ny);

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

	INFO_LOG("[PROC "<<myid<<"] - Computing tile partition: tileOverlap("<<tileOverlapX<<","<<tileOverlapY<<")");
	if(MathUtils::Compute2DGrid(ix_min,ix_max,iy_min,iy_max,Nx,Ny,m_TileSizeX,m_TileSizeY,m_TileStepSizeX,m_TileStepSizeY)<0){
		WARN_LOG("[PROC "<<myid<<"] - Failed to compute a 2D partition from input image!");
		return -1;
	}
	int nExpectedTasks= ix_min.size()*iy_min.size();
	INFO_LOG("[PROC "<<myid<<"] - #"<<nExpectedTasks<<" expected number of tasks ("<<ix_min.size()<<"x"<<iy_min.size()<<")");

	//## Compute worker tasks (check max number of tasks per worker)
	INFO_LOG("[PROC "<<myid<<"] - Computing worker task list...");
	TaskData* aTaskData= 0;
	long int workerCounter= 0;

	for(unsigned int j=0;j<iy_min.size();j++){
		for(unsigned int i=0;i<ix_min.size();i++){
			//Assign worker
			INFO_LOG("[PROC "<<myid<<"] - Assign task ("<<i<<","<<j<<") to worker no. "<<workerCounter<<"...");
				
			aTaskData= new TaskData;
			aTaskData->filename= m_InputFileName;
			aTaskData->jobId= jobId;
			aTaskData->workerId= workerCounter;
			//aTaskData->taskId= workerCounter;
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
	std::vector<int> workerIds;

	for(unsigned int i=0;i<m_taskDataPerWorkers.size();i++){
		if(m_taskDataPerWorkers[i].size()==0) continue;//no tasks present
		workerIds.push_back(i);

		//Loop over tasks present in this worker
		int nTasksInWorker= (int)(m_taskDataPerWorkers[i].size()); 
		//std::stringstream ss;	
		//ss<<"[PROC "<<myid<<"] - Worker no. "<<i<<", ";

		for(int j=0;j<nTasksInWorker;j++){
			long int ix_min= m_taskDataPerWorkers[i][j]->ix_min;
			long int ix_max= m_taskDataPerWorkers[i][j]->ix_max;
			long int iy_min= m_taskDataPerWorkers[i][j]->iy_min;
			long int iy_max= m_taskDataPerWorkers[i][j]->iy_max;
			
			//Find first neighbors among tasks inside the same worker
			for(int k=j+1;k<nTasksInWorker;k++){
				if(j==k) continue;
				long int next_ix_min= m_taskDataPerWorkers[i][k]->ix_min;
				long int next_ix_max= m_taskDataPerWorkers[i][k]->ix_max;
				long int next_iy_min= m_taskDataPerWorkers[i][k]->iy_min;
				long int next_iy_max= m_taskDataPerWorkers[i][k]->iy_max;
				
				bool isAdjacentInX= (ix_max==next_ix_min-1 || ix_min==next_ix_max+1 || (ix_min==next_ix_min && ix_max==next_ix_max));
				bool isAdjacentInY= (iy_max==next_iy_min-1 || iy_min==next_iy_max+1 || (iy_min==next_iy_min && iy_max==next_iy_max));
				bool isOverlappingInX= ( (next_ix_min>=ix_min && next_ix_min<=ix_max) || (next_ix_max>=ix_min && next_ix_max<=ix_max) );
				bool isOverlappingInY= ( (next_iy_min>=iy_min && next_iy_min<=iy_max) || (next_iy_max>=iy_min && next_iy_max<=iy_max) );
				bool isAdjacent= isAdjacentInX && isAdjacentInY;
				bool isOverlapping= isOverlappingInX && isOverlappingInY;
 
				/*
				bool isAdjacentInX= false;
				bool isAdjacentInY= false;
				if(m_UseTileOverlap){
					
				}
				else{
					isAdjacentInX= (ix_max==next_ix_min-1 || ix_min==next_ix_max+1 || (ix_min==next_ix_min && ix_max==next_ix_max) );
					isAdjacentInY= (iy_max==next_iy_min-1 || iy_min==next_iy_max+1 || (iy_min==next_iy_min && iy_max==next_iy_max) );

					isNeighborInX= (ix_max==next_ix_min-1 || ix_min==next_ix_max+1) && (iy_min==next_iy_min && iy_max==next_iy_max);
					isNeighborInY= (iy_max==next_iy_min-1 || iy_min==next_iy_max+1) && (ix_min==next_ix_min && ix_max==next_ix_max);
				}
				*/

				std::stringstream ss;	
				ss<<"[PROC "<<myid<<"] - Worker no. "<<i<<", Task "<<j<<"["<<ix_min<<","<<ix_max<<"] ["<<iy_min<<","<<iy_max<<"], NextTask "<<k<<"["<<next_ix_min<<","<<next_ix_max<<"] ["<<next_iy_min<<","<<next_iy_max<<"] ==> isAdjacentInX? "<<isAdjacentInX<<", isAdjacentInY? "<<isAdjacentInY;
				if(myid==0) INFO_LOG(ss.str());


				if(isAdjacent || isOverlapping) {
					(m_taskDataPerWorkers[i][j]->neighborTaskId).push_back(k);
					(m_taskDataPerWorkers[i][k]->neighborTaskId).push_back(j);
					(m_taskDataPerWorkers[i][j]->neighborWorkerId).push_back(i);
					(m_taskDataPerWorkers[i][k]->neighborWorkerId).push_back(i);
				}

			}//end loop next task in worker


			//Find neighbors across workers
			for(unsigned int s=i+1;s<m_taskDataPerWorkers.size();s++){
				for(unsigned int t=0;t<m_taskDataPerWorkers[s].size();t++){
					
					long int next_ix_min= m_taskDataPerWorkers[s][t]->ix_min;
					long int next_ix_max= m_taskDataPerWorkers[s][t]->ix_max;
					long int next_iy_min= m_taskDataPerWorkers[s][t]->iy_min;
					long int next_iy_max= m_taskDataPerWorkers[s][t]->iy_max;
				
					bool isAdjacentInX= (ix_max==next_ix_min-1 || ix_min==next_ix_max+1 || (ix_min==next_ix_min && ix_max==next_ix_max));
					bool isAdjacentInY= (iy_max==next_iy_min-1 || iy_min==next_iy_max+1 || (iy_min==next_iy_min && iy_max==next_iy_max));
					bool isOverlappingInX= ( (next_ix_min>=ix_min && next_ix_min<=ix_max) || (next_ix_max>=ix_min && next_ix_max<=ix_max) );
					bool isOverlappingInY= ( (next_iy_min>=iy_min && next_iy_min<=iy_max) || (next_iy_max>=iy_min && next_iy_max<=iy_max) );
					bool isAdjacent= isAdjacentInX && isAdjacentInY;
					bool isOverlapping= isOverlappingInX && isOverlappingInY;
 

					/*
					if(m_UseTileOverlap){
						isAdjacentInX= ( (next_ix_min>=ix_min && next_ix_min<=ix_max) || (next_ix_max>=ix_min && next_ix_max<=ix_max) );
						isAdjacentInY= ( (next_iy_min>=iy_min && next_iy_min<=iy_max) || (next_iy_max>=iy_min && next_iy_max<=iy_max) );
					}
					else{
						isAdjacentInX= (ix_max==next_ix_min-1 || ix_min==next_ix_max+1 || (ix_min==next_ix_min && ix_max==next_ix_max));
						isAdjacentInY= (iy_max==next_iy_min-1 || iy_min==next_iy_max+1 || (iy_min==next_iy_min && iy_max==next_iy_max));
					}
					*/

					std::stringstream ss;	
					ss<<"[PROC "<<myid<<"] - Worker no. "<<i<<", Task "<<j<<", NextWorker no. "<<s<<", NextTask "<<t<<"["<<next_ix_min<<","<<next_ix_max<<"] ["<<next_iy_min<<","<<next_iy_max<<"] ==> isAdjacentInX? "<<isAdjacentInX<<", isAdjacentInY? "<<isAdjacentInY;
					if(myid==0) INFO_LOG(ss.str());

					if(isAdjacent || isOverlapping) {
						(m_taskDataPerWorkers[i][j]->neighborTaskId).push_back(t);
						(m_taskDataPerWorkers[s][t]->neighborTaskId).push_back(j);
						(m_taskDataPerWorkers[i][j]->neighborWorkerId).push_back(s);
						(m_taskDataPerWorkers[s][t]->neighborWorkerId).push_back(i);
					}
				
				}//end loop tasks in next worker
			}//end loop workers 

		}//end loop tasks
	}//end loop workers

	int nWorkers= (int)workerIds.size();
	std::stringstream ss;
	ss<<"[PROC "<<myid<<"] - # "<<nWorkers<<" workers {";
	for(int i=0;i<nWorkers;i++){
		ss<<workerIds[i]<<",";
	}
	ss<<"}";
	INFO_LOG(ss.str());
	
	//Get main processor group
	MPI_Comm_group(MPI_COMM_WORLD, &m_WorldGroup);
	
	// Construct a group containing all of the workers (proc with tasks assigned)
	MPI_Group_incl(m_WorldGroup, nWorkers, workerIds.data() , &m_WorkerGroup);

	// Create a new communicator based on the group
	int commTag= 10;
	MPI_Comm_create_group(MPI_COMM_WORLD, m_WorkerGroup, commTag, &m_WorkerComm);

	worker_ranks = -1, nworkers = -1;
	// If this rank isn't in the new communicator, it will be
	// MPI_COMM_NULL. Using MPI_COMM_NULL for MPI_Comm_rank or
	// MPI_Comm_size is erroneous
	if (m_WorkerComm!=MPI_COMM_NULL) {
    MPI_Comm_rank(m_WorkerComm, &worker_ranks);
    MPI_Comm_size(m_WorkerComm, &nworkers);
	}
	else {
		WARN_LOG("[PROC "<<myid<<"] - Worker MPI communicator is null (this processor has no tasks and was not inserted in the worker group)!");
	}
	INFO_LOG("[PROC "<<myid<<"] - WORLD RANK/SIZE: "<<myid<<"/"<<nproc<<" WORKER RANK/SIZE: "<<worker_ranks<<"/"<<nworkers);


	//Print
	if(myid==0){
	for(unsigned int i=0;i<m_taskDataPerWorkers.size();i++){
		if(m_taskDataPerWorkers[i].size()==0) continue;//no tasks present

		
		//INFO_LOG("== Worker no. "<<i<<"==");
		for(unsigned int j=0;j<m_taskDataPerWorkers[i].size();j++){
			std::stringstream ss;	
			ss<<"[PROC "<<myid<<"] - Worker no. "<<i<<", ";

			long int ix_min= m_taskDataPerWorkers[i][j]->ix_min;
			long int ix_max= m_taskDataPerWorkers[i][j]->ix_max;
			long int iy_min= m_taskDataPerWorkers[i][j]->iy_min;
			long int iy_max= m_taskDataPerWorkers[i][j]->iy_max;
			ss<<"Task no. "<<j<<"["<<ix_min<<","<<ix_max<<"] ["<<iy_min<<","<<iy_max<<"] neighbors{";
			//INFO_LOG("Task no. "<<j<<"["<<ix_min<<","<<ix_max<<"] ["<<iy_min<<","<<iy_max<<"] neighbors{");	
			for(unsigned int k=0;k<m_taskDataPerWorkers[i][j]->neighborTaskId.size();k++){
				long int neighborWorkerId= m_taskDataPerWorkers[i][j]->neighborWorkerId[k];
				long int neighborTaskId= m_taskDataPerWorkers[i][j]->neighborTaskId[k];
				long int next_ix_min= m_taskDataPerWorkers[neighborWorkerId][neighborTaskId]->ix_min;
				long int next_ix_max= m_taskDataPerWorkers[neighborWorkerId][neighborTaskId]->ix_max;
				long int next_iy_min= m_taskDataPerWorkers[neighborWorkerId][neighborTaskId]->iy_min;
				long int next_iy_max= m_taskDataPerWorkers[neighborWorkerId][neighborTaskId]->iy_max;

				//INFO_LOG("("<<neighborWorkerId<<","<<neighborTaskId<<") ["<<next_ix_min<<","<<next_ix_max<<"] ["<<next_iy_min<<","<<next_iy_max<<"]");
				ss<<"("<<neighborWorkerId<<","<<neighborTaskId<<") ["<<next_ix_min<<","<<next_ix_max<<"] ["<<next_iy_min<<","<<next_iy_max<<"], ";
			}	
			ss<<"}";
			INFO_LOG(ss.str());	
		}//end loop tasks
		
	}//end loop workers
	}//close if


	bool hasTooManyTasks= false;
	int maxNTasksPerWorker_default= 10000;
	for(unsigned int i=0;i<m_taskDataPerWorkers.size();i++){
		long int nTasksPerWorker= (long int)m_taskDataPerWorkers[i].size();
		if(nTasksPerWorker>maxNTasksPerWorker_default){
			hasTooManyTasks= true;
			break;
		}
	}

	if(hasTooManyTasks){
		WARN_LOG("[PROC "<<myid<<"] - Too many tasks per worker (thr="<<maxNTasksPerWorker_default<<")");
		return -1;
	}

	return 0;

}//close PrepareWorkerTasks()


int SourceFinderMPI::Run(){

	double globalTimerStart= MPI_Wtime();

	//## Init options (done by all processors)
	double initTimerStart= MPI_Wtime();
	INFO_LOG("[PROC "<<myid<<"] - Initializing source finder...");
	if(Init()<0){
		ERROR_LOG("[PROC "<<myid<<"] - Initialization failed!");
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
		INFO_LOG("[PROC "<<myid<<"] - Reading input image ["<<ix_min<<","<<ix_max<<"] ["<<iy_min<<","<<iy_max<<"]...");
		Image* taskImg= ReadImage(ix_min,ix_max,iy_min,iy_max);
		if(!taskImg){
			ERROR_LOG("Reading of input image failed, skip to next task...");
			continue;
		}

		//## Set image physical boundary in task data
		long int Nx= taskImg->GetNx();
		long int Ny= taskImg->GetNy();
		double xmin= taskImg->GetXmin();//taskImg->GetXaxis()->GetBinCenter(1);
		double xmax= taskImg->GetXmax();//taskImg->GetXaxis()->GetBinCenter(Nx);
		double ymin= taskImg->GetYmin();//taskImg->GetYaxis()->GetBinCenter(1);
		double ymax= taskImg->GetYmax();//taskImg->GetYaxis()->GetBinCenter(Ny);
		m_taskDataPerWorkers[myid][j]->x_min= xmin;
		m_taskDataPerWorkers[myid][j]->x_max= xmax;
		m_taskDataPerWorkers[myid][j]->y_min= ymin;
		m_taskDataPerWorkers[myid][j]->y_max= ymax;
		
		//## Find compact sources
		INFO_LOG("[PROC "<<myid<<"] - Searching compact sources...");
		if(m_SearchCompactSources && FindCompactSources(m_taskDataPerWorkers[myid][j],taskImg)<0){
			ERROR_LOG("[PROC "<<myid<<"] - Compact source search failed!");
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
	
	//## Update task info (tile physical range) from workers
	//## The updated list of task data is available in master processor
	//## (check if it is better to replace with MPI_Gather and put it available in all workers)
	if(UpdateTaskDataFromWorkers()<0){
		ERROR_LOG("[PROC "<<myid<<"] - Updating task data from workers failed!");
		return -1;
	}

	//## Find edge sources
	MPI_Barrier(MPI_COMM_WORLD);
	if(myid==0) {
		if(FindSourcesAtEdge()<0){
			ERROR_LOG("[PROC "<<myid<<"] - Finding sources at tile edges failed!");
			return -1;
		}
	}
	
	
	//## Deblend sources
	//if(fDeblendSources) DeblendSources(fInputImg);
	
	//## Draw & Store results (done by master processor)
	if(myid==0) {
		//## Draw final sources
		if(m_IsInteractiveRun) DrawSources(m_InputImg,m_SourceCollection);
	
		//## Save to file
		if(m_SaveToFile) Save();	
		//if(m_Application && m_IsInteractiveRun) m_Application->Run();
	}//close if

	return 0;

}//close Run()


int SourceFinderMPI::FindSourcesAtEdge(){

	
	//## Find if sources (both compact and extended) are at tile edge
	//## Those found at the edge are removed from the list and added to the edge list for further processing
	double xmin_s, xmax_s, ymin_s, ymax_s;
	
	//Loop over workers
	for(unsigned int i=0;i<m_taskDataPerWorkers.size();i++){
		if(m_taskDataPerWorkers[i].size()==0) continue;//no tasks present

		//Loop over tasks per worker
		for(unsigned int j=0;j<m_taskDataPerWorkers[i].size();j++){
	
			//Loop over compact sources found
			std::vector<Source*> sources_not_at_edges;
			int nSources= (int)((m_taskDataPerWorkers[i][j]->sources).size());

			for(int k=0;k<nSources;k++){
				//Get source coordinate range
				(m_taskDataPerWorkers[i][j]->sources)[k]->GetSourceRange(xmin_s,xmax_s,ymin_s,ymax_s);

				//Check if source is at the edge of its tile
				long int xmin_tile= (m_taskDataPerWorkers[i][j])->ix_min;
				long int xmax_tile= (m_taskDataPerWorkers[i][j])->ix_max;
				long int ymin_tile= (m_taskDataPerWorkers[i][j])->iy_min;
				long int ymax_tile= (m_taskDataPerWorkers[i][j])->iy_max; 
				bool isAtTileEdgeX= (xmin_s==xmin_tile || xmax_s==xmax_tile);
				bool isAtTileEdgeY= (ymin_s==ymin_tile || ymax_s==ymax_tile);
				bool isAtTileEdge= (isAtTileEdgeX || isAtTileEdgeY);
				INFO_LOG("[PROC "<<myid<<"] - workerId="<<m_taskDataPerWorkers[i][j]->workerId<<", check if compact source no. "<<k<<"(x["<<xmin_s<<","<<xmax_s<<"] y["<<ymin_s<<","<<ymax_s<<"]) is at edge of its tile (x["<<xmin_tile<<","<<xmax_tile<<"] y["<<ymin_tile<<","<<ymax_tile<<"]), isAtTileEdgeX="<<isAtTileEdgeX<<", isAtTileEdgeY="<<isAtTileEdgeY<<", isAtTileEdge="<<isAtTileEdge);

				INFO_LOG("[PROC "<<myid<<"] - workerId="<<m_taskDataPerWorkers[i][j]->workerId<<", check if compact source no. "<<k<<"(x["<<xmin_s<<","<<xmax_s<<"] y["<<ymin_s<<","<<ymax_s<<"]) is inside neighbour tile...");
		
				//Check if source is inside neighbour tile, e.g. is in overlapping area
				//bool isAtEdge= false;
				bool isInOverlapArea= false;
				for(unsigned int l=0;l<(m_taskDataPerWorkers[i][j]->neighborWorkerId).size();l++){	
					long int neighborTaskId= (m_taskDataPerWorkers[i][j]->neighborTaskId)[l];
					long int neighborWorkerId= (m_taskDataPerWorkers[i][j]->neighborWorkerId)[l];
			
					long int xmin= (m_taskDataPerWorkers[neighborWorkerId][neighborTaskId])->ix_min;
					long int xmax= (m_taskDataPerWorkers[neighborWorkerId][neighborTaskId])->ix_max;
					long int ymin= (m_taskDataPerWorkers[neighborWorkerId][neighborTaskId])->iy_min;
					long int ymax= (m_taskDataPerWorkers[neighborWorkerId][neighborTaskId])->iy_max;

					//bool isAtEdgeX= ( (xmin_s<=xmax && xmin_s>=xmin) || (xmax_s<=xmax && xmax_s>=xmin) );
					//bool isAtEdgeY= ( (ymin_s<=ymax && ymin_s>=ymin) || (ymax_s<=ymax && ymax_s>=ymin) );
					//INFO_LOG("[PROC "<<myid<<"] - neighborWorkerId="<<neighborWorkerId<<", neighborTaskId="<<neighborTaskId<<", check if inside neighbor tile no. "<<j<<"(x["<<xmin<<","<<xmax<<"] y["<<ymin<<","<<ymax<<"]), isAtEdgeX="<<isAtEdgeX<<", isAtEdgeY="<<isAtEdgeY);
					//if( isAtEdgeX && isAtEdgeY ){
					//	isAtEdge= true;
					//	break;
					//}
					bool isOverlappingX= ( (xmin_s<=xmax && xmin_s>=xmin) || (xmax_s<=xmax && xmax_s>=xmin) );
					bool isOverlappingY= ( (ymin_s<=ymax && ymin_s>=ymin) || (ymax_s<=ymax && ymax_s>=ymin) );

					INFO_LOG("[PROC "<<myid<<"] - neighborWorkerId="<<neighborWorkerId<<", neighborTaskId="<<neighborTaskId<<", check if inside neighbor tile no. "<<j<<"(x["<<xmin<<","<<xmax<<"] y["<<ymin<<","<<ymax<<"]), isOverlappingX="<<isOverlappingX<<", isOverlappingY="<<isOverlappingY);
					if( isOverlappingX && isOverlappingY ){
						isInOverlapArea= true;
						break;
					}
				}//end loop neighbors

				//Tag source at edge is is located at the border of its tile or if it is located inside an overlapping area with another neighbor tile
				bool isAtEdge= (isAtTileEdge || isInOverlapArea);
			

				//Set edge flag in source
				if(isAtEdge) {
					(m_taskDataPerWorkers[i][j]->sources)[k]->SetEdgeFlag(true);
					(m_taskDataPerWorkers[i][j]->sources_edge).push_back( (m_taskDataPerWorkers[i][j]->sources)[k] );
				}
				else {
					(m_taskDataPerWorkers[i][j]->sources)[k]->SetEdgeFlag(false);
					sources_not_at_edges.push_back( (m_taskDataPerWorkers[i][j]->sources)[k] );
				}
			}//end loop sources	

			int nEdgeSources= (int)((m_taskDataPerWorkers[i][j]->sources_edge).size());
			INFO_LOG("[PROC "<<myid<<"] - #"<<nEdgeSources<<"/"<<nSources<<" compact sources are found at tile edges...");
	
			//Clear initial vector (DO NOT CLEAR MEMORY!) and fill with selection (then reset selection)
			(m_taskDataPerWorkers[i][j]->sources).clear();
			(m_taskDataPerWorkers[i][j]->sources).insert( 
				(m_taskDataPerWorkers[i][j]->sources).end(),
				sources_not_at_edges.begin(),
				sources_not_at_edges.end()
			);
			sources_not_at_edges.clear();

		}//end loop tasks per worker
	}//end loop workers

	return 0;

}//close FindSourcesAtEdge()



int SourceFinderMPI::UpdateTaskDataFromWorkers(){

	//## Put a barrier and collect all sources from workers in the master processor
	MPI_Barrier(MPI_COMM_WORLD);

	int MSG_TAG= 1;
	if (myid == 0) {//Receive data from the workers

		for (int i=1; i<nproc; i++) {
			//Check if this processor has tasks assigned, otherwise skip!
			if(m_taskDataPerWorkers[i].size()==0){
				INFO_LOG("[PROC "<<myid<<"] - No tasks assigned to worker no. "<<i<<", nothing to be collected, skip to next worker...");
				continue;
			}
  
			//## Probe for an incoming message from process zero
			INFO_LOG("[PROC "<<myid<<"] - Probing for message from proc "<<i);
    	MPI_Status status;

			if(MPI_Probe(i, MSG_TAG, MPI_COMM_WORLD, &status)==MPI_SUCCESS){
				INFO_LOG("[PROC "<<myid<<"] - a message has been found with the probe, with tag " << status.MPI_TAG << ", source " << status.MPI_SOURCE);

    		//## When probe returns, the status object has the size and other
    		//## attributes of the incoming message. Get the message size
				INFO_LOG("[PROC "<<myid<<"] - Getting size of message... ");
    	
				int rcvMsgSize= 0;
    		MPI_Get_count(&status, MPI_CHAR, &rcvMsgSize);

				//## Allocate a buffer to hold the incoming numbers
				INFO_LOG("[PROC "<<myid<<"] - Allocating a message of size "<<rcvMsgSize);
				if(rcvMsgSize<=0){
					ERROR_LOG("[PROC "<<myid<<"] - rcvMsg size is negative/null!");
					continue;
				}
				char* recvBuffer= (char*)malloc(rcvMsgSize);
    		
    		//## Now receive the message with the allocated buffer
    		//MPI_Recv(recvBuffer, rcvMsgSize, MPI_CHAR, MPI_ANY_SOURCE, MSG_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(recvBuffer, rcvMsgSize, MPI_CHAR, i, MSG_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    		INFO_LOG("[PROC "<<myid<<"] - Received a message of size "<<rcvMsgSize<<") from process "<<i);

				//## Update task data with received worker data	
				bool isTaskCollectionPreAllocated= true;
				if(Serializer::CharArrayToTaskDataCollection(m_taskDataPerWorkers[i],recvBuffer,rcvMsgSize,isTaskCollectionPreAllocated)<0 ){
					ERROR_LOG("[PROC "<<myid<<"] - Failed to decode recv message into task data list!");
    			if(recvBuffer) free(recvBuffer);
				}

				//## Free received buffer
    		if(recvBuffer) free(recvBuffer);
			}//close if
			else{
				ERROR_LOG("[PROC "<<myid<<"] - Message probing failed!");
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
			ERROR_LOG("[PROC "<<myid<<"] - Failed to encode task data to protobuf!");
			return -1;
		}

		//## Send buffer to master processor	
		INFO_LOG("[PROC "<<myid<<"] - Sending msg: "<<msg<<" (size="<<msg_size<<")");
		MPI_Send((void*)(msg),msg_size, MPI_CHAR, MASTER_ID, MSG_TAG, MPI_COMM_WORLD);

		//## Free buffer
		free(msg);
	}//close else

	MPI_Barrier(MPI_COMM_WORLD);


	//## Update sources in list
	if (myid == 0) {

		//Print task data
		INFO_LOG("[PROC "<<myid<<"] - Printing task data...");
		for(unsigned int i=0;i<m_taskDataPerWorkers.size();i++){
			if(m_taskDataPerWorkers[i].size()==0) continue;//no tasks present
			
			std::stringstream ss;
			ss<<"[PROC "<<myid<<"] - Worker no. "<<i<<", ";
			for(unsigned int j=0;j<m_taskDataPerWorkers[i].size();j++){
				long int ix_min= m_taskDataPerWorkers[i][j]->ix_min;
				long int ix_max= m_taskDataPerWorkers[i][j]->ix_max;
				long int iy_min= m_taskDataPerWorkers[i][j]->iy_min;
				long int iy_max= m_taskDataPerWorkers[i][j]->iy_max;
				double x_min= m_taskDataPerWorkers[i][j]->x_min;
				double x_max= m_taskDataPerWorkers[i][j]->x_max;
				double y_min= m_taskDataPerWorkers[i][j]->y_min;
				double y_max= m_taskDataPerWorkers[i][j]->y_max;

				m_xmin= x_min;
				m_xmax= x_max;
				m_ymin= y_min;
				m_ymax= y_max;
				if(m_TaskInfoTree) m_TaskInfoTree->Fill();

				ss<<"Task no. "<<j<<", PixelRange["<<ix_min<<","<<ix_max<<"] ["<<iy_min<<","<<iy_max<<"], PhysCoordRange["<<x_min<<","<<x_max<<"] ["<<y_min<<","<<y_max<<"], ";
			}//end loop tasks
			INFO_LOG(ss.str());			
	
		}//end loop workers

		for(unsigned int i=0;i<m_taskDataPerWorkers.size();i++){
			for(unsigned int j=0;j<m_taskDataPerWorkers[i].size();j++){
				m_CompactSources.insert(m_CompactSources.end(),(m_taskDataPerWorkers[i][j]->sources).begin(),(m_taskDataPerWorkers[i][j]->sources).end());
				m_CompactSources.insert(m_CompactSources.end(),(m_taskDataPerWorkers[i][j]->sources_edge).begin(),(m_taskDataPerWorkers[i][j]->sources_edge).end());

				m_ExtendedSources.insert(m_ExtendedSources.end(),(m_taskDataPerWorkers[i][j]->ext_sources).begin(),(m_taskDataPerWorkers[i][j]->ext_sources).end());
				m_ExtendedSources.insert(m_ExtendedSources.end(),(m_taskDataPerWorkers[i][j]->ext_sources_edge).begin(),(m_taskDataPerWorkers[i][j]->ext_sources_edge).end());
			}//end loop tasks
		}//end loop workers
		m_SourceCollection.insert(m_SourceCollection.end(),m_CompactSources.begin(),m_CompactSources.end());
		m_SourceCollection.insert(m_SourceCollection.end(),m_ExtendedSources.begin(),m_ExtendedSources.end());

		INFO_LOG("[PROC "<<myid<<"] - #"<<m_SourceCollection.size()<<" sources found in total (#"<<m_CompactSources.size()<<" compact, #"<<m_ExtendedSources.size()<<" extended) ...");
		
	}//close if

	return 0;

}//close UpdateTaskDataFromWorkers()

/*
int SourceFinderMPI::UpdateTaskDataFromWorkers(){
	
	//## Returns if no tasks are present for this proc
	if(m_taskDataPerWorkers[myid].size()==0){
		INFO_LOG("Proc #"<<myid<<": Nothing to be broadcasted (no tasks assigned to this worker), simply return!");
		return 0;
	}

	MPI_Barrier(m_WorkerComm);

	//## If this proc has task data, encode them for sending
	long int msg_size= 0;
	char* msg= Serializer::TaskDataCollectionToCharArray(msg_size,m_taskDataPerWorkers[myid]);
	if(!msg){
		ERROR_LOG("Proc #"<<myid<<": Failed to encode task data to protobuf!");
		return -1;
	}
	INFO_LOG("Proc #"<<myid<<": Broadcasting a buffer of size "<<msg_size<<" to other processors...");
	

	//## Gather task data sizes in advance
	long long int* msg_sizes = (long long int*)malloc(nworkers);
	MPI_Allgather(&msg_size, 1, MPI_LONG_LONG_INT, msg_sizes, 1, MPI_LONG_LONG_INT, m_WorkerComm);
		
	int maxSize= -1;
	for (int i=0;i<nworkers;i++) {	
		if(msg_sizes[i]>maxSize) maxSize= msg_sizes[i];
		INFO_LOG("Proc #"<<myid<<": msg_sizes["<<i<<"]="<<msg_sizes[i]);
	}

	char* sndBuffer= (char*)malloc(maxSize);
	for(int i=0;i<maxSize;i++) sndBuffer[i]= '#';
	memcpy ( sndBuffer, msg, msg_size );
	
	//## Allocate recv buffer and gather task data from all workers
	MPI_Barrier(m_WorkerComm);

	char** recvBuffer = (char**)malloc(nworkers);
	for(int i=0;i<nworkers;i++){
		recvBuffer[i] = (char*)malloc(maxSize);
  	//recvBuffer[i] = (char*)malloc(msg_sizes[i]);
	}

	INFO_LOG("Proc #"<<myid<<": Gather task data from all workers...");
	//MPI_Allgather( msg, msg_size, MPI_CHAR, recvBuffer, maxSize, MPI_CHAR, m_WorkerComm);
	MPI_Allgather( sndBuffer, maxSize, MPI_CHAR, recvBuffer, maxSize, MPI_CHAR, m_WorkerComm);

	
	//## Update task data with received worker data		
	INFO_LOG("Proc #"<<myid<<": Update task data with received worker data...");
	bool isTaskCollectionPreAllocated= true;
	bool isDecodingFailed= false;
	for(int i=0;i<nworkers;i++){
		if(Serializer::CharArrayToTaskDataCollection(m_taskDataPerWorkers[i],recvBuffer[i],msg_sizes[i],isTaskCollectionPreAllocated)<0 ){
			ERROR_LOG("Proc #"<<myid<<": Failed to decode recv message into task data list!");
    	isDecodingFailed= true;
		}
	}
	

	//## Free buffers
	INFO_LOG("Proc #"<<myid<<": Freeing buffers...");
	if(msg_sizes) free(msg_sizes);
	for(int i=0;i<nworkers;i++){
  	if(recvBuffer[i]) free(recvBuffer[i]);
	}
	free(recvBuffer);
	
	
	if(isDecodingFailed){
		return -1;
	}
	

	//## All processors call broadcast except those that have not assigned no tasks
	//MPI_Request requests[nproc];

	//for (int i=0;i<nproc;i++) {	
		//Skip myself...
		//if(myid==i) continue;

		//Check if this processor has tasks assigned, otherwise skip!
		//if(m_taskDataPerWorkers[i].size()==0){
		//	INFO_LOG("No tasks assigned to worker no. "<<i<<", nothing to be collected, skip to next worker...");
		//	continue;
		//}

		//## Read data broadcasted from the i-th process
		//int bcast_source= i;

		
		//Read buffer size
		//INFO_LOG("Proc #"<<myid<<": Getting size of message from process "<<bcast_source<<"... ");
		//int rcvMsgSize= 0;
		//MPI_Bcast(&rcvMsgSize, 1, MPI_INT, bcast_source, MPI_COMM_WORLD);
		

		// Allocate a buffer to hold the incoming data
		//long int rcvMsgSize= msg_sizes[i];
		//INFO_LOG("Proc #"<<myid<<": Allocating a message of size "<<rcvMsgSize<<" from process "<<bcast_source<<"... ");
		//if(rcvMsgSize<=0){
		//	ERROR_LOG("Proc #"<<myid<<": rcvMsg size is negative!");
		//	continue;
		//}
		//char* recvBuffer= (char*)malloc(rcvMsgSize);

		//Receive data buffer		
		//MPI_Ibcast(recvBuffer, rcvMsgSize, MPI_CHAR, bcast_source, MPI_COMM_WORLD,&request[i]);
		//INFO_LOG("Proc #"<<myid<<": Received a message (msg="<<recvBuffer<<", size="<<rcvMsgSize<<") from process "<<bcast_source<<"... ");

		//## Update task data with received worker data	
		//bool isTaskCollectionPreAllocated= true;
		//if(Serializer::CharArrayToTaskDataCollection(m_taskDataPerWorkers[bcast_source],recvBuffer,rcvMsgSize,isTaskCollectionPreAllocated)<0 ){
		//	ERROR_LOG("Proc #"<<myid<<": Failed to decode recv message from process "<<bcast_source<<" into task data list!");
    //	if(recvBuffer) free(recvBuffer);
		//}

		//## Free received buffer
    //if(recvBuffer) free(recvBuffer);
	//}//end loop processes
	
	//MPI_Status status[nproc];
	//int ierr= MPI_Waitall(nproc, requests, status); 
  //if(ierr!=MPI_SUCCESS){
	//	ERROR_LOG("MPI_Waitall failed rank "<<myid);
	//	return -1;
	//}
	


	
	//## All processors call broadcast except those that have not assigned no tasks
	//if(m_taskDataPerWorkers[myid].size()==0){
	//	INFO_LOG("Proc #"<<myid<<": Nothing to be broadcasted (no tasks assigned to this worker), simply return!");
	//	return 0;
	//}

	//## Put a barrier and gather all task data from workers in the master processor
	//MPI_Barrier(MPI_COMM_WORLD);

	//## Encode taskData in protobuf and send first data size to other processes
	//## MPI_Probe is not going to work with bcast
	//long int msg_size= 0;
	//char* msg= Serializer::TaskDataCollectionToCharArray(msg_size,m_taskDataPerWorkers[myid]);
	//if(!msg){
	//	ERROR_LOG("Proc #"<<myid<<": Failed to encode task data to protobuf!");
	//	return -1;
	//}

	//INFO_LOG("Proc #"<<myid<<": Broadcasting msg: "<<msg<<" (size="<<msg_size<<") to other processors...");
	
  //MPI_Bcast(&msg_size, 1, MPI_INT, myid, MPI_COMM_WORLD);
	//MPI_Bcast(msg, msg_size, MPI_CHAR, myid, MPI_COMM_WORLD);
	
	//for (int i=0;i<nproc;i++) {	
		//Skip myself...
	//	if(myid==i) continue;

		//Check if this processor has tasks assigned, otherwise skip!
	//	if(m_taskDataPerWorkers[i].size()==0){
	//		INFO_LOG("No tasks assigned to worker no. "<<i<<", nothing to be collected, skip to next worker...");
	//		continue;
	//	}

		//## Read data broadcasted from the i-th process
	//	int bcast_source= i;

		//Read buffer size
	//	INFO_LOG("Proc #"<<myid<<": Getting size of message from process "<<bcast_source<<"... ");
	//	int rcvMsgSize= 0;
	//	MPI_Bcast(&rcvMsgSize, 1, MPI_INT, bcast_source, MPI_COMM_WORLD);
		
		// Allocate a buffer to hold the incoming data
	//	INFO_LOG("Proc #"<<myid<<": Allocating a message of size "<<rcvMsgSize<<" from process "<<bcast_source<<"... ");
	//	if(rcvMsgSize<=0){
	//		ERROR_LOG("Proc #"<<myid<<": rcvMsg size is negative!");
	//		continue;
	//	}
	//	char* recvBuffer= (char*)malloc(rcvMsgSize);

		//Receive data buffer		
	//	MPI_Bcast(recvBuffer, rcvMsgSize, MPI_CHAR, bcast_source, MPI_COMM_WORLD);
	//	INFO_LOG("Proc #"<<myid<<": Received a message (msg="<<recvBuffer<<", size="<<rcvMsgSize<<") from process "<<bcast_source<<"... ");

		//## Update task data with received worker data	
	//	bool isTaskCollectionPreAllocated= true;
	//	if(Serializer::CharArrayToTaskDataCollection(m_taskDataPerWorkers[bcast_source],recvBuffer,rcvMsgSize,isTaskCollectionPreAllocated)<0 ){
	//		ERROR_LOG("Proc #"<<myid<<": Failed to decode recv message from process "<<bcast_source<<" into task data list!");
  //  	if(recvBuffer) free(recvBuffer);
	//	}

		//## Free received buffer
  //  if(recvBuffer) free(recvBuffer);

	//}//end loop processes
	

	//## Synchronize again before getting the results
  MPI_Barrier(MPI_COMM_WORLD);

	if(msg) free(msg);

	//Print task data
	INFO_LOG("Proc #"<<myid<<": Printing task data...");
	for(unsigned int i=0;i<m_taskDataPerWorkers.size();i++){
		if(m_taskDataPerWorkers[i].size()==0) continue;//no tasks present

		INFO_LOG("== Worker no. "<<i<<"==");
		for(unsigned int j=0;j<m_taskDataPerWorkers[i].size()-1;j++){
			long int ix_min= m_taskDataPerWorkers[i][j]->ix_min;
			long int ix_max= m_taskDataPerWorkers[i][j]->ix_max;
			long int iy_min= m_taskDataPerWorkers[i][j]->iy_min;
			long int iy_max= m_taskDataPerWorkers[i][j]->iy_max;
			double x_min= m_taskDataPerWorkers[i][j]->x_min;
			double x_max= m_taskDataPerWorkers[i][j]->x_max;
			double y_min= m_taskDataPerWorkers[i][j]->y_min;
			double y_max= m_taskDataPerWorkers[i][j]->y_max;
			INFO_LOG("Task no. "<<j<<"PixelRange["<<ix_min<<","<<ix_max<<"] ["<<iy_min<<","<<iy_max<<"] PhysCoordRange["<<x_min<<","<<x_max<<"] ["<<y_min<<","<<y_max<<"]");

		}//end loop tasks
	}//end loop workers
	
	return 0;

}//close UpdateTaskDataFromWorkers()
*/

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


int SourceFinderMPI::FindSources(std::vector<Source*>& sources,Image* inputImg,double seedThr,double mergeThr){

	if(!inputImg) return -1;

	//## Compute stats and bkg
	INFO_LOG("Computing image stats/bkg...");
	ImgBkgData* bkgData= ComputeStatsAndBkg(inputImg);	
	if(!bkgData){
		ERROR_LOG("Failed to compute stats/bkg info!");
		return -1;
	}

	//## Compute significance map
	Image* significanceMap= inputImg->GetSignificanceMap(bkgData,m_UseLocalBkg);
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


int SourceFinderMPI::FindCompactSources(TaskData* taskData, Image* img){

	//## Check img
	if(!img || !taskData){
		ERROR_LOG("Null ptr to input img or task data!");
		return -1;
	}
	
	//## Compute stats and bkg
	INFO_LOG("Computing image stats/bkg...");	
	ImgBkgData* bkgData= ComputeStatsAndBkg(img);	
	if(!bkgData){
		ERROR_LOG("Failed to compute stats/bkg info!");
		return -1;
	}

	//## Compute significance map
	Image* significanceMap= img->GetSignificanceMap(bkgData,m_UseLocalBkg);
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
	Image* residualImg= 0;
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
	Image* inputImg= 0;
	Image* smoothedImg= 0;	
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
	else if(m_ExtendedSearchMethod==eWaveletTransform){
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

int SourceFinderMPI::FindExtendedSources_HClust(Image*){

	//## Compute saliency
	m_SaliencyImg= m_ResidualImg->GetMultiResoSaliencyMap(
		m_SaliencyResoMin,m_SaliencyResoMax,m_SaliencyResoStep,
		m_spBeta,m_spMinArea,m_SaliencyNNFactor,m_SaliencyUseRobustPars,m_SaliencyDissExpFalloffPar,m_SaliencySpatialDistRegPar,
		m_SaliencyMultiResoCombThrFactor,
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
	Image* signalMarkerImg= m_SaliencyImg->GetBinarizedImage(signalThr,fgValue,false);
	Image* bkgMarkerImg= m_SaliencyImg->GetBinarizedImage(bkgThr,fgValue,true);
	
	//## Compute the Superpixel partition
	//SLICData* slicData= SLIC::SPGenerator(this,int regionSize,double regParam, int minRegionSize, bool useLogScaleMapping, Image* edgeImg);

	
	//## Tag the superpixel partition
	//...	

	//## Run the segmentation
	//...
	//...

	//## Clear-up

	return 0;

}//close FindExtendedSources_HClust()

int SourceFinderMPI::FindExtendedSources_ChanVese(Image* inputImg){

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

int SourceFinderMPI::FindExtendedSources_WT(Image* inputImg){

	if(!inputImg){
		//cerr<<"SourceFinder::FindExtendedSources_WT(): ERROR: Null ptr to input image given!"<<endl;
		ERROR_LOG("Null ptr to input image given!");
		return -1;
	}
	
	//## Find extended sources in the W3, W5 scales of the residual image where ONLY POINT-LIKE SOURCES are removed
	//cout<<"SourceFinder::FindExtendedSources_WT(): INFO: Find extended sources in the residual image WT-"<<m_wtScaleExtended<<"  scale ..."<<endl;
	INFO_LOG("Find extended sources in the residual image WT-"<<m_wtScaleExtended<<"  scale ...");
	std::vector<Image*> wt_extended= inputImg->GetWaveletDecomposition(m_wtScaleExtended);
	

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



Image* SourceFinderMPI::ReadImage(long int ix_min,long int ix_max,long int iy_min,long int iy_max){

	//## Check file
	FileInfo info;
	bool match_extension= false;
	if(!SysUtils::CheckFile(m_InputFileName,info,match_extension,"")){
		ERROR_LOG("Invalid input file name specified (filename="<<m_InputFileName<<"), invalid file path?!");
		return 0;
	}
	m_InputFileExtension= info.extension;


	Image* img= 0;

	//=== ROOT reading ===
	if(m_InputFileExtension==".root"){// Read image from ROOT file
		TFile* inputFile = new TFile(m_InputFileName.c_str(),"READ");
		if(!inputFile || inputFile->IsZombie()){
			ERROR_LOG("Cannot open input file "<<m_InputFileName<<"!");
			return 0;
		}
		
		Image* fullImg= (Image*)inputFile->Get(m_InputImgName.c_str());
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
		img= new Image;
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


ImgBkgData* SourceFinderMPI::ComputeStatsAndBkg(Image* img){

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
	if(img->ComputeStats(computeRobustStats,skipNegativePix,forceRecomputing)<0){
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
		double Nx= static_cast<double>(img->GetNx());
		double Ny= static_cast<double>(img->GetNy());
		boxSizeX= m_BoxSizeX*Nx;
		boxSizeY= m_BoxSizeY*Ny;
	}

	double gridSizeX= m_GridSizeX*boxSizeX;
	double gridSizeY= m_GridSizeY*boxSizeY;

	//## Compute Bkg
	ImgBkgData* bkgData= img->ComputeBkg(m_BkgEstimator,m_UseLocalBkg,boxSizeX,boxSizeY,gridSizeX,gridSizeY,m_Use2ndPassInLocalBkg,m_SkipOutliersInLocalBkg,m_SeedThr,m_MergeThr,m_NMinPix);

	if(!bkgData) {
		ERROR_LOG("Bkg computing failed!");
		return 0;
	}
		
	return bkgData;

}//close ComputeStatsAndBkg()


int SourceFinderMPI::DrawSources(Image* image,std::vector<Source*>& sources){

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

	//Save task info to file?
	if(m_TaskInfoTree){
		m_TaskInfoTree->Write();
	}

	//cout<<"SourceFinder::Save(): INFO: Closing output file..."<<endl;	
	DEBUG_LOG("Closing output file...");
	m_OutputFile->Close();

	//cout<<"SourceFinder::Save(): INFO: End save to file"<<endl;
	INFO_LOG("End save to file");

	return 0;

}//close Save()

}//close namespace
