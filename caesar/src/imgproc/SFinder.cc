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
* @file SourceFinder.cc
* @class SourceFinder
* @brief Source finder class
*
* Class to perform source finding 
* @author S. Riggi
* @date 20/01/2015
*/

#include <SFinder.h>
#include <BlobFinder.h>
#include <Image.h>
#include <Source.h>
#include <Contour.h>
#include <ConfigParser.h>
#include <BkgData.h>
#include <CodeUtils.h>
#include <MathUtils.h>
#include <SysUtils.h>
#include <Logger.h>
#include <Consts.h>

#include <SLIC.h>
#include <SLICSegmenter.h>
#include <ChanVeseSegmenter.h>
#include <LRACSegmenter.h>

#include <TaskData.h>
#include <Serializer.h>
#include <Graph.h>

#include <TObject.h>
#include <TMatrixD.h>
#include <TFile.h>
#include <TCanvas.h>

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
#include <chrono>

#ifdef MPI_ENABLED
#include <mpi.h>
#endif

using namespace std;

ClassImp(Caesar::SFinder)

namespace Caesar {

#define MASTER_ID 0
#define MAX_NTASKS_PER_WORKER 10000
 
SFinder::SFinder() 
{
	
}//close costructor


SFinder::~SFinder(){
	
	//Clearup
	INFO_LOG("Clearup source finder allocated data...");
	Clear();

}//close destructor


void SFinder::Clear()
{

	//## Close open ROOT file
	//## NB: When objects are written to file their memory is released, so don't delete them!
	if(m_procId==MASTER_ID){
		if( m_OutputFile && m_OutputFile->IsOpen() ) {
			DEBUG_LOG("Closing output ROOT file...");
			m_OutputFile->Close();
		}
	}

	//## Delete source tree
	if(m_procId==MASTER_ID){
		if(m_SourceTree && !m_saveSources) {
			DEBUG_LOG("Deleting source tree...");
			m_SourceTree->Delete();
		}
	}

	//## Delete images & objects not written to file
	if(m_InputImg && !m_saveInputMap){	
		DEBUG_LOG("Deleting input image...");
		delete m_InputImg;
		m_InputImg= 0;
	}
	if(m_BkgData && !m_saveBkgMap && !m_saveNoiseMap){
		DEBUG_LOG("Deleting bkg data...");
		delete m_BkgData;
		m_BkgData= 0;
	}
	if(m_ResidualImg && !m_saveResidualMap) {
		DEBUG_LOG("Deleting residual image...");
		delete m_ResidualImg;
		m_ResidualImg= 0;
	}
	if(m_SignificanceMap && !m_saveSignificanceMap) {
		delete m_SignificanceMap;
		m_SignificanceMap= 0;
	}			
	if(m_EdgeImg && !m_saveEdgenessMap){
		DEBUG_LOG("Deleting edgeness image...");
		delete m_EdgeImg;
		m_EdgeImg= 0;
	}
	if(m_LaplImg && !m_saveCurvatureMap){
		DEBUG_LOG("Deleting Laplacian image...");
		delete m_LaplImg;
		m_LaplImg= 0;
	}
	if(m_SegmImg && !m_saveSegmentedMap){
		DEBUG_LOG("Deleting Segmented image...");
		delete m_SegmImg;
		m_SegmImg= 0;
	}
	if(m_SaliencyImg && !m_saveSaliencyMap){
		DEBUG_LOG("Deleting Saliency image...");
		delete m_SaliencyImg;
		m_SaliencyImg= 0;
	}
	

	//## Delete TApplication??	
	//if(m_Application) {
	//	m_Application->Delete();
	//}

	//## Delete task data
	for(size_t i=0;i<m_taskDataPerWorkers.size();i++) {
		for(size_t j=0;j<m_taskDataPerWorkers[i].size();j++) {
			if(m_taskDataPerWorkers[i][j]){
				delete m_taskDataPerWorkers[i][j];
				m_taskDataPerWorkers[i][j]= 0;
			}
		}
		m_taskDataPerWorkers[i].clear();
	}
	m_taskDataPerWorkers.clear();


	//## Free comm & groups	
	#ifdef MPI_ENABLED
		if(m_mpiEnabled){
			MPI_Group_free(&m_WorkerGroup);
			MPI_Comm_free(&m_WorkerComm);
			MPI_Group_free(&m_WorldGroup);
		}
	#endif


}//close Clear()


void SFinder::InitOptions()
{

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
	m_saveToFile= true;
	m_saveConfig= true;
	m_saveSources= true;	
	m_saveResidualMap= true;
	m_saveInputMap= true;
	m_saveSaliencyMap= false;
	m_saveEdgenessMap= false;
	m_saveCurvatureMap= false;
	m_saveSegmentedMap= true;
	m_SourceTree= 0;
	m_saveDS9Region= true;
	m_DS9CatalogFileName= "";
	m_DS9RegionFormat= 1;
	m_PerfTree= 0;
		
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

	//Performance stats
	totTime= 0;
	initTime= 0;
	readImageTime= 0;
	compactSourceTime= 0;
	sourceSelectionTime= 0;
	imgResidualTime= 0;
	extendedSourceTime= 0;
	sourceDeblendTime= 0;	
	saveTime= 0;

	//Extended source finder
	m_EdgeImg= 0;
	m_LaplImg= 0;
	m_SegmImg= 0;

	//MPI vars (initialized after to actual value)
	//NB: When MPI is not used this should define only 1 process and 1 master
	m_nProc= 1;
	m_procId= 0;

	//Task info
	m_TaskInfoTree= 0;


	//Check is MPI run is enabled at build & runtime
	m_mpiEnabled= SysUtils::IsMPIInitialized();

}//close InitOptions()

int SFinder::Init(){

	//## Init options
	InitOptions();

	//## Configure from parser
	if(Configure()<0){
		ERROR_LOG("Failed to configure options from parser!");
		return -1;
	}

	//## Create TApplication if interactive run is selected
	if(!m_Application && m_IsInteractiveRun){
		m_Application= new TApplication("Application", 0, 0);
	}	

	//## Create output file
	//## NB: Done only by processor 0 in MPI run
	if(m_saveToFile && m_procId==MASTER_ID){
		m_OutputFile= new TFile(m_OutputFileName.c_str(),"RECREATE");	
		m_OutputFile->cd();
	
		//Init source tree
		if(m_saveSources){
			m_Source= 0;
			if(!m_SourceTree) m_SourceTree= new TTree("SourceInfo","SourceInfo");
			m_SourceTree->Branch("Source",&m_Source);
			m_SourceCollection.clear();
		}

		//Init task info tree
		if(!m_TaskInfoTree) m_TaskInfoTree= new TTree("TaskInfo","TaskInfo");
		m_TaskInfoTree->Branch("xmin",&m_xmin);
		m_TaskInfoTree->Branch("xmax",&m_xmax);
		m_TaskInfoTree->Branch("ymin",&m_ymin);
		m_TaskInfoTree->Branch("ymax",&m_ymax);
	

		//Init time performance tree
		if(!m_PerfTree) m_PerfTree= new TTree("PerformanceInfo","PerformanceInfo");
		m_PerfTree->Branch("tot",&totTime,"tot/D");
		m_PerfTree->Branch("init",&initTime,"init/D");
		m_PerfTree->Branch("read",&readImageTime,"read/D");
		m_PerfTree->Branch("sfinder",&compactSourceTime,"sfinder/D");
		m_PerfTree->Branch("sselector",&sourceSelectionTime,"sselector/D");
		m_PerfTree->Branch("imgres",&imgResidualTime,"imgres/D");
		m_PerfTree->Branch("extsfinder",&extendedSourceTime,"extsfinder/D");
		m_PerfTree->Branch("sdeblend",&sourceDeblendTime,"sdeblend/D");
		//m_PerfTree->Branch("save",&saveTime,"save/D");

	}//close if saveToFile


	//## Init and fill task data
	//## NB: Done by all processors in MPI run
	for(int i=0;i<m_nProc;i++){
		m_taskDataPerWorkers.push_back( std::vector<TaskData*>() );
	}
	if(PrepareWorkerTasks()<0){
		ERROR_LOG("Preparation of tasks per worker failed!");
		return -1;
	}

	return 0;

}//close Init()

int SFinder::Configure(){

	//Get image read options
	if(GET_OPTION_VALUE(inputFile,m_InputFileName)<0){
		ERROR_LOG("Failed to get inputFile option!");
		return -1;
	}
	GET_OPTION_VALUE(inputImage,m_InputImgName);
	GET_OPTION_VALUE(fitsHDUId,m_fitsHDUId);
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
	GET_OPTION_VALUE(saveToFile,m_saveToFile);
	GET_OPTION_VALUE(saveConfig,m_saveConfig);
	GET_OPTION_VALUE(saveDS9Region,m_saveDS9Region);
	GET_OPTION_VALUE(ds9RegionFile,m_DS9CatalogFileName);
	GET_OPTION_VALUE(DS9RegionFormat,m_DS9RegionFormat);
	GET_OPTION_VALUE(saveSources,m_saveSources);
	GET_OPTION_VALUE(isInteractiveRun,m_IsInteractiveRun);
	GET_OPTION_VALUE(saveResidualMap,m_saveResidualMap);
	GET_OPTION_VALUE(saveInputMap,m_saveInputMap);
	GET_OPTION_VALUE(saveSignificanceMap,m_saveSignificanceMap);
	GET_OPTION_VALUE(saveBkgMap,m_saveBkgMap);
	GET_OPTION_VALUE(saveNoiseMap,m_saveNoiseMap);
	GET_OPTION_VALUE(saveSaliencyMap,m_saveSaliencyMap);
	GET_OPTION_VALUE(saveEdgenessMap,m_saveEdgenessMap);
	GET_OPTION_VALUE(saveCurvatureMap,m_saveCurvatureMap);
	GET_OPTION_VALUE(saveSegmentedMap,m_saveSegmentedMap);
		
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
	
	//Get source deblending options
	GET_OPTION_VALUE(deblendSources,m_deblendSources);
	GET_OPTION_VALUE(deblendCurvThr,m_deblendCurvThr);
	GET_OPTION_VALUE(deblendComponentMinNPix,m_deblendComponentMinNPix);
	

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
	GET_OPTION_VALUE(useResidualInExtendedSearch,m_UseResidualInExtendedSearch);
	GET_OPTION_VALUE(activeContourMethod,m_activeContourMethod);
			
	//Get superpixel options
	GET_OPTION_VALUE(spSize,m_spSize);
	GET_OPTION_VALUE(spBeta,m_spBeta);
	GET_OPTION_VALUE(spMinArea,m_spMinArea);
	GET_OPTION_VALUE(spUseLogContrast,m_spUseLogContrast);

	//Chan-Vese options
	GET_OPTION_VALUE(cvNIters,m_cvNIters);
	GET_OPTION_VALUE(cvTimeStepPar,m_cvTimeStepPar);
	GET_OPTION_VALUE(cvWindowSizePar,m_cvWindowSizePar);
	GET_OPTION_VALUE(cvLambda1Par,m_cvLambda1Par);
	GET_OPTION_VALUE(cvLambda2Par,m_cvLambda2Par);
	GET_OPTION_VALUE(cvMuPar,m_cvMuPar);
	GET_OPTION_VALUE(cvNuPar,m_cvNuPar);
	GET_OPTION_VALUE(cvPPar,m_cvPPar);

	//LRAC algorithm options
	GET_OPTION_VALUE(lracNIters,m_lracNIters);	
	GET_OPTION_VALUE(lracLambdaPar,m_lracLambdaPar);
	GET_OPTION_VALUE(lracRadiusPar,m_lracRadiusPar);
	GET_OPTION_VALUE(lracEpsPar,m_lracEpsPar);

	//Hierarchical clustering options
	GET_OPTION_VALUE(spMergingEdgeModel,m_spMergingEdgeModel);
	GET_OPTION_VALUE(spMergingRegPar,m_spMergingRegPar);
	GET_OPTION_VALUE(spMergingNSegmentsToStop,m_spMergingNSegmentsToStop);
	GET_OPTION_VALUE(spMergingRatio,m_spMergingRatio);
	GET_OPTION_VALUE(spMergingMaxDissRatio,m_spMergingMaxDissRatio);
	GET_OPTION_VALUE(spMergingMaxDissRatio2ndNeighbours,m_spMergingMaxDissRatio2ndNeighbours);
	GET_OPTION_VALUE(spMergingDissThreshold,m_spMergingDissThreshold);
	GET_OPTION_VALUE(spMergingIncludeSpatialPars,m_spMergingIncludeSpatialPars);
	GET_OPTION_VALUE(spMergingUseRobustPars,m_spMergingUseRobustPars);
	GET_OPTION_VALUE(spMergingAddCurvDist,m_spMergingAddCurvDist);
	
	return 0;

}//close Configure()


int SFinder::RunTask(TaskData* taskData,bool storeData){

	//Check task data
	if(!taskData){
		ERROR_LOG("Null ptr to task data given!");
		return -1;
	}

	//Get task data
	long int ix_min= taskData->ix_min;
	long int ix_max= taskData->ix_max;
	long int iy_min= taskData->iy_min;
	long int iy_max= taskData->iy_max;

	//Define task data
	Image* taskImg= 0;
	ImgBkgData* bkgData= 0;
	Image* significanceMap= 0;
	Image* segmentedImg= 0;
	bool stopTask= false;
	int status= 0;

	//==================================
	//==   Read task input image
	//==================================
	INFO_LOG("[PROC "<<m_procId<<"] - Reading input image ["<<ix_min<<","<<ix_max<<"] ["<<iy_min<<","<<iy_max<<"]...");
	auto t0_read = chrono::steady_clock::now();	
	
	FileInfo info;
	taskImg= ReadImage(info,m_InputFileName,m_InputImgName,ix_min,ix_max,iy_min,iy_max);
	
	// Set image physical boundary in task data
	if(taskImg){
		double xmin= taskImg->GetXmin();//taskImg->GetXaxis()->GetBinCenter(1);
		double xmax= taskImg->GetXmax();//taskImg->GetXaxis()->GetBinCenter(Nx);
		double ymin= taskImg->GetYmin();//taskImg->GetYaxis()->GetBinCenter(1);
		double ymax= taskImg->GetYmax();//taskImg->GetYaxis()->GetBinCenter(Ny);
		taskData->x_min= xmin;
		taskData->x_max= xmax;
		taskData->y_min= ymin;
		taskData->y_max= ymax;
	}	
	else{
		ERROR_LOG("Reading of input image failed, skip to next task...");
		stopTask= true;
		status= -1;
	}
	auto t1_read = chrono::steady_clock::now();	
	readImageTime+= chrono::duration <double, milli> (t1_read-t0_read).count();
	

	//============================
	//== Find image stats & bkg
	//============================
	if(!stopTask){
		INFO_LOG("[PROC "<<m_procId<<"] - Computing image bkg...");
		bkgData= ComputeStatsAndBkg(taskImg);
		if(!bkgData){
			ERROR_LOG("Failed to compute bkg for input image!");
			stopTask= true;
			status= -1;
		}
	}

	//============================
	//== Find compact sources
	//============================	
	if(!stopTask && m_SearchCompactSources ){ 	
		INFO_LOG("[PROC "<<m_procId<<"] - Searching compact sources...");
		auto t0_sfinder = chrono::steady_clock::now();	

		significanceMap= FindCompactSources(taskImg,bkgData,taskData);
		if(!significanceMap){
			ERROR_LOG("[PROC "<<m_procId<<"] - Compact source search failed!");
			//stopTask= true;
			status= -1;
		}
		auto t1_sfinder = chrono::steady_clock::now();	
		compactSourceTime+= chrono::duration <double, milli> (t1_sfinder-t0_sfinder).count();
	}

	//============================
	//== Find extended sources
	//============================
	if(!stopTask && m_SearchExtendedSources){
		INFO_LOG("Searching extended sources...");
		auto t0_extsfinder = chrono::steady_clock::now();	

		segmentedImg= FindExtendedSources(taskImg,bkgData,taskData,storeData);
		if(!segmentedImg){
			ERROR_LOG("Extended source search failed!");
			status= -1;
		}
		auto t1_extsfinder = chrono::steady_clock::now();	
		extendedSourceTime+= chrono::duration <double, milli> (t1_extsfinder-t0_extsfinder).count();

	}//close if search extended sources


	//## Store data?
	if(storeData){
		if(taskImg) m_InputImg= taskImg;
		if(bkgData) m_BkgData= bkgData;
		if(significanceMap) m_SignificanceMap= significanceMap;
		if(segmentedImg) m_SegmImg= segmentedImg;
	}
	else{//clear all data
		if(taskImg){
			delete taskImg;
			taskImg= 0;
		}
		if(bkgData){
			delete bkgData;
			bkgData= 0;	
		}
		if(significanceMap){
			delete significanceMap;
			significanceMap= 0;
		}
		if(segmentedImg){
			delete segmentedImg;
			segmentedImg= 0;
		}
	}//close else

	return status;

}//close RunTask()



int SFinder::RunMP(){

	//Start timer
	auto t0 = chrono::steady_clock::now();

	//================================================
	//== Init options & data (done by all processors)
	//================================================
	INFO_LOG("[PROC "<<m_procId<<"] - Initializing source finder...");
	auto t0_init = chrono::steady_clock::now();
	if(Init()<0){
		ERROR_LOG("[PROC "<<m_procId<<"] - Initialization failed!");
		return -1;
	}
	auto t1_init = chrono::steady_clock::now();
	initTime= chrono::duration <double, milli> (t1_init-t0_init).count();

	
	//================================================
	//== DISTRIBUTE TASKS TO WORKERS 
	//================================================
	//## Start loop on tasks per worker
	int status= 0;
	bool storeData= true;
	if(m_mpiEnabled) storeData= false;
	size_t nTasks= m_taskDataPerWorkers[m_procId].size(); 
	INFO_LOG("[PROC "<<m_procId<<"] - Start processing of #"<<nTasks<<" tasks...");
	
	for(size_t j=0;j<nTasks;j++){

		//Run task
		if(RunTask(m_taskDataPerWorkers[m_procId][j],storeData)<0){
			ERROR_LOG("[PROC "<<m_procId<<"] - Failed to run task no. "<<j<<", skip to next!");
			status= -1;
			continue;
		}

	}//end loop tasks per worker
	
	if(status<0){
		WARN_LOG("One or more errors occurred in source finding tasks...");
	}

	
	//## Update task info (tile physical range) from workers
	//## The updated list of task data is available in master processor
	//## (check if it is better to replace with MPI_Gather and put it available in all workers)
	#ifdef MPI_ENABLED
	if(m_mpiEnabled && GatherTaskDataFromWorkers()<0){
		ERROR_LOG("[PROC "<<m_procId<<"] - Gathering task data from workers failed!");
		return -1;
	}
	#endif

	//## Merge sources found in each task in unique collection
	if(m_procId==MASTER_ID && MergeTaskData()<0){
		ERROR_LOG("[PROC "<<m_procId<<"] - Merging sources found in each task in unique collection failed!");
		return -1;
	}

	
	//## Find edge sources
	#ifdef MPI_ENABLED
		if(m_mpiEnabled) MPI_Barrier(MPI_COMM_WORLD);
	#endif
	if(m_procId==MASTER_ID) {
		if(FindSourcesAtEdge()<0){
			ERROR_LOG("[PROC "<<m_procId<<"] - Finding sources at tile edges failed!");
			return -1;
		}
	}
	
	
	//## Merge sources at edge (IMPLEMENT ME!!!)
	//...
	//...


	//============================
	//== Deblend sources
	//============================
	if(m_deblendSources && m_procId==MASTER_ID) {
		auto t0_sdeblend = chrono::steady_clock::now();	
		if(DeblendSources(m_SourceCollection)<0){
			ERROR_LOG("Failed to deblend sources!");
			return -1;
		}
		auto t1_sdeblend = chrono::steady_clock::now();	
		sourceDeblendTime= chrono::duration <double, milli> (t1_sdeblend-t0_sdeblend).count();
	}

	//Stop timer
	auto t1 = chrono::steady_clock::now();	
	totTime= chrono::duration <double, milli> (t1-t0).count();

	//============================
	//== Save to file
	//============================
	if(m_saveToFile && m_procId==MASTER_ID) {	
		auto t0_save = chrono::steady_clock::now();	
		Save();	
		auto t1_save = chrono::steady_clock::now();	
		saveTime= chrono::duration <double, milli> (t1_save-t0_save).count();
	}

	//===============================
	//== Print performance stats
	//===============================
	if(m_procId==MASTER_ID) PrintPerformanceStats();


	return 0;

}//close RunMP()



int SFinder::Run(){

	//Start timer
	auto t0 = chrono::steady_clock::now();	

	//===========================
	//== Init options & data
	//===========================
	INFO_LOG("Initializing source finder...");
	auto t0_init = chrono::steady_clock::now();	
	if(Init()<0){
		ERROR_LOG("Initialization failed!");
		return -1;
	}
	auto t1_init = chrono::steady_clock::now();	
	initTime= chrono::duration <double, milli> (t1_init-t0_init).count();
	
	//===========================
	//== Read image
	//===========================
	INFO_LOG("Reading input image...");
	auto t0_read = chrono::steady_clock::now();

	FileInfo info;
	if(m_ReadTile) m_InputImg= ReadImage(info,m_InputFileName,m_InputImgName,m_TileMinX,m_TileMaxX,m_TileMinY,m_TileMaxY);	
	else m_InputImg= ReadImage(info,m_InputFileName,m_InputImgName);
	
	//if(ReadImage()<0){
	if(!m_InputImg){
		ERROR_LOG("Reading of input image failed!");
		return -1;
	}
	m_InputFileExtension= info.extension;

	auto t1_read = chrono::steady_clock::now();	
	readImageTime= chrono::duration <double, milli> (t1_read-t0_read).count();

	//============================
	//== Find compact sources
	//============================
	INFO_LOG("Searching compact sources...");
	auto t0_sfinder = chrono::steady_clock::now();	
	if(m_SearchCompactSources && FindCompactSources()<0){
		ERROR_LOG("Compact source search failed!");
		return -1;
	}
	auto t1_sfinder = chrono::steady_clock::now();	
	compactSourceTime= chrono::duration <double, milli> (t1_sfinder-t0_sfinder).count();


	//============================
	//== Find extended sources
	//============================
	if(m_SearchExtendedSources){
		INFO_LOG("Searching extended sources...");

		// Find residual map
		INFO_LOG("Computing residual image ...");
		auto t0_res = chrono::steady_clock::now();	
		if(FindResidualMap()<0){
			ERROR_LOG("Residual map computation failed!");
			return -1;
		}
		auto t1_res = chrono::steady_clock::now();	
		imgResidualTime= chrono::duration <double, milli> (t1_res-t0_res).count();

		
		//Find extended sources
		auto t0_extsfinder = chrono::steady_clock::now();
		if(FindExtendedSources(m_ResidualImg)<0){
			ERROR_LOG("Extended source search failed!");
			return -1;
		}
		auto t1_extsfinder = chrono::steady_clock::now();	
		extendedSourceTime= chrono::duration <double, milli> (t1_extsfinder-t0_extsfinder).count();

	}//close if search extended sources

	
	//============================
	//== Deblend sources
	//============================
	if(m_deblendSources) {
		auto t0_sdeblend = chrono::steady_clock::now();	
		if(DeblendSources(m_SourceCollection)<0){
			ERROR_LOG("Failed to deblend sources!");
			return -1;
		}
		auto t1_sdeblend = chrono::steady_clock::now();	
		sourceDeblendTime= chrono::duration <double, milli> (t1_sdeblend-t0_sdeblend).count();
	}

	//============================
	//== Draw images
	//============================
	if(m_IsInteractiveRun) {
		DrawSources(m_InputImg,m_SourceCollection);
	}

	//Stop timer
	auto t1 = chrono::steady_clock::now();	
	totTime= chrono::duration <double, milli> (t1-t0).count();

	//============================
	//== Save to file
	//============================
	if(m_saveToFile) {	
		auto t0_save = chrono::steady_clock::now();	
		Save();	
		auto t1_save = chrono::steady_clock::now();	
		saveTime= chrono::duration <double, milli> (t1_save-t0_save).count();
	}

	//===============================
	//== Print performance stats
	//===============================
	PrintPerformanceStats();

	//==========================================
	//== Run TApplication (interactive run)
	//==========================================
	if(m_Application && m_IsInteractiveRun) {
		m_Application->Run();
	}
	
	return 0;

}//close Run()


int SFinder::FindSources(std::vector<Source*>& sources,Image* inputImg,double seedThr,double mergeThr,Image* searchedImg){

	if(!inputImg) {
		ERROR_LOG("Null ptr to given input image!");
		return -1;
	}

	Image* img= inputImg;
	if(searchedImg) img= searchedImg;

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
	INFO_LOG("Finding sources...");	
	int status= inputImg->FindCompactSource(
		sources,
		significanceMap,bkgData,
		seedThr,mergeThr,m_NMinPix,m_SearchNegativeExcess,m_MergeBelowSeed,
		m_SearchNestedSources,m_NestedBlobThrFactor
	);

	//## Clear data
	if(bkgData){
		delete bkgData;
		bkgData= 0;
	}
	if(significanceMap) {
		delete significanceMap;	
		significanceMap= 0;
	}

	//## Check status
	if(status<0) {
		ERROR_LOG("Source finding failed!");	
		return -1;
	}
	int nSources= static_cast<int>(sources.size());
	INFO_LOG("#"<<nSources<<" sources detected in input image...");	
	
	return 0;

}//close FindSources()



Image* SFinder::FindCompactSources(Image* inputImg, ImgBkgData* bkgData, TaskData* taskData){

	//## Check img
	if(!inputImg || !bkgData || !taskData){
		ERROR_LOG("Null ptr to input img and/or bkg/task data!");
		return nullptr;
	}

	//## Compute significance map
	Image* significanceMap= inputImg->GetSignificanceMap(bkgData,m_UseLocalBkg);
	if(!significanceMap){
		ERROR_LOG("Failed to compute significance map!");
		return nullptr;
	}

	//## Find sources
	INFO_LOG("Finding compact sources...");	
	int status= inputImg->FindCompactSource(
		taskData->sources,
		significanceMap,bkgData,
		m_SeedThr,m_MergeThr,m_NMinPix,m_SearchNegativeExcess,m_MergeBelowSeed,
		m_SearchNestedSources,m_NestedBlobThrFactor
	);

	if(status<0) {
		ERROR_LOG("Compact source finding failed!");
		delete significanceMap;
		significanceMap= 0;
		return nullptr;
	}


	//## Apply source selection?
	int nSources= static_cast<int>(taskData->sources.size());
	
	if(m_ApplySourceSelection && nSources>0){
		if(SelectSources(taskData->sources)<0){
			ERROR_LOG("Failed to select sources!");
			delete significanceMap;
			significanceMap= 0;
			return nullptr;
		}
		nSources= static_cast<int>(taskData->sources.size());
	}//close if source selection
		
	INFO_LOG("#"<<nSources<<" bright sources found and added to list...");

	return significanceMap;

}//close FindCompactSources()


int SFinder::FindCompactSources(){

	//## Compute stats and bkg
	INFO_LOG("Computing image stats/bkg...");	
	m_BkgData= ComputeStatsAndBkg(m_InputImg);	
	if(!m_BkgData){
		ERROR_LOG("Failed to compute stats/bkg info!");
		return -1;
	}

	//## Compute significance map
	m_SignificanceMap= m_InputImg->GetSignificanceMap(m_BkgData,m_UseLocalBkg);
	if(!m_SignificanceMap){
		ERROR_LOG("Failed to compute significance map!");
		return -1;
	}

	//## Find sources
	INFO_LOG("Finding compact sources...");	
	std::vector<Source*> sources;
	int status= m_InputImg->FindCompactSource(
		sources,
		m_SignificanceMap,m_BkgData,
		m_SeedThr,m_MergeThr,m_NMinPix,m_SearchNegativeExcess,m_MergeBelowSeed,
		m_SearchNestedSources,m_NestedBlobThrFactor
	);

	if(status<0) {
		ERROR_LOG("Compact source finding failed!");
		return -1;
	}

	//## Retrieve found sources 
	int nSources= static_cast<int>(sources.size());
	INFO_LOG("#"<<nSources<<" bright sources detected in input image...");
	if(nSources<=0) return 0;

	//## Apply source selection?
	int nSelSources= nSources;

	if(m_ApplySourceSelection){
		if(SelectSources(sources)<0){
			ERROR_LOG("Failed to select sources!");
			return -1;
		}
		nSelSources= sources.size();
	}//close if source selection
		
	//## Add detected sources to the list	
	m_SourceCollection.insert(m_SourceCollection.end(),sources.begin(),sources.end());
	m_CompactSources.insert(m_CompactSources.end(),sources.begin(),sources.end());

	INFO_LOG("#"<<nSelSources<<" bright sources to the list...");

	return 0;

}//close FindCompactSources()


Image* SFinder::FindResidualMap(Image* inputImg,ImgBkgData* bkgData,std::vector<Source*> const & sources){

	//Check input image
	if(!inputImg){
		ERROR_LOG("Null ptr to given input image!");
		return nullptr;
	}

	//Compute residual map
	Image* residualImg= 0;
	if(m_UseResidualInExtendedSearch && sources.size()>0){
		INFO_LOG("Computing residual image (#"<<sources.size()<<" sources present)...");
		residualImg= inputImg->GetSourceResidual(
			sources,
			m_DilateKernelSize,m_DilateSourceModel,m_DilatedSourceType,m_DilateNestedSources,	
			bkgData,m_UseLocalBkg,
			m_DilateRandomize,m_SeedBrightThr
		);
	}//close if
	else{
		INFO_LOG("Setting residual image to input image...");
		residualImg= inputImg->GetCloned("",true,true);
	}

	if(!residualImg){
		ERROR_LOG("Failed to compute residual map!");
		return nullptr;
	}
	
	return residualImg;

}//close FindResidualMap()


int SFinder::FindResidualMap(){

	//Compute residual map
	Image* residualImg= 0;
	if(m_UseResidualInExtendedSearch && m_CompactSources.size()>0){
		INFO_LOG("Computing residual image (#"<<m_CompactSources.size()<<" sources present)...");
		residualImg= m_InputImg->GetSourceResidual(
			m_CompactSources,
			m_DilateKernelSize,m_DilateSourceModel,m_DilatedSourceType,m_DilateNestedSources,	
			m_BkgData,m_UseLocalBkg,
			m_DilateRandomize,m_SeedBrightThr
		);
	}//close if
	else{
		INFO_LOG("Setting residual image to input image...");
		residualImg= m_InputImg->GetCloned("residualImg",true,true);
	}

	if(!residualImg){
		ERROR_LOG("Failed to compute residual map!");
		return -1;
	}
	m_ResidualImg= residualImg;
	
	//Compute bkg & noise map for residual img
	INFO_LOG("Computing residual image stats & bkg...");
	m_ResidualBkgData= ComputeStatsAndBkg(m_ResidualImg);
	if(!m_ResidualBkgData){
		ERROR_LOG("Failed to compute bkg data for residual map!");
		return -1;
	}

	return 0;

}//close FindResidualMap()


Image* SFinder::ComputeSmoothedImage(Image* inputImg,int model){

	//Check input image
	if(!inputImg){
		ERROR_LOG("Null ptr to given input image!");
		return nullptr;
	}

	//Compute smoothed image
	Image* smoothedImg= 0;	
	if(model==eGausFilter){
		smoothedImg= inputImg->GetSmoothedImage(m_GausFilterKernSize,m_GausFilterKernSize,m_GausFilterSigma,m_GausFilterSigma);
	}
	else if(model==eGuidedFilter){
		smoothedImg= inputImg->GetGuidedFilterImage(m_GuidedFilterRadius,m_GuidedFilterColorEps);
	}
	else{
		ERROR_LOG("Invalid smoothing algo ("<<model<<") selected!");
		return nullptr;
	}

	if(!smoothedImg){
		ERROR_LOG("Image smoothing failed!");
		return nullptr;
	}

	return smoothedImg;

}//close ComputeSmoothedImage()


Image* SFinder::FindExtendedSources(Image* inputImg,ImgBkgData* bkgData,TaskData* taskData,bool storeData)
{
	//## Check input image
	if(!inputImg || !taskData || !bkgData){
		ERROR_LOG("Null ptr to input image and/or bkg/task data given!");
		return nullptr;
	}

	//## Check segmentation method
	std::vector<int> segmMethods= {eHClust,eActiveContour,eWaveletTransform,eSaliencyThr};
	bool foundMethod= false;
	for(size_t i=0;i<segmMethods.size();i++){
		if(m_ExtendedSearchMethod==segmMethods[i]){
			foundMethod= true;
			break;
		}
	}
	if(!foundMethod){
		ERROR_LOG("Invalid extended source method selected (method="<<m_ExtendedSearchMethod<<")!");
		return nullptr;
	}


	//int status= 0;
	//bool stop= false;
	Image* searchImg= 0;
	
	//****************************
	//** Find residual map
	//****************************
	INFO_LOG("Computing residual image ...");
	Image* residualImg= FindResidualMap(inputImg,bkgData,taskData->sources);
	if(!residualImg){
		ERROR_LOG("Residual map computation failed!");
		return nullptr;
	}
	if(storeData) m_ResidualImg= residualImg;
	searchImg= residualImg;

	/*
	if(residualImg){
		if(storeData) m_ResidualImg= residualImg;
		searchImg= residualImg;
	}
	else{
		ERROR_LOG("Residual map computation failed!");
		return nullptr;
		//stop= true;
		//status= -1;
	}
	*/


	//Compute bkg & noise map for residual img
	INFO_LOG("Computing residual image stats & bkg...");
	ImgBkgData* residualBkgData= ComputeStatsAndBkg(residualImg);
	if(!residualBkgData){
		ERROR_LOG("Failed to compute bkg data for residual map!");
		if(residualImg && !storeData){
			delete residualImg;
			residualImg= 0;
		}
		return nullptr;
	}
	if(storeData) m_ResidualBkgData= residualBkgData;

	/*
	if(!stop && residualImg){
		INFO_LOG("Computing residual image stats & bkg...");
	 	residualBkgData= ComputeStatsAndBkg(residualImg);
		if(residualBkgData){
			if(storeData) m_ResidualBkgData= residualBkgData;
		}
		else{
			ERROR_LOG("Failed to compute bkg data for residual map!");
			stop= true;
			status= -1;		
		}	
	}
	*/

	//****************************
	//** Smooth image
	//****************************
	Image* smoothedImg= 0;
	if(m_UsePreSmoothing){
		smoothedImg= ComputeSmoothedImage(residualImg,m_SmoothFilter);
		if(!smoothedImg){
			ERROR_LOG("Failed to compute residual smoothed image!");
			if(residualImg && !storeData){
				delete residualImg;
				residualImg= 0;
			}
			if(residualBkgData && !storeData){
				delete residualBkgData;
				residualBkgData= 0;
			}
			return nullptr;
		}

		//Set extended source search image to smoothed residual image
		searchImg	= smoothedImg;
	}
	
	/*
	if(!stop && residualImg && m_UsePreSmoothing){
		smoothedImg= ComputeSmoothedImage(residualImg,m_SmoothFilter);
		if(smoothedImg){//Set extended source search image to smoothed residual image
			searchImg	= smoothedImg;
		}
		else{
			ERROR_LOG("Failed to compute residual smoothed image!");
			stop= true;
			status= -1;	
		}
	}//close if smoothing
	*/

	//****************************
	//** Run segmentation
	//****************************
	Image* segmentedImg= 0;
	if(m_ExtendedSearchMethod==eHClust){
		segmentedImg= FindExtendedSources_HClust(inputImg,residualBkgData,taskData,searchImg,storeData);
	}
	else if(m_ExtendedSearchMethod==eActiveContour){
		segmentedImg= FindExtendedSources_AC(inputImg,residualBkgData,taskData,searchImg,storeData);
	}
	else if(m_ExtendedSearchMethod==eWaveletTransform){
		segmentedImg= FindExtendedSources_WT(inputImg,taskData,searchImg);
	}
	else if(m_ExtendedSearchMethod==eSaliencyThr){
		segmentedImg= FindExtendedSources_SalThr(inputImg,residualBkgData,taskData,searchImg,storeData);
	}

	//Check if segmentation succeeded
	if(!segmentedImg){
		ERROR_LOG("Failed to run the segmentation algorithm!");	
		if(residualImg && !storeData){
			delete residualImg;
			residualImg= 0;
		}
		if(residualBkgData && !storeData){
			delete residualBkgData;
			residualBkgData= 0;
		}
		if(smoothedImg){
			delete smoothedImg;
			smoothedImg= 0;
		}
	}

	/*
	if(!stop){
		if(m_ExtendedSearchMethod==eHClust){
			segmentedImg= FindExtendedSources_HClust(inputImg,residualBkgData,taskData,searchImg,storeData);
		}
		else if(m_ExtendedSearchMethod==eActiveContour){
			segmentedImg= FindExtendedSources_AC(inputImg,residualBkgData,taskData,searchImg,storeData);
		}
		else if(m_ExtendedSearchMethod==eWaveletTransform){
			segmentedImg= FindExtendedSources_WT(inputImg,taskData,searchImg);
		}
		else if(m_ExtendedSearchMethod==eSaliencyThr){
			segmentedImg= FindExtendedSources_SalThr(inputImg,residualBkgData,taskData,searchImg,storeData);
		}
		
	
		//Check if segmentation succeeded
		if(!segmentedImg){
			ERROR_LOG("Failed to run the segmentation algorithm!");	
			status= -1;	
		}
	}//close if
	*/

	/*
	//Clear images?
	if(storeData && status==0){
		if(residualImg) m_ResidualImg= residualImg;
		if(residualBkgData) m_ResidualBkgData= residualBkgData;
	}
	else{
		if(residualImg) {
			delete residualImg;
			residualImg= 0;
		}
		if(residualBkgData){
			delete residualBkgData;
			residualBkgData= 0;
		}
		if(smoothedImg){
			delete smoothedImg;		
			smoothedImg= 0;
		}
	}

	//If failed return nullptr
	if(status<0){
		ERROR_LOG("Failures occurred in extended source finding, returning null ptr to segmentation map!");
		return nullptr;
	}
	*/


	//Clear data
	if(residualImg && !storeData) {
		delete residualImg;
		residualImg= 0;
	}
	if(residualBkgData && !storeData) {
		delete residualBkgData;
		residualBkgData= 0;
	}
	if(smoothedImg && !storeData) {
		delete smoothedImg;		
		smoothedImg= 0;
	}
	
	return segmentedImg;

}//close FindExtendedSources()


int SFinder::FindExtendedSources(Image* img){

	//## Check input image
	if(!img){
		ERROR_LOG("Null ptr to input image given!");
		return -1;
	}

	//## Set input map for extended source search
	Image* inputImg= 0;
	Image* smoothedImg= 0;	
	if(m_UsePreSmoothing){//Apply a smoothing stage?
		if(m_SmoothFilter==eGausFilter){
			smoothedImg= img->GetSmoothedImage(m_GausFilterKernSize,m_GausFilterKernSize,m_GausFilterSigma,m_GausFilterSigma);
		}
		else if(m_SmoothFilter==eGuidedFilter){
			smoothedImg= img->GetGuidedFilterImage(m_GuidedFilterRadius,m_GuidedFilterColorEps);
		}
		else{
			ERROR_LOG("Invalid smoothing algo selected!");
			return -1;
		}

		if(!smoothedImg){
			ERROR_LOG("Source residual image smoothing failed!");
			return -1;
		}
		//smoothedImg->SetName("smoothedImg");
		inputImg= smoothedImg;
	}//close if use smoothing
	else{
		inputImg= img;
	}

	//## Run the segmentation	
	int status= 0;
	if(m_ExtendedSearchMethod==eHClust){
		status= FindExtendedSources_HClust(inputImg);
	}
	else if(m_ExtendedSearchMethod==eActiveContour){
		status= FindExtendedSources_AC(inputImg);
	}
	else if(m_ExtendedSearchMethod==eWaveletTransform){
		status= FindExtendedSources_WT(inputImg);
	}
	else if(m_ExtendedSearchMethod==eSaliencyThr){
		status= FindExtendedSources_SalThr(inputImg);
	}
	else{
		ERROR_LOG("Invalid extended source method selected (method="<<m_ExtendedSearchMethod<<")!");
		if(smoothedImg){
			delete smoothedImg;
			smoothedImg= 0;
		} 
		return -1;
	}
	
	if(status<0){
		ERROR_LOG("Failed to run the segmentation algorithm!");	
		if(smoothedImg){
			delete smoothedImg;
			smoothedImg= 0;
		} 
		return -1;
	}

	return 0;

}//close FindExtendedSources()


Image* SFinder::FindExtendedSources_SalThr(Image* inputImg,ImgBkgData* bkgData,TaskData* taskData,Image* searchedImg,bool storeData){

	//Check input image
	if(!inputImg || !bkgData || !taskData){
		ERROR_LOG("Null ptr to input image and/or bkg/task data given!");
		return nullptr;
	}

	Image* img= inputImg;
	if(searchedImg) img= searchedImg;

	//==========================================
	//==    PRELIMINARY STAGES
	//==========================================
	//## Compute saliency
	Image* saliencyImg= img->GetMultiResoSaliencyMap(
		m_SaliencyResoMin,m_SaliencyResoMax,m_SaliencyResoStep,
		m_spBeta,m_spMinArea,m_SaliencyNNFactor,m_SaliencyUseRobustPars,m_SaliencyDissExpFalloffPar,m_SaliencySpatialDistRegPar,
		m_SaliencyMultiResoCombThrFactor,
		m_SaliencyUseBkgMap,m_SaliencyUseNoiseMap,bkgData,
		m_SaliencyThrFactor,m_SaliencyImgThrFactor
	);
	if(saliencyImg){
		if(storeData) m_SaliencyImg= saliencyImg;	
	}
	else{
		ERROR_LOG("Failed to compute saliency map!");
		return nullptr;
	}

	//## Get saliency map optimal threshold
	bool smoothPixelHisto= true;
	int pixelHistoNBins= 100;
	double signalThr= saliencyImg->FindOptimalGlobalThreshold(m_SaliencyThrFactor,pixelHistoNBins,smoothPixelHisto);
	
	//==========================================
	//==    FIND SOURCES
	//==========================================
	//## Find compact blobs in saliency map by simple thresholding
	INFO_LOG("Finding blobs in saliency map with threshold="<<signalThr<<"...");	
	bool findNegativeExcess= false;
	bool mergeBelowSeed= false;
	bool findNestedSources= false;
	int minNPix= m_NMinPix;
	std::vector<Source*> sources;
	int status= inputImg->FindCompactSource(	
		sources, saliencyImg, 0,
		signalThr,signalThr,minNPix,
		findNegativeExcess,mergeBelowSeed,findNestedSources
	);

	//Delete saliency map if not needed
	if(saliencyImg && !storeData){
		delete saliencyImg;
		saliencyImg= 0;
	}

	if(status<0){
		ERROR_LOG("Compact source finding with saliency map failed!");
		return nullptr;
	}

	//## Tag found sources as extended 
	int nSources= static_cast<int>( sources.size() );
	INFO_LOG("#"<<nSources<<" extended sources detected in input image by thresholding the saliency map...");
	for(size_t k=0;k<sources.size();k++) sources[k]->SetType(Source::eExtended);
	
	//## Add sources to extended sources?
	//## NB: Need to decide if to keep them in separate collections (for the moment put in the same collection)
	//(taskData->ext_sources).insert( (taskData->ext_sources).end(),sources.begin(),sources.end());	
	(taskData->sources).insert( (taskData->sources).end(),sources.begin(),sources.end());		
	
	//## Compute segmented map
	bool isBinary= true;
	bool invert= false;
	Image* segmMap= inputImg->GetSourceMask(sources,isBinary,invert);
	if(!segmMap){
		ERROR_LOG("Failed to compute segmented map!");
		return nullptr;
	}
	
	INFO_LOG("#"<<nSources<<" extended sources to the list...");

	return segmMap;

}//close FindExtendedSources_SalThr()


int SFinder::FindExtendedSources_SalThr(Image* inputImg){

	//Check input image
	if(!inputImg){
		ERROR_LOG("Null ptr to input image given!");
		return -1;
	}

	//==========================================
	//==    PRELIMINARY STAGES
	//==========================================
	//## Compute saliency
	m_SaliencyImg= inputImg->GetMultiResoSaliencyMap(
		m_SaliencyResoMin,m_SaliencyResoMax,m_SaliencyResoStep,
		m_spBeta,m_spMinArea,m_SaliencyNNFactor,m_SaliencyUseRobustPars,m_SaliencyDissExpFalloffPar,m_SaliencySpatialDistRegPar,
		m_SaliencyMultiResoCombThrFactor,
		m_SaliencyUseBkgMap,m_SaliencyUseNoiseMap,m_ResidualBkgData,
		m_SaliencyThrFactor,m_SaliencyImgThrFactor
	);
	if(!m_SaliencyImg){
		ERROR_LOG("Failed to compute saliency map!");
		return -1;
	}

	//## Get saliency map optimal threshold
	bool smoothPixelHisto= true;
	int pixelHistoNBins= 100;
	double signalThr= m_SaliencyImg->FindOptimalGlobalThreshold(m_SaliencyThrFactor,pixelHistoNBins,smoothPixelHisto);
	
	//==========================================
	//==    FIND SOURCES
	//==========================================
	//## Find compact blobs in saliency map by simple thresholding
	INFO_LOG("Finding blobs in saliency map with threshold="<<signalThr<<"...");	
	bool findNegativeExcess= false;
	bool mergeBelowSeed= false;
	bool findNestedSources= false;
	int minNPix= m_NMinPix;
	std::vector<Source*> sources;

	int status= m_InputImg->FindCompactSource(	
		sources, m_SaliencyImg, 0,
		signalThr,signalThr,minNPix,
		findNegativeExcess,mergeBelowSeed,findNestedSources
	);

	if(status<0){
		ERROR_LOG("Compact source finding with saliency map failed!");
		return -1;
	}

	//## Retrieve found sources 
	int nSources= static_cast<int>(sources.size());
	INFO_LOG("#"<<nSources<<" extended sources detected in input image by thresholding the saliency map...");
	if(nSources<=0) return 0;

	//## Apply source selection?
	int nSelSources= nSources;
	//...
	//...
	/*
	if(m_ApplySourceSelection){
		if(SelectSources(sources)<0){
			ERROR_LOG("Failed to select sources!");
			return -1;
		}
		nSelSources= sources.size();
	}//close if source selection
	*/

	//## Compute segmented map
	bool isBinary= true;
	bool invert= false;
	m_SegmImg= m_InputImg->GetSourceMask(sources,isBinary,invert);
	if(!m_SegmImg){
		ERROR_LOG("Failed to compute segmented map!");
		return -1;
	}

	//## Add detected sources to the list	
	m_SourceCollection.insert(m_SourceCollection.end(),sources.begin(),sources.end());
	m_ExtendedSources.insert(m_ExtendedSources.end(),sources.begin(),sources.end());

	INFO_LOG("#"<<nSelSources<<" extended sources to the list...");

	return 0;

}//close FindExtendedSources_SalThr()


Image* SFinder::FindExtendedSources_HClust(Image* inputImg,ImgBkgData* bkgData,TaskData* taskData,Image* searchedImg,bool storeData){

	//Check input image
	if(!inputImg || !bkgData || !taskData){
		ERROR_LOG("Null ptr to input image and/or bkg/task data given!");
		return nullptr;
	}

	Image* img= inputImg;
	if(searchedImg) img= searchedImg;

	//==========================================
	//==    PRELIMINARY STAGES
	//==========================================
	//## Compute saliency
	Image* saliencyImg= inputImg->GetMultiResoSaliencyMap(
		m_SaliencyResoMin,m_SaliencyResoMax,m_SaliencyResoStep,
		m_spBeta,m_spMinArea,m_SaliencyNNFactor,m_SaliencyUseRobustPars,m_SaliencyDissExpFalloffPar,m_SaliencySpatialDistRegPar,
		m_SaliencyMultiResoCombThrFactor,
		m_SaliencyUseBkgMap,m_SaliencyUseNoiseMap,bkgData,
		m_SaliencyThrFactor,m_SaliencyImgThrFactor
	);
	if(saliencyImg){
		if(storeData) m_SaliencyImg= saliencyImg;
 	}
	else{
		ERROR_LOG("Failed to compute saliency map!");
		return nullptr;
	}

	//## Threshold saliency map and get signal and bkg markers
	bool smoothPixelHisto= true;
	int pixelHistoNBins= 100;
	double signalThr= saliencyImg->FindOptimalGlobalThreshold(m_SaliencyThrFactor,pixelHistoNBins,smoothPixelHisto);
	double bkgThr= saliencyImg->FindMedianThreshold(m_SaliencyBkgThrFactor);
	if(TMath::IsNaN(signalThr) || fabs(signalThr)==TMath::Infinity()){
		ERROR_LOG("Invalid numeric threshold returned as threshold computation failed!");
		return nullptr;
	}

	INFO_LOG("Computing binarized saliency maps (signalThr="<<signalThr<<", bkgThr="<<bkgThr);
	double fgValue= 1;
	Image* signalMarkerImg= saliencyImg->GetBinarizedImage(signalThr,fgValue,false);
	Image* bkgMarkerImg= saliencyImg->GetBinarizedImage(bkgThr,fgValue,true);
	
	//## Delete saliency map if not needed
	if(saliencyImg && !storeData){
		delete saliencyImg;
		saliencyImg= 0;
	}

	//## Compute Laplacian filtered image
	INFO_LOG("Computing laplacian image...");
	Image* laplImg= ComputeLaplacianImage(img);
	if(m_LaplImg){
		if(storeData) m_LaplImg= laplImg;
	}
	else{
		ERROR_LOG("Failed to compute laplacian image, cannot perform extended source finding!");
		return nullptr;
	}

	//## Compute edge image	
	INFO_LOG("Computing edgeness image...");
	Image* edgeImg= ComputeEdgeImage(img,m_spMergingEdgeModel);
	if(edgeImg){
		if(storeData) m_EdgeImg= edgeImg;
	}
	else{
		ERROR_LOG("Failed to compute the edgeness image, cannot perform extended source finding!");
		if(laplImg && !storeData){
			delete laplImg;		
			laplImg= 0;
		}
		return nullptr;
	}

	//## Compute the Superpixel partition
	bool normalizeImage= true;
	SLICData* slicData_init= SLIC::SPGenerator(img,m_spSize,m_spBeta,m_spMinArea,normalizeImage,m_spUseLogContrast,laplImg,edgeImg);
	if(!slicData_init){
		ERROR_LOG("Failed to compute the initial superpixel partition, cannot perform extended source finding!");	
		if(laplImg && !storeData){
			delete laplImg;		
			laplImg= 0;
		}
		if(edgeImg && !storeData){
			delete edgeImg;		
			edgeImg= 0;
		}
		return nullptr;
	}

	//## Tag the superpixel partition
	if(SLIC::TagRegions(slicData_init->regions,bkgMarkerImg,signalMarkerImg)<0){
		ERROR_LOG("Failed to tag (signal vs bkg) the initial superpixel partition, cannot perform extended source finding!");
		delete slicData_init;
		slicData_init= 0;
		if(laplImg && !storeData){
			delete laplImg;		
			laplImg= 0;
		}
		if(edgeImg && !storeData){
			delete edgeImg;		
			edgeImg= 0;
		}
		if(bkgMarkerImg){
			delete bkgMarkerImg;
			bkgMarkerImg= 0;
		}	
		if(signalMarkerImg){
			delete signalMarkerImg;
			signalMarkerImg= 0;
		}
		return nullptr;
	}

	//## Delete signal/bkg marker images
	if(bkgMarkerImg){
		delete bkgMarkerImg;
		bkgMarkerImg= 0;
	}	
	if(signalMarkerImg){
		delete signalMarkerImg;
		signalMarkerImg= 0;
	}

	//==========================================
	//==    RUN SEGMENTATION
	//==========================================
	//## Run the segmentation
	INFO_LOG("Running the hierarchical clustering segmenter...");	
	SLICData slicData_segm;
	SLICSegmenter::FindSegmentation(
		*slicData_init, slicData_segm,
		m_spMergingRegPar, m_spMergingUse2ndNeighbours,
		m_spMergingNSegmentsToStop,m_spMergingRatio,
		m_spMergingMaxDissRatio,m_spMergingMaxDissRatio2ndNeighbours,m_spMergingDissThreshold,
		m_spMergingIncludeSpatialPars, m_spMergingUseRobustPars, m_spMergingAddCurvDist
	);

	//## Get segmentation results	
	INFO_LOG("Computing the segmented map from slic segmented data...");
	bool normalizeSegmImg= true;
	bool binarizeSegmImg= true;
	Image* segmentedImg= SLIC::GetSegmentedImage(inputImg,slicData_segm.regions,Region::eSignalTag,normalizeSegmImg,binarizeSegmImg);
	if(!segmentedImg){
		ERROR_LOG("Failed to compute the segmented image from slic segmented data!");
		delete slicData_init;
		slicData_init= 0;
		if(laplImg && !storeData){
			delete laplImg;		
			laplImg= 0;
		}
		if(edgeImg && !storeData){
			delete edgeImg;		
			edgeImg= 0;
		}
		return nullptr;
	}

	//## Clear-up
	if(slicData_init){
		delete slicData_init;
		slicData_init= 0;
	}
	if(laplImg && !storeData){
		delete laplImg;		
		laplImg= 0;
	}
	if(edgeImg && !storeData){
		delete edgeImg;		
		edgeImg= 0;
	}


	//## Finding blobs in masked image
	bool findNegativeExcess= false;
	bool mergeBelowSeed= false;
	bool findNestedSources= false;
	std::vector<Source*> sources;
	int status= inputImg->FindCompactSource(
		sources, segmentedImg,
		bkgData, fgValue, fgValue, 
		m_NMinPix, findNegativeExcess, mergeBelowSeed, findNestedSources
	);
	if(status<0){
		ERROR_LOG("Finding sources in hierarchical algorithm segmented mask failed!");
		return nullptr;
	}

	//## Tag sources as extended
	for(size_t k=0;k<sources.size();k++) sources[k]->SetType(Source::eExtended);
	
	//## Remove sources of negative excess (THIS METHOD SHOULD BE IMPROVED)
	if(inputImg->HasStats()){
		ImgStats* stats= inputImg->GetPixelStats();
		double imgMedian= stats->median;

		std::vector<size_t> sourcesToBeRemoved;				
		for(size_t k=0;k<sources.size();k++){	
			double Smedian= sources[k]->Median;
			if(Smedian<imgMedian) sourcesToBeRemoved.push_back(k);
		}
		CodeUtils::DeleteItems(sources, sourcesToBeRemoved);

	}//close if
	else {
		WARN_LOG("Input image has no stats computed (hint: you must have computed them before!), cannot remove negative excess from sources!");
	}	

	//## Add sources to extended sources
	//## NB: Need to decide if to keep them in separate collections (for the moment put in the same collection)
	//(taskData->ext_sources).insert( (taskData->ext_sources).end(),sources.begin(),sources.end());		
	(taskData->sources).insert( (taskData->sources).end(),sources.begin(),sources.end());		
	

	return segmentedImg;

}//close FindExtendedSources_HClust()


int SFinder::FindExtendedSources_HClust(Image* inputImg){

	//Check input image
	if(!inputImg){
		ERROR_LOG("Null ptr to input image given!");
		return -1;
	}

	//==========================================
	//==    PRELIMINARY STAGES
	//==========================================
	//## Compute saliency
	m_SaliencyImg= inputImg->GetMultiResoSaliencyMap(
		m_SaliencyResoMin,m_SaliencyResoMax,m_SaliencyResoStep,
		m_spBeta,m_spMinArea,m_SaliencyNNFactor,m_SaliencyUseRobustPars,m_SaliencyDissExpFalloffPar,m_SaliencySpatialDistRegPar,
		m_SaliencyMultiResoCombThrFactor,
		m_SaliencyUseBkgMap,m_SaliencyUseNoiseMap,m_ResidualBkgData,
		m_SaliencyThrFactor,m_SaliencyImgThrFactor
	);
	if(!m_SaliencyImg){
		ERROR_LOG("Failed to compute saliency map!");
		return -1;
	}

	//## Threshold saliency map and get signal and bkg markers
	bool smoothPixelHisto= true;
	int pixelHistoNBins= 100;
	double signalThr= m_SaliencyImg->FindOptimalGlobalThreshold(m_SaliencyThrFactor,pixelHistoNBins,smoothPixelHisto);
	double bkgThr= m_SaliencyImg->FindMedianThreshold(m_SaliencyBkgThrFactor);
	
	INFO_LOG("Computing binarized saliency maps (signalThr="<<signalThr<<", bkgThr="<<bkgThr);
	double fgValue= 1;
	Image* signalMarkerImg= m_SaliencyImg->GetBinarizedImage(signalThr,fgValue,false);
	Image* bkgMarkerImg= m_SaliencyImg->GetBinarizedImage(bkgThr,fgValue,true);
	
	//## Compute Laplacian filtered image
	INFO_LOG("Computing laplacian image...");
	m_LaplImg= ComputeLaplacianImage(inputImg);
	if(!m_LaplImg){
		ERROR_LOG("Failed to compute laplacian image, cannot perform extended source finding!");
		return -1;
	}

	//## Compute edge image	
	INFO_LOG("Computing edgeness image...");
	m_EdgeImg= ComputeEdgeImage(inputImg,m_spMergingEdgeModel);
	if(!m_EdgeImg){
		ERROR_LOG("Failed to compute the edgeness image, cannot perform extended source finding!");
		return -1;
	}

	//## Compute the Superpixel partition
	//SLICData* slicData_init= SLIC::SPGenerator(inputImg,m_spSize,m_spBeta,m_spMinArea,m_spUseLogContrast,m_EdgeImg);
	bool normalizeImage= true;
	SLICData* slicData_init= SLIC::SPGenerator(inputImg,m_spSize,m_spBeta,m_spMinArea,normalizeImage,m_spUseLogContrast,m_LaplImg,m_EdgeImg);
	if(!slicData_init){
		ERROR_LOG("Failed to compute the initial superpixel partition, cannot perform extended source finding!");	
		return -1;
	}

	//## Tag the superpixel partition
	if(SLIC::TagRegions(slicData_init->regions,bkgMarkerImg,signalMarkerImg)<0){
		ERROR_LOG("Failed to tag (signal vs bkg) the initial superpixel partition, cannot perform extended source finding!");
		delete slicData_init;
		slicData_init= 0;
		return -1;
	}
	
	//==========================================
	//==    RUN SEGMENTATION
	//==========================================
	//## Run the segmentation
	INFO_LOG("Running the hierarchical clustering segmenter...");	
	SLICData slicData_segm;
	SLICSegmenter::FindSegmentation(
		*slicData_init, slicData_segm,
		m_spMergingRegPar, m_spMergingUse2ndNeighbours,
		m_spMergingNSegmentsToStop,m_spMergingRatio,
		m_spMergingMaxDissRatio,m_spMergingMaxDissRatio2ndNeighbours,m_spMergingDissThreshold,
		m_spMergingIncludeSpatialPars, m_spMergingUseRobustPars, m_spMergingAddCurvDist
	);

	//## Get segmentation results	
	INFO_LOG("Computing the segmented map from slic segmented data...");
	bool normalizeSegmImg= true;
	bool binarizeSegmImg= true;
	m_SegmImg= SLIC::GetSegmentedImage(inputImg,slicData_segm.regions,Region::eSignalTag,normalizeSegmImg,binarizeSegmImg);
	if(!m_SegmImg){
		ERROR_LOG("Failed to compute the segmented image from slic segmented data!");
		delete slicData_init;
		slicData_init= 0;
		return -1;
	}

	//## Clear-up
	delete slicData_init;
	slicData_init= 0;

	return 0;

}//close FindExtendedSources_HClust()

Image* SFinder::ComputeLaplacianImage(Image* inputImg){

	//Check input image
	if(!inputImg){
		ERROR_LOG("Null ptr to given input image!");
		return nullptr;
	}

	//Compute laplacian image
	INFO_LOG("Computing Laplacian image ...");
	Image* laplImg= inputImg->GetLaplacianImage(true);
	if(!laplImg){
		ERROR_LOG("Failed to compute Laplacian image!");
		return nullptr;
	}

	//Compute laplacian image stats
	INFO_LOG("Compute Laplacian image stats...");
	if(laplImg->ComputeStats(true,false,false)<0){	
		ERROR_LOG("Failed to compute Laplacian image stats, returning nullptr!");
		delete laplImg;
		laplImg= 0;
		return nullptr;
	}
	INFO_LOG("Laplacian image stats");
	laplImg->PrintStats();

	return laplImg;

}//close ComputeLaplacianImage()

Image* SFinder::ComputeEdgeImage(Image* inputImg,int edgeModel){

	//Check input image
	if(!inputImg){
		ERROR_LOG("Null ptr to given input image!");
		return nullptr;
	}

	//Compute edge image according to desired model
	Image* edgeImg= 0;
	if(edgeModel == eKirschEdge){
		INFO_LOG("Computing edge image using a Kirsch model...");
		edgeImg= inputImg->GetKirschImage();	
	}
	else if(edgeModel == eChanVeseEdge){
		INFO_LOG("Computing edge image using a Chan-Vese contour model...");	
		bool returnContourImg= true;
		edgeImg= ChanVeseSegmenter::FindSegmentation (
			inputImg, 0, returnContourImg,
			m_cvTimeStepPar, m_cvWindowSizePar, m_cvLambda1Par, m_cvLambda2Par, m_cvMuPar, m_cvNuPar, m_cvPPar, m_cvNIters
		);
	}
	else {
		ERROR_LOG("Invalid edge model selected!");
		return nullptr;
	}

	//Check if edge image computing failed
	if(!edgeImg){
		ERROR_LOG("Failed to compute edge image!");
		return nullptr;
	}

	//Compute edge image stats
	INFO_LOG("Compute edge image stats...");
	if(edgeImg->ComputeStats(true,false,false)<0){
		ERROR_LOG("Failed to compute edge image stats, returning nullptr!");
		delete edgeImg;
		edgeImg= 0;
		return nullptr;
	}

	INFO_LOG("Edgeness image stats");
	edgeImg->PrintStats();
	
	return edgeImg;

}//close ComputeEdgeImage()

Image* SFinder::FindExtendedSources_AC(Image* inputImg,ImgBkgData* bkgData,TaskData* taskData,Image* searchedImg,bool storeData){

	//## Check input image
	if(!inputImg || !bkgData || !taskData){
		ERROR_LOG("Null ptr to input image and/or to bkg/task data given!");
		return nullptr;
	}

	Image* img= inputImg;
	if(searchedImg) img= searchedImg;

	INFO_LOG("Searching extended sources with the active contour method...");

	//==========================================
	//==    PRELIMINARY STAGES
	//==========================================	
	//## Compute saliency
	INFO_LOG("Computing image saliency map...");
	Image* saliencyImg= img->GetMultiResoSaliencyMap(
		m_SaliencyResoMin,m_SaliencyResoMax,m_SaliencyResoStep,
		m_spBeta,m_spMinArea,m_SaliencyNNFactor,m_SaliencyUseRobustPars,m_SaliencyDissExpFalloffPar,m_SaliencySpatialDistRegPar,
		m_SaliencyMultiResoCombThrFactor,
		m_SaliencyUseBkgMap,m_SaliencyUseNoiseMap,bkgData,
		m_SaliencyThrFactor,m_SaliencyImgThrFactor
	);
	if(saliencyImg){
		if(storeData) m_SaliencyImg= saliencyImg;
	}
	else{
		ERROR_LOG("Failed to compute saliency map!");
		return nullptr;
	}
	
	//## Get saliency map optimal threshold
	INFO_LOG("Computing saliency map optimal threshold...");
	bool smoothPixelHisto= true;
	int pixelHistoNBins= 100;
	double signalThr= saliencyImg->FindOptimalGlobalThreshold(m_SaliencyThrFactor,pixelHistoNBins,smoothPixelHisto);
	if(TMath::IsNaN(signalThr) || fabs(signalThr)==TMath::Infinity()){
		ERROR_LOG("Invalid numeric threshold returned as threshold computation failed!");
		return nullptr;
	}	

	//## Get saliency binarized image
	INFO_LOG("Thresholding the saliency map @ thr="<<signalThr<<" and compute binarized map...");
	double fgValue= 1;
	Image* signalMarkerImg= saliencyImg->GetBinarizedImage(signalThr,fgValue,false);
	if(!signalMarkerImg){
		ERROR_LOG("Failed to get saliency binarized map!");
		if(saliencyImg && !storeData){
			delete saliencyImg;
			saliencyImg= 0;
		}
		return nullptr;
	}
		
	//Delete saliency if not needed
	if(saliencyImg && !storeData){
		delete saliencyImg;
		saliencyImg= 0;
	}

	//## If binarized mage is empty (e.g. only background) do not run contour algorithm
	if(signalMarkerImg->GetMaximum()<=0){
		WARN_LOG("No signal objects detected in saliency map (only background), will not run active contour (NB: no extended sources detected in this image!)");
		delete signalMarkerImg;
		signalMarkerImg= 0;
		return nullptr;
	}

	//==========================================
	//==    RUN ACTIVE CONTOUR SEGMENTATION
	//==========================================

	//## Compute segmented image
	Image* segmentedImg= 0;
	bool returnContourImg= false;
	if(m_activeContourMethod==eChanVeseAC){//Standard ChanVese algorithm
		segmentedImg= ChanVeseSegmenter::FindSegmentation (
			img, signalMarkerImg, returnContourImg,
			m_cvTimeStepPar,m_cvWindowSizePar,m_cvLambda1Par,m_cvLambda2Par,m_cvMuPar,m_cvNuPar,m_cvPPar,m_cvNIters
		);
	}
	else if(m_activeContourMethod==eLRAC){//LRAC algorithm (with Chan-vese energy)
		segmentedImg= LRACSegmenter::FindSegmentation (
			img, signalMarkerImg,
			m_lracNIters,m_lracLambdaPar,m_lracRadiusPar,m_lracEpsPar
		);
	}
	else{
		ERROR_LOG("Invalid active contour method specified ("<<m_activeContourMethod<<")!");
		delete signalMarkerImg;
		signalMarkerImg= 0;
		return nullptr;
	}

	//Delete signal mask
	if(signalMarkerImg){
		delete signalMarkerImg;
		signalMarkerImg= 0;
	}

	if(!segmentedImg){
		ERROR_LOG("Failed to compute Active Contour image segmentation!");
		return nullptr;
	}
	
	//## Finding blobs in masked image
	bool findNegativeExcess= false;
	bool mergeBelowSeed= false;
	bool findNestedSources= false;
	std::vector<Source*> sources;
	int status= inputImg->FindCompactSource(
		sources, segmentedImg,
		bkgData, fgValue, fgValue, 
		m_NMinPix, findNegativeExcess, mergeBelowSeed, findNestedSources
	);
	if(status<0){
		ERROR_LOG("Finding sources in active contour segmented mask failed!");
		return nullptr;
	}

	//## Tag sources as extended
	for(size_t k=0;k<sources.size();k++) sources[k]->SetType(Source::eExtended);
	
	//## Remove sources of negative excess (because Chan-Vese detects them) (THIS METHOD SHOULD BE IMPROVED)
	if(inputImg->HasStats()){
		ImgStats* stats= inputImg->GetPixelStats();
		double imgMedian= stats->median;

		std::vector<size_t> sourcesToBeRemoved;				
		for(size_t k=0;k<sources.size();k++){	
			double Smedian= sources[k]->Median;
			if(Smedian<imgMedian) sourcesToBeRemoved.push_back(k);
		}
		CodeUtils::DeleteItems(sources, sourcesToBeRemoved);

	}//close if
	else {
		WARN_LOG("Input image has no stats computed (hint: you must have computed them before!), cannot remove negative excess from sources!");
	}	

	//## Add sources to extended sources
	//## NB: Need to decide if to keep them in separate collections (for the moment put in the same collection)
	//(taskData->ext_sources).insert( (taskData->ext_sources).end(),sources.begin(),sources.end());		
	(taskData->sources).insert( (taskData->sources).end(),sources.begin(),sources.end());		
	
	return segmentedImg;

}//close FindExtendedSources_AC()


int SFinder::FindExtendedSources_AC(Image* inputImg){

	//## Check input image
	if(!inputImg){
		ERROR_LOG("Null ptr to input image given!");
		return -1;
	}

	INFO_LOG("Searching extended sources with the active contour method...");

	//==========================================
	//==    PRELIMINARY STAGES
	//==========================================
	//## Compute saliency
	m_SaliencyImg= inputImg->GetMultiResoSaliencyMap(
		m_SaliencyResoMin,m_SaliencyResoMax,m_SaliencyResoStep,
		m_spBeta,m_spMinArea,m_SaliencyNNFactor,m_SaliencyUseRobustPars,m_SaliencyDissExpFalloffPar,m_SaliencySpatialDistRegPar,
		m_SaliencyMultiResoCombThrFactor,
		m_SaliencyUseBkgMap,m_SaliencyUseNoiseMap,m_ResidualBkgData,
		m_SaliencyThrFactor,m_SaliencyImgThrFactor
	);
	if(!m_SaliencyImg){
		ERROR_LOG("Failed to compute saliency map!");
		return -1;
	}

	//## Get saliency map optimal threshold
	bool smoothPixelHisto= true;
	int pixelHistoNBins= 100;
	double signalThr= m_SaliencyImg->FindOptimalGlobalThreshold(m_SaliencyThrFactor,pixelHistoNBins,smoothPixelHisto);
	
	//## Get saliency binarized image
	double fgValue= 1;
	Image* signalMarkerImg= m_SaliencyImg->GetBinarizedImage(signalThr,fgValue,false);
	if(!signalMarkerImg){
		ERROR_LOG("Failed to get saliency binarized map!");
		return -1;
	}

	//## If binarized mage is empty (e.g. only background) do not run contour algorithm
	if(signalMarkerImg->GetMaximum()<=0){
		WARN_LOG("No signal objects detected in saliency map (only background), will not run active contour (NB: no extended sources detected in this image!)");
		return -1;
	}

	//==========================================
	//==    RUN ACTIVE CONTOUR SEGMENTATION
	//==========================================

	//## Compute segmented image
	bool returnContourImg= false;
	if(m_activeContourMethod==eChanVeseAC){//Standard ChanVese algorithm
		m_SegmImg= ChanVeseSegmenter::FindSegmentation (
			inputImg, signalMarkerImg, returnContourImg,
			m_cvTimeStepPar,m_cvWindowSizePar,m_cvLambda1Par,m_cvLambda2Par,m_cvMuPar,m_cvNuPar,m_cvPPar,m_cvNIters
		);
	}
	else if(m_activeContourMethod==eLRAC){//LRAC algorithm (with Chan-vese energy)
		m_SegmImg= LRACSegmenter::FindSegmentation (
			inputImg, signalMarkerImg,
			m_lracNIters,m_lracLambdaPar,m_lracRadiusPar,m_lracEpsPar
		);
	}
	else{
		ERROR_LOG("Invalid active contour method specified ("<<m_activeContourMethod<<")!");
		delete signalMarkerImg;
		signalMarkerImg= 0;
		return -1;
	}

	if(!m_SegmImg){
		ERROR_LOG("Failed to compute ChanVese image segmentation!");
		delete signalMarkerImg;
		signalMarkerImg= 0;
		return -1;
	}
	
	//## Finding blobs in masked image
	bool findNegativeExcess= false;
	bool mergeBelowSeed= false;
	bool findNestedSources= false;
	std::vector<Source*> sources;
	int status= m_InputImg->FindCompactSource(
		sources, m_SegmImg,
		m_BkgData, fgValue, fgValue, 
		m_NMinPix, findNegativeExcess, mergeBelowSeed, findNestedSources
	);
	if(status<0){
		ERROR_LOG("Finding sources in Chan-Vese segmented mask failed!");
		return -1;
	}

	//## Tag sources as extended
	for(size_t k=0;k<sources.size();k++) sources[k]->SetType(Source::eExtended);
	
	//## Remove sources of negative excess (because Chan-Vese detects them) (THIS METHOD SHOULD BE IMPROVED)
	if(m_InputImg->HasStats()){
		ImgStats* stats= m_InputImg->GetPixelStats();
		double imgMedian= stats->median;

		std::vector<size_t> sourcesToBeRemoved;				
		for(size_t k=0;k<sources.size();k++){	
			double Smedian= sources[k]->Median;
			if(Smedian<imgMedian) sourcesToBeRemoved.push_back(k);
		}
		CodeUtils::DeleteItems(sources, sourcesToBeRemoved);

	}//close if
	else {
		WARN_LOG("Input image has no stats computed (hint: you must have computed them before!), cannot remove negative excess from sources!");
	}	

	//## Add sources to extended sources
	m_ExtendedSources.insert(m_ExtendedSources.end(),sources.begin(),sources.end());		
	m_SourceCollection.insert(m_SourceCollection.end(),sources.begin(),sources.end());

	return 0;

}//close FindExtendedSources_AC()





Image* SFinder::FindExtendedSources_WT(Image* inputImg,TaskData* taskData,Image* searchedImg){

	//## Check input image
	if(!inputImg || !taskData){
		ERROR_LOG("Null ptr to input image and/or task data given!");
		return nullptr;
	}
	Image* img= inputImg;
	if(searchedImg) img= searchedImg;
	
	//## Find extended sources in the scales of the residual image (with POINT-LIKE SOURCES removed)
	INFO_LOG("Find extended sources in the residual image WT-"<<m_wtScaleExtended<<"  scale ...");
	std::vector<Image*> wt_extended= img->GetWaveletDecomposition(m_wtScaleExtended);
	
	std::vector<Source*> sources;
	int status= FindSources(
		sources,
		inputImg,
		m_SeedThr,m_MergeThr,
		wt_extended[m_wtScaleExtended]
	);

	//## Clear-up
	for(size_t i=0;i<wt_extended.size();i++){
		if(wt_extended[i]) {
			delete wt_extended[i];
			wt_extended[i]= 0;
		}
	}

	if(status<0){
		ERROR_LOG("Extended source finding failed!");
		return nullptr;
	}
	
	
	//## Tag sources as extended
	int nSources= static_cast<int>( sources.size() );		
	INFO_LOG("#"<<nSources<<" found...");

	for(size_t i=0;i<sources.size();i++){
		sources[i]->SetType(Source::eExtended);
	}

	//## Add sources to extended sources
	//## NB: Need to decide if to keep them in separate collections (for the moment put in the same collection)
	//(taskData->ext_sources).insert( (taskData->ext_sources).end(),sources.begin(),sources.end());		
	(taskData->sources).insert( (taskData->sources).end(),sources.begin(),sources.end());		
	

	//## Compute segmented map
	bool isBinary= false;
	bool invert= false;
	Image* segmMap= inputImg->GetSourceMask(sources,isBinary,invert);
	if(!segmMap){
		ERROR_LOG("Failed to compute segmented map!");
		return nullptr;
	}

	return segmMap;

}//close FindExtendedSources_WT()



int SFinder::FindExtendedSources_WT(Image* inputImg){

	//## Check input image
	if(!inputImg){
		ERROR_LOG("Null ptr to input image given!");
		return -1;
	}
	
	//## Find extended sources in the scales of the residual image (with POINT-LIKE SOURCES removed)
	INFO_LOG("Find extended sources in the residual image WT-"<<m_wtScaleExtended<<"  scale ...");
	std::vector<Image*> wt_extended= inputImg->GetWaveletDecomposition(m_wtScaleExtended);
	
	std::vector<Source*> sources;
	int status= FindSources(
		sources, 
		wt_extended[m_wtScaleExtended],
		m_SeedThr,m_MergeThr
	);


	//## Clear-up
	for(size_t i=0;i<wt_extended.size();i++){
		if(wt_extended[i]) {
			delete wt_extended[i];
			wt_extended[i]= 0;
		}
	}

	if(status<0){
		ERROR_LOG("Extended source finding failed!");	
		return -1;
	}

	//## Tag sources as extended
	int nSources= (int)sources.size();		
	INFO_LOG("#"<<nSources<<" found...");

	for(size_t i=0;i<sources.size();i++){
		sources[i]->SetType(Source::eExtended);
	}

	//## Add sources to extended sources
	m_ExtendedSources.insert(m_ExtendedSources.end(),sources.begin(),sources.end());		
	m_SourceCollection.insert(m_SourceCollection.end(),sources.begin(),sources.end());		
		
	//## Compute segmented map
	bool isBinary= false;
	bool invert= false;
	m_SegmImg= inputImg->GetSourceMask(sources,isBinary,invert);
	if(!m_SegmImg){
		ERROR_LOG("Failed to compute segmented map!");
		return -1;
	}

	return 0;

}//close FindExtendedSources_WT()

int SFinder::SelectSources(std::vector<Source*>& sources){

	//## Apply source selection?
	int nSources= static_cast<int>(sources.size());
	if(nSources<=0) return 0;
	
	int nSelSources= 0;
	std::vector<Source*> sources_sel;

	for(int i=0;i<nSources;i++){	
		std::string sourceName= sources[i]->GetName();
		int sourceId= sources[i]->Id;
		long int NPix= sources[i]->NPix;
		double X0= sources[i]->X0;
		double Y0= sources[i]->Y0;

		//Is bad source (i.e. line-like blob, etc...)?
		if(!IsGoodSource(sources[i])) {
			INFO_LOG("Source no. "<<i<<" (name="<<sourceName<<",id="<<sourceId<<", n="<<NPix<<"("<<X0<<","<<Y0<<")) tagged as bad source, skipped!");
			sources[i]->SetGoodSourceFlag(false);
			continue;
		}
			
		//Is point-like source?
		if( IsPointLikeSource(sources[i]) ){
			INFO_LOG("Source no. "<<i<<" (name="<<sourceName<<",id="<<sourceId<<", n="<<NPix<<"("<<X0<<","<<Y0<<")) tagged as a point-like source ...");
			sources[i]->SetType(Source::ePointLike);
		}

		//Tag nested sources
		std::vector<Source*> nestedSources= sources[i]->GetNestedSources();
		for(size_t j=0;j<nestedSources.size();j++){
			std::string nestedSourceName= nestedSources[j]->GetName();
			int nestedSourceId= nestedSources[j]->Id;
			long int nestedNPix= nestedSources[j]->NPix;
			double nestedX0= nestedSources[j]->X0;
			double nestedY0= nestedSources[j]->Y0;

			if(!IsGoodSource(nestedSources[j])) {
				INFO_LOG("Source no. "<<i<<": nested source no. "<<j<<" (name="<<nestedSourceName<<",id="<<nestedSourceId<<", n="<<nestedNPix<<"("<<nestedX0<<","<<nestedY0<<")) tagged as bad source, skipped!");
				nestedSources[j]->SetGoodSourceFlag(false);
			}
			if( IsPointLikeSource(nestedSources[j]) ){
				INFO_LOG("Source no. "<<i<<": nested source no. "<<j<<" (name="<<nestedSourceName<<",id="<<nestedSourceId<<", n="<<nestedNPix<<"("<<nestedX0<<","<<nestedY0<<")) tagged as a point-like source ...");
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

bool SFinder::IsGoodSource(Source* aSource){
	
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
		INFO_LOG("BoundingBox cut not passed (BoundingBoxMin="<<BoundingBoxMin<<"<"<<m_SourceMinBoundingBox<<")");
		return false;
	}

	//## Add other check here ...
	//...
	//...

	return true;

}//close IsGoodSource()

bool SFinder::IsPointLikeSource(Source* aSource){

	if(!aSource) return false;
	if(!aSource->HasParameters()) {
		WARN_LOG("No parameters are available for this source (did you compute them?)...point-like check cannot be performed!");
		return true;
	}

	std::string sourceName= aSource->GetName();
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

int SFinder::DeblendSources(std::vector<Source*>& sources){

	//Check given source list
	if(sources.empty()){
		WARN_LOG("Empty source list, nothing to be deblended!");
		return 0;
	}

	//## Loop over image sources and perform deblending stage for non-extended sources
	for(size_t i=0;i<sources.size();i++){

		//Skip non point-like sources
		int sourceType= sources[i]->Type;
		if(sourceType!=Source::ePointLike) continue;
	
		//Deblend mother source
		//if(sources[i]->Deblend(m_deblendCurvThr,deblendComponentMinNPix)<0) {//TO BE IMPLEMENTED!!!
		//	WARN_LOG("Failed to deblend source no. "<<i<<"!");
		//}

		//Deblend nested sources
		std::vector<Source*> nestedSources= sources[i]->GetNestedSources();	
		for(size_t j=0;j<nestedSources.size();j++){
			//if(nestedSources[j] && nestedSources[j]->Deblend(m_deblendCurvThr,deblendComponentMinNPix)<0){//TO BE IMPLEMENTED!!!
			//	WARN_LOG("Failed to deblend nested source no. "<<j<<" of source no. "<<i<<"...");
			//}
		}//end loop nested sources
	}//end loop sources
	
	return 0;

}//close DeblendSources()

/*
int SFinder::ReadImage(){

	//## Check file
	FileInfo info;
	bool match_extension= false;
	if(!SysUtils::CheckFile(m_InputFileName,info,match_extension,"")){
		ERROR_LOG("Invalid input file name specified (filename="<<m_InputFileName<<"), invalid file path?!");
		return -1;
	}
	m_InputFileExtension= info.extension;
	

	//=== ROOT reading ===
	if(m_InputFileExtension==".root"){// Read image from ROOT file
		TFile* inputFile = new TFile(m_InputFileName.c_str(),"READ");
		if(!inputFile || inputFile->IsZombie()){
			ERROR_LOG("Cannot open input file "<<m_InputFileName<<"!");
			return -1;
		}
		m_InputImg=  (Image*)inputFile->Get(m_InputImgName.c_str());
		if(!m_InputImg){
			ERROR_LOG("Cannot get image from input file "<<m_InputFileName<<"!");
			return -1;
		}
	}//close if

	//=== FITS reading ===
	else if(m_InputFileExtension==".fits"){// Read image from FITS file
		m_InputImg= new Image;
		int status= 0;
		if(m_ReadTile) status= m_InputImg->ReadFITS(m_InputFileName,m_TileMinX,m_TileMaxX,m_TileMinY,m_TileMaxY); 
		else status= m_InputImg->ReadFITS(m_InputFileName);

		if(status<0){
			ERROR_LOG("Failed to read image from input file "<<m_InputFileName<<"!");
			if(m_InputImg) m_InputImg->Delete();
			return -1;
		}
	}//close else if

	//== Invalid extension ==
	else{
		ERROR_LOG("Invalid file extension detected (ext="<<m_InputFileExtension<<")!");
		return -1;
	}
	m_InputImg->SetName("img");
	
	return 0;

}//close ReadImage()
*/

Image* SFinder::ReadImage(FileInfo& info,std::string filename,std::string imgname,long int ix_min,long int ix_max,long int iy_min,long int iy_max)
{
	//## Check file
	bool match_extension= false;
	if(!SysUtils::CheckFile(filename,info,match_extension,"")){
		ERROR_LOG("Invalid input file name specified (filename="<<filename<<"), invalid file path?!");
		return nullptr;
	}
	
	//## Read image from file
	bool readTile= (ix_min!=-1 && ix_max!=-1 && iy_min!=-1 && iy_max!=-1);
	Image* img= 0;
	
	//=== ROOT reading ===
	if(info.extension==".root"){// Read image from ROOT file
		TFile* inputFile = new TFile(filename.c_str(),"READ");
		if(!inputFile || inputFile->IsZombie()){
			ERROR_LOG("Cannot open input file "<<filename<<"!");
			return nullptr;
		}
		
		if(readTile){
			//Read full image
			Image* fullImg= (Image*)inputFile->Get(imgname.c_str());
			if(!fullImg){
				ERROR_LOG("Cannot get image "<<imgname<<" from input file "<<filename<<"!");
				return nullptr;
			}
			
			//Read tile
			img= fullImg->GetTile(ix_min,ix_max,iy_min,iy_max);	
			if(!img){
				ERROR_LOG("Failed to read image tile [xmin,xmax]=["<<ix_min<<","<<ix_max<<"], [ymin,ymax]=["<<iy_min<<","<<iy_max<<"]");
				delete fullImg;
				fullImg= 0;
				return nullptr;
			}				
		}//close if read tile
		else{
			img= (Image*)inputFile->Get(imgname.c_str());
			if(!img){
				ERROR_LOG("Cannot get image "<<imgname<<" from input file "<<filename<<"!");
				return nullptr;
			}	
		}

	}//close if

	//=== FITS reading ===
	else if(info.extension==".fits"){// Read image from FITS file
		img= new Image;

		int status= 0;
		if(readTile) {
			INFO_LOG("Reading image tile (file="<<filename<<", hdu="<<m_fitsHDUId<<", range[xmin,xmax]=["<<ix_min<<","<<ix_max<<"], [ymin,ymax]=["<<iy_min<<","<<iy_max<<"])");
			status= img->ReadFITS(filename,m_fitsHDUId,ix_min,ix_max,iy_min,iy_max); 
		}
		else {
			INFO_LOG("Reading image (file="<<filename<<", hdu="<<m_fitsHDUId<<")");
			status= img->ReadFITS(filename,m_fitsHDUId);
		}

		if(status<0){
			ERROR_LOG("Failed to read image from input file "<<filename<<"!");
			if(img) {
				delete img;
				img= 0;
			}
			return nullptr;
		}
	}//close else if FITS reading

	//== Invalid extension ==
	else{
		ERROR_LOG("Invalid file extension detected (ext="<<info.extension<<")!");
		return nullptr;
	}
	
	return img;

}//close ReadImage()




ImgBkgData* SFinder::ComputeStatsAndBkg(Image* img){

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
		INFO_LOG("Setting bkg boxes to ("<<boxSizeX<<","<<boxSizeY<<") pixels ...");	
	}

	double gridSizeX= m_GridSizeX*boxSizeX;
	double gridSizeY= m_GridSizeY*boxSizeY;
	INFO_LOG("Setting grid size to ("<<gridSizeX<<","<<gridSizeY<<") pixels ...");

	//## Compute Bkg
	ImgBkgData* bkgData= img->ComputeBkg (
		m_BkgEstimator,
		m_UseLocalBkg,boxSizeX,boxSizeY,gridSizeX,gridSizeY,
		m_Use2ndPassInLocalBkg,
		m_SkipOutliersInLocalBkg,m_SeedThr,m_MergeThr,m_NMinPix
	);

	if(!bkgData) {
		ERROR_LOG("Bkg computing failed!");
		return 0;
	}
		
	return bkgData;

}//close ComputeStatsAndBkg()


int SFinder::DrawSources(Image* image,std::vector<Source*>& sources){

	/*
	INFO_LOG("Drawing sources...");
	if(!image) return -1;

	bool useCurrentCanvas= false;
	bool drawFull= false;
	image->Plot(sources,useCurrentCanvas,drawFull,eRAINBOW,true);
	*/

	return 0;

}//close DrawSources()


int SFinder::SaveDS9RegionFile(){

	//## Open file
	FILE* fout= fopen(m_DS9CatalogFileName.c_str(),"w");

	//## Saving DS9 file region
	DEBUG_LOG("Saving DS9 region header...");
	fprintf(fout,"global color=red font=\"helvetica 12 normal\" edit=1 move=1 delete=1 include=1\n");
	fprintf(fout,"image\n");

	DEBUG_LOG("Saving "<<m_SourceCollection.size()<<" sources to file...");

	for(unsigned int k=0;k<m_SourceCollection.size();k++){
		DEBUG_LOG("Dumping DS9 region info for source no. "<<k<<" ...");
		std::string regionInfo= "";
		if(m_DS9RegionFormat==ePolygonRegion) regionInfo= m_SourceCollection[k]->GetDS9Region(true);
		else if(m_DS9RegionFormat==eEllipseRegion) regionInfo= m_SourceCollection[k]->GetDS9EllipseRegion(true);
		else {
			WARN_LOG("Invalid DS9RegionType given ("<<m_DS9RegionFormat<<")");
			return -1;
		}

		fprintf(fout,"%s\n",regionInfo.c_str());
	  	
	}//end loop sources
		
	DEBUG_LOG("Closing DS9 file region...");
	fclose(fout);

	return 0;

}//close SaveDS9RegionFile()


int SFinder::Save(){

	INFO_LOG("Storing results to file & catalog...");

	//Save DS9 regions?
	if(m_saveDS9Region && SaveDS9RegionFile()<0){
		WARN_LOG("Failed to save sources to DS9 region file!");
	}

	//Check ROOT output file
	if(!m_OutputFile) {
		WARN_LOG("Null ptr to output file, nothing will be saved in ROOT file!");
		return -1;
	}
	m_OutputFile->cd();

	//Save source tree?
	if(m_saveSources){
		DEBUG_LOG("Filling source ROOT TTree...");
		for(size_t k=0;k<m_SourceCollection.size();k++){
			m_Source= m_SourceCollection[k];
			m_SourceTree->Fill();
		}
		DEBUG_LOG("Writing tree to file...");
		m_SourceTree->Write();
	}
	
	//Save config?
	if(m_saveConfig){
		TTree* ConfigTree= ConfigParser::Instance().GetConfigTree("ConfigInfo");
		if(ConfigTree) ConfigTree->Write();
	}	

	//Save performance stats	
	if(m_PerfTree){
		m_PerfTree->Fill();
		m_PerfTree->Write();
	}

	//Save input image to file?
	if(m_saveInputMap && m_InputImg){
		m_InputImg->SetNameTitle("img","img");
		m_InputImg->Write();
	}
	
	//Save residual map?
	if(m_saveResidualMap && m_ResidualImg){
		m_ResidualImg->SetNameTitle("img_residual","img_residual");
		m_ResidualImg->Write();
	}
	
	//Save significance map?
	if(m_saveSignificanceMap && m_SignificanceMap){
		m_SignificanceMap->SetNameTitle("img_significance","img_significance");
		m_SignificanceMap->Write();
	}

	//Save bkg & noise maps
	if(m_saveBkgMap && m_BkgData && m_BkgData->BkgMap){
		(m_BkgData->BkgMap)->SetNameTitle("img_bkg","img_bkg");
		(m_BkgData->BkgMap)->Write();
	}
	if(m_saveNoiseMap && m_BkgData && m_BkgData->NoiseMap){
		(m_BkgData->NoiseMap)->SetNameTitle("img_rms","img_rms");
		(m_BkgData->NoiseMap)->Write();
	}

	//Save saliency map
	if(m_saveSaliencyMap && m_SaliencyImg){
		m_SaliencyImg->SetNameTitle("img_saliency","img_saliency");
		m_SaliencyImg->Write();
	}

	//Save Laplacian
	if(m_saveCurvatureMap && m_LaplImg){
		m_LaplImg->SetNameTitle("img_lapl","img_lapl");
		m_LaplImg->Write();
	}

	//Save Edgeness
	if(m_saveEdgenessMap && m_EdgeImg){
		m_EdgeImg->SetNameTitle("img_edge","img_edge");
		m_EdgeImg->Write();
	}

	//Save segmented map
	if(m_saveSegmentedMap && m_SegmImg){
		m_SegmImg->SetNameTitle("img_segm","img_segm");
		m_SegmImg->Write();
	}

	//## Close ROOT output file
	DEBUG_LOG("Closing output file...");
	if(m_OutputFile && m_OutputFile->IsOpen()) m_OutputFile->Close();

	INFO_LOG("End save to file");

	return 0;

}//close Save()


void SFinder::PrintPerformanceStats()
{
	INFO_LOG("===========================");
	INFO_LOG("===   PERFORMANCE INFO  ===");
	INFO_LOG("===========================");
	INFO_LOG("tot (ms)= "<<totTime);
	INFO_LOG("init (ms)= "<<initTime<<" ["<<initTime/totTime*100.<<"%]");
	INFO_LOG("read image (ms)= "<<readImageTime<<" ["<<readImageTime/totTime*100.<<"%]");
	INFO_LOG("source finding (ms)= "<<compactSourceTime<<" ["<<compactSourceTime/totTime*100.<<"%]");
	INFO_LOG("source selection (ms)= "<<sourceSelectionTime<<" ["<<sourceSelectionTime/totTime*100.<<"%]");
	INFO_LOG("img residual (ms)= "<<imgResidualTime<<" ["<<imgResidualTime/totTime*100.<<"%]");
	INFO_LOG("ext source finding (ms)= "<<extendedSourceTime<<" ["<<extendedSourceTime/totTime*100.<<"%]");
	INFO_LOG("source deblending (ms)= "<<sourceDeblendTime<<" ["<<sourceDeblendTime/totTime*100.<<"%]");
	INFO_LOG("save (ms)= "<<saveTime<<" ["<<saveTime/totTime*100.<<"%]");
	INFO_LOG("===========================");

}//close PrintPerformanceStats()


int SFinder::PrepareWorkerTasks()
{
	//## Generate a uuid for this job
	std::string jobId= CodeUtils::GenerateUUID();
	DEBUG_LOG("[PROC "<<m_procId<<"] - Generated jobId: "<<jobId);
	
	//## Get input image size
	FileInfo info;
	bool match_extension= false;
	if(!SysUtils::CheckFile(m_InputFileName,info,match_extension,"")){
		ERROR_LOG("[PROC "<<m_procId<<"] - Invalid input file name specified (filename="<<m_InputFileName<<"), invalid file path?!");
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
				ERROR_LOG("[PROC "<<m_procId<<"] - Failed to open input file image "<<m_InputFileName<<" and get image size!");
				return -1;
			}
			Image* inputImg= (Image*)inputFile->Get(m_InputImgName.c_str());
			if(!inputImg) {
				ERROR_LOG("[PROC "<<m_procId<<"] - Failed to open input file image "<<m_InputFileName<<" and get image size!");
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
				ERROR_LOG("[PROC "<<m_procId<<"] - Failed to open input file image "<<m_InputFileName<<" and get image size!");
				return -1;
			}
			INFO_LOG("[PROC "<<m_procId<<"] - Input image opened with success: size="<<Nx<<"x"<<Ny);
		}
	}//close else if		
	else{
		ERROR_LOG("[PROC "<<m_procId<<"] - Invalid/unsupported file extension ("<<info.extension<<") detected!");
		return -1;
	}

	INFO_LOG("[PROC "<<m_procId<<"] - Image size: "<<Nx<<"x"<<Ny);

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

	INFO_LOG("[PROC "<<m_procId<<"] - Computing tile partition: tileSize("<<m_TileSizeX<<","<<m_TileSizeY<<"), tileOverlap("<<tileOverlapX<<","<<tileOverlapY<<")");
	if(MathUtils::Compute2DGrid(ix_min,ix_max,iy_min,iy_max,Nx,Ny,m_TileSizeX,m_TileSizeY,m_TileStepSizeX,m_TileStepSizeY)<0){
		WARN_LOG("[PROC "<<m_procId<<"] - Failed to compute a 2D partition from input image!");
		return -1;
	}
	int nExpectedTasks= ix_min.size()*iy_min.size();
	INFO_LOG("[PROC "<<m_procId<<"] - #"<<nExpectedTasks<<" expected number of tasks ("<<ix_min.size()<<"x"<<iy_min.size()<<")");

	//## Compute worker tasks (check max number of tasks per worker)
	INFO_LOG("[PROC "<<m_procId<<"] - Computing worker task list...");
	TaskData* aTaskData= 0;
	long int workerCounter= 0;

	for(size_t j=0;j<iy_min.size();j++){
		for(size_t i=0;i<ix_min.size();i++){
			//Assign worker
			INFO_LOG("[PROC "<<m_procId<<"] - Assign task ("<<i<<","<<j<<") to worker no. "<<workerCounter<<"...");
				
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

			if(workerCounter>=m_nProc-1) workerCounter= 0;
			else workerCounter++;
		}//end loop x
	}//end loop y

	
	//Fill neighbor task list
	std::vector<int> workerIds;

	for(size_t i=0;i<m_taskDataPerWorkers.size();i++){
		if(m_taskDataPerWorkers[i].size()==0) continue;//no tasks present

		//Add only processors with tasks
		workerIds.push_back(i);
		
		//Loop over tasks present in this worker
		int nTasksInWorker= static_cast<int>(m_taskDataPerWorkers[i].size()); 
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
 
				std::stringstream ss;	
				ss<<"[PROC "<<m_procId<<"] - Worker no. "<<i<<", Task "<<j<<"["<<ix_min<<","<<ix_max<<"] ["<<iy_min<<","<<iy_max<<"], NextTask "<<k<<"["<<next_ix_min<<","<<next_ix_max<<"] ["<<next_iy_min<<","<<next_iy_max<<"] ==> isAdjacentInX? "<<isAdjacentInX<<", isAdjacentInY? "<<isAdjacentInY;
				if(m_procId==MASTER_ID) INFO_LOG(ss.str());

				if(isAdjacent || isOverlapping) {
					(m_taskDataPerWorkers[i][j]->neighborTaskId).push_back(k);
					(m_taskDataPerWorkers[i][k]->neighborTaskId).push_back(j);
					(m_taskDataPerWorkers[i][j]->neighborWorkerId).push_back(i);
					(m_taskDataPerWorkers[i][k]->neighborWorkerId).push_back(i);
				}

			}//end loop next task in worker


			//Find neighbors across workers
			for(size_t s=i+1;s<m_taskDataPerWorkers.size();s++){
				for(size_t t=0;t<m_taskDataPerWorkers[s].size();t++){
					
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

					std::stringstream ss;	
					ss<<"[PROC "<<m_procId<<"] - Worker no. "<<i<<", Task "<<j<<", NextWorker no. "<<s<<", NextTask "<<t<<"["<<next_ix_min<<","<<next_ix_max<<"] ["<<next_iy_min<<","<<next_iy_max<<"] ==> isAdjacentInX? "<<isAdjacentInX<<", isAdjacentInY? "<<isAdjacentInY;
					if(m_procId==MASTER_ID) INFO_LOG(ss.str());

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

	int nWorkers= static_cast<int>(workerIds.size());
	std::stringstream ss;
	ss<<"[PROC "<<m_procId<<"] - # "<<nWorkers<<" workers {";
	for(int i=0;i<nWorkers;i++){
		ss<<workerIds[i]<<",";
	}
	ss<<"}";
	INFO_LOG(ss.str());
	

	
	//## Create a worker group (if MPI run is performed)
	m_workerRanks= -1;
	m_nWorkers= 1;

	#ifdef MPI_ENABLED
	if(m_mpiEnabled){

		//Get main processor group
		MPI_Comm_group(MPI_COMM_WORLD, &m_WorldGroup);
	
		// Construct a group containing all of the workers (proc with tasks assigned)
		MPI_Group_incl(m_WorldGroup, nWorkers, workerIds.data() , &m_WorkerGroup);

		// Create a new communicator based on the group
		int commTag= 10;
		MPI_Comm_create_group(MPI_COMM_WORLD, m_WorkerGroup, commTag, &m_WorkerComm);

		m_workerRanks = -1;
		m_nWorkers = -1;
	
		// If this rank isn't in the new communicator, it will be
		// MPI_COMM_NULL. Using MPI_COMM_NULL for MPI_Comm_rank or
		// MPI_Comm_size is erroneous
		if (m_WorkerComm!=MPI_COMM_NULL) {
    	MPI_Comm_rank(m_WorkerComm, &m_workerRanks);
    	MPI_Comm_size(m_WorkerComm, &m_nWorkers);
		}
		else {
			WARN_LOG("[PROC "<<m_procId<<"] - Worker MPI communicator is null (this processor has no tasks and was not inserted in the worker group)!");
		}

	}//close if	
	#endif
	
	INFO_LOG("[PROC "<<m_procId<<"] - WORLD RANK/SIZE: "<<m_procId<<"/"<<m_nProc<<" WORKER RANK/SIZE: "<<m_workerRanks<<"/"<<m_nWorkers);


	//Print
	if(m_procId==MASTER_ID){
		for(size_t i=0;i<m_taskDataPerWorkers.size();i++){
			if(m_taskDataPerWorkers[i].size()==0) continue;//no tasks present

			for(size_t j=0;j<m_taskDataPerWorkers[i].size();j++){
				std::stringstream ss;	
				ss<<"[PROC "<<m_procId<<"] - Worker no. "<<i<<", ";

				long int ix_min= m_taskDataPerWorkers[i][j]->ix_min;
				long int ix_max= m_taskDataPerWorkers[i][j]->ix_max;
				long int iy_min= m_taskDataPerWorkers[i][j]->iy_min;
				long int iy_max= m_taskDataPerWorkers[i][j]->iy_max;
			
				ss<<"Task no. "<<j<<"["<<ix_min<<","<<ix_max<<"] ["<<iy_min<<","<<iy_max<<"] neighbors{";
			
				for(size_t k=0;k<m_taskDataPerWorkers[i][j]->neighborTaskId.size();k++){
					long int neighborWorkerId= m_taskDataPerWorkers[i][j]->neighborWorkerId[k];
					long int neighborTaskId= m_taskDataPerWorkers[i][j]->neighborTaskId[k];
					long int next_ix_min= m_taskDataPerWorkers[neighborWorkerId][neighborTaskId]->ix_min;
					long int next_ix_max= m_taskDataPerWorkers[neighborWorkerId][neighborTaskId]->ix_max;
					long int next_iy_min= m_taskDataPerWorkers[neighborWorkerId][neighborTaskId]->iy_min;
					long int next_iy_max= m_taskDataPerWorkers[neighborWorkerId][neighborTaskId]->iy_max;

					ss<<"("<<neighborWorkerId<<","<<neighborTaskId<<") ["<<next_ix_min<<","<<next_ix_max<<"] ["<<next_iy_min<<","<<next_iy_max<<"], ";
				}	
				ss<<"}";
				INFO_LOG(ss.str());	
			}//end loop tasks
		}//end loop workers
	}//close if MASTER


	bool hasTooManyTasks= false;
	for(size_t i=0;i<m_taskDataPerWorkers.size();i++){
		long int nTasksPerWorker= static_cast<long int>(m_taskDataPerWorkers[i].size());
		if(nTasksPerWorker>MAX_NTASKS_PER_WORKER){
			hasTooManyTasks= true;
			break;
		}
	}

	if(hasTooManyTasks){
		WARN_LOG("[PROC "<<m_procId<<"] - Too many tasks per worker (thr="<<MAX_NTASKS_PER_WORKER<<")");
		return -1;
	}

	return 0;

}//close PrepareWorkerTasks()


#ifdef MPI_ENABLED
int SFinder::GatherTaskDataFromWorkers()
{
	//## Put a barrier and collect all sources from workers in the master processor
	MPI_Barrier(MPI_COMM_WORLD);

	//## Sum and average all the elapsed timers across workers 
	double initTime_sum;
	double readImageTime_sum;
	double compactSourceTime_sum;
	double sourceSelectionTime_sum;
	double imgResidualTime_sum;
	double extendedSourceTime_sum;
	MPI_Reduce(&initTime, &initTime_sum, 1, MPI_DOUBLE, MPI_SUM, MASTER_ID, m_WorkerComm);
	MPI_Reduce(&readImageTime, &readImageTime_sum, 1, MPI_DOUBLE, MPI_SUM, MASTER_ID, m_WorkerComm);
	MPI_Reduce(&compactSourceTime, &compactSourceTime_sum, 1, MPI_DOUBLE, MPI_SUM, MASTER_ID, m_WorkerComm);
	MPI_Reduce(&imgResidualTime, &imgResidualTime_sum, 1, MPI_DOUBLE, MPI_SUM, MASTER_ID, m_WorkerComm);
	MPI_Reduce(&sourceSelectionTime, &sourceSelectionTime_sum, 1, MPI_DOUBLE, MPI_SUM, MASTER_ID, m_WorkerComm);
	MPI_Reduce(&extendedSourceTime, &extendedSourceTime_sum, 1, MPI_DOUBLE, MPI_SUM, MASTER_ID, m_WorkerComm);
	INFO_LOG("[PROC "<<m_procId<<"] - cpu times (ms): {init="<<initTime<<", read="<<readImageTime<<", sourcefind="<<compactSourceTime<<", residual="<<imgResidualTime<<", sourcesel="<<sourceSelectionTime<<", extsourcefind="<<extendedSourceTime);

	if (m_procId == MASTER_ID) {
		initTime= initTime_sum;
		readImageTime= readImageTime_sum;
		compactSourceTime= compactSourceTime_sum;
		sourceSelectionTime= sourceSelectionTime_sum;
		imgResidualTime= imgResidualTime_sum;
		extendedSourceTime= extendedSourceTime_sum;
		INFO_LOG("[PROC "<<m_procId<<"] - cumulative cpu times (ms): {init="<<initTime<<", read="<<readImageTime<<", sourcefind="<<compactSourceTime<<", residual="<<imgResidualTime<<", sourcesel="<<sourceSelectionTime<<", extsourcefind="<<extendedSourceTime);
	}

	//## Merge all sources found by workers in a unique collection
	int MSG_TAG= 1;
	if (m_procId == MASTER_ID) {//Receive data from the workers

		for (int i=1; i<m_nProc; i++) {
			//Check if this processor has tasks assigned, otherwise skip!
			if(m_taskDataPerWorkers[i].size()==0){
				INFO_LOG("[PROC "<<m_procId<<"] - No tasks assigned to worker no. "<<i<<", nothing to be collected, skip to next worker...");
				continue;
			}
  
			//## Probe for an incoming message from process zero
			INFO_LOG("[PROC "<<m_procId<<"] - Probing for message from proc "<<i);
    	MPI_Status status;

			if(MPI_Probe(i, MSG_TAG, MPI_COMM_WORLD, &status)==MPI_SUCCESS){
				INFO_LOG("[PROC "<<m_procId<<"] - a message has been found with the probe, with tag " << status.MPI_TAG << ", source " << status.MPI_SOURCE);

    		//## When probe returns, the status object has the size and other
    		//## attributes of the incoming message. Get the message size
				INFO_LOG("[PROC "<<m_procId<<"] - Getting size of message... ");
    	
				int rcvMsgSize= 0;
    		MPI_Get_count(&status, MPI_CHAR, &rcvMsgSize);

				//## Allocate a buffer to hold the incoming numbers
				INFO_LOG("[PROC "<<m_procId<<"] - Allocating a message of size "<<rcvMsgSize);
				if(rcvMsgSize<=0){
					ERROR_LOG("[PROC "<<m_procId<<"] - rcvMsg size is negative/null!");
					continue;
				}
				char* recvBuffer= (char*)malloc(rcvMsgSize);
    		
    		//## Now receive the message with the allocated buffer
    		//MPI_Recv(recvBuffer, rcvMsgSize, MPI_CHAR, MPI_ANY_SOURCE, MSG_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(recvBuffer, rcvMsgSize, MPI_CHAR, i, MSG_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    		INFO_LOG("[PROC "<<m_procId<<"] - Received a message of size "<<rcvMsgSize<<") from process "<<i);

				//## Update task data with received worker data	
				bool isTaskCollectionPreAllocated= true;
				if(Serializer::CharArrayToTaskDataCollection(m_taskDataPerWorkers[i],recvBuffer,rcvMsgSize,isTaskCollectionPreAllocated)<0 ){
					ERROR_LOG("[PROC "<<m_procId<<"] - Failed to decode recv message into task data list!");
    			if(recvBuffer) free(recvBuffer);
				}

				//## Free received buffer
    		if(recvBuffer) free(recvBuffer);
			}//close if
			else{
				ERROR_LOG("[PROC "<<m_procId<<"] - Message probing failed!");
				return -1;
				//continue;
			}
		}//end loop workers
	}//close if

	else {//Send data to master
		//## First encode taskData in protobuf
		long int msg_size= 0;
		char* msg= Serializer::TaskDataCollectionToCharArray(msg_size,m_taskDataPerWorkers[m_procId]);
		if(!msg){
			ERROR_LOG("[PROC "<<m_procId<<"] - Failed to encode task data to protobuf!");
			return -1;
		}

		//## Send buffer to master processor	
		INFO_LOG("[PROC "<<m_procId<<"] - Sending msg: "<<msg<<" (size="<<msg_size<<")");
		MPI_Send((void*)(msg),msg_size, MPI_CHAR, MASTER_ID, MSG_TAG, MPI_COMM_WORLD);

		//## Free buffer
		free(msg);
	}//close else

	MPI_Barrier(MPI_COMM_WORLD);

	return 0;

}//close GatherTaskDataFromWorkers()
#endif


int SFinder::MergeTaskData()
{
	//## Update sources in list
	//## NB: Push to collection only sources NON at edge
	if (m_procId == 0) {

		//Print task data
		INFO_LOG("[PROC "<<m_procId<<"] - Printing task data...");
		for(size_t i=0;i<m_taskDataPerWorkers.size();i++){
			if(m_taskDataPerWorkers[i].size()==0) continue;//no tasks present
			
			std::stringstream ss;
			ss<<"[PROC "<<m_procId<<"] - Worker no. "<<i<<", ";
			for(size_t j=0;j<m_taskDataPerWorkers[i].size();j++){
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

		for(size_t i=0;i<m_taskDataPerWorkers.size();i++){
			for(size_t j=0;j<m_taskDataPerWorkers[i].size();j++){
				m_CompactSources.insert(m_CompactSources.end(),(m_taskDataPerWorkers[i][j]->sources).begin(),(m_taskDataPerWorkers[i][j]->sources).end());
				m_CompactSources.insert(m_CompactSources.end(),(m_taskDataPerWorkers[i][j]->sources_edge).begin(),(m_taskDataPerWorkers[i][j]->sources_edge).end());

				m_ExtendedSources.insert(m_ExtendedSources.end(),(m_taskDataPerWorkers[i][j]->ext_sources).begin(),(m_taskDataPerWorkers[i][j]->ext_sources).end());
				m_ExtendedSources.insert(m_ExtendedSources.end(),(m_taskDataPerWorkers[i][j]->ext_sources_edge).begin(),(m_taskDataPerWorkers[i][j]->ext_sources_edge).end());
			}//end loop tasks
		}//end loop workers
		m_SourceCollection.insert(m_SourceCollection.end(),m_CompactSources.begin(),m_CompactSources.end());
		m_SourceCollection.insert(m_SourceCollection.end(),m_ExtendedSources.begin(),m_ExtendedSources.end());

		INFO_LOG("[PROC "<<m_procId<<"] - #"<<m_SourceCollection.size()<<" sources found in total (#"<<m_CompactSources.size()<<" compact, #"<<m_ExtendedSources.size()<<" extended) ...");
		
	}//close if

	return 0;

}//close MergeTaskData()

int SFinder::MergeSourcesAtEdge()
{
	//## Executed only by master processor
	if(m_procId != MASTER_ID) return 0;

	//## Loop over edge sources and merge them
	//## TBD: Merge all sources regardless of their tag (compact or extended)
	struct MergedSourceInfo {
		long int source_index;
		long int worker_index;
		long int task_index;
		//std::vector<MergedSourceInfo> merged_sources;
		MergedSourceInfo(long int sindex,long int windex,long int tindex)
			: source_index(sindex), worker_index(windex), task_index(tindex)
		{}
		bool operator==(const MergedSourceInfo& obj) const {
			bool areEqual= ( 
				(obj.source_index==source_index) && 	
				(obj.worker_index==worker_index) &&
				(obj.task_index==task_index)
			);
    	return areEqual;
    }
	};//close MergedSourceInfo

	//## Fill list of edge sources to be merged and fill corresponding Graph
	std::vector<MergedSourceInfo> sourcesToBeMerged;
	Graph mergedSourceGraph;
	for(size_t i=0;i<m_taskDataPerWorkers.size();i++){
		if(m_taskDataPerWorkers[i].size()==0) continue;//no tasks present
		for(size_t j=0;j<m_taskDataPerWorkers[i].size();j++){
			int nEdgeSources= static_cast<int>((m_taskDataPerWorkers[i][j]->sources_edge).size());
			for(int k=0;k<nEdgeSources;k++){
				MergedSourceInfo mergedSourceInfo= MergedSourceInfo(k,i,j);
				sourcesToBeMerged.push_back(MergedSourceInfo(k,i,j));
				mergedSourceGraph.AddVertex();
			}
		}
	}

	//## Find adjacent edge sources
	//std::vector<MergedSourceInfo> sourcesAlreadyMerged;
	int itemPos= -1;

	for(size_t i=0;i<sourcesToBeMerged.size()-1;i++){	
		long int sindex= sourcesToBeMerged[i].source_index;
		long int windex= sourcesToBeMerged[i].worker_index;
		long int tindex= sourcesToBeMerged[i].task_index; 
		Source* source= (m_taskDataPerWorkers[windex][tindex]->sources_edge)[sindex];

		//Loop neighbors
		for(size_t j=i+1;j<sourcesToBeMerged.size();j++){	
			long int sindex_neighbor= sourcesToBeMerged[j].source_index;
			long int windex_neighbor= sourcesToBeMerged[j].worker_index;
			long int tindex_neighbor= sourcesToBeMerged[j].task_index; 
			Source* source_neighbor= (m_taskDataPerWorkers[windex_neighbor][tindex_neighbor]->sources_edge)[sindex_neighbor];

			//Check if main worker tile is physically neighbor to this
			//If not they cannot be adjacent, so skip the following check
			int itemPos= -1;
			
			if(!CodeUtils::FindItem(m_taskDataPerWorkers[windex][tindex]->neighborWorkerId,windex_neighbor,itemPos)){
				INFO_LOG("Worker id "<<windex_neighbor<<" is not physically neighbor to worker "<<windex<<", so skip the source adjacency check and go to next...");
				continue;
			}

			//Check is sources are adjacent
			bool areAdjacentSources= source->IsAdjacentSource(source_neighbor);
			if(!areAdjacentSources) continue;

			//If they are adjacent add linking in graph
			INFO_LOG("Sources (i,j)=("<<i<<" {"<<sindex<<","<<windex<<","<<tindex<<"} , "<<j<<" {"<<sindex_neighbor<<","<<windex_neighbor<<","<<tindex_neighbor<<"}) are adjacent and selected for merging...");
			mergedSourceGraph.AddEdge(i,j);

		}//end loop sources

	}//end loop sources


	//## Find all connected components in graph corresponding to 
	//## edge sources to be merged
	std::vector<std::vector<int>> connected_source_indexes;
	mergedSourceGraph.GetConnectedComponents(connected_source_indexes);
	INFO_LOG("#"<<connected_source_indexes.size()<<"/"<<sourcesToBeMerged.size()<<" edge sources will be left after merging...");
		
	//## Now merge the sources
	std::vector<MergedSourceInfo> sourcesToBeRemoved;
	bool copyPixels= false;//do not create memory for new pixels
	bool checkIfAdjacent= false;//already done before
	bool computeStatPars= false;//do not compute stats& pars at each merging
	bool computeMorphPars= false;
	bool computeRobustStats= true;
	bool forceRecomputing= false;//no need to re-compute moments (already updated in AddPixel())

	for(size_t i=0;i<connected_source_indexes.size();i++){
		if(connected_source_indexes[i].empty()) continue;

		//Get source id=0 of this component
		int index= connected_source_indexes[i][0];
		long int sindex= sourcesToBeMerged[index].source_index;
		long int windex= sourcesToBeMerged[index].worker_index;
		long int tindex= sourcesToBeMerged[index].task_index; 
		Source* source= (m_taskDataPerWorkers[windex][tindex]->sources_edge)[sindex];

		//Merge other sources in the group if any 
		int nMerged= 0;
		for(size_t j=1;j<connected_source_indexes[i].size();j++){
			int index_merged= connected_source_indexes[i][j];
			long int sindex_merged= sourcesToBeMerged[index_merged].source_index;
			long int windex_merged= sourcesToBeMerged[index_merged].worker_index;
			long int tindex_merged= sourcesToBeMerged[index_merged].task_index; 
			Source* source_merged= (m_taskDataPerWorkers[windex_merged][tindex_merged]->sources_edge)[sindex_merged];
				
			int status= source->MergeSource(source_merged,copyPixels,checkIfAdjacent,computeStatPars,computeMorphPars);
			if(status<0){
				WARN_LOG("Failed to merge sources (i,j)=("<<index<<" {"<<sindex<<","<<windex<<","<<tindex<<"} , "<<index_merged<<" {"<<sindex_merged<<","<<windex_merged<<","<<tindex_merged<<"}), skip to next...");
				continue;
			}

			//Add this source to the list of edge sources to be removed
			sourcesToBeRemoved.push_back(index_merged);

		}//end loop of sources to be merged in this component

		//If at least one was merged recompute stats & pars of merged source
		if(nMerged>0) {
			INFO_LOG("Recomputing stats & moments of merged source in merge group "<<i<<" after #"<<nMerged<<" merged source...");
			if(source->ComputeStats(computeRobustStats,forceRecomputing)<0){
				WARN_LOG("Failed to compute stats for merged source in merge group "<<i<<"...");
				continue;
			}
			if(source->ComputeMorphologyParams()<0){
				WARN_LOG("Failed to compute morph pars for merged source in merge group "<<i<<"...");
				continue;
			}
		}//close if

	}//end loop number of components

	//Remove merged source from task edge source collection
	INFO_LOG("Removing #"<<sourcesToBeRemoved.size()<<" edge sources that has been merged previously from task data collection...");
	for(size_t i=0;i<m_taskDataPerWorkers.size();i++){
		if(m_taskDataPerWorkers[i].size()==0) continue;//no tasks present

		//Loop over tasks per worker
		for(size_t j=0;j<m_taskDataPerWorkers[i].size();j++){
	
			//Loop over sources found
			std::vector<Source*> sources_merged_at_edges;
			int nSources= static_cast<int>((m_taskDataPerWorkers[i][j]->sources).size());

			for(int k=0;k<nSources;k++){
				(m_taskDataPerWorkers[i][j]->sources)[k]
			}//end 

	for(size_t i=0;i<sourcesToBeRemoved.size();i++){
		int index= sourcesToBeRemoved[i];
		long int sindex= sourcesToBeMerged[index].source_index;
		long int windex= sourcesToBeMerged[index].worker_index;
		long int tindex= sourcesToBeMerged[index].task_index;

		//Delete item
		CodeUtils::DeleteItem();
	}//end loop sources
	
	return 0;

}//close MergeSourcesAtEdge()


int SFinder::FindSourcesAtEdge()
{	
	//## Find if sources (both compact and extended) are at tile edge
	//## Those found at the edge are removed from the list and added to the edge list for further processing
	double xmin_s, xmax_s, ymin_s, ymax_s;
	
	//Loop over workers
	for(size_t i=0;i<m_taskDataPerWorkers.size();i++){
		if(m_taskDataPerWorkers[i].size()==0) continue;//no tasks present

		//Loop over tasks per worker
		for(size_t j=0;j<m_taskDataPerWorkers[i].size();j++){
	
			//Loop over sources found
			std::vector<Source*> sources_not_at_edges;
			int nSources= static_cast<int>((m_taskDataPerWorkers[i][j]->sources).size());

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
				INFO_LOG("[PROC "<<m_procId<<"] - workerId="<<m_taskDataPerWorkers[i][j]->workerId<<", check if compact source no. "<<k<<"(x["<<xmin_s<<","<<xmax_s<<"] y["<<ymin_s<<","<<ymax_s<<"]) is at edge of its tile (x["<<xmin_tile<<","<<xmax_tile<<"] y["<<ymin_tile<<","<<ymax_tile<<"]), isAtTileEdgeX="<<isAtTileEdgeX<<", isAtTileEdgeY="<<isAtTileEdgeY<<", isAtTileEdge="<<isAtTileEdge);

				INFO_LOG("[PROC "<<m_procId<<"] - workerId="<<m_taskDataPerWorkers[i][j]->workerId<<", check if compact source no. "<<k<<"(x["<<xmin_s<<","<<xmax_s<<"] y["<<ymin_s<<","<<ymax_s<<"]) is inside neighbour tile...");
		
				//Check if source is inside neighbour tile, e.g. is in overlapping area
				//bool isAtEdge= false;
				bool isInOverlapArea= false;
				for(size_t l=0;l<(m_taskDataPerWorkers[i][j]->neighborWorkerId).size();l++){	
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

					INFO_LOG("[PROC "<<m_procId<<"] - neighborWorkerId="<<neighborWorkerId<<", neighborTaskId="<<neighborTaskId<<", check if inside neighbor tile no. "<<j<<"(x["<<xmin<<","<<xmax<<"] y["<<ymin<<","<<ymax<<"]), isOverlappingX="<<isOverlappingX<<", isOverlappingY="<<isOverlappingY);
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

			int nEdgeSources= static_cast<int>((m_taskDataPerWorkers[i][j]->sources_edge).size());
			INFO_LOG("[PROC "<<m_procId<<"] - #"<<nEdgeSources<<"/"<<nSources<<" sources are found at tile edges...");
	
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

}//close namespace
