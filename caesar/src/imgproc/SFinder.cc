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
#include <Logger.h>
#include <Consts.h>

#include <SLIC.h>
#include <SLICUtils.h>
#include <SLICSegmenter.h>
#include <ChanVeseSegmenter.h>

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

using namespace std;

ClassImp(Caesar::SFinder)

namespace Caesar {

SFinder::SFinder() 
{
	
}//close costructor


SFinder::~SFinder(){
	
	//if(m_SourceTree) m_SourceTree->Delete();
	if(m_OutputFile) m_OutputFile->Close();
	//if(m_DS9CatalogFilePtr) fclose(m_DS9CatalogFilePtr);
	//if(m_Application) m_Application->Delete();

	//Delete images
	if(m_EdgeImg){
		DEBUG_LOG("Deleting edge image...");
		delete m_EdgeImg;
		m_EdgeImg= 0;
	}
	if(m_LaplImg){
		DEBUG_LOG("Deleting Laplacian image...");
		delete m_LaplImg;
		m_LaplImg= 0;
	}

}//close destructor


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
	//m_DS9CatalogFilePtr= 0;
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
	if(m_saveToFile){
		m_OutputFile= new TFile(m_OutputFileName.c_str(),"RECREATE");	
		m_OutputFile->cd();
	
		//Init source tree
		if(m_saveSources){
			m_Source= 0;
			if(!m_SourceTree) m_SourceTree= new TTree("SourceInfo","SourceInfo");
			m_SourceTree->Branch("Source",&m_Source);
			m_SourceCollection.clear();
		}

		//Init DS9 catalog
		//if(!m_DS9CatalogFilePtr) m_DS9CatalogFilePtr= fopen(m_DS9CatalogFileName.c_str(),"w");

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

	return 0;

}//close Init()

int SFinder::Configure(){

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
	if(ReadImage()<0){
		ERROR_LOG("Reading of input image failed!");
		return -1;
	}
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


int SFinder::FindSources(std::vector<Source*>& sources,Image* inputImg,double seedThr,double mergeThr){

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
	int status= inputImg->FindCompactSource(
		sources,
		significanceMap,bkgData,
		seedThr,mergeThr,m_NMinPix,m_SearchNegativeExcess,m_MergeBelowSeed,
		m_SearchNestedSources,m_NestedBlobThrFactor
	);

	if(status<0) {
		ERROR_LOG("Source finding failed!");	
		if(bkgData){
			delete bkgData;
			bkgData= 0;	
		}
		if(significanceMap) significanceMap->Delete();
		return -1;
	}
	int nSources= (int)sources.size();
	INFO_LOG("#"<<nSources<<" sources detected in input image...");	
	
	if(bkgData){
		delete bkgData;
		bkgData= 0;
	}
	if(significanceMap) significanceMap->Delete();
		
	return 0;

}//close FindSources()


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
	else if(m_ExtendedSearchMethod==eChanVese){
		status= FindExtendedSources_ChanVese(inputImg);
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
	
	//## Add detected sources to the list	
	m_SourceCollection.insert(m_SourceCollection.end(),sources.begin(),sources.end());
	m_ExtendedSources.insert(m_ExtendedSources.end(),sources.begin(),sources.end());

	INFO_LOG("#"<<nSelSources<<" extended sources to the list...");


	return 0;

}//close FindExtendedSources_SalThr()


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
	SLICData* slicData_init= SLICUtils::SPGenerator(inputImg,m_spSize,m_spBeta,m_spMinArea,normalizeImage,m_spUseLogContrast,m_LaplImg,m_EdgeImg);
	if(!slicData_init){
		ERROR_LOG("Failed to compute the initial superpixel partition, cannot perform extended source finding!");	
		return -1;
	}

	//## Tag the superpixel partition
	if(SLICUtils::TagRegions(slicData_init->regions,bkgMarkerImg,signalMarkerImg)<0){
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
	m_SegmImg= SLICUtils::GetSegmentedImage(inputImg,slicData_segm.regions,Region::eSignalTag,normalizeSegmImg,binarizeSegmImg);
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
			m_cvTimeStepPar, m_cvWindowSizePar, m_cvLambda1Par, m_cvLambda2Par, m_cvMuPar, m_cvNuPar, m_cvPPar
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


int SFinder::FindExtendedSources_ChanVese(Image* inputImg){

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

	//==========================================
	//==    RUN CHAN-VESE SEGMENTATION
	//==========================================
	//## Perform segmentation
	/*
	std::vector<Source*> sources;
	int status= inputImg->FindExtendedSource_CV(
		sources,signalMarkerImg,
		m_ResidualBkgData, m_NMinPix, m_SearchNegativeExcess,
		m_cvTimeStepPar,m_cvWindowSizePar,m_cvLambda1Par,m_cvLambda2Par,m_cvMuPar,m_cvNuPar,m_cvPPar
	);

	if(status<0){
		ERROR_LOG("ChanVese Segmentation failed!");
		return -1;
	}
	*/

	//## Compute segmented image
	bool returnContourImg= false;
	m_SegmImg= ChanVeseSegmenter::FindSegmentation (
		inputImg, signalMarkerImg, returnContourImg,
		m_cvTimeStepPar,m_cvWindowSizePar,m_cvLambda1Par,m_cvLambda2Par,m_cvMuPar,m_cvNuPar,m_cvPPar
	);
	if(!m_SegmImg){
		ERROR_LOG("Failed to compute ChanVese image segmentation!");
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

}//close FindExtendedSources_ChanVese()

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
		for(size_t j=0;j<nestedSources.size();j++){
			std::string nestedSourceName= nestedSources[j]->GetName();
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
		DEBUG_LOG("BoundingBox cut not passed (BoundingBoxMin="<<BoundingBoxMin<<"<"<<m_SourceMinBoundingBox<<")");
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
	//cout<<"SourceFinder::DrawSources(): INFO: Drawing sources..."<<endl;
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

	/*
	if(m_saveImageToFile){
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
		
	}//close if save image to file
	*/

	//cout<<"SourceFinder::Save(): INFO: Closing output file..."<<endl;	
	DEBUG_LOG("Closing output file...");
	m_OutputFile->Close();

	//cout<<"SourceFinder::Save(): INFO: End save to file"<<endl;
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


}//close namespace
