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

#include <SourceFinder.h>
#include <BlobFinder.h>
#include <Img.h>
#include <Source.h>
#include <Contour.h>
#include <ConfigParser.h>
#include <BkgData.h>
#include <CodeUtils.h>
#include <Logger.h>

#include <SLIC.h>
#include <SLICUtils.h>
#include <SLICSegmenter.h>

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

using namespace std;

ClassImp(Caesar::SourceFinder)

namespace Caesar {

SourceFinder::SourceFinder() {
	
}//close costructor


SourceFinder::~SourceFinder(){
	
	//if(m_SourceTree) m_SourceTree->Delete();
	if(m_OutputFile) m_OutputFile->Close();
	if(m_DS9CatalogFilePtr) fclose(m_DS9CatalogFilePtr);
	//if(m_Application) m_Application->Delete();
	
}//close destructor


void SourceFinder::InitOptions(){

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

int SourceFinder::Init(){

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

	return 0;

}//close Init()

int SourceFinder::Configure(){

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


int SourceFinder::Run(){

	//## Init options
	INFO_LOG("Initializing source finder...");
	if(Init()<0){
		ERROR_LOG("Initialization failed!");
		return -1;
	}
	

	//## Read input image
	INFO_LOG("Reading input image...");
	if(ReadImage()<0){
		ERROR_LOG("Reading of input image failed!");
		return -1;
	}

	//## Find compact sources
	INFO_LOG("Searching compact sources...");
	if(m_SearchCompactSources && FindCompactSources()<0){
		ERROR_LOG("Compact source search failed!");
		return -1;
	}


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

	/*
	//## Deblend sources
	if(fDeblendSources) DeblendSources(fInputImg);
	*/

	//## Draw final sources
	if(m_IsInteractiveRun) DrawSources(m_InputImg,m_SourceCollection);

	//## Save to file
	if(m_SaveToFile) Save();	
	if(m_Application && m_IsInteractiveRun) m_Application->Run();
	
	return 0;

}//close Run()


int SourceFinder::FindSources(std::vector<Source*>& sources,Img* inputImg,double seedThr,double mergeThr){

	if(!inputImg) return -1;

	//## Compute stats and bkg
	//cout<<"SourceFinder::FindSources(): INFO: Computing image stats/bkg..."<<endl;
	INFO_LOG("Computing image stats/bkg...");
	BkgData* bkgData= ComputeStatsAndBkg(inputImg);	
	if(!bkgData){
		//cerr<<"SourceFinder::FindSources(): ERROR: Failed to compute stats/bkg info!"<<endl;
		ERROR_LOG("Failed to compute stats/bkg info!");
		return -1;
	}

	//## Compute significance map
	Img* significanceMap= inputImg->GetSignificanceMap(bkgData,m_UseLocalBkg);
	if(!significanceMap){
		//cerr<<"SourceFinder::FindSources(): ERROR: Failed to compute significance map!"<<endl;
		ERROR_LOG("Failed to compute significance map!");
		return -1;
	}

	//## Find sources
	//cout<<"SourceFinder::FindSources(): INFO: Finding sources..."<<endl;
	INFO_LOG("Finding sources...");	
	int status= inputImg->FindCompactSource(sources,significanceMap,bkgData,seedThr,mergeThr,m_NMinPix,m_SearchNegativeExcess,m_MergeBelowSeed,m_SearchNestedSources,m_NestedBlobThrFactor);

	if(status<0) {
		//cerr<<"SourceFinder::FindSources(): ERROR: Source finding failed!"<<endl;
		ERROR_LOG("Source finding failed!");
		delete bkgData;
		bkgData= 0;
		significanceMap->Delete();
		return -1;
	}
	int nSources= (int)sources.size();
	//cout<<"SourceFinder::FindSources(): INFO: "<<nSources<<" sources detected in input image..."<<endl;
	INFO_LOG(nSources<<" sources detected in input image...");	
	
	delete bkgData;
	bkgData= 0;
	significanceMap->Delete();
		

	return 0;

}//close FindSources()


int SourceFinder::FindCompactSources(){

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
	int status= m_InputImg->FindCompactSource(sources,m_SignificanceMap,m_BkgData,m_SeedThr,m_MergeThr,m_NMinPix,m_SearchNegativeExcess,m_MergeBelowSeed,m_SearchNestedSources,m_NestedBlobThrFactor);
	if(status<0) {
		ERROR_LOG("Compact source finding failed!");
		//if(significanceMap) significanceMap->Delete();
		return -1;
	}

	//## Retrieve found sources 
	int nSources= (int)sources.size();
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


	//## Draw significance map
	TCanvas* SignificancePlot= new TCanvas("SignificancePlot","SignificancePlot");
	SignificancePlot->cd();
	m_SignificanceMap->Draw("COLZ");

	//## Clear-up
	//if(significanceMap) significanceMap->Delete();
	
	return 0;

}//close FindCompactSources()

int SourceFinder::FindResidualMap(){

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

int SourceFinder::FindExtendedSources(){

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

int SourceFinder::FindExtendedSources_HClust(Img*){

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

int SourceFinder::FindExtendedSources_ChanVese(Img* inputImg){

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

int SourceFinder::FindExtendedSources_WT(Img* inputImg){

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

int SourceFinder::SelectSources(std::vector<Source*>& sources){

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

bool SourceFinder::IsGoodSource(Source* aSource){
	
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

bool SourceFinder::IsPointLikeSource(Source* aSource){

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



int SourceFinder::ReadImage(){

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
		m_InputImg=  (Img*)inputFile->Get(m_InputImgName.c_str());
		if(!m_InputImg){
			ERROR_LOG("Cannot get image from input file "<<m_InputFileName<<"!");
			return -1;
		}
	}//close if

	//=== FITS reading ===
	else if(m_InputFileExtension==".fits"){// Read image from FITS file
		m_InputImg= new Img;
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
	m_InputImg->SetNameTitle("img","img");
	
	return 0;

}//close ReadImage()


BkgData* SourceFinder::ComputeStatsAndBkg(Img* img){

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


int SourceFinder::DrawSources(Img* image,std::vector<Source*>& sources){

	//cout<<"SourceFinder::DrawSources(): INFO: Drawing sources..."<<endl;
	INFO_LOG("Drawing sources...");
	if(!image) return -1;

	bool useCurrentCanvas= false;
	bool drawFull= false;
	image->Plot(sources,useCurrentCanvas,drawFull,eRAINBOW,true);
	
	return 0;

}//close DrawSources()


int SourceFinder::Save(){

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
