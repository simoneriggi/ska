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
#include <ConfigParser.h>

#include <TObject.h>
#include <TMatrixD.h>
#include <TFile.h>

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

	//Read options
	m_ReadTile= false;
	m_TileMinX= 0;
	m_TileMaxX= 0;
	m_TileMinY= 0;
	m_TileMaxY= 0;
	

}//close InitOptions()

int SourceFinder::Init(){

	//## Init options
	InitOptions();

	//## Configure from parser
	if(Configure()<0){
		cerr<<"SourceFinder::Init(): ERROR: Failed to configure options from parser!"<<endl;
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

int SourceFinder::Run(){

	//## Init options
	if(Init()<0){
		cerr<<"SourceFinder::Run(): ERROR: Initialization failed!"<<endl;
		return -1;
	}
	

	//## Read input image
	if(ReadImage()<0){
		cerr<<"SourceFinder::Run(): ERROR: Reading of input image failed!"<<endl;
		return -1;
	}

	//## Compute stats & bkg for input image
	m_BkgData= ComputeStatsAndBkg(m_InputImg);	
	if(!m_BkgData){
		cerr<<"SourceFinder::Run(): ERROR: Failed to compute image bkg!"<<endl;
		return -1;
	}

	/*
	//## Find bright sources
	if(fSearchBrightSources && FindBrightSource()<0){
		cerr<<"SourceFinder::Run(): ERROR: Bright source search failed!"<<endl;
		return -1;
	}

	//## Find faint sources
	if(fSearchFaintSources && FindFaintSource()<0){
		cerr<<"SourceFinder::Run(): ERROR: Faint source search failed!"<<endl;
		return -1;
	}

	//## Find extended sources
	if(fSearchExtendedSources && FindExtendedSource()<0){
		cerr<<"SourceFinder::Run(): ERROR: Extended source search failed!"<<endl;
		return -1;
	}

	//## Deblend sources
	if(fDeblendSources) DeblendSources(fInputImg);

	//## Draw final sources
	if(fDrawSources) DrawSources(fInputImg,false);
	*/

	//## Save to file
	if(m_SaveToFile) Save();
	
	if(m_Application && m_IsInteractiveRun) m_Application->Run();
	
	return 0;

}//close Run()

int SourceFinder::Configure(){

	//Get image read options
	GET_OPTION_VALUE(inputFile,m_InputFileName);
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
	GET_OPTION_VALUE(minNPix,m_NMinPix);
	GET_OPTION_VALUE(seedBrightThr,m_SeedBrightThr);	
	GET_OPTION_VALUE(seedThr,m_SeedThr);
	GET_OPTION_VALUE(mergeThr,m_MergeThr);


	return 0;

}//close Configure()


int SourceFinder::ReadImage(){

	//## Check file
	FileInfo info;
	bool match_extension= false;
	if(!SysUtils::CheckFile(m_InputFileName,info,match_extension,"")){
		cerr<<"SourceFinder::ReadImage(): ERROR: Invalid input file name specified (invalid file path?)!"<<endl;
		return -1;
	}
	m_InputFileExtension= info.extension;
	

	//=== ROOT reading ===
	if(m_InputFileExtension=="root"){// Read image from ROOT file
		TFile* inputFile = new TFile(m_InputFileName.c_str(),"READ");
		if(!inputFile || inputFile->IsZombie()){
			cerr<<"SourceFinder::ReadImage(): ERROR: Cannot open input file "<<m_InputFileName<<"!"<<endl;
			return -1;
		}
		m_InputImg=  (Img*)inputFile->Get(m_InputImgName.c_str());
		if(!m_InputImg){
			cerr<<"SourceFinder::ReadImage(): ERROR: Cannot get image from input file "<<m_InputFileName<<"!"<<endl;
			return -1;
		}
	}//close if

	//=== FITS reading ===
	else if(m_InputFileExtension=="fits"){// Read image from FITS file
		m_InputImg= new Img;
		int status= 0;
		if(m_ReadTile) status= m_InputImg->ReadFITS(m_InputFileName,m_TileMinX,m_TileMaxX,m_TileMinY,m_TileMaxY); 
		else status= m_InputImg->ReadFITS(m_InputFileName);

		if(status<0){
			cerr<<"SourceFinder::ReadImage(): ERROR: Failed to read image from input file "<<m_InputFileName<<"!"<<endl;
			if(m_InputImg) m_InputImg->Delete();
			return -1;
		}
	}//close else if

	//== Invalid extension ==
	else{
		cerr<<"SourceFinder::ReadImage(): ERROR: Invalid file extension detected!"<<endl;
		return -1;
	}
	m_InputImg->SetNameTitle("img","img");
	
	return 0;

}//close ReadImage()


BkgData* SourceFinder::ComputeStatsAndBkg(Img* img){

	//## Check input img
	if(!img){
		cerr<<"SourceFinder::ComputeStatsAndBkg(): ERROR: Null ptr to input image given!"<<endl;
		return 0;
	}

	//## Compute stats
	cout<<"SourceFinder::ComputeStatsAndBkg(): INFO: Computing image stats..."<<endl;	
	bool computeRobustStats= true;
	bool skipNegativePix= false;
	bool forceRecomputing= false;
	if(!img->ComputeStats(computeRobustStats,skipNegativePix,forceRecomputing)<0){
		cerr<<"SourceFinder::ComputeStatsAndBkg(): ERROR: Stats computing failed!"<<endl;
		return 0;
	}
	img->PrintStats();		

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
		cout<<"SourceFinder::ComputeStatsAndBkg(): INFO: Setting bkg boxes as ("<<m_BoxSizeX<<","<<m_BoxSizeY<<") x beam (beam="<<nPixelsInBeam<<" pixels) ..."<<endl;
		boxSizeX= nPixelsInBeam*m_BoxSizeX;
		boxSizeY= nPixelsInBeam*m_BoxSizeY;
	}
	else{
		cout<<"SourceFinder::ComputeStatsAndBkg(): INFO: Beam information is not available or its usage has been turned off, using image fractions..."<<endl;
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
		cerr<<"SourceFinder::ComputeStatsAndBkg(): ERROR: Bkg computing failed!"<<endl;
		return 0;
	}
		
	return bkgData;

}//close ComputeStatsAndBkg()


int SourceFinder::Save(){

	cout<<"SourceFinder::Save(): INFO: Storing results to file & catalog..."<<endl;

	//Save DS9 regions?
	if(m_SaveDS9Region && m_DS9CatalogFilePtr){
		cout<<"SourceFinder::Save(): INFO: Saving DS9 region header..."<<endl;
	
		fprintf(m_DS9CatalogFilePtr,"global color=red font=\"helvetica 12 normal\" edit=1 move=1 delete=1 include=1\n");
		fprintf(m_DS9CatalogFilePtr,"image\n");

		cout<<"SourceFinder::Save(): INFO: Saving "<<m_SourceCollection.size()<<" sources to file..."<<endl;
		for(unsigned int k=0;k<m_SourceCollection.size();k++){
			cout<<"SourceFinder::Save(): INFO: Dumping DS9 region info for source no. "<<k<<" ..."<<endl;
			std::string regionInfo= "";
			if(m_DS9RegionFormat==1) regionInfo= m_SourceCollection[k]->GetDS9Region(true);
			else if(m_DS9RegionFormat==2) regionInfo= m_SourceCollection[k]->GetDS9EllipseRegion(true);
			else continue;

			fprintf(m_DS9CatalogFilePtr,"%s\n",regionInfo.c_str());
	  
			cout<<"SourceFinder::Save(): INFO: Set source ptr to source "<<k<<"..."<<endl;	
			m_Source= m_SourceCollection[k];
			cout<<"SourceFinder::Save(): INFO: Saving source no. "<<k<<" to tree..."<<endl;	
			m_SourceTree->Fill();
		}//end loop sources

		cout<<"SourceFinder::Save(): INFO: Closing DS9 file region..."<<endl;	
		fclose(m_DS9CatalogFilePtr);
	}//close if SaveDS9Region()


	//Check ROOT output file
	if(!m_OutputFile) {
		cerr<<"SourceFinder::Save(): WARN: Null ptr to output file, nothing will be saved in ROOT file!"<<endl;
		return -1;
	}
	m_OutputFile->cd();

	//Save source tree?
	if(m_SaveSources){
		cout<<"SourceFinder::Save(): INFO: Writing tree to file..."<<endl;	
		m_SourceTree->Write();
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

	cout<<"SourceFinder::Save(): INFO: Closing output file..."<<endl;	
	m_OutputFile->Close();

	cout<<"SourceFinder::Save(): INFO: End save to file"<<endl;

	return 0;

}//close Save()

}//close namespace
