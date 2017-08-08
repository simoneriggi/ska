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
#include <Img.h>
#include <SysUtils.h>
#include <FITSReader.h>
#include <BkgData.h>
#include <BkgFinder.h>
#include <ConfigParser.h>

#include <RInside.h>

#include <TFile.h>

#include <iostream>
#include <vector>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

using namespace std;
using namespace Caesar;

void Usage(char* exeName){
	cout<<"=========== USAGE ==========="<<endl;
	cout<<"Usage: "<<exeName<<" [options]"<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
  cout<<"-h, --help \t Show help message and exit"<<endl;
	cout<<"-c, --config \t Config file containing option settings"<<endl;
	cout<<"=============================="<<endl;
}//close Usage()

static const struct option options_tab[] = {
  { "help", no_argument, 0, 'h' },
	{ "config", required_argument, 0, 'c' },
  {(char*)0, (int)0, (int*)0, (int)0}
};

int ReadImage();
int ComputeBkg();
int ApplySmoothing();
int OpenOutputFile();
int FindSaliency();
int Clear();
int Save();

//Options
//--> Main options
Img* inputImg= 0;
TFile* inputFile= 0;
std::string inputFileName= "";
std::string imageName= "";
TFile* outputFile= 0;	
std::string outputFileName;
Img* saliencyImg= 0;	
bool saveToFile;
bool saveConfig;
bool saveResidualMap;
bool saveBkgMap;
bool saveNoiseMap;
bool saveSignificanceMap;
bool saveInputMap;
bool saveSaliencyMap;
bool saveSources;
bool saveToFITSFile;
std::string residualMapFITSFile;
std::string inputMapFITSFile;
std::string saliencyMapFITSFile;
std::string bkgMapFITSFile;
std::string noiseMapFITSFile;
std::string significanceMapFITSFile;

//--> bkg options
BkgData* bkgData= 0;
Img* significanceMap= 0;
double boxSizeX, boxSizeY;
double gridSizeX, gridSizeY;
int bkgEstimator;
bool useLocalBkg;

int main(int argc, char *argv[]){

	if(argc<2){
		cerr<<"ERROR: Invalid number of arguments...see macro usage!"<<endl;
		Usage(argv[0]);
		exit(1);
	}
	int c = 0;
  int option_index = 0;
	
	//Parse options
	std::string configFileName= "";

	while((c = getopt_long(argc, argv, "hc:",options_tab, &option_index)) != -1) {
    
    switch (c) {
			case 0 : 
			{
				break;
			}
			case 'h':
			{
      	Usage(argv[0]);	
				exit(0);
			}
    	case 'c':	
			{
				configFileName= std::string(optarg);	
				break;	
			}
			default:
			{
      	Usage(argv[0]);	
				exit(0);
			}
    }//close switch
	}//close while
 
	//=======================
	//== Read config options 
	//=======================
	if(ConfigParser::Instance().Parse(configFileName)<0){
		cerr<<"ERROR: Failed to parse config options!"<<endl;
		return -1;
	}
	PRINT_OPTIONS();

	
	//=======================
	//== Open out file
	//=======================
	if(OpenOutputFile()<0){
		cerr<<"ERROR: Failed to open output file!"<<endl;
		return -1;	
	}
	
	//=======================
	//== Read image
	//=======================
	if(ReadImage()<0){
		cerr<<"ERROR: Failed to read image from file!"<<endl;
		return -1;
	}
	
	//=======================
	//== Background finder
	//=======================
	cout<<"INFO: Starting background finder ..."<<endl;
	if(ComputeBkg()<0){
		cerr<<"ERROR: Failed to compute bkg!"<<endl;	
		Clear();
		return -1;
	}
	significanceMap= inputImg->GetSignificanceMap(bkgData,useLocalBkg);
	if(!significanceMap){
		cerr<<"ERROR: Failed to compute significance map!"<<endl;	
		Clear();
		return -1;
	}

	//=======================
	//== Apply smoothing
	//=======================
	bool usePreSmoothing;
	if(GET_OPTION_VALUE(usePreSmoothing,usePreSmoothing)<0){
		cerr<<"ERROR: Failed to get usePreSmoothing option!"<<endl;
		return -1;
	}	
	if(usePreSmoothing){
		if(ApplySmoothing()<0){
			cerr<<"ERROR: Failed to perform image smoothing!"<<endl;
		}
	}

	//=======================
	//== Saliency filtering
	//=======================
	if(FindSaliency()<0){
		cerr<<"ERROR: Failed to compute saliency map!"<<endl;	
		Clear();
		return -1;
	}		

	//=======================
	//== Save results 
	//=======================
	Save();
	

	//=======================
	//== Clear
	//=======================
	Clear();
	
	cout<<"INFO: End saliency finder"<<endl;
	
	return 0;

}//close main

int FindSaliency(){
	
	//## Get options
	int saliencyResoMin= 0;
	int saliencyResoMax = 0; 
	int saliencyResoStep= 0;
	if(GET_OPTION_VALUE(saliencyResoMin,saliencyResoMin)<0){
		cerr<<"FindSaliency(): ERROR: Failed to get saliencyResoMin option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(saliencyResoMax,saliencyResoMax)<0){
		cerr<<"FindSaliency(): ERROR: Failed to get saliencyResoMax option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(saliencyResoStep,saliencyResoStep)<0){
		cerr<<"FindSaliency(): ERROR: Failed to get saliencyResoStep option!"<<endl;
		return -1;
	}

	double spBeta;
	int spMinArea;
	if(GET_OPTION_VALUE(spBeta,spBeta)<0){
		cerr<<"FindSaliency(): ERROR: Failed to get spBeta option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(spMinArea,spMinArea)<0){
		cerr<<"FindSaliency(): ERROR: Failed to get spMinArea option!"<<endl;
		return -1;
	}

	double saliencyNNFactor;
	double saliencySpatialRegFactor;
	bool saliencyUseRobustPars;
	double saliencyMultiResoCombThrFactor;
	double saliencyDissExpFalloffPar;
	double saliencySpatialDistRegPar;
	if(GET_OPTION_VALUE(saliencyNNFactor,saliencyNNFactor)<0){
		cerr<<"FindSaliency(): ERROR: Failed to get saliencyNNFactor option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(saliencyUseRobustPars,saliencyUseRobustPars)<0){
		cerr<<"FindSaliency(): ERROR: Failed to get saliencyUseRobustPars option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(saliencyMultiResoCombThrFactor,saliencyMultiResoCombThrFactor)<0){
		cerr<<"FindSaliency(): ERROR: Failed to get saliencyMultiResoCombThrFactor option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(saliencySpatialRegFactor,saliencySpatialRegFactor)<0){
		cerr<<"FindSaliency(): ERROR: Failed to get saliencySpatialRegFactor option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(saliencyDissExpFalloffPar,saliencyDissExpFalloffPar)<0){
		cerr<<"FindSaliency(): ERROR: Failed to get saliencyDissExpFalloffPar option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(saliencySpatialDistRegPar,saliencySpatialDistRegPar)<0){
		cerr<<"FindSaliency(): ERROR: Failed to get saliencySpatialDistRegPar option!"<<endl;
		return -1;
	}

	bool saliencyUseBkgMap, saliencyUseNoiseMap;
	if(GET_OPTION_VALUE(saliencyUseBkgMap,saliencyUseBkgMap)<0){
		cerr<<"FindSaliency(): ERROR: Failed to get saliencyUseBkgMap option!"<<endl;
		return -1;
	} 
	if(GET_OPTION_VALUE(saliencyUseNoiseMap,saliencyUseNoiseMap)<0){
		cerr<<"FindSaliency(): ERROR: Failed to get saliencyUseNoiseMap option!"<<endl;
		return -1;
	} 

	double saliencyThrFactor, saliencyImgThrFactor;
	if(GET_OPTION_VALUE(saliencyThrFactor,saliencyThrFactor)<0){
		cerr<<"FindSaliency(): ERROR: Failed to get saliencyThrFactor option!"<<endl;
		return -1;
	} 
	if(GET_OPTION_VALUE(saliencyImgThrFactor,saliencyImgThrFactor)<0){
		cerr<<"FindSaliency(): ERROR: Failed to get saliencyImgThrFactor option!"<<endl;
		return -1;
	}
		
	bool saliencyUseCurvInDiss;	
	if(GET_OPTION_VALUE(saliencyUseCurvInDiss,saliencyUseCurvInDiss)<0){
		cerr<<"FindSaliency(): ERROR: Failed to get saliencyUseCurvInDiss option!"<<endl;
		return -1;
	}
		
	cout<<"INFO: saliencyReso("<<saliencyResoMin<<","<<saliencyResoMax<<","<<saliencyResoStep<<") spBeta="<<spBeta<<". spMinArea="<<spMinArea<<", saliencyUseRobustPars="<<saliencyUseRobustPars<<", saliencyUseCurvInDiss="<<saliencyUseCurvInDiss<<" saliencyMultiResoCombThrFactor="<<saliencyMultiResoCombThrFactor<<" saliencyUseBkgMap="<<saliencyUseBkgMap<<" saliencyUseNoiseMap="<<saliencyUseNoiseMap<<" saliencyThrFactor="<<saliencyThrFactor<<", saliencyImgThrFactor="<<saliencyImgThrFactor<<endl;

	//## Compute saliency
	saliencyImg= inputImg->GetMultiResoSaliencyMap(
		saliencyResoMin,saliencyResoMax,saliencyResoStep,
		spBeta,spMinArea,saliencyNNFactor,saliencyUseRobustPars,saliencyDissExpFalloffPar,saliencySpatialDistRegPar, 
		saliencyMultiResoCombThrFactor,
  	saliencyUseBkgMap,saliencyUseNoiseMap,bkgData,
		saliencyThrFactor,saliencyImgThrFactor
	);

	if(!saliencyImg){
		cerr<<"FindSaliency(): ERROR: Failed to compute saliency map!"<<endl;
		return -1;
	}

	return 0;

}//FindSaliency()


int ApplySmoothing(){
	
	//## Get options
	int smoothFilter;
	int gausFilterKernSize;
	double gausFilterSigma;
	double guidedFilterRadius, guidedFilterColorEps;
	if(GET_OPTION_VALUE(smoothFilter,smoothFilter)<0){
		cerr<<"ApplySmoothing(): ERROR: Failed to get smoothFilter option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(gausFilterKernSize,gausFilterKernSize)<0){
		cerr<<"ApplySmoothing(): ERROR: Failed to get gausFilterKernSize option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(gausFilterSigma,gausFilterSigma)<0){
		cerr<<"ApplySmoothing(): ERROR: Failed to get gausFilterSigma option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(guidedFilterRadius,guidedFilterRadius)<0){
		cerr<<"ApplySmoothing(): ERROR: Failed to get guidedFilterRadius option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(guidedFilterColorEps,guidedFilterColorEps)<0){
		cerr<<"ApplySmoothing(): ERROR: Failed to get guidedFilterColorEps option!"<<endl;
		return -1;
	}	

	//## Apply a smoothing stage?
	Img* smoothedImg= 0;
	if(smoothFilter==Img::eGaus){
		smoothedImg= inputImg->GetSmoothedImage(gausFilterKernSize,gausFilterKernSize,gausFilterSigma,gausFilterSigma);
	}
	else if(smoothFilter==Img::eGuided){
		smoothedImg= inputImg->GetGuidedFilterImage(guidedFilterRadius,guidedFilterColorEps);
	}
	else{
		cerr<<"ApplySmoothing(): ERROR: Invalid smoothing algo selected!"<<endl;
		return -1;
	}

	if(!smoothedImg){
		cerr<<"ApplySmoothing(): ERROR: Computation of smoothed map failed!"<<endl;
		return -1;
	}

	// Compute stats
	cout<<"INFO: Computing input image stats..."<<endl;
	if(!smoothedImg->ComputeStats(true,false,false)<0){
		cerr<<"ReadImage(): ERROR: Stats computing failed!"<<endl;
		smoothedImg->Delete();
		return -1;
	}
	
	//## Replace input image with smoothed map
	TString imgName= inputImg->GetName();
	inputImg->Delete();
	smoothedImg->SetNameTitle(imgName,imgName);
	inputImg= smoothedImg;
	
	return 0;

}//close ApplySmoothing()



int ComputeBkg(){

	//## Get bkg options
	if(GET_OPTION_VALUE(boxSizeX,boxSizeX)<0 || GET_OPTION_VALUE(boxSizeY,boxSizeY)<0){
		cerr<<"ComputeBkg(): ERROR: Failed to get boxSize option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(gridSizeX,gridSizeX)<0 || GET_OPTION_VALUE(gridSizeY,gridSizeY)<0){
		cerr<<"ComputeBkg(): ERROR: Failed to get gridSize option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(bkgEstimator,bkgEstimator)<0){
		cerr<<"ComputeBkg(): ERROR: Failed to get bkgEstimator option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(useLocalBkg,useLocalBkg)<0){
		cerr<<"ComputeBkg(): ERROR: Failed to get useLocalBkg option!"<<endl;
		return -1;
	}

	//## Check grid & box size
	int Nx= inputImg->GetNbinsX();
	int Ny= inputImg->GetNbinsY();
	if(boxSizeX>=Nx || boxSizeY>=Ny || gridSizeX>=Nx || gridSizeY>=Ny){
		cerr<<"ComputeBkg(): ERROR: Box/grid size are too large compared to image size ("<<Nx<<","<<Ny<<")"<<endl;
		return -1;
	}

	//Compute bkg & noise maps
	bkgData= inputImg->ComputeBkg(bkgEstimator,useLocalBkg,boxSizeX,boxSizeY,gridSizeX,gridSizeY);
	if(!bkgData){
		cerr<<"ComputeBkg(): ERROR: Failed to compute bkg data!"<<endl;
		return -1;
	}

	return 0;

}//close ComputeBkg()


int ReadImage(){

	//## Get options
	if(GET_OPTION_VALUE(inputFile,inputFileName)<0){
		cerr<<"ReadImage(): ERROR: Failed to get inputFile option!"<<endl;
		return -1;
	}	
	if(GET_OPTION_VALUE(inputImage,imageName)<0){
		cerr<<"ReadImage(): ERROR: Failed to get inputImage option!"<<endl;
		return -1;
	}
	
	//## Check given input file and get info
	Caesar::FileInfo info;
	if(!Caesar::SysUtils::CheckFile(inputFileName,info,false)){
		cerr<<"ReadImage(): ERROR: Invalid input file ("<<inputFileName<<") specified!"<<endl;
		return -1;
	}
	std::string file_extension= info.extension;
	if(file_extension!= ".fits" && file_extension!=".root") {
		cerr<<"ReadImage(): ERROR: Invalid file extension ("<<file_extension<<")...nothing to be done!"<<endl;
		return -1;
	}

	//## Read image
	//===== ROOT reading =====
	if(file_extension==".root"){// Read image from ROOT file
		cout<<"INFO: Reading ROOT input file "<<inputFileName<<"..."<<endl;
		inputFile = new TFile(inputFileName.c_str(),"READ");
		if(!inputFile || inputFile->IsZombie()){
			cerr<<"ReadImage(): ERROR: Cannot open input file "<<inputFileName<<"!"<<endl;
			return -1;
		}
		inputImg=  (Img*)inputFile->Get(imageName.c_str());
		if(!inputImg){
			cerr<<"ReadImage(): ERROR: Cannot get image from input file "<<inputFileName<<"!"<<endl;
			return -1;
		}
	}//close if

	//===== FITS reading =====
	if(file_extension==".fits"){// Read image from FITS file
		cout<<"ReadImage(): INFO: Reading FITS input file "<<inputFileName<<"..."<<endl;
		inputImg= new Caesar::Img; 
		inputImg->SetNameTitle(imageName.c_str(),imageName.c_str());
		if(inputImg->ReadFITS(inputFileName)<0){
			cerr<<"ReadImage(): ERROR: Failed to read image from input file "<<inputFileName<<"!"<<endl;	
			return -1;
		}
	}//close else if

	if(!inputImg){
		cerr<<"ReadImage(): ERROR: Failed to read image from input file "<<inputFileName<<"!"<<endl;	
		return -1;
	}

	// Compute stats
	cout<<"INFO: Computing input image stats..."<<endl;
	if(inputImg->ComputeStats(true,false,false)<0){
		cerr<<"ReadImage(): ERROR: Stats computing failed!"<<endl;
		inputImg->Delete();
		return -1;
	}
	inputImg->PrintStats();	

	return 0;

}//close ReadImage()

int OpenOutputFile(){

	//Get options
	if(GET_OPTION_VALUE(outputFile,outputFileName)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get outputFile option!"<<endl;
		return -1;
	}	
	if(GET_OPTION_VALUE(saveToFile,saveToFile)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get saveToFile option!"<<endl;
		return -1;
	}	
	if(GET_OPTION_VALUE(saveConfig,saveConfig)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get saveConfig option!"<<endl;
		return -1;
	}	
	if(GET_OPTION_VALUE(saveResidualMap,saveResidualMap)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get saveResidualMap option!"<<endl;
		return -1;
	}	
	if(GET_OPTION_VALUE(saveBkgMap,saveBkgMap)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get saveBkgMap option!"<<endl;
		return -1;
	}	
	if(GET_OPTION_VALUE(saveNoiseMap,saveNoiseMap)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get saveNoiseMap option!"<<endl;
		return -1;
	}		
	if(GET_OPTION_VALUE(saveSignificanceMap,saveSignificanceMap)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get saveSignificanceMap option!"<<endl;
		return -1;
	}	
	if(GET_OPTION_VALUE(saveInputMap,saveInputMap)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get saveInputMap option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(saveSaliencyMap,saveSaliencyMap)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get saveSaliencyMap option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(saveSources,saveSources)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get saveSources option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(saveToFITSFile,saveToFITSFile)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get saveToFITSFile option!"<<endl;
		return -1;
	}

	if(GET_OPTION_VALUE(residualMapFITSFile,residualMapFITSFile)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get residualMapFITSFile option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(inputMapFITSFile,inputMapFITSFile)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get inputMapFITSFile option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(saliencyMapFITSFile,saliencyMapFITSFile)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get saliencyMapFITSFile option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(bkgMapFITSFile,bkgMapFITSFile)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get bkgMapFITSFile option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(noiseMapFITSFile,noiseMapFITSFile)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get noiseMapFITSFile option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(significanceMapFITSFile,significanceMapFITSFile)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get significanceMapFITSFile option!"<<endl;
		return -1;
	}

	if(saveToFile){
		outputFile= new TFile(outputFileName.c_str(),"RECREATE");
		if(!outputFile || !outputFile->IsOpen()){
			cerr<<"OpenOutputFile(): ERROR: Failed to open output file!"<<endl;
			return -1;
		}
	}
	
	return 0;

}//close OpenOutputFile()


int Clear(){

	if(inputImg) inputImg->Delete();
	if(bkgData){ 
		delete bkgData;
		bkgData= 0;
	}
	if(saliencyImg) saliencyImg->Delete();

	return 0;

}//close Clear()

int Save(){

	//## Save to ROOT?
	if(saveToFile && outputFile){
		outputFile->cd();
		if(saveConfig){
			TTree* configTree= ConfigParser::Instance().GetConfigTree();
			if(configTree) configTree->Write();
		}
		if(saveInputMap && inputImg) inputImg->Write();
		if(saveBkgMap && bkgData->BkgMap) {
			(bkgData->BkgMap)->Write();
		}
		if(saveNoiseMap && bkgData->NoiseMap) {
			(bkgData->NoiseMap)->Write();
		}
		if(saveSignificanceMap && significanceMap) {	
			significanceMap->Write();
		}
		if(saveSaliencyMap && saliencyImg) saliencyImg->Write();
		outputFile->Close();
	}

	//## Save to FITS?
	if(saveToFITSFile){
		if(saveInputMap && inputImg) inputImg->WriteFITS(inputMapFITSFile);
		if(saveBkgMap && bkgData->BkgMap) {
			(bkgData->BkgMap)->WriteFITS(bkgMapFITSFile);
		}
		if(saveNoiseMap && bkgData->NoiseMap) {
			(bkgData->NoiseMap)->WriteFITS(noiseMapFITSFile);
		}
		if(saveSignificanceMap && significanceMap) {	
			significanceMap->WriteFITS(significanceMapFITSFile);
		}
		if(saveSaliencyMap && saliencyImg) saliencyImg->WriteFITS(saliencyMapFITSFile);
	}

	return 0;
}//close Save()

