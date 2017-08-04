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
#include <ImgFITSReader.h>
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
//BkgData* bkgData= 0;
Img* significanceMap= 0;
double boxSizeX, boxSizeY;
double gridSizeX, gridSizeY;
int bkgEstimator;
bool useLocalBkg;

enum SmoothingFilter {eGaus=1,eGuided=2,eBilateral=3};

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
 

	
	//## Create RInside instance
	RInside* fR= new RInside;
	
	//=======================
	//== Read config options 
	//=======================
	//Check config file
	if(configFileName==""){
		cerr<<"ERROR: Invalid or empty config filename, see program usage!"<<endl;
		Usage(argv[0]);
		exit(1);
	}
	cout<<"INFO: Reading config from file "<<configFileName<<"..."<<endl;
	
	ConfigParser parser(configFileName);
	parser.ReadConfig();
	parser.Print();

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
	cout<<"INFO: Reading input file..."<<endl;
	if(ReadImage()<0){
		cerr<<"ERROR: Failed to read image from file!"<<endl;
		return -1;
	}
	
	//=======================
	//== Apply smoothing
	//=======================
	cout<<"INFO: Smoothing the input map? "<<ConfigParser::fUsePreSmoothing<<endl;
	if(ConfigParser::fUsePreSmoothing && ApplySmoothing()<0){
		cerr<<"ERROR: Failed to perform image smoothing!"<<endl;
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


	//=======================
	//== Saliency filtering
	//=======================
	cout<<"INFO: Computing saliency map ..."<<endl;
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
	int saliencyResoMin= ConfigParser::fSaliencyMinReso;
	int saliencyResoMax = ConfigParser::fSaliencyMaxReso; 
	int saliencyResoStep= ConfigParser::fSaliencyResoStepSize;
	double saliencyNNFactor= ConfigParser::fSaliencyNNFactor;
	bool saliencyUseRobustPars= ConfigParser::fSaliencyUseRobustPars;
	bool useCurvatureInSPMerging= ConfigParser::fUseCurvatureInSPMerging;
	double spBeta= ConfigParser::fSPRegularization;
	int spMinArea= ConfigParser::fSPMinArea;
	double saliencyFilterThresholdFactor= ConfigParser::fSaliencyFilterThresholdFactor;
	bool saliencyUseBkgMap= ConfigParser::fSaliencyUseBkgMap;
	bool saliencyUseNoiseMap= ConfigParser::fSaliencyUseNoiseMap;
	bool saliencyUseCurvatureMap= ConfigParser::fSaliencyUseCurvatureMap;
	int saliencyNormalizationMode= ConfigParser::fSaliencyNormalizationMode;
	double saliencyThresholdFactor= ConfigParser::fSaliencyThresholdFactor;
	double saliencyImgThresholdFactor= ConfigParser::fSaliencyImgThresholdFactor;

	double saliencyDissExpFalloffPar= ConfigParser::fSaliencyDissExpFalloffPar;
	double saliencySpatialDistRegPar= ConfigParser::fSaliencySpatialDistRegPar;

	//## Compute saliency map 
	cout<<"FindSaliency(): INFO: Computing saliency map: reso pars("<<saliencyResoMin<<","<<saliencyResoMax<<","<<saliencyResoStep<<"), knn="<<saliencyNNFactor<<", thr="<<saliencyFilterThresholdFactor<<" ..."<<endl;
	/*
	saliencyImg= inputImg->GetMultiResoSaliencyMap (
    saliencyResoMin, saliencyResoMax, saliencyResoStep,
		spBeta, spMinArea,
		saliencyNNFactor, saliencyUseRobustPars, useCurvatureInSPMerging,
		saliencyFilterThresholdFactor,
    saliencyUseCurvatureMap, saliencyUseBkgMap, saliencyUseNoiseMap,
		saliencyNormalizationMode,
		saliencyThresholdFactor,saliencyImgThresholdFactor
	);
	*/
	saliencyImg= inputImg->GetMultiResoSaliencyMap (
    saliencyResoMin, saliencyResoMax, saliencyResoStep,
		spBeta, spMinArea,
		saliencyNNFactor, saliencyUseRobustPars, saliencyDissExpFalloffPar, saliencySpatialDistRegPar,
		saliencyFilterThresholdFactor,
    saliencyUseCurvatureMap, saliencyUseBkgMap, saliencyUseNoiseMap,
		saliencyNormalizationMode,
		saliencyThresholdFactor,saliencyImgThresholdFactor
	);
	
	//## Compute saliency stats
	cout<<"FindSaliency(): INFO: Computing saliency map stats..."<<endl;
	if(saliencyImg->ComputeStats(true,false,true)<0){
		cerr<<"FindSaliency(): ERROR: Saliency stats computing failed!"<<endl;
		return -1;
	}	

	return 0;

}//FindSaliency()


int ApplySmoothing(){
		
	//## Apply a smoothing stage?
	Img* smoothedImg= 0;
	if(ConfigParser::fSmoothingAlgo==eGaus){
		smoothedImg= inputImg->Smooth(ConfigParser::fSmoothKernelSize,ConfigParser::fSmoothKernelSize,ConfigParser::fSmoothSigma,ConfigParser::fSmoothSigma);
	}
	else if(ConfigParser::fSmoothingAlgo==eGuided){
		smoothedImg= inputImg->GetGuidedFilterImage(ConfigParser::fGuidedSmoothRadius,ConfigParser::fGuidedSmoothColorEps);
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
	cout<<"ApplySmoothing::INFO: Computing smoothed image stats..."<<endl;
	if(smoothedImg->ComputeStats(true,false,true)<0){
		cerr<<"ApplySmoothing(): ERROR: Stats computing failed!"<<endl;
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
	
	//## Compute local image bkg
	if( inputImg->ComputeLocalBkg(
				(Img::LocalBkgMethod)ConfigParser::fLocalBkgMethod, (Img::BkgMethod)ConfigParser::fBkgEstimator,
        ConfigParser::fBoxSizeX,ConfigParser::fBoxSizeY,
        ConfigParser::fGridSizeX,ConfigParser::fGridSizeY
			)<0
	) {
		cerr<<"ComputeBkg(): ERROR: Failed to compute local bkg!"<<endl;
		return -1;
	}

	return 0;

}//close ComputeBkg()


int ReadImage(){

	//Get input file name
	std::string inputFileName= ConfigParser::fInputFileName;
	if (inputFileName == "") {
		cerr<<"ReadImage(): ERROR: Empty input file name!"<<endl;
  	return -1;
  }
	std::string inputFileExtension= inputFileName.substr(inputFileName.find_last_of(".") + 1);
	if(inputFileExtension!= "fits" && inputFileExtension!="root") {
		cerr<<"ReadImage(): ERROR: Invalid file extension ("<<inputFileExtension<<")...nothing to be done!"<<endl;
		return -1;
	}

	//--> ROOT reading
	if(inputFileExtension=="root"){// Read image from ROOT file
		TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
		if(!inputFile || inputFile->IsZombie()){
			cerr<<"ReadImage(): ERROR: Cannot open input file "<<inputFileName<<"!"<<endl;
			return -1;
		}
		inputImg=  (Img*)inputFile->Get( (ConfigParser::fROOTInputImgName).c_str() );
		if(!inputImg){
			cerr<<"ReadImage(): ERROR: Cannot get image from input file "<<inputFileName<<"!"<<endl;
			return -1;
		}
	}//close if

	//--> FITS reading
	else if(inputFileExtension=="fits"){// Read image from FITS file
		ImgFITSReader reader(inputFileName.c_str());	
		reader.ReadHeader();

		ImgFITSReader::FITSHeader header= reader.GetHeaderInfo();
		int Nx= header.Nx;
		int Ny= header.Ny;
		double Bmaj= header.Bmaj;
		double Bmin= header.Bmin;
		double dX= header.dX;
		double dY= header.dY;
		int nBeamPix= fabs(Bmaj/dX);
		cout<<"ReadImage(): INFO: Header info: Nx="<<Nx<<" Ny="<<Ny<<" nBeamPix="<<nBeamPix<<endl;
	
		if(inputImg) inputImg->Delete();
		inputImg= new Img;
		reader.Read(*inputImg);
	}//close else if

	//--> Invalid extension
	else{
		cerr<<"ReadImage(): ERROR: Invalid file extension detected!"<<endl;
		return -1;
	}
	inputImg->SetNameTitle("img","img");	

	//## Compute stats
	cout<<"ReadImage(): INFO: Computing input image stats..."<<endl;
	if(inputImg->ComputeStats(true,false,true)<0){
		cerr<<"ReadImage(): ERROR: Stats computing failed!"<<endl;
		return -1;
	}
	inputImg->DumpStats();		


	return 0;

}//close ReadImage()

int OpenOutputFile(){

	//Check if save to file option is ON
	if(!ConfigParser::fSaveToFile) return 0;

	//Open output file
	outputFile= new TFile((ConfigParser::fOutputFileName).c_str(),"RECREATE");
	if(!outputFile || !outputFile->IsOpen()){
		cerr<<"OpenOutputFile(): ERROR: Failed to open output file!"<<endl;
		return -1;
	}
	
	return 0;

}//close OpenOutputFile()


int Clear(){

	//Clear images
	if(inputImg) inputImg->Delete();
	if(saliencyImg) saliencyImg->Delete();

	return 0;

}//close Clear()

int Save(){

	//## Save to ROOT?
	if(ConfigParser::fSaveToFile && outputFile && outputFile->IsOpen()){
		outputFile->cd();
		
		if(inputImg) inputImg->Write();
		if(saliencyImg) {		
			saliencyImg->SetNameTitle("saliencyMap","saliencyMap");
			saliencyImg->Write();
		}
		outputFile->Close();
	}

	return 0;

}//close Save()

