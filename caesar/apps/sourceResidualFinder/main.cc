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
#include <Contour.h>
#include <Source.h>

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
int OpenOutputFile();
int FindSources();
int SelectSources();
bool IsGoodSource(Source* aSource);
bool IsPointLikeSource(Source* aSource);
int ComputeSourceResidual();
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
Img* residualImg= 0;
std::vector<Source*> sources;	
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

//--> Source selection
bool applySourceSelection;
double sourceMinBoundingBox;
double psCircRatioThr, psElongThr, psEllipseAreaRatioMinThr, psEllipseAreaRatioMaxThr, psMaxNPix;

//--> bkg options
BkgData* bkgData= 0;
Img* significanceMap= 0;
double boxSizeX, boxSizeY;
double gridSizeX, gridSizeY;
int bkgEstimator;
bool useLocalBkg;
bool use2ndPassInLocalBkg;
bool skipOutliersInLocalBkg;

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


	cout<<"Create Rinstance..."<<endl;
	RInside::instancePtr();
	
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
	//== Source finding
	//=======================
	if(FindSources()<0){
		cerr<<"ERROR: Failed to find sources!"<<endl;	
		Clear();
		return -1;
	}		

	//=======================
	//== Source selection
	//=======================
	if(SelectSources()<0){
		cerr<<"ERROR: Failed to select sources!"<<endl;	
		Clear();
		return -1;
	}
	
	//=======================
	//== Source residuals
	//=======================
	if(ComputeSourceResidual()<0){
		cerr<<"ERROR: Failed to compute source residual map!"<<endl;	
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
	
	cout<<"INFO: End compact source finder"<<endl;
	
	return 0;

}//close main

int ComputeSourceResidual(){

	//Get options
	bool dilateNestedSources;
	int dilateKernelSize;
	int dilatedSourceType;
	int dilateSourceModel;
	if(GET_OPTION_VALUE(dilateNestedSources,dilateNestedSources)<0){
		cerr<<"ComputeSourceResidual(): ERROR: Failed to get dilateNestedSources option!"<<endl;
		return -1;
	}	
	if(GET_OPTION_VALUE(dilateKernelSize,dilateKernelSize)<0){
		cerr<<"ComputeSourceResidual(): ERROR: Failed to get dilateKernelSize option!"<<endl;
		return -1;
	}	
	if(GET_OPTION_VALUE(dilatedSourceType,dilatedSourceType)<0){
		cerr<<"ComputeSourceResidual(): ERROR: Failed to get dilatedSourceType option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(dilateSourceModel,dilateSourceModel)<0){
		cerr<<"ComputeSourceResidual(): ERROR: Failed to get dilateSourceModel option!"<<endl;
		return -1;
	}

	//Compute residual
	residualImg= inputImg->GetSourceResidual(sources,dilateKernelSize,dilateSourceModel,dilatedSourceType,dilateNestedSources,bkgData,useLocalBkg);
	if(!residualImg){
		cerr<<"ComputeSourceResidual(): ERROR: Failed to compute residual map!"<<endl;
		return -1;
	}

	return 0;

}//close ComputeSourceResidual()

int FindSources(){
	
	//## Get options
	double seedBrightThr, mergeThr;
	int minNPix;
	bool searchNegativeExcess;
	bool mergeBelowSeed;
	bool searchNestedSources;
	double nestedBlobThrFactor;

	if(GET_OPTION_VALUE(seedBrightThr,seedBrightThr)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get seedBrightThr option!"<<endl;
		return -1;
	}	
	if(GET_OPTION_VALUE(mergeThr,mergeThr)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get mergeThr option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(minNPix,minNPix)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get minNPix option!"<<endl;
		return -1;
	}	
	if(GET_OPTION_VALUE(searchNegativeExcess,searchNegativeExcess)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get searchNegativeExcess option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(mergeBelowSeed,mergeBelowSeed)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get mergeBelowSeed option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(searchNestedSources,searchNestedSources)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get searchNestedSources option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(nestedBlobThrFactor,nestedBlobThrFactor)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get nestedBlobThrFactor option!"<<endl;
		return -1;
	}

	sources.clear();
	int status= inputImg->FindCompactSource(sources,significanceMap,bkgData,seedBrightThr,mergeThr,minNPix,searchNegativeExcess,mergeBelowSeed,searchNestedSources,nestedBlobThrFactor);
	if(status<0){
		cerr<<"FindSources(): ERROR: Bright source search failed!"<<endl;
		return -1;
	}

	return 0;

}//close FindSources()

int SelectSources(){

	//## Get options	
	if(GET_OPTION_VALUE(applySourceSelection,applySourceSelection)<0){
		cerr<<"SelectSources(): ERROR: Failed to get applySourceSelection option!"<<endl;
		return -1;
	}
	if(!applySourceSelection){
		cout<<"SelectSources(): WARN: No source selection requested!"<<endl;
		return 0;
	}
	
	if(GET_OPTION_VALUE(sourceMinBoundingBox,sourceMinBoundingBox)<0){
		cerr<<"SelectSources(): ERROR: Failed to get sourceMinBoundingBox option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(psCircRatioThr,psCircRatioThr)<0){
		cerr<<"SelectSources(): ERROR: Failed to get psCircRatioThr option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(psElongThr,psElongThr)<0){
		cerr<<"SelectSources(): ERROR: Failed to get psElongThr option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(psEllipseAreaRatioMinThr,psEllipseAreaRatioMinThr)<0){
		cerr<<"SelectSources(): ERROR: Failed to get psEllipseAreaRatioMinThr option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(psEllipseAreaRatioMaxThr,psEllipseAreaRatioMaxThr)<0){
		cerr<<"SelectSources(): ERROR: Failed to get psEllipseAreaRatioMaxThr option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(psMaxNPix,psMaxNPix)<0){
		cerr<<"SelectSources(): ERROR: Failed to get psMaxNPix option!"<<endl;
		return -1;
	}
	

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
			cout<<"SelectSources(): INFO: Source no. "<<i<<" (name="<<sourceName<<",id="<<sourceId<<", n="<<NPix<<"("<<X0<<","<<Y0<<")) tagged as bad source, skipped!"<<endl;
			sources[i]->SetGoodSourceFlag(false);
			continue;
		}
			
		//Is point-like source?
		if( IsPointLikeSource(sources[i]) ){
			cout<<"SelectSources(): INFO: Source no. "<<i<<" (name="<<sourceName<<",id="<<sourceId<<", n="<<NPix<<"("<<X0<<","<<Y0<<")) tagged as a point-like source ..."<<endl;
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
				cout<<"SelectSources(): INFO: Source no. "<<i<<": nested source no. "<<j<<" (name="<<nestedSourceName<<",id="<<nestedSourceId<<", n="<<nestedNPix<<"("<<nestedX0<<","<<nestedY0<<")) tagged as bad source, skipped!"<<endl;
				nestedSources[j]->SetGoodSourceFlag(false);
			}
			if( IsPointLikeSource(nestedSources[j]) ){
				cout<<"SelectSources(): INFO: Source no. "<<i<<": nested source no. "<<j<<" (name="<<nestedSourceName<<",id="<<nestedSourceId<<", n="<<nestedNPix<<"("<<nestedX0<<","<<nestedY0<<")) tagged as a point-like source ..."<<endl;
				nestedSources[j]->SetType(Source::ePointLike);
			}
		}//end loop nested sources
			
		//Add source to the list	
		sources_sel.push_back(sources[i]);
		nSelSources++;
	}//end loop sources

	
	cout<<"SelectSources(): INFO: Added "<<nSelSources<<" bright sources to the selected list..."<<endl;	

	//Clear initial vector (DO NOT CLEAR MEMORY!) and fill with selection (then reset selection)
	sources.clear();
	sources.insert(sources.end(),sources_sel.begin(),sources_sel.end());
	sources_sel.clear();

	return 0;

}//SelectSources()

bool IsGoodSource(Source* aSource){
	
	if(!aSource) return false;


	//## Check for pixels 	
	if(aSource->NPix<=0 || (aSource->GetPixels()).size()<=0) return false;

	//## Check for line-like source
	if( (aSource->GetContours()).size()<=0) {
		cerr<<"IsGoodSource(): WARN: No contour stored for this source, cannot perform check!"<<endl;
		return true;
	}

	double BoundingBoxMin= ((aSource->GetContours())[0])->BoundingBoxMin;
	if(BoundingBoxMin<sourceMinBoundingBox) {
		cerr<<"IsGoodSource(): INFO: BoundingBox cut not passed (BoundingBoxMin="<<BoundingBoxMin<<"<"<<sourceMinBoundingBox<<")"<<endl;
		return false;
	}

	//## Add other check here ...
	//...
	//...

	return true;

}//close IsGoodSource()

bool IsPointLikeSource(Source* aSource){

	if(!aSource) return false;
	if(!aSource->HasParameters()) {
		cerr<<"IsPointLikeSource(): WARN: No parameters are available for this source (did you compute them?)...test cannot be performed!"<<endl;
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
		if(thisContour->CircularityRatio<psCircRatioThr) {
			cout<<"SourceFinder::IsCompactSource(): INFO: Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass CircularityRatio cut (CR="<<thisContour->CircularityRatio<<"<"<<psCircRatioThr<<")"<<endl;
			isPointLike= false;
			break;
		}
		*/

		//Test elongation (how symmetrical is the shape): 0=circle,square
		if(thisContour->Elongation>psElongThr) {
			cout<<"IsPointLikeSource: INFO: Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass Elongation cut (ELONG="<<thisContour->CircularityRatio<<">"<<psElongThr<<")"<<endl;
			isPointLike= false;
			break;	
		}

		//Test ellipse fit
		if(thisContour->EllipseAreaRatio<psEllipseAreaRatioMinThr || thisContour->EllipseAreaRatio>psEllipseAreaRatioMaxThr) {
			cout<<"IsPointLikeSource: INFO: Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass EllipseAreaRatio cut (EAR="<<thisContour->EllipseAreaRatio<<" outside range ["<<psEllipseAreaRatioMinThr<<","<<psEllipseAreaRatioMaxThr<<"])"<<endl;
			isPointLike= false;
			break;	
		}

	}//end contour loop
	
	//Check number of pixels
	cout<<"IsPointLikeSource(): INFO: Source (name="<<sourceName<<","<<"id="<<sourceId<<") (NPix="<<aSource->NPix<<">"<<psMaxNPix<<")"<<endl;
	if(aSource->NPix>psMaxNPix){
		cout<<"IsPointLikeSource(): INFO: Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass nMaxPix cut (NPix="<<aSource->NPix<<">"<<psMaxNPix<<")"<<endl;
		isPointLike= false;
	}

	if(!isPointLike) return false;

	return true;

}//close IsPointLikeSource()

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
	if(GET_OPTION_VALUE(use2ndPassInLocalBkg,use2ndPassInLocalBkg)<0){
		cerr<<"ComputeBkg(): ERROR: Failed to get use2ndPassInLocalBkg option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(skipOutliersInLocalBkg,skipOutliersInLocalBkg)<0){
		cerr<<"ComputeBkg(): ERROR: Failed to get skipOutliersInLocalBkg option!"<<endl;
		return -1;
	}

	//## Get options
	double seedBrightThr, mergeThr;
	int minNPix;
	
	if(GET_OPTION_VALUE(seedBrightThr,seedBrightThr)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get seedBrightThr option!"<<endl;
		return -1;
	}	
	if(GET_OPTION_VALUE(mergeThr,mergeThr)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get mergeThr option!"<<endl;
		return -1;
	}
	if(GET_OPTION_VALUE(minNPix,minNPix)<0){
		cerr<<"OpenOutputFile(): ERROR: Failed to get minNPix option!"<<endl;
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
	bkgData= inputImg->ComputeBkg(bkgEstimator,useLocalBkg,boxSizeX,boxSizeY,gridSizeX,gridSizeY,use2ndPassInLocalBkg,skipOutliersInLocalBkg,seedBrightThr,mergeThr,minNPix);
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
	for(unsigned int i=0;i<sources.size();i++){
		if(sources[i]){
			delete sources[i];
			sources[i]= 0;
		}
	}
	sources.clear();
	if(residualImg) residualImg->Delete();

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
		if(saveResidualMap && residualImg) residualImg->Write();
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
		if(saveResidualMap && residualImg) residualImg->WriteFITS(residualMapFITSFile);
	}

	return 0;
}//close Save()

