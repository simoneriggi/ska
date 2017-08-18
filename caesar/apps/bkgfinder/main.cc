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
	cout<<"-i, --input \t Input file name containing image to be read (only .fits/.root supported)"<<endl;
	cout<<"-b, --boxsize \t Size of sampling box in pixels (often a multiple of the image beam size)"<<endl;
	cout<<"-g, --gridsize \t Granularity size of the interpolation grid in pixels (i.e. half the box size)"<<endl;
	cout<<"-e, --estimator \t Bkg estimator used in the sampling box (1=mean, 2=median, 3=biweight, 4=clipped median)"<<endl;
	cout<<"-f, --fitsout \t Write results in FITS files"<<endl;
	cout<<"-o, --output \t Output file name (1 file for ROOT out and multiple files if FITS out option is selected)"<<endl;
	cout<<"-I, --imgname \t Image name in input ROOT file (if non standard)"<<endl;
	cout<<"-s, --significance \t Compute and store also the significance map (along with bkg and noise maps)"<<endl;
	cout<<"=============================="<<endl;
}//close Usage()

static const struct option options_tab[] = {
  { "help", no_argument, 0, 'h' },
	{ "input", required_argument, 0, 'i' },
	{ "output", optional_argument, 0, 'o' },
	{ "fitsout", no_argument, 0, 'f' },
	{ "gridsize", required_argument, 0, 'g' },
	{ "boxsize", required_argument, 0, 'b' },
	{ "estimator", required_argument, 0, 'e' },
	{ "significance", no_argument, 0, 's' },
	{ "imgname", optional_argument, 0, 'I' },
  {(char*)0, (int)0, (int*)0, (int)0}
};




int main(int argc, char *argv[]){

	if(argc<2){
		cerr<<"ERROR: Invalid number of arguments...see macro usage!"<<endl;
		Usage(argv[0]);
		exit(1);
	}
	int c = 0;
  int option_index = 0;
	
	//Parse options
	bool useDefaultOutput= true;
	std::string inputFileName= "";
	std::string outputFileName= "";
	std::string outputFileName_bkg= "";
	std::string outputFileName_noise= "";
	std::string outputFileName_significance= "";
	bool writeToFITS= false;
	double boxSize= 0;
	double gridSize= 0;
	int bkgEstimator= 0;
	bool computeSignificance= false;
	std::string imageName= "img";
	Caesar::Img* inputImg= 0;
	TFile* outputFile= 0;

	while((c = getopt_long(argc, argv, "hfsi:o::b:g:e:I::",options_tab, &option_index)) != -1) {
    
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
    	case 'i':	
			{
				inputFileName= std::string(optarg);	
				break;	
			}
			case 'o':	
			{
				outputFileName= std::string(optarg);	
				useDefaultOutput= false;
				break;	
			}
			case 'f':
			{
      	writeToFITS= true;
				break;
			}
			case 'b':
			{
				boxSize= atof(optarg);
				break;
			}
			case 'g':
			{
				gridSize= atof(optarg);
				break;
			}
			case 'e':
			{
				bkgEstimator= atoi(optarg);
				break;
			}
			case 's':	
			{
				computeSignificance= true;
				break;
			}
			case 'I':	
			{
				imageName= std::string(optarg);	
				break;	
			}
			default:
			{
      	Usage(argv[0]);	
				exit(0);
			}
    }//close switch
	}//close while
 
	
	//Set bkg estimator
	Caesar::BkgEstimator estimator= Caesar::eMedianBkg;
	if(bkgEstimator==Caesar::eMedianBkg){
		estimator= Caesar::eMedianBkg;
	}
	else if(bkgEstimator==Caesar::eMeanBkg){
		estimator= Caesar::eMeanBkg;
	}
	else if(bkgEstimator==Caesar::eBiWeightBkg){
		estimator= Caesar::eBiWeightBkg;
	}	
	else if(bkgEstimator==Caesar::eMedianClippedBkg){
		estimator= Caesar::eMedianClippedBkg;
	}
	else{
		cerr<<"ERROR: Invalid bkg estimator specified!"<<endl;
		if(inputImg) inputImg->Delete();
		return -1;
	}

	//=======================
	//== Read image 
	//=======================
	// Check given input file and get info
	Caesar::FileInfo info;
	if(!Caesar::SysUtils::CheckFile(inputFileName,info,false)){
		cerr<<"ERROR: Invalid input file ("<<inputFileName<<") specified!"<<endl;
		return 0;
	}
	std::string file_extension= info.extension;
	if(file_extension!= ".fits" && file_extension!=".root") {
		cerr<<"ERROR: Invalid file extension ("<<file_extension<<")...nothing to be done!"<<endl;
		return 0;
	}
	
	// Set output filenames
	if(useDefaultOutput){
		std::string basefilename_wext= info.filename_wext;
		if(writeToFITS){//FITS output
			std::string outputFileNamePrefix= basefilename_wext;
			outputFileName_bkg= outputFileNamePrefix + std::string("_bkg.fits");
			outputFileName_noise= outputFileNamePrefix + std::string("_noise.fits");	
			outputFileName_significance= outputFileNamePrefix + std::string("_significance.fits");
		}
		else{//ROOT output
			outputFileName= basefilename_wext + std::string("_bkg.root");
			outputFile= new TFile(outputFileName.c_str(),"RECREATE");
		}
	}//close if
	else{
		if(writeToFITS){//FITS output
			std::string outputFileNamePrefix= outputFileName;
			outputFileName_bkg= outputFileNamePrefix + std::string("_bkg.fits");
			outputFileName_noise= outputFileNamePrefix + std::string("_noise.fits");	
			outputFileName_significance= outputFileNamePrefix + std::string("_significance.fits");
		}
		else{
			outputFile= new TFile(outputFileName.c_str(),"RECREATE");
		}
	}//close else
	


	//--> ROOT reading
	if(file_extension==".root"){// Read image from ROOT file
		cout<<"INFO: Reading ROOT input file "<<inputFileName<<"..."<<endl;
		TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
		if(!inputFile || inputFile->IsZombie()){
			cerr<<"ERROR: Cannot open input file "<<inputFileName<<"!"<<endl;
			return -1;
		}
		inputImg=  (Caesar::Img*)inputFile->Get(imageName.c_str());
		if(!inputImg){
			cerr<<"ERROR: Cannot get image from input file "<<inputFileName<<"!"<<endl;
			return -1;
		}
	}//close if

	//--> FITS reading
	if(file_extension==".fits"){// Read image from FITS file
		cout<<"INFO: Reading FITS input file "<<inputFileName<<"..."<<endl;
		inputImg= new Caesar::Img; 
		inputImg->SetNameTitle(imageName.c_str(),imageName.c_str());
		if(inputImg->ReadFITS(inputFileName)<0){
			cerr<<"ERROR: Failed to read image from input file "<<inputFileName<<"!"<<endl;	
			return -1;
		}
	}//close else if

	if(!inputImg){
		cerr<<"ERROR: Failed to read image from input file "<<inputFileName<<"!"<<endl;	
		return -1;
	}

	// Compute stats
	cout<<"INFO: Computing input image stats..."<<endl;
	if(inputImg->ComputeStats(true,false,false)<0){
		cerr<<"ReadImage(): ERROR: Stats computing failed!"<<endl;
		return -1;
	}
	inputImg->PrintStats();	

	//=======================
	//== Background finder
	//=======================
	cout<<"INFO: Starting background finder ..."<<endl;
	
	//Check grid & box size
	int Nx= inputImg->GetNbinsX();
	int Ny= inputImg->GetNbinsY();
	if(boxSize>=min(Nx,Ny) || gridSize>=min(Nx,Ny) ){
		cerr<<"ERROR: Box/grid size are too large compared to image size ("<<Nx<<","<<Ny<<")"<<endl;
		inputImg->Delete();
		return -1;
	}

	//Compute bkg & noise maps
	Caesar::BkgData* bkgData= inputImg->ComputeBkg(estimator,true,boxSize,boxSize,gridSize,gridSize);
	if(!bkgData){
		cerr<<"ERROR: Failed to compute bkg data!"<<endl;
		inputImg->Delete();
		return -1;
	}

	//=======================
	//== Output 
	//=======================
	Caesar::Img* BkgMap= bkgData->BkgMap;
	Caesar::Img* NoiseMap= bkgData->NoiseMap;
	Caesar::Img* SignificanceMap= 0;
	if(computeSignificance) SignificanceMap= inputImg->GetSignificanceMap(bkgData,true);
	
	if(writeToFITS){
		BkgMap->WriteFITS(outputFileName_bkg);
		NoiseMap->WriteFITS(outputFileName_noise);
		if(computeSignificance && SignificanceMap) SignificanceMap->WriteFITS(outputFileName_significance);
	}
	else{
		outputFile->cd();
		BkgMap->Write();	
		NoiseMap->Write();
		if(computeSignificance && SignificanceMap) SignificanceMap->Write();
		outputFile->Close();
	}

	//Clear
	if(inputImg) inputImg->Delete();
	if(bkgData){ 
		delete bkgData;
		bkgData= 0;
	}

	cout<<"INFO: End background finder"<<endl;
	
	return 0;

}//close main



