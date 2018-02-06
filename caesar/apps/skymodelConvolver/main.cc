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

#include <Image.h>
#include <Source.h>

#include <ConfigParser.h>
#include <Logger.h>
#include <CodeUtils.h>

//ROOT headers
#include <TFile.h>
#include <TTree.h>

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
	cout<<"-i, --input=[INPUT_FILE] \t Input file name containing sources to be read in ROOT TTree (.root)"<<endl;
	cout<<"-o, --output=[OUTPUT_FILE] \t Output file name "<<endl;
	cout<<"-B, --bmaj=[BMAJ] \t Bmaj in arcsec (NB: Mandatory)"<<endl;
	cout<<"-b, --bmin=[BMAJ] \t Bmin in arcsec (NB: Mandatory)"<<endl;
	cout<<"-a, --bpa=[BMAJ] \t Bpa in degrees (NB: Mandatory)"<<endl;
	cout<<"-t, --threshold=[THRESHOLD] \t Flux threshold below which pixels are removed from sources (default=0)"<<endl;
	cout<<"-T, --threshold_ext=[THRESHOLD_EXT] \t Flux threshold below which pixels are removed from extended sources (default=0)"<<endl;
	cout<<"-m, --mergeOverlappingSources \t Merge overlapping thresholds after convolution (default=no)"<<endl;
	cout<<"-s, --nsigmas=[NSIGMAS] \t Number of sigmas used in convolution gaussian kernel (default=10)"<<endl;
	cout<<"-v, --verbosity=[LEVEL] \t Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG) (default=INFO)"<<endl;
	
	cout<<"=============================="<<endl;
}//close Usage()

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "verbosity", required_argument, 0, 'v'},
	{ "input", required_argument, 0, 'i' },	
	{ "output", optional_argument, 0, 'o' },
	{ "bmaj", required_argument, 0, 'B'}, //Bmaj
	{ "bmin", required_argument, 0, 'b'}, //Bmin
	{ "bpa", required_argument, 0, 'a'}, //Pos angle
	{ "threshold", required_argument, 0, 't'}, //Flux threshold below which pixels are removed from convolved sources
	{ "threshold_ext", required_argument, 0, 'T'}, //Flux threshold below which pixels are removed from convolved extended sources
	{ "mergeOverlappingSources", no_argument, 0, 'm'},
	{ "nsigmas", required_argument, 0, 's'}, //Gaus conv kernel nsigmas (default=10)
  {(char*)0, (int)0, (int*)0, (int)0}
};


//Options
int verbosity= 4;//INFO level
std::string fileName= "";
std::string outputFileName= "skymodel_conv.root";
double Bmaj= -1;
double Bmin= -1;
double Bpa= -1;
double fluxThr= -1;
double fluxThr_ext= -1;
bool mergeOverlappingSources= false;
int nSigmas= 10;
int minPixels= 5;


//Vars
Image* img= 0;
std::vector<Source*> sources;
Image* img_conv= 0;
std::vector<Source*> sources_conv;
TFile* outputFile= 0;

//Functions
int ParseOptions(int argc, char *argv[]);
std::string GetStringLogLevel(int verbosity);
int ReadData();
int RunConvolver();
void Save();
void Init();

int main(int argc, char *argv[])
{
	//================================
	//== Parse command line options
	//================================
	if(ParseOptions(argc,argv)<0){
		ERROR_LOG("Failed to parse command line options!");
		return -1;
	}
	
	//=======================
	//== Init
	//=======================
	INFO_LOG("Initializing data...");
	Init();

	//=======================
	//== Read data
	//=======================
	INFO_LOG("Reading source data ...");
	if(ReadData()<0){
		ERROR_LOG("Reading of source data failed!");
		return -1;
	}

	//=======================
	//== Convolver
	//=======================
	INFO_LOG("Convolve skymodel sources ...");
	if(RunConvolver()<0){
		ERROR_LOG("Skymodel source convolver failed!");
		return -1;
	}

	//=======================
	//== Save
	//=======================
	INFO_LOG("Saving data to file ...");
	Save();

	INFO_LOG("End skymodel convolver");

	return 0;

}//close main



int ParseOptions(int argc, char *argv[])
{
	//## Check args
	if(argc<2){
		cerr<<"ERROR: Invalid number of arguments (n="<<argc<<") ...see program usage!"<<endl;
		Usage(argv[0]);
		exit(1);
	}

	//## Parse options
	int c = 0;
  int option_index = 0;

	while((c = getopt_long(argc, argv, "hi:o:v:t:T:B:b:a:s:m",options_tab, &option_index)) != -1) {
    
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
				fileName= std::string(optarg);	
				break;	
			}
			case 'o':	
			{
				outputFileName= std::string(optarg);	
				break;	
			}
			case 'v':	
			{
				verbosity= atoi(optarg);	
				break;	
			}
			case 't':
			{
				fluxThr= atof(optarg);
				break;
			}
			case 'T':
			{
				fluxThr_ext= atof(optarg);
				break;
			}
			case 'B':
			{
				Bmaj= atof(optarg);
				break;
			}
			case 'b':
			{
				Bmin= atof(optarg);
				break;
			}
			case 'a':
			{
				Bpa= atof(optarg);
				break;
			}
			case 'm':
			{
				mergeOverlappingSources= true;
				break;
			}
			case 's':
			{
				nSigmas= atoi(optarg);
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
	//== Init Logger 
	//=======================
	//## Set logging level
	std::string sloglevel= GetStringLogLevel(verbosity);
	LoggerManager::Instance().CreateConsoleLogger(sloglevel,"logger","System.out");
	
	//=======================
	//== Check args 
	//=======================
	if(Bmin<0 || Bmaj<0 || Bpa<0){
		ERROR_LOG("Invalid beam parameters given (hint: beam pars are mandatory and shall be >0)");
		return -1;
	}
	if(fileName==""){
		ERROR_LOG("Invalid/empty input file name given!");
		return -1;
	}
	if(nSigmas<=0){
		ERROR_LOG("Invalid nsigma arg given (hint: shall be >0)!");
		return -1;
	}

	return 0;

}//close ParseOptions()


int RunConvolver()
{
	//================================================
	// ==> Convolver steps
	//1) Read skymodel image (to get image size)	
	//2) Loop over true sources
	//    - Read source i-th
	//    - Create mask image with source i-th
	//    - Convolve mask image with given beam
	//    - Apply threshold to convolved image (e.g. set to zero all pixels below a given significance)
	//    - Find source in thresholded convolved image and add to smoothed source list
	//3) Make a mask image with all smooothed sources --> write to fits & root
	//4) Merge sources sharing pixels (if enabled) --> write to root
	//=================================================

	img_conv= img->GetCloned("",true,true);
	img_conv->Reset();

	for(size_t i=0;i<sources.size();i++){
		//Get threshold
		int type= sources[i]->Type;
		double thr= fluxThr;
		if(type==Source::eCompact || type==Source::ePointLike){
			thr= fluxThr;
		}
		else if(type==Source::eExtended || type==Source::eCompactPlusExtended){
			thr= fluxThr_ext;
		}
		else{
			thr= fluxThr;
		}
	
		//Get image source mask
		bool isBinary= false;	
		bool invert= false;
		Image* sourceImg= img->GetSourceMask({sources[i]},isBinary,invert);
		if(!sourceImg){
			ERROR_LOG("Failed to compute image mask for source "<<i+1<<"!");
			return -1;
		}

		//Convolve source image with beam
		Image* sourceImg_conv= sourceImg->GetBeamConvolvedImage(Bmaj,Bmin,Bpa,nSigmas);
		if(!sourceImg_conv){
			ERROR_LOG("Failed to convolve source image "<<i+1<<"!");
			delete sourceImg;	
			sourceImg= 0;
			return -1;
		}

		//Delete source image
		delete sourceImg;
		sourceImg= 0;

		//Find convolved source
		std::vector<Source*> csources;
		if(sourceImg_conv->FindCompactSource(csources,thr,minPixels)<0){
			ERROR_LOG("Failed to find convolved source "<<i+1<<"!");
			delete sourceImg_conv;	
			sourceImg_conv= 0;
			return -1;
		}
	
		
		if(csources.size()>1){
			ERROR_LOG("More than 1 source found in convolved image (this should not occur...)!");
			delete sourceImg_conv;
			sourceImg_conv= 0;
			return -1;
		}
		if(csources.empty()){
			WARN_LOG("Source "<<i+1<<" not found after convolution (below npix/flux threshold?) will be removed...");
			delete sourceImg_conv;
			sourceImg_conv= 0;
			continue;
		}

		//Add convolved source to list
		sources_conv.push_back(csources[0]);		

		//Add convolved image to skymodel
		bool computeStats= false;
		img_conv->Add(sourceImg_conv,1,computeStats);

		//Delete convolved source image
		delete sourceImg_conv;
		sourceImg_conv= 0;

		
	}//end loop sources

	INFO_LOG("#"<<sources_conv.size()<<" sources present after convolution...");

	//Compute stats of skymodel convolved image
	INFO_LOG("Computing stats of skymodel convolved image...");
	img_conv->ComputeStats(true);

	return 0;

}//close RunConvolver()


void Init()
{
	//Open output file
	if(!outputFile) outputFile= new TFile(outputFileName.c_str(),"RECREATE");

	
}//close Init()


std::string GetStringLogLevel(int verbosity){

	std::string slevel= "";
	if(verbosity<=0) slevel= "FATAL";
	else if(verbosity==1) slevel= "FATAL";
	else if(verbosity==2) slevel= "ERROR";
	else if(verbosity==3) slevel= "WARN";
	else if(verbosity==4) slevel= "INFO";
	else if(verbosity>5) slevel= "DEBUG";
	else slevel= "OFF";

	return slevel;

}//close GetStringLogLevel()

int ReadData()
{
	//Open file with source collection
	TFile* inputFile= new TFile(fileName.c_str(),"READ");
	if(!inputFile){
		ERROR_LOG("Failed to open file "<<fileName<<"!");
		return -1;
	}

	//Read skymodel image
	img= (Image*)inputFile->Get("img");
	if(!img){
		ERROR_LOG("Failed to read skymodel image from input file "<<fileName<<"!");
		return -1;
	}
	
	//Get access to source trees
	Source* aSource= 0;

	TTree* sourceTree= (TTree*)inputFile->Get("SourceInfo");
	if(!sourceTree || sourceTree->IsZombie()){
		ERROR_LOG("Failed to get access to source tree in file "<<fileName<<"!");	
		return -1;
	}
	sourceTree->SetBranchAddress("Source",&aSource);

	
	//Read sources
	INFO_LOG("Reading #"<<sourceTree->GetEntries()<<" sources in file "<<fileName<<"...");
	for(int i=0;i<sourceTree->GetEntries();i++){
		sourceTree->GetEntry(i);

		Source* source= new Source;
		*source= *aSource;
		sources.push_back(source);
	}//end loop sources

	INFO_LOG("#"<<sources.size()<<" sources read...");

	
	return 0;

}//close ReadData()

void Save()
{
	//Save data to file
	if(outputFile && outputFile->IsOpen()){
		outputFile->cd();		
		if(img_conv) {
			img_conv->SetName("img");
			img_conv->Write(); 
		}
		outputFile->Close();
	}

}//close Save()



