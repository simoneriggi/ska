#include <BkgFinder.h>
#include <ImgFITSReader.h>
#include <Img.h>
#include <ConfigParser.h>
#include <RInside.h>

#include <iostream>
#include <vector>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

using namespace std;

void Usage(){
	cout<<"*** USAGE ***"<<endl;
	cout<<"[EXE] --config=[CONFIG-FILE]"<<endl;
	cout<<"*************"<<endl;
}

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "config", required_argument, 0, 'c' },
  {(char*)0, (int)0, (int*)0, (int)0}
};

std::string configFileName= "";
Img* fInputImg= 0;
TFile* fOutputFile= 0;
int Init();
int ReadImage();

int main(int argc, char *argv[]){

	int c = 0;
  int option_index = 0;
	std::string inputFileName= "";

	while((c = getopt_long(argc, argv, "hc:",options_tab, &option_index)) != -1) {
    
    switch (c) {
			case 'h':
			{
      	Usage();	
				exit(0);
			}
    	case 'c':	
			{
				configFileName= std::string(optarg);	
				break;	
			}
    	default:
			{
      	Usage();	
				exit(0);
			}
    }//close switch
	}//close while
 
	//## Check config file
	if(configFileName==""){
		cerr<<"ERROR: Invalid or empty config filename, see program usage!"<<endl;
		Usage();
		exit(1);
	}

	//## Create RInside instance
	RInside* fR= new RInside;
	

	//## Read config file
	cout<<"INFO: Reading config from file "<<configFileName<<"..."<<endl;
	ConfigParser parser(configFileName);
	parser.ReadConfig();
	parser.Print();

	//## Output file init
	TFile* fOutputFile= 0;
	if(ConfigParser::fSaveToFile){
		fOutputFile= new TFile((ConfigParser::fOutputFileName).c_str(),"RECREATE");
		fOutputFile->cd();
	}

	//## Read image 
	if(ReadImage()<0){
		cerr<<"ERROR: Reading image failed!"<<endl;
		exit(1);
	}

	
	//## Background finder
	cout<<"INFO: Starting background finder ..."<<endl;
	Img* BackgroundMap= 0;
	Img* NoiseMap= 0;
	
	int status = 0;

	if(ConfigParser::fUseLocalBkg){ 		
		BkgFinder finder;
		BkgFinder::BkgMapData* bkgMapData= 0;

		if((Img::LocalBkgMethod)ConfigParser::fLocalBkgMethod==Img::eSuperpixelBkg){
			finder.SetSegmentationRegionSize(ConfigParser::fSPSize);
			finder.SetSegmentationRegularization(ConfigParser::fSPRegularization);
			finder.SetSegmentationMinRegionArea(ConfigParser::fSPMinArea);
			finder.SetSegmentationDistanceEps(ConfigParser::fSPMergingDistEps);		
			//status= finder.FindSuperpixelBkg(fInputImg, (Img::BkgMethod)ConfigParser::fBkgEstimator, ConfigParser::fBoxSize,ConfigParser::fBoxSize, ConfigParser::fGridSize, ConfigParser::fGridSize);
			bkgMapData= finder.FindSuperpixelBkg(fInputImg, (Img::BkgMethod)ConfigParser::fBkgEstimator, ConfigParser::fBoxSize,ConfigParser::fBoxSize, ConfigParser::fGridSize, ConfigParser::fGridSize);
		}
		else if((Img::LocalBkgMethod)ConfigParser::fLocalBkgMethod==Img::eGridBkg){
			//status= finder.FindGridBkg(fInputImg, (Img::BkgMethod)ConfigParser::fBkgEstimator, ConfigParser::fBoxSize,ConfigParser::fBoxSize, ConfigParser::fGridSize, ConfigParser::fGridSize);
			bkgMapData= finder.FindGridBkg(fInputImg, (Img::BkgMethod)ConfigParser::fBkgEstimator, ConfigParser::fBoxSize,ConfigParser::fBoxSize, ConfigParser::fGridSize, ConfigParser::fGridSize);
		}
		else{
			cerr<<"ERROR: Invalid local background method selected!"<<endl;
			return -1;
		}

		if(status<0){
			cerr<<"ERROR: Local background finder failed!"<<endl;
			return -1;
		}
		//BackgroundMap= finder.GetInterpolatedBackgroundMap();
		//NoiseMap= finder.GetInterpolatedRMSMap();

		BackgroundMap= bkgMapData->BkgMap;
		NoiseMap= bkgMapData->RMSMap;
		
	}//close if local bkg
	else{//compute global background
		status= fInputImg->ComputeBkg((Img::BkgMethod)ConfigParser::fBkgEstimator);
		if(status<0) {
			cerr<<"ERROR: Global background computing failed!"<<endl;
			return -1;
		}
		fInputImg->DumpBkg();
	}
	
				
	//## Save to file
	if(ConfigParser::fSaveToFile && fOutputFile){
		fOutputFile->cd();	
		if(BackgroundMap) {
			BackgroundMap->SetNameTitle("BkgMap","BkgMap");
			BackgroundMap->Write();
		}
		if(NoiseMap) {
			NoiseMap->SetNameTitle("NoiseMap","NoiseMap");
			NoiseMap->Write();
		}
		fOutputFile->Close();
	}

	cout<<"INFO: End background finder"<<endl;
	
	return 0;

}//close main




int ReadImage(){

	//Get input file name
	std::string fInputFileName= ConfigParser::fInputFileName;
	if (fInputFileName == "") {
		cerr<<"ReadImage(): ERROR: Empty input file name!"<<endl;
  	return -1;
  }
	std::string fInputFileExtension= fInputFileName.substr(fInputFileName.find_last_of(".") + 1);
	if(fInputFileExtension!= "fits" && fInputFileExtension!="root") {
		cerr<<"ReadImage(): ERROR: Invalid file extension ("<<fInputFileExtension<<")...nothing to be done!"<<endl;
		return -1;
	}

	//--> ROOT reading
	if(fInputFileExtension=="root"){// Read image from ROOT file
		TFile* inputFile = new TFile(fInputFileName.c_str(),"READ");
		if(!inputFile || inputFile->IsZombie()){
			cerr<<"ReadImage(): ERROR: Cannot open input file "<<fInputFileName<<"!"<<endl;
			return -1;
		}
		fInputImg=  (Img*)inputFile->Get("FullImage");
		if(!fInputImg){
			cerr<<"ReadImage(): ERROR: Cannot get image from input file "<<fInputFileName<<"!"<<endl;
			return -1;
		}
	}//close if

	//--> FITS reading
	else if(fInputFileExtension=="fits"){// Read image from FITS file
		ImgFITSReader reader(fInputFileName.c_str());	
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
	
		if(fInputImg) fInputImg->Delete();
		fInputImg= new Img;
		reader.Read(*fInputImg);
	}//close else if

	//--> Invalid extension
	else{
		cerr<<"ReadImage(): ERROR: Invalid file extension detected!"<<endl;
		return -1;
	}
	fInputImg->SetNameTitle("img","img");	

	//## Compute stats
	cout<<"ReadImage(): INFO: Computing input image stats..."<<endl;
	if(!fInputImg->ComputeStats(true,false,true)<0){
		cerr<<"ReadImage(): ERROR: Stats computing failed!"<<endl;
		return -1;
	}
	fInputImg->DumpStats();		

	return 0;

}//close ReadImage()


