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

//Caesar headers
#include <Image.h>
#include <Logger.h>
#include <SysUtils.h>
#include <FITSReader.h>

//R headers
#include <RInside.h>

//ROOT headers
#include <TFile.h>

//C++ headers
#include <iostream>
#include <vector>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <chrono>

using namespace std;
using namespace Caesar;

void Usage(char* exeName){
	cout<<"=========== USAGE ==========="<<endl;
	cout<<"Usage: "<<exeName<<" [options]"<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
  cout<<"-h, --help \t Show help message and exit"<<endl;
	cout<<"-i, --input \t Input file name containing image to be read (only .fits/.root supported)"<<endl;
	cout<<"-n, --nthreads \t Number of threads to be used for reading (-1=all available threads)"<<endl;
	cout<<"-v, --verbosity \t Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG)"<<endl;
	cout<<"=============================="<<endl;
}//close Usage()

static const struct option options_tab[] = {
  { "help", no_argument, 0, 'h' },
	{ "input", required_argument, 0, 'i' },
	{ "nthreads", required_argument, 0, 'n' },
  {(char*)0, (int)0, (int*)0, (int)0}
};


//Options
std::string inputFileName= "";
int verbosity= 4;//INFO level
int nthreads= -1;

//Globar vars
Caesar::Image* inputImg= 0;
std::string imageName= "img";
Caesar::FileInfo info;

//Functions
int ParseOptions(int argc, char *argv[]);
std::string GetStringLogLevel(int verbosity);
int ReadImage();
int ComputeStats();
void Clear();
void Save();

//=======================================
//===             MAIN               ====
//=======================================
int main(int argc, char *argv[]){

	auto t0 = chrono::steady_clock::now();
	
	//================================
	//== Parse command line options
	//================================
	auto t0_parse = chrono::steady_clock::now();
	if(ParseOptions(argc,argv)<0){
		ERROR_LOG("Failed to parse command line options!");
		Clear();
		return -1;
	}
	auto t1_parse = chrono::steady_clock::now();
	double dt_parse= chrono::duration <double, milli> (t1_parse-t0_parse).count();

	//=======================
	//== Read image
	//=======================
	auto t0_read = chrono::steady_clock::now();
	if(ReadImage()<0){
		ERROR_LOG("Failed to read image from file!");		
		Clear();
		return -1;
	}
	auto t1_read = chrono::steady_clock::now();
	double dt_read= chrono::duration <double, milli> (t1_read-t0_read).count();

	//=======================
	//== Compute stats
	//=======================
	auto t0_stats = chrono::steady_clock::now();
	if(ComputeStats()<0){
		ERROR_LOG("Failed to read image from file!");		
		Clear();
		return -1;
	}
	auto t1_stats = chrono::steady_clock::now();
	double dt_stats= chrono::duration <double, milli> (t1_stats-t0_stats).count();

	//=======================
	//== Clear data
	//=======================
	Clear();

	//=======================
	//== Print perf stats
	//=======================
	auto t1 = chrono::steady_clock::now();
	double dt= chrono::duration <double, milli> (t1-t0).count();

	INFO_LOG("===========================");
	INFO_LOG("===   PERFORMANCE INFO  ===");
	INFO_LOG("===========================");
	INFO_LOG("dt(ms)= "<<dt);
	INFO_LOG("dt_read(ms)= "<<dt_read<<" ["<<dt_read/dt*100.<<"%]");
	INFO_LOG("dt_stats(ms)= "<<dt_stats<<" ["<<dt_stats/dt*100.<<"%]");
	INFO_LOG("===========================");
	
	
	INFO_LOG("End image read");
	
	return 0;

}//close main


int ParseOptions(int argc, char *argv[])
{

	//## Check args
	if(argc<2){
		cerr<<"ERROR: Invalid number of arguments...see macro usage!"<<endl;
		Usage(argv[0]);
		exit(1);
	}
	
	//## Parse options
	int c = 0;
  int option_index = 0;
	while((c = getopt_long(argc, argv, "hi:n:",options_tab, &option_index)) != -1) {
    
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
			case 'n':	
			{
				nthreads= atoi(optarg);	
				break;	
			}
			default:
			{
      	Usage(argv[0]);	
				exit(0);
			}
    }//close switch
	}//close while

	//## Set logging level
	std::string sloglevel= GetStringLogLevel(verbosity);
	LoggerManager::Instance().CreateConsoleLogger(sloglevel,"logger","System.out");
	
	//## Set number of threads
	if(nthreads>0) SysUtils::SetOMPThreads(nthreads);

	
	return 0;

}//close ParseOptions()

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



int ReadImage(){

	// Check given input file and get info
	Caesar::FileInfo info;
	if(!Caesar::SysUtils::CheckFile(inputFileName,info,false)){
		ERROR_LOG("Invalid input file ("<<inputFileName<<") specified!");
		return -1;
	}
	std::string file_extension= info.extension;
	if(file_extension!= ".fits" && file_extension!=".root") {
		ERROR_LOG("Invalid file extension ("<<file_extension<<")...nothing to be done!");
		return -1;
	}

	//--> ROOT reading
	if(file_extension==".root"){// Read image from ROOT file
		INFO_LOG("Reading ROOT input file "<<inputFileName<<"...");
		TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
		if(!inputFile || inputFile->IsZombie()){
			ERROR_LOG("Cannot open input file "<<inputFileName<<"!");
			return -1;
		}
			
		inputImg=  (Caesar::Image*)inputFile->Get(imageName.c_str());
		if(!inputImg){
			ERROR_LOG("Cannot get image from input file "<<inputFileName<<"!");
			return -1;
		}
	}//close if

	//--> FITS reading
	if(file_extension==".fits"){// Read image from FITS file
		INFO_LOG("Reading FITS input file "<<inputFileName<<"...");
		
		inputImg= new Caesar::Image;
		inputImg->SetName(imageName);
		if(inputImg->ReadFITS(inputFileName)<0){
			ERROR_LOG("Failed to read image from input file "<<inputFileName<<"!");
			return -1;
		}
	}//close else if

	if(!inputImg){
		ERROR_LOG("Failed to read image from input file "<<inputFileName<<"!");
		return -1;
	}

	return 0;

}//close ReadImage()

int ComputeStats(){

	//## Compute stats
	INFO_LOG("Computing input image stats...");
	if(inputImg->ComputeStats(true,false,false)<0){
		ERROR_LOG("Stats computing failed!");
		return -1;
	}
	inputImg->PrintStats();	

	return 0;

}//close ComputeStats()


void Clear(){

	//Clear input image
	if(inputImg) inputImg->Delete();

}//close Clear()

