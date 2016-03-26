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
#include <FITSReader.h>
#include <SysUtils.h>
#include <Img.h>

#include <TFITS.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TMath.h>
#include <TF1.h>
#include <TF2.h>
#include <TF12.h>
#include <TROOT.h>
#include <TLine.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TMinuit.h>
#include <TMultiDimFit.h>
#include <TCut.h>
#include <TEntryList.h>
#include <TVector3.h>
#include <TPaveStats.h>
#include <TProfile.h>
#include <TApplication.h>
#include <TMatrixD.h>
#include <TRandom.h>
#include <TRandom3.h>

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <vector>
#include <stdexcept>

#include <vector>

#include <getopt.h>

using namespace std;

std::string FullFileName;
TStyle* myStyle;

void Usage(char* exeName){
	cout<<"=========== USAGE ==========="<<endl;
	cout<<"Usage: "<<exeName<<" [options]"<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
  cout<<"-h, --help \t Show help message and exit"<<endl;
	cout<<"-i, --input \t Input file name containing image to be read (.fits)"<<endl;
	cout<<"-x, --xmin \t Minimum x pixel id to be read"<<endl;
	cout<<"-X, --xmax \t Maximum x pixel id to be read"<<endl;
	cout<<"-y, --ymin \t Minimum y pixel id to be read"<<endl;
	cout<<"-Y, --ymax \t Maximum y pixel id to be read"<<endl;
	cout<<"-o, --output \t Output file name (ROOT format)"<<endl;
	cout<<"-I, --imgname \t Image name in input ROOT file (if non standard)"<<endl;
	cout<<"=============================="<<endl;
}//close Usage()

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "input", required_argument, 0, 'i' },
	{ "xmin", optional_argument, 0, 'x' },
	{ "xmax", optional_argument, 0, 'X' },
	{ "ymin", optional_argument, 0, 'y' },
	{ "ymax", optional_argument, 0, 'Y' },
	{ "imgname", optional_argument, 0, 'I' },
	{ "output", optional_argument, 0, 'o' },
  {(char*)0, (int)0, (int*)0, (int)0}
};

void SetGraphicsStyle();


int main(int argc, char **argv){

	if(argc<2){
		cerr<<"ERROR: Invalid number of arguments...see macro usage!"<<endl;
		Usage(argv[0]);
		exit(1);
	}

	//## Get command args
	int c = 0;
  int option_index = 0;
	std::string inputFileName= "";
	std::string outputFileName= "";
	bool readFullImage= true;
	bool useDefaultOutput= true;
	int minx= -1;
	int maxx= -1;
	int miny= -1;
	int maxy= -1;
	std::string imageName= "img";
	
	while((c = getopt_long(argc, argv, "hi:x::X::y::Y::o::I::",options_tab, &option_index)) != -1) {
    
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
			case 'x':	
			{
				minx= atof(optarg);	
				readFullImage= false;
				break;	
			}
			case 'X':	
			{
				maxx= atof(optarg);	
				readFullImage= false;
				break;	
			}
			case 'y':	
			{
				miny= atof(optarg);	
				readFullImage= false;
				break;	
			}
			case 'Y':	
			{
				maxy= atof(optarg);	
				readFullImage= false;
				break;	
			}
			case 'o':	
			{
				outputFileName= std::string(optarg);	
				useDefaultOutput= false;
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

	//## Check coords range in case 
	if(!readFullImage) {
		if(minx>=maxx || miny>=maxy){
			cerr<<"ERROR: Invalid coord range selected (x["<<minx<<","<<maxx<<"] y["<<miny<<","<<maxy<<"])"<<endl;
			exit(1);
		}
	}

	//## Check given input file and get info
	Caesar::FileInfo info;
	if(!Caesar::SysUtils::CheckFile(inputFileName,info,true,".fits")){
		cerr<<"ERROR: Invalid input file ("<<inputFileName<<") specified!"<<endl;
		exit(1);
	}
	if(useDefaultOutput){
		std::string basefilename_wext= info.filename_wext;
		outputFileName= basefilename_wext + std::string(".root");
	}
	
	//## Print options
	cout<<"*******************************"<<endl;
	cout<<"input file: "<<inputFileName<<endl;
	if(!readFullImage) cout<<"x["<<minx<<","<<maxx<<"] y["<<miny<<","<<maxy<<"]"<<endl;
	cout<<"image name: "<<imageName<<endl;
	cout<<"output file: "<<outputFileName<<endl;
	cout<<"*******************************"<<endl;
	
	
	

	//## Run program
	Caesar::Img* image= 0;
	TFile* OutputFile= 0;

	try{
		// Set graphics style
		SetGraphicsStyle();

		//Create output file
		OutputFile= new TFile(outputFileName.c_str(),"RECREATE");	
		OutputFile->cd();

		// Read image
		cout<<"INFO: Reading image from file..."<<endl;
		
		image= new Caesar::Img; 
		image->SetNameTitle(imageName.c_str(),imageName.c_str());

		int status= 0;
		if(readFullImage){
			status= image->ReadFITS(inputFileName);
		}
		else {
			status= image->ReadFITS(inputFileName,minx,maxx,miny,maxy);
		}
		
		if(status<0){
			cerr<<"ERROR: Failed to read image from FITS file!"<<endl;
			if(image) image->Delete();
			if(OutputFile && OutputFile->IsOpen()) OutputFile->Close();
			return -1;
		}

		int nX= image->GetNbinsX();
		int nY= image->GetNbinsY();
		cout<<"INFO: image info: nX="<<nX<<" nY="<<nY<<endl;
		
	

		cout<<"INFO: Writing to file..."<<endl;
		if(OutputFile && OutputFile->IsOpen()){
			if(image) image->Write();
			OutputFile->Close();
		}

	}//close try block
	catch(std::exception const & e) {
		cerr << "ERROR: Image read failed with status "<<e.what()<<endl;
		if(image) image->Delete();
		if(OutputFile && OutputFile->IsOpen()) OutputFile->Close();
		return -1;
	}

	return 0; 

}//close macro


void SetGraphicsStyle(){

	myStyle= new TStyle("myStyle","myStyle");
	myStyle->SetCanvasDefH(700); 
  myStyle->SetCanvasDefW(700); 

	myStyle->SetFrameBorderMode(0);
	myStyle->SetCanvasBorderMode(0);
  myStyle->SetPadBorderMode(0);
  myStyle->SetPadColor(0);
  myStyle->SetCanvasColor(0);
  myStyle->SetTitleFillColor(0);
  myStyle->SetTitleBorderSize(1);
  myStyle->SetStatColor(0);
  myStyle->SetStatBorderSize(1);
  myStyle->SetOptTitle(0);
  myStyle->SetOptStat(0);
  myStyle->SetOptFit(1);
	myStyle->SetOptLogx(0);
	myStyle->SetOptLogy(0);
  //myStyle->SetPalette(1,0);
  myStyle->SetTitleBorderSize(0);//border size of Title PavelLabel
  myStyle->SetTitleX(0.1f);
	myStyle->SetTitleW(0.8f);
  myStyle->SetStatY(0.975);                
  myStyle->SetStatX(0.95);                
  myStyle->SetStatW(0.2);                
  myStyle->SetStatH(0.15);                
  myStyle->SetTitleXOffset(0.8);
  myStyle->SetTitleYOffset(1.1);
  myStyle->SetMarkerStyle(8);
  myStyle->SetMarkerSize(0.4);
  myStyle->SetFuncWidth(1.);
  myStyle->SetPadTopMargin(0.1);
  myStyle->SetPadBottomMargin(0.12);
  myStyle->SetPadLeftMargin(0.15);
  myStyle->SetPadRightMargin(0.15);
  myStyle->SetTitleSize(0.06,"X");
  myStyle->SetTitleSize(0.06,"Y");
  myStyle->SetTitleSize(0.06,"Z");
  myStyle->SetTitleFont(52,"X");
  myStyle->SetTitleFont(52,"Y");
  myStyle->SetTitleFont(52,"Z");
  myStyle->SetLabelFont(42,"X");
  myStyle->SetLabelFont(42,"Y");
  myStyle->SetLabelFont(42,"Z");

  myStyle->SetErrorX(0.);
	
	gROOT->SetStyle("myStyle");

	/*
	const int ncolors= 2000;	
	const int Number= 12;
	int myTemperaturePalette[ncolors];
  double r[]    = {0.164, 0.150, 0.250, 0.450, 0.670, 0.880, 1.000, 1.000,1.000,0.970,0.850,0.650};
  double g[]    = {0.043, 0.306, 0.630, 0.853, 0.973, 1.000,1.000,0.880,0.679,0.430,0.150,0.000};
  double b[]    = {0.850,1.000,1.000,1.000,1.000,1.000,0.750,0.600,0.450,0.370,0.196,0.130};
  double stop[] = {0.,0.2,0.3,0.4,0.5,0.7,0.75,0.8,0.85,0.9,0.95,1};
  int FI = TColor::CreateGradientColorTable(Number, stop, r, g, b, ncolors);
  for (int i=0;i<ncolors;i++) {
		myTemperaturePalette[i] = FI+i;
	}
	*/

	/*
  int HotToColdPalette[ncolors];
  double r_HotToCold[Number];
  double g_HotToCold[Number];
  double b_HotToCold[Number];
  //double stop_HotToCold[] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.5,0.6,0.7,0.8,1};
	double stop_HotToCold[] = {0.,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1};

  for (int i=0;i<Number;i++) {
   r_HotToCold[i]= r[Number-i-1];
   g_HotToCold[i]= g[Number-i-1];
   b_HotToCold[i]= b[Number-i-1];
  }
  int FI_HotToCold = TColor::CreateGradientColorTable(Number, stop_HotToCold, r_HotToCold, g_HotToCold, b_HotToCold, ncolors);
  for (int i=0;i<ncolors;i++) {
  	HotToColdPalette[i] = FI_HotToCold + i;
  }
	*/


	gStyle->SetNumberContours(999);
	//gStyle->SetPalette(ncolors,myTemperaturePalette);

	//gStyle->SetPalette(53);//Black Body
	//gStyle->SetPalette(52);//Black & White
	gStyle->SetPalette(55);//Rainbow


}//close SetGraphicsStyle()





