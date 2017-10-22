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
	cout<<"-i, --input \t Input file name containing sources to be read in ROOT TTree (.root)"<<endl;
	cout<<"-I, --input_rec \t Input file name containing reconstructed/detected sources to be read in ROOT TTree (.root)"<<endl;
	cout<<"-o, --output \t Output file name "<<endl;
	cout<<"-t, --threshold \t Fraction of matching pixels to consider sources equal "<<endl;
	cout<<"-v, --verbosity \t Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG)"<<endl;
	
	cout<<"=============================="<<endl;
}//close Usage()

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "input", required_argument, 0, 'i' },
	{ "input_rec", required_argument, 0, 'I' },
	{ "verbosity", required_argument, 0, 'v'},
	{ "output", optional_argument, 0, 'o' },
	{ "threshold", required_argument, 0, 't'},	
  {(char*)0, (int)0, (int*)0, (int)0}
};


//Options
std::string fileName= "";
std::string fileName_rec= "";
int verbosity= 4;//INFO level
float matchingThreshold= 0.9;

//Globar vars
TFile* outputFile= 0;
std::string outputFileName= "MatchOutput.root";
TTree* matchedSourceInfo;
std::vector<Source*> sources;
std::vector<Source*> sources_rec;
int SourceFoundFlag;
std::string SourceName;
int SourceType;
int SourceSimType;
double X0;
double Y0;
double S;
double Smax;
double S_true;
double X0_true;
double Y0_true;
int SourceType_rec;
double S_rec;
double X0_rec;
double Y0_rec;
double Smax_rec;
double MatchFraction;
int HasFitInfo;
double S_fit;
double X0_fit;
double Y0_fit;


//Functions
int ParseOptions(int argc, char *argv[]);
std::string GetStringLogLevel(int verbosity);
int CorrelateSourceCatalogs();
int ReadSourceData();
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
	//== Read source data
	//=======================
	INFO_LOG("Reading source data from given files...");
	if(ReadSourceData()<0){
		ERROR_LOG("Reading of source data failed!");
		return -1;
	}

	//=======================
	//== Correlate catalogs
	//=======================
	INFO_LOG("Correlating source catalogs...");
	if(CorrelateSourceCatalogs()<0){
		ERROR_LOG("Correlating source data failed!");
		return -1;
	}

	//=======================
	//== Save
	//=======================
	INFO_LOG("Saving data to file ...");
	Save();

	INFO_LOG("End source catalog correlator");

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

	while((c = getopt_long(argc, argv, "hi:I:o:v:t:",options_tab, &option_index)) != -1) {
    
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
			case 'I':	
			{
				fileName_rec= std::string(optarg);	
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
				matchingThreshold= atof(optarg);
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

int ReadSourceData()
{

	//Open files
	TFile* f= new TFile(fileName.c_str(),"READ");
	if(!f){
		ERROR_LOG("Failed to open file "<<fileName<<"!");
		return -1;
	}
	
	TFile* f_rec= new TFile(fileName_rec.c_str(),"READ");
	if(!f_rec){
		ERROR_LOG("Failed to open file "<<fileName_rec<<"!");
		return -1;
	}

	//Get access to source trees
	Source* aSource= 0;

	TTree* sourceTree= (TTree*)f->Get("SourceInfo");
	if(!sourceTree || sourceTree->IsZombie()){
		ERROR_LOG("Failed to get access to source tree in file "<<fileName<<"!");	
		return -1;
	}
	sourceTree->SetBranchAddress("Source",&aSource);

	TTree* sourceTree_rec= (TTree*)f_rec->Get("SourceInfo");
	if(!sourceTree_rec || sourceTree_rec->IsZombie()){
		ERROR_LOG("Failed to get access to source tree in file "<<fileName_rec<<"!");	
		return -1;
	}
	sourceTree_rec->SetBranchAddress("Source",&aSource);

	//Read sources
	INFO_LOG("Reading #"<<sourceTree->GetEntries()<<" sources in file "<<fileName<<"...");
	for(int i=0;i<sourceTree->GetEntries();i++){
		sourceTree->GetEntry(i);

		Source* source= new Source;
		*source= *aSource;
		sources.push_back(source);
	}//end loop sources

	INFO_LOG("Reading #"<<sourceTree_rec->GetEntries()<<" sources in file "<<fileName_rec<<"...");
	for(int i=0;i<sourceTree_rec->GetEntries();i++){
		sourceTree_rec->GetEntry(i);

		Source* source= new Source;
		*source= *aSource;
		sources_rec.push_back(source);
	}//end loop sources

	return 0;

}//close ReadSourceData()

int CorrelateSourceCatalogs()
{
	//Check source sizes
	if(sources.empty() || sources_rec.empty()){
		WARN_LOG("One or both source collections are empty, nothing to correlate...");
		return 0;
	}

	struct MatchingSourceInfo{	
		MatchingSourceInfo(float f,size_t index)
			: matchedPixelFraction(f), sourceIndex(index)
		{}
		double matchedPixelFraction;
		size_t sourceIndex;
	};

	INFO_LOG("Correlating ("<<sources.size()<<","<<sources_rec.size()<<") sources...");
	long int nFoundSources= 0;

	//Loop over first catalogue (e.g. simulated/true sources)
	for(size_t i=0;i<sources.size();i++){
		
		if(i%100==0) INFO_LOG("Find cross-match for source no. "<<i<<"/"<<sources.size()<<"...");

		//Init 1st collection pars
		long int NPix= sources[i]->GetNPixels();
		std::vector<MatchingSourceInfo> matched_info;	
		SourceFoundFlag= 0;
		SourceName= std::string(sources[i]->GetName());
		SourceType= sources[i]->Type;
		SourceSimType= sources[i]->SimType;
		X0= sources[i]->X0;
		Y0= sources[i]->Y0;
		S= sources[i]->GetS();
		Smax= sources[i]->GetSmax();
		S_true= S;
		X0_true= X0;
		Y0_true= Y0;
		bool hasTrueInfo= sources[i]->HasTrueInfo();
		if(hasTrueInfo){
			S_true= sources[i]->GetTrueFlux();
			sources[i]->GetTruePos(X0_true,Y0_true);
		}	

		//Init rec collection pars
		S_rec= -1;
		Smax_rec= -1;
		X0_rec= -1;
		Y0_rec= -1;
		SourceType_rec= -1;
		MatchFraction= 0.;
		HasFitInfo= 0;
		S_fit= -1;
		X0_fit= -1;
		Y0_fit= -1;

		//Loop over second catalogue (e.g. detected/reconstructed sources)
		for(size_t j=0;j<sources_rec.size();j++){
			long int NPix_rec= sources_rec[j]->GetNPixels();
			std::string SourceName_rec= std::string(sources_rec[j]->GetName());
			X0_rec= sources_rec[j]->X0;
			Y0_rec= sources_rec[j]->Y0;

			long int NMatchingPixels= sources[i]->GetNMatchingPixels(sources_rec[j]);
			double matchingPixelFraction= (double)(NMatchingPixels)/(double)(NPix);
			DEBUG_LOG("Source "<<SourceName<<" (X0="<<X0<<", Y0="<<Y0<<", N="<<NPix<<"): finding matching with source "<<SourceName_rec<<" (X0="<<X0_rec<<", Y0="<<Y0_rec<<", N="<<NPix_rec<<"), NMatchingPixels="<<NMatchingPixels<<" f="<<matchingPixelFraction<<" (t="<<matchingThreshold<<")");
			
			if(NMatchingPixels<=0 || matchingPixelFraction<matchingThreshold) continue;

			matched_info.push_back(MatchingSourceInfo(matchingPixelFraction,j));

		}//end loop 2nd collection

		//Fill source matching stats
		if(!matched_info.empty()){
			nFoundSources++;
			SourceFoundFlag= 1;	

			//Find best match in case of multiple matchings
			long int index_best= matched_info[0].sourceIndex;
			float matchFraction_best= matched_info[0].matchedPixelFraction; 
			
			if(matched_info.size()>1){
				for(size_t k=1;k<matched_info.size();k++){
					long int index= matched_info[k].sourceIndex;
					float matchFraction= matched_info[k].matchedPixelFraction; 
					if(matchFraction>matchFraction_best) {
						index_best= index;
						matchFraction_best= matchFraction;
					}
				}//end loop matchings
			}//close if

			SourceType_rec= sources_rec[index_best]->Type;
			X0_rec= sources_rec[index_best]->X0;
			Y0_rec= sources_rec[index_best]->Y0;
			S_rec= sources_rec[index_best]->GetS();
			Smax_rec= sources_rec[index_best]->GetSmax();
			MatchFraction= matchFraction_best;

			//If source has fit info loop over fitted components to find best match
			//NB: If fit is good this should increase the positional & flux accuracy
			HasFitInfo= sources_rec[index_best]->HasFitInfo();
			if(HasFitInfo){
				SourceFitPars fitPars= sources_rec[index_best]->GetFitPars();
				double dist_best= 1.e+99;
				for(int i=0;i<fitPars.GetNComponents();i++){
					double X0_fitcomp= fitPars.GetParValue(i,"x0");	
					double Y0_fitcomp= fitPars.GetParValue(i,"y0");
					double S_fitcomp= fitPars.GetParValue(i,"A");
					double dx= X0_fitcomp-X0_true;
					double dy= Y0_fitcomp-Y0_true;
					double dist= sqrt(dx*dx+dy*dy);
					if(dist<dist_best){
						dist_best= dist;
						S_fit= S_fitcomp;
						X0_fit= X0_fitcomp;
						Y0_fit= Y0_fitcomp;
					}
				}//end loop fitted components
			}//close has fit info

		}//close if has match
	
		matchedSourceInfo->Fill();

	}//end loop collection

	INFO_LOG("#"<<nFoundSources<<"/"<<sources.size()<<" sources found...");

	return 0;

}//close CorrelateSourceCatalogs()


void Init(){

	//Open output file
	if(!outputFile) outputFile= new TFile(outputFileName.c_str(),"RECREATE");

	if(!matchedSourceInfo) matchedSourceInfo= new TTree("data","data");
	matchedSourceInfo->Branch("found",&SourceFoundFlag,"found/I");
	matchedSourceInfo->Branch("name",&SourceName);
	matchedSourceInfo->Branch("type",&SourceType,"type/I");
	matchedSourceInfo->Branch("simtype",&SourceSimType,"simtype/I");
	matchedSourceInfo->Branch("S",&S,"S/D");
	matchedSourceInfo->Branch("Smax",&Smax,"Smax/D");	
	matchedSourceInfo->Branch("X0",&X0,"X0/D");
	matchedSourceInfo->Branch("Y0",&Y0,"Y0/D");
	matchedSourceInfo->Branch("S_true",&S_true,"S_true/D");
	matchedSourceInfo->Branch("X0_true",&X0_true,"X0_true/D");
	matchedSourceInfo->Branch("Y0_true",&Y0_true,"Y0_true/D");
	matchedSourceInfo->Branch("type_rec",&SourceType_rec,"type_rec/I");
	matchedSourceInfo->Branch("S_rec",&S_rec,"S_rec/D");
	matchedSourceInfo->Branch("Smax_rec",&Smax_rec,"Smax_rec/D");
	matchedSourceInfo->Branch("X0_rec",&X0_rec,"X0_rec/D");
	matchedSourceInfo->Branch("Y0_rec",&Y0_rec,"Y0_rec/D");
	matchedSourceInfo->Branch("MatchFraction",&MatchFraction,"MatchFraction/D");
	matchedSourceInfo->Branch("HasFitInfo",&HasFitInfo,"HasFitInfo/I");
	matchedSourceInfo->Branch("S_fit",&S_fit,"S_fit/D");
	matchedSourceInfo->Branch("X0_fit",&X0_fit,"X0_fit/D");
	matchedSourceInfo->Branch("Y0_fit",&Y0_fit,"Y0_fit/D");

}//close Init()

void Save()
{
	if(outputFile && outputFile->IsOpen()){
		outputFile->cd();		
		if(matchedSourceInfo) matchedSourceInfo->Write();
		outputFile->Close();
	}

}//close Save()


