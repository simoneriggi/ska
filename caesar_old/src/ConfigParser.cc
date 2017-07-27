/**
* @file ConfigParser.cc
* @class ConfigParser
* @brief Parse the configuration file containing program parameters
* 
* @author S. Riggi
* @date 25/04/2010
*/

#include "ConfigParser.h"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

using namespace std;

//MAIN OPTIONS
std::string ConfigParser::fInputFileName;
std::string ConfigParser::fOutputFileName;
std::string ConfigParser::fDS9CatalogFileName;
bool ConfigParser::fSaveToFile;
bool ConfigParser::fSaveImageToFile;
int ConfigParser::fSaveImageType;
int ConfigParser::fDS9RegionFormat;
bool ConfigParser::fDrawSources;
bool ConfigParser::fIsInteractive;

//BACKGROUND OPTIONS
bool ConfigParser::fUseLocalBkg;
int ConfigParser::fLocalBkgMethod;
int ConfigParser::fBkgEstimator;
double ConfigParser::fBoxSize;
double ConfigParser::fGridSize;
double ConfigParser::fBoxSizeX;
double ConfigParser::fBoxSizeY;
double ConfigParser::fGridSizeX;
double ConfigParser::fGridSizeY;

	
// SOURCE FINDING OPTIONS
bool ConfigParser::fDeblendSources;
double ConfigParser::fCurvatureThreshold;
double ConfigParser::fPeakThreshold;
int ConfigParser::fSourceComponentMinNPix;
bool ConfigParser::fSearchNestedSources;
bool ConfigParser::fUseCurvatureMixture;
double ConfigParser::fCurvatureWeight;
int ConfigParser::fNPixMin;
double ConfigParser::fSeedBrightThreshold;
double ConfigParser::fSeedThreshold;
double ConfigParser::fMergeThreshold;
bool ConfigParser::fSearchBrightSources;
bool ConfigParser::fSearchExtendedSources;
bool ConfigParser::fSearchFaintSources;
int ConfigParser::fExtendedSearchMethod;
bool ConfigParser::fSearchNegativeExcess;
int ConfigParser::fWTScaleForFaintSourceSearch;
int ConfigParser::fWTScaleForExtendedSourceSearch;
bool ConfigParser::fUseResidualImageInExtendedSearch;

//SOURCE SELECTION
bool ConfigParser::fApplySourceSelection;
bool ConfigParser::fTagPointSources;
double ConfigParser::fMinBoundingBox;
double ConfigParser::fPointSourceCircRatioThr;
double ConfigParser::fPointSourceElongThr;
double ConfigParser::fPointSourceMinEllipseAreaRatio;
double ConfigParser::fPointSourceMaxEllipseAreaRatio;
int ConfigParser::fPointSourceMaxNPix;

//SOURCE RESIDUAL FINDING OPTIONS
int ConfigParser::fSourceDilateKernelSize;
bool ConfigParser::fDilateNestedSources;
int ConfigParser::fDilatedSourceType;
int ConfigParser::fDilateSourceModel;
bool ConfigParser::fRandomizeInDilate;
double ConfigParser::fRandSigmaInDilate;

//SEGMENTATION OPTIONS
int ConfigParser::fSPSize;
double ConfigParser::fSPRegularization;
int ConfigParser::fSPMinArea;
double ConfigParser::fSPMergingDistEps;
int ConfigParser::fSPMergingAlgo;
double ConfigParser::fSPMergingRatio;
double ConfigParser::fSPMergingRegularization;
bool ConfigParser::fUse2ndNeighborsInSPMerging;
int ConfigParser::fMinMergedSP;
double ConfigParser::fSPMergingDistThreshold;
bool ConfigParser::fSPUseLogContrast;
bool ConfigParser::fUsePixelRatioCut;
double ConfigParser::fPixelRatioCut;
bool ConfigParser::fTagSignificativeSP;
int ConfigParser::fSPTaggingMethod;
double ConfigParser::fSignificantSPRatio;
bool ConfigParser::fSPMergingUseAdaptingDistThreshold;
double ConfigParser::fSPMergingAdaptingDistThresholdScale;
double ConfigParser::fSPMergingMaxDissRatio;
double ConfigParser::fSPMergingMaxDissRatio2ndNeighbor;
bool ConfigParser::fUseCurvatureInSPMerging;
int ConfigParser::fSPMergingAggloMethod;
int	ConfigParser::fSPMergingMinClustSize;
double ConfigParser::fSPMergingMaxHeightQ;
int ConfigParser::fSPMergingDeepSplitLevel;
int ConfigParser::fSPMergingEdgeModel;

//SALIENCY
double ConfigParser::fSaliencyThresholdFactor;
double ConfigParser::fBkgSaliencyThresholdFactor;
double ConfigParser::fSaliencyImgThresholdFactor;
int ConfigParser::fSaliencyMinReso;
int ConfigParser::fSaliencyMaxReso;
int ConfigParser::fSaliencyResoStepSize;
bool ConfigParser::fSaliencyUseRobustPars;
bool ConfigParser::fSaliencyUseBkgMap;
bool ConfigParser::fSaliencyUseNoiseMap;
bool ConfigParser::fSaliencyUseCurvatureMap;
double ConfigParser::fSaliencyNNFactor;
double ConfigParser::fSaliencyFilterThresholdFactor;
int ConfigParser::fSaliencyNormalizationMode;

//CHAN-VESE OPTIONS
double ConfigParser::fCVTimeStep;
double ConfigParser::fCVWindowSize;
double ConfigParser::fCVLambda1Par;
double ConfigParser::fCVLambda2Par;
double ConfigParser::fCVMuPar;
double ConfigParser::fCVNuPar;
double ConfigParser::fCVPPar;

//SMOOTHING
bool ConfigParser::fUsePreSmoothing;
int ConfigParser::fSmoothingAlgo;
int ConfigParser::fSmoothKernelSize;
double ConfigParser::fSmoothSigma;
int ConfigParser::fGuidedSmoothRadius;
double ConfigParser::fGuidedSmoothColorEps;
		

ConfigParser::ConfigParser(std::string filename){

	fConfigFileName= filename;

	//MAIN OPTIONS
	fInputFileName= "";
	fOutputFileName= "Output.root";
	fDS9CatalogFileName= "DS9SourceCatalog.reg";
	fSaveToFile= true;
	fSaveImageToFile= false;
	fSaveImageType= 2;//Residual image
	fDS9RegionFormat= 1;	
	fDrawSources= true;
	fIsInteractive= true;

	//BACKGROUND OPTIONS
	fUseLocalBkg= false;
	fLocalBkgMethod= Img::eGridBkg;
	fBkgEstimator= Img::eMedianBkg;
	fBoxSize= 105;
	fGridSize= 21;
	fBoxSizeX= 105;
	fBoxSizeY= 105;
	fGridSizeX= 21;
	fGridSizeY= 21;

	// SOURCE FINDING OPTIONS
	fDeblendSources= true;
	fCurvatureThreshold= 0;
	fPeakThreshold= 3;
	fSourceComponentMinNPix= 6;
	fSearchNestedSources= true;
	fUseCurvatureMixture= false;
	fCurvatureWeight= 0.7;
	fNPixMin= 8;
	fSeedBrightThreshold= 20;
	fSeedThreshold= 5;
	fMergeThreshold= 2.6;
	fSearchBrightSources= true;
	fSearchExtendedSources= true;
	fSearchFaintSources= true;
	fExtendedSearchMethod= 1;
	fSearchNegativeExcess= false;
	fWTScaleForFaintSourceSearch= 1;
	fWTScaleForExtendedSourceSearch= 4;
	fUseResidualImageInExtendedSearch= true;

	// SOURCE SELECTION
	fApplySourceSelection= true;
	fTagPointSources= true;
	fMinBoundingBox= 2;
	fPointSourceCircRatioThr= 0.4;
	fPointSourceElongThr= 0.7;
	fPointSourceMinEllipseAreaRatio= 0.6;
	fPointSourceMaxEllipseAreaRatio= 1.4;
	fPointSourceMaxNPix= 500;

	// SOURCE RESIDUAL MAP OPTIONS
	fSourceDilateKernelSize= 5;
	fDilateNestedSources= false;
	fDilatedSourceType= -1;//all sources dilated
	fDilateSourceModel= 1;//median
	fRandomizeInDilate= false;
	fRandSigmaInDilate= 1;

	//SEGMENTATION OPTIONS
	fSPSize= 10;
	fSPRegularization= 100;
	fSPMinArea= 5;
	fSPMergingDistEps= 5;
	fSPMergingAlgo= 2;//hierarchical algo as default
	fSPMergingRatio= 0.3;
	fSPMergingRegularization= 0.1;
	fUse2ndNeighborsInSPMerging= true;
	fMinMergedSP= 1;
	fSPMergingMaxDissRatio= 1.15;
	fSPMergingMaxDissRatio2ndNeighbor= 1.15;
	fSPMergingDistThreshold= 0.1;
	fSPUseLogContrast= false;
	fUsePixelRatioCut= false;
	fPixelRatioCut= 0.5;
	fTagSignificativeSP= false;
	fSPTaggingMethod= 1;
	fSignificantSPRatio= 0.5;
	fSPMergingEdgeModel= 1;	
	fUseCurvatureInSPMerging= true;
	fSPMergingUseAdaptingDistThreshold= false;
	fSPMergingAdaptingDistThresholdScale= 500;
	fSPMergingAggloMethod= 2;
	fSPMergingMinClustSize= 3;
	fSPMergingMaxHeightQ= 0.95;
	fSPMergingDeepSplitLevel= 1;

	//SALIENCY
	fSaliencyThresholdFactor= 2;
	fBkgSaliencyThresholdFactor= 1;
	fSaliencyImgThresholdFactor= 1;
	fSaliencyMinReso= 20;
	fSaliencyMaxReso= 60;
	fSaliencyResoStepSize= 10;
	fSaliencyUseRobustPars= true;
	fSaliencyUseBkgMap= true;
	fSaliencyUseNoiseMap= true;
	fSaliencyUseCurvatureMap= false;
	fSaliencyNNFactor= 0.1;
	fSaliencyFilterThresholdFactor= 0.8;
	fSaliencyNormalizationMode= 2;//adaptive thr

	//CHAN-VESE OPTIONS
	fCVTimeStep= 0.1;
	fCVWindowSize= 1;	
	fCVLambda1Par= 1;
	fCVLambda2Par= 2;
	fCVMuPar= 0.5;
	fCVNuPar= 0;
	fCVPPar= 1;

	//SMOOTHING
	fUsePreSmoothing= false;
	fSmoothingAlgo= 2;
	fSmoothKernelSize= 5;
	fSmoothSigma= 1;
	fGuidedSmoothRadius= 12;
	fGuidedSmoothColorEps= 0.04;

}//close costructor


ConfigParser::~ConfigParser(){

}//close destructor


void ConfigParser::ReadConfig() {

	// read configuration from file
	cout<<"ConfigParser::ReadConfig(): INFO: Reading and parsing file "<<fConfigFileName.c_str()<<endl;

  ifstream in;  
  in.open(fConfigFileName.c_str());
  if(!in.good()) {
    string errMsg = "ConfigParser::ReadConfig(): ERROR: Cannot read config file " + fConfigFileName
      + " \n  ********* no configuration settings!!! ********* ";
    throw std::runtime_error(errMsg);
		exit(1);
  }

	//Start parsing the config file
	char buffer[1000];//container for a full line
  std::string descriptor;//container for config descriptor 
  in.getline(buffer,1000);//get the full line
	std::string parsedline;

	
	while(std::getline(in,parsedline)) {
		
		char first_char= *(parsedline.c_str());
		
		if(first_char!='#' && first_char!='\n' && first_char!=' '){
			stringstream line(parsedline);
			stringstream line_copy(parsedline);
     
			line_copy >> descriptor;
      
      if(descriptor!="\n" && descriptor!=""){
				//get all config parameters

				//########################
				//####  RUN CONFIG
				//########################
				if(descriptor.compare("inputFile")==0){
					line >> descriptor >> fInputFileName;			
		  	}
				else if(descriptor.compare("isInteractive")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fIsInteractive= true;
					else if(thisFlagValue.compare("F")==0) fIsInteractive= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for isInteractive, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				else if(descriptor.compare("saveToFile")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fSaveToFile= true;
					else if(thisFlagValue.compare("F")==0) fSaveToFile= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for saveToFile, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				else if(descriptor.compare("saveImageToFile")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue >> fSaveImageType;			
					if(thisFlagValue.compare("T")==0) fSaveImageToFile= true;
					else if(thisFlagValue.compare("F")==0) fSaveImageToFile= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for saveImageToFile, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				else if(descriptor.compare("outputFile")==0){
					line >> descriptor >> fOutputFileName;			
		  	}
				else if(descriptor.compare("DS9CatalogFile")==0){
					line >> descriptor >> fDS9CatalogFileName >> fDS9RegionFormat;			
		  	}
				else if(descriptor.compare("drawSources")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fDrawSources= true;
					else if(thisFlagValue.compare("F")==0) fDrawSources= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for drawSources, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				
				

				//## BKG OPTIONS
				else if(descriptor.compare("useLocalBkg")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseLocalBkg= true;
					else if(thisFlagValue.compare("F")==0) fUseLocalBkg= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for UseLocalBkg, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				else if(descriptor.compare("localBkgMethod")==0){
					line >> descriptor >> fLocalBkgMethod;
				}
				else if(descriptor.compare("bkgEstimator")==0){
					line >> descriptor >> fBkgEstimator;
				}
				//else if(descriptor.compare("localBkgBoxSize")==0){
				//	line >> descriptor >> fBoxSize;
				//}
				else if(descriptor.compare("localBkgBoxSize")==0){
					line >> descriptor >> fBoxSizeX >> fBoxSizeY;
				}
				//else if(descriptor.compare("localBkgGridSize")==0){
				//	line >> descriptor >> fGridSize;
				//}
				else if(descriptor.compare("localBkgGridSize")==0){
					line >> descriptor >> fGridSizeX >> fGridSizeY;
				}
	
				//## SOURCE FINDING OPTIONS
				else if(descriptor.compare("deblendSources")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue >> fCurvatureThreshold >> fSourceComponentMinNPix;
					if(thisFlagValue.compare("T")==0) fDeblendSources= true;
					else if(thisFlagValue.compare("F")==0) fDeblendSources= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for deblendSources, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				else if(descriptor.compare("useCurvatureMixture")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue >> fCurvatureWeight;			
					if(thisFlagValue.compare("T")==0) fUseCurvatureMixture= true;
					else if(thisFlagValue.compare("F")==0) fUseCurvatureMixture= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for UseCurvatureMixture, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				else if(descriptor.compare("peakThreshold")==0){
					line >> descriptor >> fPeakThreshold;
				}
				else if(descriptor.compare("minNPix")==0){
					line >> descriptor >> fNPixMin;
				}
				
				else if(descriptor.compare("searchNestedSources")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fSearchNestedSources= true;
					else if(thisFlagValue.compare("F")==0) fSearchNestedSources= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for SearchNestedSources, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				else if(descriptor.compare("seedBrightThr")==0){
					line >> descriptor >> fSeedBrightThreshold;
				}
				else if(descriptor.compare("seedThr")==0){
					line >> descriptor >> fSeedThreshold;
				}
				else if(descriptor.compare("mergeThr")==0){
					line >> descriptor >> fMergeThreshold;
				}
				else if(descriptor.compare("searchBrightSources")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fSearchBrightSources= true;
					else if(thisFlagValue.compare("F")==0) fSearchBrightSources= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for searchBrightSources, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				else if(descriptor.compare("searchExtendedSources")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fSearchExtendedSources= true;
					else if(thisFlagValue.compare("F")==0) fSearchExtendedSources= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for SearchExtendedSources, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				else if(descriptor.compare("searchFaintSources")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fSearchFaintSources= true;
					else if(thisFlagValue.compare("F")==0) fSearchFaintSources= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for SearchFaintSources, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				else if(descriptor.compare("extendedSearchMethod")==0){
					line >> descriptor >> fExtendedSearchMethod;
				}
				else if(descriptor.compare("searchNegativeExcess")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fSearchNegativeExcess= true;
					else if(thisFlagValue.compare("F")==0) fSearchNegativeExcess= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for SearchNegativeExcess, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				else if(descriptor.compare("wtScaleFaint")==0){
					line >> descriptor >> fWTScaleForFaintSourceSearch;
				}
				else if(descriptor.compare("wtScaleExtended")==0){
					line >> descriptor >> fWTScaleForExtendedSourceSearch;
				}	
				else if(descriptor.compare("useResidualImageInExtendedSearch")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseResidualImageInExtendedSearch= true;
					else if(thisFlagValue.compare("F")==0) fUseResidualImageInExtendedSearch= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useResidualImageInExtendedSearch, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				
				//## SOURCE SELECTION
				else if(descriptor.compare("applySourceSelection")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fApplySourceSelection= true;
					else if(thisFlagValue.compare("F")==0) fApplySourceSelection= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for applySourceSelection, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}	
				else if(descriptor.compare("tagPointSources")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fTagPointSources= true;
					else if(thisFlagValue.compare("F")==0) fTagPointSources= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for tagPointSources, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}	
				else if(descriptor.compare("sourceMinBoundingBox")==0){
					line >> descriptor >> fMinBoundingBox;
				}
				else if(descriptor.compare("pointSourceCircRatioThr")==0){
					line >> descriptor >> fPointSourceCircRatioThr;
				}
				else if(descriptor.compare("pointSourceElongThr")==0){
					line >> descriptor >> fPointSourceElongThr;
				}
				else if(descriptor.compare("pointSourceEllipseAreaRatioThr")==0){
					line >> descriptor >> fPointSourceMinEllipseAreaRatio >> fPointSourceMaxEllipseAreaRatio;
				}
				else if(descriptor.compare("pointSourceMaxNPix")==0){
					line >> descriptor >> fPointSourceMaxNPix;
				}
		
				//## SOURCE RESIDUAL MAP OPTIONS
				else if(descriptor.compare("dilateNestedSources")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fDilateNestedSources= true;
					else if(thisFlagValue.compare("F")==0) fDilateNestedSources= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for dilateNestedSources, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}	
				else if(descriptor.compare("dilateKernelSize")==0){
					line >> descriptor >> fSourceDilateKernelSize;
				}
				else if(descriptor.compare("dilatedSourceType")==0){
					line >> descriptor >> fDilatedSourceType;
				}

				else if(descriptor.compare("dilateSourceModel")==0){
					line >> descriptor >> fDilateSourceModel;
				}			
				else if(descriptor.compare("dilateRandomize")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue >> fRandSigmaInDilate;			
					if(thisFlagValue.compare("T")==0) fRandomizeInDilate= true;
					else if(thisFlagValue.compare("F")==0) fRandomizeInDilate= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for dilateRandomize, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}	

				//## SEGMENTATION ALGO OPTIONS
				else if(descriptor.compare("spInitPars")==0){
					line >> descriptor >> fSPSize >> fSPRegularization >> fSPMinArea;
				}
				else if(descriptor.compare("spMergingAlgo")==0){
					line >> descriptor >> fSPMergingAlgo;
				}
				else if(descriptor.compare("spHierMergingPars")==0){
					line >> descriptor >> fMinMergedSP >> fSPMergingRatio >> fSPMergingRegularization >> fSPMergingDistThreshold;
				}
				else if(descriptor.compare("spHierMergingMaxDissRatio")==0){
					line >> descriptor >> fSPMergingMaxDissRatio >> fSPMergingMaxDissRatio2ndNeighbor;
				}
				else if(descriptor.compare("useCurvatureInSPMerging")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUseCurvatureInSPMerging= true;
					else if(thisFlagValue.compare("F")==0) fUseCurvatureInSPMerging= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useCurvatureInSPMerging, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				else if(descriptor.compare("use2ndNeighborsInSPMerging")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUse2ndNeighborsInSPMerging= true;
					else if(thisFlagValue.compare("F")==0) fUse2ndNeighborsInSPMerging= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for use2ndNeighborsInSPMerging, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				else if(descriptor.compare("useLogContrastInSPGeneration")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fSPUseLogContrast= true;
					else if(thisFlagValue.compare("F")==0) fSPUseLogContrast= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useLogContrastInSPGeneration, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				else if(descriptor.compare("usePixelRatioCut")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue >> fPixelRatioCut;			
					if(thisFlagValue.compare("T")==0) fUsePixelRatioCut= true;
					else if(thisFlagValue.compare("F")==0) fUsePixelRatioCut= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for usePixelRatioCut, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				else if(descriptor.compare("tagSignificantSP")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue >> fSPTaggingMethod >> fSignificantSPRatio;			
					if(thisFlagValue.compare("T")==0) fTagSignificativeSP= true;
					else if(thisFlagValue.compare("F")==0) fTagSignificativeSP= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for tagSignificantSP, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				

				else if(descriptor.compare("useAdaptiveDistThreshold")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue >> fSPMergingAdaptingDistThresholdScale;			
					if(thisFlagValue.compare("T")==0) fSPMergingUseAdaptingDistThreshold= true;
					else if(thisFlagValue.compare("F")==0) fSPMergingUseAdaptingDistThreshold= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for useAdaptiveDistThreshold, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				else if(descriptor.compare("spMergingEdgeModel")==0){
					line >> descriptor >> fSPMergingEdgeModel;
				}
				

				else if(descriptor.compare("spMergingAggloMethod")==0){
					line >> descriptor >> fSPMergingAggloMethod;
				}
				else if(descriptor.compare("spMergingMinClustSize")==0){
					line >> descriptor >> fSPMergingMinClustSize;
				}
				else if(descriptor.compare("spMergingMaxHeightQ")==0){
					line >> descriptor >> fSPMergingMaxHeightQ;
				}
				else if(descriptor.compare("spMergingDeepSplitLevel")==0){
					line >> descriptor >> fSPMergingDeepSplitLevel;
				}


				//## SALIENCY MAP
				else if(descriptor.compare("saliencyThrRatio")==0){
					line >> descriptor >> fSaliencyThresholdFactor >> fBkgSaliencyThresholdFactor;
				}
				else if(descriptor.compare("saliencyImgThrRatio")==0){
					line >> descriptor >> fSaliencyImgThresholdFactor;
				}
				
				else if(descriptor.compare("saliencyResoPars")==0){
					line >> descriptor >> fSaliencyMinReso >> fSaliencyMaxReso >> fSaliencyResoStepSize;
				}
				else if(descriptor.compare("saliencyUseRobustPars")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fSaliencyUseRobustPars= true;
					else if(thisFlagValue.compare("F")==0) fSaliencyUseRobustPars= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for saliencyUseRobustPars, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				else if(descriptor.compare("saliencyUseBkgMap")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fSaliencyUseBkgMap= true;
					else if(thisFlagValue.compare("F")==0) fSaliencyUseBkgMap= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for saliencyUseBkgMap, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				else if(descriptor.compare("saliencyUseNoiseMap")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fSaliencyUseNoiseMap= true;
					else if(thisFlagValue.compare("F")==0) fSaliencyUseNoiseMap= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for saliencyUseNoiseMap, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				else if(descriptor.compare("saliencyUseCurvatureMap")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fSaliencyUseCurvatureMap= true;
					else if(thisFlagValue.compare("F")==0) fSaliencyUseCurvatureMap= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for fSaliencyUseCurvatureMap, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				else if(descriptor.compare("saliencyNNFactor")==0){
					line >> descriptor >> fSaliencyNNFactor;
				}
				else if(descriptor.compare("saliencyFilterThresholdFactor")==0){
					line >> descriptor >> fSaliencyFilterThresholdFactor;
				}
				else if(descriptor.compare("saliencyNormalizationMode")==0){
					line >> descriptor >> fSaliencyNormalizationMode;
				}


				//## CHAN-VESE SEGMENTATION ALGO OPTIONS
				else if(descriptor.compare("cvTimeStepPar")==0){
					line >> descriptor >> fCVTimeStep;
				}
				else if(descriptor.compare("cvWindowSizePar")==0){
					line >> descriptor >> fCVWindowSize;
				}
				else if(descriptor.compare("cvLambdaPar")==0){
					line >> descriptor >> fCVLambda1Par >> fCVLambda2Par;
				}
				else if(descriptor.compare("cvMuPar")==0){
					line >> descriptor >> fCVMuPar;
				}
				else if(descriptor.compare("cvNuPar")==0){
					line >> descriptor >> fCVNuPar;
				}
				else if(descriptor.compare("cvPPar")==0){
					line >> descriptor >> fCVPPar;
				}
	
				//## SMOOTHING
				else if(descriptor.compare("usePreSmoothing")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fUsePreSmoothing= true;
					else if(thisFlagValue.compare("F")==0) fUsePreSmoothing= false;
					else{
						string errMsg = "ConfigParser::ReadConfig(): ERROR: Invalid setting for usePreSmoothing, use T or F!";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}
				else if(descriptor.compare("smoothingFilter")==0){
					line >> descriptor >> fSmoothingAlgo;
				}
				else if(descriptor.compare("gausSmoothFilterPars")==0){
					line >> descriptor >> fSmoothKernelSize >> fSmoothSigma;
				}
				else if(descriptor.compare("guidedFilterPars")==0){
					line >> descriptor >> fGuidedSmoothRadius >> fGuidedSmoothColorEps;
				}
	

				else{
					//config setting not defined
					line >> descriptor;
					string errMsg = "ConfigParser::ReadConfig(): ERROR: Descriptor " + descriptor
      	                + " not defined \n  ********* bad settings!!! ********* ";
    			throw std::runtime_error(errMsg);
					exit(1);
				}//close else

			}//close if descriptor
		}//close if buffer

		if (!in.good()) break;
	}//close while

	in.close();

}//close ReadConfig()



void ConfigParser::Print(){

	//print parsed information
	cout<<"################################"<<endl;
	cout<<"###     PARSED  SETTINGS   #####"<<endl;
	cout<<"################################"<<endl;
	
}//close Print()

