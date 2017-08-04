/**
* @file ConfigParser.h
* @class ConfigParser
* @brief Parse the configuration file containing program parameters
* 
* @author S. Riggi
* @date 25/04/2010
*/
#ifndef _CONFIGPARSER_H_
#define _CONFIGPARSER_H_

#include <vector>
#include <string>

#include <Img.h>


class ConfigParser {
  

	public:
    
		/** 
		\brief Class constructor
 		*/
  	ConfigParser(std::string filename);
		/** 
		\brief Class destructor
 		*/
  	~ConfigParser();
		/** 
		\brief Read the config file, parse and set info to be used by other classes
 		*/
		void ReadConfig();
		/** 
		\brief Print parsed information
 		*/
		void Print();
	
	private:

		

	public:

		std::string fConfigFileName;
	
		//## MAIN OPTIONS
		static std::string fInputFileName;
		static std::string fOutputFileName;
		static std::string fDS9CatalogFileName;
		static std::string fROOTInputImgName;
		static bool fSaveToFile;
		static bool fSaveImageToFile;
		static int fSaveImageType;
		static int fDS9RegionFormat;
		static bool fDrawSources;
		static bool fIsInteractive;

		//## BKG OPTIONS
		static bool fUseLocalBkg;
		static int fLocalBkgMethod;
		static int fBkgEstimator;
		static double fBoxSize;
		static double fGridSize;
		static double fBoxSizeX;
		static double fBoxSizeY;
		static double fGridSizeX;
		static double fGridSizeY;
	
		//## SOURCE FINDING OPTIONS
		static bool fDeblendSources;
		static double fCurvatureThreshold;
		static double fPeakThreshold;
		static int fSourceComponentMinNPix;
		static bool fUseCurvatureMixture;
		static double fCurvatureWeight;
		static bool fSearchNestedSources;
		static int fNPixMin;
		static double fSeedBrightThreshold;
		static double fSeedThreshold;
		static double fMergeThreshold;
		static bool fSearchBrightSources;
		static bool fSearchExtendedSources;
		static bool fSearchFaintSources;
		static int fExtendedSearchMethod;	
		static bool fSearchNegativeExcess;
		static int fWTScaleForFaintSourceSearch;
		static int fWTScaleForExtendedSourceSearch;
		static bool fUseResidualImageInExtendedSearch;

		//## SOURCE SELECTION	
		static bool fApplySourceSelection;
		static double fMinBoundingBox;
		static bool fTagPointSources;
		static double fPointSourceCircRatioThr;
		static double fPointSourceElongThr;
		static double fPointSourceMinEllipseAreaRatio;
		static double fPointSourceMaxEllipseAreaRatio;
		static int fPointSourceMaxNPix;

		//## SOURCE RESIDUAL MAP OPTiONS
		static int fSourceDilateKernelSize;
		static bool fDilateNestedSources;
		static int fDilatedSourceType;
		static int fDilateSourceModel;
		static bool fRandomizeInDilate;
		static double fRandSigmaInDilate;

		//## SEGMENTATION ALGO OPTIONS
		static int fSPSize;
		static double fSPRegularization;
		static int fSPMinArea;
		static double fSPMergingDistEps;
		static int fSPMergingAlgo;
		static double fSPMergingRatio;
		static double fSPMergingRegularization;
		static bool fUse2ndNeighborsInSPMerging;
		static int fMinMergedSP;
		static double fSPMergingMaxDissRatio;
		static double fSPMergingMaxDissRatio2ndNeighbor;
		static double fSPMergingDistThreshold;
		static bool fSPMergingUseAdaptingDistThreshold;
		static double fSPMergingAdaptingDistThresholdScale;
		static double fSPMergingDist;
		static bool fSPUseLogContrast;
		static bool fUsePixelRatioCut;
		static double fPixelRatioCut;
		static bool fTagSignificativeSP;
		static int fSPTaggingMethod;
		static double fSignificantSPRatio;
		static bool fUseCurvatureInSPMerging;
		static int fSPMergingEdgeModel;

		static int fSPMergingAggloMethod;
		static int fSPMergingMinClustSize;
		static double fSPMergingMaxHeightQ;
		static int fSPMergingDeepSplitLevel;

		//## SALIENCY COMPUTATION
		static int fSaliencyMinReso;
		static int fSaliencyMaxReso;
		static int fSaliencyResoStepSize;
		static double fBkgSaliencyThresholdFactor;
		static double fSaliencyThresholdFactor;
		static double fSaliencyImgThresholdFactor;
		static bool fSaliencyUseRobustPars;
		static bool fSaliencyUseBkgMap;
		static bool fSaliencyUseNoiseMap;
		static bool fSaliencyUseCurvatureMap;
		static double fSaliencyNNFactor;
		static double fSaliencyFilterThresholdFactor;
		static int fSaliencyNormalizationMode;
	
		static double fSaliencyDissExpFalloffPar;
		static double fSaliencySpatialDistRegPar;

		//## CHAN-VESE ALGO OPTIONS
		static double fCVTimeStep;
		static double fCVWindowSize;
		static double fCVLambda1Par;
		static double fCVLambda2Par;
		static double fCVMuPar;
		static double fCVNuPar;
		static double fCVPPar;

		//## FILTER STAGE
		static bool fUsePreSmoothing;
		static int fSmoothingAlgo;
		static int fSmoothKernelSize;
		static double fSmoothSigma;
		static int fGuidedSmoothRadius;
		static double fGuidedSmoothColorEps;
		

		friend class SourceFinder;	
		friend class BkgFinder;
};


#endif
 
