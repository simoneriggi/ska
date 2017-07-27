/**
* @file SourceFinder.h
* @class SourceFinder
* @brief SourceFinder
*
* Image class
* @author S. Riggi
* @date 17/08/2015
*/

#ifndef SourceFinder_h
#define SourceFinder_h 1

#include <Img.h>
//#include <VLSlicSegmentation.h>
#include <ChanVeseSegmentation.h>
#include <SLICSegmenter.h>
#include <SLICUtils.h>

#include <RInside.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TH3D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TF1.h>
#include <TF12.h>
#include <TF2.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TStyle.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TApplication.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <iostream>
#include <time.h>
#include <ctime>

using namespace std;


class SourceFinder {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    SourceFinder();
		/** 
		\brief Class destructor
 		*/
  	~SourceFinder();

	public:

		int Run();
		
	private:
	
		enum SaveImageType {eInputImage=1,eResidualImage=2,eSegmentedImage=3};
		enum ExtendedSearchMethod {eWT=1,eSPSegmentation=2,eCVSegmentation=3};
		enum SmoothingFilter {eGaus=1,eGuided=2,eBilateral=3};

		int SetConfig();
		int Init();	
		int ReadImage();
		int FindSource(Img* img,double seedThr,double mergeThr);
		int FindBrightSource();	
		int FindFaintSource();	
		int FindExtendedSource();
		int FindExtendedSource_WT(Img*);	
		int FindExtendedSource_SPSegmentation(Img*);
		int FindExtendedSource_SPSegmentation_new(Img*);
		
		int FindExtendedSource_CVContour(Img*);
		int FindNestedSource(Source* aSource);	
		int DrawSources(Img* img,bool drawSourcePlots=false);	
		void SetGraphicsStyle();
		void Save();	
		bool IsGoodSource(Source* aSource);
		bool IsCompactPointLike(Source* aSource);
		int DeblendSources(Img* img);
		
	private:

		//## RInside& ROOT App instance
		RInside* fR;
		TApplication* fApplication;
		TStyle* myStyle;
		bool fSaveToFile;
		bool fSaveImageToFile;
		int fSaveImageType;
		bool fDrawSources;
		bool fIsInteractive;
		
		//Input data
		std::string fInputFileName;
		std::string fInputFileExtension;
		Img* fInputImg;

		//Output data
		FILE* fDS9CatalogFilePtr;
		FILE* fDS9CatalogFilePtr2;
		std::string fDS9CatalogFileName;
		int fDS9RegionFormat;
		std::string fOutputFileName;
		TFile* fOutputFile;
		TTree* fOutputTree;
		Source* fSource;
		std::vector<Source*> fSourceCollection;
		TTree* fConfigInfo;	

		//Output segmentation data
		TTree* fSPMergingInfo;
		int fLevelId;
		int fNDissEntries;	
		int fNDissEntries_merged;
		static const int MAXNENTRIES= 100000;
		double fDissList[MAXNENTRIES];
		double fDissList_merged[MAXNENTRIES];
		double fDissMin_merged;
		double fDissMax_merged;
		double fDissMin;
		double fDissMax;	
		double fDissMean;
		double fDissRMS;
		double fDissMedian;
		double fDissMedianRMS;		
		double fDissMean0;
		double fDissRMS0;
		double fDissMedian0;
		double fDissMedianRMS0;	
		int fNR;
		double fMSE;
		Img* fSegmentedImage;
		Img* fFinalSegmentedImage;
		Img* fFinalSignalSegmentedImage;
		TGraph* fSegmentationContourGraph;
		
		Img* fInitSaliencyMap;
		Img* fSaliencyMap;
		Img* fSumSaliencyMap;
		
		//Bkg options
		bool fUseLocalBackground;
		int fLocalBkgModel;
		int fBkgEstimator;
		double fBkgGridSize;
		double fBkgBoxSize;
		double fBkgGridSizeX;
		double fBkgGridSizeY;
		double fBkgBoxSizeX;
		double fBkgBoxSizeY;

		//Source options
		bool fDeblendSources;
		double fCurvatureThreshold;
		double fPeakThreshold;
		int fSourceComponentMinNPix;
		bool fUseCurvatureMixture;
		double fCurvatureWeight;
		bool fSearchNestedSources;
		int fNPixMin;
		double fSeedBrightThreshold;
		double fSeedThreshold;
		double fMergeThreshold;
		bool fSearchBrightSources;
		bool fSearchExtendedSources;
		bool fSearchFaintSources;
		int fExtendedSearchMethod;
		bool fSearchNegativeExcess;

		int fWTScaleForFaintSourceSearch;
		int fWTScaleForExtendedSourceSearch;
		bool fUseResidualImageInExtendedSearch;
		
		//Source selection
		bool fApplySourceSelection;
		double fMinBoundingBox;
		bool fTagPointSources;
		double fPointSourceCircRatioThr;
		double fPointSourceElongThr;
		double fPointSourceMinEllipseAreaRatio;
		double fPointSourceMaxEllipseAreaRatio;
		int fPointSourceMaxNPix;

		//Bright source residual options
		int fSourceDilateKernelSize;
		bool fDilateNestedSources;
		int fDilatedSourceType;
		int fDilateSourceModel;
		bool fRandomizeInDilate;
		double fRandSigmaInDilate;

		//Segmentation algo
		//VLSlicSegmentation* fSLICSegmentation;
		SLICSegmenter* fSLICSegmenter;
		int fSPSize;
		double fSPRegularization;
		int fSPMinArea;
		double fSPMergingDistEps;
		int fSPMergingAlgo;
		double fSPMergingRatio;
		double fSPMergingRegularization;
		bool fUse2ndNeighborsInSPMerging;
		int fMinMergedSP;
		double fSPMergingDistThreshold;
		bool fSPUseLogContrast;
		bool fUsePixelRatioCut;
		double fPixelRatioCut;
		bool fTagSignificativeSP;
		int fSPTaggingMethod;
		double fSignificantSPRatio;
		bool fSPMergingUseAdaptingDistThreshold;
		double fSPMergingAdaptingDistThresholdScale;		
		double fSPMergingMaxDissRatio;
		double fSPMergingMaxDissRatio2ndNeighbor;
		bool fUseCurvatureInSPMerging;
		int fSPMergingAggloMethod;
		int fSPMergingMinClustSize;
		double fSPMergingMaxHeightQ;
		int fSPMergingDeepSplitLevel;
		int fSPMergingEdgeModel;
		
		//Saliency map
		double fSaliencyThresholdFactor;
		double fBkgSaliencyThresholdFactor;
		double fSaliencyImgThresholdFactor;
		int fSaliencyMinReso;
		int fSaliencyMaxReso;
		int fSaliencyResoStepSize;
		bool fSaliencyUseRobustPars;
		bool fSaliencyUseBkgMap;
		bool fSaliencyUseNoiseMap;
		bool fSaliencyUseCurvatureMap;
		double fSaliencyNNFactor;
		double fSaliencyFilterThresholdFactor;
		int fSaliencyNormalizationMode;

		//Chan-Vese segmentation
		ChanVeseSegmentation* fCVSegmentation;
		double fCVTimeStep;
		double fCVWindowSize;
		double fCVLambda1Par;
		double fCVLambda2Par;
		double fCVMuPar;
		double fCVNuPar;
		double fCVPPar;
		
		//Smoothing
		bool fUsePreSmoothing;
		int fSmoothingAlgo;
		int fSmoothKernelSize;
		double fSmoothSigma;
		int fGuidedSmoothRadius;
		double fGuidedSmoothColorEps;
		
};

#endif


