/**
* @file SLICSegmenter.h
* @class SLICSegmenter
* @brief SLICSegmenter
*
* In this class, an over-segmentation is created of an image, provided by the
* step-size (distance between initial cluster locations) and the colour
* distance parameter.
* @author S. Riggi
* @date 15/06/2015
*/



#ifndef SLIC_SEGMENTER_H
#define SLIC_SEGMENTER_H

#include "Img.h"
#include "Region.h"
#include "ChanVeseSegmentation.h"

//#include <mathop.h>

#include <TVector2.h>
#include <TGraph.h>
#include <TText.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include <stdio.h>
#include <math.h>
#include <vector>
#include <float.h>
using namespace std;

// 2d matrices are handled by 2d vectors
#define vec2dd std::vector< std::vector<double> >
#define vec2di std::vector< std::vector<int> >
#define vec2db std::vector< std::vector<bool> >
// The number of iterations run by the clustering algorithm
#define NR_ITERATIONS 10


class SLICSegmenter {

	public:
  	/**
		\brief Class constructor
		*/
    SLICSegmenter();
		/**
		\brief Class destructor
		*/
    ~SLICSegmenter();
		
		struct SPMergingInfo {
			int levelId;
			int NR;//number of regions at current hierarchical level
			double MSE;//mean squared error (cost function)
			double DissMin;//minimum dissimilarity for this hierarchy merging
			double DissMax;//maximum dissimilarity for this hierarchy merging
			double DissMean;
			double DissRMS;
			double DissMedian;
			double DissMedianRMS;
			double DissMean0;
			double DissRMS0;
			double DissMedian0;
			double DissMedianRMS0;
		};

	public:
    
		/**
		\brief Run the image segmentation 
 		*/
		int RunSegmentation(Img* img, int regionSize=10,double regularization=100, int minRegionSize=5,bool mergeRegions=false);
	
		/**
		\brief Hierarhical merge the superpixels 
 		*/
		int SPHierarchicalMerger(std::vector<Region*>& inputRegions,std::vector< std::vector<long int> >& inputLabels,int mergerTag=1,int mergedTag=1,bool includeSpatialPars=false,int edgeModel=2);	
		
		/**
		\brief Get the pixel cluster ids
 		*/
		std::vector< std::vector<long int> > GetPixelClusterIds(){return fPixelClusterIds;}
		std::vector< std::vector<long int> > GetMergedPixelClusterIds(){return fMergedPixelClusterIds;}
		
		/**
		\brief Get the saliency map
 		*/
		Img* GetSaliencyMap(){return fSaliencyImg;}

		/**
		\brief Get the pre-merging segmented map
 		*/
		Img* GetPreMergingSegmentedMap(){return fPreMergingSegmentedImg;}

		/**
		\brief Get the list of segmented regions
 		*/
		std::vector<Region*> GetRegions() {return fRegions;}
		std::vector<Region*> GetMergedRegions() {return fMergedRegions;}

		/**
		\brief Set log contrast mapping on/off 
 		*/
		void SetLogContrastMapping(bool choice){fUseLogNormalizedImage=choice;}
		/**
		\brief Set Superpixel merging ratio parameter
 		*/
		void SetSPMergingRatio(double value){fSPMergingRatio=value;}
		/**
		\brief Set Superpixel similarity regularization
 		*/
		void SetSPMergingRegularization(double value){fSPMergingRegularization=value;}
		/**
		\brief Set Superpixel similarity regularization
 		*/
		void SetMinMergedSP(int value){fMinMergedSP=value;}
		/**
		\brief Set Superpixel similarity distance threshold
 		*/
		void SetSPMergingDistThreshold(double value){fSPMergingDistThreshold=value;}
		
		/**
		\brief Set on/off usage of pixel ratio cut
 		*/
		void UsePixelRatioCut(bool choice){fUsePixelRatioCut=choice;}	
		/**
		\brief Set Superpixel similarity distance threshold
 		*/
		void SetPixelRatioCut(double value){fPixelRatioCut=value;}

		/**
		\brief Set Superpixel similarity distance threshold
 		*/
		void SetSignificantSPRatioCut(double value){fSignificantSPRatio=value;}

		/**
		\brief Use 2nd-order neighbors in pixel merging
 		*/
		void Use2ndNeighborsInSPMerging(bool choice){fUse2ndNeighborsInSPMerging=choice;}
		

		/**
		\brief Set saliency threshold factor
 		*/
		void SetSaliencyThresholdFactor(double value){fSaliencyThresholdFactor=value;}
		void SetBkgSaliencyThresholdFactor(double value){fBkgSaliencyThresholdFactor=value;}
		void SetSaliencyImgThresholdFactor(double value){fSaliencyImgThresholdFactor=value;}
		
		/**
		\brief Set adaptive threshold threshold scale factor
 		*/
		void SetAdaptiveThresholdScale(double value){fSPMergingAdaptingDistThresholdScale=value;}
		/**
		\brief Set on/off usage of pixel ratio cut
 		*/
		void UseAdaptiveDistThreshold(bool choice){fSPMergingUseAdaptingDistThreshold=choice;}	

		/**
		\brief Set adaptive threshold threshold scale factor
 		*/
		void SetMaxDissRatio(double value){fSPMergingMaxDissRatio=value;}
		void SetMaxDissRatioFor2ndNeighbors(double value){fSPMergingMaxDissRatio_2ndNeighbor=value;}
		
		/**
		\brief Set merge model
 		*/
		void SetEdgeModel(int value){fSPMergingEdgeModel=value;}
		

		/**
		\brief Set saliency min/max/step reso
 		*/
		void SetSaliencyResoPars(int v1,int v2,int v3){
			fSaliencyMinReso=v1;
			fSaliencyMaxReso=v2;
			fSaliencyResoStepSize=v3;
		}
		
		/**
		\brief Use robust pars in saliency 
 		*/
		void UseRobustParsInSaliency(bool choice){fSaliencyUseRobustPars=choice;}	
		/**
		\brief Use bkg map in saliency
 		*/
		void UseBkgMapInSaliency(bool choice){fSaliencyUseBkgMap=choice;}	
		/**
		\brief Use noise map in saliency
 		*/
		void UseNoiseMapInSaliency(bool choice){fSaliencyUseNoiseMap=choice;}	
		/**
		\brief Use curvature map in saliency
 		*/
		void UseCurvatureMapInSaliency(bool choice){fSaliencyUseCurvatureMap=choice;}	
		/**
		\brief Set fraction of nearest-neighbors to be used in saliency map (1=all neighbors)
 		*/
		void SetNNFactorInSaliency(double value){fSaliencyNNFactor=value;}
		/**
		\brief Set filter threshold in saliency map
 		*/
		void SetFilterThresholdFactorInSaliency(double value){fSaliencyFilterThresholdFactor=value;}
		/**
		\brief Set filter threshold in saliency map
 		*/
		void SetSaliencyNormalizationMode(int value){fSaliencyNormalizationMode=value;}
		/**
		\brief Set min number of pixels in saliency map thresholding
 		*/
		void SetMinNPixInSaliencyThresholding(int value){fNPixMin=value;}

		/**
		\brief Set saliency dissimilarity exponential falloff par
 		*/
		void SetDissExpFalloffParInSaliency(double value){fSaliencyDissExpFalloffPar= value;}
		/**
		\brief Set saliency spatial dist regularization par
 		*/
		void SetSpatialDistRegParInSaliency(double value){fSaliencySpatialDistRegPar= value;}
		

		/**
		\brief Set box size in bkg map computation
 		*/
		void SetBoxSizeInLocalBkgMap(double v1,double v2){fBoxSizeX=v1;fBoxSizeY=v2;}
		/**
		\brief Set grid size in bkg map computation
 		*/
		void SetGridSizeInLocalBkgMap(double v1,double v2){fGridSizeX=v1;fGridSizeY=v2;}
		
		

		/**
		\brief Set SP merging agglo method
 		*/
		void SetSPMergingAggloMethod(int choice){fSPMergingAggloMethod=choice;}	
		/**
		\brief Set SP merging min cluster size
 		*/
		void SetSPMergingMinClustSize(int choice){fSPMergingMinClustSize=choice;}	
		/**
		\brief Set SP merging max height quantile (i.e. 0.95)
 		*/
		void SetSPMergingMaxHeightQ(double choice){fSPMergingMaxHeightQ=choice;}	
		/**
		\brief Set SP merging deep split level
 		*/
		void SetSPMergingDeepSplitLevel(int choice){fSPMergingDeepSplitLevel=choice;}	

		/**
		\brief Use curvature param in SP merging
 		*/
		void UseCurvatureInSPMerging(bool choice){fSPMergingUseCurvature=choice;}

		/**
		\brief Get merging info
 		*/
		std::vector<SPMergingInfo> GetMergingInfo(){return fSPMergingInfo;}

		/**
		\brief Set CV time step par
 		*/
		void SetCVTimeStep(double value){fCVTimeStep=value;}
		/**
		\brief Set CV window size par
 		*/
		void SetCVWindowSize(double value){fCVWindowSize=value;}
		/**
		\brief Set CV lambda 1 par
 		*/
		void SetCVLambda1Par(double value){fCVLambda1Par=value;}
		void SetCVLambda2Par(double value){fCVLambda2Par=value;}
		/**
		\brief Set CV mu par
 		*/
		void SetCVMuPar(double value){fCVMuPar=value;}
		/**
		\brief Set CV nu par
 		*/
		void SetCVNuPar(double value){fCVNuPar=value;}
		/**
		\brief Set CV p par
 		*/
		void SetCVPPar(double value){fCVPPar=value;}
		
			
		
	private:
		
		/**
  	\brief Initialize image and class data
 		*/
		int Init(Img*);
		/**
		\brief Compute the over-segmentation based on the step-size and relative weighting of the pixel and colour values.
 		*/
		int SuperpixelGenerator(int regionSize=10,double regularization=100, int minRegionSize=5);
		
		/**
		\brief Merge the superpixels 
 		*/
		int SPMaxSimilarityMerger(std::vector<Region*>& inputRegions,std::vector< std::vector<long int> >& inputLabels,int mergerTag=1,int mergedTag=1);		
		
		/**
		\brief Multi-step merger
 		*/
		int MultiStepSPMerger();
		

	private:
  	
		// The cluster assignments for each pixel
		std::vector< std::vector<long int> > fPixelClusterIds;
		std::vector< std::vector<long int> > fMergedPixelClusterIds;
		    
    // The LAB and xy values of the centers
		std::vector<Region*> fRegions;
		std::vector<Region*> fMergedRegions;
		    
		//The color lookup table
		TH1D* fLUT;

		//The curvature image
		Img* fLaplImg;
		
		//The edge-filtered image
		Img* fEdgeFilterImg;
		Img* fKirschEdgeFilterImg;
		
		//Image
		Img* fInputImg;
		Img* fImg;

		//Saliency image
		Img* fSaliencyImg;
		Img* fPreMergingSegmentedImg;
		
		//Merging info
		std::vector<SPMergingInfo> fSPMergingInfo;
		double fAbsDissMedian;
		double fAbsDissMedianRMS;

		//Options	
		bool fUseLogNormalizedImage;

		//==> Superpixel generation
		int fSPSize;
		double fSPRegularization;
		int fSPMinArea;

		//==> Superpixel Merging
		double fSPMergingRatio;
		bool fSPMergingIncludeSpatialPars;
		double fSPMergingRegularization;
		int fMinMergedSP;
		double fSPMergingDistThreshold;
		double fSPMergingMaxDissRatio;
		double fSPMergingMaxDissRatio_2ndNeighbor;
		bool fSPMergingUseCurvature;
		bool fUsePixelRatioCut;
		double fPixelRatioCut;
		bool fTagSignificativeSP;
		int fSPTaggingMethod;
		double fSignificantSPRatio;
		bool fSPMergingUseAdaptingDistThreshold;
		double fSPMergingAdaptingDistThresholdScale;

		int fSPMergingEdgeModel;
		bool fUse2ndNeighborsInSPMerging;

		//==> R Hierarchical clustering 
		int fSPMergingAggloMethod;
		int fSPMergingMinClustSize;
		double fSPMergingMaxHeightQ;
		int fSPMergingDeepSplitLevel;
 
		//==> Saliency options
		int fSaliencyMinReso;
		int fSaliencyMaxReso;
		int fSaliencyResoStepSize;
		double fSaliencyThresholdFactor;
		double fBkgSaliencyThresholdFactor;
		double fSaliencyImgThresholdFactor;
		bool fSaliencyUseRobustPars;
		bool fSaliencyUseBkgMap;
		bool fSaliencyUseNoiseMap;
		bool fSaliencyUseCurvatureMap;
		double fSaliencyNNFactor;
		double fSaliencyFilterThresholdFactor;
		int fSaliencyNormalizationMode;
		double fSaliencyDissExpFalloffPar;
		double fSaliencySpatialDistRegPar;

		int fNPixMin;
	
		//==> Bkg map options
		double fBoxSizeX;
		double fBoxSizeY;
		double fGridSizeX;
		double fGridSizeY;

		//==> Chan-Vese options
		ChanVeseSegmentation* fCVSegmentation;
		double fCVTimeStep;
		double fCVWindowSize;
		double fCVLambda1Par;
		double fCVLambda2Par;
		double fCVMuPar;
		double fCVNuPar;
		double fCVPPar;
		
};

#endif

