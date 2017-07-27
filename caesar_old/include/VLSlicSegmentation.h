/**
* @file VLSlicSegmentation.h
* @class VLSlicSegmentation
* @brief VLSlicSegmentation
*
* In this class, an over-segmentation is created of an image, provided by the
* step-size (distance between initial cluster locations) and the colour
* distance parameter.
* @author S. Riggi
* @date 15/06/2015
*/



#ifndef VL_SLIC_SEGMENTATION_H
#define VL_SLIC_SEGMENTATION_H

#include "Img.h"
#include "Region.h"

#include <mathop.h>

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


class VLSlicSegmentation {

	public:
  	/**
		\brief Class constructor
		*/
    VLSlicSegmentation();
		/**
		\brief Class destructor
		*/
    ~VLSlicSegmentation();
		
		enum SuperPixelMergeAlgoType {eDBSCAN=1,eHIER=2};
		enum SuperPixelTaggingMethod {eSaliencyThresholding=1,eMahalanobisThresholding=2};

		struct SPMergingInfo {
			int levelId;
			int NR;//number of regions at current hierarchical level
			double MSE;//mean squared error (cost function)
			double DissMin;//minimum dissimilarity for this hierarchy merging
			double DissMax;//maximum dissimilarity for this hierarchy merging
			double DissMedian;
			double DissMedianRMS;
			double DissMedian0;
			double DissMedianRMS0;
			double MergedDissMin;
			double MergedDissMax;
			std::vector<double> DissList;
			std::vector<double> MergedDissList;
		};

		struct SLICSimilarityData {
			TMatrixD* DissimilarityMatrix;
			TMatrixD* SimilarityMatrix;
			TMatrixD* AbsDissimilarityMatrix;
			std::vector<double> DList;
			double Dmin;
			double Dmax;
			double Dmedian;
			double Dmedianrms;
			~SLICSimilarityData() { 
				if(DissimilarityMatrix) DissimilarityMatrix->Delete();
				if(SimilarityMatrix) SimilarityMatrix->Delete();
				if(AbsDissimilarityMatrix) AbsDissimilarityMatrix->Delete();
				DList.clear();
			}
		};


	public:
    
		/**
		\brief Returns a graph with single pixel wide contour around the clusters.
 		*/
   	TGraph* ComputeClusterContours(Img*,std::vector< std::vector<long int> >, std::vector<Region*>);
		
		TGraph* GetClusterContours(){return fContours;}
		
		/**
		\brief Returns a new image with the pixels colored according to their cluster means
 		*/
   	Img* GetClusterColoredImage(Img* image,std::vector<Region*> regions,bool drawOnlySignificative=false,bool colorWithSaliency=false);

		/**
		\brief Get the pixel cluster ids
 		*/
		std::vector< std::vector<long int> > GetPixelClusterIds(){return fPixelClusterIds;}
		/**
		\brief Get the list of connected regions per each region
 		*/
		std::vector< std::vector<int> > GetConnectedClusterList(){return fConnectedClusterList;}
		
		/**
		\brief Get the list of segmented regions
 		*/
		std::vector<Region*> GetRegions() {return fRegions;}

		/**
		\brief Run the image segmentation 
 		*/
		int RunSegmentation(Img* img, int regionSize=10,double regularization=100, int minRegionSize=5, bool mergeRegions=false, int MergeAlgoType=2);

		/**
		\brief Find the background region
 		*/
		Region* FindBackgroundRegion(double CL=0.975,bool includeSpatialPar=false,bool includeCurvPar=true,bool useOnlyPosExcess=true,bool useRobustPars=false);
		Region* FindSimpleBackgroundRegion(int nPixThreshold=1000);


		Img* GetSegmentationInputImg(){return fImg;}
		std::vector<TText*> GetRegionTextLabels(){return fRegionTextList;}

		/**
		\brief Get the computed background region
 		*/
		Img* GetBackgroundImage(){return fBackgroundImg;}

		
		/**
  	\brief Merge regions according to a hierarchical clustering method
 		*/
		int HierarchicalMergeRegions();
		/**
  	\brief Merge regions according to a DBSCAN clustering method
 		*/
		void MergeRegions(int minNpts=1,double EpsColor=5);
		
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
		\brief Set on/off usage of 2nd-order neighbors in superpixel merging 
 		*/
		void Use2ndNeighborsInSPMerging(bool choice){fUse2ndNeighborsInSPMerging=choice;}
		/**
		\brief Set Superpixel similarity regularization
 		*/
		void SetMinMergedSP(int value){fMinMergedSP=value;}
		/**
		\brief Set Superpixel similarity distance threshold
 		*/
		void SetSPMergingDistThreshold(double value){fSPMergingDistThreshold=value;}
		/**
		\brief Set on/off usage of 2nd-order neighbors in superpixel merging 
 		*/
		void UseRobustParamsInSPMerging(bool choice){fUseRobustParamsInSPMerging=choice;}
		
		/**
		\brief Set on/off usage of pixel ratio cut
 		*/
		void UsePixelRatioCut(bool choice){fUsePixelRatioCut=choice;}	
		/**
		\brief Set Superpixel similarity distance threshold
 		*/
		void SetPixelRatioCut(double value){fPixelRatioCut=value;}

		/**
		\brief Set on/off usage of significant SP tagging
 		*/
		void TagSignificantSP(bool choice){fTagSignificativeSP=choice;}
		/**
		\brief Set significant SP tagging method
 		*/
		void SetSignificantSPTaggingMethod(int choice){fSPTaggingMethod=choice;}
		/**
		\brief Set Superpixel similarity distance threshold
 		*/
		void SetSignificantSPRatioCut(double value){fSignificantSPRatio=value;}

		/**
		\brief Set saliency threshold factor
 		*/
		void SetSaliencyThresholdFactor(double value){fSaliencyThresholdFactor=value;}

		/**
		\brief Set adaptive threshold threshold scale factor
 		*/
		void SetAdaptiveThresholdScale(double value){fSPMergingAdaptingDistThresholdScale=value;}
		/**
		\brief Set on/off usage of pixel ratio cut
 		*/
		void UseAdaptiveDistThreshold(bool choice){fSPMergingUseAdaptingDistThreshold=choice;}	

		/**
		\brief Get merging info
 		*/
		std::vector<SPMergingInfo> GetMergingInfo(){return fSPMergingInfo;}

		/**
  	\brief Get saliency map
 		*/
		Img* GetSaliencyMap(){return fSaliencyImg;}
		Img* GetSumSaliencyMap(){return fSumSaliencyImg;}
		/**
  	\brief Compute saliency map
 		*/
		Img* ComputeSaliencyMap(std::vector<Region*> regions,int knn=100);
		
		
	private:
		
		/**
  	\brief Initialize image and class data
 		*/
		int Init(Img*);

		/**
  	\brief Compute the over-segmentation based on the step-size and relative weighting of the pixel and colour values.
 		*/
		int RunSLICSegmentation(int regionSize,double regularization,int minRegionSize);
		
		
		/**
  	\brief Find region nearest-neighbour
 		*/
		std::vector<int> findNeighbors(int regionId, double epsColor);
		void findNeighbors(std::vector<Region*> regions,std::map<int,int> mapping,std::vector< std::vector<long int> > labels,bool addSecondLevelNeighbors=true,bool useRobustParams=false);
		void findNeighbors_v2(std::vector<Region*> regions,std::map<int,int> mapping,std::vector< std::vector<long int> > labels,bool addSecondLevelNeighbors=true,bool useRobustParams=false);

		/**
  	\brief Compute a similarity matrix between regions
 		*/
		//std::pair<TMatrixD,TMatrixD> ComputeRegionSimilarity(std::vector<Region*> regions,std::map<int,int> mapping,double beta=0.1);
		VLSlicSegmentation::SLICSimilarityData* ComputeRegionSimilarity(std::vector<Region*> regions,std::map<int,int> mapping,double beta=0.1);
		
		/**
  	\brief Compute Google page rank of regions
 		*/
		std::vector<double> ComputePageRank(TMatrixD W,double alpha,double tol2);
		/**
  	\brief Apply region thresholding on the basis of saliency information
 		*/
		int ApplySaliencyThresholding(std::vector<Region*> regions,double saliencyThresholdFactor=2);
		/**
  	\brief Apply region thresholding on the basis of Mahalanobis distance information
 		*/
		int ApplyMahalanobisThresholding(std::vector<Region*> regions,double CL=0.975,bool includeSpatialPar=false,bool includeCurvPar=true,bool useOnlyPosExcess=true,bool useRobustPars=false);
		/**
  	\brief Tag regions (significative, salient)
 		*/
		int TagRegions(std::vector<Region*> regions,double saliencyThresholdFactor=2,double CL=0.975,bool includeSpatialPar=false,bool includeCurvPar=true);
		
	private:
  	
		// The cluster assignments and distance values for each pixel.
		std::vector< std::vector<long int> > fPixelClusterIds;
		    
    // The LAB and xy values of the centers
		std::vector<Region*> fRegions;
		std::vector< std::vector<int> > fConnectedClusterList;//list of connected region ids per each region
		
		//The cluster contours
		TGraph* fContours;
		    
		//The color lookup table
		TH1D* fLUT;

		//The curvature image
		Img* fLaplImg;
		Img::StatsData* fLaplImgStats;

		//The edge-filtered image
		Img* fEdgeFilterImg;
		Img::StatsData* fEdgeFilterImgStats;

		//Image
		bool fUseLogNormalizedImage;
		Img* fInputImg;
		Img* fImg;
		Img::StatsData* fImgStats;
    std::vector<TText*> fRegionTextList;
		TText* fRegionText;
		Img* fBackgroundImg;

		//Merging info
		std::vector<SPMergingInfo> fSPMergingInfo;
		Img* fSaliencyImg;
		Img* fSumSaliencyImg;
		
		//Options	
		double fSPMergingRatio;
		double fSPMergingRegularization;
		bool fUse2ndNeighborsInSPMerging;
		int fMinMergedSP;
		double fSPMergingDistThreshold;
		bool fUseRobustParamsInSPMerging;
		bool fUsePixelRatioCut;
		double fPixelRatioCut;
		bool fTagSignificativeSP;
		int fSPTaggingMethod;
		double fSignificantSPRatio;
		double fSaliencyThresholdFactor;
		bool fSPMergingUseAdaptingDistThreshold;
		double fSPMergingAdaptingDistThresholdScale;
};

#endif

