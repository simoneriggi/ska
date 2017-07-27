/**
* @file SLICUtils.h
* @class SLICUtils
* @brief SLICUtils
*
* In this class, an over-segmentation is created of an image, provided by the
* step-size (distance between initial cluster locations) and the colour
* distance parameter.
* @author S. Riggi
* @date 15/06/2015
*/



#ifndef SLIC_UTILS_H
#define SLIC_UTILS_H

#include "Img.h"
#include "Region.h"
#include "Contour.h"

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


class SLICUtils {

	public:
  	/**
		\brief Class constructor
		*/
    SLICUtils();
		/**
		\brief Class destructor
		*/
    ~SLICUtils();
		
		typedef std::map<int,std::vector<int>> SLICBoundaryPixMap;//key: neighbour id, value: list of shared pix ids
		typedef SLICBoundaryPixMap::iterator SLICBoundaryPixMapIterator;
		typedef std::vector< std::vector<int> > SLICConnectedRegions;

		struct SLICContourData {
			Contour* contour;
			std::vector<SLICBoundaryPixMap> boundaryData;
			SLICConnectedRegions connectedRegionIds; 
			~SLICContourData() { 
				cout<<"SLICUtils::~SLICContourData(): Deleting contour..."<<endl;
				if(contour) {
					delete contour;
					contour= 0;
				}		
				cout<<"SLICUtils::~SLICContourData(): Clearing map/vectors..."<<endl;
				connectedRegionIds.clear();
				boundaryData.clear();
			}
		};

		struct SLICSimilarityData {
			TMatrixD* DissimilarityMatrix;
			TMatrixD* AdjacencyMatrix;
			TMatrixD* AbsDissimilarityMatrix;
			TMatrixD* SaliencyDissimilarityMatrix;
			TMatrixD* NeighborMatrix;
			std::vector< std::vector<int> > DissimilaritySortIndexMatrix;
			std::vector<double> saliencyList;
			double Dmin;
			double Dmax;
			double Dmedian;
			double Dmedianrms;
			double Emin;
			double Emax;
			SLICSimilarityData() {
				DissimilarityMatrix= 0;
				AdjacencyMatrix= 0;
				AbsDissimilarityMatrix= 0;
				SaliencyDissimilarityMatrix= 0;
				NeighborMatrix= 0;
			}
			~SLICSimilarityData() { 
				cout<<"SLICUtils::~SLICSimilarityData(): Deleting struct data..."<<endl;
				if(DissimilarityMatrix) DissimilarityMatrix->Delete();
				if(AdjacencyMatrix) AdjacencyMatrix->Delete();
				if(AbsDissimilarityMatrix) AbsDissimilarityMatrix->Delete();
				if(SaliencyDissimilarityMatrix) SaliencyDissimilarityMatrix->Delete();
				if(NeighborMatrix) NeighborMatrix->Delete();
				cout<<"SLICUtils::~SLICSimilarityData(): Clearing vector..."<<endl;
				DissimilaritySortIndexMatrix.clear();
				saliencyList.clear();
			}
		};

	public:
    
		static SLICContourData* ComputeBoundaryContours(Img* image,std::vector< std::vector<long int> > labels, std::vector<Region*> regions);
		static Img* GetSegmentedImage(Img* img,std::vector<Region*> regions,bool drawOnlySignificant=false,bool drawOnlySalient=false,bool normalize=false,bool binarize=false);
		static Img* GetSegmentedImage(Img* img,std::vector<Region*> regions,int selectedTag=-1,bool normalize=false,bool binarize=false);
		
		static int TagRegions(Img* image,Img* saliencyMap,std::vector<Region*> regions,double saliencyThresholdFactor=5,double saliencyBkgThresholdFactor=0.5,double CL=0.975,bool includeSpatialPar=false,bool includeCurvPar=true,int knn=100);
		static int TagRegions_last(Img* image,Img* saliencyMap,std::vector<Region*> regions,double saliencyThresholdFactor,double saliencyBkgThresholdFactor,int minNPix);


		static Img* GetSaliencyMap(Img* image,std::vector<Region*> regions,int knn);
		static Img* GetSaliencyMap(Img* image,TMatrixD* dissMatrix,std::vector<Region*> regions,int knn);

		static SLICSimilarityData* ComputeRegionSimilarity(Img* edgeImage,SLICUtils::SLICContourData* contourData, std::vector<Region*> regions,double beta=0.5,bool includeSpatialDist=false,int mergedTag=-1);
		static SLICSimilarityData* ComputeDissimilarityMatrix(Img* edgeImage,SLICUtils::SLICContourData* contourData, std::vector<Region*> regions,double beta=0.5,bool includeSpatialDist=false);

		static std::vector<TText*> GetRegionTextLabels(std::vector<Region*> regions,int selectedTag=-1);

		static int TagOccludedRegions(Img* image,std::vector< std::vector<long int> > labels,std::vector<Region*> regions,int mainRegionTag,int secondaryRegionTag);

	private:
		
		
		
	private:
  	
		
};

#endif

