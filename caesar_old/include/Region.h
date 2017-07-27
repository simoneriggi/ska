/**
* @file Region.h
* @class Region
* @brief Region
*
* Image region
* @author S. Riggi
* @date 22/06/2015
*/



#ifndef Region_h
#define Region_h 1


#include "Contour.h"
#include "Source.h"

#include <TH1D.h>
#include <TH2D.h>
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
#include <TGraph.h>

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

using namespace std;

class Img;

class Region {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    Region();
		//Region(const Region &obj);

		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~Region();

		enum RegionImageMode {eBinarizedMap=0,eFluxMap=1,eMeanMap=2,eRMSMap=3};
		
		enum RegionTag {eBkgTag=0,eSignalTag=1,eUntagged=2};

		struct Pixel {	
			int id;//global bin id of reference image	
			double S;//pixel intensity
			TVector3 color;//color vector	
			double x;//pixel x coordinate
			double y;//pixel y coordinate
			int ix;//pixel id x
			int iy;//pixel id y
			double S_curv;//curvature estimator
			double S_edge;//edge estimator
			bool isOnEdge;//flag marking if pixel is found on region contour
			double distanceToEdge;//distance to the edge (=0 for edge pixels)		
		};


		struct EdgeInfo{
			double x;
			double y;
			int labelId;
			double S_curv;
			double S_edge;
		};

		struct RobustStats{
			double median;
			double mad;
			double median_curv;
			double mad_curv;
			double entropy;
		};
		
		struct Parameters {
			TVectorD* pars;
			TVectorD* robustPars;
			TVectorD* spatialPars;
			Parameters() {
				pars= 0;
				robustPars= 0;
				spatialPars= 0;
			}
			~Parameters() { 
				if(pars) pars->Delete();
				if(robustPars) robustPars->Delete();
				if(spatialPars) spatialPars->Delete();	
			}
		};//close Parameters()

		struct NeighborInfo {
			int id;//neighbor region id	
			double D;//distance from this region to neighbor Dij= D(i,i+j) 
			int order;//neighbor order (1st=1,2nd=2)
			double H;//entropy of merged region E(i+j)	
			double Edgeness;//edgeness estimator for shared boundaries
			bool operator==(const NeighborInfo& m) const {
        return ( (m.id==id) && (m.D==D) );
    	}
		};
		struct FindNeighborId {
    	int id;
    	FindNeighborId(int id) : id(id) {}
    	bool operator () ( const NeighborInfo& l) const {
        return id == l.id;
    	}
		};
		typedef std::vector<Pixel> PixelCollection;
		typedef std::map<int,Pixel> PixelMap;


	public:

		//Setters & getters	
		void SetId(int id){fId=id;}
		void SetNPixels(int value){fNPix=value;}	
		std::pair<double,double> GetDistance(Region* aRegion,bool useRobustParams=false,bool normalizeParams=true,bool addCurvDist=true);
		std::pair<double,double> GetAsymmDistance(Region* aRegion,bool useRobustParams=false,bool normalizeParams=true,bool addSpatialDist=false,bool addCurvDist=true);
		double GetDissimilarity(Region* aRegion,bool useRobustParams=false,bool normalizeParams=true,bool addCurvDist=true);
		std::vector<double> GetEntropyDistance(Region* aRegion);
		TVectorD GetParamVector(bool includeSpatialPar=false,bool includeCurvPar=true,bool useRobustPars=false);
		Region::Parameters* GetParams(bool includeCurvPar=true);

		//Pixel manipulation
		PixelCollection GetPixels(){return fPixelCollection;}
		void AddPixel(Pixel pixel);
		int ComputeParameters(bool computeContours=false,bool computeRobustStats=true,bool forceRecomputing=false);

		void AddRegion(Region* aRegion,bool addPixels=true);
		void AddNeighbor(int neighborId){fNeighbourRegions.push_back(neighborId);}
		int AddNeighborInfo(Region* aNeighborRegion,int order=1,bool useRobustParams=true,bool normalizeParams=true,bool addCurvDist=true);
		void AddNeighborInfo(NeighborInfo aNeighborInfo);
		void SortNeighborInfo(bool isAscending=true,bool useTotDistance=false,double beta=0);

		int ComputeKLD(Region* aRegion,double& KLD);
		
		/**
		* \brief Dump region info
		*/
		void Dump(){
			cout<<"*** REGION NO. "<<fId<<" ***"<<endl;
			cout<<"N= "<<fNPix<<" M1="<<fM1<<", M2="<<fM2<<", M3="<<fM3<<", M4="<<fM4<<endl;
			cout<<"X0="<<fX0<<" Y0="<<fY0<<" Mean="<<fMean<<" RMS="<<fRMS<<" Skewness="<<fSkewness<<", Kurtosis="<<fKurtosis<<" Median="<<fMedian<<" MedianRMS="<<fMedianRMS<<endl;
			cout<<"Smin/Smax="<<fSmin<<"/"<<fSmax<<" Xmin/Xmax="<<fXmin<<"/"<<fXmax<<", Ymin/Ymax="<<fYmin<<"/"<<fYmax<<endl;
			//cout<<"Image RMS="<<fImageRMS<<endl;
			cout<<"****************************"<<endl;
		}

		/**
		* \brief Create and return a source from this region
		*/
		Source* GetSource();

		/**
		* \brief Get region image
		*/
		Img* GetImage(RegionImageMode mode);
		/**
		* \brief Draw region
		*/
		void Draw(RegionImageMode mode=eFluxMap);


		void AddSubRegionId(int id){fSubRegionIds.push_back(id);}


	private:

		int ComputeRobustStats(RobustStats& rstats);

		static bool compareNeighborByAscendingDist(const NeighborInfo &a, const NeighborInfo &b) {
   		return ( (a.D) < (b.D) );
		}
		static bool compareNeighborByDescendingDist(const NeighborInfo &a, const NeighborInfo &b) {
   		return ( (a.D) > (b.D) );
		}
		static bool compareNeighborByAscendingTotDist(const NeighborInfo &a, const NeighborInfo &b,double beta) {
   		return ( (a.D+beta*a.Edgeness) < (b.D+beta*b.Edgeness) );
		}
		static bool compareNeighborByDescendingTotDist(const NeighborInfo &a, const NeighborInfo &b,double beta) {
   		return ( (a.D+beta*a.Edgeness) > (b.D+beta*b.Edgeness) );
		}

		static bool compareNeighborDist(const NeighborInfo &a, const NeighborInfo &b,bool isAscending,double beta) {	
			bool res= false;
			if(isAscending) res= ( (a.D+beta*a.Edgeness) < (b.D+beta*b.Edgeness) );
			else res= ( (a.D+beta*a.Edgeness) > (b.D+beta*b.Edgeness) );
   		return res;
		}
		
		void UpdateMoments(Pixel pixel);	
		void ResetMoments();

	public:

		long int fId;//Region id
		int fNPix;//Number of pixels in region
		
		int fImageSizeX;
		int fImageSizeY;
		double fImageMinX;
		double fImageMaxX;
		double fImageMinY;
		double fImageMaxY;
		double fImageMinS;
		double fImageMaxS;
		double fImageMinScurv;
		double fImageMaxScurv;	
		double fImageMinSedge;
		double fImageMaxSedge;
		double fImageRMS;
		
		//Pixel intensity moments
		double fM1;//1st moment
		double fM2;//2nd moment
		double fM3;//3rd moment
		double fM4;//4th moment
		double fMean;//mean = M1/N
		double fRMS;
		double fKurtosis;
		double fSkewness;
		double fMedian;
		double fMedianRMS;
		double fMedian_curv;
		double fMedianRMS_curv;
		
		double fS;//sum of pixel signals
		double fSmax;//max of pixel signals
		double fSmin;//min of pixel signals
		double fSxx;
		double fSyy;
		double fSxy;

		double fS_curv;//sum of pixel curvature
		double fM1_curv;//1st moment
		double fM2_curv;//2nd moment
		double fMean_curv;
		double fRMS_curv;

		double fS_edge;//sum of edge estimator
		double fH;//sample Entropy
		double fParamH;//entropy of moment vector (mean,rms,...)

		int fPixIdmax;//id of pixel with max signal
		int fPixIdmin;//id of pixel with min signal
		double fWx0;//Signal-weighted X position average
		double fWy0;//Signal-weighted Y position average
		double fX0;//X position average
		double fY0;//Y position average

		double fXmin;
		double fXmax;
		double fYmin;
		double fYmax;
		int fIx_min;
		int fIx_max;
		int fIy_min;
		int fIy_max;
	
		TVector3 fColorSum;
		TVector3 fColorM1;
		TVector3 fColorM2;
		TVector3 fColor;
	
		PixelCollection fPixelCollection;

		std::vector<int> fNeighbourRegions;
		std::vector<NeighborInfo> fNeighbourRegionInfo;
		std::vector<int> fSubRegionIds;
		

		bool fIsSignificative;
		double fMahalanobisDistance;
		bool fIsSalient;
		double fSaliency;
		int fTag;

		double fAbsDissMin;
		double fAbsDissMax;
		bool fIsOccluded;
		
		

		Contour* fContour;
		std::vector<Contour*> fContourCollection;
		bool fHasParameters;
		std::vector<int> fContourPointIndex;
		//std::vector<cv::Point2f> fAroundContourPoints;
		std::vector<EdgeInfo> fAroundContourPoints;
		//std::map<int,std::vector<cv::Point2f>> fAroundContourPointsMap;

	private:
	
		friend class sorter;	
		
};

class sorter {
	
public:
	sorter(bool isAscending,double beta) : beta_(beta),isAscending_(isAscending) {}
  bool operator()(const Region::NeighborInfo &a, const Region::NeighborInfo &b) const {
  	return Region::compareNeighborDist(a, b, isAscending_, beta_);
  }
	bool isAscending_;
	double beta_;
};


#ifdef __MAKECINT__
#pragma link C++ class Region+;
#pragma link C++ class vector<Region>+;
#pragma link C++ class vector<Region*>+;
#endif


#endif
