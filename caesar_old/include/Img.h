/**
* @file Img.h
* @class Img
* @brief Img
*
* Image class
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef Img_h
#define Img_h 1

#include "Utils.h"
#include "Region.h"
#include "Source.h"
#include "GradientFilter.h"
#include "Contour.h"

//#include <wcs.h>


#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
//#include <opencv2/ximgproc.hpp>
//#include <opencv2/ximgproc/edge_filter.hpp>
#include <opencv2/saliency/saliencyBaseClasses.hpp>
#include <opencv2/saliency.hpp>

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

class Utils;
class TH2D;
class TH2F;
class Img;
typedef std::vector<Img> ImgCollection;
typedef std::vector<Img*> ImgPtrCollection;		


//class Img : public TH2D {
class Img : public TH2F {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    Img();
		Img(const Img& img);
		Img(const char* name, const char* title, int nbinsx, const float* xbins, int nbinsy, const float* ybins);
		Img(const char* name, const char* title, int nbinsx, float xlow, float xup, int nbinsy, float ylow, float yup);

		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~Img();

		
		//## Img metadata
		struct MetaData {
	    int Nx;
			int Ny;
			int Cx;
			int Cy;
			double Xc;
			double Yc;
			double dX;
			double dY;
			double RotX;
			double RotY;
			std::string CoordTypeX;
			std::string CoordTypeY;
			std::string BUnit;
			double Bmaj;
			double Bmin;
			double Bpa;
			double Epoch;
		};

    
		struct StatsData {
			int n;
			double min;
			double max;
			double mean;
			double meanErr;
			double rms;
			double rmsErr;
      double skewness;	
			double skewnessErr;			
			double kurtosis;
      double kurtosisErr;
			double median;
			double medianRMS;
			double bwLocation;
			double bwScale;
			double bwLocationIter;
			double bwScaleIter;
			double clippedMedian;
			double clippedRMS;
			double fitmean;
			double fitmeanErr;
			double fitrms;
			double fitrmsErr;
		};

		struct BkgData {
			int ix_min;
			int iy_min;
			int ix_max;
			int iy_max;
			int npix;
			bool isReliableBkg;
			double bkgLevel;
			double bkgRMS;
		};

		//## Bkg finder options
		enum BkgMethod {eMeanBkg=1,eMedianBkg=2,eRobustBkg=3,eSimpleRobustBkg=4,eDummyRobustBkg=5,eBiWeightBkg=6,eMedianClippedBkg=7};
		enum LocalBkgMethod {eGridBkg=1,eSuperpixelBkg=2};

		struct BkgOptions {			
			BkgMethod method;
			bool skipNegativePixels;
			double boxSizeX;
			double boxSizeY;	
		};

		enum KernelType {eSOBEL_HORIZ=1,eSOBEL_VERT=2,eSCHARR_HORIZ=3,eSCHARR_VERT=4,eLAPLACE=5,eLoG=6,eNormLoG=7,eKIRSCH=8};
		enum DilationModel {eBkg=1,eSourceMean=2,eSourceMedian=3};
		
	public:
	
		/**
		* \brief Draw plain image
		*/
		TCanvas* DrawPlain(bool useLogScale=false);

		//Fill image
		/**
		* \brief FillPixel
		*/
		int FillPixel(double x,double y,double w, bool includeNegativePixelsInStats=true);
		void SetPixel(int ix,int iy,double w, bool includeNegativePixelsInStats=true);
		
		int ReadFile(std::string filename);
		int ReadFITSFile(std::string filename,int ix_min=-1,int ix_max=-1,int iy_min=-1,int iy_max=-1);
		
		/**
		* \brief Get size of synthetic beam
		*/
		double GetBeamSize(){
			if(!fHasMetaData) return 0;
			double beamSize= sqrt(fabs(metadata.Bmaj*metadata.Bmin));
			return beamSize;	
		}
		/**
		* \brief Get pixel scale
		*/
		int GetBeamSizeInPixel(){
			if(!fHasMetaData) return 0;
			double beamSize= GetBeamSize();
			double pixScale= sqrt(fabs(metadata.dX*metadata.dY));
			int npixInBeam= int(ceil(beamSize/pixScale));
			return npixInBeam;
		}

		//Filter functions
    /**
		* \brief Get image wavelength decomposition
		*/
		ImgPtrCollection GetWaveletDecomposition(int nScales);
		/**
		* \brief Smooth image
		*/
		Img* Smooth(int size_x=3,int size_y=3,double sigma_x=1,double sigma_y=1);

		/**
		* \brief Compute saliency map (multi parameter version)
		*/	
		Img* GetSaliencyMap_MultiParVersion(int reso=20,double regFactor=0.01,int minRegionSize=5,double knnFactor=0.1,bool useRobust=false,bool addCurvDist=true,bool interpolate=false,double expFalloffPar=100,double distanceRegPar=1);
		Img* GetSaliencyMap_MultiParVersion(std::vector<Region*>& regions,double knnFactor=0.1,bool useRobust=false,bool addCurvDist=true,bool interpolate=false,double expFalloffPar=100,double distanceRegPar=1);
	
		/**
		* \brief Compute saliency map (single color parameter version)
		*/
		Img* GetSaliencyMap(int reso=20,double regFactor=1,int minRegionSize=10,double knnFactor=0.2,bool useRobustPars=false,double expFalloffPar=100,double distanceRegPar=1,bool interpolate=false);
		Img* GetSaliencyMap(std::vector<Region*>& regions,double knnFactor=0.2,bool useRobustPars=false,double expFalloffPar=100,double distanceRegPar=1,bool interpolate=false);

		/**
		* \brief Compute smoothed saliency map (Parazzi et al method)
		*/
		Img* GetSmoothedSaliencyMap(std::vector<Region*> regions,double sigmaS=0.5,double sigmaC=0.5,double saliencyKFactor=6,bool useRobust=false,bool addCurvDist=true);
		Img* GetSmoothedSaliencyMap(int reso=20,double regFactor=0.01,int minRegionSize=5,double sigmaS=0.5,double sigmaC=0.5,double saliencyKFactor=6,bool useRobust=false,bool addCurvDist=true);
		
		/**
		* \brief Compute saliency map according to Luo et al method (ref. Luo et al, EUVIP, Jun 2013, Paris, France. 2013)
		*/	
		Img* GetSaliencyMap_LuoMethod(int reso=20,double regFactor=1,int minRegionSize=10,double expFalloffPar=1,double distanceRegPar=1);
		Img* GetSaliencyMap_LuoMethod(std::vector<Region*>& regions,double expFalloffPar=1,double distanceRegPar=1);
		
		/**
		* \brief Compute saliency map according to Zhang et al method (ref. )
		*/
		Img* GetSaliencyMap_GofermanMethod(int reso=20,double regFactor=1,int minRegionSize=10,double knnFactor=0.1,double distanceRegPar=1,double expFalloffPar=1);
		Img* GetSaliencyMap_GofermanMethod(std::vector<Region*>& regions,double knnFactor=0.1,double distanceRegPar=1,double expFalloffPar=1);
	
		/**
		* \brief Compute saliency map according to spectral method (ref. )
		*/
		Img* GetSaliencyMap_SpectralRes();
		
		/**
		* \brief Compute multi-resolution saliency map (ref. )
		*/
		//Img* GetMultiResoSaliencyMap(int resoMin=20,int resoMax=60,int resoStep=10,double beta=0.01,int minRegionSize=1,double knnFactor=0.1,bool useRobust=false,bool addCurvDist=true,double thr=0.3,bool addCurvMap=false,bool addBkgMap=true,bool addNoiseMap=true,int normalizationAcrossResoMode=2,double medianThrFactor=2,double medianImgThrFactor=1);

		Img* GetMultiResoSaliencyMap(int resoMin=20,int resoMax=60,int resoStep=10,double beta=1,int minRegionSize=10,double knnFactor=0.2,bool useRobustPars=false,double expFalloffPar=100, double distanceRegPar=1, double thr=0.3,bool addCurvMap=false,bool addBkgMap=true,bool addNoiseMap=true,int normalizationAcrossResoMode=2,double medianThrFactor=2,double medianImgThrFactor=1);
		
		/**
		* \brief Compute multi-resolution smoothed saliency map (ref. )
		*/
		Img* GetMultiResoSmoothedSaliencyMap(int resoMin=20,int resoMax=60,int resoStep=10,double beta=0.01,int minRegionSize=1,double sigmaS=0.5,double sigmaC=0.5,double saliencyKFactor=6,bool useRobust=false,bool addCurvDist=true,double thr=0.3,bool addCurvMap=false,bool addBkgMap=true,bool addNoiseMap=true);
	
		/**
		* \brief Get multi-scale LoG convolved image
		*/
		Img* GetMultiScaleBlobFilterMap(int kernelSize=7,double sigmaMin=1,double sigmaMax=5,double sigmaStep=1,double significantScaleFractionThr=0.5,int normalizationAcrossScaleMode=3,double thrFactor=2,double imgThrFactor=1);
		Img* GetMultiScaleBlobMask(int kernelSize=7,double sigmaMin=1,double sigmaMax=5,double sigmaStep=1,int thrModel=1,double thrFactor=2);

    /**
		* \brief Get matrix from image
		*/
		TMatrixD* GetMatrix();
		/**
		* \brief Get normalized image
		*/
		Img* GetNormalizedImage(int normmin=1,int normmax=256, bool skipEmptyBins=true);
		Img* GetLogNormalizedImage(int normmin=1,int normmax=256,bool skipEmptyBins=true);
		Img* GetSigmoidNormalizedImage(int normmin=1,int normmax=256,double x0=130, double sigma=1, bool skipEmptyBins=true);
		/**
		* \brief Get binarized image
		*/
		Img* GetBinarized(double threshold,double bkgValue=0,double fgValue=1,bool isLowerThreshold=false);
		/**
		* \brief Get image gradients along x and y
		*/
		std::vector<Img*> GetGradientImages(std::string KernelTyle="SCHARR");
		/**
		* \brief Get gradient image
		*/
		Img* GetGradientImage(std::string KernelTyle="SCHARR");
		/**
		* \brief Get laplacian image
		*/
		Img* GetLaplacianImage(bool invert=false);
		/**
		* \brief Get laplacian weighted image
		*/
		Img* GetLaplacianWeightedImage(double f=0.7,std::string filtType="LAPL",int size=3,double scale=1);
		/**
		* \brief Get laplacian of gaussian image
		*/
		Img* GetLoGImage(bool invert=false);	
		/**
		* \brief Get scale-normalized laplacian of gaussian image
		*/
		Img* GetNormLoGImage(int size=3,double scale=1,bool invert=false);
		/**
		* \brief Get Kirsch convoluted image
		*/
		Img* GetKirschImage();
		/**
		* \brief Get bilateral filter image
		*/
		Img* GetBilateralFilterImage(int KernelSize=5, double sigmaColor=10, double sigmaSpace=10);
		/**
		* \brief Get guided filter image
		*/
		Img* GetGuidedFilterImage(int radius=12,double eps=0.04);
		/**
		* \brief Find edges
		*/
		Img* FindEdges(std::vector<int>& edgePixelIds,std::string filter="LoG",int size=3,double scale=1);
		/**
		* \brief Find edges with Canny method
		*/
		Img* FindCannyEdges(std::vector<int>& edgePixelIds,int KernelSize=5,double lowThreshold=150,double thresholdRatio=1.5);
		
		/**
		* \brief Returns gabor filtered image
		*/
		Img* GetGaborFilteredImage(int KernelSize=31,int nThetaSteps=18,double sigma=10,double lambda=20,double gamma=1,double psi=0);

		/**
		* \brief Get linear Hough Transform
		*/
		std::vector<TLine*> LinearHoughTransform(double edgeThreshold=100,double accumulatorReso=1,double thetaReso=1,int threshold=30, double minLineLength=0, double maxLineGap=0);

		/**
		* \brief Get circle Hough Transform
		*/
		std::vector<TEllipse*> HoughTransform(int accumulatorReso=1,int minCircleDistance=10,int minCircleRadius=5,int maxCircleRadius=1000,double edgeThreshold=100,double accumulatorThreshold=30);
		/**
		* \brief Get skeleton image
		*/
		Img* GetSkeleton(int KernelSize=3,double threshold=127,int maxIter=100);

		/**
		* \brief Returns image dilated with given kernel
		*/
		Img* Dilate(TGraph& peaks,int KernSize=5);
		/**
		* \brief Find peaks in image with multiple dilation operation
		*/
		TGraph* FindPeaks(int peakShiftTolerance=1);
		std::vector<int> FindPeakIds(int peakShiftTolerance=1);
		
		TGraph* FindMorphPeakIds(int KernMin=3,int KernMax=7,int KernStep=2,int peakShiftTolerance=1,bool isValley=false);

		/**
		* \brief Apply morphology operation on input image 
		*/
		Img* MorphFilter(std::string operation,int KernSize=5,int structElementType=1);

		

		/**
		* \brief Convert image to OpenCV mat float image
		*/
		cv::Mat ImgToMat(std::string encoding="64");

		/**
		* \brief Find image contour
		*/
		std::vector<Contour*> FindContour();
		
		//Metadata functions
    /**
		* \brief Set image metadata information
		*/
		void SetMetaData(MetaData data) {metadata= data;fHasMetaData=true;}
    /**
		* \brief Get image metadata information
		*/
		MetaData GetMetaData() {return metadata;}
		 /**
		* \brief Has metadata information
		*/
		bool HasMetaData(){return fHasMetaData;}

   	/**
		* \brief Get world coordinates
		*/
		bool GetWorldCoord(int ix,int iy,double& xpos, double& ypos);
		

		// Stats functions
		/**
		* \brief Compute stats information 
		*/
		int ComputeStats(bool computeRobustStats,bool skipNegativePixels=false,bool forceRecomputing=false);
		

    /**
		* \brief Dump stats information 
		*/
		void DumpStats(){
			if(!this->HasStats()) return;
			cout<<"*** IMG STATS ***"<<endl;
			cout<<"N="<<fPixelStats->n<<" min/max="<<fPixelStats->min<<"/"<<fPixelStats->max<<endl;
			cout<<"M1="<<fM1<<" M2="<<fM2<<" M3="<<fM3<<" M4="<<fM4<<endl;
			cout<<"Mean: "<<fPixelStats->mean<<" +- "<<fPixelStats->meanErr<<endl;
			cout<<"RMS: "<<fPixelStats->rms<<" +- "<<fPixelStats->rmsErr<<endl;
			cout<<"Skewness: "<<fPixelStats->skewness<<" +- "<<fPixelStats->skewnessErr<<endl;
      cout<<"Kurtosis: "<<fPixelStats->kurtosis<<" +- "<<fPixelStats->kurtosisErr<<endl;
			cout<<"Median: "<<fPixelStats->median<<" +- "<<fPixelStats->medianRMS<<endl;
			cout<<"BiweightLocation: "<<fPixelStats->bwLocation<<" +- "<<fPixelStats->bwScale<<endl;
			cout<<"*****************"<<endl;
		}//close Img::DumpStats()
		/**
		* \brief Get stats information
		*/
		StatsData* GetPixelStats(){return fPixelStats;}
		/**
		* \brief Check if stats has been computed 
		*/
		bool HasStats(){return (fHasStats && fPixelStats);}
		TH1D* GetPixelHisto(){return fPixelHisto;}

		/**
		* \brief Compute pixel histo
		*/
		TH1D* GetPixelHisto(int nbins,bool normalize=false);
		/**
		* \brief Find Otsu threshold
		*/
		double FindOtsuThreshold(int nbins=100);

		/**
		* \brief Find valley threshold
		*/
		double FindValleyThreshold(int nbins=100,bool smooth=true);

		/**
		* \brief Init tile size (sizing vectors, ...)
		*/
		void InitTiles(int TileSizeX,int TileSizeY);
		/**
		* \brief Set tile size (used for background calculation)
		*/
		void SetTileSize(int TileSizeX,int TileSizeY){
			//fTileSizeX= TileSizeX;
			//fTileSizeY= TileSizeY;
			InitTiles(TileSizeX,TileSizeY);
		}		

		/**
		* \brief Get the absolute tile index given pixel coordinate
		*/
		int GetTileIndex(int ix,int iy){
			if (ix<0 || iy<0 || ix>=this->GetNbinsX() || iy>=this->GetNbinsY()) return -1;
			int tileIdX= ix/fTileSizeX;
			int tileIdY= iy/fTileSizeY;
			int index= tileIdX + fNTilesX*tileIdY;
      return index;
		}	

		// Tile manipulation functions
    /**
		* \brief Extract a sub-image/tile from the image according to given range
		*/
		Img* GetTile(int tileId);
    Img* GetTile(int ix_min,int ix_max,int iy_min,int iy_max);
	
		/**
		* \brief Check if bkg data has been computed 
		*/
		bool HasBkgData(){return (fHasBkgData && (fBkgData.size()>0));}
		/**
		* \brief Get computed bkg data
		*/
		std::vector<BkgData*> GetBkgData(){return fBkgData;}
		/**
		* \brief Set bkg data
		*/
		void SetBkgData(std::vector<BkgData*> vect){
			//fBkgData=vect;
			fBkgData.clear();
			for(unsigned int i=0;i<vect.size();i++){
				BkgData* thisBkgData= new BkgData;
				*thisBkgData= *vect[i];
				fBkgData.push_back(thisBkgData);
			}
			fHasBkgData=true;
		}

		
		/**
		* \brief Get computed bkg map
		*/
		TH2D* GetBkgLevelMap(){return fBackgroundLevelMap;}
		/**
		* \brief Set bkg map
		*/
		void SetBkgLevelMap(TH2D* map){fBackgroundLevelMap= map;}
		/**
		* \brief Get computed bkg rms
		*/
		TH2D* GetBkgRMSMap(){return fBackgroundRMSMap;}
		/**
		* \brief Set bkg rms
		*/
		void SetBkgRMSMap(TH2D* map){fBackgroundRMSMap= map;}
		
	
		/**
		* \brief Get computed bkg map
		*/
		Img* GetInterpolatedBkgLevelMap(){return fInterpolatedBackgroundLevelMap;}
		/**
		* \brief Set computed bkg map
		*/
		void SetInterpolatedBkgLevelMap(Img* map){
			if(!map) return;
			//fInterpolatedBackgroundLevelMap=map;
			if(fInterpolatedBackgroundLevelMap) {
				delete fInterpolatedBackgroundLevelMap;
				fInterpolatedBackgroundLevelMap= 0;
			}
			fInterpolatedBackgroundLevelMap= new Img;
			*fInterpolatedBackgroundLevelMap= *map;
		}
		/**
		* \brief Get computed bkg rms
		*/
		Img* GetInterpolatedBkgRMSMap(){return fInterpolatedBackgroundRMSMap;}
		/**
		* \brief Set computed bkg rms
		*/
		void SetInterpolatedBkgRMSMap(Img* map){
			if(!map) return;
			if(fInterpolatedBackgroundRMSMap) {
				delete fInterpolatedBackgroundRMSMap;
				fInterpolatedBackgroundRMSMap= 0;
			}
			fInterpolatedBackgroundRMSMap= new Img;
			*fInterpolatedBackgroundRMSMap= *map;
		}
		
		
		/**
		* \brief Reset background
		*/
		void ResetBkg(){
			fHasBkgData= false;
			fBkgData.clear();
			fBkgData.resize(0);
			if(fInterpolatedBackgroundLevelMap) {
				delete fInterpolatedBackgroundLevelMap;
				fInterpolatedBackgroundLevelMap= 0;
			}
			if(fInterpolatedBackgroundRMSMap) {
				delete fInterpolatedBackgroundRMSMap;
				fInterpolatedBackgroundRMSMap= 0;
			}
		}
	

		/**
		* \brief Compute global background info
		*/
		int ComputeBkg(Img::BkgMethod method=eRobustBkg);
		/**
		* \brief Compute local background info
		*/
		int ComputeLocalBkg(Img::LocalBkgMethod method=eGridBkg,Img::BkgMethod estimator=eMedianBkg, int boxSizeX=100, int boxSizeY=100, double boxSlideOffsetX=1, double boxSlideOffsetY=1,int SPSize=10,double SPRegularization=100,int SPMinArea=5,bool useTwoPass=true);
		/**
		* \brief Dump bkg information 
		*/
		void DumpBkg(){
			if(!this->HasBkgData()) return;
			cout<<"*** IMG BKG ***"<<endl;
			for(unsigned int j=0;j<fBkgData.size();j++){
				cout<<"Tile "<<j+1<<": bkg="<<fBkgData[j]->bkgLevel<<", rms="<<fBkgData[j]->bkgRMS<<endl;				
			}//end loop tiles
			cout<<"*****************"<<endl;
		}//close Img::DumpBkg()
		/**
		* \brief Get bkg in pixel
		*/
		int GetPixelBkg(int ix,int iy,BkgData& bkgInfo,bool useLocalBackground=false);
		
		//######################
		//## Segmentation
		//######################
		std::vector<Region*> FindSegmentation(int RegionSize,double regularization, int minRegionSize,double eps, TGraph& contours,Img& coloredImg);

		/**
		* \brief Perform Superpixel segmentation
		*/
		Img* FindVLSLICSegmentation(std::vector<Region*>& regions,int RegionSize=10,double regularization=100,int minRegionArea=5,bool useLogContrast=false,bool mergeRegions=true,int algoType=2);
		
		/**
		* \brief Perform ChanVese segmentation
		*/
		Img* FindCVSegmentation(double dt=0.1,double h=1,double lambda1=1.0,double lambda2=2.0,double mu=0.5,double nu=0,double p=1);

		/**
		* \brief Find flood fill pixel list
		*/
		std::vector<int> FloodFill(int seedPixelId,double mergeThr,bool mergeBelowSeed=false);

		/**
		* \brief Find blobs of connected pixels around a seed and above specified thresholds
		*/
		std::vector<Source*> FindBlobs(Img* significanceMap,double seedThr=5,double mergeThr=2.5,int minPixels=10,bool findNegativeExcess=false,bool useLocalBackground=false,bool mergeBelowSeed=false);

		/**
		* \brief Find blobs of connected pixels around a seed and above specified thresholds
		*/
		Source* FindSeededBlob(int seedPixelId,double mergeThr,bool mergeBelowSeed=false);
	
		/**
		* \brief Find multi-scale blobs
		*/
		std::vector<Source*> FindMultiScaleBlobs(int kernelFactor=6,double sigmaMin=1,double sigmaMax=1,double sigmaStep=1,int thrModel=2,double thrFactor=1.5,double seedThr=1,double mergeThr=1,int minNPix=5);
		std::vector<Source*> FindMultiScaleBlobs(Img* blobMask,double seedThr,double mergeThr,int minNPix);

		/**
		* \brief Find compact source using a Flood-Fill method
		*/ 
		int FindCompactSource(Img* significanceMap,double seedThreshold=5, double mergeThreshold=2.5,int minPixels=10,bool findNegativeExcess=false,bool findNestedSources=true,bool useLocalBackground=false,bool mergeBelowSeed=false,double peakThreshold=3);
		/**
		* \brief Find extended sources using a segmentation method
		*/
		int FindExtendedSource(int RegionSize=10,double regularization=100,int minRegionSize=10,double mergeEps=5,double threshold=0);
		
		/**
		* \brief Get source map
		*/
		Img* GetSourceMap(bool isBinary=false);
		/**
		* \brief Get source mask
		*/
		Img* GetSourceMask(std::vector<Source*> sources,bool isBinary);
		/**
		* \brief Get source mask
		*/
		Img* GetSourceMask(std::vector<Source*> sources,int mode);
		/**
		* \brief Get masked image
		*/
		Img* GetMask(Img* mask,bool isBinary=false);

		/**
		* \brief Replace source pixels in image with random background, dilating region around the source with the specified kernel size
		*/
		int DilateSource(Source* source,bool useLocalBackground=false,int KernSize=3,bool skipToNested=false,int sourceType=-1,int dilateModel=eSourceMedian,bool randomize=false,double randSigma=1);

		/**
		* \brief Returns a residual image obtained by dilating given sources with a random background
		*/
		Img* GetSourceResidual(bool useLocalBackground=false,int KernSize=3,bool skipToNested=false,int sourceType=-1,int dilateModel=eSourceMedian,bool randomize=false,double randSigma=1);
		/**
		* \brief Get image sources
		*/
		std::vector<Source*> GetSources(){return fSourceCollection;}
		/**
		* \brief Add source to list
		*/	
		void AddSource(Source* s){
			if(!s) return;
			Source* aSource= new Source;
			*aSource= *s;
			fSourceCollection.push_back(aSource);
			fHasSources= true;
		}
		/**
		* \brief Set sources
		*/
		void SetSources(std::vector<Source*> list){
			ResetSources();
			for(unsigned int i=0;i<list.size();i++){
				Source* aSource= new Source;
				*aSource= *list[i];
				fSourceCollection.push_back(aSource);
			}
			fHasSources= true;
		}
		/**
		* \brief Check if image has sources
		*/
		bool HasSources(){return (fHasSources && (fSourceCollection.size()>0));}
		/**
		* \brief Reset sources
		*/
		void ResetSources(){			
			for(unsigned int i=0;i<fSourceCollection.size();i++){
				if(fSourceCollection[i]) {
					delete fSourceCollection[i];
					fSourceCollection[i]= 0;
				}
			}
			fSourceCollection.clear();
			fHasSources= false;
		}

		/**
		* \brief Compute significance map
		*/
		Img* GetSignificanceMap(bool useLocalBackground=false);
		/**
		* \brief Get significance map from bkg & noise maps
		*/
		Img* GetSignificanceMap(Img* BkgMap, Img* NoiseMap);

		/**
		* \brief Get number of pixels in image
		*/	
		int GetNPixels(){return fNpix;}
		double GetPixelMin(){return fPixelMin;}
		double GetPixelMax(){return fPixelMax;}


		/**
		* \brief Compute Zernike moments
		*/
		std::vector<double> GetZernikeMoments(int order=6,double radius=-1,int method=1);
		/**
		* \brief Compute Hu moments
		*/
		std::vector<double> GetHuMoments();


	private:
		/**
		* \brief Init class attributes
		*/
		void Init();
	
		/**
		* \brief Compute image stats parameters from moments 
		*/
		void ComputeStatsParams(bool computeRobustStats=true,bool skipNegativePixels=false);
		/**
		* \brief Update moments
		*/
		void UpdateMoments(int ix,int iy,double w);
		/**
		* \brief Reset stats and/or moments
		*/
		void ResetStats(bool resetMoments=false);

		
  public:
    std::string errflag;

	private:

		//Metadata
		MetaData metadata;		
		bool fHasMetaData;
	
    //Tile region data
    int fTileSizeX;
		int fTileSizeY;
    int fNTilesX;
		int fNTilesY;	
		
		//Stats data
		int fNpix;//npixels	
		double fM1;//1st moments
		double fM2;//2nd moment
		double fM3;//3rd moment
		double fM4;//4th moment
		double fPixelMin;
		double fPixelMax;
		bool fHasStats;	
		StatsData* fPixelStats;
		TH1D* fPixelHisto;

		//Bkg data
		bool fHasBkgData;
		std::vector<BkgData*> fBkgData;
		TH2D* fBackgroundLevelMap;
		TH2D* fBackgroundRMSMap;
		Img* fInterpolatedBackgroundLevelMap;
		Img* fInterpolatedBackgroundRMSMap;
		
		//Sources
		bool fHasSources;
		std::vector<Source*> fSourceCollection;

		ClassDef(Img,1)

};


#ifdef __MAKECINT__
#pragma link C++ class Img+; 
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<Img*>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<Img>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ struct Img::MetaData+;
//#pragma link C++ class Img::MetaData+;
#endif

#ifdef __MAKECINT__
#pragma link C++ struct Img::StatsData+;
//#pragma link C++ class Img::StatsData+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<Img::StatsData>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ struct Img::BkgData+;
#pragma link C++ struct Img::BkgData*+;
#pragma link C++ class std::vector<Img::BkgData>+;
#pragma link C++ class std::vector<Img::BkgData*>+;
#endif

/*
#ifdef __MAKECINT__
#pragma link C++ class WorldCoor+;
#pragma link C++ class poly+;
#pragma link C++ class wcsprm+;
#pragma link C++ class linprm+;
#pragma link C++ class celprm+;
#pragma link C++ class prjprm+;
#pragma link C++ class IRAFsurface+;
#pragma link C++ class Distort+;
#endif
*/
//#ifdef __MAKECINT__
//#pragma link C class WorldCoor*;
//#endif


#endif


