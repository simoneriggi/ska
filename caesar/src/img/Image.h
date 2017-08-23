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

/**
* @file Image.h
* @class Image
* @brief Image
*
* Image class
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _IMAGE_h
#define _IMAGE_h 1

#include <Consts.h>
#include <ImgStats.h>
#include <ImgMetaData.h>

#include <FITSReader.h>
#include <AstroUtils.h>
#include <StatsUtils.h>
#include <BkgFinder.h>
#include <MorphFilter.h>
#include <GraphicsUtils.h>
#include <CodeUtils.h>
#include <Logger.h>

//WCSTOOLS
#include <wcs.h>


//ROOT
#include <TH2F.h>
#include <TMath.h>
#include <TMatrixD.h>
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
#include <iostream>
#include <time.h>
#include <ctime>

using namespace std;

namespace cv {
	class Mat;
}

namespace Caesar{


//class Image : public TObject {
class Image : public TNamed {

  public:
		
    /** 
		\brief Class constructor: initialize structures.
 		*/
		Image();
    Image(const Image& img);		
		Image(long int nbinsx,long int nbinsy,float xlow=-1,float ylow=-1,std::string name="");	
		Image(long int nbinsx,long int nbinsy,std::vector<float>const& pixels,float xlow=-1,float ylow=-1,std::string name="");
		Image(long int nbinsx,long int nbinsy,float w,float xlow=-1,float ylow=-1,std::string name="");
		Image& operator=(const Image &img);
		virtual void Copy(TObject &hnew) const;

		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~Image();

		

	public:
				
		//======================================
		//==       SETTERS/GETTERS method  
		//======================================
		/**
		* \brief Set image name
		*/
		//void SetName(std::string name){m_name= name;}
		void SetName(std::string name){TNamed::SetName(name.c_str());}
		/**
		* \brief Get image name
		*/
		//std::string GetName(){return m_name;}
		std::string GetName(){return std::string(TNamed::GetName());}

		/**
		* \brief Get npixels
		*/
		long int GetNPixels(){return m_Nx*m_Ny;}

		/**
		* \brief Return pixel list
		*/
		const std::vector<float>& GetPixels() const {return m_pixels;}

		/**
		* \brief Get pixel vector size
		*/
		long int GetPixelDataSize(){
			return (long int)(m_pixels.size());
		}
		/**
		* \brief Has pixel data
		*/
		bool HasPixelData() {
			long int n= GetPixelDataSize();
			return (n>0 && n==GetNPixels());
		}

		/**
		* \brief Set image size (NB: this cleares pixel vector and allocated space!!!)
		*/
		int SetSize(long int Nx,long int Ny,float xlow=0,float ylow=0);

		/**
		* \brief Get size
		*/
		void GetSize(long int& Nx,long int& Ny){Nx=m_Nx; Ny=m_Ny;}
		/**
		* \brief Get size x
		*/
		long int GetNx(){return m_Nx;}
		/**
		* \brief Get size y
		*/
		long int GetNy(){return m_Ny;}

		/**
		* \brief Get xmin
		*/
		float GetXmin(){return m_Xmin;}

		/**
		* \brief Get xmax
		*/
		float GetXmax(){return m_Xmin + m_Nx-1;}

		/**
		* \brief Get ymin
		*/
		float GetYmin(){return m_Ymin;}

		/**
		* \brief Get ymax
		*/
		float GetYmax(){return m_Ymin + m_Ny-1;}
		/**
		* \brief Get range
		*/
		void GetRange(float& xmin,float& xmax,float& ymin,float& ymax){
			xmin= GetXmin();
			xmax= GetXmax();
			ymin= GetYmin();
			ymax= GetYmax();
		}

		//================================
		//==       Utils method  
		//================================
		float GetX(long int binx){
			return GetXmin() + binx;
		}
		float GetY(long int biny){
			return GetYmin() + biny;
		}
		long int GetBin(long int binx,long int biny){
			return binx + biny*m_Nx;
		}	
		long int GetBinX(long int gbin){
			return gbin % m_Nx;
		}
		long int GetBinY(long int gbin){
			long binx= GetBinX(gbin);
			return ((gbin-binx)/m_Nx)%m_Ny;
		}	

		/**
		* \brief Get pixel value
		*/
		float GetPixelValue(long int gbin){
			if(gbin<0 || gbin>=(long int)(m_pixels.size())) {
				WARN_LOG("Invalid pixel bin ("<<gbin<<") requested to be accessed (hint: check if image size is not initialized or given bin exceed size), returning nan");
				return std::numeric_limits<float>::quiet_NaN();
			}
			return m_pixels[gbin];
		}

		/**
		* \brief Get pixel value
		*/
		float GetPixelValue(long int ix,long int iy){
			long int gbin= GetBin(ix,iy);	
			return GetPixelValue(gbin);
		}

		/**
		* \brief Set pixel value
		*/
		int SetPixelValue(long int gbin,double w){
			if(gbin<0 || gbin>=(long int)(m_pixels.size())) {
				WARN_LOG("Invalid pixel bin ("<<gbin<<") requested to be accessed (hint: check if image size is not initialized or given bin exceed size), will not set value!");
				return -1;
			}
			m_pixels[gbin]= w;
			return 0;
		}
			
		/**
		* \brief Set pixel value
		*/
		int SetPixelValue(long int ix,long int iy,double w){
			long int gbin= GetBin(ix,iy);	
			return SetPixelValue(gbin,w);
		}

		/**
		* \brief Get bin content (equivalent to GetPixelValue)
		*/
		float GetBinContent(long int ix,long int iy){
			long int gbin= GetBin(ix,iy);	
			return GetPixelValue(gbin);
		}
		/**
		* \brief Get bin content (equivalent to GetPixelValue)
		*/
		float GetBinContent(long int gbin){
			return GetPixelValue(gbin);
		}
	
		/**
		* \brief Set bin content (equivalent to SetPixelValue)
		*/
		int SetBinContent(long int ix,long int iy,double w){
			long int gbin= GetBin(ix,iy);	
			return SetPixelValue(gbin,w);
		}
		/**
		* \brief Set bin content (equivalent to SetPixelValue)
		*/
		int SetBinContent(long int gbin,double w){
			return SetPixelValue(gbin,w);
		}

		/**
		* \brief Clone image (return a new image)
		*/
		Image* GetCloned(std::string name,bool copyMetaData=true,bool resetStats=true){
			Image* clone= new Image;
			*clone= *this;
			clone->SetName(name.c_str());			
			if(!copyMetaData) clone->ClearMetaData();
			if(resetStats) clone->ResetImgStats(true,true);
			return clone;
		}
			
		/**
		* \brief Reset image (zero all pixels and clear stats, size & meta-data preserved)
		*/
		void Reset(){
			//Set all pixels to 0
			std::fill(m_pixels.begin(), m_pixels.end(), 0);

			//Reset and delete stats
			this->ResetImgStats(true,true);
		}
		
		/**
		* \brief Get pixel minimum
		*/
		float GetMinimum(){
			return m_StatMoments.minVal;
		}

		/**
		* \brief Get pixel maximum
		*/
		float GetMaximum(){
			return m_StatMoments.maxVal;
		}

		/**
		* \brief Get pixel world coordinate
		*/
		/*
		int GetWorldCoord(int ix,int iy,double& xpos, double& ypos) {
			return Caesar::AstroUtils::PixelToWCSCoords(this,ix,iy,xpos,ypos);
		}
		*/
		
		/**
		* \brief Check if given bin id is within image range
		*/
		bool HasBin(long int gbin){
			return ( gbin>=0 && gbin<(long int)(m_pixels.size()) );
		}
		bool HasBin(long int binIdX,long int binIdY){
			long int binId= GetBin(binIdX,binIdY);
			return HasBin(binId);
		}
		/*
		bool HasBin(double x,double y){
			int binId= this->FindBin(x,y);
			return HasBin(binId);
		}
		bool HasBinX(int binIdX){
			int Nx= this->GetNbinsX();
			return (binIdX>=1 && binIdX<=Nx );
		}
		bool HasBinX(double x){
			int binIdX= this->GetXaxis()->FindBin(x);
			return HasBinX(binIdX);
		}
		bool HasBinY(int binIdY){
			int Ny= this->GetNbinsY(); 
			return (binIdY>=1 && binIdY<=Ny );
		}
		bool HasBinY(double y){
			int binIdY= this->GetYaxis()->FindBin(y);
			return HasBinY(binIdY);
		}
		*/

		/**
		* \brief Check if bin content is within given value range
		*/
		bool IsBinContentInRange(long int binIdX,long int binIdY,double minThr,double maxThr){	
			if(!HasBin(binIdX,binIdY)) return false;
			double w= this->GetBinContent(binIdX,binIdY);
			return (w>=minThr && w<=maxThr);
		}

		/**
		* \brief Check if bin is at edge
		*/
		bool IsEdgeBin(long int binIdX,long int binIdY){
			return ( binIdX==0 || binIdX==(m_Nx-1) || binIdY==0 || binIdY==(m_Ny-1) );
		}
		/**
		* \brief Check if bin is at edge
		*/
		bool IsEdgeBin(long int gbin){
			long int binIdX= GetBinX(gbin);
			long int binIdY= GetBinY(gbin);
			return IsEdgeBin(binIdX,binIdY);
		}

		//================================
		//==       Fill methods
		//================================
		/**
		* \brief Check filled pixel
		*/
		int CheckFillPixel(long int gbin,double w);

		/**
		* \brief Fill pixels (to be used to compute stats at fill time). NB: use only in single-thread otherwise computed moments are not correct (+ race conditions)!!!
		*/
		int FillPixel(long int ix,long int iy,double w,bool useNegativePixInStats=true);
		int FillPixel(long int gbin,double w,bool useNegativePixInStats=true);
		/**
		* \brief Fill pixel (multithreaded version). This updated the moments per thread. To get the cumulative moments use the StatsUtils
		*/
		#ifdef OPENMP_ENABLED
		int FillPixelMT(Caesar::StatMoments<double>& moments,long int ix,long int iy,double w,bool useNegativePixInStats=true);
		int FillPixelMT(Caesar::StatMoments<double>& moments,long int gbin,double w,bool useNegativePixInStats=true);
		#endif		

		/**
		* \brief Fill pixels (to be used to compute stats at fill time)
		*/
		void FillFromMat(cv::Mat&,bool useNegativePixInStats=true);

		

		//================================
		//==       READ/WRITE methods
		//================================
		/**
		* \brief Read image from FITS file
		*/
		int ReadFITS(std::string filename,int hdu_id=1,int ix_min=-1,int ix_max=-1,int iy_min=-1,int iy_max=-1);
		/**
		* \brief Write image to FITS file
		*/
		int WriteFITS(std::string outfilename);
		/**
		* \brief Get image subregion or tile
		*/
		Image* GetTile(long int ix_min,long int ix_max,long int iy_min,long int iy_max,std::string imgname="");
		/**
		* \brief Read image from an image file
		*/
		//int ReadImageFile(std::string filename);

		//================================
		//==       METADATA methods
		//================================
    /**
		* \brief Set image metadata information
		*/
		void CopyMetaData(ImgMetaData* data) {
			ClearMetaData();
			if(!data) return;
			if(!m_MetaData) m_MetaData= new ImgMetaData;
			*m_MetaData = *data;
			m_HasMetaData= true;
		}
		int SetMetaData(ImgMetaData* data) {
			if(!data){
				ERROR_LOG("Null ptr to given metadata!");
				return -1;
			}
			m_MetaData= data;
			m_HasMetaData= true;
			return 0;
		}
    /**
		* \brief Get image metadata information
		*/
		ImgMetaData* GetMetaData() {return m_MetaData;}
		/**
		* \brief Has metadata information
		*/
		bool HasMetaData(){return (m_HasMetaData && m_MetaData);}
		

		//================================
		//==       STATS methods
		//================================
		/**
		* \brief Update moments (multithreaded version)
		*/
		int ComputeMoments(bool skipNegativePixels=false);
		/**
		* \brief Get stats information
		*/
		ImgStats* GetPixelStats(){return m_Stats;}
		/**
		* \brief Check if stats has been computed 
		*/
		bool HasStats(){return (m_HasStats && m_Stats);}
		/**
		* \brief Compute stats information 
		*/
		int ComputeStats(bool computeRobustStats,bool skipNegativePixels=false,bool forceRecomputing=false);
		/**
		* \brief Print stats information 
		*/
		void PrintStats(){
			if(!HasStats()) return;
			m_Stats->Print();
		}//close PrintStats()
		/**
		* \brief Print stats information 
		*/
		void LogStats(std::string level="INFO"){
			if(!HasStats()) return;
			m_Stats->Log(level);
		}//close LogStats()

		/**
		* \brief Set stat moments
		*/
		void SetMoments(Caesar::StatMoments<double> moments){
			m_StatMoments= moments;
			this->ResetImgStats(false,true);
		}
		
		/**
		* \brief Get stat moments
		*/
		Caesar::StatMoments<double> GetMoments(){return m_StatMoments;}

		/**
		* \brief Get tile mean stats
		*/
		int GetTileMeanStats(float& mean,float& stddev,long int& npix,long int ix_min,long int ix_max,long int iy_min,long int iy_max,bool skipNegativePixels=false); 
		
		/**
		* \brief Get tile median stats
		*/
		int GetTileMedianStats(float& median,float& mad,long int& npix,long int ix_min,long int ix_max,long int iy_min,long int iy_max,bool skipNegativePixels=false); 
		
		/**
		* \brief Get tile median stats
		*/
		int GetTileClippedStats(Caesar::ClippedStats<float>& clipped_stats,long int& npix,double clipsigma,long int ix_min,long int ix_max,long int iy_min,long int iy_max,bool skipNegativePixels=false); 
		
		/**
		* \brief Get tile biweight stats
		*/
		int GetTileBiWeightStats(float& bwLocation,float& bwScale,long int& npix,double C,long int ix_min,long int ix_max,long int iy_min,long int iy_max,bool skipNegativePixels=false); 
		
		//================================
		//==       BKG METHODS
		//================================
		/**
		* \brief Compute local bkg
		*/
		ImgBkgData* ComputeBkg(int estimator=eMedianBkg,bool computeLocalBkg=true,int boxSizeX=100,int boxSizeY=100, double gridStepSizeX=10, double gridStepSizeY=10, bool use2ndPass=true,bool skipOutliers=false,double seedThr=5,double mergeThr=2.6,int minPixels=10);
		/**
		* \brief Compute significance map
		*/
		Image* GetSignificanceMap(ImgBkgData* bkgData,bool useLocalBkg=false);

		//================================
		//==    THRESHOLDING METHODS
		//================================
		//== Thresholding methods ==
		/**
		* \brief Compute pixel histo
		*/
		TH1D* GetPixelHisto(int nbins=100,bool normalize=false);

		/**
		* \brief Find optimal global threshold
		*/
		double FindOptimalGlobalThreshold(int nbins=100);

		/**
		* \brief Find Otsu threshold
		*/
		double FindOtsuThreshold(int nbins=100);
		/**
		* \brief Find valley threshold
		*/
		double FindValleyThreshold(int nbins=100,bool smooth=true);
		/**
		* \brief Find median global threshold
		*/
		//double FindMedianThreshold(double thrFactor);
		/**
		* \brief Find optimal global threshold
		*/
		//double FindOptimalGlobalThreshold(double thrFactor,int nbins=100,bool smooth=true);
		/**
		* \brief Get binarized image
		*/
		Image* GetBinarizedImage(double threshold,double fgValue=1,bool isLowerThreshold=false);


		//=========================================
		//==    SOURCE FINDING METHODS           ==
		//=========================================
		/**
		* \brief Find compact sources
		*/
		int FindCompactSource(std::vector<Source*>&,Image* floodImg=0,ImgBkgData* bkgData=0,double seedThr=5,double mergeThr=2.6,int minPixels=10,bool findNegativeExcess=false,bool mergeBelowSeed=false,bool findNestedSources=false,double nestedBlobThrFactor=1);
		/**
		* \brief Find nested sources
		*/
		int	FindNestedSource(std::vector<Source*>& sources,ImgBkgData* bkgData=0,int minPixels=5,double nestedBlobThreshold=1);

		/**
		* \brief Find extended sources with ChanVese method
		*/
		int FindExtendedSource_CV(std::vector<Source*>&,ImgBkgData* bkgData=0,int minPixels=10,bool findNegativeExcess=false,double dt=0.1,double h=1,double lambda1=1.0,double lambda2=2.0,double mu=0.5,double nu=0,double p=1);

		/**
		* \brief Find extended sources with Hierarchical Clustering method
		*/
		//int FindExtendedSource_HClust(std::vector<Source*>&,Img* saliencyImg,Img* edgeImg);

		
		//=========================================
		//==     FILTER METHODS                  ==
		//=========================================
		/**
		* \brief Get masked image
		*/
		Image* GetMask(Image* mask,bool isBinary=false);
		
		/**
		* \brief Get source masked image
		*/
		Image* GetSourceMask(std::vector<Source*>const& sources,bool isBinary=false,bool invert=false);

		/**
		* \brief Returns a residual image obtained by dilating given sources with a random background
		*/
		Image* GetSourceResidual(std::vector<Source*>const& sources,int KernSize=5,int dilateModel=MorphFilter::eDilateWithBkg,int dilateSourceType=-1,bool skipToNested=false,ImgBkgData* bkgData=0,bool useLocalBkg=false,bool randomize=false,double zThr=20);


		/**
		* \brief Scale image by a factor 'c'
		*/
		int Scale(double c);

		/**
		* \brief Add an image to this (this= this + img*c)
		*/
		int Add(Image* img,double c=1,bool computeStats=false);

		/**
		* \brief Get normalized image
		*/
		Image* GetNormalizedImage(std::string normScale="LINEAR",int normmin=1,int normmax=256, bool skipEmptyBins=true);

		/**
		* \brief Get image laplacian
		*/
		Image* GetLaplacianImage(bool invert=false);

		/**
		* \brief Get guided filter image
		*/
		Image* GetGuidedFilterImage(int radius=12,double eps=0.04);
		/**
		* \brief Smooth image
		*/
		Image* GetSmoothedImage(int size_x=3,int size_y=3,double sigma_x=1,double sigma_y=1);

		/**
		* \brief Get single-reso saliency map
		*/
		Image* GetSaliencyMap(int reso=20,double regFactor=1,int minRegionSize=10,double knnFactor=1,bool useRobust=false,double expFalloffPar=100,double distanceRegPar=1);

		/**
		* \brief Get multi-reso saliency map
		*/
		Image* GetMultiResoSaliencyMap(int resoMin=20,int resoMax=60,int resoStep=10,double beta=1,int minRegionSize=10,double knnFactor=0.2,bool useRobustPars=false,double expFalloffPar=100,double distanceRegPar=1,double salientMultiplicityThrFactor=0.7,bool addBkgMap=true,bool addNoiseMap=true,ImgBkgData* bkgData=0,double saliencyThrFactor=2,double imgThrFactor=1);

		/**
		* \brief Get image wavelength decomposition
		*/
		std::vector<Image*> GetWaveletDecomposition(int nScales);
		
		/**
		* \brief Get image gradient
		*/
		Image* GetGradientImage(bool invert=false);

		/**
		* \brief Get kirsch image
		*/
		Image* GetKirschImage();
				
		/**
		* \brief Get laplacian of gaussian image
		*/
		Image* GetLoGImage(bool invert=false);	
		
		/**
		* \brief Get scale-normalized laplacian of gaussian image
		*/
		Image* GetNormLoGImage(int size=3,double scale=1,bool invert=false);
		
		

		//=========================================
		//==   CONVERT METHODS
		//=========================================
		/**
		* \brief Export to ROOT TH2
		*/
		TH2D* GetHisto2D(std::string histoname="");

		/**
		* \brief Convert image to OpenCV mat float image
		*/
		cv::Mat GetOpenCVMat(std::string encoding="64");	
		
		/**
		* \brief Get ROOT matrix from image
		*/
		TMatrixD* GetMatrix();
		

		
		//==================================
		//==    DRAW METHODS             ===
		//==================================
		/**
		* \brief Draw image (plain image, no WCS, no sources)
		*/
		int Draw(int palette=Caesar::eRAINBOW,bool drawFull=false,bool useCurrentCanvas=false,std::string units="");
		/**
		* \brief Draw image and sources (no WCS axis)
		*/
		int Draw(std::vector<Source*>const& sources,int palette=Caesar::eRAINBOW,bool drawFull=false,bool useCurrentCanvas=false,std::string units="");

		/**
		* \brief Draw image
		*/
		//int Plot(std::vector<Source*>const&,bool useCurrentCanvas=true,bool drawFull=false,int paletteStyle=Caesar::eRAINBOW,bool drawColorPalette=true,bool putWCAxis=false,int coordSystem=-1,std::string units="Jy/beam");
		/**
		* \brief Set draw range
		*/
		//void SetDrawRange(double zmin,double zmax);


		

	protected:
		/**
		* \brief Init data
		*/
		void Init();
		/**
		* \brief Reset stats and/or moments
		*/
		void ResetImgStats(bool resetMoments=false,bool clearStats=false);
		/**
		* \brief Clear stats
		*/
		void ClearImgStats(){
			if(!m_Stats) return;
			delete m_Stats;
			m_Stats= 0;
			m_HasStats= false;
		}
		
		/**
		* \brief Update moments
		*/
		//void UpdateMoments(double w);	
		
		/**
		* \brief Compute image stats parameters from moments 
		*/
		void ComputeStatsParams(bool computeRobustStats=true,bool skipNegativePixels=false);
		/**
		* \brief Clear stats
		*/
		void ClearMetaData(){
			if(!m_MetaData) return;
			delete m_MetaData;
			m_MetaData= 0;
			m_HasMetaData= false;
		}

		/**
		* \brief Get tile pixels
		*/
		int GetTilePixels(std::vector<float>& pixels,long int ix_min,long int ix_max,long int iy_min,long int iy_max,bool skipNegativePixels=false);

		
	protected:
		//Image name
		std::string m_name;

		//Pixel list
		std::vector<float> m_pixels;

		//Image size
		long int m_Nx;
		long int m_Ny;

		//X Y min coordinates (normally (0,0))
		float m_Xmin;
		float m_Ymin;
	
		//Meta-data
		bool m_HasMetaData;
		ImgMetaData* m_MetaData;
		
		//Stats data
		bool m_HasStats;
		ImgStats* m_Stats;
		StatMoments<double> m_StatMoments;//stat moments & min/max

		/*
		long long int m_Npix;//npixels	
		double m_M1;//1st moments
		double m_M2;//2nd moment
		double m_M3;//3rd moment
		double m_M4;//4th moment
		double m_PixelMin;
		double m_PixelMax;
		*/

		ClassDef(Image,1)

};

#ifdef __MAKECINT__
#pragma link C++ enum Image::ImageType+;
//#pragma link C++ class StatsData+;
//#pragma link C++ class MetaData+;
#pragma link C++ class Image+;
#pragma link C++ class vector<Image>+;
#pragma link C++ class vector<Image*>+;
#endif

}//close namespace

#endif
