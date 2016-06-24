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

#include <FITSReader.h>
#include <AstroUtils.h>
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

class Source;


enum WCSType {
	eJ2000= 0,
	eB1950= 1,
	eGALACTIC= 2,
	eECLIPTIC= 3,
	eALTAZ= 4,
	eLINEAR= 5
};

enum ImgFilters {
	eGausFilter= 1,
	eGuidedFilter= 2,
	eWaveletFilter= 3,
	eLoGFilter= 4
};

class ImgMetaData : public TObject {

	public:
		ImgMetaData(){
			Nx= 0; Ny= 0;
			Cx= 0; Cy= 0;
			dX= 0; dY= 0;
			RotX= 0; RotY= 0;
			CoordTypeX= ""; CoordTypeY= "";
			BUnit= "";
			Bmaj= 0; Bmin= 0; Bpa= 0;
			Epoch= 0;
			m_wcsType= "";
		}
		virtual ~ImgMetaData(){};

	public: 
		/**
		* \brief Set cards from FITS file
		*/
		void SetFITSCards(Caesar::FITSFileInfo fits_info){
			Nx= (fits_info.header).Nx;
			Ny= (fits_info.header).Ny;
			Cx= (fits_info.header).Cx;
			Cy= (fits_info.header).Cy;
			Xc= (fits_info.header).Xc;
			Yc= (fits_info.header).Yc;
			CoordTypeX= (fits_info.header).CoordTypeX;
			CoordTypeY= (fits_info.header).CoordTypeY;
			BUnit= (fits_info.header).BUnit;
			Bmaj= (fits_info.header).Bmaj;
			Bmin= (fits_info.header).Bmin;
			Bpa= (fits_info.header).Bpa;
			dX= (fits_info.header).dX;
			dY= (fits_info.header).dY;
			RotX= (fits_info.header).RotX;
			RotY= (fits_info.header).RotY;
			Epoch= (fits_info.header).Epoch;
		}//close SetFITSCards()

		/**
		* \brief Get current WCS type
		*/
		std::string GetWCSType(){return m_wcsType;}
		
		/**
		* \brief Get world coordinate system 
		*/
		WorldCoor* GetWorldCoord(int coordSystem=-1){

			//Compute the wcs from vars
			WorldCoor* wcs= wcskinit(Nx,Ny,(char*)CoordTypeX.c_str(),(char*)CoordTypeY.c_str(),Cx,Cy,Xc,Yc,NULL,dX,dY,RotY,(int)(Epoch),Epoch);
			std::string wcsType= std::string(getwcsout(wcs));
			cout<<"ImgMetaData::GetWorldCoord(): INFO: wcsType="<<wcsType<<endl;

			//Convert wcs to desired type
			char* flag = (char*)("");
			if(coordSystem==eGALACTIC)
				flag = (char*)("GALACTIC");	
			else if(coordSystem==eJ2000)
				flag = (char*)("FK5");
			else if(coordSystem==eB1950)
				flag = (char*)("FK4");
			else if(coordSystem==-1 && m_wcsType!="")					
				flag = (char*)(m_wcsType.c_str());
			
			if(strcmp(flag,"")!=0) {
				wcsoutinit (wcs,flag);
				m_wcsType= std::string(flag);
			}
			
			wcsType= std::string(getwcsout(wcs));
			cout<<"ImgMetaData::GetWorldCoord(): INFO: wcsType="<<wcsType<<endl;

			return wcs;
		}//close GetWorldCoord()

		/**
		* \brief Get size of synthetic beam
		*/
		double GetBeamSize(){
			double beamSize= sqrt(fabs(Bmaj*Bmin));
			return beamSize;	
		}
		/**
		* \brief Get pixel scale
		*/
		int GetBeamSizeInPixel(){
			double beamSize= GetBeamSize();
			double pixScale= sqrt(fabs(dX*dY));
			int npixInBeam= int(ceil(beamSize/pixScale));
			return npixInBeam;
		}
		/**
		* \brief Get flux correction from beam
		*/
		double GetBeamFluxIntegral(){
			double fx= Bmaj;
			double fy= Bmin;
			double A= TMath::Pi()*fx*fy/(4*log(2));//2d gaussian area with FWHM=fx,fy
			return A;
		}

	public:
		//Image size
		int Nx;
		int Ny;

		//Reference pixel id
		int Cx;
		int Cy;

		//Reference pixel coords
		double Xc;
		double Yc;

		//Pixel size
		double dX;
		double dY;

		//System rotation info
		double RotX;
		double RotY;

		//Type of astro coords
		std::string CoordTypeX;
		std::string CoordTypeY;

		//Units
		std::string BUnit;

		//Beam info
		double Bmaj;
		double Bmin;
		double Bpa;

		//Obs Epoch
		double Epoch;

	private:
		std::string m_wcsType;
		
	ClassDef(ImgMetaData,3)

};

#ifdef __MAKECINT__
#pragma link C++ class ImgMetaData+;
//#pragma link C++ enum WCSType+;
#endif


class ImgStats : public TObject {

	public:
		ImgStats(){
			Reset();		
		};
		virtual ~ImgStats(){};

	public:
		void Reset(){
			n= 0;
			min= 0;
			max= 0;
			mean= 0;
			meanErr= 0;	
			rms= 0;
			rmsErr= 0;
			skewness= 0;
			skewnessErr= 0;
			median= 0;
			medianRMS= 0;
			bwLocation= 0;
			bwScale= 0;
			bwLocationIter= 0;
			bwScaleIter= 0;
			clippedMedian= 0;
			clippedRMS= 0;
		}
		void Log(std::string level="INFO"){
			LOG(level,GetPrintable());
		}
		void Print(){
			cout<<"*** IMG STATS ***"<<endl;
			cout<<"N="<<n<<" min/max="<<min<<"/"<<max<<endl;
			cout<<"Mean: "<<mean<<" +- "<<meanErr<<endl;
			cout<<"RMS: "<<rms<<" +- "<<rmsErr<<endl;
			cout<<"Skewness: "<<skewness<<" +- "<<skewnessErr<<endl;
      cout<<"Median: "<<median<<", MAD: "<<medianRMS<<endl;
			cout<<"BiWeight Location: "<<bwLocation<<", Scale: "<<bwScale<<endl;		
			cout<<"Clipped Median: "<<clippedMedian<<" MAD: "<<clippedRMS<<endl;
			cout<<"*****************"<<endl;
		}
		std::string GetPrintable(){
			std::stringstream ss;
			ss<<"IMG STATS: ";
			ss<<"N="<<n<<" min/max="<<min<<"/"<<max<<", ";
			ss<<"Mean: "<<mean<<" +- "<<meanErr<<", ";
			ss<<"RMS: "<<rms<<" +- "<<rmsErr<<", ";
			ss<<"Skewness: "<<skewness<<" +- "<<skewnessErr<<", ";
			ss<<"Median: "<<median<<", MAD: "<<medianRMS<<", ";
			ss<<"BiWeight Location: "<<bwLocation<<", Scale: "<<bwScale<<", ";
			ss<<"Clipped Median: "<<clippedMedian<<" MAD: "<<clippedRMS;
			return ss.str();
		}

	public:			
		int n;
		double min;
		double max;
		double mean;
		double meanErr;
		double rms;
		double rmsErr;
  	double skewness;	
		double skewnessErr;			
		double median;
		double medianRMS;
		double bwLocation;
		double bwScale;
		double bwLocationIter;
		double bwScaleIter;
		double clippedMedian;
		double clippedRMS;
		
	ClassDef(ImgStats,2)

};

#ifdef __MAKECINT__
#pragma link C++ class ImgStats+;
#endif


class Img : public TH2F {

  public:
		
    /** 
		\brief Class constructor: initialize structures.
 		*/
		Img();
    Img(const Img& img);
		Img(const char* name, const char* title, int nbinsx, const float* xbins, int nbinsy, const float* ybins);
		Img(const char* name, const char* title, int nbinsx, float xlow, float xup, int nbinsy, float ylow, float yup);
		Img& operator=(const Img &img);
		virtual void Copy(TObject &hnew) const;

		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~Img();

		enum ImgType {
			eUnknown=0,
			eFluxMap=1,
			eMeanFluxMap=2,
			eSignificanceMap=3,
			eBkgMap=4,
			eNoiseMap=5,
			eBinaryMap=6,
			eResidualMap=7,
			ePullMap=8,
			eCurvatureMap=9
		};
		enum ImgSmoothFilter {
			eGaus=1,
			eGuided=2,
			eWT=3,
			eGradient=4,
			eLaplacian=5,
			eLoG=6,	
			eNormLoG=7
		};

	public:
				
		//== Fill methods ==
		/**
		* \brief Fill pixels (to be used to compute stats at fill time)
		*/
		int FillPixel(double x,double y,double w,bool useNegativePixInStats=true);
		/**
		* \brief Fill pixels (to be used to compute stats at fill time)
		*/
		void FillFromMat(cv::Mat,bool useNegativePixInStats=true);

		/**
		* \brief Check if given bin id is within image range
		*/
		bool HasBin(int binId){
			return (!this->IsBinOverflow(binId) && !this->IsBinUnderflow(binId) );
		}
		bool HasBin(int binIdX,int binIdY){
			int binId= this->GetBin(binIdX,binIdY);
			return HasBin(binId);
		}
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
		bool IsBinContentInRange(int binIdX,int binIdY,double minThr,double maxThr){	
			if(!HasBin(binIdX,binIdY)) return false;
			double w= this->GetBinContent(binIdX,binIdY);
			return (w>=minThr && w<=maxThr);
		}

		//== Read methods ==
		/**
		* \brief Read image from FITS file
		*/
		int ReadFITS(std::string filename,int ix_min=-1,int ix_max=-1,int iy_min=-1,int iy_max=-1);
		/**
		* \brief Write image to FITS file
		*/
		int WriteFITS(std::string outfilename);
		/**
		* \brief Get image subregion or tile
		*/
		Img* GetTile(int ix_min,int ix_max,int iy_min,int iy_max);

		
		//== Stats methods ==
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

		//== Meta-Data methods ==
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
		/**
		* \brief Get pixel world coordinate
		*/
		int GetWorldCoord(int ix,int iy,double& xpos, double& ypos) {
			return Caesar::AstroUtils::PixelToWCSCoords(this,ix,iy,xpos,ypos);
		}
		

		//== Bkg finder methods ==
		/**
		* \brief Compute local bkg
		*/
		BkgData* ComputeBkg(int estimator=BkgFinder::eMedianBkg,bool computeLocalBkg=true,int boxSizeX=100,int boxSizeY=100, double gridStepSizeX=10, double gridStepSizeY=10, bool use2ndPass=true,bool skipOutliers=false,double seedThr=5,double mergeThr=2.6,int minPixels=10);
		/**
		* \brief Compute significance map
		*/
		Img* GetSignificanceMap(BkgData* bkgData,bool useLocalBkg=false);

		//== Thresholding methods ==
		/**
		* \brief Compute pixel histo
		*/
		TH1D* GetPixelHisto(int nbins=100,bool normalize=false);
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
		double FindMedianThreshold(double thrFactor);
		/**
		* \brief Find optimal global threshold
		*/
		double FindOptimalGlobalThreshold(double thrFactor,int nbins=100,bool smooth=true);
		/**
		* \brief Get binarized image
		*/
		Img* GetBinarizedImage(double threshold,double fgValue=1,bool isLowerThreshold=false);


		//== Source finding methods ==
		/**
		* \brief Find compact sources
		*/
		int FindCompactSource(std::vector<Source*>&,Img* floodImg=0,BkgData* bkgData=0,double seedThr=5,double mergeThr=2.6,int minPixels=10,bool findNegativeExcess=false,bool mergeBelowSeed=false,bool findNestedSources=false,double nestedBlobThrFactor=1);
		/**
		* \brief Find nested sources
		*/
		int	FindNestedSource(std::vector<Source*>& sources,BkgData* bkgData=0,int minPixels=5,double nestedBlobThreshold=1);

		/**
		* \brief Find extended sources with ChanVese method
		*/
		int FindExtendedSource_CV(std::vector<Source*>&,BkgData* bkgData=0,int minPixels=10,bool findNegativeExcess=false,double dt=0.1,double h=1,double lambda1=1.0,double lambda2=2.0,double mu=0.5,double nu=0,double p=1);

		/**
		* \brief Find extended sources with Hierarchical Clustering method
		*/
		int FindExtendedSource_HClust(std::vector<Source*>&,Img* saliencyImg,Img* edgeImg);


		//== Filter methods ==
		/**
		* \brief Get guided filter image
		*/
		Img* GetGuidedFilterImage(int radius=12,double eps=0.04);
		/**
		* \brief Get image wavelength decomposition
		*/
		std::vector<Img*> GetWaveletDecomposition(int nScales);
		/**
		* \brief Get kirsch image
		*/
		Img* GetKirschImage();
		/**
		* \brief Get image gradient
		*/
		Img* GetGradientImage(bool invert=false);
		/**
		* \brief Get image laplacian
		*/
		Img* GetLaplacianImage(bool invert=false);
		/**
		* \brief Get laplacian of gaussian image
		*/
		Img* GetLoGImage(bool invert=false);	
		/**
		* \brief Get scale-normalized laplacian of gaussian image
		*/
		Img* GetNormLoGImage(int size=3,double scale=1,bool invert=false);
		/**
		* \brief Smooth image
		*/
		Img* GetSmoothedImage(int size_x=3,int size_y=3,double sigma_x=1,double sigma_y=1);
		/**
		* \brief Smooth image
		*/
		Img* GetSaliencyMap(int resoMin=20,int resoMax=60,int resoStep=10,double beta=1,int minRegionSize=10,double knnFactor=0.2,double spatialRegFactor=6,bool useRobust=true,bool addCurvDist=true,double salientMultiplicityThrFactor=0.7,bool addBkgMap=true,bool addNoiseMap=true,BkgData* bkgData=0,double saliencyThrFactor=2,double imgThrFactor=1);

		//== Convert/Scale util methods
		/**
		* \brief Get matrix from image
		*/
		TMatrixD* GetMatrix();
		/**
		* \brief Convert image to OpenCV mat float image
		*/
		cv::Mat ImgToMat(std::string encoding="64");	
		/**
		* \brief Get normalized image
		*/
		Img* GetNormalizedImage(std::string normScale="LINEAR",int normmin=1,int normmax=256, bool skipEmptyBins=true);
		/**
		* \brief Clone image
		*/
		Img* GetCloned(std::string name,bool copyMetaData=true,bool resetStats=true){
			//cout<<"Img::GetCloned(): Start..."<<endl;
			//Img* clone= (Img*)this->Clone(name.c_str());
			Img* clone= new Img;
			*clone= *this;
			clone->SetNameTitle(name.c_str(),name.c_str());			
			//cout<<"Img::GetCloned(): After clone..."<<endl;
			//if(copyMetaData) clone->CopyMetaData(m_MetaData);	
			if(!copyMetaData) clone->ClearMetaData();
			//cout<<"Img::GetCloned(): After copy meta data..."<<endl;			
			if(resetStats) clone->ResetImgStats(true,true);
			//cout<<"Img::GetCloned(): After reset stats..."<<endl;
			return clone;
		}//close GetCloned()
		/**
		* \brief Get masked image
		*/
		Img* GetMask(Img* mask,bool isBinary=false);
		/**
		* \brief Get source masked image
		*/
		Img* GetSourceMask(std::vector<Source*>const& sources,bool isBinary=false,bool invert=false);
		/**
		* \brief Returns a residual image obtained by dilating given sources with a random background
		*/
		Img* GetSourceResidual(std::vector<Source*>const& sources,int KernSize=5,int dilateModel=MorphFilter::eDilateWithBkg,int dilateSourceType=-1,bool skipToNested=false,BkgData* bkgData=0,bool useLocalBkg=false,bool randomize=false,double zThr=20);

		//== Drawing methods
		/**
		* \brief Draw image
		*/
		int Plot(std::vector<Source*>const&,bool useCurrentCanvas=true,bool drawFull=false,int paletteStyle=Caesar::eRAINBOW,bool drawColorPalette=true,bool putWCAxis=false,int coordSystem=-1,std::string units="Jy/beam");
		/**
		* \brief Set draw range
		*/
		void SetDrawRange(double zmin,double zmax);

	private:
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
		void UpdateMoments(int ix,int iy,double w);
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

		

	private:

		//Meta-data
		bool m_HasMetaData;
		ImgMetaData* m_MetaData;
		
		//Stats data
		bool m_HasStats;
		ImgStats* m_Stats;
		
		int m_Npix;//npixels	
		double m_M1;//1st moments
		double m_M2;//2nd moment
		double m_M3;//3rd moment
		double m_M4;//4th moment
		double m_PixelMin;
		double m_PixelMax;
		

		ClassDef(Img,1)

};

#ifdef __MAKECINT__
#pragma link C++ enum Img::ImgType+;
#pragma link C++ class StatsData+;
#pragma link C++ class MetaData+;
#pragma link C++ class Img+;
#pragma link C++ class vector<Img>+;
#pragma link C++ class vector<Img*>+;
#endif

}//close namespace

#endif
