/**
* @file Source.h
* @class Source
* @brief Source
*
* Source
* @author S. Riggi
* @date 20/01/2015
*/



#ifndef Source_h
#define Source_h 1

#include <Contour.h>
//#include <SLICSegmenter.h>
//#include <SLICUtils.h>

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

#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <complex>

class Img;

using namespace std;

class Source {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    Source();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~Source();

		enum SourceType {eCompact=1,ePointLike=2,eExtended=3,eExtendedSegm=4,eUnknown=5};
		enum PixelType {eSeed=1,eNormal=2,eHalo=3};
		enum SourceFlag {eReal=1,eCandidate=2,eFake=3};
		enum SourceImageMode {eBinarizedMap=0,eFluxMap=1,eSignificanceMap=2,eSourceSignificanceMap=3,eFluxCurvMixtureMap=4,eFluxCurvMixtureSignificanceMap=5,eFluxCurvMap=6,eMeanFluxMap=7};
		
		struct Pixel {	
			double S;//pixel intensity
			PixelType Type;//pixel flag
			double Z;//pixel significance in reference image
			double Curv;//pixel curvature value in reference image
			double BkgLevel;//estimate of pixel background level
			double NoiseLevel;//estimate of pixel noise level
			int id;//global bin id of reference image
			int ix;//pixel id x
			int iy;//pixel id y
			double x;//pixel x coordinate
			double y;//pixel y coordinate
		};

		struct FitInfo {
			double A;
			double AErr;
			double Cx;
			double Cy;
			double CxErr;
			double CyErr;
			double sigmaX;
			double sigmaY;
			double sigmaXErr;
			double sigmaYErr;
			double theta;
			double thetaErr;
			double ellMaj;
			double ellMin;
		};

		typedef std::vector<Pixel> PixelCollection;
		typedef std::map<int,Pixel> PixelMap;

		typedef std::vector<Source> SourceCollection;
		typedef std::vector<Source*> SourcePtrCollection;


	public:

		//Setters & getters	
		void SetType(SourceType choice){fType=choice;}
		void SetFlag(SourceFlag choice){fFlag=choice;}
		void SetId(int id){fId=id;}
		void SetName(std::string name){fName=name;}

		void SetNPixels(int value){fNPix=value;}	
		bool IsAtEdge(){return fHasPixelsAtEdge;}	
		void SetEdgeFlag(bool choice){fHasPixelsAtEdge=choice;}

		//Pixel manipulation
		PixelCollection GetPixels(){return fPixelCollection;}
		void AddPixel(Pixel pixel);
		int ComputeStats();
		int ComputeMorphologyParams();
		int ComputeZernikeMoments(int order=6);

		std::vector<Contour*> GetContours(){return fContourCollection;}

		/**
		* \brief Add source fit
		*/
		void AddFitInfo(FitInfo* info){
			if(!info) return;
			fFitInfo.push_back(info);
			fHasFitInfo= true;
		}

		/**
		* \brief Check if this is a good source
		*/
		bool IsGoodSource();
		/**
		* \brief Check if this is a compact/point-source candidate
		*/
		bool IsCompactPointLike();

		/**
		* \brief Check if the sources has multicomponent inside
		*/
		int Deblend(double curvThr=0,int componentMinNPix=6);

		/**
		* \brief Return source image
		*/
		Img* GetImage(SourceImageMode mode=eFluxMap,double curvMixtureFract=0);

		/**
		* \brief Find sub source inside this source
		*/
		int FindNestedSource(double seedThr=5,double mergeThr=2.5,int nPixMin=10,bool findPointLike=true);
		int FindNestedSource(double curvThr=0,int nPixMin=6,double peakThr=5);
		/**
		* \brief Draw source
		*/
		void Draw(SourceImageMode mode=eFluxMap,double curvMixtureFract=0,bool drawNested=false,bool drawOnlyNested=false,bool drawContour=true,bool drawEllipse=false,bool drawBoundingBox=false,bool drawLegends=false);

		/**
		* \brief Dump source info
		*/
		void Dump(bool dumpContourInfo=false,bool dumpNestedSourceInfo=false){
			cout<<"*** SOURCE ID: "<<fId<<" ***"<<endl;
			cout<<"Name: "<<fName<<" Type: "<<fType<<" HasNestedSources? "<<fHasNestedSources<<" AtEdge? "<<fHasPixelsAtEdge<<endl;
			cout<<"C("<<fX0<<","<<fY0<<" Cw("<<fWx0<<","<<fWy0<<") Stot="<<fS<< " Smin/Smax="<<fSmin<<"/"<<fSmax<<endl;
			cout<<"N= "<<fNPix<<" M1="<<fM1<<", M2="<<fM2<<", M3="<<fM3<<", M4="<<fM4<<endl;
			cout<<"Mean="<<fMean<<" RMS="<<fRMS<<" Skewness="<<fSkewness<<", Kurtosis="<<fKurtosis<<" Median="<<fMedian<<" MedianRMS="<<fMedianRMS<<endl;
			cout<<"Xmin/Xmax="<<fXmin<<"/"<<fXmax<<", Ymin/Ymax="<<fYmin<<"/"<<fYmax<<endl;
			if(dumpContourInfo){
				for(unsigned int i=0;i<fContourCollection.size();i++){
					if(fContourCollection[i]) fContourCollection[i]->Dump();
				}
			}
			if(dumpNestedSourceInfo && fHasNestedSources){
				for(unsigned int i=0;i<fNestedSourceCollection.size();i++){
					if(fNestedSourceCollection[i]) fNestedSourceCollection[i]->Dump(dumpContourInfo,dumpNestedSourceInfo);
				}
			}
			cout<<"****************************"<<endl;
		}

		/**
		* \brief Dump ellipse source info
		*/
		std::string DumpDS9EllipseRegionInfo(bool dumpNestedSourceInfo=false){
			//ellipse x y radius radius angle
			std::stringstream sstream;
			sstream<<"ellipse ";
			for(unsigned int i=0; i<fContourCollection.size(); i++){ 
				if(!fContourCollection[i]->HasEllipseFit) continue;
				double EllX= fContourCollection[i]->EllipseCenter.x;
				double EllY= fContourCollection[i]->EllipseCenter.y;
				double EllMajAxis= fContourCollection[i]->EllipseMajAxis;
				double EllMinAxis= fContourCollection[i]->EllipseMinAxis;
				double EllRotAxis= fContourCollection[i]->EllipseRotAngle-90;//-90 comes from DS9 strange format!
				cv::Point2f BBoxCenter= fContourCollection[i]->BoundingBoxCenter;
				double BBoxMinAxis=  fContourCollection[i]->BoundingBoxMin;	
				double BBoxMajAxis= fContourCollection[i]->BoundingBoxMaj;
				double BBoxAngle= fContourCollection[i]->BoundingBoxAngle;	
				sstream<<EllX+1<<" "<<EllY+1<<" "<<(EllMajAxis/2)<<" "<<(EllMinAxis/2)<<" "<<(EllRotAxis)<<" ";
				//sstream<<BBoxCenter.x+1<<" "<<BBoxCenter.y+1<<" "<<(BBoxMajAxis/2)<<" "<<(BBoxMinAxis/2)<<" "<<(BBoxAngle+90)<<" ";
			}
			sstream<<"# text={S"<<fId<<"}";

			if(dumpNestedSourceInfo && fHasNestedSources){			
				sstream<<endl;
				for(unsigned int k=0;k<fNestedSourceCollection.size();k++){
					sstream<<"ellipse ";
					std::vector<Contour*> nestedContours= fNestedSourceCollection[k]->fContourCollection;
					for(unsigned int i=0; i<nestedContours.size(); i++){ 
						if(!nestedContours[i]->HasEllipseFit) continue;
						double EllX= nestedContours[i]->EllipseCenter.x;
						double EllY= nestedContours[i]->EllipseCenter.y;
						double EllMajAxis= nestedContours[i]->EllipseMajAxis;
						double EllMinAxis= nestedContours[i]->EllipseMinAxis;
						double EllRotAxis= nestedContours[i]->EllipseRotAngle-90;//-90 comes from DS9 strange format!
						cv::Point2f BBoxCenter= fContourCollection[i]->BoundingBoxCenter;
						double BBoxMinAxis=  fContourCollection[i]->BoundingBoxMin;	
						double BBoxMajAxis= fContourCollection[i]->BoundingBoxMaj;
						double BBoxAngle= fContourCollection[i]->BoundingBoxAngle;	
						sstream<<EllX+1<<" "<<EllY+1<<" "<<(EllMajAxis/2)<<" "<<(EllMinAxis/2)<<" "<<(EllRotAxis)<<" ";
						//sstream<<BBoxCenter.x+1<<" "<<BBoxCenter.y+1<<" "<<(BBoxMajAxis/2)<<" "<<(BBoxMinAxis/2)<<" "<<(BBoxAngle)<<" ";
					}//end loop contours
					sstream<<"# text={S"<<fId<<"_Nest"<<k<<"}";
					if(k!=fNestedSourceCollection.size()-1) sstream<<endl;
				}//end loop nested sources
			}//close dumpNestedSourceInfo

			return sstream.str();
		}//close DumpDS9EllipseRegionInfo()


		std::string DumpDS9RegionInfo(bool dumpNestedSourceInfo=false){
			//global color=red dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 image
			std::stringstream sstream;
			sstream<<"polygon ";
			for(unsigned int i=0; i<fContourCollection.size(); i++){ 
				int nPoints= fContourCollection[i]->GetN();
				for(int j=0;j<nPoints;j++){
					cv::Point2f* contPnt= fContourCollection[i]->GetPoint(j);
					if(!contPnt) continue;
					sstream<<(int)contPnt->x+1<<" "<<(int)contPnt->y+1<<" ";
				}
			}
			//sstream<<"# text={S"<<fName<<"-"<<fId<<"}";
			sstream<<"# text={S"<<fId<<"}";

			if(dumpNestedSourceInfo && fHasNestedSources){			
				sstream<<endl;
				for(unsigned int k=0;k<fNestedSourceCollection.size();k++){
					sstream<<"polygon ";
					//std::string nestedSourceRegionInfo= fNestedSourceCollection[k]->DumpDS9RegionInfo(false);
					//sstream<<nestedSourceRegionInfo;
					std::vector<Contour*> nestedContours= fNestedSourceCollection[k]->fContourCollection;
					for(unsigned int i=0; i<nestedContours.size(); i++){ 
						int nPoints= nestedContours[i]->GetN();
						for(int j=0;j<nPoints;j++){
							cv::Point2f* contPnt= nestedContours[i]->GetPoint(j);
							if(!contPnt) continue;
							sstream<<(int)contPnt->x+1<<" "<<(int)contPnt->y+1<<" ";
						}
					}//end loop contours
					sstream<<"# text={S"<<fId<<"_Nest"<<k<<"}";
					if(k!=fNestedSourceCollection.size()-1) sstream<<endl;
				}//end loop nested sources
			}//close dumpNestedSourceInfo

			return sstream.str();
		}

		/**
		* \brief Check if source overlap with a given source
		*/
		bool IsOverlapping(Source* aSource,double& overlapArea,double overlapThreshold=0.9);
		/**
		* \brief Is source inside given source
		*/
		bool IsInsideSource(Source* aSource){
			if(!aSource) return false;
			bool isInsideX= (fXmin>=aSource->fXmin && fXmax<=aSource->fXmax);
			bool isInsideY= (fYmin>=aSource->fYmin && fYmax<=aSource->fYmax);
			bool isInside= (isInsideX && isInsideY);			
			return isInside;
		}
		/**
		* \brief Check if a given point is inside the source
		*/
		bool IsPointInsideSource(double x,double y,double distTolerance=2);

		/**
		* \brief Check if source has computed parameters
		*/
		bool HasParameters(){return fHasParameters;}

		/**
		* \brief Reset nested sources
		*/
		void ResetNestedSources(){
			for(unsigned int k=0;k<fNestedSourceCollection.size();k++){
				if(fNestedSourceCollection[k]) {
					delete fNestedSourceCollection[k];
					fNestedSourceCollection[k]= 0;
				}
			}
			fNestedSourceCollection.clear();
			fNestedSourceCollection.resize(0);
			if(fNestedSource){
				delete fNestedSource;
				fNestedSource= 0;
			}
			fHasNestedSources= false;
		}//close ResetNestedSources()

		/**
		* \brief Add nested sources
		*/
		void AddNestedSource(Source* aNestedSource){
			if(!aNestedSource) return;
			int nNestedSources= fNestedSourceCollection.size();
			int nestedId= nNestedSources+1;
			TString nestedName= Form("%s_N%d",fName.c_str(),nestedId);
			fNestedSource= new Source;
			*(fNestedSource)= *(aNestedSource);//copy source
			fNestedSource->fName= std::string(nestedName);
			fNestedSource->fId= nestedId;
			fNestedSource->fType= aNestedSource->fType;
			fNestedSource->fDepthLevel= this->fDepthLevel+1;
			fNestedSourceCollection.push_back(fNestedSource);
			fHasNestedSources= true;
		}	

		/**
		* \brief Find nested source with hierarchical clustering
		*/
		int FindNestedClusterRegions(int regionSize=10,double regularization=100, int minRegionSize=5,double SPMergingRatio=0.3,double SPMergingRegularization=0.1,double MaxDissRatio=1.15);

	private:

		
	public:

		std::string fName;
		long int fId;
		int fFlag;
		int fType;
		int fNPix;	
		bool fIsGoodSource;
		bool fIsSignificant;
		bool fIsSalient;

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

		double fS;//sum of pixel signals
		double fSmax;//max of pixel signals
		double fSmin;//min of pixel signals
		double fSxx;
		double fSyy;
		double fSxy;
		double fFluxCorrection;//correction factor to go from Jy/beam to Jy

		double fCurvMin;//min of pixel curvature
		double fCurvMax;//max of pixel curvature

		int fPixIdmax;//id of pixel with max signal
		int fPixIdmin;//id of pixel with min signal
		double fF;//sum of pixel fluxes
		double fFlux;//estimated source flux (i.e. 2D gaussian fit)
		double fWx0;//Signal-weighted X position average
		double fWy0;//Signal-weighted Y position average
		double fX0;//X position average
		double fY0;//Y position average

		//Bounding box
		double fXmin;
		double fXmax;
		double fYmin;
		double fYmax;
		int fIx_min;
		int fIx_max;
		int fIy_min;
		int fIy_max;
	
		Contour* fContour;
		std::vector<Contour*> fContourCollection;

		bool fHasPixelsAtEdge;
		PixelCollection fPixelCollection;
		bool fHasParameters;
		//std::map<int,int> fPixelIndexMap;//map key: , map value: 
	
		//Nested source info
		int fDepthLevel;
		bool fHasNestedSources;
		Source* fNestedSource;
		std::vector<Source*> fNestedSourceCollection;
		
		//Fit info	
		bool fHasFitInfo;
		std::vector<FitInfo*> fFitInfo;
		double fFitChi2;
		double fFitNDF;
		double fFitNFreePars;

		//2D morphological pars
		cv::Moments moments;
		double HuMoments[7];
		std::vector<double> ZMMoments;
		

};

#ifdef __MAKECINT__
#pragma link C++ class Source+;
#pragma link C++ class vector<Source>+;
#pragma link C++ class vector<Source*>+;
#endif


#endif
