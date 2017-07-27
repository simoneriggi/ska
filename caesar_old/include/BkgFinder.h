/**
* @file BkgFinder.h
* @class BkgFinder
* @brief BkgFinder
*
* BkgFinder
* @author S. Riggi
* @date 20/01/2015
*/



#ifndef BkgFinder_h
#define BkgFinder_h 1

#include <Img.h>

#include <RInside.h>

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
#include <TGraph2D.h>


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


class BkgFinder {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    BkgFinder();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~BkgFinder();

		struct BkgMapData {	
			std::vector<Img::BkgData*> bkgData;
			Img* BkgMap;
			Img* RMSMap;
			BkgMapData() {
				BkgMap= 0;
				RMSMap= 0;
				bkgData.clear();
			}
			~BkgMapData() {
				if(BkgMap){
					delete BkgMap;
					BkgMap= 0;
				}
				if(RMSMap){
					delete RMSMap;
					RMSMap= 0;
				}
				for(unsigned int i=0;i<bkgData.size();i++){
					if(bkgData[i]){
						delete bkgData[i];
						bkgData[i]= 0;
					}
				}//end loop 
				bkgData.clear();
			}//close destructor 
		};

	public:
	
		/**
		* \brief Estimate bkg from an image (tile collection)
		*/
		BkgFinder::BkgMapData* FindBkg(Img* img,Img::BkgMethod method, int boxSizeX, int boxSizeY, double boxSlideOffsetX=1, double boxSlideOffsetY=1);
		BkgFinder::BkgMapData* FindGridBkg(Img* img,Img::BkgMethod method, int boxSizeX, int boxSizeY, double gridStepSizeX=1, double gridStepSizeY=1);
		BkgFinder::BkgMapData* FindSuperpixelBkg(Img* img,Img::BkgMethod method, int boxSizeX, int boxSizeY, double gridStepSizeX=1, double gridStepSizeY=1);

		/**
		* \brief Compute the background estimate for a tile image
		*/
		Img::BkgData* ComputeBkg(Img* img, Img::BkgMethod method);

		void SetSegmentationRegionSize(int s){fRegionSize=s;}
		void SetSegmentationRegularization(double s){fRegularization=s;}
		void SetSegmentationMinRegionArea(int s){fMinRegionArea=s;}
		void SetSegmentationDistanceEps(double s){fColorEps=s;}
	

	private:
		
		

		void Init();

		Img::BkgData* GetMeanBkg(Img*);
		Img::BkgData* GetMedianBkg(Img*);
		Img::BkgData* GetBiWeightBkg(Img*);
		Img::BkgData* GetMedianClippedBkg(Img*);
		Img::BkgData* GetRobustBkg(Img*,int regionSize=10,double regularization=100,int minRegionArea=5,double colorEps=5);
		Img::BkgData* GetSimpleRobustBkg(Img*,int regionSize=10,double regularization=100,int minRegionArea=5,double colorEps=5);
		Img::BkgData* GetDummyRobustBkg(Img*);
		std::vector<double> linspace(double first, double last, int len);

	private:

		RInside* fR;

		int nTiles;
		int nTilesX;
		int nTilesY;

		int npix;
		double mean;
		double rms;
		double median;
		double medianRMS;
		double bkgLevel;		
		double bkgRMS;
		int tileId;
		int tileMinX;
		int tileMaxX;
		int tileMinY;
		int tileMaxY;

		int fRegionSize;
		double fRegularization;
		int fMinRegionArea;
		double fColorEps;

			
		ClassDef(BkgFinder,1)

};

#ifdef __MAKECINT__
#pragma link C++ class BkgFinder+; 
#endif

#endif
