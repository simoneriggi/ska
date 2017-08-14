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
* @file BkgFinder.h
* @class BkgFinder
* @brief BkgFinder
*
* Class for computing local background data in images
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef BkgFinder_h
#define BkgFinder_h 1

#include <Consts.h>
#include <TObject.h>

namespace Caesar {

class Img;
class Image;
class BkgSampleData;
class BkgData;
class ImgBkgData;

class BkgFinder : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    BkgFinder();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~BkgFinder();

		//enum BkgEstimator {eMeanBkg=1,eMedianBkg=2,eBiWeightBkg=3,eMedianClippedBkg= 4};
		//enum BkgMethod {eGridBkg=1,eSuperpixelBkg=2};

	public:
	
		//=========================================
		//==  NEW IMAGE METHODS 
		//=========================================
		/**
		* \brief Find image background
		*/
		static ImgBkgData* FindBkg(Image* img,int estimator=eMedianBkg,bool computeLocalBkg=true,int boxSizeX=100,int boxSizeY=100, double gridStepSizeX=10, double gridStepSizeY=10, bool use2ndPass=true,bool skipOutliers=false,double seedThr=5,double mergeThr=2.6,int minPixels=10);
	
		
		//=========================================
		//==  OLD IMAGE METHODS 
		//=========================================
		/**
		* \brief Find image background
		*/
		static BkgData* FindBkg(Img* img,int estimator=eMedianBkg,bool computeLocalBkg=true,int boxSizeX=100,int boxSizeY=100, double gridStepSizeX=10, double gridStepSizeY=10, bool use2ndPass=true,bool skipOutliers=false,double seedThr=5,double mergeThr=2.6,int minPixels=10);
		
		
	private:

		//=========================================
		//==  NEW IMAGE METHODS 
		//=========================================
		/**
		* \brief Compute global background
		*/
		static int ComputeGlobalBkg(ImgBkgData* bkgData,Image* img,int estimator);
		/**
		* \brief Compute background of sample/tile image
		*/
		static int ComputeSampleBkg(BkgSampleData& bkgSampleData,Image* img,int estimator,long int ix_min=-1,long int ix_max=-1,long int iy_min=-1,long int iy_max=-1);
		/**
		* \brief Find local grid background
		*/
		static int FindLocalGridBkg(ImgBkgData* bkgData,Image* img,int estimator,long int boxSizeX,long int boxSizeY,double gridStepSizeX,double gridStepSizeY,bool use2ndPass);
		/**
		* \brief Compute local grid background
		*/
		static int ComputeLocalGridBkg(ImgBkgData* bkgData,Image* img,int estimator,long int boxSizeX,long int boxSizeY,double gridStepSizeX,double gridStepSizeY);
		
		//=========================================
		//==  OLD IMAGE METHODS 
		//=========================================
		/**
		* \brief Compute global background
		*/
		static int ComputeGlobalBkg(BkgData* bkgData,Img* img,int estimator);
		/**
		* \brief Compute background of sample/tile image
		*/
		static int ComputeSampleBkg(BkgSampleData& bkgSampleData,Img* img,int estimator);
		/**
		* \brief Find local grid background
		*/
		static int FindLocalGridBkg(BkgData* bkgData,Img* img,int estimator,int boxSizeX, int boxSizeY, double gridStepSizeX,double gridStepSizeY,bool use2ndPass);
		/**
		* \brief Compute local grid background
		*/
		static int ComputeLocalGridBkg(BkgData* bkgData,Img* img,int estimator,int boxSizeX,int boxSizeY,double gridStepSizeX,double gridStepSizeY);
		
	private:

		ClassDef(BkgFinder,1)

};

#ifdef __MAKECINT__
#pragma link C++ class BkgFinder+; 
#endif

}//close namespace 


#endif


