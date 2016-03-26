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
* @file BlobFinder.h
* @class BlobFinder
* @brief Blob finder class
*
* Class to perform blob finding 
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef BlobFinder_h
#define BlobFinder_h 1

#include <TObject.h>
#include <TMatrixD.h>

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


namespace Caesar {


class Img;
class BkgData;
class Blob;
class Source;


class BlobFinder : public TObject {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		BlobFinder();
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~BlobFinder();

	public:
		static int FloodFill(Img* img,std::vector<int>& clusterPixelIds,int seedPixelId,double floodMinThr,double floodMaxThr);
		
		template <class T>
			static int FindBlobs(Img* inputImg,std::vector<T*>& blobs,Img* floodImg=0,BkgData* bkgData=0,double seedThr=5,double mergeThr=2.6,int minPixels=10,bool findNegativeExcess=false,bool mergeBelowSeed=false);

		static Img* GetMultiScaleBlobMask(Img* img,int kernelFactor,double sigmaMin,double sigmaMax,double sigmaStep,int thrModel=2,double thrFactor=1);

	public:
		
		ClassDef(BlobFinder,1)

};//close BlobFinder()



#ifdef __MAKECINT__
#pragma link C++ class BlobFinder+;
#pragma link C++ function FindBlob<Blob>(Caesar::Img*,std::vector<Blob*>&,Caesar::Img*,Caesar::BkgData*,double,double,int,bool,bool);
#pragma link C++ function FindBlob<Source>(Caesar::Img*,std::vector<Source*>&,Caesar::Img*,Caesar::BkgData*,double,double,int,bool,bool);
#endif

}//close namespace

#endif

