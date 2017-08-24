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
* @file SLICUtils.h
* @class SLICData
* @brief SLIC utils class
*
* Superpixel utility functions
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _SLIC_UTILS_h
#define _SLIC_UTILS_h 1

#include <Contour.h>

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

class Image;
class Region;
class SLICContourData;
class SLICSimilarityData;
class SLICNeighborData;
class SLICNeighborCollection;
class SLICData;

class SLICUtils : public TObject {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SLICUtils();
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~SLICUtils();

	public:
		
		/** 
		\brief Compute superpixel boundary contours
 		*/
		static SLICContourData* ComputeBoundaryContours(SLICData* slicData);

		/** 
		\brief Compute region similarities
 		*/
		static SLICSimilarityData* ComputeRegionSimilarity(SLICData* slicData,std::vector<SLICNeighborCollection>& neighbors,double beta=0.5);

		/** 
		\brief Find superpixel neighbors
 		*/
		static int FindNeighbors(std::vector<SLICNeighborCollection>& neighbors,SLICData* slicData,SLICContourData* contourData,bool get2ndNeighbors=true,int selectedTag=-1,bool includeSpatialDist=false);

		/** 
		\brief Compute segmented image given a list of tagged regions
 		*/
		static Image* GetSegmentedImage(Image* img,std::vector<Region*>const& regions,int selectedTag=-1,bool normalize=false,bool binarize=false);

		/** 
		\brief Count number of regions per tag
 		*/
		static int CountTaggedRegions(std::vector<Region*>const& regions,int& NSig,int& NBkg,int& NUntagged);		

		/** 
		\brief Tag regions into signal/bkg according to signal & bkg marker images
 		*/
		static int TagRegions(std::vector<Region*>& regions,Image* binaryMap_bkg,Image* binaryMap_signal);

	public:
		
		ClassDef(SLICUtils,1)

};//close SLICUtils()



#ifdef __MAKECINT__
#pragma link C++ class SLICUtils+;
#endif

}//close namespace

#endif

