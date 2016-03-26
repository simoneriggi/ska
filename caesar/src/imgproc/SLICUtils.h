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

#ifndef SLICUtils_h
#define SLICUtils_h 1

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


class Img;
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
		
		static SLICContourData* ComputeBoundaryContours(SLICData* slicData);
		//static SLICSimilarityData* ComputeRegionSimilarity(SLICData* slicData,SLICContourData* contourData,double beta=0.5,bool includeSpatialDist=false,int mergedTag=-1);
		static SLICSimilarityData* ComputeRegionSimilarity(SLICData* slicData,std::vector<SLICNeighborCollection>& neighbors,double beta=0.5);

		static int FindNeighbors(std::vector<SLICNeighborCollection>& neighbors,SLICData* slicData,SLICContourData* contourData,bool get2ndNeighbors=true,int selectedTag=-1,bool includeSpatialDist=false);

		static Img* GetSegmentedImage(Img* img,std::vector<Region*>const& regions,int selectedTag=-1,bool normalize=false,bool binarize=false);
		static int CountTaggedRegions(std::vector<Region*>const& regions,int& NSig,int& NBkg,int& NUntagged);		
		static int TagRegions(std::vector<Region*>& regions,Img* binaryMap_bkg,Img* binaryMap_signal);

	public:
		
		ClassDef(SLICUtils,1)

};//close SLICUtils()



#ifdef __MAKECINT__
#pragma link C++ class SLICUtils+;
#endif

}//close namespace

#endif

