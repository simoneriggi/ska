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
* @file SaliencyFilter.h
* @class SaliencyFilter
* @brief Class implementing saliency filtering
*
* Saliency Filter
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef SaliencyFilter_h
#define SaliencyFilter_h 1

//ROOT
#include <TObject.h>

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


namespace Caesar{

class Region;
class Img;
class BkgData;

class SaliencyFilter : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    SaliencyFilter();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~SaliencyFilter();

	public:
		/**
		* \brief Compute saliency map for one resolution
		*/
		//static Img* ComputeSaliencyMap(Img* img,int reso,double regFactor,int minRegionSize,double knnFactor,double spatialRegFactor,bool useRobust,bool addCurvDist);
		static Img* ComputeSaliencyMap(Img* img,int reso=20,double regFactor=1,int minRegionSize=10,double knnFactor=1,bool useRobust=false,double expFalloffPar=100,double distanceRegPar=1);
		
		/**
		* \brief Compute multi resolution saliency map 
		*/
		//static Img* ComputeMultiResoSaliencyMap(Img* img,int resoMin=20,int resoMax=60,int resoStep=10,double beta=1,int minRegionSize=10,double knnFactor=0.2,double spatialRegFactor=6,bool useRobust=true,bool addCurvDist=true,double salientMultiplicityThrFactor=0.7,bool addBkgMap=true,bool addNoiseMap=true,BkgData* bkgData=0,double saliencyThrFactor=2,double imgThrFactor=1);
		static Img* ComputeMultiResoSaliencyMap(Img* img,int resoMin=20,int resoMax=60,int resoStep=10,double beta=1,int minRegionSize=10,double knnFactor=1,bool useRobustPars=false,double expFalloffPar=100,double distanceRegPar=1,double salientMultiplicityThrFactor=0.7,bool addBkgMap=true,bool addNoiseMap=true,BkgData* bkgData=0,double saliencyThrFactor=2,double imgThrFactor=1);

	private:
		
		//static Img* ComputeSaliencyMap(Img* img,std::vector<Region*>const& regions,double knnFactor,double spatialRegFactor,bool useRobust,bool addCurvDist);
		static Img* ComputeSaliencyMap(Img* img,std::vector<Region*>const& regions,double knnFactor,bool useRobust,double expFalloffPar,double distanceRegPar);

	private:

	ClassDef(SaliencyFilter,1)

};

#ifdef __MAKECINT__
#pragma link C++ class SaliencyFilter+;
#endif

}//close namespace

#endif

