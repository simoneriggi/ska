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
* @file MorphFilter.h
* @class MorphFilter
* @brief Class implementing morphological filtering
*
* Morphological Filter
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef MorphFilter_h
#define MorphFilter_h 1

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

using namespace std;


namespace Caesar{

class Source;
class BkgData;
class Img;

class MorphFilter : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    MorphFilter();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~MorphFilter();

		enum DilationModel {eDilateWithBkg=1,eDilateWithSourceMedian=2};

	public:
	
		static Img* Dilate(Img* img,int KernSize,bool returnPeakImg=false);
		static int DilateAroundSources(Img* img,std::vector<Source*>const& sources,int KernSize=5,int dilateModel=eDilateWithBkg,int dilateSourceType=-1,bool skipToNested=false,BkgData* bkgData=0,bool useLocalBkg=false,bool randomize=false);

	private:
		static int DilateAroundSource(Img* img,Source* source,int KernSize,int dilateModel,int dilateSourceType,bool skipToNested,BkgData* bkgData,bool useLocalBkg,bool randomize);
		static int FindDilatedSourcePixels(Img* img,Source* source,int KernSize,std::vector<int>& pixelsToBeDilated);
		
	private:

	ClassDef(MorphFilter,1)

};

#ifdef __MAKECINT__
#pragma link C++ class MorphFilter+;
#endif

}//close namespace

#endif
