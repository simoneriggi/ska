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
* @file SLIC.h
* @class SLIC
* @brief SLIC generator class
*
* Superpixel generator
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _SLIC_h
#define _SLIC_h 1

#include <SLICData.h>
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


namespace Caesar {

class Image;

class SLIC : public TObject {

  public:
		
    /** 
		\brief Class constructor: initialize structures.
 		*/
		SLIC();
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~SLIC();


	public:


		/** 
		\brief Generate a superpixel partition given the passed options
 		*/
		static SLICData* SPGenerator(Image* img,int regionSize=10,double regParam=1,int minRegionSize=10,bool useLogScaleMapping=false,Image* edgeImg=0);

		
	private:
		/** 
		\brief Initialize the superpixel data structure
 		*/
		static int SetSPData(SLICData* slicData,Image* img,bool useLogScaleMapping,Image* edgeImg);

		
	private:

		ClassDef(SLIC,1)

};

#ifdef __MAKECINT__
#pragma link C++ class SLIC+;
#endif

}//close namespace

#endif

