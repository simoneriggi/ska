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
* @file AstroUtils.h
* @class AstroUtils
* @brief Utility functions for astronomical tasks
*
* Utility functions for astronomical tasks
* @author S. Riggi
* @date 15/01/2016
*/


#ifndef AstroUtils_h
#define AstroUtils_h 1

#include <wcs.h>

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
#include <time.h>
#include <ctime>

using namespace std;


namespace Caesar {

class Img;
class ImgMetaData;
class ImgStatsData;

class AstroUtils : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    AstroUtils();
		/**
		* \brief Class destructor: free allocated memory
		*/
   	virtual ~AstroUtils();

		
	public:

		static int PixelToWCSCoords(Caesar::Img* image,WorldCoor* wcs,double ix,double iy,double& xpos, double& ypos);
		static int PixelToWCSCoords(Caesar::Img* image,double ix,double iy,double& xpos, double& ypos,int coordSystem=-1);
		
		
	private:
	
		ClassDef(AstroUtils,1)
};

#ifdef __MAKECINT__
#pragma link C++ class AstroUtils+;
#endif	

}//close namespace


#endif 
