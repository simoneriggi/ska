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
* @file Pixel.h
* @class Pixel
* @brief Pixel data class
*
* Pixel class
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef Pixel_h
#define Pixel_h 1

#include <TObject.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

#include <functional>
#include <utility>
#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <iostream>
#include <time.h>
#include <ctime>

using namespace std;


namespace Caesar {

class Blob;

class Pixel : public TObject {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		Pixel();
		/**
		* \brief Copy constructor
		*/
		Pixel(const Pixel& pixel);
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~Pixel();

		/**
		* \brief Assignment Operator
		*/
		Pixel& operator=(const Pixel& pixel);
		/**
		* \brief Equality Operator
		*/
		bool operator==(const Pixel& pixel) const {
  		return (this->x == pixel.x && this->y == pixel.y);
		}
		/**
		* \brief < Operator
		*/
		bool operator<(const Pixel& pixel) const {
			return std::tie (this->x,this->y) < std::tie (pixel.x,pixel.y);
  	}

		/**
		* \brief Copy method
		*/
		void Copy(TObject& pixel) const;

		enum PixelType {eNormal=1,eSeed=2,eHalo=3};

	public: 
		void SetPhysCoords(double xx,double yy){x=xx;y=yy;}
		void SetCoords(long int i,long int j){ix=i;iy=j;}
		void SetBkg(double bkg,double noise){bkgLevel=bkg;noiseLevel=noise;}
		void SetCurv(double val){S_curv=val;}
		void SetEdge(double val){S_edge=val;}
		std::pair<double,double> GetBkg(){return std::make_pair(bkgLevel,noiseLevel);}
		double GetCurv(){return S_curv;}
		double GetEdge(){return S_edge;}
		

	private:
		void Init();
		void UpdateMoments(Pixel* pixel);
		void ResetMoments();
		void ResetStats();
		void ResetPixels();

	public:
		long int id;//global bin id of reference image	
		PixelType type;//pixel flag
		double S;//pixel intensity
		double x;//pixel x coordinate
		double y;//pixel y coordinate
		long int ix;//pixel id x
		long int iy;//pixel id y
		bool isOnEdge;//flag marking if pixel is found on region contour
		double distanceToEdge;//distance to the edge (=0 for edge pixels)		

	private:
		//Set private to force external class to explicitly call setter methods
		double S_curv;//curvature estimator
		double S_edge;//edge estimator
					
		double bkgLevel;//estimate of pixel background level
		double noiseLevel;//estimate of pixel noise level

	friend class Blob;

	ClassDef(Pixel,1)

};
typedef std::vector<Pixel*> PixelCollection;
typedef std::map<int,Pixel*> PixelMap;

struct PixelMatcher {
	bool operator()(const Pixel* lhs, const Pixel* rhs) const { 
		return std::tie (lhs->x,lhs->y) < std::tie (rhs->x,rhs->y);
	}
	static bool AreAdjacent(const Pixel* lhs, const Pixel* rhs) {
		double distX= rhs->x - lhs->x;
		double distY= rhs->y - lhs->y;
		bool areAdjacent= (fabs(distX)<=1 && fabs(distY)<=1);
  	return areAdjacent;
	}
};//close PixelMatcher()



#ifdef __MAKECINT__
#pragma link C++ class Pixel+;
#pragma link C++ class vector<Pixel*>+;
#pragma link C++ class map<int,Pixel*>+;
#pragma link C++ enum Pixel::PixelType+;
#endif

}//close namespace


#endif

