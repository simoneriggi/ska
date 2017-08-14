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
* @file BkgData.h
* @class BkgData
* @brief BkgData
*
* Class for storing bkg data
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef BkgData_h
#define BkgData_h 1

#include <Img.h>
#include <Image.h>
#include <CodeUtils.h>
#include <Logger.h>

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

using namespace std;


namespace Caesar {

class BkgSampleData : public TObject{

	public:

		BkgSampleData(){
			ix_min= 0; ix_max= 0;
			iy_min= 0; iy_max= 0;
			npix= 0;
			isReliable= true;
			bkgLevel= 0;
			bkgRMS= 0;
		}
		virtual ~BkgSampleData(){};

	public:
		void CopyBkgData(BkgSampleData aBkgSample){
			bkgLevel= aBkgSample.bkgLevel;
			bkgRMS= aBkgSample.bkgRMS;
		}
		void Log(std::string level="INFO"){
			LOG(level,GetPrintable());
		}
		void Print(){
			cout<<"== BKG SAMPLE DATA NO. "<<id<<" =="<<endl;
			cout<<"N="<<npix<<" xrange("<<ix_min<<","<<ix_max<<") yrange("<<iy_min<<","<<iy_max<<")"<<endl;
			cout<<"bkgLevel="<<bkgLevel<<" bkgRMS="<<bkgRMS<<endl;
			cout<<"=================================="<<endl;
		}
		std::string GetPrintable(){
			std::stringstream ss;
			ss<<"BkgSample no. "<<id<<": ";
			ss<<"N="<<npix<<", xrange("<<ix_min<<","<<ix_max<<"), yrange("<<iy_min<<","<<iy_max<<"), ";
			ss<<"bkgLevel="<<bkgLevel<<", bkgRMS="<<bkgRMS;
			return ss.str();
		}

	public:	
		int id;
		int ix_min;
		int iy_min;
		int ix_max;
		int iy_max;
		int npix;
		bool isReliable;
		double bkgLevel;
		double bkgRMS;	

	ClassDef(BkgSampleData,1)
};

//=====================================
//==   NEW BKG DATA
//=====================================

class ImgBkgData : public TObject {

	public:
		ImgBkgData();
		virtual ~ImgBkgData();

	public:
		void ClearSamplings(){
			BkgSamplings.clear();
		}
		void ClearBkgMap(){
			if(!BkgMap) return;
			delete BkgMap;
			BkgMap= 0;
		}
		void ClearNoiseMap(){
			if(!NoiseMap) return;
			delete NoiseMap;
			NoiseMap= 0;
		}

		void Clear(){
			ClearSamplings();
			ClearBkgMap();
			ClearNoiseMap();
		}
		void CopyBkgMap(Caesar::Image* aMap);		
		void CopyNoiseMap(Caesar::Image* aMap);
		bool HasLocalBkg(){return (BkgMap && NoiseMap);}
		
	public:
		std::vector<BkgSampleData> BkgSamplings;
		Image* BkgMap;//the interpolated bkg map
		Image* NoiseMap;//the interpolated noise map
		double gBkg;
		double gNoise;

	ClassDef(ImgBkgData,1)

};//close class ImgBkgData


//=====================================
//==   OLD BKG DATA
//=====================================

class BkgData : public TObject {

	public:
		BkgData();
		virtual ~BkgData();

	public:
		void ClearSamplings(){
			BkgSamplings.clear();
		}
		void ClearBkgMap(){
			if(!BkgMap) return;
			delete BkgMap;
			BkgMap= 0;
		}
		void ClearNoiseMap(){
			if(!NoiseMap) return;
			delete NoiseMap;
			NoiseMap= 0;
		}

		void Clear(){
			ClearSamplings();
			ClearBkgMap();
			ClearNoiseMap();
		}
		void CopyBkgMap(Caesar::Img* aMap);		
		void CopyNoiseMap(Caesar::Img* aMap);
		bool HasLocalBkg(){return (BkgMap && NoiseMap);}
		
	public:
		std::vector<BkgSampleData> BkgSamplings;
		Img* BkgMap;//the interpolated bkg map
		Img* NoiseMap;//the interpolated noise map
		double gBkg;
		double gNoise;

	ClassDef(BkgData,1)

};//close class

#ifdef __MAKECINT__
#pragma link C++ class BkgSampleData+;
#pragma link C++ class vector<BkgSampleData>+;
#pragma link C++ class BkgData+;
#pragma link C++ class ImgBkgData+;
#endif


}//close namespace 


#endif


