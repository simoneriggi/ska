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
* @file BkgData.cc
* @class BkgData
* @brief BkgData
*
* Class for storing bkg data
* @author S. Riggi
* @date 20/01/2015
*/

#include <BkgData.h>
#include <Img.h>

#include <TObject.h>
#include <TString.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
using namespace std;


ClassImp(Caesar::BkgSampleData)
ClassImp(Caesar::BkgData)

namespace Caesar {

BkgData::BkgData() {
	BkgSamplings.clear();
	BkgMap= 0;
	NoiseMap= 0;
	gBkg= 0;
	gNoise= 0;
}//close costructor

BkgData::~BkgData() {
	Clear();
}//close destructor

void BkgData::CopyBkgMap(Caesar::Img* aMap){
	if(!aMap) return;
	TString mapName= "bkgMap";
	if(BkgMap) {
		mapName= BkgMap->GetName();
		delete BkgMap;
		BkgMap= 0;
	}
	BkgMap= aMap->GetCloned(std::string(mapName),true,true);
}//close CopyBkgMap()
		
void BkgData::CopyNoiseMap(Caesar::Img* aMap){
	if(!aMap) return;
	TString mapName= "noiseMap";
	if(NoiseMap) {
		mapName= NoiseMap->GetName();	
		delete NoiseMap;
		NoiseMap= 0;
	}
	NoiseMap= aMap->GetCloned(std::string(mapName),true,true);

}//close CopyNoiseMap()
		
}//close namespace

