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
* @file SLICData.cc
* @class SLICData
* @brief SLIC data class
*
* Superpixel data
* @author S. Riggi
* @date 20/01/2015
*/

#include <SLICData.h>
#include <Region.h>
#include <Img.h>

#include <TObject.h>
#include <TMatrixD.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <numeric>
#include <time.h>
#include <ctime>
#include <queue>

using namespace std;

ClassImp(Caesar::SLICData)
ClassImp(Caesar::SLICContourData)
ClassImp(Caesar::SLICNeighborData)
ClassImp(Caesar::SLICSimilarityData)

namespace Caesar {

SLICData::SLICData() {
	inputImg= 0;
	edgeImg= 0;
	laplImg= 0;
	regions.clear();
	labels.clear();
}//close costructor

SLICData::~SLICData(){
	Clear();
}//close destructor

void SLICData::Clear(){

	ClearImages();
	ClearRegions();

}//close ClearAll() 


void SLICData::ClearImages(){
	if(inputImg) inputImg->Delete();
	if(edgeImg) edgeImg->Delete();
	if(laplImg) laplImg->Delete();
}//close ClearImages()

void SLICData::ClearRegions(){
	
	for(unsigned int i=0;i<regions.size();i++){
		if(regions[i]){
			delete regions[i];
			regions[i]= 0;
		}
	}
	regions.clear();
	labels.clear();

}//close ClearRegions()

}//close namespace


