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
* @file Pixel.cc
* @class Pixel
* @brief Pixel data class
*
* Pixel class
* @author S. Riggi
* @date 20/01/2015
*/

#include <Pixel.h>

#include <TObject.h>

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

ClassImp(Caesar::Pixel)

namespace Caesar {

Pixel::Pixel() {
	
	Init();
	
}//close costructor


Pixel::~Pixel(){
	

}//close destructor


void Pixel::Init(){

	//Init values
	id= -1;
	type= eNormal;
	S= 0; 
	S_curv= 0; 
	S_edge= 0;
	x= -1; 
	y=-1; 
	ix= -1; 
	iy= -1; 
	isOnEdge= false; 
	distanceToEdge= std::numeric_limits<double>::infinity();
	bkgLevel= 0;
	noiseLevel= 0;

}//close Init()

}//close namespace

