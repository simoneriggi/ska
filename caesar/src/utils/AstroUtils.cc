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
* @file AstroUtils.cc
* @class AstroUtils
* @brief Utility functions for astronomical tasks
*
* Utility functions for astronomical tasks
* @author S. Riggi
* @date 15/01/2016
*/


#include <AstroUtils.h>
#include <Img.h>

#include <TObject.h>



#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <ctime>

using namespace std;

ClassImp(Caesar::AstroUtils)

namespace Caesar {

AstroUtils::AstroUtils(){

}

AstroUtils::~AstroUtils(){

}



int AstroUtils::PixelToWCSCoords(Caesar::Img* image,WorldCoor* wcs,double ix,double iy,double& xpos, double& ypos) {

	//Check pixel values in input
	if(!image){
		cerr<<"AstroUtils::PixelToWCSCoords(): ERROR: Null image ptr given!"<<endl;
		return -1;	
	}

	//Get image range
	int Nx= image->GetNbinsX();
	int Ny= image->GetNbinsY();
	double xmin= image->GetXaxis()->GetXmin();
	double ymin= image->GetYaxis()->GetXmin();
	double xmax= image->GetXaxis()->GetXmax();
	double ymax= image->GetYaxis()->GetXmax();

	if(ix<xmin || iy<ymin || ix>xmax || iy>ymax ){
		cerr<<"AstroUtils::PixelToWCSCoords(): ERROR: Invalid pix range selected (ix="<<ix<<", iy="<<iy<<")"<<endl;
		return -1;	
	}

	//Check WCS
	if(!wcs){
		cerr<<"AstroUtils::PixelToWCSCoords(): ERROR: Null ptr to given WCS!"<<endl;
		return -1;
	}

	//Convert coords
	pix2wcs (wcs,ix,iy,&xpos, &ypos);

	return 0;

}//close PixelToWCSCoords()


int AstroUtils::PixelToWCSCoords(Caesar::Img* image,double ix,double iy,double& xpos, double& ypos,int coordSystem) {

	//Check pixel values in input
	if(!image){
		cerr<<"AstroUtils::PixelToWCSCoords(): ERROR: Null image ptr given!"<<endl;
		return -1;	
	}

	//Get image range
	int Nx= image->GetNbinsX();
	int Ny= image->GetNbinsY();
	double xmin= image->GetXaxis()->GetXmin();
	double ymin= image->GetYaxis()->GetXmin();
	double xmax= image->GetXaxis()->GetXmax();
	double ymax= image->GetYaxis()->GetXmax();

	if(ix<xmin || iy<ymin || ix>xmax || iy>ymax ){
		cerr<<"AstroUtils::PixelToWCSCoords(): ERROR: Invalid pix range selected (ix="<<ix<<", iy="<<iy<<")"<<endl;
		return -1;	
	}

	
	//Check image meta-data
	if(!image->HasMetaData() ){
    cerr<<"AstroUtils::PixelToWCSCoords(): WARNING: No metadata available in image!"<<endl;
		return -1;
	}
	Caesar::ImgMetaData* metadata= image->GetMetaData();	
	
	//Get the coord system
	WorldCoor* wcs= metadata->GetWorldCoord(coordSystem);
	if(!wcs){
		cerr<<"AstroUtils::PixelToWCSCoords: WARNING: Failed to get WorldCoord system from metadata!"<<endl;
		return -1;
	}

	//Convert coords
	pix2wcs (wcs,ix,iy,&xpos, &ypos);

	//Clear up
	delete wcs;
	wcs= 0;

	return 0;
		
}//close PixelToWCSCoords()







}//close namespace



