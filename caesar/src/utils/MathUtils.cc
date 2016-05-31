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
* @file MathUtils.cc
* @class MathUtils
* @brief Utility functions for math tasks
*
* Utility functions for math tasks
* @author S. Riggi
* @date 15/01/2016
*/


#include <MathUtils.h>
#include <CodeUtils.h>
#include <Logger.h>

#include <TMath.h>

#include <linterp.h>

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

ClassImp(Caesar::MathUtils)

namespace Caesar {

MathUtils::MathUtils(){

}

MathUtils::~MathUtils(){

}

int MathUtils::BiLinearInterpolation(std::vector<double>const& sampled_gridX, std::vector<double>const& sampled_gridY,std::vector<double>const& sampledZ,std::vector<double>const& interp_gridX,std::vector<double>const& interp_gridY,std::vector<double>& interpZ){

	//## Check args
	long int nSamplesX= (long int)sampled_gridX.size();
	long int nSamplesY= (long int)sampled_gridY.size();
	long int nSamplesZ= (long int)sampledZ.size();
	if(nSamplesX<=0 || nSamplesY<=0){
		//cerr<<"MathUtils::BiLinearInterpolation(): ERROR: Invalid sample grid size given!"<<endl;
		ERROR_LOG("Invalid sample grid size given!");
		return -1;
	}
	long int num_elements = nSamplesX*nSamplesY;
	if(nSamplesZ!=num_elements){
		//cerr<<"MathUtils::BiLinearInterpolation(): ERROR: Invalid sample Z given (it must be equal to Nx x Ny and indexed as ix*Ny+iy)!"<<endl;
		ERROR_LOG("Invalid sample Z given (it must be equal to Nx x Ny and indexed as ix*Ny+iy)!");
		return -1;
	}

	long int nInterpX= (long int)interp_gridX.size();
	long int nInterpY= (long int)interp_gridY.size();
  long int num_interp_elements = nInterpX*nInterpY;
	if(nInterpX<=0 || nInterpY<=0){
		//cerr<<"MathUtils::BiLinearInterpolation(): ERROR: Invalid interpolation grid given (size must be >0 in both directions)!"<<endl;
		ERROR_LOG("Invalid interpolation grid given (size must be >0 in both directions)!");
		return -1;
	}
	interpZ.clear();

	//## Perform the 2D interpolation
	try {		
	
		// Construct the grid in each dimension (note that we will pass in a sequence of iterators pointing to the beginning of each grid)
		//cout<<"MathUtils::BiLinearInterpolation(): INFO: Build 2D grid for interpolation (nSamplesX="<<nSamplesX<<", nSamplesY="<<nSamplesY<<")..."<<endl;
		DEBUG_LOG("Build 2D grid for interpolation (nSamplesX="<<nSamplesX<<", nSamplesY="<<nSamplesY<<")...");

  	std::vector< std::vector<double>::const_iterator > grid_iter_list;
  	grid_iter_list.push_back(sampled_gridX.begin());
  	grid_iter_list.push_back(sampled_gridY.begin());
  
  	// the size of the grid in each dimension
  	array<int,2> grid_sizes;
  	grid_sizes[0] = nSamplesX;
  	grid_sizes[1] = nSamplesY;
  
  	// construct the interpolator. the last two arguments are pointers to the underlying data
		//cout<<"MathUtils::BiLinearInterpolation(): INFO: Build the bkg interpolator..."<<endl;
		DEBUG_LOG("Build the bkg interpolator...");
  	InterpMultilinear<2, double> interpolator_ML(grid_iter_list.begin(), grid_sizes.begin(), sampledZ.data(), sampledZ.data() + num_elements);
		
		
		//Construct interpolated grid
		// interpolate multiple values: create sequences for each coordinate
		//cout<<"MathUtils::BiLinearInterpolation(): INFO: Build the interpolated grid ("<<nInterpX<<","<<nInterpY<<") x("<<interp_gridX[0]<<","<<interp_gridX[nInterpX-1]<<") y("<<interp_gridY[0]<<","<<interp_gridY[nInterpY-1]<<")..."<<endl;
		DEBUG_LOG("Build the interpolated grid ("<<nInterpX<<","<<nInterpY<<") x("<<interp_gridX[0]<<","<<interp_gridX[nInterpX-1]<<") y("<<interp_gridY[0]<<","<<interp_gridY[nInterpY-1]<<")...");
  	
  	//std::vector<double> interp_gridX = CodeUtils::linspace(xlim[0],xlim[1], Nx);
		//std::vector<double> interp_gridY = CodeUtils::linspace(ylim[0],ylim[1], Ny);
		

  	std::vector<double> interp_x(num_interp_elements);
  	std::vector<double> interp_y(num_interp_elements);
  	for (unsigned int i=0; i<interp_gridX.size(); i++) {
    	for (unsigned int j=0; j<interp_gridY.size(); j++) {
				long int gBin= i*interp_gridY.size() + j;
	  		interp_x[gBin] = interp_gridX[i];
	  		interp_y[gBin] = interp_gridY[j];
			}
  	}
 	 	interpZ.assign(num_interp_elements,0);

		// pass in a sequence of iterators, one for each coordinate
  	std::vector< std::vector<double>::iterator > interp_list;
  	interp_list.push_back(interp_x.begin());
  	interp_list.push_back(interp_y.begin());
  
		//Interpolate sequence
		//cout<<"MathUtils::BiLinearInterpolation(): INFO: Run the interpolation on grid..."<<endl;
		DEBUG_LOG("Run the interpolation on grid...");
		interpolator_ML.interp_vec(num_interp_elements, interp_list.begin(), interp_list.end(), interpZ.begin());
  	
	}//close try block
	catch( std::exception &ex ) {
		//cerr << "MathUtils::BiLinearInterpolation(): ERROR: Exception detected in interpolation: " << ex.what() << endl;
		ERROR_LOG("Exception detected in interpolation (err=" << ex.what()<<")");
		return -1;
  } 
	catch(...) { 
		//cerr << "MathUtils::BiLinearInterpolation(): ERROR: C++ exception (unknown reason) in interpolation!" << endl;
		ERROR_LOG("Unknown exception caught in interpolation!");
		return -1;
  }		
	
	return 0;

}//close BiLinearInterpolation()


}//close namespace



