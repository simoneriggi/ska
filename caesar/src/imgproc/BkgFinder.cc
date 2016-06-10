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
* @file BkgFinder.cc
* @class BkgFinder
* @brief BkgFinder
*
* Class for computing local background data in images
* @author S. Riggi
* @date 20/01/2015
*/

#include <BkgFinder.h>
#include <BkgData.h>
#include <Source.h>
#include <Img.h>
#include <CodeUtils.h>
#include <MathUtils.h>
#include <Logger.h>


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
using namespace std;

ClassImp(Caesar::BkgFinder)

namespace Caesar {


BkgFinder::BkgFinder() {

}//close costructor


BkgFinder::~BkgFinder() {

}//close destructor


BkgData* BkgFinder::FindBkg(Img* img,int estimator,bool computeLocalBkg,int boxSizeX,int boxSizeY, double gridStepSizeX,double gridStepSizeY,bool use2ndPass,bool skipOutliers,double seedThr,double mergeThr,int minPixels)
{
	//## Check input image
	if(!img) {
		ERROR_LOG("Null ptr to given image!");
		return 0;
	}
		
	//## Compute stats
	if(!img->HasStats()){//computing stats
		if( img->ComputeStats(true,false,false)<0 ){
			ERROR_LOG("Failed to compute stats!");
			return 0;	
		}
	}//close if

	//## Init 
	BkgData* bkgData= new BkgData;

	//## Compute global bkg
	INFO_LOG("Computing global bkg...");
	if(ComputeGlobalBkg(bkgData,img,estimator)<0){
		ERROR_LOG("Failed to compute global bkg!");
		delete bkgData;
		bkgData= 0;
		return 0;
	}

	//## Compute local bkg?		
	if(computeLocalBkg){
		INFO_LOG("Computing local bkg ...");
		int status= FindLocalGridBkg(bkgData,img,estimator,boxSizeX,boxSizeY,gridStepSizeX,gridStepSizeY,use2ndPass);
		if(status<0){
			ERROR_LOG("Failed to compute local grid bkg!");
			delete bkgData;
			bkgData= 0;
			return 0;
		}
		DEBUG_LOG("Local bkg computation completed!");
	}//close if computeLocalBkg

	//## Skip outliers?
	//## Search and exclude significant blobs (both positive & negative excesses) 
	//## using the first estimate bkg/noise
	if(skipOutliers){
		INFO_LOG("Improving bkg estimate by skipping outliers ...");

		//Get significance map
		Img* significanceMap= img->GetSignificanceMap(bkgData,computeLocalBkg);
		if(!significanceMap){
			ERROR_LOG("Failed to compute the significance map (needed to exclude blobs)!");
			delete bkgData;
			bkgData= 0;
			return 0;
		}	

		//Find blobs
		DEBUG_LOG("Finding compact blobs to be tagged as outliers...");
		std::vector<Source*> blobs;
		bool findNegativeExcess= true;
		bool mergeBelowSeed= false;
		bool findNestedSources= false;
		int status= img->FindCompactSource(blobs,significanceMap,bkgData,seedThr,mergeThr,minPixels,findNegativeExcess,mergeBelowSeed,findNestedSources);
		if(status<0){
			ERROR_LOG("Failed to find significant blobs!");
			delete bkgData;
			bkgData= 0;
			significanceMap->Delete();
			return 0;
		}

		//Find image without outliers (set to zero)
		DEBUG_LOG("Computing image without outliers (set to zero)...");
		Img* img_wOutliers= img->GetSourceMask(blobs,false,true);//invert mask
		if(!img_wOutliers){
			ERROR_LOG("Failed to compute image with blob outliers subtracted!");
			delete bkgData;
			bkgData= 0;
			significanceMap->Delete();
			for(unsigned int i=0;i<blobs.size();i++){
				if(blobs[i]){
					delete blobs[i];
					blobs[i]= 0;
				}
			}
			blobs.clear();
			return 0;
		}//close if

		//Recompute bkg on residual map (using this function recursively)
		//Do not skip outliers this time!
		DEBUG_LOG("Recomputing bkg on residual map...");
		BkgData* robustBkgData= FindBkg(img_wOutliers,estimator,computeLocalBkg,boxSizeX,boxSizeY,gridStepSizeX,gridStepSizeY,use2ndPass,false);
		if(!robustBkgData){
			ERROR_LOG("Failed to compute bkg over image with blob outliers subtracted!");
			delete bkgData;
			bkgData= 0;
			significanceMap->Delete();
			for(unsigned int i=0;i<blobs.size();i++){
				if(blobs[i]){
					delete blobs[i];
					blobs[i]= 0;
				}
			}
			blobs.clear();
			return 0;
		}

		//Override main bkgData with robust estimates
		for(unsigned int i=0;i<(robustBkgData->BkgSamplings).size();i++) {
			(bkgData->BkgSamplings)[i]= (robustBkgData->BkgSamplings)[i];
		}
		bkgData->CopyBkgMap(robustBkgData->BkgMap);//copy new bkg map (delete previous)
		bkgData->CopyNoiseMap(robustBkgData->NoiseMap);//copy new noise map (delete previous)
		bkgData->gBkg= robustBkgData->gBkg;
		bkgData->gNoise= robustBkgData->gNoise;
				
		//Delete stuff
		significanceMap->Delete();
		for(unsigned int i=0;i<blobs.size();i++){
			if(blobs[i]){
				delete blobs[i];
				blobs[i]= 0;
			}
		}
		blobs.clear();
		delete robustBkgData;
		robustBkgData= 0;

	}//close if skip outliers

	return bkgData;

}//close FindBkg()


int BkgFinder::FindLocalGridBkg(BkgData* bkgData,Img* img,int estimator,int boxSizeX, int boxSizeY, double gridStepSizeX, double gridStepSizeY,bool use2ndPass){

	//## Compute bkg data
	INFO_LOG("Computing local bkg (1st pass)...");
	if(ComputeLocalGridBkg(bkgData,img, estimator, boxSizeX, boxSizeY, gridStepSizeX,gridStepSizeY)<0){	
		ERROR_LOG("Computation of local background failed for this image!");
		return -1;
	}

	//## Improve rms by recomputing stuff from residual map 
	if(use2ndPass){
		INFO_LOG("Improving rms estimation with a 2nd pass...");

		TString residualMapName= Form("%s_residual",img->GetName());
		Img* residualMap= img->GetCloned(std::string(residualMapName),true,true);
		residualMap->Add(bkgData->BkgMap,-1);//subtract the bkg level model

		//Compute bkg for residual map
		BkgData* bkgData_residual= new BkgData;
		if(ComputeLocalGridBkg(bkgData_residual,residualMap,estimator,boxSizeX, boxSizeY, gridStepSizeX,gridStepSizeY)<0){
			ERROR_LOG("Computation of local background failed @ 2nd pass!");
			delete residualMap;
			residualMap= 0;
			delete bkgData_residual;
			bkgData_residual= 0;
			return -1;
		}
	
		//Update rms data in image data
		DEBUG_LOG("Update bkg RMS sampling data in image data...");
		for(unsigned int i=0;i<(bkgData_residual->BkgSamplings).size();i++) {
			(bkgData->BkgSamplings)[i].bkgRMS= (bkgData_residual->BkgSamplings)[i].bkgRMS;
		}

		DEBUG_LOG("Copying new noise map to bkg data...");
		bkgData->CopyNoiseMap(bkgData_residual->NoiseMap);//copy new noise map (delete previous)

		//Delete stuff
		DEBUG_LOG("Deleting allocated data...");
		delete residualMap;
		residualMap= 0;
		delete bkgData_residual;
		bkgData_residual= 0;
		DEBUG_LOG("End local bkg computation");
	}//close if

	return 0;

}//close FindLocalGridBkg()


int BkgFinder::ComputeLocalGridBkg(BkgData* bkgData,Img* img,int estimator,int boxSizeX, int boxSizeY, double gridStepSizeX, double gridStepSizeY){

	//## Check input image
	if(!img) {
		ERROR_LOG("Null ptr to given image!");
		return -1;
	}

	//## Check options
	int Nx= img->GetNbinsX();
	int Ny= img->GetNbinsY();
	if(boxSizeX<=0 || boxSizeX>=Nx || boxSizeY<=0 || boxSizeY>=Ny) {
		ERROR_LOG("Invalid box size given (too small or larger than image size)!");
		return -1;
	}
	if(gridStepSizeX<=0 || gridStepSizeY<=0 ){
		ERROR_LOG("Invalid grid step size given (null or negative)!");
		return -1;
	}
	
	int TileSizeX= boxSizeX;
	int TileSizeY= boxSizeY;
	double xlim[2]= {img->GetXaxis()->GetXmin(),img->GetXaxis()->GetXmax()};
	double ylim[2]= {img->GetYaxis()->GetXmin(),img->GetYaxis()->GetXmax()};
	
	DEBUG_LOG("N("<<Nx<<","<<Ny<<") TileSize=("<<TileSizeX<<","<<TileSizeY<<") GridStepSize("<<gridStepSizeX<<","<<gridStepSizeY<<")");
	
	//## Initialize and count number of rows & cols
	int indexX_start= TileSizeX/2;
	int indexY_start= TileSizeY/2;
	int indexX_end= Nx-1-TileSizeX/2;
	int indexY_end= Ny-1-TileSizeY/2;
	int indexX= indexX_start;
	int indexY= indexY_start;
	int nTiles= 0;
	int nTilesX= 0;
	int nTilesY= 0;
	std::vector<int> ixList;
	std::vector<int> iyList;
	DEBUG_LOG("Counting number of tiles from ("<<indexX<<","<<indexY<<") up to ("<<indexX_end<<","<<indexY_end<<")...");

	while(indexY<=indexY_end){
		iyList.push_back(indexY);
		nTilesY++;
		indexY+= gridStepSizeY;
	}//end while loop

	while(indexX<=indexX_end){
		ixList.push_back(indexX);	
		nTilesX++;
		indexX+= gridStepSizeX;
	}
	nTiles= nTilesX*nTilesY;

	DEBUG_LOG("nTilesX="<<nTilesX<<"("<<ixList.size()<<") nTilesY="<<nTilesY<<"("<<iyList.size()<<") nTiles="<<nTiles);

	//## Loop over all number of tiles and compute bkg info for each one
	DEBUG_LOG("Computing tile background...");	

	//Vectors to store sampling bkg x & y coordinates
	std::vector<double> sampledGridX;
	std::vector<double> sampledGridY;
	
	//Matrix where to store bkg and rms info (used to rescue bad tiles)
	int bkgRescueTileSize= 3;
	TMatrixD* sampledBkgLevel= 0;
	TMatrixD* sampledRMS= 0;
		
	//Init bkg data
	try {	
		//Init vectors
		sampledGridX.assign(nTilesX,0);
		sampledGridY.assign(nTilesY,0);

		//Allocate matrix
		sampledBkgLevel= new TMatrixD(nTilesX,nTilesY);
		sampledBkgLevel->Zero();
		sampledRMS= new TMatrixD(nTilesX,nTilesY);
		sampledRMS->Zero();
		
		int counter= 0;
	
		for(int i=0;i<nTilesX;i++){
			int ix= ixList[i];
			int ix_min= ix-TileSizeX/2;
			int ix_max= ix_min+TileSizeX-1;
			double x= img->GetXaxis()->GetBinCenter(ix+1);
			sampledGridX[i]= x;

			for(int j=0;j<nTilesY;j++){
				counter++;
				int iy= iyList[j];
				int iy_min= iy-TileSizeY/2;
				int iy_max= iy_min+TileSizeY-1;
				double y= img->GetYaxis()->GetBinCenter(iy+1);
				sampledGridY[j]= y;

				//## Get tile from image	
				//cout<<"BkgFinder::FindGridBkg(): INFO: Retrieving tile no. "<<counter<<" C("<<ix<<","<<iy<<" Cxy("<<x<<","<<y<<") ix("<<ix_min<<","<<ix_max<<") iy("<<iy_min<<","<<iy_max<<") x("<<x_min<<","<<x_max<<") y("<<y_min<<","<<y_max<<")"<<endl;

				Img* TileImg= img->GetTile(ix_min,ix_max,iy_min,iy_max);
				if(!TileImg) {
					//cerr<<"BkgFinder::ComputeLocalGridBkg(): ERROR: Cannot get tile from image!"<<endl; 
					throw std::runtime_error("Cannot get tile from image!");
				}
			
				//## Compute bkg for this tile
				bool isValidSampling= true;
				BkgSampleData tileBkgData;
				tileBkgData.id= nTiles;
				if(ComputeSampleBkg(tileBkgData,TileImg,estimator)<0){
					//cerr<<"BkgFinder::ComputeLocalGridBkg(): WARN: Background estimation failed for tile no. "<<nTiles<<"!"<<endl;
					WARN_LOG("Background estimation failed for tile no. "<<nTiles<<"!");
					isValidSampling= false;
				}

				//## Delete tile image
				if(TileImg) TileImg->Delete();

				//## Try to assign a valid bkg estimate using neighbors	
				if(!isValidSampling){
					int nGoodPreviousSamples= 0;
					for(int s=1;s<=bkgRescueTileSize;s++){
						if(j-s>=0 && (*sampledBkgLevel)(i,j-s)!=0 && (*sampledRMS)(i,j-s)!=0){
							(*sampledBkgLevel)(i,j)+= (*sampledBkgLevel)(i,j-s);
							(*sampledRMS)(i,j)+= (*sampledRMS)(i,j-s);
							nGoodPreviousSamples++;
						}
						if(i-s>=0 && (*sampledBkgLevel)(i-s,j)!=0 && (*sampledRMS)(i-s,j)!=0){
							(*sampledBkgLevel)(i,j)+= (*sampledBkgLevel)(i-s,j);
							(*sampledRMS)(i,j)+= (*sampledRMS)(i-s,j);
							nGoodPreviousSamples++;
						}
					}
					if(nGoodPreviousSamples>0) {
						(*sampledBkgLevel)(i,j)/= (double)nGoodPreviousSamples;
						(*sampledRMS)(i,j)/= (double)nGoodPreviousSamples;
					}	
					else{
						(*sampledBkgLevel)(i,j)= 0;
						(*sampledRMS)(i,j)= 0;
					}
				}//close if !isValidSampling


				//## Fill bkg samplings
				(bkgData->BkgSamplings).push_back(tileBkgData);
				(*sampledBkgLevel)(i,j)= tileBkgData.bkgLevel;
				(*sampledRMS)(i,j)= tileBkgData.bkgRMS;			
			}//end loop Y
		}//end loop X
	}//close try block
	catch( std::exception &ex ) {
		//cerr << "BkgFinder::ComputeLocalGridBkg(): ERROR: Exception detected: " << ex.what() << endl;
		ERROR_LOG("Exception detected in bkg tile computation (err=" << ex.what()<<")");
		return -1;
  } 
	catch(...) {
		//cerr << "BkgFinder::ComputeLocalGridBkg(): ERROR: C++ exception (unknown reason)" << endl;
		ERROR_LOG("Unknown exception catched in bkg tile computation!"); 
		return -1;
  }	

	
	//## Perform the 2D interpolation
	INFO_LOG("Start bkg 2d interpolation...");
	try {		
		
		// Construct the grid in each dimension (note that we will pass in a sequence of iterators pointing to the beginning of each grid)
		//cout<<"BkgFinder::ComputeLocalGridBkg(): INFO: Build 2D grid for interpolation (nTilesX="<<nTilesX<<", nTilesY="<<nTilesY<<")..."<<endl;
		DEBUG_LOG("Build 2D grid for interpolation (nTilesX="<<nTilesX<<", nTilesY="<<nTilesY<<")...");

		long int num_elements= nTilesX*nTilesY;
		std::vector<double> fbkg_values(num_elements,0);
		std::vector<double> frms_values(num_elements,0);
		std::vector<double> interp_gridX = CodeUtils::linspace(xlim[0],xlim[1], Nx);
		std::vector<double> interp_gridY = CodeUtils::linspace(ylim[0],ylim[1], Ny);
		std::vector<double> interpolatedBkg;
  	std::vector<double> interpolatedRMS;
  	for (int i=0; i<nTilesX; i++) {
    	for (int j=0; j<nTilesY; j++) {
				long int gBin= i*nTilesY + j;
				fbkg_values[gBin] = (*sampledBkgLevel)(i,j);
				frms_values[gBin] = (*sampledRMS)(i,j);
			}
  	}

		//cout<<"BkgFinder::ComputeLocalGridBkg(): INFO: Interpolating bkg ...."<<endl;
		DEBUG_LOG("Interpolating bkg map ...");
		int status= MathUtils::BiLinearInterpolation(sampledGridX,sampledGridY,fbkg_values,interp_gridX,interp_gridY,interpolatedBkg);
		if(status<0){
			//cerr<<"BkgFinder::ComputeLocalGridBkg(): INFO: Failed to interpolate bkg ...."<<endl;
			throw std::runtime_error("Failed to interpolate bkg!");
		}

		//cout<<"BkgFinder::ComputeLocalGridBkg(): INFO: Interpolating noise ...."<<endl;
		DEBUG_LOG("Interpolating noise map ...");
		status= MathUtils::BiLinearInterpolation(sampledGridX,sampledGridY,frms_values,interp_gridX,interp_gridY,interpolatedRMS);
		if(status<0){
			//cerr<<"BkgFinder::ComputeLocalGridBkg(): INFO: Failed to interpolate noise ...."<<endl;
			throw std::runtime_error("Failed to interpolate noise!");
		}
		
		//## Fill images
		//cout<<"BkgFinder::ComputeLocalGridBkg(): INFO: Init bkg image..."<<endl;
		DEBUG_LOG("Init bkg image...");
		TString bkgImgName= Form("%s_bkg",img->GetName());
		bkgData->BkgMap= img->GetCloned(std::string(bkgImgName),true,true);
		(bkgData->BkgMap)->Reset();
		
		//cout<<"BkgFinder::ComputeLocalGridBkg(): INFO: Init noise image..."<<endl;
		DEBUG_LOG("Init noise image...");
		TString noiseImgName= Form("%s_noise",img->GetName());
		bkgData->NoiseMap= img->GetCloned(std::string(noiseImgName),true,true);
		(bkgData->NoiseMap)->Reset();
		
		//cout<<"BkgFinder::ComputeLocalGridBkg(): INFO: Filling bkg/noise images..."<<endl;
		DEBUG_LOG("Filling bkg/noise images...");
		for (unsigned int i=0; i<interp_gridX.size(); i++) {
    	for (unsigned int j=0; j<interp_gridY.size(); j++) {
				long int gBin= i*interp_gridY.size() + j;
				double thisBkgValue= interpolatedBkg[gBin];
	  		double thisRMSValue= interpolatedRMS[gBin];
				if(thisBkgValue==0 || thisRMSValue<=0 ){
					//cout<<"BkgFinder::ComputeLocalGridBkg(): WARN: Interpolated value is zero (bkg="<<thisBkgValue<<", rms="<<thisRMSValue<<")"<<endl;
					WARN_LOG("Interpolated value is zero (bkg="<<thisBkgValue<<", rms="<<thisRMSValue<<")");
				}

				double binX= img->GetXaxis()->GetBinCenter(i+1);
				double binY= img->GetYaxis()->GetBinCenter(j+1);				
				(bkgData->BkgMap)->FillPixel(binX,binY,thisBkgValue);
				(bkgData->NoiseMap)->FillPixel(binX,binY,thisRMSValue);
			}//end loop bins Y
  	}//end loop bins X

		//Delete stuff
		//cout<<"BkgFinder::ComputeLocalGridBkg(): INFO: Deleting sample bkg/noise..."<<endl;
		DEBUG_LOG("Deleting sample bkg/noise...");
		if(sampledBkgLevel) sampledBkgLevel->Delete();
		if(sampledRMS) sampledRMS->Delete();
		
		
	}//close try block
	catch( std::exception &ex ) {
		//cerr << "BkgFinder::ComputeLocalGridBkg(): ERROR: Exception detected in interpolation: " << ex.what() << endl;		
		ERROR_LOG("Exception detected in interpolation (err="<<ex.what()<<")");
		if(sampledBkgLevel) sampledBkgLevel->Delete();
		if(sampledRMS) sampledRMS->Delete();
		return -1;
  } 
	catch(...) { 
		//cerr << "BkgFinder::ComputeLocalGridBkg(): ERROR: C++ exception (unknown reason) in interpolation!" << endl;
		ERROR_LOG("Unknown exception caught in interpolation!");
		if(sampledBkgLevel) sampledBkgLevel->Delete();
		if(sampledRMS) sampledRMS->Delete();
		return -1;
  }		

	//cout<<"BkgFinder::ComputeLocalGridBkg(): INFO: End local bkg computation"<<endl;
	DEBUG_LOG("End local bkg computation");

	return 0;

}//close BkgFinder::ComputeLocalGridBkg()


int BkgFinder::ComputeGlobalBkg(BkgData* bkgData,Img* img,int estimator){

	//## Compute bkg
	BkgSampleData bkgSampleData;
	if(ComputeSampleBkg(bkgSampleData,img,estimator)<0){
		//cerr<<"BkgFinder::ComputeGlobalBkg(): ERROR: Failed to compute bkg estimator!"<<endl;
		ERROR_LOG("Failed to compute bkg estimator!");
		return -1;
	}			
		
	bkgData->gBkg= bkgSampleData.bkgLevel;
	bkgData->gNoise= bkgSampleData.bkgRMS;

	return 0;

}//close ComputeGlobalBkg()



int BkgFinder::ComputeSampleBkg(BkgSampleData& bkgSampleData,Img* img,int estimator){

	//## Check input image
	if(!img) {
		//cerr<<"BkgFinder::ComputeSampleBkg(): ERROR: Null ptr to given image!"<<endl;
		ERROR_LOG("Null ptr to given image!");
		return -1;
	}
		
	//## Compute and get stats for this tile
	if(!img->HasStats()){//computing stats
		if( img->ComputeStats(true,false,false)<0 ){	
			//cerr<<"BkgFinder::ComputeSampleBkg(): ERROR: Failed to compute stats!"<<endl;
			ERROR_LOG("Failed to compute stats!");	
			return -1;	
		}
	}//close if
	Caesar::ImgStats* stats= img->GetPixelStats();
	
	//## Fill bkg sampling data
	bkgSampleData.npix= stats->n;
	bkgSampleData.ix_min= img->GetXaxis()->GetBinCenter(1);
	bkgSampleData.ix_max= img->GetXaxis()->GetBinCenter(img->GetNbinsX());
	bkgSampleData.iy_min= img->GetYaxis()->GetBinCenter(1);
	bkgSampleData.iy_max= img->GetYaxis()->GetBinCenter(img->GetNbinsY());

	// == mean bkg ==
	if(estimator == BkgEstimator::eMeanBkg){
		bkgSampleData.bkgLevel= stats->mean;
		bkgSampleData.bkgRMS= stats->rms;
	}
		
	//== median bkg ==
	else if(estimator == BkgEstimator::eMedianBkg){
		bkgSampleData.bkgLevel= stats->median;
		bkgSampleData.bkgRMS= stats->medianRMS;
	}

	//== BIWEIGHT bkg ==
	else if(estimator == BkgEstimator::eBiWeightBkg){
		bkgSampleData.bkgLevel= stats->bwLocation;
		bkgSampleData.bkgRMS= stats->bwScale;
	}
	//== CLIPPED MEDIAN bkg ==
	else if(estimator == BkgEstimator::eMedianClippedBkg){
		bkgSampleData.bkgLevel= stats->clippedMedian;
		bkgSampleData.bkgRMS= stats->clippedRMS;
	}

	//invalid estimator
	else{
		//cerr<<"BkgFinder::ComputeSampleBkg(): ERROR: Invalid bkg estimator selected !"<<endl;
		ERROR_LOG("Invalid bkg estimator selected!");
		return -1;
	}

	//## Integrity checks
	if(bkgSampleData.npix<10){
		//cout<<"BkgFinder::ComputeSampleBkg(): WARN: Too few pixels available (n="<<bkgSampleData.npix<<") for bkg estimation, unreliable sampled bkg!"<<endl;
		WARN_LOG("Too few pixels available (n="<<bkgSampleData.npix<<") for bkg estimation, unreliable sampled bkg!");	
		bkgSampleData.isReliable= false;
	}

	if( bkgSampleData.bkgLevel==0 || TMath::IsNaN(bkgSampleData.bkgLevel) || fabs(bkgSampleData.bkgLevel)==TMath::Infinity() ){
		//cout<<"BkgFinder::ComputeSampleBkg(): WARN: Invalid bkg level estimated (bkgLevel="<<bkgSampleData.bkgLevel<<"), unreliable sampled bkg!"<<endl;
		WARN_LOG("Invalid bkg level estimated (bkgLevel="<<bkgSampleData.bkgLevel<<"), unreliable sampled bkg!");
		bkgSampleData.isReliable= false;
	}
	if( bkgSampleData.bkgRMS==0 || TMath::IsNaN(bkgSampleData.bkgRMS) || fabs(bkgSampleData.bkgRMS)==TMath::Infinity() ){
		//cout<<"BkgFinder::ComputeSampleBkg(): WARN: Invalid noise level estimated (bkgRMS="<<bkgSampleData.bkgRMS<<"), unreliable sampled bkg!"<<endl;
		WARN_LOG("Invalid noise level estimated (bkgRMS="<<bkgSampleData.bkgRMS<<"), unreliable sampled bkg!");
		bkgSampleData.isReliable= false;
	}

	return 0;

}//close BkgFinder::ComputeSampleBkg()

}//close namespace

