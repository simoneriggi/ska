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
#include <Image.h>
#include <CodeUtils.h>
#include <MathUtils.h>
#include <Logger.h>
#include <Consts.h>

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


BkgFinder::BkgFinder() 
{

}//close costructor


BkgFinder::~BkgFinder() 
{

}//close destructor


//======================================================
//==            NEW IMAGE METHODS
//======================================================
ImgBkgData* BkgFinder::FindBkg(Image* img,int estimator,bool computeLocalBkg,int boxSizeX,int boxSizeY, double gridStepSizeX,double gridStepSizeY,bool use2ndPass,bool skipOutliers,double seedThr,double mergeThr,int minPixels)
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
	ImgBkgData* bkgData= new ImgBkgData;

	
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
		INFO_LOG("Local bkg computation completed!");
	}//close if computeLocalBkg

	
	//## Skip outliers?
	//## Search and exclude significant blobs (both positive & negative excesses) 
	//## using the first estimate bkg/noise
	if(skipOutliers){
		INFO_LOG("Improving bkg estimate by skipping outliers ...");

		//Get significance map
		INFO_LOG("Computing the significance map ...");
		Image* significanceMap= img->GetSignificanceMap(bkgData,computeLocalBkg);
		if(!significanceMap){
			ERROR_LOG("Failed to compute the significance map (needed to exclude blobs)!");
			delete bkgData;
			bkgData= 0;
			return 0;
		}	

		//Find blobs
		INFO_LOG("Finding compact blobs to be tagged as outliers...");
		std::vector<Source*> blobs;
		bool findNegativeExcess= true;
		bool mergeBelowSeed= false;
		bool findNestedSources= false;
		int status= img->FindCompactSource(blobs,significanceMap,bkgData,seedThr,mergeThr,minPixels,findNegativeExcess,mergeBelowSeed,findNestedSources);
		if(status<0){
			ERROR_LOG("Failed to find significant blobs!");
			delete bkgData;
			bkgData= 0;
			delete significanceMap;
			significanceMap= 0;
			return 0;
		}

		//Find image without outliers (set to zero)
		INFO_LOG("Computing image without outliers (set to zero)...");
		Image* img_wOutliers= img->GetSourceMask(blobs,false,true);//invert mask
		if(!img_wOutliers){
			ERROR_LOG("Failed to compute image with blob outliers subtracted!");
			delete bkgData;
			bkgData= 0;
			delete significanceMap;
			significanceMap= 0;
			for(unsigned int i=0;i<blobs.size();i++){
				if(blobs[i]){
					delete blobs[i];
					blobs[i]= 0;
				}
			}
			blobs.clear();
			return 0;
		}//close if

		INFO_LOG("img without outliers info: min/max="<<img_wOutliers->GetMinimum()<<"/"<<img_wOutliers->GetMaximum());

		//Recompute bkg on residual map (using this function recursively)
		//Do not skip outliers this time!
		INFO_LOG("Recomputing bkg on residual map...");
		ImgBkgData* robustBkgData= FindBkg(img_wOutliers,estimator,computeLocalBkg,boxSizeX,boxSizeY,gridStepSizeX,gridStepSizeY,use2ndPass,false);
		if(!robustBkgData){
			ERROR_LOG("Failed to compute bkg over image with blob outliers subtracted!");
			delete bkgData;
			bkgData= 0;
			delete significanceMap;
			significanceMap= 0;
			for(unsigned int i=0;i<blobs.size();i++){
				if(blobs[i]){
					delete blobs[i];
					blobs[i]= 0;
				}
			}
			blobs.clear();
			return 0;
		}
		INFO_LOG("Robust bkg map min/max="<<(robustBkgData->BkgMap)->GetMinimum()<<"/"<<(robustBkgData->BkgMap)->GetMaximum()<<", rms map min/max="<<(robustBkgData->NoiseMap)->GetMinimum()<<"/"<<(robustBkgData->NoiseMap)->GetMaximum());


		//Override main bkgData with robust estimates
		for(unsigned int i=0;i<(robustBkgData->BkgSamplings).size();i++) {
			(bkgData->BkgSamplings)[i]= (robustBkgData->BkgSamplings)[i];
		}
		bkgData->CopyBkgMap(robustBkgData->BkgMap);//copy new bkg map (delete previous)
		bkgData->CopyNoiseMap(robustBkgData->NoiseMap);//copy new noise map (delete previous)
		bkgData->gBkg= robustBkgData->gBkg;
		bkgData->gNoise= robustBkgData->gNoise;
				
		//Delete stuff
		delete significanceMap;
		significanceMap= 0;
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


int BkgFinder::FindLocalGridBkg(ImgBkgData* bkgData,Image* img,int estimator,long int boxSizeX,long int boxSizeY, double gridStepSizeX, double gridStepSizeY,bool use2ndPass){

	//## Compute bkg data
	INFO_LOG("Computing local bkg (1st pass)...");
	if(ComputeLocalGridBkg(bkgData,img, estimator, boxSizeX, boxSizeY, gridStepSizeX,gridStepSizeY)<0){	
		ERROR_LOG("Computation of local background failed for this image!");
		return -1;
	}

	
	//## Improve rms by recomputing stuff from residual map 
	if(use2ndPass){
		INFO_LOG("Improving rms estimation with a 2nd pass...");

		TString residualMapName= Form("%s_residual",img->GetName().c_str());
		Image* residualMap= img->GetCloned(std::string(residualMapName),true,true);
		residualMap->Add(bkgData->BkgMap,-1);//subtract the bkg level model

		//Compute bkg for residual map
		ImgBkgData* bkgData_residual= new ImgBkgData;
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



int BkgFinder::ComputeLocalGridBkg(ImgBkgData* bkgData,Image* img,int estimator,long int boxSizeX,long int boxSizeY,double gridStepSizeX, double gridStepSizeY){

	//## Check input image
	if(!img) {
		ERROR_LOG("Null ptr to given image!");
		return -1;
	}

	//## Check options
	long int Nx= img->GetNx();
	long int Ny= img->GetNy();
	if(boxSizeX<=0 || boxSizeX>=Nx || boxSizeY<=0 || boxSizeY>=Ny) {
		ERROR_LOG("Invalid box size ("<<boxSizeX<<","<<boxSizeY<<") given (too small or larger than image size)!");
		return -1;
	}
	if(gridStepSizeX<=0 || gridStepSizeY<=0 ){
		ERROR_LOG("Invalid grid step size ("<<gridStepSizeX<<","<<gridStepSizeY<<") given (null or negative)!");
		return -1;
	}
	
	long int TileSizeX= boxSizeX;
	long int TileSizeY= boxSizeY;
	double xlim[2]= {img->GetXmin(),img->GetXmax()};
	double ylim[2]= {img->GetYmin(),img->GetYmax()};
	
	DEBUG_LOG("N("<<Nx<<","<<Ny<<") TileSize=("<<TileSizeX<<","<<TileSizeY<<") GridStepSize("<<gridStepSizeX<<","<<gridStepSizeY<<")");
	
	//## Initialize and count number of rows & cols
	long int indexX_start= TileSizeX/2;
	long int indexY_start= TileSizeY/2;
	long int indexX_end= Nx-1-TileSizeX/2;
	long int indexY_end= Ny-1-TileSizeY/2;
	long int indexX= indexX_start;
	long int indexY= indexY_start;
	long int nTiles= 0;
	long int nTilesX= 0;
	long int nTilesY= 0;
	std::vector<long int> ixList;
	std::vector<long int> iyList;
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
	sampledGridX.assign(nTilesX,0);
	sampledGridY.assign(nTilesY,0);

	long int num_elements= nTilesX*nTilesY;
	std::vector<double> fbkg_values(num_elements,0);
	std::vector<double> frms_values(num_elements,0);

	BkgSampleData initBkgData;
	(bkgData->BkgSamplings).clear();
	(bkgData->BkgSamplings).assign(num_elements,initBkgData);
	
	//Matrix where to store bkg and rms info (used to rescue bad tiles)
	int bkgRescueTileSize= 3;
	
	//Loop over tiles to sample bkg & noise locally
	#ifdef OPENMP_ENABLED
	#pragma omp parallel for collapse(2)
	#endif
	for(long int i=0;i<nTilesX;i++){
		for(long int j=0;j<nTilesY;j++){

			long int index= i*nTilesY + j;

			//## Set grid X info (moved in nested loop because of collapse openmp statements)
			long int ix= ixList[i];
			long int ix_min= ix-TileSizeX/2;
			long int ix_max= ix_min+TileSizeX-1;
			double x= img->GetX(ix);
			sampledGridX[i]= x;

			//## Set grid y info
			long int iy= iyList[j];
			long int iy_min= iy-TileSizeY/2;
			long int iy_max= iy_min+TileSizeY-1;
			double y= img->GetY(iy);
			sampledGridY[j]= y;

			//## Compute bkg for this tile
			bool isValidSampling= true;
			BkgSampleData tileBkgData;
			tileBkgData.id= index;
			if(ComputeSampleBkg(tileBkgData,img,estimator,ix_min,ix_max,iy_min,iy_max)<0){
				WARN_LOG("Background estimation failed for tile (i,j)=("<<i<<","<<j<<")!");
				isValidSampling= false;
			}

			//## Store sampled bkg info
			(bkgData->BkgSamplings)[index]= tileBkgData;

			//## If sampling failed try to assign a valid bkg estimate using neighbors	
			if(isValidSampling){
				fbkg_values[index]= tileBkgData.bkgLevel;
				frms_values[index]= tileBkgData.bkgRMS;		
			}
			else{//failed sampling

				#ifdef OPENMP_ENABLED
				#pragma omp critical 
				#endif
				{
					int nGoodPreviousSamples= 0;
					double recovered_bkg= 0;
					double recovered_rms= 0;

					for(int s=1;s<=bkgRescueTileSize;s++){
						int index_y= i*nTilesY + (j-s);
						int index_x= (i-s)*nTilesY + j;
						double bkg_y= fbkg_values[index_y];
						double rms_y= frms_values[index_y];
						double bkg_x= fbkg_values[index_x];
						double rms_x= frms_values[index_x];
					
						bool isGoodSampling_y= (j-s>=0 && bkg_y!=0 && rms_y!=0);
						bool isGoodSampling_x= (i-s>=0 && bkg_x!=0 && rms_x!=0);

						if(isGoodSampling_y){
							recovered_bkg+= bkg_y;
							recovered_rms+= rms_y;
							nGoodPreviousSamples++;
						}
						if(isGoodSampling_x){
							recovered_bkg+= bkg_x;
							recovered_rms+= rms_x;
							nGoodPreviousSamples++;
						}	
					}//end loop previous tiles

					if(nGoodPreviousSamples>0) {
						recovered_bkg/= (double)nGoodPreviousSamples;
						recovered_rms/= (double)nGoodPreviousSamples;
						fbkg_values[index]= recovered_bkg;
						frms_values[index]= recovered_rms;		
					}	
					else{//giving up!
						fbkg_values[index]= 0;
						frms_values[index]= 0;
					}
				}//close else isValidSampling
			}//close critical region

		}//end loop Y
	}//end loop X
	
	
	//## Perform the 2D interpolation
	INFO_LOG("Start bkg 2d interpolation...");
	std::vector<double> interp_gridX;
	std::vector<double> interp_gridY;
	
	//Construct the interpolation grid
	try{
		interp_gridX = CodeUtils::linspace(xlim[0],xlim[1], Nx);
		interp_gridY = CodeUtils::linspace(ylim[0],ylim[1], Ny);
	}
	catch( std::exception& e ) {
		ERROR_LOG("C++ exception while creating the interpolation grid (err="<<e.what()<<")");
		return -1;
  } 
		

	//If OPENMP if defined split the bkg noise interpolation in two separate concurrent thread	
	std::vector<double> interpolatedBkg;
	std::vector<double> interpolatedRMS;
	bool splitWork= false;
	int errflag= 0;
	#ifdef OPENMP_ENABLED
		splitWork= true;
	#endif

	#ifdef OPENMP_ENABLED
	#pragma omp parallel num_threads(2) reduction(+: errflag)
	#endif
	{
		int thread_id= omp_get_thread_num();

		//Perform the bkg map interpolation	
		if( !splitWork || (splitWork && thread_id==0) ){
			DEBUG_LOG("Interpolating bkg map ...");
			
			try{
				int status= MathUtils::BiLinearInterpolation(sampledGridX,sampledGridY,fbkg_values,interp_gridX,interp_gridY,interpolatedBkg);
				if(status<0){
					ERROR_LOG("Failed to interpolate the bkg map!");
					errflag++;
				}
			}
			catch( std::exception& e ) {
				ERROR_LOG("C++ exception while interpolating the bkg map (err="<<e.what()<<")");
				errflag++;
  		}
			catch(...) { 
				ERROR_LOG("Unknown exception caught while interpolating the bkg map!");
				errflag++;
  		}
		}//close if thread 0
			
		//Perform the noise map interpolation	
		if( !splitWork || (splitWork && thread_id==1) ){
						
			DEBUG_LOG("Interpolating noise map ...");
			try{
				int status= MathUtils::BiLinearInterpolation(sampledGridX,sampledGridY,frms_values,interp_gridX,interp_gridY,interpolatedRMS);
				if(status<0){
					ERROR_LOG("Failed to interpolate noise!");
					errflag++;
				}
			}
			catch( std::exception& e ) {
				ERROR_LOG("C++ exception while interpolating the noise map (err="<<e.what()<<")");
				errflag++;
  		}
			catch(...) { 
				ERROR_LOG("Unknown exception caught while interpolating the noise map!");
				errflag++;
  		}
		}//close if thread 1

	}//close parallel section
	
	//Check errors in interpolation
	if(errflag>0){
		ERROR_LOG("One or more failures occurred in bkg/noise map interpolation...");
		return -1;
	}
	
	
	//## Fill images
	DEBUG_LOG("Init bkg image...");
	TString bkgImgName= Form("%s_bkg",img->GetName().c_str());
	(bkgData->BkgMap)= img->GetCloned(std::string(bkgImgName),true,true);
	(bkgData->BkgMap)->Reset();
		
	DEBUG_LOG("Init noise image...");
	TString noiseImgName= Form("%s_noise",img->GetName().c_str());
	(bkgData->NoiseMap)= img->GetCloned(std::string(noiseImgName),true,true);
	(bkgData->NoiseMap)->Reset();
		
	INFO_LOG("Filling bkg/noise images after interpolation ...");
	#ifdef OPENMP_ENABLED
		Caesar::StatMoments<double> bkg_moments_t;		
		Caesar::StatMoments<double> rms_moments_t;	
		std::vector<Caesar::StatMoments<double>> bkg_parallel_moments;
		std::vector<Caesar::StatMoments<double>> rms_parallel_moments;
	
		//#pragma omp declare reduction (merge : std::vector<Caesar::StatMoments<double>> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
		#pragma omp parallel private(bkg_moments_t,rms_moments_t) //reduction(merge: bkg_parallel_moments,rms_parallel_moments)
		{

			int thread_id= omp_get_thread_num();
			int nthreads= SysUtils::GetOMPThreads();

			#pragma omp single
   		{
     		bkg_parallel_moments.assign(nthreads,Caesar::StatMoments<double>());
				rms_parallel_moments.assign(nthreads,Caesar::StatMoments<double>());
   		}

			#pragma omp for collapse(2) 
			for (size_t i=0; i<interp_gridX.size(); i++) {
   			for (size_t j=0; j<interp_gridY.size(); j++) {
					long int gBin= i*interp_gridY.size() + j;
					double thisBkgValue= interpolatedBkg[gBin];
	  			double thisRMSValue= interpolatedRMS[gBin];
					if(thisBkgValue==0 || thisRMSValue<=0 ){
						WARN_LOG("Interpolated value is zero (bkg="<<thisBkgValue<<", rms="<<thisRMSValue<<")");
					}
	
					(bkgData->BkgMap)->FillPixelMT(bkg_moments_t,i,j,thisBkgValue);
					(bkgData->NoiseMap)->FillPixelMT(rms_moments_t,i,j,thisRMSValue);
				}//end loop bins Y
  		}//end loop bins X
		
			//Fill parallel moments per thread
			//bkg_parallel_moments.push_back(bkg_moments_t);
			//rms_parallel_moments.push_back(rms_moments_t);
			bkg_parallel_moments[thread_id]= bkg_moments_t;
			rms_parallel_moments[thread_id]= rms_moments_t;

			DEBUG_LOG("Thread id="<<thread_id<<" bkg min/max="<<bkg_parallel_moments[thread_id].minVal<<"/"<<bkg_parallel_moments[thread_id].maxVal);
			DEBUG_LOG("Thread id="<<thread_id<<" rms min/max="<<rms_parallel_moments[thread_id].minVal<<"/"<<rms_parallel_moments[thread_id].maxVal);
			
		}//close parallel section

		//Update moments from parallel estimates
		//- Bkg map moments
		Caesar::StatMoments<double> bkg_moments;
		if(Caesar::StatsUtils::ComputeMomentsFromParallel(bkg_moments,bkg_parallel_moments)<0){
			ERROR_LOG("Failed to compute cumulative bkg map moments from parallel estimates (NB: image will have wrong moments!)");
			return -1;
		}
		//DEBUG_LOG("bkg moments aggregated");
		//bkg_moments.Print();
		(bkgData->BkgMap)->SetMoments(bkg_moments);
		
		//DEBUG_LOG("bkg moments aggregated (after set)");
		//((bkgData->BkgMap)->GetMoments()).Print();

		//- RMS map moments
		Caesar::StatMoments<double> rms_moments;
		if(Caesar::StatsUtils::ComputeMomentsFromParallel(rms_moments,rms_parallel_moments)<0){
			ERROR_LOG("Failed to compute cumulative rms map moments from parallel estimates (NB: image will have wrong moments!)");
			return -1;
		}
		//DEBUG_LOG("rms moments aggregated");
		//rms_moments.Print();
		(bkgData->NoiseMap)->SetMoments(rms_moments);

		//DEBUG_LOG("rms moments aggregated (after set)");
		//((bkgData->NoiseMap)->GetMoments()).Print();

	#else
		for (size_t i=0; i<interp_gridX.size(); i++) {
   		for (size_t j=0; j<interp_gridY.size(); j++) {
				long int gBin= i*interp_gridY.size() + j;
				double thisBkgValue= interpolatedBkg[gBin];
	  		double thisRMSValue= interpolatedRMS[gBin];
				if(thisBkgValue==0 || thisRMSValue<=0 ){
					WARN_LOG("Interpolated value is zero (bkg="<<thisBkgValue<<", rms="<<thisRMSValue<<")");
				}
	
				(bkgData->BkgMap)->FillPixel(i,j,thisBkgValue);
				(bkgData->NoiseMap)->FillPixel(i,j,thisRMSValue);
			}//end loop bins Y
  	}//end loop bins X
	#endif
	
	DEBUG_LOG("End local bkg computation");

	return 0;

}//close ComputeLocalGridBkg()


int BkgFinder::ComputeGlobalBkg(ImgBkgData* bkgData,Image* img,int estimator){

	//## Compute bkg
	BkgSampleData bkgSampleData;
	if(ComputeSampleBkg(bkgSampleData,img,estimator)<0){
		ERROR_LOG("Failed to compute bkg estimator!");
		return -1;
	}			
		
	bkgData->gBkg= bkgSampleData.bkgLevel;
	bkgData->gNoise= bkgSampleData.bkgRMS;

	return 0;

}//close ComputeGlobalBkg()



int BkgFinder::ComputeSampleBkg(BkgSampleData& bkgSampleData,Image* img,int estimator,long int ix_min,long int ix_max,long int iy_min,long int iy_max){

	//## Check input image
	if(!img) {
		ERROR_LOG("Null ptr to given image!");
		return -1;
	}
	
	//## Check if read tile
	bool isTile= (ix_min!=-1 && ix_max!=-1 && iy_min!=-1 && iy_max!=-1);
	if(!isTile){//full image
		ix_min= 0;
		ix_max= img->GetNx()-1;
		iy_min= 0;
		iy_max= img->GetNy()-1;	
	}

	//## Fill bkg sampling data
	bkgSampleData.ix_min= ix_min;
	bkgSampleData.ix_max= ix_max;
	bkgSampleData.iy_min= iy_min;
	bkgSampleData.iy_max= iy_max;

	bool skipNegativePixels= false;
	int status= 0;
	long int npix= 0;

	// == mean bkg ==
	if(estimator == BkgEstimator::eMeanBkg){
		float mean, rms;
		status= img->GetTileMeanStats(mean,rms,npix,ix_min,ix_max,iy_min,iy_max,skipNegativePixels);
		if(status==0){
			bkgSampleData.bkgLevel= mean;
			bkgSampleData.bkgRMS= rms;
			bkgSampleData.npix= npix;
		}
	}
		
	//== median bkg ==
	else if(estimator == BkgEstimator::eMedianBkg){
		float median, medianRMS;
		status= img->GetTileMedianStats(median,medianRMS,npix,ix_min,ix_max,iy_min,iy_max,skipNegativePixels);
		if(status==0){
			bkgSampleData.bkgLevel= median;
			bkgSampleData.bkgRMS= medianRMS;
			bkgSampleData.npix= npix;
		}
	}

	//== BIWEIGHT bkg ==
	else if(estimator == BkgEstimator::eBiWeightBkg){
		float bwLocation, bwScale;
		double C= 6;
		status= img->GetTileBiWeightStats(bwLocation,bwScale,npix,C,ix_min,ix_max,iy_min,iy_max,skipNegativePixels);	
		if(status==0){
			bkgSampleData.bkgLevel= bwLocation;
			bkgSampleData.bkgRMS= bwScale;
			bkgSampleData.npix= npix;
		}
	}
	//== CLIPPED MEDIAN bkg ==
	else if(estimator == BkgEstimator::eMedianClippedBkg){
		ClippedStats<float> clipped_stats;
		double clipSigma= 3;
		status= img->GetTileClippedStats(clipped_stats,npix,clipSigma,ix_min,ix_max,iy_min,iy_max,skipNegativePixels);
		if(status==0){
			bkgSampleData.bkgLevel= clipped_stats.median;
			bkgSampleData.bkgRMS= clipped_stats.stddev;
			bkgSampleData.npix= npix;
		}
	}

	//invalid estimator
	else{
		ERROR_LOG("Invalid bkg estimator selected!");
		return -1;
	}
	
	//Check if stats computing failed
	if(status<0){
		ERROR_LOG("Failed to compute stats for tile [("<<ix_min<<","<<ix_max<<"), ("<<iy_min<<","<<iy_max<<")]");
		return -1;
	}

	//## Integrity checks
	
	if(bkgSampleData.npix<10){
		WARN_LOG("Too few pixels available (n="<<bkgSampleData.npix<<") for bkg estimation, unreliable sampled bkg!");	
		bkgSampleData.isReliable= false;
	}

	if( bkgSampleData.bkgLevel==0 || TMath::IsNaN(bkgSampleData.bkgLevel) || fabs(bkgSampleData.bkgLevel)==TMath::Infinity() ){
		WARN_LOG("Invalid bkg level estimated (bkgLevel="<<bkgSampleData.bkgLevel<<"), unreliable sampled bkg!");
		bkgSampleData.isReliable= false;
	}
	if( bkgSampleData.bkgRMS==0 || TMath::IsNaN(bkgSampleData.bkgRMS) || fabs(bkgSampleData.bkgRMS)==TMath::Infinity() ){
		WARN_LOG("Invalid noise level estimated (bkgRMS="<<bkgSampleData.bkgRMS<<"), unreliable sampled bkg!");
		bkgSampleData.isReliable= false;
	}

	return 0;

}//close BkgFinder::ComputeSampleBkg()










//======================================================
//==            OLD IMAGE METHODS
//======================================================

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
	long int Nx= img->GetNbinsX();
	long int Ny= img->GetNbinsY();
	if(boxSizeX<=0 || boxSizeX>=Nx || boxSizeY<=0 || boxSizeY>=Ny) {
		ERROR_LOG("Invalid box size ("<<boxSizeX<<","<<boxSizeY<<") given (too small or larger than image size)!");
		return -1;
	}
	if(gridStepSizeX<=0 || gridStepSizeY<=0 ){
		ERROR_LOG("Invalid grid step size ("<<gridStepSizeX<<","<<gridStepSizeY<<") given (null or negative)!");
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
		ERROR_LOG("Null ptr to given image!");
		return -1;
	}
		
	//## Compute and get stats for this tile
	if(!img->HasStats()){//computing stats
		if( img->ComputeStats(true,false,false)<0 ){	
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

