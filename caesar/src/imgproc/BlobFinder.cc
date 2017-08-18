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
* @file BlobFinder.cc
* @class BlobFinder
* @brief Blob finder class
*
* Class to perform blob finding 
* @author S. Riggi
* @date 20/01/2015
*/


#include <BlobFinder.h>
#include <Img.h>
#include <Image.h>
#include <BkgData.h>
#include <Blob.h>
#include <Source.h>
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
#include <numeric>
#include <time.h>
#include <ctime>
#include <queue>

using namespace std;

ClassImp(Caesar::BlobFinder)

namespace Caesar {

BlobFinder::BlobFinder() {

	
}//close costructor

BlobFinder::~BlobFinder(){
	

}//close destructor

//=========================================
//==  NEW IMAGE METHODS 
//=========================================
template <class T>
int BlobFinder::FindBlobs(Image* inputImg,std::vector<T*>& blobs,Image* floodImg,ImgBkgData* bkgData,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool mergeBelowSeed){

	//## Check input img
	if(!inputImg){
		ERROR_LOG("Null ptr to given input image!");
		return -1;
	}

	//## Check if the flood map is provided otherwise set to the input map
	//## NB: In source search it should be the significance map
	//## NB2: It could be used to fill blobs in the input map, conditional to another map (i.e. a binary mask)
	if(!floodImg) floodImg= inputImg;

	//## Check if bkg data are provided and if local bkg is available
	//## If local bkg is available use it, otherwise use global bkg
	bool hasBkgData= false;
	bool hasLocalBkg= false;
	if(bkgData){
		hasBkgData= true;
		hasLocalBkg= bkgData->HasLocalBkg();
	}

	//## Set flood threshold
	//Example: merge=4, seed=5  [4,5] o [4,+inf]
	// merge=-4 seed=-5         [-5,-4] o [-inf,-4]            [-inf,4]
	double floodMinThr= mergeThr;
	double floodMaxThr= std::numeric_limits<double>::infinity();
	double floodMinThr_inv= -std::numeric_limits<double>::infinity();
	double floodMaxThr_inv= -mergeThr;
	if(mergeBelowSeed) {
		floodMaxThr= seedThr;
		floodMinThr_inv= seedThr;
	}

	DEBUG_LOG("Flood thr("<<floodMinThr<<","<<floodMaxThr<<") Flood inv thr("<<floodMinThr_inv<<","<<floodMaxThr_inv<<")");
	
	//## Find seed pixels (above seed threshold)	
	std::vector<long int> pixelSeeds;	
	std::vector<bool> isNegativeExcessSeed;
	long int Nx= inputImg->GetNx();
	long int Ny= inputImg->GetNy();
	long int Ntot= Nx*Ny;
	
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){	
			int gBin= inputImg->GetBin(i+1,j+1);	
			double Z= floodImg->GetBinContent(i+1,j+1);
			bool isNegative= false;
			if(fabs(Z)>=seedThr) {
				if(Z<0) isNegative= true;
				pixelSeeds.push_back(gBin);	
				isNegativeExcessSeed.push_back(isNegative);
			}
		}//end loop y
	}//end loop x
	
	DEBUG_LOG("#"<<pixelSeeds.size()<<" seeds found ...");

	//## Perform cluster finding starting from detected seeds
	int nBlobs= 0;
	T* aBlob= 0;
	Pixel* aPixel= 0;
	std::vector<bool> isAddedInCluster(Ntot,false);
	
	for(unsigned int k=0;k<pixelSeeds.size();k++){
		long int seedPixelId= pixelSeeds[k];	
		long int binX= inputImg->GetBinX(seedPixelId);
		long int binY= inputImg->GetBinY(seedPixelId);
		
		//Check if this seed bin has been already assigned to a cluster
		if(isAddedInCluster[seedPixelId]) continue;
		
		//Skip negative excess seed if not requested
		if(!findNegativeExcess && isNegativeExcessSeed[k]) continue;
		
		//Compute flooded pixels
		std::vector<long int> clusterPixelIds;
		int status= 0;
		if(isNegativeExcessSeed[k]){
			status= FloodFill(floodImg,clusterPixelIds,seedPixelId,floodMinThr_inv,floodMaxThr_inv);
		}
		else {
			status= FloodFill(floodImg,clusterPixelIds,seedPixelId,floodMinThr,floodMaxThr);
		}
		if(status<0) {
			WARN_LOG("Flood fill failed, skip seed!");
			continue;
		}

		//Append cluster pixels to a blob object
		size_t nClusterPixels= clusterPixelIds.size();
		if(nClusterPixels==0 || (int)nClusterPixels<minPixels) {//skip small blobs
			DEBUG_LOG("Blob pixels found @ (x,y)=("<<binX<<","<<binY<<") (N="<<nClusterPixels<<") below npix threshold (thr="<<minPixels<<"), skip blob!");
			continue;
		}
		DEBUG_LOG("Blob found @ (x,y)=("<<binX<<","<<binY<<") (N="<<nClusterPixels<<")");
		
		nBlobs++;	
		
		DEBUG_LOG("Adding new blob (# "<<nBlobs<<") to list...");
		TString blobName= Form("%s_blobId%d",inputImg->GetName().c_str(),nBlobs);
		aBlob= new T;
		aBlob->SetId(nBlobs);	
		aBlob->SetName(std::string(blobName));
		
		for(size_t l=0;l<nClusterPixels;l++){
			long int clusterPixelId= clusterPixelIds[l];	
			if(isAddedInCluster[clusterPixelId]) continue;
			isAddedInCluster[clusterPixelId]= true;//do not forget to add to list of taken pixels!

			long int clusterPixelIdX= inputImg->GetBinX(clusterPixelId);
			long int clusterPixelIdY= inputImg->GetBinY(clusterPixelId);
			double S= inputImg->GetBinContent(clusterPixelId);			
			double Z= floodImg->GetBinContent(clusterPixelId);

			double x= inputImg->GetX(clusterPixelIdX);
			double y= inputImg->GetY(clusterPixelIdY);
			long int ix= clusterPixelIdX;
			long int iy= clusterPixelIdY;
			
			aPixel= new Pixel;
			aPixel->S= S;
			if(fabs(Z)>=seedThr) aPixel->type= Pixel::eSeed;
			else aPixel->type= Pixel::eNormal;
			aPixel->id= clusterPixelId;
			aPixel->SetPhysCoords(x,y);
			aPixel->SetCoords(ix,iy);
			
			if( inputImg->IsEdgeBin(clusterPixelIdX,clusterPixelIdY) ) {
				aBlob->SetEdgeFlag(true);
			}

			//Set bkg data if available
			if(hasBkgData){
				double bkgLevel= bkgData->gBkg;
				double noiseLevel= bkgData->gNoise;
				if(hasLocalBkg){
					bkgLevel= (bkgData->BkgMap)->GetBinContent(clusterPixelId);
					noiseLevel= (bkgData->NoiseMap)->GetBinContent(clusterPixelId);
				}
				aPixel->SetBkg(bkgLevel,noiseLevel);
			}//close if
	
			aBlob->AddPixel(aPixel);
		}//end loop cluster pixels

		//## Check if blobs has pixels
		if(!aBlob->HasPixels()){//no pixel...delete blob!
			delete aBlob;
			aBlob= 0;
			continue;
		}

		//## Compute stats
		DEBUG_LOG("Computing blob stats...");
		aBlob->ComputeStats();
		
		//## Compute morphology parameters
		DEBUG_LOG("Computing blob morphology params...");
		aBlob->ComputeMorphologyParams();

		//## Add blob to list
		blobs.push_back(aBlob);
		
	}//end loop seeds

	INFO_LOG("#"<<blobs.size()<<" blobs found!");

	return 0;

}//close BlobFinder::FindBlobs()
template int BlobFinder::FindBlobs<Blob>(Image* img,std::vector<Blob*>& blobs,Image*,ImgBkgData*,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool mergeBelowSeed);
template int BlobFinder::FindBlobs<Source>(Image* img,std::vector<Source*>& blobs,Image*,ImgBkgData*,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool mergeBelowSeed);


int BlobFinder::FloodFill(Image* img,std::vector<long int>& clusterPixelIds,long int seedPixelId,double floodMinThr,double floodMaxThr){
	
	//Init
	clusterPixelIds.clear();

	//Check image and given seed id
	if(!img){
		ERROR_LOG("Null ptr to image given!");
		return -1;
	}
	if(!img->HasBin(seedPixelId)){//check if given seed actually exists
		ERROR_LOG("Given seed id is outside image range!");
		return -1;
	}

	//Check given flood range
	double seedSignal= img->GetPixelValue(seedPixelId);
	if(seedSignal<floodMinThr || seedSignal>floodMaxThr){
		WARN_LOG("Given flood threshold range does not contain seed, no blobs detected!");
		return -1;
	}
	
	//Add seed to queue and loop over queue
	std::queue<long int> pixelQueue;
	pixelQueue.push(seedPixelId);
	
	int Ntot= img->GetNPixels();
	std::vector<bool> isAddedInQueue(Ntot,false);	
	std::vector<bool> isAddedInCluster(Ntot,false);

	while(!pixelQueue.empty()){

		//Take first pixel in queue, process it and then remove from the queue
		long int gBinId= pixelQueue.front();
		long int binIdX= img->GetBinX(gBinId);
		long int binIdY= img->GetBinY(gBinId);
		pixelQueue.pop();

		//Loop on row pixels above threshold
		while (img->IsBinContentInRange(binIdX-1,binIdY,floodMinThr,floodMaxThr)){
    	binIdX--;
    }//close while loop
    
		bool spanUp = false;
    bool spanDown = false;
		 
		while (img->IsBinContentInRange(binIdX,binIdY,floodMinThr,floodMaxThr)) {
   		long int gBinId_cluster= img->GetBin(binIdX,binIdY);
			if(!isAddedInCluster[gBinId_cluster]) {
				clusterPixelIds.push_back(gBinId_cluster);
				isAddedInCluster[gBinId_cluster]= true;
			}
			
			//Search up pixels
			long int gBinId_up= img->GetBin(binIdX,binIdY+1);

			if (!spanUp && img->IsBinContentInRange(binIdX,binIdY+1,floodMinThr,floodMaxThr)) {
      	if(!isAddedInQueue[gBinId_up]) {
					pixelQueue.push(gBinId_up);
					isAddedInQueue[gBinId_up]= true;
					spanUp = true;
				} 
			}//close if
			else if (spanUp && !img->IsBinContentInRange(binIdX,binIdY+1,floodMinThr,floodMaxThr)){
				spanUp = false;
      }

			//Search down pixel
			long int gBinId_down= img->GetBin(binIdX,binIdY-1);

			if (!spanDown && img->IsBinContentInRange(binIdX,binIdY-1,floodMinThr,floodMaxThr)) {
     		if(!isAddedInQueue[gBinId_down]) {
					pixelQueue.push(gBinId_down);
					isAddedInQueue[gBinId_down]= true;
					spanDown = true;
				} 
      }//close if 
			else if (spanDown && !img->IsBinContentInRange(binIdX,binIdY-1,floodMinThr,floodMaxThr)) {
				spanDown = false;
      }
      binIdX++;
		}//end while loop
	}//end queue loop
	
	//Append cluster pixels to a source object
	DEBUG_LOG("#"<<clusterPixelIds.size()<<" cluster pixels found around given seed "<<seedPixelId);
	
	return 0;

}//close BlobFinder::FloodFill()

//=========================================
//==  OLD IMAGE METHODS 
//=========================================
int BlobFinder::FloodFill(Img* img,std::vector<int>& clusterPixelIds,int seedPixelId,double floodMinThr,double floodMaxThr){
	
	//Init
	clusterPixelIds.clear();

	//Check image and given seed id
	if(!img){
		ERROR_LOG("Null ptr to image given!");
		return -1;
	}
	if(!img->HasBin(seedPixelId)){//check if given seed actually exists
		ERROR_LOG("Given seed id is outside image range!");
		return -1;
	}

	//Check given flood range
	double seedSignal= img->GetBinContent(seedPixelId);
	if(seedSignal<floodMinThr || seedSignal>floodMaxThr){
		WARN_LOG("Given flood threshold range does not contain seed, no blobs detected!");
		return -1;
	}
	
	//Add seed to queue and loop over queue
	std::queue<int> pixelQueue;
	pixelQueue.push(seedPixelId);
	
	int Ntot= img->GetNcells();//this accounts for underflow/overflow bins (e.g. Ntot= (Nx+2)*(Ny+2))
	std::vector<bool> isAddedInQueue(Ntot,false);	
	std::vector<bool> isAddedInCluster(Ntot,false);

	while(!pixelQueue.empty()){

		//Take first pixel in queue, process it and then remove from the queue
		int gBinId= pixelQueue.front(); 
		int binIdX, binIdY, binIdZ; 
		img->GetBinXYZ(gBinId,binIdX,binIdY,binIdZ);
		pixelQueue.pop();

		//Loop on row pixels above threshold
		while (img->IsBinContentInRange(binIdX-1,binIdY,floodMinThr,floodMaxThr)){
    	binIdX--;
    }//close while loop
    
		bool spanUp = false;
    bool spanDown = false;
		 
		while (img->IsBinContentInRange(binIdX,binIdY,floodMinThr,floodMaxThr)) {
   		int gBinId_cluster= img->GetBin(binIdX,binIdY);
			if(!isAddedInCluster[gBinId_cluster]) {
				clusterPixelIds.push_back(gBinId_cluster);
				isAddedInCluster[gBinId_cluster]= true;
			}
			
			//Search up pixels
			int gBinId_up= img->GetBin(binIdX,binIdY+1);

			if (!spanUp && img->IsBinContentInRange(binIdX,binIdY+1,floodMinThr,floodMaxThr)) {
      	if(!isAddedInQueue[gBinId_up]) {
					pixelQueue.push(gBinId_up);
					isAddedInQueue[gBinId_up]= true;
					spanUp = true;
				} 
			}//close if
			else if (spanUp && !img->IsBinContentInRange(binIdX,binIdY+1,floodMinThr,floodMaxThr)){
				spanUp = false;
      }

			//Search down pixel
			int gBinId_down= img->GetBin(binIdX,binIdY-1);

			if (!spanDown && img->IsBinContentInRange(binIdX,binIdY-1,floodMinThr,floodMaxThr)) {
     		if(!isAddedInQueue[gBinId_down]) {
					pixelQueue.push(gBinId_down);
					isAddedInQueue[gBinId_down]= true;
					spanDown = true;
				} 
      }//close if 
			else if (spanDown && !img->IsBinContentInRange(binIdX,binIdY-1,floodMinThr,floodMaxThr)) {
				spanDown = false;
      }
      binIdX++;
		}//end while loop
	}//end queue loop
	
	//Append cluster pixels to a source object
	//int nClusterPixels= (int)clusterPixelIds.size();
	//cout<<"BlobFinder::FloodFill(): INFO: #"<<nClusterPixels<<" cluster pixels found around given seed "<<seedPixelId<<endl;
	
	return 0;

}//close BlobFinder::FloodFill()

template <class T>
int BlobFinder::FindBlobs(Img* inputImg,std::vector<T*>& blobs,Img* floodImg,BkgData* bkgData,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool mergeBelowSeed){

	//## Check input img
	if(!inputImg){
		ERROR_LOG("Null ptr to given input image!");
		return -1;
	}

	//## Check if the flood map is provided otherwise set to the input map
	//## NB: In source search it should be the significance map
	//## NB2: It could be used to fill blobs in the input map, conditional to another map (i.e. a binary mask)
	if(!floodImg) floodImg= inputImg;

	//## Check if bkg data are provided and if local bkg is available
	//## If local bkg is available use it, otherwise use global bkg
	bool hasBkgData= false;
	bool hasLocalBkg= false;
	if(bkgData){
		hasBkgData= true;
		hasLocalBkg= bkgData->HasLocalBkg();
	}

	//## Set flood threshold
	//Example: merge=4, seed=5  [4,5] o [4,+inf]
	// merge=-4 seed=-5         [-5,-4] o [-inf,-4]            [-inf,4]
	double floodMinThr= mergeThr;
	double floodMaxThr= std::numeric_limits<double>::infinity();
	double floodMinThr_inv= -std::numeric_limits<double>::infinity();
	double floodMaxThr_inv= -mergeThr;
	if(mergeBelowSeed) {
		floodMaxThr= seedThr;
		floodMinThr_inv= seedThr;
	}

	DEBUG_LOG("Flood thr("<<floodMinThr<<","<<floodMaxThr<<") Flood inv thr("<<floodMinThr_inv<<","<<floodMaxThr_inv<<")");
	
	//## Find seed pixels (above seed threshold)	
	std::vector<int> pixelSeeds;	
	std::vector<bool> isNegativeExcessSeed;
	int Nx= inputImg->GetNbinsX();
	int Ny= inputImg->GetNbinsY();
	int Ntot= inputImg->GetNcells();
	
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){	
			int gBin= inputImg->GetBin(i+1,j+1);	
			double Z= floodImg->GetBinContent(i+1,j+1);
			bool isNegative= false;
			if(fabs(Z)>=seedThr) {
				if(Z<0) isNegative= true;
				pixelSeeds.push_back(gBin);	
				isNegativeExcessSeed.push_back(isNegative);
			}
		}//end loop y
	}//end loop x
	
	DEBUG_LOG("#"<<pixelSeeds.size()<<" seeds found ...");

	//## Perform cluster finding starting from detected seeds
	int nBlobs= 0;
	T* aBlob= 0;
	Pixel* aPixel= 0;
	std::vector<bool> isAddedInCluster(Ntot,false);
	
	for(unsigned int k=0;k<pixelSeeds.size();k++){
		int seedPixelId= pixelSeeds[k];
		int binX, binY, binZ;
		inputImg->GetBinXYZ(seedPixelId,binX,binY,binZ);
		
		//Check if this seed bin has been already assigned to a cluster
		if(isAddedInCluster[seedPixelId]) continue;
		
		//Skip negative excess seed if not requested
		if(!findNegativeExcess && isNegativeExcessSeed[k]) continue;
		
		//Compute flooded pixels
		std::vector<int> clusterPixelIds;
		int status= 0;
		if(isNegativeExcessSeed[k]){
			status= FloodFill(floodImg,clusterPixelIds,seedPixelId,floodMinThr_inv,floodMaxThr_inv);
		}
		else {
			status= FloodFill(floodImg,clusterPixelIds,seedPixelId,floodMinThr,floodMaxThr);
		}
		if(status<0) {
			WARN_LOG("Flood fill failed, skip seed!");
			continue;
		}

		//Append cluster pixels to a blob object
		int nClusterPixels= (int)clusterPixelIds.size();
		if(nClusterPixels==0 || nClusterPixels<minPixels) {//skip small blobs
			DEBUG_LOG("Blob pixels found @ (x,y)=("<<binX<<","<<binY<<") (N="<<nClusterPixels<<") below npix threshold (thr="<<minPixels<<"), skip blob!");
			continue;
		}
		DEBUG_LOG("Blob found @ (x,y)=("<<binX<<","<<binY<<") (N="<<nClusterPixels<<")");
		
		nBlobs++;	
		
		DEBUG_LOG("Adding new blob (# "<<nBlobs<<") to list...");
		TString blobName= Form("%s_blobId%d",std::string(inputImg->GetName()).c_str(),nBlobs);
		aBlob= new T;
		aBlob->SetId(nBlobs);	
		aBlob->SetName(std::string(blobName));
		
		for(int l=0;l<nClusterPixels;l++){
			int clusterPixelId= clusterPixelIds[l];	
			if(isAddedInCluster[clusterPixelId]) continue;
			isAddedInCluster[clusterPixelId]= true;//do not forget to add to list of taken pixels!
			int clusterPixelIdX, clusterPixelIdY, clusterPixelIdZ;
			inputImg->GetBinXYZ(clusterPixelId,clusterPixelIdX,clusterPixelIdY,clusterPixelIdZ);

			double S= inputImg->GetBinContent(clusterPixelId);			
			double Z= floodImg->GetBinContent(clusterPixelId);
			double x= inputImg->GetXaxis()->GetBinCenter(clusterPixelIdX);
			double y= inputImg->GetYaxis()->GetBinCenter(clusterPixelIdY);
			int ix= clusterPixelIdX-1;
			int iy= clusterPixelIdY-1;
			
			aPixel= new Pixel;
			aPixel->S= S;
			if(fabs(Z)>=seedThr) aPixel->type= Pixel::eSeed;
			else aPixel->type= Pixel::eNormal;
			aPixel->id= clusterPixelId;
			aPixel->SetPhysCoords(x,y);
			aPixel->SetCoords(ix,iy);
			
			if( clusterPixelIdX<=1 || clusterPixelIdY<=1 || clusterPixelIdX>=Nx || clusterPixelIdY>=Ny) 
				aBlob->SetEdgeFlag(true);

			//Set bkg data if available
			if(hasBkgData){
				double bkgLevel= bkgData->gBkg;
				double noiseLevel= bkgData->gNoise;
				if(hasLocalBkg){
					bkgLevel= (bkgData->BkgMap)->GetBinContent(clusterPixelId);
					noiseLevel= (bkgData->NoiseMap)->GetBinContent(clusterPixelId);
				}
				aPixel->SetBkg(bkgLevel,noiseLevel);
			}//close if
	
			aBlob->AddPixel(aPixel);
		}//end loop cluster pixels

		//## Check if blobs has pixels
		if(!aBlob->HasPixels()){//no pixel...delete blob!
			delete aBlob;
			aBlob= 0;
			continue;
		}

		//## Compute stats
		DEBUG_LOG("Computing blob stats...");
		aBlob->ComputeStats();
		
		//## Compute morphology parameters
		DEBUG_LOG("Computing blob morphology params...");
		aBlob->ComputeMorphologyParams();

		//## Add blob to list
		blobs.push_back(aBlob);
		
	}//end loop seeds

	INFO_LOG("#"<<blobs.size()<<" blobs found!");

	return 0;

}//close BlobFinder::FindBlobs()
template int BlobFinder::FindBlobs<Blob>(Img* img,std::vector<Blob*>& blobs,Img*,BkgData*,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool mergeBelowSeed);
template int BlobFinder::FindBlobs<Source>(Img* img,std::vector<Source*>& blobs,Img*,BkgData*,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool mergeBelowSeed);


Img* BlobFinder::GetMultiScaleBlobMask(Img* img,int kernelFactor,double sigmaMin,double sigmaMax,double sigmaStep,int thrModel,double thrFactor){

	//## Check imge
	if(!img){
		ERROR_LOG("Null ptr to given image!");
		return 0;
	}

	//## Init scales
	int nScales= (sigmaMax-sigmaMin)/sigmaStep + 1;
	
	TString imgName= Form("%s_blobMask",std::string(img->GetName()).c_str());	
	Img* blobMask= img->GetCloned(std::string(imgName),true,true);
	blobMask->Reset();

	int nbins= 100;
	std::vector<Img*> filterMaps;
	std::vector<double> thresholdLevels;
	
	for(int i=0;i<nScales;i++){
		double sigma= sigmaMin + i*sigmaStep;
		int kernelSize= kernelFactor*sigma;	
		if(kernelSize%2==0) kernelSize++;
		INFO_LOG("Computing LoG map @ scale "<<sigma<<" (step="<<sigmaStep<<", kernsize="<<kernelSize<<")");

		//Compute LoG filter
		Img* filterMap= img->GetNormLoGImage(kernelSize,sigma,true);
		filterMap->ComputeStats(true,false,true);
		filterMaps.push_back(filterMap);

		//Compute threshold levels
		ImgStats* imgStats= filterMap->GetPixelStats();	
		double median= imgStats->median;
		double medianRMS= imgStats->medianRMS;
		double medianThr= thrFactor*median;
		double medianRMSThr= thrFactor*medianRMS;
		double otsuThr= filterMap->FindOtsuThreshold(nbins);
		double valleyThr= filterMap->FindValleyThreshold(nbins,true);
		double optimalThr= std::max(std::min(otsuThr,valleyThr),medianThr);
		double thrLevel= medianRMSThr;
		if(thrModel==1) thrLevel= optimalThr;
		else if(thrModel==2) thrLevel= medianRMSThr;
		else thrLevel= medianRMSThr;
		thresholdLevels.push_back(thrLevel);	
	}//end loop reso
	
	//Find blobs across scales
	for(int i=0;i<blobMask->GetNbinsX();i++){
		double x= img->GetXaxis()->GetBinCenter(i+1);
		for(int j=0;j<blobMask->GetNbinsY();j++){
			double y= img->GetYaxis()->GetBinCenter(j+1);
			double binContent= img->GetBinContent(i+1,j+1);
			if(binContent==0) continue;

			double wsum= 0;
			int counter= 0;
			for(unsigned int k=0;k<filterMaps.size();k++){
				double w= filterMaps[k]->GetBinContent(i+1,j+1);
				if(w<thresholdLevels[k]) continue;
				wsum+= w;
				counter++;
			}//end loop scales

			blobMask->FillPixel(x,y,counter);			
		}//end loop y
	}//end loop x

	//Clear 
	for(unsigned int k=0;k<filterMaps.size();k++){
		if(filterMaps[k]) filterMaps[k]->Delete();		
	}
	filterMaps.clear();
	return blobMask;

}//close GetMultiScaleBlobMask()

}//close namespace
