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
* @file SLIC.cc
* @class SLIC
* @brief SLIC generator class
*
* Superpixel generator
* @author S. Riggi
* @date 20/01/2015
*/

#include <SLIC.h>
#include <SLICData.h>
#include <Region.h>
#include <Image.h>
#include <CodeUtils.h>
#include <Logger.h>

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

ClassImp(Caesar::SLIC)

namespace Caesar {


SLIC::SLIC() {

}//close costructor

SLIC::~SLIC(){

}//close destructor



int SLIC::SetSPData(SLICData* slicData,Image* img,bool useLogScaleMapping,Image* edgeImage){

	//## Check SLICData and clear existing data
	if(!slicData){
		ERROR_LOG("Null ptr to given SLIC data!");
		return -1;
	}
	DEBUG_LOG("Clearing SLIC data...");
	slicData->Clear();
	
	//## Compute image stats (needed to normalize it)
	if(!img->HasStats()){
		INFO_LOG("Input image has no stats computed, computing them...");
		if(img->ComputeStats(true,false,true)<0) {
			ERROR_LOG("Failed to compute input image stats!");
			return -1;
		}
	}
	
	//## Compute and set input image for SLIC generation (normalized to range [1-256])
	double normmin= 1;
	double normmax= 256;
	INFO_LOG("Normalize input image in range ["<<normmin<<","<<normmax<<"...");
	Image* inputImg= 0;
	if(useLogScaleMapping) {
		inputImg= img->GetNormalizedImage("LOG",normmin,normmax,true);
		//slicData->inputImg= img->GetNormalizedImage("LOG",normmin,normmax,true);
	}
	else {
		inputImg= img->GetNormalizedImage("LINEAR",normmin,normmax,true);
		//slicData->inputImg= img->GetNormalizedImage("LINEAR",normmin,normmax,true);
	}
	if(!slicData->inputImg){
		ERROR_LOG("Failed to compute input norm image!");
		return -1;
	}

	
	//if(!((slicData->inputImg)->HasStats())){
	if( !inputImg->HasStats() ){	
		INFO_LOG("Computing norm image stats...");
		//if((slicData->inputImg)->ComputeStats(true,false,false)<0) {
		if(inputImg->ComputeStats(true,false,false)<0) {
			ERROR_LOG("Failed to compute norm image stats!");
			return -1;
		}
	}//close if hasStats
	

	//## Compute and set laplacian image
	INFO_LOG("Computing laplacian image...");
	//slicData->laplImg= img->GetLaplacianImage(true);
	Image* laplImg= img->GetLaplacianImage(true);
	//if(!slicData->laplImg){
	if(!laplImg){
		ERROR_LOG("Failed to compute lapl image!");
		return -1;
	}
	//if(!((slicData->laplImg)->HasStats())){
	if(!laplImg->HasStats()){
		INFO_LOG("Laplacian image has no stats computed, computing them...");
		//if((slicData->laplImg)->ComputeStats(true,false,false)<0) {	
		if(laplImg->ComputeStats(true,false,false)<0) {
			ERROR_LOG("Failed to compute lapl image stats!");
			return -1;
		}
	}
	
	//## Set edge image (if given in input)
	if(edgeImage){
		//slicData->edgeImg= edgeImage;
		//if(!(slicData->edgeImg)->HasStats()){
		if(!edgeImage->HasStats()){
			INFO_LOG("Edge image has no stats computed, computing them...");
			//if((slicData->edgeImg)->ComputeStats(true,false,false)<0) {
			if(edgeImage->ComputeStats(true,false,false)<0) {
				ERROR_LOG("Failed to compute edge image stats!");
				return -1;
			}
		}
	}//close if is edgeImg

	/*
	//## Init pixel labels
	for(long int i=0;i<img->GetNx();i++){
		(slicData->labels).push_back( std::vector<long int>() );
		for(long int j=0;j<img->GetNy();j++){
			(slicData->labels)[i].push_back(0);
		}//end loop bins Y
	}//end loop bins X
	*/

	//## Finally set SLIC data
	if(slicData->SetData(inputImg,laplImg,edgeImage)<0){
		ERROR_LOG("Failed to set data in SLIC class (hint: check image sizes)!");
		return -1;
	}

	return 0;

}//close SetSPData()


SLICData* SLIC::SPGenerator(Image* img,int regionSize,double regParam, int minRegionSize, bool useLogScaleMapping, Image* edgeImage){

	//## Check inputs
	if(!img){
		ERROR_LOG("Null ptr to given image!");
		return nullptr;
	}
	long long Nx= img->GetNx();
	long long Ny= img->GetNy();
	if(Nx<=0 || Ny<=0 || regionSize<=0 || regParam<0){
		ERROR_LOG("Invalid image size ("<<Nx<<"x"<<Ny<<") or input options (regionSize="<<regionSize<<", regParam="<<regParam<<")");
		return nullptr;
	}

	//## Initialize SLIC storage data
	INFO_LOG("Initialize SLIC storage data...");
	SLICData* slicData= new SLICData;
	if(SetSPData(slicData,img,useLogScaleMapping,edgeImage)<0){
		ERROR_LOG("Initialization failed!");
		delete slicData;
		slicData= 0;
		return nullptr;
	}

	//## Set image range data
	/*
	long int nx= (slicData->inputImg)->GetNx();
	long int ny= (slicData->inputImg)->GetNy();
	float Xmin= (slicData->inputImg)->GetXmin();
	float Ymin= (slicData->inputImg)->GetYmin();
	float Smin= (slicData->inputImg)->GetMinimum();
	float Smax= (slicData->inputImg)->GetMaximum();
	ImgRange imgRange= ImgRange(nx,ny,Smin,Smax,Xmin,Ymin);
	*/

	//## Init SLIC algo data
	INFO_LOG("Allocating SLIC algo data...");
	Image* normImg= slicData->inputImg;
	const unsigned long long maxNumIterations = 100;
	const long long numRegionsX = (unsigned long long) ceil( (double)Nx/regionSize) ;
  const long long numRegionsY = (unsigned long long) ceil( (double)Ny/regionSize) ;
  const long long numRegions = numRegionsX * numRegionsY ;
  const long long numPixels = Nx * Ny;
	
	float* edgeMap = 0;
	edgeMap= (float*)calloc(numPixels, sizeof(float));
	float* centers = 0;
	centers= (float*)malloc(sizeof(float)*3*numRegions);
	unsigned int* masses = 0;
	masses= (unsigned int*)malloc(sizeof(unsigned int) * numPixels);
	if(!edgeMap || !centers || !masses){
		ERROR_LOG("Failed to allocate memory for the algorithm!");
		delete slicData;
		slicData= 0;
		if(masses) free(masses);
  	if(centers) free(centers);
  	if(edgeMap) free(edgeMap);
		return nullptr;
	}

	/*
	unsigned int* segmentation = 0;
	try{
		segmentation= new unsigned int[Nx*Ny];//vector where the segmentation is stored
	}
	catch(...){
		ERROR_LOG("Failed to allocate memory for the algorithm!");
		delete slicData;
		slicData= 0;	
		if(masses) free(masses);
  	if(centers) free(centers);
  	if(edgeMap) free(edgeMap);
		return nullptr;
	}
	*/

	std::vector<long int> segmentation(numPixels,0);
	

	//## Compute edge map (gradient strength)
	INFO_LOG("Compute edge map (gradient strength)...");	
	#ifdef OPENMP_ENABLED
	#pragma omp parallel for collapse(2)
	#endif
	for (long long j=1;j<Ny-1;++j) {
  	for (long long i=1;i<Nx-1;++i) {
			long int iy= j;
			long int ix= i;
			long int index= i+j*Nx;
			double w_left= normImg->GetPixelValue(ix-1,iy);
			double w_right= normImg->GetPixelValue(ix+1,iy);
			double w_up= normImg->GetPixelValue(ix,iy+1);
			double w_down= normImg->GetPixelValue(ix,iy-1);
			double grad= (w_left-w_right)*(w_left-w_right) + (w_up-w_down)*(w_up-w_down);
			edgeMap[index]= grad;
    }//end loop x
  }//end loop y
  
	//## Initialize K-means centers	
	INFO_LOG("Initializing K-means centers...");

 	#ifdef OPENMP_ENABLED
	#pragma omp parallel for
	#endif
	for(long long index=0; index<(signed)numRegions; index++){
		long long u= index % numRegionsX;
		long long v= ((index-u)/numRegionsX) % numRegionsY;
		
		long long centerx = 0;
    long long centery = 0;
    double minEdgeValue = std::numeric_limits<double>::infinity();

		long long x = (long long)round(regionSize * (u + 0.5));
    long long y = (long long)round(regionSize * (v + 0.5));
    x = max(min(x, Nx-1),0LL);
    y = max(min(y, Ny-1),0LL);

		// Search in a 3x3 neighbourhood the smallest edge response
    for (long long yp = max(0LL, y-1); yp <= min(Ny-1, y+1) ; ++yp) {
      for (long long xp = max(0LL, x-1); xp <= min(Nx-1, x+1) ; ++xp) {
        double thisEdgeValue = edgeMap[xp+yp*Nx];
        if (thisEdgeValue < minEdgeValue) {
        	minEdgeValue = thisEdgeValue ;
          centerx = xp ;
          centery = yp ;
        }
      }//end loop xp
    }//end loop yp
		
		// Initialize the new center at this location
		long long array_index= 3*index;
		centers[array_index] = (float) centerx ;
    centers[array_index+1] = (float) centery ;
    centers[array_index+2] = normImg->GetPixelValue(centerx,centery);
		//DEBUG_LOG("(u,v)("<<u<<","<<v<<") (x,y)=("<<x<<","<<y<<"), array_index="<<array_index);
	}
	
	/*
	//=== OLD VERSION NON PARALLELIZED ===
	long long i = 0;
  for (long long v=0; v<(signed)numRegionsY; ++v) {
    for (long long u=0; u<(signed)numRegionsX; ++u) {
      long long centerx = 0;
      long long centery = 0;
      double minEdgeValue = std::numeric_limits<double>::infinity();

      long long x = (long long)round(regionSize * (u + 0.5));
      long long y = (long long)round(regionSize * (v + 0.5));
      x = max(min(x, Nx-1),0LL);
      y = max(min(y, Ny-1),0LL);

      // Search in a 3x3 neighbourhood the smallest edge response
      for (long long yp = max(0LL, y-1); yp <= min(Ny-1, y+1) ; ++yp) {
        for (long long xp = max(0LL, x-1); xp <= min(Nx-1, x+1) ; ++xp) {
          double thisEdgeValue = edgeMap[xp+yp*Nx];
          if (thisEdgeValue < minEdgeValue) {
            minEdgeValue = thisEdgeValue ;
            centerx = xp ;
            centery = yp ;
          }
        }//end loop xp
      }//end loop yp

      // Initialize the new center at this location
			long long index= u+v*numRegionsX;
			INFO_LOG("(u,v)("<<u<<","<<v<<") (x,y)=("<<x<<","<<y<<"), i="<<i<<" index="<<index);
      centers[i++] = (float) centerx ;
      centers[i++] = (float) centery ;
      centers[i++] = normImg->GetPixelValue(centerx,centery);
	
    }//end loop X
  }//end loop Y
	//============================================
	*/

	//## Run k-means iterations
	INFO_LOG("Running "<<maxNumIterations<<" iterations...");
 	double previousEnergy = std::numeric_limits<double>::infinity();
  double startingEnergy= 0;

  for (unsigned long long iter=0;iter < maxNumIterations ; ++iter) {
    double factor = regParam/(regionSize*regionSize);
    double energy = 0 ;

    // assign pixels to centers
		#ifdef OPENMP_ENABLED
		#pragma omp parallel for collapse(2) reduction(+: energy) 
		#endif
    for (long long y = 0 ; y < (signed)Ny ; ++y) {
      for (long long x = 0 ; x < (signed)Nx ; ++x) {
        long long u = floor((double)x/regionSize - 0.5);
        long long v = floor((double)y/regionSize - 0.5);
				float z = normImg->GetPixelValue(x,y);

        double minDistance = std::numeric_limits<double>::infinity();
        for (long long vp = max(0LL, v); vp <= min(numRegionsY-1, v+1) ; ++vp) {
          for (long long up = max(0LL, u); up <= min(numRegionsX-1, u+1) ; ++up) {
            long long region = up  + vp * numRegionsX;
            double centerx = centers[3*region + 0];
            double centery = centers[3*region + 1];
						double centerz = centers[3*region + 2];
						//float z = normImg->GetPixelValue(x,y);//moved out from nested loop

            double spatial = (x - centerx) * (x - centerx) + (y - centery) * (y - centery);
            double appearance = (z - centerz) * (z - centerz);
            double distance= appearance + factor * spatial;
            if (distance<minDistance) {
              minDistance = distance;
              //segmentation[x + y * Nx] = (unsigned int)region;
							segmentation[x + y * Nx] = static_cast<long int>(region);
            }
          }//end loop up
        }//end loop vp
        energy += minDistance;
      }//end loop X
    }//end loop Y

    // check energy termination conditions
    if (iter == 0) {
      startingEnergy = energy ;
    } 
		else {
      if ((previousEnergy - energy) < 1.e-5 * (startingEnergy - energy)) {
        break ;
      }
    }
    previousEnergy = energy;

    // Recompute centers
    memset(masses, 0, sizeof(unsigned int) * Nx * Ny);
    memset(centers, 0, sizeof(float) * 3*numRegions);

    for (long long y = 0 ; y<Ny ; ++y) {
      for (long long x = 0 ; x<Nx ; ++x) {
        long long pixel = x + y * Nx;
        long long region = segmentation[pixel] ;
        masses[region]++;
        centers[region * 3 + 0] += x;
        centers[region * 3 + 1] += y;
				centers[region * 3 + 2] += normImg->GetPixelValue(x,y);
      }//end loop x
    }//end loop y

    for (long long region = 0 ; region < (signed)numRegions ; ++region) {
      double mass = max((double)masses[region], 1e-8) ;
			for (long long i = 3*region; i<3*(region + 1); ++i) {
        centers[i]/= mass;
			}
    }//end loop region

  }//end loop iterations


	//Free stuff
  if(masses) free(masses);
  if(centers) free(centers);
  if(edgeMap) free(edgeMap);

	//## Eliminate small regions
  INFO_LOG("Eliminating small regions...");
  unsigned int* cleaned = (unsigned int*)calloc(numPixels, sizeof(unsigned int));
  unsigned long long* segment = (unsigned long long*)malloc(sizeof(unsigned long long) * numPixels);
	if(!cleaned || !segment){	
		if(cleaned) free(cleaned);
		if(segment) free(segment);
		delete slicData;
		slicData= 0;
		return nullptr;	
	}

 	const long long dx [] = {+1, -1,  0,  0} ;
  const long long dy [] = { 0,  0, +1, -1} ;
 
  for (long long pixel=0; pixel < (signed)numPixels ; ++pixel) {
		if (cleaned[pixel]) continue;
		
    unsigned int label = segmentation[pixel];
    unsigned long long numExpanded = 0 ;
    unsigned long long segmentSize = 0 ;
    segment[segmentSize++] = pixel;

		// Find an adjacent label, for possible use later
    // Find cleanedLabel as the label of an already cleaned region neighbour of this pixel
		unsigned int cleanedLabel = label + 1;
    cleaned[pixel] = label + 1;
    long long x = pixel % Nx;
    long long y = pixel / Nx;
    for (long long direction = 0 ; direction < 4 ; ++direction) {
    	long long xp = x + dx[direction] ;
      long long yp = y + dy[direction] ;
      long long neighbor = xp + yp * Nx ;
      if (0<=xp && xp<(signed)Nx && 0<=yp && yp<(signed)Ny && cleaned[neighbor]) {
      	cleanedLabel = cleaned[neighbor] ;
      }
    }//end loop direction

    // Expand the segment	
		while (numExpanded < segmentSize) {
			long long open = segment[numExpanded++] ;
      x = open % Nx ;
      y = open / Nx ;
      for (long long direction = 0 ; direction < 4 ; ++direction) {
      	long long xp = x + dx[direction] ;
        long long yp = y + dy[direction] ;
        long long neighbor = xp + yp * Nx ;
        if (0<=xp && xp<(signed)Nx && 0<=yp && yp<(signed)Ny && cleaned[neighbor]==0 && segmentation[neighbor]==label) {
        	cleaned[neighbor] = label + 1 ;
          segment[segmentSize++] = neighbor ;
        }//close if
      }//end loop for
    }//end while loop

    //Change label to cleanedLabel if the segment is too small
		if ((signed)segmentSize < minRegionSize) {
			while (segmentSize > 0) {
      	cleaned[segment[--segmentSize]] = cleanedLabel ;
      }//end while loop
    }//close if 

  }//end pixel loop

	//## Restore base 0 indexing of the regions
  //for (long long pixel=0; pixel<(signed)numPixels; ++pixel) cleaned[pixel]--;
  //memcpy(segmentation, cleaned, numPixels * sizeof(unsigned int)) ;
		
	for (long long pixel=0; pixel<(signed)numPixels; ++pixel) {
		cleaned[pixel]--;
		segmentation[pixel]= cleaned[pixel];
	}

	//## Delete data
	DEBUG_LOG("Free data...");
	if(cleaned) free(cleaned);
  if(segment) free(segment);

	
	//## Allocate regions
	INFO_LOG("Allocating "<<numRegions<<" regions list...");
	std::vector<Region*> regions;
	Region* aRegion= 0;
	for(long long k=0;k<numRegions;k++) {		
		aRegion= new Region();
		regions.push_back(aRegion);
	}

	
	//## Fill regions
	INFO_LOG("Filling regions...");
	Pixel* aPixel= 0;

	for (long int ix=0; ix<Nx; ix++) {
  	for (long int iy=0;iy<Ny; iy++) {

			//long int id= ix + iy * Nx;
			long int id= img->GetBin(ix,iy);
			double x= img->GetX(ix);
			double y= img->GetY(iy);
			double S= img->GetPixelValue(ix,iy);
			double S_curv= slicData->GetScurv(ix,iy);//(slicData->laplImg)->GetPixelValue(ix,iy);
			double S_edge= slicData->GetSedge(ix,iy);//if(slicData->edgeImg) S_edge= (slicData->edgeImg)->GetPixelValue(ix,iy);

			long int label= static_cast<long int>(segmentation[id]);
			
			//slicData->SetPixelLabel(ix,iy,label);
			//(slicData->labels)[ix][iy]= label;//TO BE REMOVED AFTER CHECKS!

			if(label<0 || label>=(signed)numRegions){
				WARN_LOG("Skip label "<<label<<" for pixel ("<<ix<<","<<iy<<")...");
				continue;
			}

			//Create and fill a pixel
			aPixel= new Pixel;
			aPixel->id= id;
			aPixel->SetPhysCoords(x,y);
			aPixel->SetCoords(ix,iy);
			aPixel->S= S;
			aPixel->SetCurv(S_curv);
			aPixel->SetEdge(S_edge);
		
			//Add pixel to region
			//slicData->SetRegionId(label,label);
			//if(S!=0) slicData->AddPixelToRegion(label,aPixel);
		
			regions[label]->Id= label;
			if(S!=0) regions[label]->AddPixel(aPixel);		

			//(slicData->regions)[label]->Id= label;
			//if(pixelVal!=0) (slicData->regions)[pixelLabel]->AddPixel(aPixel);
    }//end loop Ny
  }//end loop Nx

	//Clear segmentation list
	//if(segmentation) delete [] segmentation; 	
	
	//## Set regions and labels to slic data
	INFO_LOG("Setting regions and pixel labels in slic data...");
	slicData->SetRegions(regions);
	
	
	//## Remove regions without pixels (important in case of empty image zones)
	INFO_LOG("Removing regions without pixels...");
	slicData->RemoveEmptyRegions();

	/*
	std::vector<size_t> regionsToBeDeleted;
	std::map<long int,long int> goodRegionIdMap;
	for(size_t k=0;k<(slicData->regions).size();k++) {	
		long int regionId= slicData->GetRegionId(k);//(slicData->regions)[k]->Id;
		long int nPix= slicData->GetRegionSize(k);//(slicData->regions)[k]->NPix;
		goodRegionIdMap.insert( std::pair<long int,long int>(regionId,1) );
		if(nPix<=0) {
			regionsToBeDeleted.push_back(k);
			goodRegionIdMap[regionId]= -1;
		}
	}
	for (long int ix=0; ix<Nx; ix++) {
		for (long int iy=0;iy<Ny; iy++) {
			long int index= ix + iy * Nx;
			long int labelId= slicData->GetPixelLabel(index);
			//long int labelId= (slicData->labels)[ix][iy];
			long int corrFactor= goodRegionIdMap[labelId];
			
			slicData->ScalePixelLabel(index,corrFactor);
			(slicData->labels)[ix][iy]*= corrFactor; 
		}
	}
	slicData->DeleteRegions(regionsToBeDeleted);//CodeUtils::DeleteItems(slicData->regions, regionsToBeDeleted);
	*/
	//##########################################################################

	INFO_LOG("#"<<slicData->GetNRegions()<<" regions present after cleanup...");

	//## Compute region parameters
	INFO_LOG("Computing region parameters...");
	if(slicData->ComputeRegionParameters()<0){
		WARN_LOG("One/more errors occurred while computing region parameters!");
	}

	/*
	for(size_t k=0;k<(slicData->regions).size();k++) {	
		long int nPix= (slicData->regions)[k]->NPix;
		if(nPix<=0) continue;
		(slicData->regions)[k]->ComputeStats(true,false);
	}//end loop regions
	*/

	INFO_LOG("End superpixel generation...");

	return slicData;

}//close SLIC::SuperpixelGenerator()




}//close namespace

