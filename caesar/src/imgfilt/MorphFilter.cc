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
* @file MorphFilter.cc
* @class MorphFilter
* @brief Class implementing morphological filtering
*
* Morphological Filter
* @author S. Riggi
* @date 20/01/2015
*/


#include <MorphFilter.h>
#include <Image.h>
#include <MathUtils.h>
#include <Source.h>
#include <Pixel.h>
#include <BkgData.h>
#include <Logger.h>

#include <TObject.h>
#include <TRInterface.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>

using namespace std;

ClassImp(Caesar::MorphFilter)

namespace Caesar {

MorphFilter::MorphFilter() {

}//close costructor


MorphFilter::~MorphFilter(){

}//close destructor

//================================================
//==      NEW IMAGE METHODS
//================================================
/*
Image* MorphFilter::Dilate(Image* img,int KernSize,bool returnPeakImg){
		
	if(!img){
		ERROR_LOG("Null ptr to given image!");
		return nullptr;
	}
	long int Nx= img->GetNx();
	long int Ny= img->GetNy();

	//## Convert image to OpenCV format
	cv::Mat mat= img->GetOpenCVMat("64");

	//## Init dilation options
	cv::Size kernel_size(KernSize,KernSize);
	cv::Mat element= cv::getStructuringElement(cv::MORPH_RECT, kernel_size, cv::Point(-1,-1));
	
	//## Dilate image
	cv::Mat mat_dilated;
	int iterations= 1;
	cv::dilate(mat, mat_dilated, element, cv::Point(-1,-1),iterations,cv::BORDER_CONSTANT);

	//## Compare original and dilated image
	cv::Mat mat_cmp = cv::Mat::zeros(Ny,Nx,CV_8UC1);
	cv::compare(mat, mat_dilated, mat_cmp, cv::CMP_EQ);
	
	//## Convert back dilated image 
	TString imgName= Form("%s_Dilated_kernSize%d",img->GetName().c_str(),KernSize);	
	Img* DilatedImg= img->GetCloned(std::string(imgName),true,true);
	DilatedImg->Reset();
	
	imgName= Form("%s_DilatePeak_kernSize%d",img->GetName().c_str(),KernSize);	
	Img* peakImg= img->GetCloned(std::string(imgName),true,true);
	peakImg->Reset();
	
	int npeaks= 0;
	
	for(long int j=0;j<Ny;j++){
		int rowId= Ny-1-j;
		double y= img->GetYaxis()->GetBinCenter(j+1);

		for(long int i=0;i<Nx;i++){
			int colId= i;
			double x= img->GetXaxis()->GetBinCenter(i+1);
			
			double binContent= img->GetBinContent(i+1,j+1);
			if(binContent==0) continue;
			double matrixElement= mat_dilated.at<double>(rowId,colId);
				
			DilatedImg->FillPixel(x,y,matrixElement);

			float mat_comparison= (float)mat_cmp.at<unsigned char>(rowId,colId);
			if(mat_comparison==255){
				//## Check surrounding pixels (do not tag as peak if the surrounding 3x3 is flat)
				bool isFlatArea= true;
				for(int ix=i-1;ix<i+1;ix++){
					for(int iy=j-1;iy<j+1;iy++){
						if(ix==i && iy==j) continue;
						double w= img->GetBinContent(ix+1,iy+1);
						if(w!=binContent) {
							isFlatArea= false;
							break;
						}
					}//end loop kernel y
				}//end loop kernel x

				if(!isFlatArea){
					peakImg->FillPixel(x,y,1);
					npeaks++;
					DEBUG_LOG("Peaks #"<<npeaks<<" detected @ ("<<x<<","<<y<<")");
				}
			}
		}//end loop x
	}//end loop y
	
	Img* returnImg= 0;
	if(returnPeakImg){
		DilatedImg->Delete();
		returnImg= peakImg;
	}
	else{
		peakImg->Delete();
		returnImg= DilatedImg;
	}

	return returnImg;

}//close Dilate()
*/

int MorphFilter::FindDilatedSourcePixels(Image* img,Source* source,int KernSize,std::vector<long int>& pixelsToBeDilated){

	//## Check input source
	if(!source) {
		ERROR_LOG("Null ptr to input source given!");
		return -1;
	}
	
	//## Init dilation kernel
	if(KernSize%2==0){
		ERROR_LOG("KernSize argument should be an odd number!");
		return -1;
	}
	int dilateSize= KernSize/2;
	cv::Mat element= cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(KernSize,KernSize));
	
	//PixelCollection sourcePixels= source->GetPixels();
	int nDilatedPixels= 0;
	long int Nx= img->GetNx();
	long int Ny= img->GetNy();

	for(int l=0;l<source->GetNPixels();l++){
		Pixel* thisPixel= source->GetPixel(l);
		long int id= thisPixel->id;
		long int binx= img->GetBinX(id);
		long int biny= img->GetBinY(id);
		long int ix= binx;
		long int iy= biny;			

		//Find if pixel was already added to the list of dilated pixels
		std::vector<long int>::iterator it= std::find(pixelsToBeDilated.begin(),pixelsToBeDilated.end(),id);
		if( (it==pixelsToBeDilated.end() || pixelsToBeDilated.empty()) ) pixelsToBeDilated.push_back(id);
			
		for(int tx=-dilateSize;tx<=dilateSize;tx++){
			long int binx_next= binx + tx;
			long int colId= tx + dilateSize;

			for(int ty=-dilateSize;ty<=dilateSize;ty++){	
				long int biny_next= biny+ty;
				long int rowId= tx + dilateSize;

				if(ix+tx<Nx && ix+tx>=0 && iy+ty<Ny && iy+ty>=0){
					double kernValue= (double)element.at<char>(rowId,colId);
					long int gBinId= img->GetBin(binx_next,biny_next);
					it= std::find(pixelsToBeDilated.begin(),pixelsToBeDilated.end(),gBinId);
					if( kernValue>0 && (it==pixelsToBeDilated.end() || pixelsToBeDilated.empty()) ) {
						pixelsToBeDilated.push_back(gBinId);
						nDilatedPixels++;
					}
				}//close if
			}//end loop kernel
		}//end loop kernel
	}//end loop pixels		
	
	//## Replace selected pixels
	DEBUG_LOG("#"<<nDilatedPixels<<" pixels to be dilated...");
		
	return 0;

}//close FindDilatedSourcePixels()


int MorphFilter::DilateAroundSource(Image* img,Source* source,int KernSize,int dilateModel,int dilateSourceType,bool skipToNested,ImgBkgData* bkgData,bool useLocalBkg,bool randomize,double zThr){

	//## Check input source
	if(!source) {
		ERROR_LOG("Null ptr to input source given!");
		return -1;
	}
	bool hasNestedSources= source->HasNestedSources();
	bool hasStats= source->HasStats();
	if(!hasStats){
		WARN_LOG("No stats computed for input source...computing!");
		source->ComputeStats(true,true);
	}
	int sourceType= source->Type;
	double sourceMedian= source->Median;
	double sourceMedianRMS= source->MedianRMS;
	
	//## Skip faint sources
	//Get pixel seeds
	bool isFaintSource= false;
	std::vector<int> seedPixelIndexes= source->GetSeedPixelIndexes();
	DEBUG_LOG("#"<<seedPixelIndexes.size()<<" seed pixels...");

	double Zmax= -1.e+99;

	#ifdef OPENMP_ENABLED
	#pragma omp parallel for reduction(max: Zmax)
	#endif
	for(size_t i=0;i<seedPixelIndexes.size();i++){
		int index= seedPixelIndexes[i];
		Pixel* aPixel= source->GetPixel(index);
		if(!aPixel) continue;
		long int id= aPixel->id;
		std::pair<double,double> bkgInfo= aPixel->GetBkg();
		double bkgLevel= bkgInfo.first;
		double noiseLevel= bkgInfo.second;
		double w= img->GetBinContent(id);
		double Z= (w-bkgLevel)/noiseLevel;
		if(fabs(Z)>Zmax) Zmax= Z;
	}

	if(Zmax<zThr){
		INFO_LOG("Source is below significance threshold for dilation (Z="<<Zmax<<"<"<<zThr<<"), skip it!");
		return 0;	
	}
	
	//## Check R interface
	double sigmaTrunc= 1;//trunc random gaussian to +-sigmaTrunc	

	DEBUG_LOG("Retrieve RInterface instance...");
	ROOT::R::TRInterface& fR= ROOT::R::TRInterface::Instance();
	std::string randomGenCmd= std::string("rtruncnorm(1, a=-sigmaTrunc, b=sigmaTrunc, mean = 0, sd = 1)");
	try{
		fR.Execute("library(\"truncnorm\");");
		fR["sigmaTrunc"]= sigmaTrunc;
	}
	catch( std::exception &ex ) {
		ERROR_LOG("C++ exception catched while loading R library truncnorm (err=" << ex.what() <<")");
		return -1;
  } 
	catch(...) { 
		ERROR_LOG("Unknown exception catched while loading R library truncnorm!");
		return -1;
  }	


	//## Find pixels to be dilated
	//== MOTHER SOURCE
	std::vector<long int> pixelsToBeDilated;
	DEBUG_LOG("hasNestedSources="<<hasNestedSources);
	if(!hasNestedSources && (dilateSourceType==-1 || sourceType==dilateSourceType) ){
		DEBUG_LOG("Dilating mother sources...");
		FindDilatedSourcePixels(img,source,KernSize,pixelsToBeDilated);	
	}
	//== NESTED SOURCES
	if(skipToNested && hasNestedSources){
		DEBUG_LOG("Dilating nested sources...");
		std::vector<Source*> nestedSources= source->GetNestedSources();
		for(unsigned int k=0;k<nestedSources.size();k++){
			int nestedSourceType= nestedSources[k]->Type;
			if(dilateSourceType==-1 || nestedSourceType==sourceType){
				FindDilatedSourcePixels(img,nestedSources[k],KernSize,pixelsToBeDilated);
			}
		}//end loop nested sources
	}//close if dilate nested

	//## Replace dilated pixels with model		
	if(dilateModel==eDilateWithSourceMedian){
		double BkgRealization= sourceMedian;
		double BkgRMS= sourceMedianRMS;
		if(randomize){
			#ifdef OPENMP_ENABLED
			#pragma omp parallel for
			#endif
			for(size_t l=0;l<pixelsToBeDilated.size();l++){
				long int id= pixelsToBeDilated[l];			
				double r= fR.Eval(randomGenCmd.c_str());
				double bkg= BkgRealization + r*BkgRMS;
				//BkgRealization+= r*BkgRMS;
				img->SetPixelValue(id,bkg);
			}//end loop pixels 	
		}
		else{
			#ifdef OPENMP_ENABLED
			#pragma omp parallel for
			#endif
			for(size_t l=0;l<pixelsToBeDilated.size();l++){
				long int id= pixelsToBeDilated[l];			
				img->SetPixelValue(id,BkgRealization);
			}//end loop pixels 
		}
	}//close if
  else if(dilateModel==eDilateWithBkg){
		if(useLocalBkg){
			if(randomize){
				#ifdef OPENMP_ENABLED
				#pragma omp parallel for
				#endif
				for(size_t l=0;l<pixelsToBeDilated.size();l++){
					long int id= pixelsToBeDilated[l];			
					double BkgRealization= (bkgData->BkgMap)->GetPixelValue(id);
					double BkgRMS= (bkgData->NoiseMap)->GetPixelValue(id);
					double r= fR.Eval(randomGenCmd.c_str());
					//BkgRealization+= r*BkgRMS;
					double bkg= BkgRealization + r*BkgRMS;
					img->SetPixelValue(id,bkg);
				}//end loop pixels
			}
			else{
				#ifdef OPENMP_ENABLED
				#pragma omp parallel for
				#endif
				for(size_t l=0;l<pixelsToBeDilated.size();l++){
					long int id= pixelsToBeDilated[l];			
					double BkgRealization= (bkgData->BkgMap)->GetPixelValue(id);
					img->SetPixelValue(id,BkgRealization);
				}//end loop pixels
			}
		}//close if
		else{
			double BkgRealization= bkgData->gBkg;
			double BkgRMS= bkgData->gNoise;	
			if(randomize){
				#ifdef OPENMP_ENABLED
				#pragma omp parallel for
				#endif
				for(size_t l=0;l<pixelsToBeDilated.size();l++){
					long int id= pixelsToBeDilated[l];			
					double r= fR.Eval(randomGenCmd.c_str());
					double bkg= BkgRealization + r*BkgRMS;
					//BkgRealization+= r*BkgRMS;
					img->SetPixelValue(id,bkg);
				}//end loop pixels
			}
			else{
				#ifdef OPENMP_ENABLED
				#pragma omp parallel for
				#endif
				for(size_t l=0;l<pixelsToBeDilated.size();l++){
					long int id= pixelsToBeDilated[l];			
					img->SetPixelValue(id,BkgRealization);
				}//end loop pixels
			}
		}//close else
	}//close else if

	//Force recomputation of stats if present, otherwise recompute only moments
	bool skipNegativePixels= false;
	bool computeRobustStats= true;	
	bool forceRecomputing= true;
	int status= 0;
	if(img->HasStats()) status= img->ComputeStats(computeRobustStats,skipNegativePixels,forceRecomputing);
	else status= img->ComputeMoments(skipNegativePixels);
		
	if(status<0){
		WARN_LOG("Failed to recompute moments/stats after source dilation!");
		return -1;
	}

	return 0;

}//close DilateAroundSource()


int MorphFilter::DilateAroundSources(Image* img,std::vector<Source*>const& sources,int KernSize,int dilateModel,int dilateSourceType,bool skipToNested,ImgBkgData* bkgData,bool useLocalBkg,bool randomize,double zThr){
	
	//## Check input image
	if(!img){
		ERROR_LOG("Null ptr to given image!");
		return -1;
	}
	
	//## Check bkg data
	if(dilateModel==eDilateWithBkg){
	 	if(!bkgData){
			ERROR_LOG("Selected to use bkg dilation but null ptr to bkg data!");
			return -1;
		}
		if(useLocalBkg && !bkgData->HasLocalBkg()){
			ERROR_LOG("Selected to use local bkg but no local bkg data are available!");
			return -1;
		}
	}//close if

	//## Check source list
	if(sources.size()<=0){
		WARN_LOG("Source list empty, nothing to be dilated!");
		return 0;
	}

	//## Start dilating sources
	for(unsigned int k=0;k<sources.size();k++){	
		int status= DilateAroundSource(img,sources[k],KernSize,dilateModel,dilateSourceType,skipToNested,bkgData,useLocalBkg,randomize,zThr);
		if(status<0){
			WARN_LOG("Source dilation failed for source no. "<<k<<" ...");
		}
	}//end loop sources

	return 0;

}//close DilateAroundSources()


}//close namespace
