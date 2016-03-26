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
#include <Img.h>
#include <MathUtils.h>
#include <Source.h>
#include <Pixel.h>
#include <BkgData.h>

#include <TObject.h>
#include <TRInterface.h>

//#include <Rcpp.h>
//#include <RInside.h>
//using namespace Rcpp;

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

Img* MorphFilter::Dilate(Img* img,int KernSize,bool returnPeakImg){
		
	if(!img){
		cerr<<"MorphFilter::Dilate(): Null ptr to given image!"<<endl;
		return 0;
	}
	int Nx= img->GetNbinsX();
	int Ny= img->GetNbinsY();

	//## Convert image to OpenCV format
	cv::Mat mat= img->ImgToMat("64");

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
	TString imgName= Form("%s_Dilated_kernSize%d",std::string(img->GetName()).c_str(),KernSize);	
	Img* DilatedImg= img->GetCloned(std::string(imgName),true,true);
	DilatedImg->Reset();
	
	imgName= Form("%s_DilatePeak_kernSize%d",std::string(img->GetName()).c_str(),KernSize);	
	Img* peakImg= img->GetCloned(std::string(imgName),true,true);
	peakImg->Reset();
	
	int npeaks= 0;
	
	for(int j=0;j<Ny;j++){
		int rowId= Ny-1-j;
		double y= img->GetYaxis()->GetBinCenter(j+1);

		for(int i=0;i<Nx;i++){
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
					//peaks.push_back(binId);
					peakImg->FillPixel(x,y,1);
					npeaks++;
					cout<<"Img::Dilate(): INFO: Peaks #"<<npeaks<<" detected @ ("<<x<<","<<y<<")"<<endl;
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

}//close MorphFilter::Dilate()


int MorphFilter::FindDilatedSourcePixels(Img* img,Source* source,int KernSize,std::vector<int>& pixelsToBeDilated){

	//## Check input source
	if(!source) {
		cerr<<"MorphFilter::DilateAroundSource(): ERROR: Null ptr to input source given!"<<endl;
		return -1;
	}
	
	//## Init dilation kernel
	if(KernSize%2==0){
		cerr << "MorphFilter::DilateAroundSource(): ERROR: KernSize should be an odd number!" << endl;
		return -1;
	}
	int dilateSize= KernSize/2;
	cv::Mat element= cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(KernSize,KernSize));
	
	//PixelCollection sourcePixels= source->GetPixels();
	int nDilatedPixels= 0;
	int Nx= img->GetNbinsX();
	int Ny= img->GetNbinsY();

	//for(unsigned int l=0;l<sourcePixels.size();l++){
	for(int l=0;l<source->GetNPixels();l++){
		//Pixel* thisPixel= sourcePixels[l];
		Pixel* thisPixel= source->GetPixel(l);
		int id= thisPixel->id;
		int binx, biny, binz;
		img->GetBinXYZ(id,binx,biny,binz);
		int ix= binx-1;
		int iy= biny-1;			

		std::vector<int>::iterator it= std::find(pixelsToBeDilated.begin(),pixelsToBeDilated.end(),id);
		if( (it==pixelsToBeDilated.end() || pixelsToBeDilated.empty()) ) pixelsToBeDilated.push_back(id);
			
		for(int tx=-dilateSize;tx<=dilateSize;tx++){
			int binx_next= binx+tx;
			int colId= tx + dilateSize;
			for(int ty=-dilateSize;ty<=dilateSize;ty++){	
				int biny_next= biny+ty;
				int rowId= tx + dilateSize;
				if(ix+tx<Nx && ix+tx>=0 && iy+ty<Ny && iy+ty>=0){
					double kernValue= (double)element.at<char>(rowId,colId);
					int gBinId= img->GetBin(binx_next,biny_next);
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
	cout<<"MorphFilter::FindDilatedSourcePixels(): INFO: #"<<nDilatedPixels<<" pixels to be dilated..."<<endl;
		
	return 0;

}//close FindDilatedSourcePixels()

int MorphFilter::DilateAroundSource(Img* img,Source* source,int KernSize,int dilateModel,int dilateSourceType,bool skipToNested,BkgData* bkgData,bool useLocalBkg,bool randomize){

	//## Check input source
	if(!source) {
		cerr<<"MorphFilter::DilateAroundSource(): ERROR: Null ptr to input source given!"<<endl;
		return -1;
	}
	bool hasNestedSources= source->HasNestedSources();
	bool hasStats= source->HasStats();
	if(!hasStats){
		cerr<<"MorphFilter::DilateAroundSource(): WARN: No stats computed for input source...computing!"<<endl;
		source->ComputeStats(true,true);
	}
	int sourceType= source->Type;
	double sourceMedian= source->Median;
	double sourceMedianRMS= source->MedianRMS;
	
	//## Check R interface
	double sigmaTrunc= 1;//trunc random gaussian to +-sigmaTrunc	

	
	//RInside* fR= 0;
	cout<<"MorphFilter::DilateAroundSource(): INFO: Retrieve RInterface instance..."<<endl;
	ROOT::R::TRInterface& fR= ROOT::R::TRInterface::Instance();
	cout<<"MorphFilter::DilateAroundSource(): INFO: done!"<<endl;
	std::string randomGenCmd= std::string("rtruncnorm(1, a=-sigmaTrunc, b=sigmaTrunc, mean = 0, sd = 1)");
	try{
		/*
		//Get instance to RInside
		fR= RInside::instancePtr();
		if(!fR){
			cerr<<"MorphFilter::DilateAroundSource(): WARN: Cannot retrieve RInside instance (did you create one in your main application?)!"<<endl;
			fR= new RInside;		
		} 
		(*fR)["sigmaTrunc"]= sigmaTrunc;	
		fR->parseEvalQ( std::string("library(\"truncnorm\");") );
		*/
		fR.Execute("library(\"truncnorm\");");
		fR["sigmaTrunc"]= sigmaTrunc;

	}//close try block
	
	catch( std::exception &ex ) {
		cerr << "MorphFilter::DilateAroundSource(): ERROR: Exception catched: " << ex.what() << endl;
		return -1;
  } 
	catch(...) { 
		cerr << "MorphFilter::DilateAroundSource(): ERROR: C++ exception (unknown reason)" << endl;
		return -1;
  }	


	//## Find pixels to be dilated
	//== MOTHER SOURCE
	std::vector<int> pixelsToBeDilated;
	cout<<"MorphFilter::DilateAroundSource(): INFO: hasNestedSources="<<hasNestedSources<<endl;
	if(!hasNestedSources && (dilateSourceType==-1 || sourceType==dilateSourceType) ){
		cout<<"MorphFilter::DilateAroundSource(): INFO: Dilating mother sources..."<<endl;
		FindDilatedSourcePixels(img,source,KernSize,pixelsToBeDilated);	
	}
	//== NESTED SOURCES
	if(skipToNested && hasNestedSources){
		cout<<"MorphFilter::DilateAroundSource(): INFO: Dilating nested sources..."<<endl;
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
			for(unsigned int l=0;l<pixelsToBeDilated.size();l++){
				int id= pixelsToBeDilated[l];			
				double r= fR.Eval(randomGenCmd.c_str());
				//double r= Rcpp::as<double>( fR->parseEval(randomGenCmd) );
				BkgRealization+= r*BkgRMS;
				img->SetBinContent(id,BkgRealization);
			}//end loop pixels 	
		}
		else{
			for(unsigned int l=0;l<pixelsToBeDilated.size();l++){
				int id= pixelsToBeDilated[l];			
				img->SetBinContent(id,BkgRealization);
			}//end loop pixels 
		}
	}//close if
  else if(dilateModel==eDilateWithBkg){
		if(useLocalBkg){
			if(randomize){
				for(unsigned int l=0;l<pixelsToBeDilated.size();l++){
					int id= pixelsToBeDilated[l];			
					double BkgRealization= (bkgData->BkgMap)->GetBinContent(id);
					double BkgRMS= (bkgData->NoiseMap)->GetBinContent(id);
					//double r= Rcpp::as<double>( fR->parseEval(randomGenCmd) );
					double r= fR.Eval(randomGenCmd.c_str());
					BkgRealization+= r*BkgRMS;
					img->SetBinContent(id,BkgRealization);
				}//end loop pixels
			}
			else{
				for(unsigned int l=0;l<pixelsToBeDilated.size();l++){
					int id= pixelsToBeDilated[l];			
					double BkgRealization= (bkgData->BkgMap)->GetBinContent(id);
					img->SetBinContent(id,BkgRealization);
				}//end loop pixels
			}
		}//close if
		else{
			double BkgRealization= bkgData->gBkg;
			double BkgRMS= bkgData->gNoise;	
			if(randomize){
				for(unsigned int l=0;l<pixelsToBeDilated.size();l++){
					int id= pixelsToBeDilated[l];			
					//double r= Rcpp::as<double>( fR->parseEval(randomGenCmd) );
					double r= fR.Eval(randomGenCmd.c_str());
					BkgRealization+= r*BkgRMS;
					img->SetBinContent(id,BkgRealization);
				}//end loop pixels
			}
			else{
				for(unsigned int l=0;l<pixelsToBeDilated.size();l++){
					int id= pixelsToBeDilated[l];			
					img->SetBinContent(id,BkgRealization);
				}//end loop pixels
			}
		}//close else
	}//close else if

	return 0;

}//close DilateAroundSource()


int MorphFilter::DilateAroundSources(Img* img,std::vector<Source*>const& sources,int KernSize,int dilateModel,int dilateSourceType,bool skipToNested,BkgData* bkgData,bool useLocalBkg,bool randomize){
	
	//## Check input image
	if(!img){
		cerr<<"MorphFilter::DilateAroundSource(): ERROR: Null ptr to given image!"<<endl;
		return -1;
	}
	
	//## Check bkg data
	if(dilateModel==eDilateWithBkg){
	 	if(!bkgData){
			cerr<<"MorphFilter::DilateAroundSources(): ERROR: Selected to use bkg dilation but null ptr to bkg data!"<<endl;
			return -1;
		}
		if(useLocalBkg && !bkgData->HasLocalBkg()){
			cerr<<"MorphFilter::DilateAroundSources(): ERROR: Selected to use local bkg but no local bkg data are available!"<<endl;
			return -1;
		}
	}//close if

	//## Check source list
	if(sources.size()<=0){
		cerr<<"MorphFilter::DilateAroundSources(): WARN: Source list empty, nothing to be dilated!"<<endl;
		return 0;
	}

	//## Start dilating sources
	for(unsigned int k=0;k<sources.size();k++){	
		int status= DilateAroundSource(img,sources[k],KernSize,dilateModel,dilateSourceType,skipToNested,bkgData,useLocalBkg,randomize);
		if(status<0){
			cerr<<"MorphFilter::DilateAroundSources(): WARN: Source dilation failed for source no. "<<k<<" ..."<<endl;				
		}
	}//end loop sources

	return 0;

}//close DilateAroundSources()

}//close namespace
