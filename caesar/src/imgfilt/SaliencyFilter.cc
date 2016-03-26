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
* @file SaliencyFilter.cc
* @class SaliencyFilter
* @brief Class implementing saliency filtering
*
* Saliency Filter
* @author S. Riggi
* @date 20/01/2015
*/

#include <SaliencyFilter.h>
#include <Img.h>
#include <Region.h>
#include <SLIC.h>
#include <BkgData.h>
#include <CodeUtils.h>

#include <TObject.h>
#include <TVectorD.h>
#include <TGraph2D.h>


#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>

using namespace std;

ClassImp(Caesar::SaliencyFilter)

namespace Caesar {

SaliencyFilter::SaliencyFilter() {

}//close costructor


SaliencyFilter::~SaliencyFilter(){

}//close destructor


Img* SaliencyFilter::ComputeSaliencyMap(Img* img,int reso,double regFactor,int minRegionSize,double knnFactor,double spatialRegFactor,bool useRobust,bool addCurvDist){

	//## Normalize image
	double NormMin= 1;
	double NormMax= 256;
	Img* img_norm= img->GetNormalizedImage("LINEAR",NormMin,NormMax);
	if(!img_norm){
		cerr<<"SaliencyFilter::ComputeSaliencyMap(): ERROR: Failed to normalize input image!"<<endl;
		return 0;
	}

	//## Compute segmentation in superpixels
	bool useLogScaleMapping= false;
	//SLICData* slicData= SLIC::SPGenerator(img,reso,regFactor,minRegionSize,useLogScaleMapping,0);	//pass null edgeImg (not needed here)
	SLICData* slicData= SLIC::SPGenerator(img_norm,reso,regFactor,minRegionSize,useLogScaleMapping,0);	//pass null edgeImg (not needed here)
	if(!slicData){
		cerr<<"SaliencyFilter::ComputeSaliencyMap(): ERROR: Superpixel segmentation failed!"<<endl;
		return 0;
	}

	//Get results
	std::vector<Region*> regions= (slicData->regions);
	int nRegions= (int)regions.size();
	if(nRegions<=0) {
		cerr<<"SaliencyFilter::ComputeSaliencyMap(): ERROR: Superpixel segmentation returned no regions!"<<endl;
		delete slicData;
		slicData= 0;
		return 0;
	}

	//## Compute saliency map
	//Img* saliencyMap= ComputeSaliencyMap(img,regions,knnFactor,useRobust,addCurvDist);
	Img* saliencyMap= ComputeSaliencyMap(img_norm,regions,knnFactor,spatialRegFactor,useRobust,addCurvDist);
	if(!saliencyMap){
		cerr<<"SaliencyFilter::ComputeSaliencyMap(): ERROR: Failed to compute saliency map!"<<endl;
		delete slicData;
		slicData= 0;
		return 0;
	}

	img_norm->Delete();
	delete slicData;
	slicData= 0;
	
	return saliencyMap;

}//close ComputeSaliencyMap()


Img* SaliencyFilter::ComputeSaliencyMap(Img* img,std::vector<Region*>const& regions,double knnFactor,double spatialRegFactor,bool useRobust,bool addCurvDist){

	//## Check regions
	int nRegions= (int)regions.size();
	if(nRegions<=0) {
		cerr<<"SaliencyFilter::ComputeSaliencyMap(): ERROR: No regions given!"<<endl;
		return 0;
	}
	
	//## Select number of nearest neighbor regions to be used in saliency computation
	int knn_min= 10;//in number of regions
	int knn_chosen= (int)(std::round(knnFactor*nRegions));//percentage of available regions
	int knn= std::max(knn_chosen,knn_min);
	int KNN= knn;
	if(knn>nRegions || knn<0) KNN= nRegions;
	cout<<"KNN="<<KNN<<endl;
	
	//## Compute region parameters (including robust stats)
	int nPars= 2;
	int nPars_robust= 2;
	if(addCurvDist){
		nPars+= 2;
		nPars_robust+= 2;
	}
	int nPars_spatial= 2;
	TMatrixD params(nRegions,nPars);
	params.Zero();
	TMatrixD params_robust(nRegions,nPars_robust);
	params_robust.Zero();
	TMatrixD params_spatial_noNorm(nRegions,nPars_spatial);
	params_spatial_noNorm.Zero();
	TMatrixD params_spatial(nRegions,nPars_spatial);
	params_spatial.Zero();

	cout<<"SaliencyFilter::ComputeSaliencyMap(): INFO: Compute region pars..."<<endl;
	for(int i=0;i<nRegions;i++){
		//Compute stats pars
		if(!regions[i]->HasStats()){
			regions[i]->ComputeStats(true,false);
		}
		//regions[i]->Print();

		//Get params vectors
		Region::RegionPars* regionPars= regions[i]->GetParams(addCurvDist);
		if(!regionPars){
			cerr<<"SaliencyFilter::ComputeSaliencyMap(): ERROR: Failed to get pars for region no. "<<i<<"!"<<endl;
			return 0;
		}
		
		//Fill param matrix
		params[i]= *(regionPars->pars);
		params_robust[i]= *(regionPars->robustPars);
		params_spatial[i]= *(regionPars->spatialPars);
		params_spatial_noNorm[i]= *(regionPars->spatialPars);

		//Clear up
		if(regionPars){
			delete regionPars;
			regionPars= 0;
		}
	}//end loop regions

	//## Find min & max
	cout<<"SaliencyFilter::ComputeSaliencyMap(): INFO: Compute param min/max..."<<endl;
	TVectorD params_min(nPars);
	TVectorD params_max(nPars);
	params_min.Zero();
	params_max.Zero();
	for(int j=0;j<nPars;j++){
		TVectorD v = TMatrixDColumn(params,j);
		params_min(j)= v.Min();
		params_max(j)= v.Max();
	}//end loop pars

	TVectorD params_robust_min(nPars_robust);
	TVectorD params_robust_max(nPars_robust);
	params_robust_min.Zero();
	params_robust_max.Zero();
	for(int j=0;j<nPars_robust;j++){
		TVectorD v = TMatrixDColumn(params_robust,j);
		params_robust_min(j)= v.Min();
		params_robust_max(j)= v.Max();
	}//end loop pars

	TVectorD params_spatial_min(nPars_spatial);
	TVectorD params_spatial_max(nPars_spatial);
	params_spatial_min.Zero();
	params_spatial_max.Zero();
	for(int j=0;j<nPars_spatial;j++){
		TVectorD v = TMatrixDColumn(params_spatial,j);
		params_spatial_min(j)= v.Min();
		params_spatial_max(j)= v.Max();
	}//end loop pars

	//## Normalize parameters to [0,1]	
	cout<<"SaliencyFilter::ComputeSaliencyMap(): INFO: Normalize region pars.."<<endl;
	double NormMin= 0;
	double NormMax= 1;
	//double NormMin= 1;
	//double NormMax= 256;
	for(int i=0;i<nRegions;i++){
		for(int j=0;j<nPars;j++){
			double parValue= params(i,j); 
			double parValue_min= params_min(j);
			double parValue_max= params_max(j);
			double parValue_norm= NormMin + (NormMax-NormMin)*(parValue-parValue_min)/(parValue_max-parValue_min);
			params(i,j)= parValue;
			//params(i,j)= parValue_norm;
		}//end loop pars
		for(int j=0;j<nPars_robust;j++){
			double parValue= params_robust(i,j); 
			double parValue_min= params_robust_min(j);
			double parValue_max= params_robust_max(j);
			double parValue_norm= NormMin + (NormMax-NormMin)*(parValue-parValue_min)/(parValue_max-parValue_min);
			params_robust(i,j)= parValue;
			//params_robust(i,j)= parValue_norm;
		}//end loop pars
		for(int j=0;j<nPars_spatial;j++){
			double parValue= params_spatial(i,j); 
			double parValue_min= params_spatial_min(j);
			double parValue_max= params_spatial_max(j);
			double parValue_norm= NormMin + (NormMax-NormMin)*(parValue-parValue_min)/(parValue_max-parValue_min);
			//params_spatial(i,j)= parValue;
			params_spatial(i,j)= parValue_norm;
		}//end loop pars
	}//end loop regions

	//## DEBUG!!!
	//params_robust.Print();

	//## Compute mutual region distances
	TMatrixD* DissimilarityMatrix= new TMatrixD(nRegions,nRegions);
	DissimilarityMatrix->Zero();
	TMatrixD* DissimilarityMatrix_robust= new TMatrixD(nRegions,nRegions);
	DissimilarityMatrix_robust->Zero();

	for(int i=0;i<nRegions-1;i++){
		TVectorD x_i = TMatrixDRow(params,i);
		TVectorD xrobust_i = TMatrixDRow(params_robust,i);
		TVectorD xspatial_i = TMatrixDRow(params_spatial,i);
		for(int j=i+1;j<nRegions;j++){
			TVectorD x_j = TMatrixDRow(params,j);
			TVectorD xrobust_j = TMatrixDRow(params_robust,j);
			TVectorD xspatial_j = TMatrixDRow(params_spatial,j);

			double dist2= (x_i-x_j).Norm2Sqr();
			double dist2_robust= (xrobust_i-xrobust_j).Norm2Sqr();
			double dist2_spatial= (xspatial_i-xspatial_j).Norm2Sqr();

			double dist= sqrt(dist2);
			double dist_robust= sqrt(dist2_robust);
			double dist_spatial= sqrt(dist2_spatial);

			//double dissimilarity= sqrt(dist2)/(1+sqrt(dist2_spatial));
			//double dissimilarity_robust= sqrt(dist2_robust)/(1+sqrt(dist2_spatial));

			double dissimilarity= dist * exp(-spatialRegFactor*dist_spatial);
			double dissimilarity_robust= dist_robust * exp(-spatialRegFactor*dist_spatial);

			DissimilarityMatrix->operator()(i,j)= dissimilarity;
			DissimilarityMatrix->operator()(j,i)= dissimilarity;
			DissimilarityMatrix_robust->operator()(i,j)= dissimilarity_robust;
			DissimilarityMatrix_robust->operator()(j,i)= dissimilarity_robust;
		}//end loop regions
	}//end loop regions

	
	TString imgName= Form("%s_saliency",img->GetName());
	Img* saliencyImg= (Img*)img->GetCloned(std::string(imgName),true,true);
	saliencyImg->SetNameTitle(imgName,imgName);
	saliencyImg->Reset();

	//## Compute saliency (METHOD1)
	TGraph2D* interpolationGraph= new TGraph2D;

	for(int i=0;i<nRegions;i++){
		
		//Sort dissimilarities
		std::vector<double> dissList;
		std::vector<double> dissList_robust;

		for(int j=0;j<nRegions;j++){
			double D= DissimilarityMatrix->operator()(i,j);
			double D_robust= DissimilarityMatrix_robust->operator()(i,j);
			dissList.push_back(D);
			dissList_robust.push_back(D_robust);
		}//end loop regions

		std::vector<double> sorted;
		std::vector<size_t> sort_index;//sorting index
		CodeUtils::sort( dissList,sorted,sort_index);	

		std::vector<double> sorted_robust;
		std::vector<size_t> sort_index_robust;//sorting index
		CodeUtils::sort( dissList_robust,sorted_robust,sort_index_robust);	
		
		//Compute saliency over k-th neighbors
		double sum= 0;
		double sum_robust= 0;
		for(int k=0;k<KNN;k++){
			size_t index= sort_index[k];
			double diss= dissList[index];
			size_t index_robust= sort_index_robust[k];
			double diss_robust= dissList_robust[index_robust];
			sum+= diss;
			sum_robust+= diss_robust;
		}//end loop 
		double Saliency= 0;
		if(useRobust) Saliency= 1.-exp(-sum_robust/(double)(KNN));
		else Saliency= 1.-exp(-sum/(double)(KNN));
				
		//Fill interpolation graph
		double Cx= regions[i]->X0;	
		double Cy= regions[i]->Y0;	
		interpolationGraph->SetPoint(i,Cx,Cy,Saliency);

		//Fill image
		//PixelCollection pixels= regions[i]->GetPixels();
		//for(unsigned int j=0;j<pixels.size();j++){//loop on pixels inside region
		for(int j=0;j<regions[i]->GetNPixels();j++){//loop on pixels inside region
			Pixel* thisPixel= regions[i]->GetPixel(j);
			double thisPixelX= thisPixel->x;
			double thisPixelY= thisPixel->y;
			//double thisPixelX= pixels[j]->x;
			//double thisPixelY= pixels[j]->y;
			saliencyImg->FillPixel(thisPixelX,thisPixelY,Saliency);
		}//end loop pixels in region
	}//end loop regions
	
	//Clear matrix
	if(DissimilarityMatrix) DissimilarityMatrix->Delete();
	if(DissimilarityMatrix_robust) DissimilarityMatrix_robust->Delete();

	
	//## Interpolate?
	//if(interpolate){
		cout<<"Img::GetSaliencyMap(): INFO: Interpolating saliency map..."<<endl;
		for(int i=0;i<saliencyImg->GetNbinsX();i++){	
			double x= saliencyImg->GetXaxis()->GetBinCenter(i+1);				
			for(int j=0;j<saliencyImg->GetNbinsY();j++){
				double y= saliencyImg->GetYaxis()->GetBinCenter(j+1);
				double interpSaliency= interpolationGraph->Interpolate(x,y);
				if(interpSaliency!=0){
					saliencyImg->SetBinContent(i+1,j+1,interpSaliency);
				}
			}//end loop bins Y
		}//end loop binsX
	//}//close if
	if(interpolationGraph) interpolationGraph->Delete();
	

	return saliencyImg;

}//close ComputeSaliencyMap()


Img* SaliencyFilter::ComputeMultiResoSaliencyMap(Img* img,int resoMin,int resoMax,int resoStep,double beta,int minRegionSize,double knnFactor,double spatialRegFactor,bool useRobust,bool addCurvDist,double salientMultiplicityThrFactor,bool addBkgMap,bool addNoiseMap,BkgData* bkgData,double saliencyThrFactor,double imgThrFactor){
	
	//## Check input img
	if(!img){
		cerr<<"SaliencyFilter::ComputeMultiResoSaliencyMap(): ERROR: Null ptr to given input image!"<<endl;
		return 0;
	}
	
	//## Check bkg data	
	if(!bkgData && (addBkgMap || addNoiseMap)){
		cerr<<"SaliencyFilter::ComputeMultiResoSaliencyMap(): ERROR: Selected to use bkgdata in saliency computation but no bkg data given!"<<endl;
		return 0;
	}

	//## Compute img stats
	if(!img->HasStats()){
		img->ComputeStats(true,false,false);
	}
	ImgStats* imgStats= img->GetPixelStats();
	if(!imgStats){
		cerr<<"SaliencyFilter::ComputeMultiResoSaliencyMap(): ERROR: No stats available for this image!"<<endl;
		return 0;
	}
	double imgMedian= imgStats->median;
	double imgMedianThr= imgThrFactor*imgMedian;
	int nReso= (resoMax-resoMin)/resoStep + 1;
	
	TString imgName= Form("%s_saliencyMean",img->GetName());
	Img* saliencyImg_mean= (Img*)img->GetCloned(std::string(imgName),true,true);
	saliencyImg_mean->SetNameTitle(imgName,imgName);
	saliencyImg_mean->Reset();

	//double NormMin= 0;
	//double NormMax= 1;
	double NormMin= 1;
	double NormMax= 256;
	std::vector<Img*> salMaps;
	std::vector<double> salMapsThresholds;
	int nbins= 100;
	
	for(int i=0;i<nReso;i++){
		int reso= resoMin + i*resoStep;
		cout<<"SaliencyFilter::ComputeMultiResoSaliencyMap(): INFO: Computing saliency map @ reso "<<reso<<" (step="<<resoStep<<")"<<endl;

		//Compute saliency map @ current reso and normalize
		Img* salMap= ComputeSaliencyMap(img,reso,beta,minRegionSize,knnFactor,spatialRegFactor,useRobust,addCurvDist);
		Img* salMap_norm= salMap->GetNormalizedImage("LINEAR",NormMin,NormMax);
		saliencyImg_mean->Add(salMap_norm);		
		if(salMap) salMap->Delete();
		salMaps.push_back(salMap_norm);

		//Compute stats		
		salMap_norm->ComputeStats(true,false,false);
		double salMedian= (salMap_norm->GetPixelStats())->median;
		double medianThr= saliencyThrFactor*salMedian;
		double otsuThr= salMap_norm->FindOtsuThreshold(nbins);
		double valleyThr= salMap_norm->FindValleyThreshold(nbins,true);
		//double salThr= std::max(medianThr,otsuThr);
		//double salThr= std::max(std::min(otsuThr,valleyThr),medianThr);
		double salThr= medianThr;
		salMapsThresholds.push_back(salThr);
	}//end loop reso
	
	//Normalize final saliency
	saliencyImg_mean->Scale(1./(double)nReso);
	double minSaliency= saliencyImg_mean->GetMinimum();
	
	cout<<"SaliencyFilter::ComputeMultiResoSaliencyMap(): INFO: Normalize saliency sum over reso..."<<endl;
	imgName= Form("%s_saliencyCombined",img->GetName());
	Img* saliencyImg= (Img*)saliencyImg_mean->GetCloned(std::string(imgName),true,true);
	saliencyImg->SetNameTitle(imgName,imgName);
	saliencyImg->Reset();
	saliencyImg_mean->Delete();
	
	//## Normalize saliency (using adaptive threshold)
	int salientMultiplicityThr= (int)std::round(salientMultiplicityThrFactor*nReso);

	for(int i=0;i<saliencyImg->GetNbinsX();i++){
		for(int j=0;j<saliencyImg->GetNbinsY();j++){
			double imgBinContent= img->GetBinContent(i+1,j+1);	
			if(imgBinContent==0) continue;
			bool isPositiveExcess= (imgBinContent>imgMedianThr);
			double wmin= 1.e+99;
			double wmax= -1.e+99;
			int saliencyMultiplicity= 0;
			for(unsigned int k=0;k<salMaps.size();k++){
				double thisw= salMaps[k]->GetBinContent(i+1,j+1);
				double thisThreshold= salMapsThresholds[k];
				if(thisw>thisThreshold) saliencyMultiplicity++;
				if(thisw<wmin) wmin= thisw;
				if(thisw>=wmax) wmax= thisw;
			}//end loop multi reso
			
			//cout<<"bin ("<<i<<","<<j<<") saliencyMultiplicity="<<saliencyMultiplicity<<" (salientMultiplicityThr="<<salientMultiplicityThr<<")"<<endl;
			if(saliencyMultiplicity>=salientMultiplicityThr){
				if(isPositiveExcess) saliencyImg->SetBinContent(i+1,j+1,wmax);
				else saliencyImg->SetBinContent(i+1,j+1,minSaliency);
			}
			else {
				saliencyImg->SetBinContent(i+1,j+1,wmin);
			}
		}//end loop bins Y
	}//end loop bins X
	

	//Clear map list
	for(unsigned int k=0;k<salMaps.size();k++){
		if(salMaps[k]) salMaps[k]->Delete();		
	}
	salMaps.clear();
	
	//Normalize bkg and noise maps
	if(addBkgMap){
		cout<<"SaliencyFilter::ComputeMultiResoSaliencyMap(): INFO: Normalize bkg map..."<<endl;
		Img* bkgImg= (bkgData->BkgMap)->GetNormalizedImage("LINEAR",NormMin,NormMax);
		if(bkgImg){
			saliencyImg->Add(bkgImg);	
			bkgImg->Delete();
		}
		else{
			cerr<<"SaliencyFilter::ComputeMultiResoSaliencyMap(): WARN: Failed to normalize bkg map, cannot add it to saliency computation!"<<endl;
		}
	}//close if

	if(addNoiseMap){
		Img* noiseImg= (bkgData->NoiseMap)->GetNormalizedImage("LINEAR",NormMin,NormMax);
		if(noiseImg){
			saliencyImg->Add(noiseImg);
			noiseImg->Delete();	
		}
		else{
			cerr<<"SaliencyFilter::ComputeMultiResoSaliencyMap(): WARN: Failed to normalize noise map, cannot add it to saliency computation!"<<endl;			
		}
	}//close if


	//Normalize final map
	cout<<"SaliencyFilter::ComputeMultiResoSaliencyMap(): INFO: Normalize final maps..."<<endl;
	imgName= Form("%s_saliencyMultiReso",img->GetName());
	Img* saliencyMap= saliencyImg->GetNormalizedImage("LINEAR",NormMin,NormMax);
	saliencyMap->SetNameTitle(imgName,imgName);
	if(saliencyImg) saliencyImg->Delete();

	for(int i=0;i<saliencyMap->GetNbinsX();i++){
		for(int j=0;j<saliencyMap->GetNbinsY();j++){
			double imgBinContent= img->GetBinContent(i+1,j+1);
			if(imgBinContent==0){
				saliencyMap->SetBinContent(i+1,j+1,0);
			}
		}
	}		
	saliencyMap->ComputeStats(true,false,true);

	return saliencyMap;

}//close GetMultiResoSaliencyMap()


}//close namespace
