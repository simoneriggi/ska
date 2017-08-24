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
* @file SLICUtils.cc
* @class SLICData
* @brief SLIC utils class
*
* Superpixel utility functions
* @author S. Riggi
* @date 20/01/2015
*/

#include <SLICUtils.h>
#include <Region.h>
#include <Contour.h>
#include <SLICData.h>
#include <StatsUtils.h>
#include <CodeUtils.h>
#include <Image.h>


#include <TObject.h>
#include <TMatrixD.h>
#include <TVector2.h>

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

ClassImp(Caesar::SLICUtils)

namespace Caesar {

SLICUtils::SLICUtils() {

	
}//close costructor

SLICUtils::~SLICUtils(){
	

}//close destructor

Image* SLICUtils::GetSegmentedImage(Image* image,std::vector<Region*>const& regions,int selectedTag,bool normalize,bool binarize){

	//## Check input
	if(!image) {
		ERROR_LOG("Null ptr to input image given!");
		return 0;
	}

	int nRegions= (int)regions.size();
	if(nRegions<=0) {
		ERROR_LOG("No regions available, nothing to be done!");
		return 0;
	}

	//## Use image normalization?
	double A= 0;
	double B= 1;
	if(image->HasStats() && normalize){
		ImgStats* stats= image->GetPixelStats();
		double NormMin= image->GetMinimum();
		double NormMax= image->GetMaximum();
		A= NormMin - (NormMax-NormMin)*stats->min/(stats->max-stats->min);
		B= (NormMax-NormMin)/(stats->max-stats->min);
	}	
	
	//## Fill image with region means
	TString imgName= Form("%s_segmented",image->GetName().c_str());
	Image* segmentedImg= (Image*)image->GetCloned(std::string(imgName),true,true);
	segmentedImg->Reset();
	
  for(int i=0;i<nRegions;i++){//loop on regions
		int regionTag= regions[i]->Tag;
		int nPixelsInRegion= (int)regions[i]->NPix;
		DEBUG_LOG("Region no. "<<i<<" N="<<nPixelsInRegion<<" regionTag="<<regionTag);
		
		if(selectedTag!=-1 && regionTag!=selectedTag) continue;
		
		double Mean= regions[i]->Mean;
		double colorValue= A+B*Mean;
		if(binarize) colorValue= 1;

		for(int j=0;j<nPixelsInRegion;j++){//loop on pixels inside region
			double thisPixelX= (regions[i]->GetPixel(j))->x;
			double thisPixelY= (regions[i]->GetPixel(j))->y;
			segmentedImg->Fill(thisPixelX,thisPixelY,colorValue);
		}//end loop pixels in region	
	}//end loop regions
	
	return segmentedImg;

}//close GetSegmentedImage()


SLICContourData* SLICUtils::ComputeBoundaryContours(SLICData* slicData) {
  	
	//## Check input
	if(!slicData || !slicData->inputImg) return 0;
	int nRegions= slicData->GetNRegions();
	if(nRegions<=0) return 0;

	//## Init data
	cout<<"SLICUtils::ComputeBoundaryContours(): INFO: Computing contours from NR="<<nRegions<<" regions..."<<endl;
	long int Nx= (slicData->inputImg)->GetNx();
	long int Ny= (slicData->inputImg)->GetNy();
	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};
		
	std::map<int,int> mapping;
	std::vector< std::vector<bool> > isNeighborLinkTaken;
	SLICContourData* contourData= new SLICContourData;
	(contourData->contour)= new Contour;
	
	for(int k=0;k<nRegions;k++){
		int regionId= (slicData->regions)[k]->Id;
		mapping.insert( std::pair<int,int>(regionId,k) );
		contourData->connectedRegionIds.push_back( std::vector<int>() );
		contourData->boundaryData.push_back( SLICBoundaryPixMap() );
		isNeighborLinkTaken.push_back( std::vector<bool>() );
		for(int s=0;s<nRegions;s++) isNeighborLinkTaken[k].push_back(false);
	}//end loop regions 

  //## Loop over all the pixels
  for (long int i=0; i<Nx; i++) {
  	for (long int j=0; j<Ny; j++) {
			//int regionId= (slicData->pixelLabels)->operator()(i,j);
			int regionId= (slicData->labels)[i][j];
			if(regionId<0) continue;
			int regionIndex= mapping[regionId];			
   		int nr_p = 0;
			std::map<int,int> connectedRegionCounterMap;
			connectedRegionCounterMap.clear();
			std::map<int,std::vector<int>> connectedRegionBoundaryPixelListMap;
			connectedRegionBoundaryPixelListMap.clear();
            
      // Compare the pixel to its 8 neighbours
      for (int k=0; k<8; k++) {
				int x = i + dx8[k];
				int y = j + dy8[k];
				         
        if (x >= 0 && x < Nx && y >= 0 && y < Ny) {
					//int regionId_neighbor= (slicData->pixelLabels)->operator()(x,y);
       		int regionId_neighbor= (slicData->labels)[x][j];
       	
        	if (regionId!=regionId_neighbor && regionId_neighbor>=0) {	
          	nr_p++;
						connectedRegionCounterMap[regionId_neighbor]++;
						long int gBin_neighbor= (slicData->inputImg)->GetBin(x,y);
						connectedRegionBoundaryPixelListMap[regionId_neighbor];
						connectedRegionBoundaryPixelListMap[regionId_neighbor].push_back(gBin_neighbor);
          }
        }//close if
      }//end loop neighbours
            
      // Add the pixel to the contour list of corresponding region
     	if (nr_p>= 2) {
				long int gBin= (slicData->inputImg)->GetBin(i,j);
				double binX= (slicData->inputImg)->GetX(i);
				double binY= (slicData->inputImg)->GetY(j);	
				(contourData->contour)->AddPoint(TVector2(binX,binY));
      	
				//Add neighbor ids and boundary pixel ids
				std::map<int,std::vector<int>>::iterator counterListIt = connectedRegionBoundaryPixelListMap.begin();
				for (counterListIt=connectedRegionBoundaryPixelListMap.begin(); counterListIt!=connectedRegionBoundaryPixelListMap.end(); ++counterListIt){
					int neighborId= counterListIt->first;
					std::vector<int> neighborPixIds= counterListIt->second;
					int counts= (int)neighborPixIds.size();	
					int neighborIndex= mapping[neighborId];

					if(counts>=2) {		
						std::vector<int>::iterator it = std::find (contourData->connectedRegionIds[regionIndex].begin(), contourData->connectedRegionIds[regionIndex].end(), neighborIndex);
						
						if( contourData->connectedRegionIds[regionIndex].empty() || it==contourData->connectedRegionIds[regionIndex].end() ) {
							contourData->connectedRegionIds[regionIndex].push_back(neighborIndex);
						}//close if	

						(contourData->boundaryData[regionIndex])[neighborId];//add connection with neighbor if not existing in the map
						((contourData->boundaryData[regionIndex])[neighborId]).push_back(gBin);//add this contour pixel
						for(int t=0;t<counts;t++){
							int gBin_neighbor= neighborPixIds[t];
							it = std::find ( ((contourData->boundaryData[regionIndex])[neighborId]).begin(), ((contourData->boundaryData[regionIndex])[neighborId]).end(), gBin_neighbor);
							if( ((contourData->boundaryData[regionIndex])[neighborId]).empty() || it==((contourData->boundaryData[regionIndex])[neighborId]).end() ) {
								((contourData->boundaryData[regionIndex])[neighborId]).push_back(gBin_neighbor);
							}//close if
						}//end loop counts
						
					}//close if counts>2
				}//end loop map counter iterator
      }//close if nr_p>2
    }//end loop image Ny
  }//end loop image Nx
	
	return contourData;

}//close SLICUtils::ComputeBoundaryContours()


int SLICUtils::FindNeighbors(std::vector<SLICNeighborCollection>& neighbors,SLICData* slicData,SLICContourData* contourData,bool get2ndNeighbors,int selectedTag,bool includeSpatialDist){
	
	//## Check input
	if(!slicData || !contourData){
		cerr<<"SLICUtils::FindNeighbors(): ERROR: Null ptr to SLIC data given!"<<endl;
		return -1;
	}
	if(!contourData) return -1;
	if(!slicData->edgeImg) return -1;

	int nRegions= slicData->GetNRegions();
	if(nRegions<=0){
		cerr<<"SLICUtils::FindNeighbors(): WARN: No regions available, nothing to be done!"<<endl;
		return 0;
	}
	
	//## Fill list of 1st and 2-nd neighbors for each region
	//## Include only selected tags 	
	SLICConnectedRegions connectedRegionIds= contourData->connectedRegionIds;
	std::vector<SLICBoundaryPixMap> boundaryData= contourData->boundaryData;
	std::vector< std::vector<long int> > neighborIndexList;
	double Emax= 1;

	bool useRobustParams= false;
	bool normalizeParams= true;
	bool addSpatialDist= false;
	bool addCurvDist= true;

	for(int i=0; i<nRegions; i++) {
		int regionId= (slicData->regions)[i]->Id;			
		neighborIndexList.push_back( std::vector<long int>() );
		SLICNeighborCollection aNeighborCollection;
		neighbors.push_back(aNeighborCollection);

		//Fill 1st-order neighbors
		for(unsigned int j=0;j<connectedRegionIds[i].size();j++){//loop over 1st neighbors
			int neighborIndex= connectedRegionIds[i][j];
			int neighborTag= (slicData->regions)[neighborIndex]->Tag;
			int neighborId= (slicData->regions)[neighborIndex]->Id;

			if(selectedTag==-1 || neighborTag==selectedTag){
				//Find symmetric dissimilarity
				double Dsym, Dsym_space;
				int status= (slicData->regions)[i]->GetDistance(Dsym, Dsym_space,(slicData->regions)[neighborIndex],useRobustParams,normalizeParams,addCurvDist);
				if(status<0){
					cerr<<"SLICUtils::FindNeighbors(): WARN: Symm distance calculation failed!"<<endl;
				}

				//Find asymm dissimilarities
				double Diss, DissNeighbor;
				status= (slicData->regions)[i]->GetAsymmDistance(Diss,DissNeighbor,(slicData->regions)[neighborIndex],useRobustParams,normalizeParams,includeSpatialDist,addCurvDist);
				if(status<0){
					cerr<<"SLICUtils::FindNeighbors(): WARN: Asymm distance calculation failed!"<<endl;
				}
		
				//Find edgeness
				SLICBoundaryPixMapIterator it= boundaryData[i].find(neighborId);
				double E= 0;
				if(it!=boundaryData[i].end()){//region is among neighbors, compute edgeness terms
					std::vector<int> sharedPixelIds= (boundaryData[i])[neighborId];
					int nBoundaryPixels= (int)sharedPixelIds.size();			
					for(int t=0;t<nBoundaryPixels;t++) {
						int gBin= sharedPixelIds[t];
						double S_edge= (slicData->edgeImg)->GetBinContent(gBin);
						E+= S_edge;
					}	
					if(nBoundaryPixels>0) E/= (double)nBoundaryPixels;  
				}//close if
	
				SLICBoundaryPixMapIterator it_neighbor= boundaryData[neighborIndex].find(regionId);
				double E_neighbor= 0;

				if(it_neighbor!=boundaryData[neighborIndex].end()){//region is among neighbors, compute edgeness terms
					std::vector<int> sharedPixelIds= (boundaryData[neighborIndex])[regionId];
					int nBoundaryPixels= (int)sharedPixelIds.size();			
					for(int t=0;t<nBoundaryPixels;t++) {
						int gBin= sharedPixelIds[t];
						double S_edge= (slicData->edgeImg)->GetBinContent(gBin);
						E_neighbor+= S_edge;
					}	
					if(nBoundaryPixels>0) E_neighbor/= (double)nBoundaryPixels;  
				}//close if

				SLICNeighborData neighborData;
				neighborData.Order= 1;
				neighborData.Id= neighborId;
				neighborData.Index= neighborIndex;
				neighborData.Tag= neighborTag;
				neighborData.D= Diss;
				neighborData.D_n= DissNeighbor;
				neighborData.E= E;	
				neighborData.E_n= E_neighbor;
				neighborData.Dsym= Dsym;

				neighborIndexList[i].push_back(neighborIndex);
				neighbors[i].Add(neighborData);
			}//close if selected region tag

			//Fill 2nd-order neighbors?
			if(!get2ndNeighbors) continue;

			for(unsigned int t=0;t<connectedRegionIds[neighborIndex].size();t++){//loop over 2nd neighbors
				int neighborIndex_2nd= connectedRegionIds[neighborIndex][t];
				int neighborTag_2nd= (slicData->regions)[neighborIndex_2nd]->Tag;
				int neighborId_2nd= (slicData->regions)[neighborIndex_2nd]->Id;
			
				if( (selectedTag!=-1 || neighborTag_2nd==selectedTag) && neighborId_2nd!=regionId ){
					std::vector<long int>::iterator it= std::find(neighborIndexList[i].begin(),neighborIndexList[i].end(),neighborIndex_2nd);
						
					if( neighborIndexList[i].size()==0 || it==neighborIndexList[i].end() ){
						//Find symm diss for 2nd neighbors
						double Dsym_2nd, Dsym_space_2nd;
						int status= (slicData->regions)[i]->GetDistance(Dsym_2nd, Dsym_space_2nd,(slicData->regions)[neighborIndex_2nd],useRobustParams,normalizeParams,addCurvDist);
						if(status<0){
							cerr<<"SLICUtils::FindNeighbors(): WARN: Symm distance calculation failed!"<<endl;
						}
	
						//Find asymm diss for 2nd neighbors
						double Diss_2nd, DissNeighbor_2nd;
						status= (slicData->regions)[i]->GetAsymmDistance(Diss_2nd,DissNeighbor_2nd,(slicData->regions)[neighborIndex_2nd],useRobustParams,normalizeParams,includeSpatialDist,addCurvDist);
						if(status<0){
							cerr<<"SLICUtils::FindNeighbors(): WARN: Asymm distance calculation failed!"<<endl;
						}
						SLICNeighborData neighborData_2nd;
						neighborData_2nd.Order= 2;
						neighborData_2nd.Id= neighborId_2nd;
						neighborData_2nd.Index= neighborIndex_2nd;
						neighborData_2nd.Tag= neighborTag_2nd;
						neighborData_2nd.D= Diss_2nd;
						neighborData_2nd.D_n= DissNeighbor_2nd;
						neighborData_2nd.E= Emax;	
						neighborData_2nd.E_n= Emax;
						neighborData_2nd.Dsym= Dsym_2nd;
						neighbors[i].Add(neighborData_2nd);

						neighborIndexList[i].push_back(neighborIndex_2nd);
					}//close if item not found in collection	
				}//close if
			}//end loop 2nd-neighbors
		}//end loop 1st-neighbors
	}//end loop regions

	neighborIndexList.clear();

	return 0;

}//close FindNeighbors()


SLICSimilarityData* SLICUtils::ComputeRegionSimilarity(SLICData* slicData,std::vector<SLICNeighborCollection>& neighbors,double beta){

	//## Check input data
	if(!slicData) return 0;
	int nRegions= slicData->GetNRegions();
	if(nRegions<=0) return 0;
	if(!slicData->edgeImg) return 0;
	ImgStats* edgeImgStats= (slicData->edgeImg)->GetPixelStats();
	if(!edgeImgStats){
		cerr<<"SLICUtils::ComputeRegionSimilarity(): ERROR: No stats for edge image...return!"<<endl;
		return 0;
	}
	double EdgeImgMin= edgeImgStats->min;
	double EdgeImgMax= edgeImgStats->max;
	if(EdgeImgMin==EdgeImgMax) EdgeImgMin= 0;
	double Emin_norm= 0;
	double Emax_norm= 1;
	double Dmin_norm= 0;
	double Dmax_norm= 1;
	cout<<"SLICUtils::ComputeRegionSimilarity(): INFO: Compute region similarities (nRegions="<<nRegions<<")"<<endl;

	//## Init matrix
	//TMatrixD* DissimilarityMatrix= new TMatrixD(nRegions,nRegions);
	//DissimilarityMatrix->Zero();

	TMatrixD* AdjacencyMatrix= new TMatrixD(nRegions,nRegions);
	AdjacencyMatrix->Zero();

	//TMatrixD* EdgenessMatrix= new TMatrixD(nRegions,nRegions);
	//EdgenessMatrix->Zero();
	//(*EdgenessMatrix)+= Emax_norm;

	//## Fill matrix
	double Dmin= 1.e+99;
	double Dmax= -1.e+99;
	double Emin= 1.e+99;
	double Emax= -1.e+99;
	/*
	for(unsigned int i=0;i<neighbors.size();i++){
		for(unsigned int j=0;j<neighbors[i].size();j++){
			int neighborIndex= neighbors[i][j].Index;
			int neighborOrder= neighbors[i][j].Order;
			double D= neighbors[i][j].D;
			double E= neighbors[i][j].E;

			(*DissimilarityMatrix)(i,j)= D;
			(*EdgenessMatrix)(i,j)= E;
			if(D>Dmax) Dmax= D;
			if(D<Dmin) Dmin= D;
			if(E>Emax) Emax= E;
			if(E<Emin) Emin= E;
		}//end loop neighbor regions
	}//end loop regions
	*/

	std::vector<double> Dlist;

	for(unsigned int i=0;i<neighbors.size();i++){
		std::vector<SLICNeighborData> thisNeighbors= neighbors[i].GetNeighbors();
		for(unsigned int j=0;j<thisNeighbors.size();j++){//loop over neighbor list for this region
			
			int neighborIndex= thisNeighbors[j].Index;
			double D= thisNeighbors[j].D;
			double E= thisNeighbors[j].E;

			//(*DissimilarityMatrix)(i,neighborIndex)= D;
			//(*EdgenessMatrix)(i,neighborIndex)= E;
			if(D>Dmax) Dmax= D;
			if(D<Dmin) Dmin= D;
			if(E>Emax) Emax= E;
			if(E<Emin) Emin= E;
			Dlist.push_back(D);
		}//end loop neighbor regions
	}//end loop regions

	double Dmedian= StatsUtils::GetMedianFast<double>(Dlist);
	double Dmad= StatsUtils::GetMADFast(Dlist,Dmedian);
	double Dmedianrms= Dmad*1.4826;
	

	//## Normalize matrix
	const double SMALL_NUMBER= 0.00000000001;
	
	/*
	for(unsigned int i=0;i<neighbors.size();i++){
		for(unsigned int j=0;j<neighbors[i].size();j++){
			int neighborIndex= neighbors[i][j].Index;
			int neighborOrder= neighbors[i][j].Order;
			double D= neighbors[i][j].D;
			double E= neighbors[i][j].E;

			double D_norm= Dmin_norm + (Dmax_norm-Dmin_norm)*(D-Dmin)/(Dmax-Dmin);
			double E_norm= Emin_norm + (Emax_norm-Emin_norm)*(E-Emin)/(Emax-Emin);
			double Dtot= (1-beta)*D_norm + beta*E_norm;
			//double Dtot= D_norm + beta*E_norm;
			Dtot+= SMALL_NUMBER;//to avoid dividing by zero!

			neighbors[i][j].Dtot= Dtot;
			(*DissimilarityMatrix)(i,neighborIndex)= Dtot;
			(*EdgenessMatrix)(i,neighborIndex)= E_norm;
			if(neighborOrder>0) (*AdjacencyMatrix)(i,neighborIndex)= 1./Dtot;	
		}//end loop neighbor regions
	}//end loop regions
	*/

	for(unsigned int i=0;i<neighbors.size();i++){
		std::vector<SLICNeighborData> thisNeighbors= neighbors[i].GetNeighbors();
		for(unsigned int j=0;j<thisNeighbors.size();j++){//loop over neighbor list for this region
			int neighborIndex= thisNeighbors[j].Index;
			int neighborOrder= thisNeighbors[j].Order;
			double D= thisNeighbors[j].D;
			double D_n= thisNeighbors[j].D_n;
			double E= thisNeighbors[j].E;
			double E_n= thisNeighbors[j].E_n;

			double D_norm= Dmin_norm + (Dmax_norm-Dmin_norm)*(D-Dmin)/(Dmax-Dmin);
			double D_n_norm= Dmin_norm + (Dmax_norm-Dmin_norm)*(D_n-Dmin)/(Dmax-Dmin);
			
			double E_norm= Emin_norm + (Emax_norm-Emin_norm)*(E-Emin)/(Emax-Emin);
			double E_n_norm= Emin_norm + (Emax_norm-Emin_norm)*(E_n-Emin)/(Emax-Emin);
			double Dtot= (1-beta)*D_norm + beta*E_norm;
			double Dtot_n= (1-beta)*D_n_norm + beta*E_n_norm;
			//double Dtot= D_norm + beta*E_norm;
			Dtot+= SMALL_NUMBER;//to avoid dividing by zero!
			Dtot_n+= SMALL_NUMBER;//to avoid dividing by zero!

			(neighbors[i].GetNeighbor(j))->Dtot= Dtot;//Set Dtot in current neighbor
			(neighbors[i].GetNeighbor(j))->Dtot_n= Dtot_n;//Set Dtot in current neighbor

			//(*DissimilarityMatrix)(i,neighborIndex)= Dtot;
			if(neighborOrder>0) (*AdjacencyMatrix)(i,neighborIndex)= 1./Dtot;	
		}//end loop neighbor regions
	}//end loop regions
		
	//## Normalize adjacency matrix by rows
	for(int i=0;i<AdjacencyMatrix->GetNrows();i++){
		double sum= 0;	
		for(int j=0;j<AdjacencyMatrix->GetNcols();j++) {
			sum+= (*AdjacencyMatrix)(i,j);
		}
		if(sum!=0) {
			for(int j=0;j<AdjacencyMatrix->GetNcols();j++) (*AdjacencyMatrix)(i,j)/= sum;
		}
	}//end loop rows		
			
				
	//## Return data
	SLICSimilarityData* SimilarityData= new SLICSimilarityData;
	//SimilarityData->DissimilarityMatrix= DissimilarityMatrix;
	SimilarityData->AdjacencyMatrix= AdjacencyMatrix;  
	SimilarityData->Dmin= Dmin;
	SimilarityData->Dmax= Dmax;
	SimilarityData->Dmedian= Dmedian;
	SimilarityData->Dmedianrms= Dmedianrms;
	SimilarityData->Emin= Emin;
	SimilarityData->Emax= Emax;
	
	//EdgenessMatrix->Delete();

	return SimilarityData;

}//close ComputeRegionSimilarity()




/*
SLICSimilarityData* SLICUtils::ComputeRegionSimilarity(SLICData* slicData,SLICContourData* contourData,double beta,bool includeSpatialDist,int mergedTag){

	//## Check input data
	if(!slicData) return 0;
	int nRegions= slicData->GetNRegions();
	if(nRegions<=0) return 0;
	if(!contourData) return 0;
	if(!slicData->edgeImg) return 0;
	ImgStats* edgeImgStats= (slicData->edgeImg)->GetPixelStats();
	if(!edgeImgStats){
		cerr<<"SLICUtils::ComputeRegionSimilarity(): ERROR: No stats for edge image...return!"<<endl;
		return 0;
	}
	double EdgeImgMin= edgeImgStats->min;
	double EdgeImgMax= edgeImgStats->max;
	if(EdgeImgMin==EdgeImgMax) EdgeImgMin= 0;
	double Emin_norm= 0;
	double Emax_norm= 1;
	double Dmin_norm= 0;
	double Dmax_norm= 1;
	cout<<"SLICUtils::ComputeRegionSimilarity(): INFO: Compute region similarities (nRegions="<<nRegions<<")"<<endl;

	//## Init data
	SLICSimilarityData* SimilarityData= new SLICSimilarityData;
	(SimilarityData->DissimilarityMatrix)= new TMatrixD(nRegions,nRegions);
	(SimilarityData->DissimilarityMatrix)->Zero();

	(SimilarityData->AdjacencyMatrix)= new TMatrixD(nRegions,nRegions);
	(SimilarityData->AdjacencyMatrix)->Zero();

	(SimilarityData->AbsDissimilarityMatrix)= new TMatrixD(nRegions,nRegions);
	(SimilarityData->AbsDissimilarityMatrix)->Zero();

	(SimilarityData->NeighborMatrix)= new TMatrixD(nRegions,nRegions);
	(SimilarityData->NeighborMatrix)->Zero();
	(*SimilarityData->NeighborMatrix)+= 999;
	for(int i=0;i<nRegions;i++) (*(SimilarityData->NeighborMatrix))(i,i)= 0;
	
	TMatrixD* EdgenessMatrix= new TMatrixD(nRegions,nRegions);
	EdgenessMatrix->Zero();
	
	std::vector<SLICBoundaryPixMap> boundaryData= contourData->boundaryData;
	SLICConnectedRegions connectedRegionIds= contourData->connectedRegionIds;

	std::vector< std::vector<double> > DissMatrix;
	std::vector<double> AbsDissList;
	for(int i=0;i<nRegions;i++){
		DissMatrix.push_back( std::vector<double>() );
		(SimilarityData->DissimilaritySortIndexMatrix).push_back( std::vector<int>() );
		for(int j=0;j<nRegions;j++) {
			DissMatrix[i].push_back(1.e+99);
			EdgenessMatrix->operator()(i,j)= Emax_norm;
		}
	}

	
	//## Consider only 1-st and 2-nd neighbors
	std::vector< std::vector<int> > neighborIndexList_2nd;
	std::vector<double> EdgenessList;
	std::vector<double> DissList;

	for(int i=0; i<nRegions; i++) {
		int regionId= (slicData->regions)[i]->Id;

		//Fill list of 2-nd neighbors	
		neighborIndexList_2nd.push_back( std::vector<int>() );		
		
		for(unsigned int j=0;j<connectedRegionIds[i].size();j++){//loop over 1st neighbors
			int neighborIndex= connectedRegionIds[i][j];

			for(unsigned int t=0;t<connectedRegionIds[neighborIndex].size();t++){//loop over 2nd neighbors
				int neighborIndex_2nd= connectedRegionIds[neighborIndex][t];
				//int neighborTag_2nd= (slicData->regions)[neighborIndex_2nd]->Tag;
				//if(mergedTag!=-1 && neighborTag_2nd!=mergedTag) continue;//skip neighbors if tag is different from the requested

				std::vector<int>::iterator it= std::find(connectedRegionIds[i].begin(),connectedRegionIds[i].end(),neighborIndex_2nd);
				std::vector<int>::iterator it2= std::find(neighborIndexList_2nd[i].begin(),neighborIndexList_2nd[i].end(),neighborIndex_2nd);
		
				if( it==connectedRegionIds[i].end() && (it2==neighborIndexList_2nd[i].end() || neighborIndexList_2nd[i].size()==0) ){
					neighborIndexList_2nd[i].push_back(neighborIndex_2nd);
				}
			}//end loop 2nd-neighbors
		}//end loop 1st-neighbors

		//Loop over 1st-neighbors
		for(unsigned int j=0;j<connectedRegionIds[i].size();j++){
			int neighborIndex= connectedRegionIds[i][j];
			int neighborId= (slicData->regions)[neighborIndex]->Id;
			int neighborTag= (slicData->regions)[neighborIndex]->Tag;	

			//Set neighbor id
			(*(SimilarityData->NeighborMatrix))(i,neighborIndex)= 1;
			(*(SimilarityData->NeighborMatrix))(neighborIndex,i)= 1;
			
			if(mergedTag!=-1 && neighborTag!=mergedTag) continue;//skip neighbors if tag is different from the requested

			//Compute dissimilarity
			double Diss, DissNeighbor;
			int status= (slicData->regions)[i]->GetAsymmDistance(Diss,DissNeighbor,(slicData->regions)[neighborIndex],false,true,includeSpatialDist);	
			if(status<0){
				cerr<<"SLICUtils::ComputeRegionSimilarity(): WARN: Asymm distance calculation failed!"<<endl;
			}
			(*(SimilarityData->DissimilarityMatrix))(i,neighborIndex)= Diss;
			(*(SimilarityData->DissimilarityMatrix))(neighborIndex,i)= DissNeighbor;

			DissList.push_back(Diss);
			DissList.push_back(DissNeighbor);

			//Compute edge term
			//Find if region j-th is among neighbors of i-th
			SLICBoundaryPixMapIterator it= boundaryData[i].find(neighborId);
			
			if(it!=boundaryData[i].end()){//region is among neighbors, compute edgeness terms
				std::vector<int> sharedPixelIds= (boundaryData[i])[neighborId];
				int nBoundaryPixels= (int)sharedPixelIds.size();
				double E= 0;
				for(int t=0;t<nBoundaryPixels;t++) {
					int gBin= sharedPixelIds[t];
					double S_edge= (slicData->edgeImg)->GetBinContent(gBin);
					//double S_edge_norm= Emin_norm + (Emax_norm-Emin_norm)*(S_edge-EdgeImgMin)/(EdgeImgMax-EdgeImgMin);
					//E+= S_edge_norm;
					E+= S_edge;
				}	
				if(nBoundaryPixels>0) E/= (double)nBoundaryPixels; 
						
				std::vector<int> sharedPixelIds_neighbors= (boundaryData[neighborIndex])[regionId];
				int nBoundaryPixels_neighbors= (int)sharedPixelIds_neighbors.size();
				double E_neighbor= 0;
				for(int t=0;t<nBoundaryPixels_neighbors;t++) {
					int gBin= sharedPixelIds_neighbors[t];
					double S_edge= (slicData->edgeImg)->GetBinContent(gBin);
					//double S_edge_norm= Emin_norm + (Emax_norm-Emin_norm)*(S_edge-EdgeImgMin)/(EdgeImgMax-EdgeImgMin);
					//E_neighbor+= S_edge_norm;
					E_neighbor+= S_edge;
				}	
				if(nBoundaryPixels_neighbors>0) E_neighbor/= (double)nBoundaryPixels_neighbors; 

				(*EdgenessMatrix)(i,neighborIndex)= E;
				(*EdgenessMatrix)(neighborIndex,i)= E_neighbor;
				EdgenessList.push_back(E);
				EdgenessList.push_back(E_neighbor);
			}//close if
		}//end loop neighbors list

		//Loop over 2nd-neighbors
		for(unsigned int j=0;j<neighborIndexList_2nd[i].size();j++){
			int neighborIndex_2nd= neighborIndexList_2nd[i][j];
			//int neighborId_2nd= (slicData->regions)[neighborIndex_2nd]->Id;
			int neighborTag_2nd= (slicData->regions)[neighborIndex_2nd]->Tag;	

			//Set neighbor id
			(*(SimilarityData->NeighborMatrix))(i,neighborIndex_2nd)= 2;
			(*(SimilarityData->NeighborMatrix))(neighborIndex_2nd,i)= 2;
			
			if(mergedTag!=-1 && neighborTag_2nd!=mergedTag) continue;//skip neighbors if tag is different from the requested

			//Compute dissimilarity
			double Diss, DissNeighbor;
			int status= (slicData->regions)[i]->GetAsymmDistance(Diss,DissNeighbor,(slicData->regions)[neighborIndex_2nd],false,true,includeSpatialDist);
			if(status<0){
				cerr<<"SLICUtils::ComputeRegionSimilarity(): WARN: Asymm distance calculation failed!"<<endl;
			}
			(*(SimilarityData->DissimilarityMatrix))(i,neighborIndex_2nd)= Diss;
			(*(SimilarityData->DissimilarityMatrix))(neighborIndex_2nd,i)= DissNeighbor;
			DissList.push_back(Diss);
			DissList.push_back(DissNeighbor);

			(*EdgenessMatrix)(i,neighborIndex_2nd)= Emax_norm;
			(*EdgenessMatrix)(neighborIndex_2nd,i)= Emax_norm;
		}//end loop 2nd neighbors

	}//end loop regions
	

	//## Normalize dissimilarity matrix in [0,1]	
	std::sort(DissList.begin(),DissList.end());
	double Dmin= 0;
	double Dmax= DissList[DissList.size()-1];
	std::sort(EdgenessList.begin(),EdgenessList.end());
	double Emin= EdgeImgMin;
	double Emax= EdgeImgMax;
	cout<<"SLICUtils::ComputeRegionSimilarity(): INFO: Normalizing dissimilarity & edgeness in ["<<Dmin_norm<<","<<Dmax_norm<<"], Dmin/Dmax="<<Dmin<<"/"<<Dmax<<" Emin/Emax="<<Emin<<"/"<<Emax<<endl;

	if(Dmax<=Dmin || Emax<=Emin){
		cerr<<"SLICUtils::ComputeRegionSimilarity(): ERROR: Invalid normalization const for dissimilatiry and/or Edgeness!"<<endl;
		return 0;
	}
	
	const double SMALL_NUMBER= 0.00000000001;
	double AbsDmin= 1.e+99;
	double AbsDmax= -1.e+99;

	for(int i=0; i<nRegions; i++) {
		//int regionId= (slicData->regions)[i]->Id;

		//Fill 1-st neighbors
		//cout<<"SLICUtils::ComputeRegionSimilarity(): INFO: Region no. "<<i<<": filling 1st neighbor info..."<<endl;
		for(unsigned int j=0;j<connectedRegionIds[i].size();j++){
			int neighborIndex= connectedRegionIds[i][j];
			//int neighborId= (slicData->regions)[neighborIndex]->Id;
			int neighborTag= (slicData->regions)[neighborIndex]->Tag;
			if(mergedTag!=-1 && neighborTag!=mergedTag) continue;//skip neighbors if tag is different from the requested

			double D= (*(SimilarityData->DissimilarityMatrix))(i,neighborIndex);
			double D_norm= Dmin_norm + (Dmax_norm-Dmin_norm)*(D-Dmin)/(Dmax-Dmin);
			
			double E= (*EdgenessMatrix)(i,neighborIndex);
			double E_norm= Emin_norm + (Emax_norm-Emin_norm)*(E-Emin)/(Emax-Emin);

			double Dtot= (1-beta)*D_norm + beta*E_norm;
			//double Dtot= D_norm + beta*E_norm;			
			Dtot+= SMALL_NUMBER;//to avoid dividing by zero!

			(*(SimilarityData->DissimilarityMatrix))(i,neighborIndex)= Dtot;	
			(*(SimilarityData->AdjacencyMatrix))(i,neighborIndex)= 1./Dtot;
			(*(SimilarityData->AbsDissimilarityMatrix))(i,neighborIndex)= D;
			(*EdgenessMatrix)(i,neighborIndex)= E_norm;
			AbsDissList.push_back( D );
			DissMatrix[i][neighborIndex]= Dtot;

			if(D>AbsDmax) AbsDmax= D;
			if(D<AbsDmin) AbsDmin= D;
			//cout<<"Neighbor "<<neighborIndex<<": D="<<D<<" Dtot="<<Dtot<<endl;
		}//end loop neighbors

		//Fill 2-nd neighbors
		for(unsigned int j=0;j<neighborIndexList_2nd[i].size();j++){
			int neighborIndex_2nd= neighborIndexList_2nd[i][j];
			//int neighborId_2nd= (slicData->regions)[neighborIndex_2nd]->Id;
			int neighborTag_2nd= (slicData->regions)[neighborIndex_2nd]->Tag;			
			if(mergedTag!=-1 && neighborTag_2nd!=mergedTag) continue;//skip neighbors if tag is different from the requested

			double D= (*(SimilarityData->DissimilarityMatrix))(i,neighborIndex_2nd);
			double D_norm= Dmin_norm + (Dmax_norm-Dmin_norm)*(D-Dmin)/(Dmax-Dmin);
			
			double E_norm= Emax_norm;
			double Dtot= (1-beta)*D_norm + beta*E_norm;
			//double Dtot= D_norm + beta*E_norm;			
			Dtot+= SMALL_NUMBER;//to avoid dividing by zero!

			(*(SimilarityData->DissimilarityMatrix))(i,neighborIndex_2nd)= Dtot;	
			(*(SimilarityData->AdjacencyMatrix))(i,neighborIndex_2nd)= 1./Dtot;
			(*(SimilarityData->AbsDissimilarityMatrix))(i,neighborIndex_2nd)= D;
			(*EdgenessMatrix)(i,neighborIndex_2nd)= E_norm;
			AbsDissList.push_back( D );
			DissMatrix[i][neighborIndex_2nd]= Dtot;

			if(D>AbsDmax) AbsDmax= D;
			if(D<AbsDmin) AbsDmin= D;

			//cout<<"Neighbor "<<neighborIndex_2nd<<": D="<<D<<" Dtot="<<Dtot<<endl;
		}//end loop 2nd neighbors

	}//end loop regions
	

	//std::sort(AbsDissList.begin(),AbsDissList.end());
	//double AbsDmin= AbsDissList[0];
	//double AbsDmax= AbsDissList[AbsDissList.size()-1];
	double AbsDmedian= StatsUtils::GetMedianFast<double>(AbsDissList);
	double AbsDmad= StatsUtils::GetMADFast(AbsDissList,AbsDmedian);
	double AbsDmedianrms= AbsDmad*1.4826;
	cout<<"SLICUtils::ComputeRegionSimilarity(): INFO: AbsDissList size="<<AbsDissList.size()<<", AbsDmin/AbsDmax="<<AbsDmin<<"/"<<AbsDmax<<" AbsDmedian="<<AbsDmedian<<" AbsDmedianrms="<<AbsDmedianrms<<endl;

	SimilarityData->Dmin= AbsDmin;
	SimilarityData->Dmax= AbsDmax;
	SimilarityData->Dmedian= AbsDmedian;
	SimilarityData->Dmedianrms= AbsDmedianrms;
	SimilarityData->Emin= Emin;
	SimilarityData->Emax= Emax;

	//Normalize similarity matrix by rows
	for(int i=0;i<(SimilarityData->AdjacencyMatrix)->GetNrows();i++){
		double sum= 0;	
		for(int j=0;j<(SimilarityData->AdjacencyMatrix)->GetNcols();j++) {
			sum+= (*(SimilarityData->AdjacencyMatrix))(i,j);
		}
		if(sum!=0) {
			for(int j=0;j<(SimilarityData->AdjacencyMatrix)->GetNcols();j++) (*(SimilarityData->AdjacencyMatrix))(i,j)/= sum;
		}
	}//end loop rows

	//## Sort dissimilarity per row	and store sort index
	for(unsigned int i=0;i<DissMatrix.size();i++){
		std::vector<size_t> sort_index;//sorting index
		std::vector<double> sorted;
		CodeUtils::sort( DissMatrix[i],sorted,sort_index);
		for(unsigned int k=0;k<sort_index.size();k++){
			int index= sort_index[k];
			double diss= DissMatrix[i][index];
			//if(diss<=1.e-12 || i==index) continue;
			if(diss<=0 || (int)(i)==index) continue;
			
			//(SimilarityData->DissimilaritySortIndexMatrix)[i][k]= index;
			(SimilarityData->DissimilaritySortIndexMatrix)[i].push_back(index);
		}
	}//end loop rows

	//Clear up
	EdgenessMatrix->Delete();

	return SimilarityData;

}//close ComputeRegionSimilarity()
*/


int SLICUtils::TagRegions(std::vector<Region*>& regions,Image* binaryMap_bkg,Image* binaryMap_signal){

	if(!binaryMap_bkg || !binaryMap_signal ) {
		ERROR_LOG("No binary maps provided!");
		return -1;
	}
	
	int nRegions= (int)regions.size();
	if(nRegions<=0) {
		ERROR_LOG("No regions available, nothing to be tagged!");
		return -1;
	}
	
	DEBUG_LOG("Tag regions...");
	int nBkgReg= 0;
	int nSignalReg= 0;
	int nUntaggedReg= 0;

	for(int i=0;i<nRegions;i++){	
		int nPix= regions[i]->NPix;
		int nBkg= 0;
		int nUntagged= 0;
		int nSignal= 0;
		if(nPix<=0) {
			regions[i]->Tag= Region::eUntagged;
			continue;
		}

		for(int j=0;j<nPix;j++){//loop on pixels inside region
			long int pixId= (regions[i]->GetPixel(j))->id;
			bool isSignal= (binaryMap_signal->GetBinContent(pixId)>0);
			bool isBkg= (binaryMap_bkg->GetBinContent(pixId)>0);
			if(isBkg && !isSignal) nBkg++;
			else if(isSignal && !isBkg) nSignal++;
			else nUntagged++;
		}//end loop pixels in region
		
		//Tag using a majority criterion
		regions[i]->Tag= Region::eUntagged;
		if(nSignal>nBkg && nSignal>nUntagged) {
			regions[i]->Tag= Region::eSignalTag;
			nSignalReg++;
		}
		else if(nBkg>nSignal && nBkg>nUntagged) {
			regions[i]->Tag= Region::eBkgTag;
			nBkgReg++;
		}
		else if(nUntagged>nSignal && nUntagged>nBkg) {
			regions[i]->Tag= Region::eUntagged;
			nUntaggedReg++;
		}
	}//end loop regions

	INFO_LOG("(nS,nB,nU)=("<<nSignalReg<<","<<nBkgReg<<","<<nUntaggedReg<<")");

	return 0;

}//close TagRegions()

int SLICUtils::CountTaggedRegions(std::vector<Region*>const& regions,int& nSig,int& nBkg,int& nUntagged){

	//Check regions
	int nRegions= (int)regions.size();
	if(nRegions<=0){
		nSig= 0;
		nBkg= 0;
		nUntagged= 0;
		ERROR_LOG("No regions given!");
		return -1;
	}

	nBkg= 0;
	nSig= 0;
	nUntagged= 0;
	for(unsigned int i=0;i<regions.size();i++) {
		int tag= regions[i]->Tag;
		if(tag==Region::eBkgTag) nBkg++;
		else if(tag==Region::eSignalTag) nSig++;
		else if(tag==Region::eUntagged) nUntagged++;
	}//end loop regions
			
	return 0;

}//close CountTaggedRegions()



}//close namespace


