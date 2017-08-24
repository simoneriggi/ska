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
* @file SLICSegmenter.cc
* @class SLICSegmenter
* @brief SLICSegmenter
*
* In this class, an over-segmentation is created of an image, provided by the
* step-size (distance between initial cluster locations) and the colour
* distance parameter.
* @author S. Riggi
* @date 15/06/2015
*/

#include <SLICSegmenter.h>
#include <Region.h>
#include <SLICData.h>
#include <SLICUtils.h>
#include <StatsUtils.h>
#include <CodeUtils.h>
#include <Image.h>

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


ClassImp(Caesar::SLICSegmenter)


namespace Caesar {


SLICSegmenter::SLICSegmenter() {
	
}//close constructor

SLICSegmenter::~SLICSegmenter(){
	
}//close destructor


int SLICSegmenter::CheckData(SLICData const& slicData){

	//Check given slic input data
	if(!slicData.inputImg){
		ERROR_LOG("Null ptr to input image in slicData!");
		return -1;
	}
	if(!slicData.edgeImg){
		ERROR_LOG("Null ptr to edge image in slicData!");
		return -1;
	}	
	//if(!slicData.pixelLabels){
	//	ERROR_LOG("Null ptr to pixel labels in slicData!");
	//	return -1;
	//}

	int nRegions= slicData.GetNRegions();
	if(nRegions<=1){
		WARN_LOG("No (or too few) regions (NR="<<nRegions<<" in slicData!");
		return -1;
	}

	return 0;

}//close Init()


int SLICSegmenter::FindSegmentation(SLICData const& slicData,SLICData& segmSlicData,double SPMergingRegularization,bool SPMergingIncludeSpatialPars,bool use2ndNeighborsInSPMerging,int minMergedSP,double SPMergingRatio, double SPMergingMaxDissRatio,double SPMergingMaxDissRatio_2ndNeighbor,double SPMergingDissThreshold){

	//## Check input data
	INFO_LOG("Checking initial slic data given...");
	if(CheckData(slicData)<0){
		ERROR_LOG("Invalid slicData given (missing edge image information or empty region list?)!");
		return -1;
	}
	
	//## Create mapping between region id and index
	std::map<int,int> regionIdMap_original= slicData.GetRegionIdMap();
	
	//## Copy given slic data
	SLICData* slicData_segmMultiStep= new SLICData;//for 1st stage segmentation
	*slicData_segmMultiStep= slicData;
	
	segmSlicData.Clear();//for final hierarchical merging
	segmSlicData= slicData;
	
	//## Adaptively merge superpixels using max similarity measure
	INFO_LOG("Adaptively merge superpixels using max similarity measure...");
	int status= MultiStepSPMerger(slicData,*slicData_segmMultiStep,SPMergingRegularization,SPMergingIncludeSpatialPars,use2ndNeighborsInSPMerging);
	if(status<0){
		ERROR_LOG("Multi-step superpixel merger failed!");
		return -1;
	}

	//## Check first stage segmentation
	int nSig= 0;
	int nBkg= 0;
	int nUntagged= 0;
	INFO_LOG("Counting tagged regions after adaptive merging stage...");
	SLICUtils::CountTaggedRegions(slicData_segmMultiStep->regions,nSig,nBkg,nUntagged);
	if(nSig<=0 || slicData_segmMultiStep->GetNRegions()<=1){
		WARN_LOG("No signal regions (or no regions) found after adaptive merging, stop processing.");
		return 0;
	}


	//## Create a copy of original regions and change tag according to segmented regions
	INFO_LOG("Create a copy of original regions and change tag according to segmented regions...");
	for(int i=0;i<segmSlicData.GetNRegions();i++) {
		(segmSlicData.regions)[i]->Tag= Region::eBkgTag;
	}//end loop original regions

	for(int i=0;i<slicData_segmMultiStep->GetNRegions();i++) {
		Region* thisRegion= slicData_segmMultiStep->GetRegion(i);
		int regionId= thisRegion->Id;
		int regionTag= thisRegion->Tag;
		int regionIndex_original= regionIdMap_original[regionId]; 

		Region* theOriginalRegion= segmSlicData.GetRegion(regionIndex_original);
		theOriginalRegion->Tag= regionTag;//set new tag
		
		int nSubRegions= (slicData_segmMultiStep->regions)[i]->GetNSubRegions();
		for(int j=0;j<nSubRegions;j++){
			int subregionId= (slicData_segmMultiStep->regions)[i]->GetSubRegionId(j);
			int subregionIndex_original= regionIdMap_original[subregionId]; 
			(segmSlicData.regions)[subregionIndex_original]->Tag= regionTag;
		}	
	}//end loop regions
	
	INFO_LOG("Remove slic data computed in adaptive merging stage...");
	if(slicData_segmMultiStep){
		delete slicData_segmMultiStep;
		slicData_segmMultiStep= 0;
	}
	
	//## Hierarchically merge signal regions
	INFO_LOG("Hierarchically merge signal regions (#"<<segmSlicData.GetNRegions()<<" regions present)...");
	if(segmSlicData.GetNRegions()>2 && nSig>0){
		bool includeSpatialPars= true;
		int nRegions_beforeMerging= segmSlicData.GetNRegions();
		int mergingStatus= SPHierarchicalMerger(segmSlicData,Region::eSignalTag,Region::eSignalTag,minMergedSP,SPMergingRatio,SPMergingRegularization,includeSpatialPars,use2ndNeighborsInSPMerging,SPMergingMaxDissRatio,SPMergingMaxDissRatio_2ndNeighbor,SPMergingDissThreshold);
		if(mergingStatus<0){
			ERROR_LOG("Merging of signal regions failed!");
			return -1;
		}	
		int nRegions_afterMerging= segmSlicData.GetNRegions();
		INFO_LOG("Merged signal regions (NR="<<nRegions_afterMerging-nRegions_beforeMerging<<")");
	}//close if
	

	
	//## Finally add all background regions together
	INFO_LOG("Finally add all background regions together...");
	int bkgLabel= -1;
	Region* bkgRegion= new Region;
	bkgRegion->Id= 0;
	//std::vector<Region*> signalPlusMergedBkgRegions;
	bool isLabelSet= false;
	bool copyPixels= true;
	int nBkgRegions= 0;
	int nSignalRegions= 0;
	std::map<int,int> relabelMap;
	std::vector<int> regionsToBeDeleted; 

	for(int i=0;i<segmSlicData.GetNRegions();i++) {
		Region* thisRegion= segmSlicData.GetRegion(i);
		int tag= thisRegion->Tag;
		int regionId= thisRegion->Id;
		relabelMap.insert( std::pair<int,int>(regionId,regionId) );

		if(tag==Region::eBkgTag) {
			if(!isLabelSet){
				bkgLabel= regionId;
				bkgRegion->Tag= Region::eBkgTag;
			
				long int imageSizeX, imageSizeY;
				thisRegion->GetImageSize(imageSizeX,imageSizeY);

				double imageXmin, imageXmax, imageYmin, imageYmax;
				thisRegion->GetImageRange(imageXmin,imageXmax,imageYmin,imageYmax);
	
				double Smin, Smax;
				thisRegion->GetImageSRange(Smin,Smax);

				double ScurvMin, ScurvMax;
				thisRegion->GetImageScurvRange(ScurvMin,ScurvMax);
			
				double SedgeMin, SedgeMax;	
				thisRegion->GetImageSedgeRange(SedgeMin,SedgeMax);

				bkgRegion->SetImageSize(imageSizeX,imageSizeY);
				bkgRegion->SetImageRange(imageXmin,imageXmax,imageYmin,imageYmax);
				bkgRegion->SetImageSRange(Smin,Smax);
				bkgRegion->SetImageScurvRange(ScurvMin,ScurvMax);
				bkgRegion->SetImageSedgeRange(SedgeMin,SedgeMax);
				bkgRegion->SetImageRMS(thisRegion->GetImageRMS());
	
				isLabelSet= true;
			}//close if
			else {
				regionsToBeDeleted.push_back(i);
			}

			relabelMap[regionId]= bkgLabel;
		
			bkgRegion->AddRegion(thisRegion,true,copyPixels);		
			nBkgRegions++;
		}//close if bkg regions
		else if(tag==Region::eSignalTag){
			nSignalRegions++;
		}
		else{
			WARN_LOG("Untagged region present (this should not occur CHECK!!!!)");
		}
	}//end loop regions


	bkgRegion->Id= bkgLabel;
	
	//Remove merged regions and compute region pars	
	INFO_LOG("Remove merged regions and compute region pars	...");
	CodeUtils::DeleteItems(segmSlicData.regions, regionsToBeDeleted);
	for(int i=0;i<segmSlicData.GetNRegions();i++){
		segmSlicData.regions[i]->ComputeStats(true,false);
	}//end loop merged regions

	
	//Update pixel labels	
	INFO_LOG("Update pixel labels	...");
	for(size_t i=0;i<(segmSlicData.labels).size();i++) {
		for(size_t j=0;j<(segmSlicData.labels)[i].size();j++) {
			int oldLabel= (segmSlicData.labels)[i][j];
			int newLabel= relabelMap[oldLabel];
			(segmSlicData.labels)[i][j]= newLabel;
		}
	}
	
	return 0;

}//close FindSegmentation()


int SLICSegmenter::MultiStepSPMerger(SLICData const& slicData,SLICData& segmSlicData,double SPMergingRegularization,bool SPMergingIncludeSpatialPars,bool use2ndNeighborsInSPMerging){
	
	//## Check if regions have been tagged
	//## Return if no signal-tagged regions is present
	int nSig= 0;
	int nBkg= 0;
	int nUntagged= 0;
	int status= SLICUtils::CountTaggedRegions(segmSlicData.regions,nSig,nBkg,nUntagged);
	if(status<0 || (nSig==0 && nBkg==0) ){
		ERROR_LOG("Regions are not tagged (you must tag regions before running segmentation)!");
		return -1;
	}	
	if(nSig==0){//Tag all untagged regions as bkg and end!		
		WARN_LOG("No signal-tagged regions available, tagging all regions as background...");
		for(int i=0;i<segmSlicData.GetNRegions();i++) {
			int tag= (segmSlicData.regions)[i]->Tag;
			if(tag==Region::eUntagged) (segmSlicData.regions)[i]->Tag= Region::eBkgTag;
		}//end loop regions	
		return 0;
	}//close if
	 
	//## Run region merging
	bool stopAlgo= false;
	int nTotMergedRegions= 0;
	int stageCounter= 0;

	while(!stopAlgo){

		//## Update tagging stats
		nTotMergedRegions= 0;
		if(stageCounter>0){
			SLICUtils::CountTaggedRegions(segmSlicData.regions,nSig,nBkg,nUntagged);
		}
		stageCounter++;
		INFO_LOG("Stage no. "<<stageCounter<<" NR="<<segmSlicData.GetNRegions()<<": (nBkg,nSig,nUntagged)=("<<nBkg<<","<<nSig<<","<<nUntagged<<")");
		
		//## == STAGE 1==
		//## Merge non-tagged regions with bkg-tagged regions
		if(nBkg>0 && nUntagged>0) {//start 1st stage if there are untagged regions present	
			INFO_LOG("=================");
			INFO_LOG("==   STAGE 1   ==");
			INFO_LOG("=================");
			INFO_LOG("Start 1st stage...");		
			int nregions_preStage1= segmSlicData.GetNRegions();
			SPMaxSimilarityMerger(segmSlicData,Region::eBkgTag,Region::eUntagged,SPMergingRegularization,SPMergingIncludeSpatialPars, use2ndNeighborsInSPMerging);
			int nregions_postStage1= segmSlicData.GetNRegions();
			int nMergedRegions_1stStage= nregions_preStage1-nregions_postStage1;
			nTotMergedRegions+= nMergedRegions_1stStage;
		}//close if	

		//Update tagging stats
		SLICUtils::CountTaggedRegions(segmSlicData.regions,nSig,nBkg,nUntagged);
		INFO_LOG("After 1st stage: NR="<<segmSlicData.GetNRegions()<<": (nBkg,nSig,nUntagged)=("<<nBkg<<","<<nSig<<","<<nUntagged<<")");
		
		//## == STAGE 2 ==
		//## Merge non-tagged regions adaptively
		if(nUntagged>0){//start 2nd stage if there are untagged regions available
			INFO_LOG("=================");
			INFO_LOG("==   STAGE 2   ==");
			INFO_LOG("=================");
			INFO_LOG("Start 2nd stage...");
			int nregions_preStage2= segmSlicData.GetNRegions();
			SPMaxSimilarityMerger(segmSlicData,Region::eUntagged,Region::eUntagged,SPMergingRegularization,SPMergingIncludeSpatialPars, use2ndNeighborsInSPMerging);
			int nregions_postStage2= segmSlicData.GetNRegions();
			int nMergedRegions_2ndStage= nregions_preStage2-nregions_postStage2;
			nTotMergedRegions+= nMergedRegions_2ndStage;	
		}//close if
		
		//Update tagging stats
		SLICUtils::CountTaggedRegions(segmSlicData.regions,nSig,nBkg,nUntagged);
		INFO_LOG("After 2nd stage: NR="<<segmSlicData.GetNRegions()<<": (nBkg,nSig,nUntagged)=("<<nBkg<<","<<nSig<<","<<nUntagged<<")");
	
		//## Check end condition
		if(nUntagged<=0) {
			INFO_LOG("NR="<<segmSlicData.GetNRegions()<<": No untagged regions left, algorithm end!");
			stopAlgo= true;
		}
		if(nTotMergedRegions==0){
			INFO_LOG("NR="<<segmSlicData.GetNRegions()<<": No regions merged in all stage, mark remaining as signal and end algorithm!");
			for(int i=0;i<segmSlicData.GetNRegions();i++) {
				int tag= (segmSlicData.regions)[i]->Tag;
				if(tag==Region::eUntagged) (segmSlicData.regions)[i]->Tag= Region::eSignalTag;
			}//end loop regions
			stopAlgo= true;
		}//close if

	}//end while loop 

	return 0;

}//close MultiStepSPMerger()



int SLICSegmenter::SPMaxSimilarityMerger(SLICData& segmSlicData,int mergerTag,int mergedTag,double SPMergingRegularization,bool SPMergingIncludeSpatialPars,bool use2ndNeighborsInSPMerging){

	//## Check regions
	int nRegions= segmSlicData.GetNRegions();
	if(nRegions<=0) {
		WARN_LOG("No regions available, nothing to be merged!");
		return 0;
	}
	
	//## Create the mapping of regionId and vector index
	INFO_LOG("Create the mapping of regionId and region list index...");
	std::map<int,int> regionIdMap= segmSlicData.GetRegionIdMap();

	/*
	//std::map<int,int> regionIdMap_top;	
	//std::vector< std::vector<int> > mergedRegionList;
	
	for(int k=0;k<nRegions;k++){
		int regionId= (segmSlicData.regions)[k]->Id;
		regionIdMap.insert( std::pair<int,int>(regionId,k) );
		//regionIdMap_top.insert( std::pair<int,int>(regionId,k) );
		//mergedRegionList.push_back( std::vector<int>() );
	}//end loop regions
	*/

	//std::map<int,int> regionIdMap_initialSegm= slicData.GetRegionIdMap();

	/*
	for(unsigned int k=0;k<(slicData.regions).size();k++){
		int regionId= (slicData.regions)[k]->Id;
		regionIdMap_initialSegm.insert( std::pair<int,int>(regionId,k) );
	}//end loop regions
	*/

	//## Compute region contour info (neighbors, ...)
	INFO_LOG("Finding region neighbors (NR="<<nRegions<<") ...");
	SLICContourData* contourData= SLICUtils::ComputeBoundaryContours(&segmSlicData);
	if(!contourData){
		ERROR_LOG("Failed to compute the slic contour data!");
		return -1;
	}
	SLICConnectedRegions connectedRegionIds= contourData->connectedRegionIds;

	//## Find neighbors list
	//SLICNeighbors neighbors;
	SLICNeighborCollections neighbors;
	int selectedTag= -1;//mergedTag
	int status= SLICUtils::FindNeighbors(neighbors,&segmSlicData,contourData,use2ndNeighborsInSPMerging,selectedTag,false);
	if(status<0){
		ERROR_LOG("Failed to find neighbors list!");
		return -1;
	}
		
	if(contourData){
		delete contourData;
		contourData= 0;
	}

	/*
	//## Fill list of 2-nd neighbors (if requested to be used in max similarity merging)
	if(use2ndNeighborsInSPMerging){
		std::vector< std::vector<int> > neighborIndexList_2nd;
		for(int i=0; i<nRegions; i++) {
			//int regionId= (segmSlicData.regions)[i]->Id;
			neighborIndexList_2nd.push_back( std::vector<int>() );		
		
			for(unsigned int j=0;j<connectedRegionIds[i].size();j++){//loop over 1st neighbors
				int neighborIndex= connectedRegionIds[i][j];
				for(unsigned int t=0;t<connectedRegionIds[neighborIndex].size();t++){//loop over 2nd neighbors
					int neighborIndex_2nd= connectedRegionIds[neighborIndex][t];	
					std::vector<int>::iterator it= std::find(connectedRegionIds[i].begin(),connectedRegionIds[i].end(),neighborIndex_2nd);
					std::vector<int>::iterator it2= std::find(neighborIndexList_2nd[i].begin(),neighborIndexList_2nd[i].end(),neighborIndex_2nd);
					if( it==connectedRegionIds[i].end() && (it2==neighborIndexList_2nd[i].end() || neighborIndexList_2nd[i].size()==0) ){
						neighborIndexList_2nd[i].push_back(neighborIndex_2nd);
					}
				}//end loop 2nd-neighbors
			}//end loop 1st-neighbors
			connectedRegionIds[i].insert(connectedRegionIds[i].end(),neighborIndexList_2nd[i].begin(),neighborIndexList_2nd[i].end());
		}//end loop regions
	}//close if use 2nd neighbors
	*/

	//## Compute similarity matrix	
	INFO_LOG("Compute region similarity matrix...");
	//SLICSimilarityData* similarityData= SLICUtils::ComputeRegionSimilarity(&segmSlicData,contourData,SPMergingRegularization,SPMergingIncludeSpatialPars,mergedTag);
	SLICSimilarityData* similarityData= SLICUtils::ComputeRegionSimilarity(&segmSlicData,neighbors,SPMergingRegularization);
	TMatrixD* AdjacencyMatrix= similarityData->AdjacencyMatrix;
	//TMatrixD* DissimilarityMatrix= similarityData->DissimilarityMatrix; 
		
	
	//## Compute page rank of segments and sort
	INFO_LOG("Compute page rank ...");
	std::vector<double> ranks;
	status= StatsUtils::ComputePageRank(ranks,AdjacencyMatrix->T());//pass transpose of adjacency matrix
	if(status<0 || ranks.size()<=0){
		WARN_LOG("PageRank failed, cannot perform region merging!");
		return -1;
	}
	std::vector<size_t> sort_index;//sorting index
	std::vector<double> ranks_sorted;
	CodeUtils::sort_descending(ranks,ranks_sorted,sort_index);
		
	
	//## Loop over sorted ranks and select regions to be merged
	INFO_LOG("Start merging loop ...");
	std::vector< std::vector<int> > regionsToBeMerged;
	for(unsigned int i=0;i<(segmSlicData.regions).size();i++) regionsToBeMerged.push_back( std::vector<int>() );
	std::vector<int> regionsToBeDeleted;
	std::vector<int> regionsIdToBeDeleted;
	int nMergedRegions= 0;

	for(unsigned int i=0;i<sort_index.size();i++){
		size_t index= sort_index[i];//region index in regions list
		int regionId= (segmSlicData.regions)[index]->Id;
		int regionTag= (segmSlicData.regions)[index]->Tag;
		
		if(regionTag!=mergerTag && mergerTag!=-1) {
			INFO_LOG("Skip region no. "<<i<<" (id="<<regionId<<", tag="<<regionTag<<")...");
			continue;
		}
			
		//Check if this seed was not already merged by a previous (best ranked) region
		std::vector<int>::iterator seedfinderIt= std::find(regionsToBeDeleted.begin(),regionsToBeDeleted.end(),index);
		if(seedfinderIt!=regionsToBeDeleted.end()){
			INFO_LOG("Seed ranked region (id="<<regionId<<") was already selected for merging in a previous node, skip this merging!");
			continue;
		}

	
		//Loop over untagged neighbors and find best merging
		int nGoodNeighbors= 0;
		std::vector<SLICNeighborData> thisNeighbors= neighbors[index].GetNeighbors();

		//for(unsigned int j=0;j<connectedRegionIds[index].size();j++){//loop over 1st-neighbors
		for(unsigned int j=0;j<thisNeighbors.size();j++){//loop over neighbors of region with ranked index

			//int neighborIndex= connectedRegionIds[index][j];
			int neighborIndex= thisNeighbors[j].Index;
			int neighborId= (segmSlicData.regions)[neighborIndex]->Id;
			int neighborTag= (segmSlicData.regions)[neighborIndex]->Tag;

			if(mergerTag==Region::eBkgTag && neighborTag==mergerTag){//for bkg-untagged merging look at neighbors with tag different than bkg
				continue;
			}
			if(mergerTag==Region::eUntagged && neighborTag!=mergedTag){//for untagged-untagged merging look at neighbors with 'untagged' tag
				continue;
			}
			nGoodNeighbors++;
			
			//double Delta_ij= (*DissimilarityMatrix)(index,neighborIndex);
			double Delta_ij= thisNeighbors[j].Dsym;
	
			//Loop over neighbors of this neighbors (2nd neighbors)
			int closerIndex= -1;
			double closerDiss= 1.e+99;	
			std::vector<SLICNeighborData> thisNeighbors_2nd= neighbors[neighborIndex].GetNeighbors();

			for(unsigned int k=0;k<thisNeighbors_2nd.size();k++){
			//for(unsigned int k=0;k<connectedRegionIds[neighborIndex].size();k++){
				//int neighborIndex_2nd= connectedRegionIds[neighborIndex][k];
				//double Delta_jk= (*DissimilarityMatrix)(neighborIndex,neighborIndex_2nd);		
				
				int neighborIndex_2nd= thisNeighbors_2nd[k].Index;
				double Delta_jk= thisNeighbors_2nd[k].Dsym; 

				if(Delta_jk<closerDiss && Delta_jk>0){	
					closerDiss= Delta_jk;	
					closerIndex= neighborIndex_2nd;
				}
			}//end loop 2nd neighbors

			//If diss Dij is the maximum among Djk(s) merge ij!
			if(closerIndex!=(int)index) continue;


			//Check if this closer region was not already selected to be merged to another node previously	
			std::vector<int>::iterator finderIt= std::find(regionsToBeDeleted.begin(),regionsToBeDeleted.end(),neighborIndex);				
			if(finderIt!=regionsToBeDeleted.end()){	
				INFO_LOG("Closer neighbor no. "<<neighborIndex<<" (id="<<neighborId<<") was already selected for merging in a previous node, skip this merging!");
				continue;
			}

			//Check if this neighbor was not previously selected as a merger
			if(regionsToBeMerged[neighborIndex].size()>0){
				INFO_LOG("Closer neighbor no. "<<neighborIndex<<" (id="<<neighborId<<") was previously selected as a merger, skip this merging!");
				continue;
			}

			//Merging selected!
			INFO_LOG("New merging ("<<index<<"-"<<neighborIndex<<"), Delta_ij="<<Delta_ij<<" closerDiss="<<closerDiss);
			//int regionIndex_A= regionIdMap_top[regionId];
			//int regionIndex_B= regionIdMap_top[neighborId];		

			regionsToBeMerged[index].push_back(neighborIndex);
			regionsToBeDeleted.push_back(neighborIndex);
			regionsIdToBeDeleted.push_back(neighborId);		
			//mergedRegionList[regionIndex_A].push_back(regionIndex_B);
			nMergedRegions++;

		}//end loop neighbors
				
	}//end loop ranks

	//## Delete contour data
	/*
	if(contourData){
		delete contourData;
		contourData= 0;
	}
	*/
	if(similarityData){
		delete similarityData;
		similarityData= 0;
	}
		
	//## Merge regions and update region and label list
	if(nMergedRegions>0){
		INFO_LOG("Merge the selected regions...");
		std::map<int,int> newLabelMap;

		for(unsigned int k=0;k<(segmSlicData.regions).size();k++){
			int regionId= (segmSlicData.regions)[k]->Id;	
			int regionTag= (segmSlicData.regions)[k]->Tag;
			newLabelMap.insert( std::pair<int,int>(regionId,regionId) );
			for(unsigned int j=0;j<regionsToBeMerged[k].size();j++){
				int mergedRegionIndex= regionsToBeMerged[k][j];
				int mergedRegionId= (segmSlicData.regions)[mergedRegionIndex]->Id;
				int mergedRegionTag= (segmSlicData.regions)[mergedRegionIndex]->Tag;
				newLabelMap[mergedRegionId]= regionId;

				//Change tag in original segmentation
				//int mergedRegionIndex_originalSegm= (slicData.regionIdMap)[mergedRegionId];//Original regions
				//int mergedRegionId_originalSegm= (slicData.regions)[mergedRegionIndex_originalSegm]->Id;
				//cout<<"SLICSegmenter::SPMaxSimilarityMerger(): INFO: Region no. "<<k<<" (id="<<regionId<<",tag="<<regionTag<<") : merging region id="<<mergedRegionId<<", tag="<<mergedRegionTag<<" (index="<<mergedRegionIndex<<",orig index="<<mergedRegionIndex_originalSegm<<", orig id="<<mergedRegionId_originalSegm<<")"<<endl;		

				INFO_LOG("Region no. "<<k<<" (id="<<regionId<<",tag="<<regionTag<<") : merging region id="<<mergedRegionId<<", tag="<<mergedRegionTag<<" (index="<<mergedRegionIndex<<")");			

				(segmSlicData.regions)[k]->AddRegion((segmSlicData.regions)[mergedRegionIndex],true,true);
				//(segmSlicData.regions)[k]->AddSubRegionId(mergedRegionId);//NEED TO CHECK!!!
			
			}//end loop regions to be merged
		}//end loop regions
	
		//## Delete aggregated region from region list and index map
		INFO_LOG("Deleting regions aggregated in this step from the main list...");
		CodeUtils::DeleteItems((segmSlicData.regions), regionsToBeDeleted);
		for(size_t k=0;k<regionsIdToBeDeleted.size();k++) regionIdMap.erase(regionsIdToBeDeleted[k]);

		//## Update map and recompute parameters & contours (will be used by nearest neighbors search)
		INFO_LOG("Updating region parameters & contours...");
		for(unsigned int k=0;k<(segmSlicData.regions).size();k++){
			(segmSlicData.regions)[k]->ComputeStats(false,true);
			int regionId= (segmSlicData.regions)[k]->Id;
			regionIdMap[regionId]= k;
		}//end loop regions

		//## Update pixel labels
		INFO_LOG("Updating pixel labels...");
		for(unsigned int i=0;i<(segmSlicData.labels).size();i++) {
			for(unsigned int j=0;j<(segmSlicData.labels)[i].size();j++) {
				int oldLabel= (segmSlicData.labels)[i][j];
				int newLabel= newLabelMap[oldLabel];
				(segmSlicData.labels)[i][j]= newLabel;
			}//end loop
		}//end loop 
	}//close if merge regions
		
	INFO_LOG(nMergedRegions<<"/"<<segmSlicData.GetNRegions()<<" regions merged at this stage...");

	return 0;

}//close SPMaxSimilarityMerger()


int SLICSegmenter::SPHierarchicalMerger(SLICData& slicData,int mergerTag,int mergedTag,int minMergedSP,double SPMergingRatio,double SPMergingRegularization,bool includeSpatialPars, bool use2ndNeighborsInSPMerging,double SPMergingMaxDissRatio,double SPMergingMaxDissRatio_2ndNeighbor,double SPMergingDissThreshold){

	//## Check regions
	int nRegions= slicData.GetNRegions();
	if(nRegions<=0) {
		WARN_LOG("No regions available, nothing to be merged!");
		return 0;
	}

	//## Init region list to be passed to algorithm
	//cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Copying initial superpixel partition to tmp partition..."<<endl;
	//std::vector<Region*> regions;
	//regions.assign((slicData.regions).begin(),(slicData.regions).end());
	
	//## Create the mapping of regionId and vector index
	INFO_LOG("Create the mapping of regionId and region list index...");
	std::map<int,int> regionIdMap= slicData.GetRegionIdMap();

	/*
	//std::map<int,int> regionIdMap_top;	
	//std::vector<std::vector<int>> mergedRegionList;
	for(int k=0;k<nRegions;k++){
		int regionId= (slicData.regions)[k]->Id;
		regionIdMap.insert( std::pair<int,int>(regionId,k) );
		//regionIdMap_top.insert( std::pair<int,int>(regionId,k) );
		//mergedRegionList.push_back( std::vector<int>() );
	}//end loop regions
	*/
	/*
	//Init pixel labels: copy initial pixel labels to tmp list
	std::vector< std::vector<long int> > labels;
	for(unsigned int i=0;i<(slicData.labels).size();i++) {
		labels.push_back( std::vector<long int>() );
		labels[i].assign((slicData.labels)[i].begin(),(slicData.labels)[i].end());
	}
	*/

	//## Run a hierarchical clustering till all segments are merged in one
	int hierarchyLevel= 0;
	double DissMedian0= 0;
	int nMergedRegionsInHierarchyLevel= 0;
	int itemPos= -1;
	
	INFO_LOG("Starting hierarchical merging...");

	while( (slicData.GetNRegions())>minMergedSP ){
		nMergedRegionsInHierarchyLevel= 0;
		
		//## Compute region id map
		std::map<int,int> regionIdMap= slicData.GetRegionIdMap();

		//## Compute region contour info (neighbors, ...)
		INFO_LOG("Finding region neighbors at hierarchy level "<<hierarchyLevel<<", NR="<<slicData.GetNRegions()<<" ...");
		SLICContourData* contourData= SLICUtils::ComputeBoundaryContours(&slicData);
		SLICConnectedRegions connectedRegionIds= contourData->connectedRegionIds;

		//## Find neighbors list
		SLICNeighborCollections neighbors;
		int selectedTag= mergedTag;
		int status= SLICUtils::FindNeighbors(neighbors,&slicData,contourData,use2ndNeighborsInSPMerging,selectedTag,includeSpatialPars);
		if(status<0){
			ERROR_LOG("Failed to find neighbors list!");
			return -1;
		}
		if(contourData){
			delete contourData;
			contourData= 0;
		}
		/*
		std::vector< std::vector<int> > neighborIndexList_1st;
		std::vector< std::vector<int> > neighborIndexList_2nd;

		for(int i=0; i<slicData->GetNRegions(); i++) {
			int regionId= (slicData.regions)[i]->Id;
			//Fill list of 1st and 2-nd neighbors	
			neighborIndexList_1st.push_back( std::vector<int>() );
			neighborIndexList_2nd.push_back( std::vector<int>() );
		
			for(unsigned int j=0;j<connectedRegionIds[i].size();j++){//loop over 1st neighbors
				int neighborIndex= connectedRegionIds[i][j];
				int neighborTag= regions[neighborIndex]->fTag;
				int neighborId= regions[neighborIndex]->fId;
				if(mergedTag!=-1 && neighborTag==mergedTag){
					neighborIndexList_1st[i].push_back(neighborIndex);
				}

				for(unsigned int t=0;t<connectedRegionIds[neighborIndex].size();t++){//loop over 2nd neighbors
					int neighborIndex_2nd= connectedRegionIds[neighborIndex][t];
					int neighborTag_2nd= regions[neighborIndex_2nd]->fTag;
					int neighborId_2nd= regions[neighborIndex_2nd]->fId;
					if(mergedTag!=-1 && neighborTag_2nd==mergedTag && neighborId_2nd!=regionId){
						std::vector<int>::iterator it= std::find(connectedRegionIds[i].begin(),connectedRegionIds[i].end(),neighborIndex_2nd);
						std::vector<int>::iterator it2= std::find(neighborIndexList_2nd[i].begin(),neighborIndexList_2nd[i].end(),neighborIndex_2nd);
		
						if( it==connectedRegionIds[i].end() && (it2==neighborIndexList_2nd[i].end() || neighborIndexList_2nd[i].size()==0) ){
							neighborIndexList_2nd[i].push_back(neighborIndex_2nd);
						}	
					}//close if
				}//end loop 2nd-neighbors
			}//end loop 1st-neighbors
		}//end loop regions
		*/


		
		//## Compute similarity matrix	
		INFO_LOG("Compute region similarity at hierarchy level "<<hierarchyLevel<<" ...");
		SLICSimilarityData* similarityData= SLICUtils::ComputeRegionSimilarity(&slicData,neighbors,SPMergingRegularization);
		TMatrixD* AdjacencyMatrix= similarityData->AdjacencyMatrix;
		double DissMedian= similarityData->Dmedian;
 		if(hierarchyLevel==0) {
			DissMedian0= DissMedian;
		}

		/*
		Img* edgeImg= fEdgeFilterImg;
		if(edgeModel==1) edgeImg= fKirschEdgeFilterImg;
		SLICUtils::SLICSimilarityData* similarityData= SLICUtils::ComputeRegionSimilarity(edgeImg,contourData,regions,fSPMergingRegularization,includeSpatialPars,mergedTag);
		TMatrixD* AdjacencyMatrix= similarityData->AdjacencyMatrix;
		TMatrixD* DissimilarityMatrix= similarityData->DissimilarityMatrix; 
		TMatrixD* AbsDissimilarityMatrix= similarityData->AbsDissimilarityMatrix; 
		TMatrixD* NeighborMatrix= similarityData->NeighborMatrix; 
		std::vector< std::vector<int> > DissSortIndexes= similarityData->DissimilaritySortIndexMatrix;
		AbsDissMedian= similarityData->Dmedian;
 		AbsDissMedianRMS= similarityData->Dmedianrms;
 		AbsDissMin= similarityData->Dmin;
		AbsDissMax= similarityData->Dmax;

		if(hierarchyLevel==0) {
			AbsDissMedian0= AbsDissMedian;
			AbsDissMedianRMS0= AbsDissMedianRMS;
		}
		
		double MSE= 0.;
		for(unsigned int k=0;k<regions.size();k++){		
			double M2= regions[k]->fM2;
			MSE+= M2;
		}//end loop regions
		MSE/= N;
		spMergeInfo.MSE= MSE;
		*/

		//## Compute page rank of segments and sort
		INFO_LOG("Compute page rank at hierarchy level "<<hierarchyLevel<<" ...");
		std::vector<double> ranks;
		status= StatsUtils::ComputePageRank(ranks,AdjacencyMatrix->T());//pass transpose of adjacency matrix
		if(status<0 || ranks.size()<=0){
			WARN_LOG("PageRank failed, cannot perform region merging!");
			return -1;
		}
		
		std::vector<size_t> sort_index;//sorting index
		std::vector<double> ranks_sorted;
		CodeUtils::sort_descending(ranks,ranks_sorted,sort_index);
		
		//## Count max number of mergeable regions
		int nMaxMergeableRegions= 0;
		for(unsigned int i=0;i<sort_index.size();i++){
			size_t index= sort_index[i];//region index in regions list
			int regionTag= slicData.regions[index]->Tag;
			if(regionTag!=mergerTag && mergerTag!=-1) continue;
			nMaxMergeableRegions++;
		}//end loop regions
		int maxRegionsToMerge= std::round(nMaxMergeableRegions*SPMergingRatio);

		//## Loop over sorted ranks and select regions to be merged
		int nMergedRegions= 0;
		int nMergeableRegions= 0;
		std::vector< std::vector<int> > regionsToBeMerged;
		for(int i=0;i<slicData.GetNRegions();i++) regionsToBeMerged.push_back( std::vector<int>() );
		std::vector<int> regionsToBeDeleted;
		std::vector<int> regionsIdToBeDeleted;

		for(unsigned int i=0;i<sort_index.size();i++){
			size_t index= sort_index[i];//region index in regions list
			//double thisRank= ranks[index];//region rank
			int regionId= (slicData.regions)[index]->Id;
			int regionTag= (slicData.regions)[index]->Tag;
			//int mapIndex= regionIdMap[regionId];
			
			if(regionTag!=mergerTag && mergerTag!=-1) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Skip ranked seed region (id="<<regionId<<",tag="<<regionTag<<") as tag is different from mergerTag!"<<endl;
				continue;
			}

			//Stop merging above nmerge threshold						
			if(nMergeableRegions+1>maxRegionsToMerge) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Maximum number of regions mergeable reached for this hierarchy level ("<<nMergeableRegions+1<<">"<<maxRegionsToMerge<<")"<<endl;
				break;
			}
			
			//Check if this seed was not already merged by a previous (best ranked) region
			if(!CodeUtils::FindItem(regionsIdToBeDeleted,regionId,itemPos)){
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Seed ranked region (id="<<regionId<<") was already selected for merging in a previous node, skip this merging!"<<endl;
				continue;
			}

			/*
			std::vector<int>::iterator seedfinderIt= std::find(regionsIdToBeDeleted.begin(),regionsIdToBeDeleted.end(),regionId);
			if(seedfinderIt!=regionsIdToBeDeleted.end()){
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Seed ranked region (id="<<regionId<<") was already selected for merging in a previous node, skip this merging!"<<endl;
				continue;
			}
			*/

			//Chech if the seed region has any neighbors (according to the selected merged tag)
			int NN= (int)neighbors[index].GetN();
			if(NN<=0){
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Skip ranked seed region (id="<<regionId<<",tag="<<regionTag<<") as no neighbors are available!"<<endl;
				continue;
			}
			/*
			int NN_1st= (int)neighborIndexList_1st[index].size();
			int NN_2nd= (int)neighborIndexList_2nd[index].size();	
			if(NN_1st<=0 && NN_2nd<=0){
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Skip ranked seed region (id="<<regionId<<",tag="<<regionTag<<") as no neighbors are available!"<<endl;
				continue;
			}
			*/

			//Get closest neighbor
			/*
			if(DissSortIndexes[index].size()<=0) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Skip seed ranked region (id="<<regionId<<") as no neighbors are available!"<<endl;
				continue;
			}
			int closerNeighborIndex= DissSortIndexes[index][0];
			int closerNeighborId= (slicData.regions)[closerNeighborIndex]->Id;
			int closerNeighborTag= (slicData.regions)[closerNeighborIndex]->Tag;
			*/

			int closerNeighborPosInCollection= neighbors[index].FindCloserByDissTot();//this is not the region index in original vector!!!
			SLICNeighborData* closerNeighbor= neighbors[index].GetNeighbor(closerNeighborPosInCollection);
			int closerNeighborId= closerNeighbor->Id;
			int closerNeighborTag= closerNeighbor->Tag; 
			int closerNeighborOrder= closerNeighbor->Order;
			int closerNeighborIndex= closerNeighbor->Index;//regionIdMap[closerNeighborId];
			
			if(closerNeighborId==regionId) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Skip neighbor (id="<<closerNeighborId<<") as equal to seed region!"<<endl;
				continue;//skip same region
			}
			if(closerNeighborTag!=mergedTag && mergedTag!=-1) {	
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Skip ranked seed region (id="<<regionId<<",tag="<<regionTag<<") as neighbor (id="<<closerNeighborId<<",tag="<<closerNeighborTag<<") tag is different from mergedTag!"<<endl;
				continue;
			}
			nMergeableRegions++;
			
			
			double Delta_ij= closerNeighbor->Dtot;
			double Delta_ji= closerNeighbor->Dtot_n;
			double AbsDelta_ij= closerNeighbor->D;
			
			//Find current region index among best neighbors of closer
			int N_half= ceil(neighbors[closerNeighborIndex].GetN()/2.);
			if(N_half>=2){
				bool isSeedAmongNeighborCloser= neighbors[closerNeighborIndex].IsIdAmongNClosersByTotDiss(regionId,N_half);
				if(!isSeedAmongNeighborCloser){
					cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Neighbor region id="<<closerNeighborId<<" rejected for merging with region "<<regionId<<" since this is not among the closest neighbors of the selected neighbor!"<<endl;
					continue;
				}
			}

			/*
			double Delta_ij= (*DissimilarityMatrix)(index,closerNeighborIndex);//include edgeness
			double Delta_ji= (*DissimilarityMatrix)(closerNeighborIndex,index);//include edgeness
			double AbsDelta_ij= (*AbsDissimilarityMatrix)(index,closerNeighborIndex);//without edgeness
			double AbsDelta_ji= (*AbsDissimilarityMatrix)(closerNeighborIndex,index);//without edgeness
			int closerNeighborness= (*NeighborMatrix)(index,closerNeighborIndex);
			
			//Get neighbors of the closer region
			std::vector<int> closerNeighborNNIndexes;
			closerNeighborNNIndexes.assign(
				DissSortIndexes[closerNeighborIndex].begin(),
				DissSortIndexes[closerNeighborIndex].begin() + ceil(DissSortIndexes[closerNeighborIndex].size()/2.)
			);

			//Find current region index among best neighbors of closer
			std::vector<int>::iterator it= std::find(closerNeighborNNIndexes.begin(),closerNeighborNNIndexes.end(),index);
			if(it==closerNeighborNNIndexes.end() ) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Neighbor region id="<<closerNeighborId<<" rejected for merging with region "<<regionId<<" since this is not among the closest neighbors of selected neighbor!"<<endl;
				continue;
			}
			*/
			
			//## Check similarity difference
			//## 1st neighbors are merged (flood-fill)
			//## 2nd neighbors are merged if their similarities are not too different
			/*
			if(closerNeighborness==1 && (Delta_ji>Delta_ij*SPMergingMaxDissRatio) ) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Neighbor region id="<<closerNeighborId<<" (1st neighbor) rejected for merging with region "<<regionId<<"(Delta_ji="<<Delta_ji<<">"<<fSPMergingMaxDissRatio<<"*Delta_ij="<<Delta_ij<<endl;
				continue;
			}
			else if(closerNeighborness==2 && (Delta_ji>Delta_ij*SPMergingMaxDissRatio_2ndNeighbor) ) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Neighbor region id="<<closerNeighborId<<" (2nd neighbor) rejected for merging with region "<<regionId<<"(Delta_ji="<<Delta_ji<<">"<<fSPMergingMaxDissRatio_2ndNeighbor<<"*Delta_ij="<<Delta_ij<<endl;
				continue;
			}
			*/

			//Check if this closer region was not already selected as a seed merger previously
			if( regionsToBeMerged[closerNeighborIndex].size()>0 ){	
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Closer neighbor (id="<<closerNeighborId<<") selected for merging was before selected as a primary merger, skip this merging!"<<endl;
				continue;
			}

			//Check if this closer region was not already selected to be merged to another node previously	
			if(!CodeUtils::FindItem(regionsIdToBeDeleted,closerNeighborId,itemPos)){
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Closer neighbor (id="<<closerNeighborId<<") was already selected for merging in a previous node, skip this merging!"<<endl;
				continue;
			}

			/*
			std::vector<int>::iterator finderIt= std::find(regionsIdToBeDeleted.begin(),regionsIdToBeDeleted.end(),closerNeighborId);
			if(finderIt!=regionsIdToBeDeleted.end()){
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Closer neighbor (id="<<closerNeighborId<<") was already selected for merging in a previous node, skip this merging!"<<endl;
				continue;
			}
			*/

			
			//## Apply dissimilarity threshold			
			if(use2ndNeighborsInSPMerging && closerNeighborOrder==2){
				double DissThreshold= SPMergingDissThreshold*DissMedian0;
			
				if(Delta_ji>Delta_ij*SPMergingMaxDissRatio_2ndNeighbor){
					cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Neighbor region id="<<closerNeighborId<<" (2nd neighbor) rejected for merging with region "<<regionId<<"(Delta_ji="<<Delta_ji<<">"<<SPMergingMaxDissRatio_2ndNeighbor<<"*Delta_ij="<<Delta_ij<<endl;
					continue;
				}
				if(AbsDelta_ij>DissThreshold){
					cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Region id="<<regionId<<": closer neighbor (id="<<closerNeighborId<<") cannot be merged as dissimilarity is too large (AbsDiss="<<AbsDelta_ij<<", Diss="<<Delta_ij<<">"<<DissThreshold<<")"<<endl;
					continue;
				}	
			}//close if


			/*
			if(closerNeighborness==2 && (Delta_ji>Delta_ij*SPMergingMaxDissRatio_2ndNeighbor) ) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Neighbor region id="<<closerNeighborId<<" (2nd neighbor) rejected for merging with region "<<regionId<<"(Delta_ji="<<Delta_ji<<">"<<fSPMergingMaxDissRatio_2ndNeighbor<<"*Delta_ij="<<Delta_ij<<endl;
				continue;
			}

			if(closerNeighborness==2 && AbsDelta_ij>DissThreshold){
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Region id="<<regionId<<": closer neighbor (id="<<closerNeighborId<<") cannot be merged as dissimilarity is too large (AbsDiss="<<AbsDelta_ij<<", Diss="<<Delta_ij<<">"<<DissThreshold<<")"<<endl;
				continue;
			}
			*/

			cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Region (id="<<regionId<<") merging neighbor (id="<<closerNeighborId<<"): order="<<closerNeighborOrder<<", Delta_ij="<<Delta_ij<<", Delta_ji="<<Delta_ji<<" ratio="<<Delta_ji/Delta_ij<<endl;

			
			//int regionIndex_A= regionIdMap_top[regionId];
			//int regionIndex_B= regionIdMap_top[closerNeighborId];
			
			regionsToBeMerged[index].push_back(closerNeighborIndex);
			regionsToBeDeleted.push_back(closerNeighborIndex);
			regionsIdToBeDeleted.push_back(closerNeighborId);		
			//mergedRegionList[regionIndex_A].push_back(regionIndex_B);
	
			nMergedRegions++;
		}//end loop ranks

		
		//## Delete contour data
		/*
		if(contourData){
			delete contourData;
			contourData= 0;
		}
		*/
		if(similarityData){
			delete similarityData;
			similarityData= 0;
		}
		
		
		//## Merge regions
		cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Merge the selected regions..."<<endl;
		std::map<int,int> newLabelMap;

		for(int k=0;k<slicData.GetNRegions();k++){
			int regionId= (slicData.regions)[k]->Id;
			newLabelMap.insert( std::pair<int,int>(regionId,regionId) );
			for(unsigned int j=0;j<regionsToBeMerged[k].size();j++){
				int mergedRegionIndex= regionsToBeMerged[k][j];
				int mergedRegionId= (slicData.regions)[mergedRegionIndex]->Id;
				newLabelMap[mergedRegionId]= regionId;
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Region no. "<<k<<" (id="<<regionId<<") : merging region id="<<mergedRegionId<<" (index="<<mergedRegionIndex<<")"<<endl;		
				(slicData.regions)[k]->AddRegion((slicData.regions)[mergedRegionIndex],true,true);
			}//end loop regions to be merged
		}//end loop regions
	
		//## Delete aggregated region from region list and index map
		cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Deleting regions aggregated in this step from the main list..."<<endl;
		CodeUtils::DeleteItems((slicData.regions), regionsToBeDeleted);
		for(size_t k=0;k<regionsIdToBeDeleted.size();k++) regionIdMap.erase(regionsIdToBeDeleted[k]);

		//## Update map and recompute parameters & contours (will be used by nearest neighbors search)
		cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Updating region parameters & contours..."<<endl;
		for(int k=0;k<slicData.GetNRegions();k++){
			(slicData.regions)[k]->ComputeStats(false,true);
			int regionId= (slicData.regions)[k]->Id;
			regionIdMap[regionId]= k;
		}//end loop regions
	
		//## Update pixel labels
		cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Updating pixel labels..."<<endl;	
		for(unsigned int i=0;i<(slicData.labels).size();i++) {
			for(unsigned int j=0;j<(slicData.labels)[i].size();j++) {
				int oldLabel= (slicData.labels)[i][j];
				int newLabel= newLabelMap[oldLabel];
				(slicData.labels)[i][j]= newLabel;
			}
		}
	
		nMergedRegionsInHierarchyLevel= nMergedRegions;
		hierarchyLevel++;

		if(nMergedRegionsInHierarchyLevel==0){
			cerr<<"SLICSegmenter::SPHierarchicalMerger(): WARN: No regions merged in this stage, exit to avoid stuck!"<<endl;	
			break;
		}
		cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: "<<nMergedRegionsInHierarchyLevel<<"/"<<slicData.GetNRegions()<<" regions aggregated at this level hierarchy..."<<endl;
		

	}//end while loop


	cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: "<<hierarchyLevel<<" hierarchy levels aggregated: N="<<slicData.GetNRegions()<<" regions left"<<endl;


	/*
	//## Final tag of composite regions according to a majority criterion, that is region tagged as bkg if the majority of sub-regions are non-significative
	cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Tagging significant regions..."<<endl;
	for(int i=0;i<regions.size();i++){
		int regionId= regions[i]->fId;
		int regionTag= regions[i]->fTag;
		int regionIndex_top= regionIdMap_top[regionId];//region index at the top of the hierarchy
		std::vector<int> subregionIds= regions[i]->fSubRegionIds;

		bool isSignificant= inputRegions[regionIndex_top]->fIsSignificative;
		bool isSalient= inputRegions[regionIndex_top]->fIsSalient;
		double salientFraction= 0;
		int nSubRegions= 1+subregionIds.size();
		if(isSalient) salientFraction++;
		bool isAnySignificant= false;
		if(isSignificant) isAnySignificant= true;

		for(unsigned int j=0;j<subregionIds.size();j++){
			int subRegionId= subregionIds[j];
			int subRegionIndex_top= regionIdMap_top[subRegionId];
			//bool subIsSignificant= fRegions[subRegionIndex_top]->fIsSignificative;
			//bool subIsSalient= fRegions[subRegionIndex_top]->fIsSalient;	
			bool subIsSignificant= inputRegions[subRegionIndex_top]->fIsSignificative;
			bool subIsSalient= inputRegions[subRegionIndex_top]->fIsSalient;
			if(subIsSalient) salientFraction++;
			if(subIsSignificant) isAnySignificant= true;
		}
		salientFraction/= (double)nSubRegions;
			
		if( (salientFraction>fSignificantSPRatio) || isAnySignificant ){
			cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Region id="<<regionId<<" tagged as significant as (significanceRatio="<<salientFraction<<">"<<fSignificantSPRatio<<")..."<<endl;
			regions[i]->fIsSignificative= true;
		}
		else {
			regions[i]->fIsSignificative= false;
		}	
	}//end loop regions
		
	//## Tag background?
	if(fUsePixelRatioCut){
		//Set region with maximum number of pixels
		for(int i=0;i<regions.size();i++){
			int regionId= regions[i]->fId;	
			int npix= regions[i]->fNPix;
			double pixelRatio= (double)npix/(double)N;
			if(pixelRatio>fPixelRatioCut) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Region id="<<regionId<<" has too many pixels (ratio="<<pixelRatio<<">"<<fPixelRatioCut<<") and will be tagged as bkg..."<<endl;
				regions[i]->fIsSignificative= false;
			}
		}
	}//close if
	
	//## Copy final results to global data (fRegions, fPixelClusterIds)
	inputRegions.clear();
	inputRegions.assign(regions.begin(),regions.end());
	for(unsigned int i=0;i<labels.size();i++) {
		for(unsigned int j=0;j<labels[i].size();j++) {
			inputLabels[i][j]= labels[i][j];
		}
	}
	*/	
	
	return 0;

}//close SPHierarchicalMerger()

}//close namespace


