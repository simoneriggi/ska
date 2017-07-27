#include <SLICUtils.h>
#include <OutlierDetector.h>

#include <TColor.h>
#include <TMath.h>
#include <TStyle.h>
#include <TColor.h>
#include <TPad.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TText.h>
#include <TPolyLine.h>
#include <Math/QuantFuncMathCore.h>

#include <mathop.h>

#include <chrono>
using namespace std;
using namespace std::chrono;

#define SMALL_NUMBER 0.00000000001;

SLICUtils::SLICUtils() {

}//close constructor


SLICUtils::~SLICUtils() {	

}//close destructor

Img* SLICUtils::GetSegmentedImage(Img* image,std::vector<Region*> regions,int selectedTag,bool normalize,bool binarize){

	//## Check input
	if(!image) {
		cerr<<"SLICUtils::GetSegmentedImage(): ERROR: Null ptr to input image given!"<<endl;
		return 0;
	}

	int nRegions= (int)regions.size();
	if(nRegions<=0) {
		cerr<<"SLICUtils::GetSegmentedImage(): ERROR: No regions available, nothing to be done!"<<endl;
		return 0;
	}

	//## Use image normalization?
	double A= 0;
	double B= 1;
	if(image->HasStats() && normalize){
		Img::StatsData* stats= image->GetPixelStats();
		double NormMin= image->GetMinimum();
		double NormMax= image->GetMaximum();
		A= NormMin - (NormMax-NormMin)*stats->min/(stats->max-stats->min);
		B= (NormMax-NormMin)/(stats->max-stats->min);
	}	
	
	//## Fill image with region means
	Img* segmentedImg= (Img*)image->Clone("segmentedImg");
	segmentedImg->Reset();
	
  for(int i=0;i<nRegions;i++){//loop on regions
		int regionId= regions[i]->fId;	
		int regionTag= regions[i]->fTag;

		if(regionTag!=-1 && regionTag!=selectedTag) continue;
		
		int nPixelsInRegion= (int)regions[i]->fNPix;
		double Mean= regions[i]->fMean;
		double colorValue= A+B*Mean;
		if(binarize) colorValue= 1;

		for(int j=0;j<nPixelsInRegion;j++){//loop on pixels inside region
			int thisPixelId= (regions[i]->fPixelCollection)[j].id;
			segmentedImg->SetBinContent(thisPixelId,colorValue);
		}//end loop pixels in region	
	}//end loop regions
	
	return segmentedImg;

}//close GetSegmentedImage()



Img* SLICUtils::GetSegmentedImage(Img* image,std::vector<Region*> regions,bool drawOnlySignificant,bool drawOnlySalient,bool normalize,bool binarize){

	//## Check input
	if(!image) {
		cerr<<"SLICUtils::GetSegmentedImage(): ERROR: Null ptr to input image given!"<<endl;
		return 0;
	}

	int nRegions= (int)regions.size();
	if(nRegions<=0) {
		cerr<<"SLICUtils::GetSegmentedImage(): ERROR: No regions available, nothing to be done!"<<endl;
		return 0;
	}

	//## Use image normalization?
	double A= 0;
	double B= 1;
	if(image->HasStats() && normalize){
		Img::StatsData* stats= image->GetPixelStats();
		double NormMin= image->GetMinimum();
		double NormMax= image->GetMaximum();
		A= NormMin - (NormMax-NormMin)*stats->min/(stats->max-stats->min);
		B= (NormMax-NormMin)/(stats->max-stats->min);
	}	
	
	//## Fill image with region means
	Img* segmentedImg= (Img*)image->Clone("segmentedImg");
	segmentedImg->Reset();
	
  for(int i=0;i<nRegions;i++){//loop on regions
		int regionId= regions[i]->fId;
		bool isSignificative= regions[i]->fIsSignificative;
		bool isSalient= regions[i]->fIsSalient;
		cout<<"SLICUtils::GetSegmentedImage(): INFO: Region id="<<regionId<<" isSignificative?"<<isSignificative<<endl;
		regions[i]->Dump();

		if(drawOnlySignificant && !isSignificative) continue;
		if(drawOnlySalient && !isSalient) continue;
		
		int nPixelsInRegion= (int)regions[i]->fNPix;
		double Mean= regions[i]->fMean;
		double colorValue= A+B*Mean;
		if(binarize) colorValue= 1;

		for(int j=0;j<nPixelsInRegion;j++){//loop on pixels inside region
			int thisPixelId= (regions[i]->fPixelCollection)[j].id;
			segmentedImg->SetBinContent(thisPixelId,colorValue);
		}//end loop pixels in region	
	}//end loop regions
	
	return segmentedImg;
	
}//close GetSegmentedImage()


int SLICUtils::TagRegions_last(Img* image,Img* saliencyMap,std::vector<Region*> regions,double saliencyThresholdFactor,double saliencyBkgThresholdFactor,int minNPix){

	cout<<"SLICUtils::TagRegions(): INFO: Tag regions..."<<endl;
	if(!image ) return -1;

	if(!image->HasStats()){
		cerr<<"SLICUtils::TagRegions(): ERROR: Image has no stats!"<<endl;
		return -1;
	}
	Img::StatsData* imageStats= image->GetPixelStats();
	double mean= imageStats->mean; 
	double rms= imageStats->rms;

	int nRegions= (int)regions.size();
	if(nRegions<=0) {
		cerr<<"SLICUtils::TagRegions(): ERROR: No regions available, nothing to be tagged!"<<endl;
		return -1;
	}
	
	//## Compute saliency significance map
	if(!saliencyMap){
		cerr<<"SLICUtils::TagRegions(): ERROR: Cannot get the saliency map"<<endl;
		return -1;
	}	
	if(!saliencyMap->HasStats()){
		cerr<<"SLICUtils::TagRegions(): ERROR: Saliency has no stats computed!"<<endl;
		return -1;
	}
	if(!saliencyMap->HasBkgData()){
		cerr<<"SLICUtils::TagRegions(): ERROR: Saliency has no bkg data computed!"<<endl;
		return -1;
	}
	//if(!saliencyMap->ComputeStats(true,false,true)<0){
	//	cerr<<"SLICUtils::TagRegions(): ERROR: Stats computing failed!"<<endl;
	//	return -1;
	//}
	Img::StatsData* saliencyStats= saliencyMap->GetPixelStats();
	double saliencyMedian= saliencyStats->median;
	double saliencyMedianRMS= saliencyStats->medianRMS;
	
	int nHistoBins= 100;
	double saliencyOtsuThr= saliencyMap->FindOtsuThreshold(nHistoBins);
	double saliencyValleyThr= saliencyMap->FindValleyThreshold(nHistoBins,true);
	double saliencyMedianThr= saliencyThresholdFactor*saliencyMedian;
	double saliencyBkgThr= saliencyBkgThresholdFactor*saliencyMedian;
	//double saliencyThr= std::max(saliencyOtsuThr,saliencyMedianThr);
	double saliencyThr= std::max(std::min(saliencyOtsuThr,saliencyValleyThr),saliencyMedianThr);


	//Get signal binary map
	bool findNegativeExcess= false;
	bool findNestedSources= false;
	bool useLocalBackground= false; 
	bool mergeBelowSeed= false;
	double seedThr= saliencyThresholdFactor;//5;
	double mergeThr= 2.6;
	double peakThr= 3;
	Img* saliencySignificanceMap= saliencyMap->GetSignificanceMap(useLocalBackground);
	saliencyMap->ResetSources();
	saliencyMap->FindCompactSource(saliencyMap,saliencyThr,saliencyThr,minNPix,findNegativeExcess,findNestedSources,useLocalBackground,mergeBelowSeed,peakThr);
	//saliencyMap->FindCompactSource(saliencySignificanceMap,seedThr,mergeThr,minNPix,findNegativeExcess,findNestedSources,useLocalBackground,mergeBelowSeed,peakThr);
	Img* saliencyMap_signalBinarized= saliencyMap->GetSourceMap(true);

	saliencySignificanceMap->Delete();

	//Get bkg binary map
	Img* saliencyMap_bkgBinarized= saliencyMap->GetBinarized(saliencyBkgThr,0,1,true);

	cout<<"SLICUtils::TagRegions(): INFO: <saliency>="<<saliencyMedian<<" saliencyMedianThr="<<saliencyMedianThr<<" saliencyOtsuThr="<<saliencyOtsuThr<<" saliencyThr="<<saliencyThr<<" saliencyBkgThr="<<saliencyBkgThr<<endl;

	int nBkgReg= 0;
	int nSignalReg= 0;
	int nUntaggedReg= 0;

	for(int i=0;i<nRegions;i++){	
		int nPix= regions[i]->fNPix;
		int nBkg= 0;
		int nUntagged= 0;
		int nSignal= 0;
		if(nPix<=0) {
			regions[i]->fTag= Region::eUntagged;
			regions[i]->fSaliency= 0;
			continue;
		}
		
		double regionMean= regions[i]->fMean;
		double averageSaliency= 0;
		for(int j=0;j<nPix;j++){//loop on pixels inside region
			int pixId= (regions[i]->fPixelCollection)[j].id;
			double saliency= saliencyMap->GetBinContent(pixId);	
			averageSaliency+= saliency;
			 
			bool isSignal= (saliencyMap_signalBinarized && saliencyMap_signalBinarized->GetBinContent(pixId)>0);
			bool isBkg= (saliencyMap_bkgBinarized && saliencyMap_bkgBinarized->GetBinContent(pixId)>0);

			if(isBkg && !isSignal) nBkg++;
			else if(isSignal && !isBkg) nSignal++;
			else nUntagged++;
		}//end loop pixels in region
		averageSaliency/= (double)nPix;
		regions[i]->fSaliency= averageSaliency;

		//Tag using a majority criterion
		regions[i]->fTag= Region::eUntagged;
		if(nSignal>nBkg && nSignal>nUntagged && regionMean>mean) {
			regions[i]->fTag= Region::eSignalTag;
			nSignalReg++;
		}
		else if(nBkg>nSignal && nBkg>nUntagged) {
			regions[i]->fTag= Region::eBkgTag;
			nBkgReg++;
		}
		else if(nUntagged>nSignal && nUntagged>nBkg) {
			regions[i]->fTag= Region::eUntagged;
			nUntaggedReg++;
		}
	}//end loop regions

	cout<<"TagRegions_last(): INFO: (nS,nB,nU)=("<<nSignalReg<<","<<nBkgReg<<","<<nUntaggedReg<<")"<<endl;

	return 0;

}//close TagRegions_last()


int SLICUtils::TagRegions(Img* image,Img* saliencyMap,std::vector<Region*> regions,double saliencyThresholdFactor,double saliencyBkgThresholdFactor,double CL,bool includeSpatialPar,bool includeCurvPar,int knn){

	cout<<"SLICUtils::TagRegions(): INFO: Tag regions..."<<endl;
	if(!image ) return -1;

	int nRegions= (int)regions.size();
	if(nRegions<=0) {
		cerr<<"SLICUtils::TagRegions(): ERROR: No regions available, nothing to be tagged!"<<endl;
		return -1;
	}
	
	//## Compute saliency
	if(!saliencyMap){
		cerr<<"SLICUtils::TagRegions(): ERROR: Cannot get the saliency map"<<endl;
		return -1;
	}
	if(!saliencyMap->ComputeStats(true,false,true)<0){
		cerr<<"SLICUtils::TagRegions(): ERROR: Stats computing failed!"<<endl;
		if(saliencyMap) saliencyMap->Delete();
		return -1;
	}
	Img::StatsData* saliencyStats= saliencyMap->GetPixelStats();
	double saliencyMedian= saliencyStats->median;
	double saliencyMedianRMS= saliencyStats->medianRMS;
	double saliencyThr= saliencyMedian + saliencyThresholdFactor*saliencyMedianRMS;
	double saliencyBkgThr= saliencyMedian + saliencyBkgThresholdFactor*saliencyMedianRMS;
	//double saliencyThr= saliencyThresholdFactor*saliencyMedian;
	//double saliencyBkgThr= saliencyBkgThresholdFactor*saliencyMedian;
	cout<<"SLICUtils::TagRegions(): INFO: <saliency>="<<saliencyMedian<<" saliencyThr="<<saliencyThr<<" saliencyBkgThr="<<saliencyBkgThr<<endl;

	
	for(int i=0;i<nRegions;i++){	
		int nPix= regions[i]->fNPix;
		if(nPix<=0) regions[i]->fSaliency= 0;
		int binId= (regions[i]->fPixelCollection)[0].id;
		double Saliency= saliencyMap->GetBinContent(binId);
		regions[i]->fSaliency= Saliency;
		//(*saliency_data_matrix)(i,0)= Saliency;
	}//end loop regions
	
	
	//## Compute robust mean/covariance/Mahalanobis distance
	int nPars= 2;
	if(includeSpatialPar) nPars+= 2;
	if(includeCurvPar) nPars+= 2;
	TMatrixD* data_matrix= new TMatrixD(nRegions,nPars);
	for(int i=0;i<nRegions;i++){
		TVectorD regionPars= regions[i]->GetParamVector(includeSpatialPar,includeCurvPar);//only color params (no centroids!)
		for(int j=0;j<regionPars.GetNoElements();j++)	{
			(*data_matrix)(i,j)= regionPars(j);	
		}//end loop par dim
	}//end loop regions
		
	//## Perform outlier analysis and retrieve results
	OutlierDetector outlierDetector;
	if( outlierDetector.FindOutliers(data_matrix)<0 ) {
		cerr<<"SLICUtils::TagRegions(): ERROR: Outlier detection failed!"<<endl;
		if(data_matrix) data_matrix->Delete();
		return -1;
	}

	std::vector<int> data_flags= outlierDetector.GetDataFlags();
	std::vector<int> outliers_ids= outlierDetector.GetOutliers();
	std::vector<double> MHDs= outlierDetector.GetDistances();
	TMatrixD* robustCov= outlierDetector.GetRobustCov();
	TVectorD* robustMean= outlierDetector.GetRobustMean();
	double cutoffValue= outlierDetector.GetCutoff();
		
	//## Tag significative regions according to Mahalanobis distance
	//## CL: 0.975 (conventional) 0.997 (3sigma) 0.954 (2sigma), 0.683 (1sigma)
	//## Tag salient regions
	double chi2quantile= ROOT::Math::chisquared_quantile(CL,nPars);
	double outlierCut= chi2quantile;
	double saliency_chi2quantile= ROOT::Math::chisquared_quantile(CL,1);
	double saliency_outlierCut= saliency_chi2quantile;
	//double saliencyThr= saliencyThresholdFactor*saliencyAvg;
	//cout<<"SLICUtils::TagRegions(): INFO: <saliency>="<<saliencyAvg<<" saliencyThr="<<saliencyThr<<endl;

	for(int i=0;i<nRegions;i++){
		int regionId= regions[i]->fId;
		int outlierFlag= data_flags[i];	
		double MHD= MHDs[i];
		//double saliencyMHD= saliency_MHDs[i];

		bool isOutlier= (MHD>outlierCut);
		bool isPositiveOutlier= ( (*data_matrix)(i,0)>(*robustMean)(0) );
		double S= regions[i]->fSaliency;
		bool isSalient= (S>saliencyThr);	
		bool isBelowBkgSalientThr= (S<saliencyBkgThr);
		regions[i]->fMahalanobisDistance= MHD;
		
		if(isSalient && isPositiveOutlier) regions[i]->fIsSalient= true;
		else regions[i]->fIsSalient= false;
		if(isOutlier && isPositiveOutlier) regions[i]->fIsSignificative= true;
		else regions[i]->fIsSignificative= false;

		//Tag saliency
		if( isPositiveOutlier && (S>saliencyThr) ) regions[i]->fTag= Region::eSignalTag;
		else if(S<saliencyBkgThr) regions[i]->fTag= Region::eBkgTag;
		else regions[i]->fTag= Region::eUntagged;

		//cout<<"SLICUtils::TagRegions(): INFO: Region id="<<regionId<<": MHD="<<MHD<<" (cut="<<outlierCut<<",cutoffValue="<<cutoffValue<<") isOutlier?"<<isOutlier<<" isPositiveOutlier? "<<isPositiveOutlier<<" outlierFlag="<<outlierFlag<<" isSalient?"<<isSalient<<endl;		
	}//end loop regions

	//Clear-up
	if(data_matrix) data_matrix->Delete();

	return 0;

}//close TagRegions()


Img* SLICUtils::GetSaliencyMap(Img* image,TMatrixD* dissMatrix,std::vector<Region*> regions,int knn){

	//## Compute saliency
	cout<<"SLICUtils::GetSaliencyMap(): INFO: Applying saliency thresholding..."<<endl;
	int nRegions= (int)regions.size();
	if(nRegions<=0 || !image || !dissMatrix) return 0;
	
	//## Compute saliency
	TString imgName= Form("%s_saliency",image->GetName());
	Img* saliencyImg= (Img*)image->Clone(imgName);
	saliencyImg->SetNameTitle(imgName,imgName);
	saliencyImg->Reset();

	std::vector<double> regionSaliency;

	if(knn<nRegions){ 
	
		for(int i=0;i<nRegions;i++){
			std::vector<double> dissimilarityList;
			for(int j=0;j<nRegions;j++){
				double diss= (*dissMatrix)(i,j);
				dissimilarityList.push_back(diss);
			}//end loop regions

			std::vector<double> sorted;
			std::vector<size_t> sort_index;//sorting index
			Utils::sort( dissimilarityList,sorted,sort_index);		
					
			//Compute saliency over k-th neighbors
			double sum= 0;
			double Saliency= 0;
			for(int k=0;k<knn;k++){
				size_t index= sort_index[k];
				if(index==i) continue;
				double diss= (*dissMatrix)(i,index);
				sum+= diss;
			}
			if(sum!=0) Saliency= 1.-exp(-sum/(double)(knn-1));
			regionSaliency.push_back(Saliency);		
		}//end loop regions

	}//close if
	else{
		for(int i=0;i<nRegions;i++){
			double Saliency= 0;
			double sum= 0;
			for(int j=0;j<nRegions;j++){
				if(i==j) continue;
				double diss= (*dissMatrix)(i,j);
				sum+= diss;
			}//end loop regions
			if(sum!=0) Saliency= 1.-exp(-sum/(double)(nRegions-1));
			regionSaliency.push_back(Saliency);
		}//end loop regions
	}//close else


	//## Fill image
	for(int i=0;i<nRegions;i++){
		double Saliency= regionSaliency[i];
		int nPixelsInRegion= (int)regions[i]->fNPix;
		for(int j=0;j<nPixelsInRegion;j++){//loop on pixels inside region
			int thisPixelId= (regions[i]->fPixelCollection)[j].id;
			saliencyImg->SetBinContent(thisPixelId,Saliency);
		}//end loop pixels in region	
	}//end loop regions

	return saliencyImg;

}//close GetSaliencyMap()


Img* SLICUtils::GetSaliencyMap(Img* image,std::vector<Region*> regions,int knn){

	//## Compute saliency
	cout<<"SLICUtils::GetSaliencyMap(): INFO: Applying saliency thresholding..."<<endl;
	int nRegions= (int)regions.size();
	if(nRegions<=0 || !image) return 0;
	
	//## Compute saliency
	TString imgName= Form("%s_saliency",image->GetName());
	Img* saliencyImg= (Img*)image->Clone(imgName);
	saliencyImg->SetNameTitle(imgName,imgName);
	saliencyImg->Reset();

	std::vector<double> regionSaliency;

	if(knn<nRegions){ 
		for(int i=0;i<nRegions;i++){
			double Saliency= 0;
			std::vector<double> distList;
			for(int j=0;j<nRegions;j++){
				std::pair<double,double> dists= regions[i]->GetDistance(regions[j],false,true);
				double dist_appearance= dists.first;
				distList.push_back(dist_appearance);
			}//end loop regions

			std::vector<double> sorted;
			std::vector<size_t> sort_index;//sorting index
			Utils::sort( distList,sorted,sort_index);		
			
			//Compute saliency over k-th neighbors
			double sum= 0;
			for(int k=0;k<knn;k++){
				size_t index= sort_index[k];
				double diss= regions[i]->GetDissimilarity(regions[index],false,true);
				sum+= diss;
			}	
			if(sum!=0) Saliency= 1.-exp(-sum/(double)(knn));
			regionSaliency.push_back(Saliency);		
		}//end loop regions

	}//close if
	else{
		for(int i=0;i<nRegions;i++){
			double Saliency= 0;
			double sum= 0;
			for(int j=0;j<nRegions;j++){
				if(i==j) continue;
				double diss= regions[i]->GetDissimilarity(regions[j],false,true);
				sum+= diss;
			}//end loop regions
			if(sum!=0) Saliency= 1.-exp(-sum/(double)(nRegions-1));
			regionSaliency.push_back(Saliency);
		}//end loop regions
	}//close else


	//## Fill image
	for(int i=0;i<nRegions;i++){
		double Saliency= regionSaliency[i];
		int nPixelsInRegion= (int)regions[i]->fNPix;
		for(int j=0;j<nPixelsInRegion;j++){//loop on pixels inside region
			int thisPixelId= (regions[i]->fPixelCollection)[j].id;
			saliencyImg->SetBinContent(thisPixelId,Saliency);
		}//end loop pixels in region	
	}//end loop regions

	return saliencyImg;

}//close GetSaliencyMap()


SLICUtils::SLICContourData* SLICUtils::ComputeBoundaryContours(Img* image,std::vector< std::vector<long int> > labels, std::vector<Region*> regions) {
  	
	//## Check input
	if(!image) return 0;
	int nRegions= (int)regions.size();
	if(nRegions<=0) return 0;

	//## Init data
	cout<<"SLICUtils::GetClusterContours(): INFO: Computing contours from NR="<<regions.size()<<" regions..."<<endl;
	int Nx= image->GetNbinsX();
	int Ny= image->GetNbinsY();
	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};
		
	int nContourPts= 0;
	std::map<int,int> mapping;
	std::vector< std::vector<bool> > isNeighborLinkTaken;
	SLICContourData* contourData= new SLICContourData;
	(contourData->contour)= new Contour;
	
	for(int k=0;k<nRegions;k++){
		int regionId= regions[k]->fId;
		mapping.insert( std::pair<int,int>(regionId,k) );
		contourData->connectedRegionIds.push_back( std::vector<int>() );
		contourData->boundaryData.push_back( SLICBoundaryPixMap() );
		isNeighborLinkTaken.push_back( std::vector<bool>() );
		for(int s=0;s<nRegions;s++) isNeighborLinkTaken[k].push_back(false);
	}//end loop regions 

  //## Loop over all the pixels
  for (int i=0; i<Nx; i++) {
  	for (int j=0; j<Ny; j++) {
			int regionId= labels[i][j];
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
					int regionId_neighbor= labels[x][y];
       	
        	if (regionId!=regionId_neighbor && regionId_neighbor>=0) {	
          	nr_p++;
						connectedRegionCounterMap[regionId_neighbor]++;
						int gBin_neighbor= image->GetBin(x+1,y+1);
						connectedRegionBoundaryPixelListMap[regionId_neighbor];
						connectedRegionBoundaryPixelListMap[regionId_neighbor].push_back(gBin_neighbor);
          }
        }//close if
      }//end loop neighbours
            
      // Add the pixel to the contour list of corresponding region
     	if (nr_p>= 2) {
				int gBin= image->GetBin(i+1,j+1);
				double binX= image->GetXaxis()->GetBinCenter(i+1);
				double binY= image->GetYaxis()->GetBinCenter(j+1);	
				(contourData->contour)->AddPoint(cv::Point2f(binX,binY));
      	
				//Add neighbor ids and boundary pixel ids
				std::map<int,std::vector<int>>::iterator counterListIt = connectedRegionBoundaryPixelListMap.begin();
				for (counterListIt=connectedRegionBoundaryPixelListMap.begin(); counterListIt!=connectedRegionBoundaryPixelListMap.end(); ++counterListIt){
					int neighborId= counterListIt->first;
					std::vector<int> neighborPixIds= counterListIt->second;
					int counts= (int)neighborPixIds.size();	
					int neighborIndex= mapping[neighborId];

					if(counts>=2) {
						
						//std::vector<int>::iterator it = std::find (contourData->connectedRegionIds[regionIndex].begin(), contourData->connectedRegionIds[regionIndex].end(), neighborId);	
						std::vector<int>::iterator it = std::find (contourData->connectedRegionIds[regionIndex].begin(), contourData->connectedRegionIds[regionIndex].end(), neighborIndex);
						
						if( contourData->connectedRegionIds[regionIndex].empty() || it==contourData->connectedRegionIds[regionIndex].end() ) {
							//contourData->connectedRegionIds[regionIndex].push_back(neighborId);
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


SLICUtils::SLICSimilarityData* SLICUtils::ComputeDissimilarityMatrix(Img* edgeImage,SLICUtils::SLICContourData* contourData, std::vector<Region*> regions,double beta,bool includeSpatialDist){

	//## Check input data
	int nRegions= (int)regions.size();
	if(nRegions<=0) return 0;
	if(!contourData) return 0;
	if(!edgeImage) return 0;
	Img::StatsData* edgeImgStats= edgeImage->GetPixelStats();
	if(!edgeImgStats){
		cerr<<"SLICUtils::ComputeDissimilarityMatrix(): ERROR: No stats for edge image...return!"<<endl;
		return 0;
	}
	double EdgeImgMin= edgeImgStats->min;
	double EdgeImgMax= edgeImgStats->max;
	double Emin_norm= 0;
	double Emax_norm= 1;
	double Dmin_norm= 0;
	double Dmax_norm= 1;
	double DSmin_norm= 0;
	double DSmax_norm= 1;
	cout<<"SLICUtils::ComputeDissimilarityMatrix(): INFO: Init matrix (nRegions="<<nRegions<<")"<<endl;

	TMatrixD* DissimilarityMatrix= new TMatrixD(nRegions,nRegions);
	DissimilarityMatrix->Zero();//fill with 0s

	TMatrixD* SaliencyDissimilarityMatrix= new TMatrixD(nRegions,nRegions);
	SaliencyDissimilarityMatrix->Zero();//fill with 0s

	double sigmaP= 0.25;
	double sigmaC= 20/255.;
	TMatrixD* SaliencySpatialWeightMatrix= new TMatrixD(nRegions,nRegions);
	SaliencySpatialWeightMatrix->Zero();//fill with 0s
	TMatrixD* SaliencyColorWeightMatrix= new TMatrixD(nRegions,nRegions);
	SaliencyColorWeightMatrix->Zero();//fill with 0s
	
	for(int i=0;i<nRegions;i++) {
		(*SaliencySpatialWeightMatrix)(i,i)= 1;
 		(*SaliencyColorWeightMatrix)(i,i)= 1;
 	}

	TMatrixD* EdgenessMatrix= new TMatrixD(nRegions,nRegions);
	EdgenessMatrix->Zero();
	(*EdgenessMatrix)+= EdgeImgMax;//fill with max edges
	for(int i=0;i<EdgenessMatrix->GetNrows();i++) (*EdgenessMatrix)(i,i)= 0;//fill diagonal with 0s

	TMatrixD* TotDissimilarityMatrix= new TMatrixD(nRegions,nRegions);
	TotDissimilarityMatrix->Zero();//fill with 0s

	TMatrixD* SpatialDissimilarityMatrix= new TMatrixD(nRegions,nRegions);
	SpatialDissimilarityMatrix->Zero();
	
	std::vector<SLICBoundaryPixMap> boundaryData= contourData->boundaryData;
	SLICConnectedRegions connectedRegionIds= contourData->connectedRegionIds;

	
	//## Compute dissimilarity and edgeness matrix	
	cout<<"SLICUtils::ComputeDissimilarityMatrix(): INFO: Compute region dissimilarity and edgeness (nRegions="<<nRegions<<")"<<endl;

	std::vector<double> EdgenessList;
	std::vector<double> DissList;
	std::vector<double> SpatialDissList;

	for(int i=0; i<nRegions-1; i++) {
		int regionId= regions[i]->fId;

		for(int j=i+1; j<nRegions; j++) {				
			int neighborId= regions[j]->fId;

			//Compute appearance term
			std::pair<double,double> dists= regions[i]->GetDistance(regions[j],false,true);
			double dist_color= dists.first;
			double dist_space= dists.second;
			double Diss= dist_color;
			if(includeSpatialDist) Diss+= dist_space;
			(*DissimilarityMatrix)(i,j)= Diss;
			(*DissimilarityMatrix)(j,i)= Diss;

			(*SpatialDissimilarityMatrix)(i,j)= dist_space;
			(*SpatialDissimilarityMatrix)(j,i)= dist_space;

			DissList.push_back(Diss);
			SpatialDissList.push_back(dist_space);
			
			//Compute edge term
			//Find if region j-th is among neighbors of i-th
			SLICBoundaryPixMapIterator it= boundaryData[i].find(neighborId);
			
			if(it!=boundaryData[i].end()){//region is among neighbors, compute edgeness terms
				std::vector<int> sharedPixelIds= (boundaryData[i])[neighborId];
				int nBoundaryPixels= (int)sharedPixelIds.size();
				double E= 0;
				for(int t=0;t<nBoundaryPixels;t++) {
					int gBin= sharedPixelIds[t];
					double S_edge= edgeImage->GetBinContent(gBin);
					double S_edge_norm= Emin_norm + (Emax_norm-Emin_norm)*(S_edge-EdgeImgMin)/(EdgeImgMax-EdgeImgMin);
					E+= S_edge;
				}	
				if(nBoundaryPixels>0) E/= (double)nBoundaryPixels; 
						
				std::vector<int> sharedPixelIds_neighbors= (boundaryData[j])[regionId];
				int nBoundaryPixels_neighbors= (int)sharedPixelIds_neighbors.size();
				double E_neighbor= 0;
				for(int t=0;t<nBoundaryPixels_neighbors;t++) {
					int gBin= sharedPixelIds_neighbors[t];
					double S_edge= edgeImage->GetBinContent(gBin);
					double S_edge_norm= Emin_norm + (Emax_norm-Emin_norm)*(S_edge-EdgeImgMin)/(EdgeImgMax-EdgeImgMin);
					E_neighbor+= S_edge;
				}	
				if(nBoundaryPixels_neighbors>0) E_neighbor/= (double)nBoundaryPixels_neighbors; 

				double Etot= 0.5*(E+E_neighbor);
				(*EdgenessMatrix)(i,j)= Etot;
				(*EdgenessMatrix)(j,i)= Etot;
				EdgenessList.push_back(Etot);
			}//close else
			else{//regions are not direct neighbors
				(*EdgenessMatrix)(i,j)= EdgeImgMax;//set to maximum edge term
				(*EdgenessMatrix)(j,i)= EdgeImgMax;
			}
		}//end loop regions
	}//end loop regions

	std::sort(DissList.begin(),DissList.end());
	double Dmin= 0;
	double Dmax= DissList[DissList.size()-1];
	std::sort(EdgenessList.begin(),EdgenessList.end());
	double Emin= 0;
	double Emax= EdgeImgMax;
	std::sort(SpatialDissList.begin(),SpatialDissList.end());
	double DSmin= 0;
	double DSmax= SpatialDissList[SpatialDissList.size()-1];
	
	cout<<"SLICUtils::ComputeDissimilarityMatrix(): INFO: Normalizing dissimilarity & edgeness in ["<<Dmin_norm<<","<<Dmax_norm<<"], Dmin/Dmax="<<Dmin<<"/"<<Dmax<<" Emin/Emax="<<Emin<<"/"<<Emax<<" DSmin/DSmax="<<DSmin<<"/"<<DSmax<<endl;

	//## Normalize the dissimilarities and compute tot diss
	for(int i=0; i<nRegions-1; i++) {
		int regionId= regions[i]->fId;

		for(int j=i+1; j<nRegions; j++) {
			int neighborId= regions[j]->fId;
		
			double D= (*DissimilarityMatrix)(i,j);
			double D_norm= Dmin_norm + (Dmax_norm-Dmin_norm)*(D-Dmin)/(Dmax-Dmin);
			(*DissimilarityMatrix)(i,j)= D_norm;
			(*DissimilarityMatrix)(j,i)= D_norm;
			
			double E= (*EdgenessMatrix)(i,j);
			double E_norm= Emin_norm + (Emax_norm-Emin_norm)*(E-Emin)/(Emax-Emin);
			//cout<<"(i,j)=("<<i<<","<<j<<") E="<<E<<" E_norm="<<E_norm<<" Emin/Emax="<<Emin<<"/"<<Emax<<endl;
			(*EdgenessMatrix)(i,j)= E_norm;
			(*EdgenessMatrix)(j,i)= E_norm;

			double DS= (*SpatialDissimilarityMatrix)(i,j);
			double DS_norm= DSmin_norm + (DSmax_norm-DSmin_norm)*(DS-DSmin)/(DSmax-DSmin);
			(*SpatialDissimilarityMatrix)(i,j)= DS_norm;
			(*SpatialDissimilarityMatrix)(j,i)= DS_norm;
			
			double Dtot= (1-beta)*D_norm + beta*E_norm;
			//double Dtot= D_norm + beta*E_norm;			
			//Dtot+= SMALL_NUMBER;//to avoid dividing by zero!
			(*TotDissimilarityMatrix)(i,j)= Dtot;
			(*TotDissimilarityMatrix)(j,i)= Dtot;

			double Dsal= Dtot/(1.+DS_norm);
			(*SaliencyDissimilarityMatrix)(i,j)= Dsal;
			(*SaliencyDissimilarityMatrix)(j,i)= Dsal;

			double w= exp(-DS_norm/(2.*sigmaP*sigmaP));
			(*SaliencySpatialWeightMatrix)(i,j)= w;
			(*SaliencySpatialWeightMatrix)(j,i)= w;
	
			double wc= exp(-D_norm/(2.*sigmaC*sigmaC));
			(*SaliencyColorWeightMatrix)(i,j)= wc;
			(*SaliencyColorWeightMatrix)(j,i)= wc;
						
		}//end loop regions
	}//end loop regions

	
	//Normalize similarity spatial W by rows so that sum_j(Wij)=1	
	//## Compute D & U
	std::vector<double> DVector;
	std::vector<double> UVector;
		
	for(int i=0;i<SaliencySpatialWeightMatrix->GetNrows();i++){
		double sum= 0;	
		double sumc= 0;	
	
		//Normalize weights to unit sum
		for(int j=0;j<SaliencySpatialWeightMatrix->GetNcols();j++) {
			sum+= (*SaliencySpatialWeightMatrix)(i,j);
			sumc+= (*SaliencyColorWeightMatrix)(i,j);
		}
		if(sum==0 || sumc==0) {
			cerr<<"SLICUtils::ComputeDissimilarityMatrix(): ERROR: Null sum in saliency weights!"<<endl;
			return 0;
		}

		//Compute weighted color position
		TVector2 mu_i= TVector2(0,0); 
		for(int j=0;j<SaliencySpatialWeightMatrix->GetNcols();j++) {
			(*SaliencySpatialWeightMatrix)(i,j)/= sum;
			(*SaliencyColorWeightMatrix)(i,j)/= sum;

			TVector2 Pj= TVector2(regions[j]->fX0,regions[j]->fY0);
			double wc= (*SaliencyColorWeightMatrix)(i,j);
			mu_i+= wc*Pj;
		}

		//Compute U & D terms
		double Pdist_min= 0;
		double Pdist_max= sqrt(pow(regions[i]->fImageMaxX-regions[i]->fImageMinX,2) + pow(regions[i]->fImageMaxY-regions[i]->fImageMinY,2));
		double salD= 0;
		double salU= 0;
		double saliency= 0;
		for(int j=0;j<SaliencySpatialWeightMatrix->GetNcols();j++) {
			double wc= (*SaliencyColorWeightMatrix)(i,j);
			TVector2 Pj= TVector2(regions[j]->fX0,regions[j]->fY0);
			double Pdist= (Pj-mu_i).Mod();
			double Pdist_norm= DSmin_norm + (DSmax_norm-DSmin_norm)*(Pdist-Pdist_min)/(Pdist_max-Pdist_min);
			double Pdist_norm_weighted= Pdist_norm*wc;
			salD+= Pdist_norm_weighted;
			
			double d= (*DissimilarityMatrix)(i,j);//already normalized
			double ws= (*SaliencySpatialWeightMatrix)(i,j);//already normalized to unit sum
			double u= d*ws;
			salU+= u;
		}
		DVector.push_back(salD);	
		UVector.push_back(salU);	

	}//end loop rows

	//Normalize D and U to [0,1]
	std::vector<double> saliencyU;
	saliencyU.assign(UVector.begin(),UVector.end());	
	std::vector<double> saliencyD;
	saliencyD.assign(DVector.begin(),DVector.end());
	std::vector<double> saliencyList;

	std::sort(DVector.begin(),DVector.end());
	std::sort(UVector.begin(),UVector.end());
	double minSalU= UVector[0];
	double maxSalU= UVector[UVector.size()-1];
	double minSalD= DVector[0];
	double maxSalD= DVector[DVector.size()-1];
	double normMin= 0;
	double normMax= 1;
	double saliencyKFactor= 6;
	
	for(unsigned int i=0;i<saliencyU.size();i++){
		double salU= saliencyU[i];
		double salU_norm= normMin + (normMax-normMin)*(salU-minSalU)/(maxSalU-minSalU);
		double salD= saliencyD[i];
		double salD_norm= normMin + (normMax-normMin)*(salD-minSalD)/(maxSalD-minSalD);
		salU_norm+= SMALL_NUMBER;
		salD_norm+= SMALL_NUMBER;
		saliencyU[i]= salU_norm;
		saliencyD[i]= salD_norm;
		double Saliency= salU_norm*exp(-saliencyKFactor*salD_norm);
		saliencyList.push_back(Saliency);
	}//end loop regions



	/*
	//## Find adjacency matrix (Consider only 1-st and 2-nd neighbors)
	TMatrixD AdjMatrix(nRegions,nRegions);
	AdjMatrix.Zero();
	for(int i=0;i<AdjMatrix.GetNrows();i++) AdjMatrix(i,i)= 1;
 
	std::vector< std::vector<int> > neighborIndexList;
	for(int i=0; i<nRegions; i++) {
		int regionId= regions[i]->fId;

		neighborIndexList.push_back( std::vector<int>() );		
		
		//Fill list of 1-st neighbors	
		neighborIndexList[i].assign(connectedRegionIds[i].begin(),connectedRegionIds[i].end());

		//Fill list of 2-nd neighbors			
		for(unsigned int j=0;j<connectedRegionIds[i].size();j++){//loop over 1st neighbors
			int neighborIndex= connectedRegionIds[i][j];

			for(unsigned int t=0;t<connectedRegionIds[neighborIndex].size();t++){//loop over 2nd neighbors
				int neighborIndex_2nd= connectedRegionIds[neighborIndex][t];
				std::vector<int>::iterator it= std::find(neighborIndexList[i].begin(),neighborIndexList[i].end(),neighborIndex_2nd);
		
				if( it==neighborIndexList[i].end() || neighborIndexList[i].size()==0 ){
					neighborIndexList[i].push_back(neighborIndex_2nd);
				}
			}//end loop 2nd-neighbors
		}//end loop 1st-neighbors

		//Loop over all neighbors and fill adj matrix
		for(unsigned int j=0;j<neighborIndexList[i].size();j++){
			int neighborIndex= neighborIndexList[i][j];
			AdjMatrix(i,neighborIndex)= 1;
			AdjMatrix(neighborIndex,i)= 1;		
		}//end loop all neighbors
	}//end loop regions


	//## Set non-neighbors elements to max value
	double DtotMax= TotDissimilarityMatrix->Max();
	cout<<"DtotMax="<<DtotMax<<endl;
	for(int i=0; i<nRegions-1; i++) {
		int regionId= regions[i]->fId;

		for(int j=i+1; j<nRegions; j++) {
			int neighborId= regions[j]->fId;
			if(AdjMatrix(i,j)==0) {
				(*TotDissimilarityMatrix)(i,j)= DtotMax;
				(*TotDissimilarityMatrix)(j,i)= DtotMax;
			}
		}//end loop regions
	}//end loop regions
	*/

	/*
	cout<<"==Dissimilarity Matrix  =="<<endl;
	DissimilarityMatrix->Print();
	cout<<"======================="<<endl;

	cout<<"== Edgeness Matrix =="<<endl;
	EdgenessMatrix->Print();
	cout<<"======================="<<endl;

	//cout<<"== Adjacency Matrix =="<<endl;
	//AdjMatrix.Print();
	//cout<<"======================="<<endl;

	cout<<"== Tot Dissimilarity Matrix =="<<endl;
	TotDissimilarityMatrix->Print();
	cout<<"======================="<<endl;
	
	cout<<"== Saliency Dissimilarity Matrix  =="<<endl;
	SaliencyDissimilarityMatrix->Print();
	cout<<"======================="<<endl;
	*/

	//clear up
	TotDissimilarityMatrix->Delete();
	//DissimilarityMatrix->Delete();
	EdgenessMatrix->Delete();
	SpatialDissimilarityMatrix->Delete();
	SaliencySpatialWeightMatrix->Delete();
	SaliencyColorWeightMatrix->Delete();

	//return data
	SLICUtils::SLICSimilarityData* similarityData= new SLICUtils::SLICSimilarityData;
	similarityData->DissimilarityMatrix= DissimilarityMatrix;
	similarityData->SaliencyDissimilarityMatrix= SaliencyDissimilarityMatrix;
	similarityData->saliencyList= saliencyList;
	
	return similarityData;

}//close ComputeDissimilarityMatrix()




SLICUtils::SLICSimilarityData* SLICUtils::ComputeRegionSimilarity(Img* edgeImage,SLICUtils::SLICContourData* contourData, std::vector<Region*> regions,double beta,bool includeSpatialDist,int mergedTag){

	//## Check input data
	int nRegions= (int)regions.size();
	if(nRegions<=0) return 0;
	if(!contourData) return 0;
	if(!edgeImage) return 0;
	Img::StatsData* edgeImgStats= edgeImage->GetPixelStats();
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
			//DissMatrix[i].push_back(0);
			DissMatrix[i].push_back(1.e+99);
			(*EdgenessMatrix)(i,j)= Emax_norm;
		}
	}

	

	/*
	//## Compute dissimilarity matrix	
	for(int i=0; i<nRegions-1; i++) {
		int regionId= regions[i]->fId;

		for(int j=i+1; j<nRegions; j++) {
			int neighborId= regions[j]->fId;
			std::pair<double,double> dists= regions[i]->GetAsymmDistance(regions[j],false,true,includeSpatialDist);
				
			//Appearance distance
			double Diss_appearance= dists.first;//appeareance dist
			double DissNeighbor_appearance= dists.second;//appeareance dist

			//Spatial distance (TO BE COMPUTED!!!)
			double Diss_spatial= 0;//spatial dist
			double DissNeighbor_spatial= 0;//spatial dist
			
			//TOT distance
			double Diss= Diss_appearance;//total diss
			double DissNeighbor= DissNeighbor_appearance;//total diss

			if(includeSpatialDist) {
				Diss+= Diss_spatial;
				DissNeighbor+= DissNeighbor_spatial;
			}
			(*(SimilarityData->DissimilarityMatrix))(i,j)= Diss;
			(*(SimilarityData->DissimilarityMatrix))(j,i)= DissNeighbor;

			//Compute edge term
			//Find if region j-th is among neighbors of i-th
			SLICBoundaryPixMapIterator it= boundaryData[i].find(neighborId);
			double E= Emin_norm ;
			double E_neighbor= Emin_norm;
			if(it==boundaryData[i].end()){//not found among neighbors, set edgeness to maximum value allowed
				E= Emax_norm;
				E_neighbor= Emax_norm;
			}
			else{//region is among neighbors, compute edgeness terms
				std::vector<int> sharedPixelIds= (boundaryData[i])[neighborId];
				int nBoundaryPixels= (int)sharedPixelIds.size();
				for(int t=0;t<nBoundaryPixels;t++) {
					int gBin= sharedPixelIds[t];
					double S_edge= edgeImage->GetBinContent(gBin);
					double S_edge_norm= Emin_norm + (Emax_norm-Emin_norm)*(S_edge-EdgeImgMin)/(EdgeImgMax-EdgeImgMin);
					E+= S_edge_norm;
				}	
				if(nBoundaryPixels>0) E/= (double)nBoundaryPixels; 
						
				std::vector<int> sharedPixelIds_neighbors= (boundaryData[j])[regionId];
				int nBoundaryPixels_neighbors= (int)sharedPixelIds_neighbors.size();
				for(int t=0;t<nBoundaryPixels_neighbors;t++) {
					int gBin= sharedPixelIds_neighbors[t];
					double S_edge= edgeImage->GetBinContent(gBin);
					double S_edge_norm= Emin_norm + (Emax_norm-Emin_norm)*(S_edge-EdgeImgMin)/(EdgeImgMax-EdgeImgMin);
					E_neighbor+= S_edge_norm;
				}	
				if(nBoundaryPixels_neighbors>0) E_neighbor/= (double)nBoundaryPixels_neighbors; 
			}//close else

			(*EdgenessMatrix)(i,j)= E;
			(*EdgenessMatrix)(j,i)= E_neighbor;
		}//end loop regions
	}//end loop regions
	*/

	//## Consider only 1-st and 2-nd neighbors
	std::vector< std::vector<int> > neighborIndexList_2nd;
	std::vector<double> EdgenessList;
	std::vector<double> DissList;

	for(int i=0; i<nRegions; i++) {
		int regionId= regions[i]->fId;

		//Fill list of 2-nd neighbors	
		neighborIndexList_2nd.push_back( std::vector<int>() );		
		
		for(unsigned int j=0;j<connectedRegionIds[i].size();j++){//loop over 1st neighbors
			int neighborIndex= connectedRegionIds[i][j];

			for(unsigned int t=0;t<connectedRegionIds[neighborIndex].size();t++){//loop over 2nd neighbors
				int neighborIndex_2nd= connectedRegionIds[neighborIndex][t];
				int neighborTag_2nd= regions[neighborIndex_2nd]->fTag;
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
			int neighborId= regions[neighborIndex]->fId;
			int neighborTag= regions[neighborIndex]->fTag;	

			//Set neighbor id
			(*(SimilarityData->NeighborMatrix))(i,neighborIndex)= 1;
			(*(SimilarityData->NeighborMatrix))(neighborIndex,i)= 1;
			
			if(mergedTag!=-1 && neighborTag!=mergedTag) continue;//skip neighbors if tag is different from the requested

			//Compute dissimilarity
			std::pair<double,double> dists= regions[i]->GetAsymmDistance(regions[neighborIndex],false,true,includeSpatialDist);
			double Diss= dists.first;
			double DissNeighbor= dists.second;
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
					double S_edge= edgeImage->GetBinContent(gBin);
					double S_edge_norm= Emin_norm + (Emax_norm-Emin_norm)*(S_edge-EdgeImgMin)/(EdgeImgMax-EdgeImgMin);
					//E+= S_edge_norm;
					E+= S_edge;
				}	
				if(nBoundaryPixels>0) E/= (double)nBoundaryPixels; 
						
				std::vector<int> sharedPixelIds_neighbors= (boundaryData[neighborIndex])[regionId];
				int nBoundaryPixels_neighbors= (int)sharedPixelIds_neighbors.size();
				double E_neighbor= 0;
				for(int t=0;t<nBoundaryPixels_neighbors;t++) {
					int gBin= sharedPixelIds_neighbors[t];
					double S_edge= edgeImage->GetBinContent(gBin);
					double S_edge_norm= Emin_norm + (Emax_norm-Emin_norm)*(S_edge-EdgeImgMin)/(EdgeImgMax-EdgeImgMin);
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
			int neighborId_2nd= regions[neighborIndex_2nd]->fId;
			int neighborTag_2nd= regions[neighborIndex_2nd]->fTag;	

			//Set neighbor id
			(*(SimilarityData->NeighborMatrix))(i,neighborIndex_2nd)= 2;
			(*(SimilarityData->NeighborMatrix))(neighborIndex_2nd,i)= 2;
			
			if(mergedTag!=-1 && neighborTag_2nd!=mergedTag) continue;//skip neighbors if tag is different from the requested

			//Compute dissimilarity
			std::pair<double,double> dists= regions[i]->GetAsymmDistance(regions[neighborIndex_2nd],false,true,includeSpatialDist);
			double Diss= dists.first;
			double DissNeighbor= dists.second;
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
	

	for(int i=0; i<nRegions; i++) {
		int regionId= regions[i]->fId;

		//Fill 1-st neighbors
		//cout<<"SLICUtils::ComputeRegionSimilarity(): INFO: Region no. "<<i<<": filling 1st neighbor info..."<<endl;
		for(unsigned int j=0;j<connectedRegionIds[i].size();j++){
			int neighborIndex= connectedRegionIds[i][j];
			int neighborId= regions[neighborIndex]->fId;
			int neighborTag= regions[neighborIndex]->fTag;
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

			//cout<<"Neighbor "<<neighborIndex<<": D="<<D<<" Dtot="<<Dtot<<endl;
		}//end loop neighbors

		//Fill 2-nd neighbors
		for(unsigned int j=0;j<neighborIndexList_2nd[i].size();j++){
			int neighborIndex_2nd= neighborIndexList_2nd[i][j];
			int neighborId_2nd= regions[neighborIndex_2nd]->fId;
			int neighborTag_2nd= regions[neighborIndex_2nd]->fTag;			
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

			//cout<<"Neighbor "<<neighborIndex_2nd<<": D="<<D<<" Dtot="<<Dtot<<endl;
		}//end loop 2nd neighbors

	}//end loop regions
	

	/*
	for(int i=0;i<(SimilarityData->DissimilarityMatrix)->GetNrows();i++){
		for(int j=0;j<(SimilarityData->DissimilarityMatrix)->GetNcols();j++){
			double D= (*(SimilarityData->DissimilarityMatrix))(i,j);
			double D_norm= Dmin_norm + (Dmax_norm-Dmin_norm)*(D-Dmin)/(Dmax-Dmin);
			D2sum+= (D-Dmean)*(D-Dmean);

			double E= (*EdgenessMatrix)(i,j);
			double E_norm= Emin_norm + (Emax_norm-Emin_norm)*(E-Emin)/(Emax-Emin);
			
			//double Dtot= (1-beta)*D_norm + beta*E_norm;
			double Dtot= D_norm + beta*E_norm;			
			Dtot+= SMALL_NUMBER;//to avoid dividing by zero! 

			if(i==j){
				(*(SimilarityData->DissimilarityMatrix))(i,j)= 0;
				(*(SimilarityData->AdjacencyMatrix))(i,j)= 0;
				(*(SimilarityData->AbsDissimilarityMatrix))(i,j)= 0;
				(*EdgenessMatrix)(i,j)= 0;
			}
			else{
	
				if(D!=0){
					(*(SimilarityData->DissimilarityMatrix))(i,j)= Dtot;	
					(*(SimilarityData->AdjacencyMatrix))(i,j)= 1./Dtot;
					(*(SimilarityData->AbsDissimilarityMatrix))(i,j)= D;
					(*EdgenessMatrix)(i,j)= E_norm;
					AbsDissList.push_back( (*(SimilarityData->DissimilarityMatrix))(i,j) );
				}
				else{
					(*(SimilarityData->DissimilarityMatrix))(i,j)= 0;	
					(*(SimilarityData->AdjacencyMatrix))(i,j)= 0;
					(*(SimilarityData->AbsDissimilarityMatrix))(i,j)= D;
					(*EdgenessMatrix)(i,j)= Emax_norm;
				}
			}
			DissMatrix[i][j]= (*(SimilarityData->DissimilarityMatrix))(i,j);	
		}//end loop cols
	}//end loop rows
	*/


	std::sort(AbsDissList.begin(),AbsDissList.end());
	double AbsDmin= AbsDissList[0];
	double AbsDmax= AbsDissList[AbsDissList.size()-1];
	double AbsDmedian= Utils::GetMedian(AbsDissList,true);
	double AbsDmad= Utils::GetMAD(AbsDissList,AbsDmedian);
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
	for(int i=0;i<DissMatrix.size();i++){
		std::vector<size_t> sort_index;//sorting index
		std::vector<double> sorted;
		Utils::sort( DissMatrix[i],sorted,sort_index);
		for(unsigned int k=0;k<sort_index.size();k++){
			size_t index= sort_index[k];
			double diss= DissMatrix[i][index];
			//if(diss<=1.e-12 || i==index) continue;
			if(diss<=0 || i==index) continue;
			
			//(SimilarityData->DissimilaritySortIndexMatrix)[i][k]= index;
			(SimilarityData->DissimilaritySortIndexMatrix)[i].push_back(index);
		}
	}//end loop rows

	
	//cout<<"== Final Dissimilarity Matrix (tot diss) =="<<endl;
	//(SimilarityData->DissimilarityMatrix)->Print();
	//cout<<"======================="<<endl;

	/*
	cout<<"== Adjacency Matrix =="<<endl;
	(SimilarityData->AdjacencyMatrix)->Print();
	cout<<"======================="<<endl;
	
	cout<<"== Edgeness Matrix =="<<endl;
	EdgenessMatrix->Print();
	image
	cout<<"======================="<<endl;
	*/

	//Clear up
	EdgenessMatrix->Delete();

	return SimilarityData;

}//close ComputeRegionSimilarity()



int SLICUtils::TagOccludedRegions(Img* image,std::vector< std::vector<long int> > labels,std::vector<Region*> regions,int mainRegionTag,int secondaryRegionTag){

	if(!image) return -1;

	//## Compute region contour info (neighbors, ...)
	cout<<"SLICSegmenter::SuperpixelMaxSimilarityMerger(): INFO: Finding region neighbors (NR="<<regions.size()<<") ..."<<endl;
	SLICUtils::SLICContourData* contourData= SLICUtils::ComputeBoundaryContours(image,labels,regions);
	SLICUtils::SLICConnectedRegions connectedRegionIds= contourData->connectedRegionIds;
	
	//## Decomment if you want to use 2nd neighbors in max similarity merging
	//Fill list of 2-nd neighbors	
	std::vector< std::vector<int> > neighborIndexList_2nd;
	for(unsigned int i=0; i<regions.size(); i++) {
		int regionId= regions[i]->fId;
		int regionTag= regions[i]->fTag;
		regions[i]->fIsOccluded= false;
		
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

	for(unsigned int i=0;i<regions.size();i++) {
		int regionId= regions[i]->fId;
		int regionTag= regions[i]->fTag;
		regions[i]->fIsOccluded= false;
		if(regionTag!=mainRegionTag) continue;

		bool hasSurroundingRegions= false;
		for(unsigned int j=0;j<connectedRegionIds[i].size();j++){
			int neighborIndex= connectedRegionIds[i][j];
			int neighborId= regions[neighborIndex]->fId;
			int neighborTag= regions[neighborIndex]->fTag;
			if(neighborTag==secondaryRegionTag){
				hasSurroundingRegions= true;
				break;
			}
		}//end loop all neighbors

		if(!hasSurroundingRegions) regions[i]->fIsOccluded= true;
	}//end loop regions

	if(contourData){
		delete contourData;
		contourData= 0;
	}	
	
	return 0;

}//close TagOccludedRegions()


std::vector<TText*> SLICUtils::GetRegionTextLabels(std::vector<Region*> regions,int selectedTag){

	std::vector<TText*> regionTextList;
	regionTextList.clear();

	int nRegions= (int)regions.size();
	if(nRegions<=0){
		cerr<<"SLICUtils::GetRegionTextLabels(): WARN: No regions available, nothing to compute!"<<endl;
		return regionTextList;
	}

	TText* regionText= 0;

	for(int k=0;k<nRegions;k++){
		int nPix= regions[k]->fNPix;
		int regionId= regions[k]->fId;
		int regionTag= regions[k]->fTag;
		if(nPix<=0) continue;
		if(selectedTag!=-1 && regionTag!=selectedTag) continue;

		double Cx= regions[k]->fX0;
		double Cy= regions[k]->fY0;
		double textX= Cx;
		double textY= Cy;
		regionText= new TText(textX,textY,Form("%d",regionId));
		regionText->SetTextSize(0.015);
		regionText->SetTextColor(kBlack);
		regionTextList.push_back(regionText);	
	}//end loop new regions	
	
	return regionTextList;
	
}//close GetRegionTextLabels()


