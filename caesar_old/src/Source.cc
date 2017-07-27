/**
* @file Source.cc
* @class Source
* @brief Source
*
* Source class
* @author S. Riggi
* @date 20/01/2015
*/

#include <Source.h>
#include <Img.h>
#include <Contour.h>
#include <EllipseUtils.h>
#include <SourceFitter.h>
#include <ZernikeMoments.h>

#include "SLICSegmenter.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TGraph2D.h>
#include <TF1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TCut.h>
#include <TEventList.h>
#include <TMath.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TColor.h>
#include <TEllipse.h>
#include <TPolyLine.h>
#include <TPaveText.h>
#include <TSpectrum2.h>
#include <TMatrixDEigen.h>

#include <Math/WrappedTF1.h>
#include <Math/GSLIntegrator.h>
#include <Math/GSLMinimizer.h>
#include <Math/Functor.h>
#include <Math/WrappedFunction.h>
#include <Math/WrappedParamFunction.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>
#include <Math/SpecFunc.h>
#include <Math/DistFunc.h>



#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
using namespace cv;


#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>

using namespace std;


ClassImp(Source)

Source::Source(){

	fName= "";
	fId= -999;
	fFlag= eCandidate;
	fType= eUnknown;
	fNPix= 0;
	fIsGoodSource= true;
	fS= 0;	
	fSmin= 1.e+99;
	fSmax= -1.e+99;
	fCurvMin= 1.e+99;
	fCurvMax= -1.e+99;
	fSxx= 0;
	fSyy= 0;
	fSxy= 0;
	fPixIdmax= -1;
	fPixIdmin= -1;
	fF= 0;
	fFlux= 0;
	fWx0= 0;
	fWy0= 0;
	fX0= 0;
	fY0= 0;
	fFluxCorrection= 1;
	
	fNPix= 0;

	fM1= 0;
	fM2= 0;
	fM3= 0;
	fM4= 0;
	fMean= 0;
	fRMS= 0;
	fKurtosis= 0;
	fSkewness= 0;
	fMedian= 0;
	fMedianRMS= 0;

	fXmin= 1.e+99;
	fXmax= -1.e+99;
	fYmin= 1.e+99;
	fYmax= -1.e+99;
	
	fIx_min= 1.e+99;
	fIx_max= -1.e+99;
	fIy_min= 1.e+99;
	fIy_max= -1.e+99;

	fContour= 0;
	fContourCollection.clear();
	fContourCollection.resize(0);	

	fHasPixelsAtEdge= false;
	fPixelCollection.clear();
	fHasParameters= false;

	fDepthLevel= 0;
	fHasNestedSources= false;
	fNestedSource= 0;
	fNestedSourceCollection.clear();

	//Moments= ...
	for(int k=0;k<7;k++) HuMoments[k]= 0;
	ZMMoments.clear();

}//close costructor


Source::~Source(){

	/*
	for(int i=0;i<fContourCollection.size();i++){
		if(fContourCollection[i]) {
			fContourCollection[i]->Delete();		
			fContourCollection[i]= 0;
		}
	}
	if(fContour) {
		fContour->Delete();
		fContour= 0;
	}
	fContourCollection.clear();
	fContourCollection.resize(0);
	*/

}//close destructor


void Source::AddPixel(Pixel pixel){

	//Append pixel to list
	fPixelCollection.push_back(pixel);

	//Update moments
	double w= pixel.S;
	double curv= pixel.Curv;
	double x= pixel.x;	
	double y= pixel.y;
	int ix= pixel.ix;
	int iy= pixel.iy;
	int id= pixel.id;
	
	if(w<fSmin) {
		fSmin= w;
		fPixIdmin= id;
	}
	if(w>fSmax) {
		fSmax= w;
		fPixIdmax= id;
	}

	if(curv<fCurvMin) fCurvMin= curv;
	if(curv>fCurvMax) fCurvMax= curv;

	fS+= w;
	fSxx+= w*x*x;
	fSyy+= w*y*y;
	fSxy+= w*x*y;

	if(x<fXmin) fXmin= x;
	if(x>fXmax) fXmax= x;
	if(y<fYmin) fYmin= y;
	if(y>fYmax) fYmax= y;
	if(ix<fIx_min) fIx_min= ix;
	if(ix>fIx_max) fIx_max= ix;
	if(iy<fIy_min) fIy_min= iy;
	if(iy>fIy_max) fIy_max= iy;

	fWx0+= w*x;
	fWy0+= w*y;
	fX0+= x;
	fY0+= y;

	fNPix++;
  double delta = w - fM1;
  double delta_n = delta/fNPix;
  double delta_n2 = delta_n * delta_n;
  double f = delta * delta_n * (fNPix-1);
  fM1+= delta_n;
  fM4+= f * delta_n2 * (fNPix*fNPix - 3*fNPix + 3) + 6 * delta_n2 * fM2 - 4 * delta_n * fM3;
  fM3+= f * delta_n * (fNPix - 2) - 3 * delta_n * fM2;
  fM2+= f;	


}//close AddPixel()


int Source::ComputeStats(){

	if(fNPix<=0 || fPixelCollection.size()<=0) return -1;
		
	fX0/= (double)(fNPix);
	fY0/= (double)(fNPix);
	fWx0/= fS;
	fWy0/= fS;
	fSxx/= fS;
	fSyy/= fS;
	fSxy/= fS;

	//fF= fS/;//need beam correction to compute flux
	
	fMean= fM1;
	fRMS= 0;
	if(fNPix>2) {
		fRMS= sqrt(fM2/(fNPix-1));
	}

	fSkewness= 0;
	fKurtosis= 0;
	if(fM2!=0) {
  	fSkewness= sqrt(fNPix)*fM3/pow(fM2,1.5);//need to adjust for finite population?
  	fKurtosis= fNPix*fM4/(fM2*fM2);//also given without -3
	}
	
	//## Compute robust stats (median, MAD, ...)	
	std::vector<double> Pixels;
	Pixels.clear();
	Pixels.resize(0);
	
	for(unsigned int i=0;i<fPixelCollection.size();i++){
		double pixelValue= fPixelCollection[i].S;
		if( pixelValue==0 ) continue;
		Pixels.push_back(pixelValue);
	}

	//Sort and compute median for all image	
	std::sort(Pixels.begin(),Pixels.end());
	fMedian= Utils::GetMedian(Pixels,true);
	
	//Compute MAD = median(|x_i-median|)
	//cout<<"Img::ComputeStatsParams(): INFO: Computing MAD for full image..."<<endl;
	std::vector<double> MADs;
	for(unsigned j=0;j<Pixels.size();j++){
		double mad= fabs(Pixels[j]-fMedian);
		MADs.push_back(mad);
	}
	std::sort(MADs.begin(),MADs.end());
	double medianMAD= Utils::GetMedian(MADs,true);
	fMedianRMS= medianMAD*1.4826;//0.6744888;

	return 0;

}//close ComputeStats()


bool Source::IsGoodSource(){
	
	//## Check for pixels 	
	if(fNPix<=0 || fPixelCollection.size()<=0) return false;

	//## Check for line-like source
	if(fContourCollection.size()<=0) {
		cerr<<"Source::IsGoodSource(): WARN: No contour stored for this source, cannot perform check!"<<endl;
		return true;
	}

	double BoundingBoxMin= fContourCollection[0]->BoundingBoxMin;
	if(BoundingBoxMin<2) {
		cerr<<"Source::IsGoodSource(): INFO: BoundingBox cut not passed (BoundingBoxMin="<<BoundingBoxMin<<"<2)"<<endl;
		return false;
	}

	//## Add other check here ...
	//...
	//...

	return true;

}//close Source::isGoodSource()


bool Source::IsCompactPointLike(){

	if(!fHasParameters) {
		cerr<<"Source::IsCompactSource(): WARN: No parameters are available for this source (did you compute them?)...test cannot be performed!"<<endl;
		return true;
	}

	//Loop over contours and check if all of them have circular features
	bool isPointLike= true;
	for(int i=0;i<fContourCollection.size();i++){
		Contour* thisContour= fContourCollection[i];

		//Test circularity ratio: 1= circle
		if(thisContour->CircularityRatio<0.4) {
			cout<<"Source::IsCompactSource(): INFO: Source does not pass CircularityRatio cut (CR="<<thisContour->CircularityRatio<<"<0.4)"<<endl;
			isPointLike= false;
			break;
		}

		//Test elongation (how symmetrical is the shape): 0=circle,square
		if(thisContour->Elongation>0.7) {
			cout<<"Source::IsCompactSource(): INFO: Source does not pass Elongation cut (ELONG="<<thisContour->CircularityRatio<<">0.4)"<<endl;
			isPointLike= false;
			break;	
		}

		//Test ellipse fit
		if(thisContour->EllipseAreaRatio<0.6 || thisContour->EllipseAreaRatio>1.4) {
			cout<<"Source::IsCompactSource(): INFO: Source does not pass EllipseAreaRatio cut (EAR="<<thisContour->EllipseAreaRatio<<"<0.6)"<<endl;
			isPointLike= false;
			break;	
		}

	}//end contour loop
	
	//Check number of pixels
	if(fNPix>500){
		cout<<"Source::IsCompactSource(): INFO: Source does not pass nMaxPix cut (NPix="<<fNPix<<">500)"<<endl;
		isPointLike= false;
	}

	if(!isPointLike) return false;

	return true;

}//close Source::IsCompactSource()


int Source::ComputeMorphologyParams(){

	cout<<"Source::ComputeMorphologyParams(): INFO: Computing morphology parameters..."<<endl;
	if(fNPix<=0 || fPixelCollection.size()<=0) return -1;
		
	//######################################
	//## Find the source bounding box
	//######################################
	double xRange[2]= {fXmin,fXmax};
	double yRange[2]= {fYmin,fYmax};	
	int ixRange[2]= {fIx_min,fIx_max};
	int iyRange[2]= {fIy_min,fIy_max};
	
	//Bounding box in (x,y) coordinates
	int boundingBoxX[2];
	int boundingBoxY[2];
	int deltaPix= 50;
	boundingBoxX[0]= xRange[0]-deltaPix;
	boundingBoxX[1]= xRange[1]+deltaPix;
	boundingBoxY[0]= yRange[0]-deltaPix;
	boundingBoxY[1]= yRange[1]+deltaPix;
	int nBoxX= boundingBoxX[1]-boundingBoxX[0]+1;
	int nBoxY= boundingBoxY[1]-boundingBoxY[0]+1;

	//Bounding box in (ix,iy) coordinates
	int boundingBoxIX[2];
	int boundingBoxIY[2];	
	boundingBoxIX[0]= ixRange[0]-deltaPix;
	boundingBoxIX[1]= ixRange[1]+deltaPix;
	boundingBoxIY[0]= iyRange[0]-deltaPix;
	boundingBoxIY[1]= iyRange[1]+deltaPix;
	int nBoxIX= boundingBoxIX[1]-boundingBoxIX[0]+1;
	int nBoxIY= boundingBoxIY[1]-boundingBoxIY[0]+1;

	//cout<<"Source::ComputeMorphologyParams(): INFO: xRange("<<xRange[0]<<","<<xRange[1]<<"), yRange("<<yRange[0]<<","<<yRange[1]<<")"<<endl;
	//cout<<"Source::ComputeMorphologyParams(): INFO: boundingBoxX("<<boundingBoxX[0]<<","<<boundingBoxX[1]<<"), boundingBoxY("<<boundingBoxY[0]<<","<<boundingBoxY[1]<<")"<<"  nBoxX="<<nBoxX<<", nBoxY="<<nBoxY<<endl;
	//cout<<"Source::ComputeMorphologyParams(): INFO: boundingBoxIX("<<boundingBoxIX[0]<<","<<boundingBoxIX[1]<<"), boundingBoxIY("<<boundingBoxIY[0]<<","<<boundingBoxIY[1]<<")"<<"  nBoxIX="<<nBoxIX<<", nBoxIY="<<nBoxIY<<endl;


	//## Fill image and binarized image
	//TString imgName= Form("OriginalImg_%s",fName.c_str());
	//Img* originalImg= new Img(imgName,imgName,nBoxX,boundingBoxX[0]-0.5,boundingBoxX[1]+0.5,nBoxY,boundingBoxY[0]-0.5,boundingBoxY[1]+0.5);
	//originalImg->Sumw2();

	cv::Mat binarizedImg = cv::Mat::zeros(nBoxIY, nBoxIX, CV_8UC1);
	cv::Mat rasterImg = cv::Mat::zeros(nBoxIY, nBoxIX, CV_64FC1);

	for(unsigned int k=0;k<fPixelCollection.size();k++){
		Pixel thisPixel= fPixelCollection[k];
		double thisS= thisPixel.S;	
		double thisZ= thisPixel.Z;	
		int thisId= thisPixel.id;
		double thisX= thisPixel.x;
		double thisY= thisPixel.y;
		int ix= thisPixel.ix;
		int iy= thisPixel.iy;
		//originalImg->FillPixel(thisX,thisY,thisZ);
		
		ix-= boundingBoxIX[0];
		iy-= boundingBoxIY[0];
		int rowId= nBoxIY-1-iy;
		int colId= nBoxIX-1-ix;
		binarizedImg.at<uchar>(rowId, colId, 0) = 1;
		rasterImg.at<double>(rowId, colId, 0) = thisS;		
	}//end loop pixels


	//## Compute contour
	cv::Mat binarizedImg_clone= binarizedImg.clone();//image will be modified by findContours! 
	std::vector<std::vector<cv::Point>> contours; // Vector for storing contour
  std::vector<cv::Vec4i> hierarchy;
	cv::findContours( binarizedImg_clone, contours, hierarchy,CV_RETR_EXTERNAL,CV_CHAIN_APPROX_NONE, Point(0,0) ); // Find only the external contours in the image
  //cv::findContours( binarizedImg_clone, contours, hierarchy,CV_RETR_CCOMP, CV_CHAIN_APPROX_NONE ); // Find the contours in the image organizing in a 2-level hierarchy
	//cv::findContours( binarizedImg_clone, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_NONE, Point(0,0) );
	
	for(unsigned int i=0; i<contours.size(); i++){ // iterate through each contour
		int nContourPts= (int)contours[i].size();
		if(nContourPts<=0) continue;

		//Create and fill contour
		fContour= new Contour;

		cout<<"Source::ComputeMorphologyParams(): INFO: Contour no. "<<i+1<<": (";
		
		for(int j=0;j<nContourPts;j++){
			int contx= contours[i][j].x;
			int conty= contours[i][j].y;
			int rowId= nBoxIY-1-conty;
			int colId= nBoxIX-1-contx;
			int contx_transf= colId + boundingBoxX[0];
			int conty_transf= rowId + boundingBoxY[0];
			fContour->AddPoint(cv::Point2f(contx_transf,conty_transf));
			cout<<"("<<contx_transf<<","<<conty_transf<<"), ";
		}//end loop points in contour
		cout<<")"<<endl;
		

		//Compute contour parameters
		if(fContour->ComputeParameters()<0){
			cerr<<"Source::ComputeMorphologyParams(): WARN: Failed to compute parameters for contour no. "<<i<<"!"<<endl;
			//delete fContour;
			//continue;
		}
		
		//Add contour to the list
		fContourCollection.push_back(fContour);	

	}//end loop contours


	//## Compute HuMoments
	moments= cv::moments(rasterImg, false);
	cv::HuMoments(moments, HuMoments);

	//## Compute zernike moments
	int order= 6;
	ComputeZernikeMoments(order);
	
	fHasParameters= true;

	return 0;


}//close Source::ComputeMorphologyParams()


int Source::ComputeZernikeMoments(int order){

	//## Get source image
	Img* fluxMap= GetImage(Source::eFluxMap);
	if(!fluxMap) return -1;
	
	//## Compute Zernike moments
	double radius= -1;
	ZMMoments.clear();
	ZMMoments= ZernikeMoments::GetZernike2D_Direct(fluxMap, order,radius);
	//ZMMoments= ZernikeMoments::GetZernike2D(fluxMap, order,radius);

	fluxMap->Delete();

	return 0;

}//close 

bool Source::IsPointInsideSource(double x,double y,double distTolerance){

	//Check contour
	if(fContourCollection.size()<=0){
		cerr<<"Source::IsPointInsideSource(): WARN: No contour available, check cannot be performed!"<<endl;
		return false;
	}
	
	//Loop over external contours (for disjoint regions, i.e. those merged, there are more than 1 contour)
	bool isInsideAnyContour= false;
	bool measureDist= true;
	
	for(unsigned int s=0;s<fContourCollection.size();s++){
		std::vector<cv::Point2f> contour_points= fContourCollection[s]->GetPoints();
		double distanceToEdge= cv::pointPolygonTest(contour_points, cv::Point2f(x,y), measureDist);

		//If outside this contour, skip and test next
		if(distanceToEdge<0){//outside contour
			continue;
		}
		else if(distanceToEdge==0) {//on contour edge
			continue;
		}
		else if(distTolerance>0 && distanceToEdge>distTolerance) {//inside contour
			isInsideAnyContour= true;
			break;						
		}//close if inside contour

	}//end loop contours

	return isInsideAnyContour;

}//close IsPointInsideSource()

int Source::FindNestedSource(double curvThr,int nPixMin,double peakThr){

	cout<<"Source::FindNestedSource(): INFO: Searching for nested sources..."<<endl;
	
	//Reset current nested source collection
	fNestedSourceCollection.clear();
	fNestedSourceCollection.resize(0);
	fNestedSource= 0;
	fHasNestedSources= false;

	//## Get source/curv/significance images
	Img* sourceImg= GetImage(eFluxMap);
	Img* significanceImg= GetImage(eSignificanceMap);
	Img* curvImg= GetImage(Source::eFluxCurvMap);
	if(!sourceImg || !significanceImg || !curvImg){
		cerr<<"Source::FindNestedSource(): ERROR: Cannot get source/significance or curv map!"<<endl;
		if(sourceImg) sourceImg->Delete();
		if(significanceImg) significanceImg->Delete();
		if(curvImg) curvImg->Delete();
		return -1;
	}

	//## Find peaks in curvature map
	int tol= 1;//1 pixel-tolerance
	TGraph* peaks= curvImg->FindPeaks(tol);
	if(!peaks) {
		cout<<"Source::FindNestedSource(): INFO: No peaks detected!"<<endl;
		sourceImg->Delete();
		significanceImg->Delete();
		curvImg->Delete();
		return 0;
	}
	int nPeaks= peaks->GetN();
	cout<<"Source::FindNestedSource(): INFO: #"<<nPeaks<<" peaks detected!"<<endl;

	//## Find blobs image by thresholding the curvature map and then by flood-fill
	Img* curvImg_binary= curvImg->GetBinarized(curvThr);
	Img* maskedImg= sourceImg->GetMask(curvImg_binary,false);

	for(int k=0;k<nPeaks;k++){
		double x, y;
		peaks->GetPoint(k,x,y);
		int peakBinId= sourceImg->FindBin(x,y);
		double peakZ= significanceImg->GetBinContent(peakBinId);
		if(sourceImg->IsBinOverflow(peakBinId) || sourceImg->IsBinUnderflow(peakBinId) ) continue;
		if(peakZ<peakThr) continue;//skip peaks below threshold
		curvImg_binary->AddBinContent(peakBinId,1);
	}	
	double fgValue= 1;	
	bool findNegativeExcess= false;
	bool findNestedSources= false;
	bool useLocalBackground= false;
	bool mergeBelowSeed= false;//true!
	//std::vector<int> nestedBlobsIds= curvImg_binary->FloodFill(fgValue+1,fgValue,mergeBelowSeed);

	int status= maskedImg->FindCompactSource(curvImg_binary,fgValue+1,fgValue,nPixMin,findNegativeExcess,findNestedSources,useLocalBackground,mergeBelowSeed);
	if(status<0){
		cerr<<"Source::FindNestedSource(): WARN: Blob flood-filling failed!"<<endl;
		sourceImg->Delete();
		significanceImg->Delete();
		curvImg->Delete();
		peaks->Delete();	
		curvImg_binary->Delete();
		maskedImg->Delete();
		return -1;
	}

	//## Set nested sources
	std::vector<Source*> NestedSources= maskedImg->GetSources();
	int nNestedSources= (int)NestedSources.size();
	if(nNestedSources<=0){
		cout<<"Source::FindNestedSource(): INFO: No nested source found!"<<endl;
		sourceImg->Delete();
		significanceImg->Delete();
		curvImg->Delete();
		peaks->Delete();	
		curvImg_binary->Delete();
		maskedImg->Delete();
		return 0;
	}
	cout<<"Source::FindNestedSource(): INFO: #"<<nNestedSources<<" nested sources found..."<<endl;	
	
	for(int i=0;i<nNestedSources;i++) {
		fNestedSource= new Source;
		*(fNestedSource)= *(NestedSources[i]);
		fNestedSource->fType= NestedSources[i]->fType;
		fNestedSource->fDepthLevel= this->fDepthLevel+1;
		fNestedSourceCollection.push_back(fNestedSource);	
	}//end loop nested sources

	fHasNestedSources= true;

	//## Clear-up
	sourceImg->Delete();
	significanceImg->Delete();
	curvImg->Delete();
	peaks->Delete();	
	curvImg_binary->Delete();
	maskedImg->Delete();
	
	return 0;

}//close FindNestedSource()

int Source::FindNestedSource(double seedThr,double mergeThr,int nPixMin,bool findPointLike){

	cout<<"Source::FindNestedSource(): INFO: Searching for nested sources..."<<endl;

	//Reset current nested source collection
	fNestedSourceCollection.clear();
	fNestedSourceCollection.resize(0);
	fNestedSource= 0;
	fHasNestedSources= false;

	//## Get source image
	Img* sourceImg= GetImage(eFluxMap);
	Img* sourceCurvatureImg= GetImage(eFluxCurvMixtureMap,0.7);

	//## Compute source image stats
	int status= sourceCurvatureImg->ComputeStats(true,false,true);//recomputing stats
	if(status<0) {
		cerr<<"Source::FindNestedSource(): ERROR: Stats computing failed!"<<endl;
		return -1;
	}
	sourceCurvatureImg->DumpStats();

	//## Compute source image bkg
	//if(fNPix>10000) status= sourceImg->ComputeBkg(Img::eRobustBkg); 
	//else status= sourceImg->ComputeBkg(Img::eMedianBkg); 	
	status= sourceCurvatureImg->ComputeBkg(Img::eMedianBkg); 
	if(status<0) {
		cerr<<"Source::FindNestedSource(): ERROR: Bkg computing failed!"<<endl;
		return -1;
	}
	sourceCurvatureImg->DumpBkg();
	this->Dump();

	//## Compute source significance map	
	bool useLocalBackground= false;
	Img* significanceMap= sourceCurvatureImg->GetSignificanceMap(useLocalBackground);
	if(!significanceMap) {
		cerr<<"Source::FindNestedSource(): ERROR: Failed to compute significance map!"<<endl;
		return -1;
	}

	//## Compute nested sources
	bool findNegativeExcess= false;
	bool findNestedSources= false;
	bool useCurvatureMixture= true;
	double curvatureMixtureFraction= 0.7;
	status= sourceImg->FindCompactSource(significanceMap,seedThr,mergeThr,nPixMin,findNegativeExcess,findNestedSources);
	//status= sourceImg->FindCompactSource(significanceMap,5,3,8,findNegativeExcess,findNestedSources);
	if(status<0) {
		cerr<<"Source::FindNestedSource(): ERROR: Source finding failed!"<<endl;
		return -1;
	}


	//Get source list
	std::vector<Source*> NestedSources= sourceImg->GetSources();
	cout<<"Source::FindNestedSource(): INFO: "<<NestedSources.size()<<" nested source available..."<<endl;	

	for(unsigned int i=0;i<NestedSources.size();i++) {
		if(findPointLike && NestedSources[i]->fType!=ePointLike ) {
			cout<<"Source::FindNestedSource(): INFO: Nested blob no. "<<i<<" skipped by point-like cut..."<<endl;
			continue;
		}
		fNestedSource= new Source;
		*(fNestedSource)= *(NestedSources[i]);
		fNestedSource->fType= NestedSources[i]->fType;
		fNestedSource->fDepthLevel= this->fDepthLevel+1;
		fNestedSourceCollection.push_back(fNestedSource);	
	}//end loop sources

	cout<<"Source::FindNestedSource(): INFO: "<<fNestedSourceCollection.size()<<" nested sources found after selection cuts!"<<endl;
	if(fNestedSourceCollection.size()>0) fHasNestedSources= true;

	//## Clear-up
	if(sourceImg) sourceImg->Delete();
	if(sourceCurvatureImg) sourceCurvatureImg->Delete();
	if(significanceMap) significanceMap->Delete();
	
	/*
	for(int i=0;i<NestedSources.size();i++){
		if(NestedSources[i]) {
			delete NestedSources[i];
		}
	}
	*/
	return 0;

}//close FindNestedSource()


int Source::FindNestedClusterRegions(int regionSize,double regularization,int minRegionSize,double SPMergingRatio,double SPMergingRegularization,double MaxDissRatio){

	//Get source image
	Img* sourceImg= GetImage(Source::eFluxMap,0);
		
	//Parameters
	int MinMergedSP= 1;
	bool SPMergingIncludeSpatialPars= false;
	double PixelRatioCut= 0.4;
	bool UseCurvatureInSPMerging= true;
	double SignificantSPRatio= 0.5;
	int SPMergingEdgeModel= 1;
	double SPMergingDistThreshold= 3.5;

	//Compute superpixel partition
	SLICSegmenter slic;
	slic.SetLogContrastMapping(false);
	slic.SetSPMergingRatio(SPMergingRatio);
	slic.SetSPMergingRegularization(SPMergingRegularization);
	slic.SetMinMergedSP(MinMergedSP);
	slic.SetSPMergingDistThreshold(SPMergingDistThreshold);
	slic.UsePixelRatioCut(true);
	slic.SetPixelRatioCut(PixelRatioCut);
	slic.SetSignificantSPRatioCut(SignificantSPRatio);
	slic.SetMaxDissRatio(MaxDissRatio);
	slic.SetMaxDissRatioFor2ndNeighbors(MaxDissRatio);		
	slic.UseCurvatureInSPMerging(UseCurvatureInSPMerging);
	slic.SetEdgeModel(SPMergingEdgeModel);

	int status= slic.RunSegmentation(sourceImg, regionSize,regularization, minRegionSize, false);
	if(status<0){
		cerr<<"Source::FindNestedClusterRegions(): ERROR: Superpixel partition failed!"<<endl;
		sourceImg->Delete();
		return -1;
	}
	
	std::vector<Region*> regions= slic.GetRegions();
	std::vector< std::vector<long int> > labels= slic.GetPixelClusterIds();

	//Set all regions as signal
	int nRegions= (int)regions.size();
	for(int i=0;i<nRegions;i++){
		regions[i]->fTag= Region::eSignalTag;
	}//end loop regions

	//Hierarchical merge regions
	status= slic.SPHierarchicalMerger(regions,labels,-1,-1);		
	if(status<0){
		cerr<<"Source::FindNestedClusterRegions(): ERROR: Superpixel merging failed!"<<endl;
		sourceImg->Delete();
		return -1;
	}

	std::vector<Region*> regions_merged= slic.GetRegions();
	int nMergedRegions= (int)regions_merged.size();
	cout<<"Source::FindNestedClusterRegions(): INFO: #"<<nMergedRegions<<" regions merged out of "<<nRegions<<" ..."<<endl;
	if(nMergedRegions==nRegions || nMergedRegions==1){
		cerr<<"Source::FindNestedClusterRegions(): WARN: No regions merged or only 1 merged regions at the end!"<<endl;
		return -1;
	}

	for(int i=0;i<nMergedRegions;i++){
		Source* aSource= regions_merged[i]->GetSource();
		if(aSource) this->AddNestedSource(aSource);
	}//end loop merged regions

	return 0;

}//close FindNestedClusterRegions()


int Source::Deblend(double curvThr,int componentMinNPix){
	int status= SourceFitter::FitSource(this,curvThr,componentMinNPix);
	return status;
}//close Deblend()



Img* Source::GetImage(SourceImageMode mode,double curvMixtureFract){

	//Bounding box in (x,y) coordinates
	double xRange[2]= {fXmin,fXmax};
	double yRange[2]= {fYmin,fYmax};	
	int boundingBoxX[2];
	int boundingBoxY[2];
	int deltaPix= 50;
	boundingBoxX[0]= xRange[0]-deltaPix;
	boundingBoxX[1]= xRange[1]+deltaPix;
	boundingBoxY[0]= yRange[0]-deltaPix;
	boundingBoxY[1]= yRange[1]+deltaPix;
	int nBoxX= boundingBoxX[1]-boundingBoxX[0]+1;
	int nBoxY= boundingBoxY[1]-boundingBoxY[0]+1;

	//## Fill image and binarized image
	TString imgName= Form("SourceImg_%s_mode%d",fName.c_str(),mode);
	Img* sourceImg= new Img(imgName,imgName,nBoxX,boundingBoxX[0]-0.5,boundingBoxX[1]+0.5,nBoxY,boundingBoxY[0]-0.5,boundingBoxY[1]+0.5);
	Img* sourceSignificanceMap= 0;

	if(mode==eBinarizedMap){
		for(unsigned int k=0;k<fPixelCollection.size();k++){
			Pixel thisPixel= fPixelCollection[k];
			double thisX= thisPixel.x;
			double thisY= thisPixel.y;
			sourceImg->FillPixel(thisX,thisY,1);
		}//end loop pixels
	}//close if
	else if(mode==eFluxMap){
		for(unsigned int k=0;k<fPixelCollection.size();k++){
			Pixel thisPixel= fPixelCollection[k];
			double thisX= thisPixel.x;
			double thisY= thisPixel.y;
			double thisS= thisPixel.S;	
			double thisZ= thisPixel.Z;	
			sourceImg->FillPixel(thisX,thisY,thisS);
			//sourceImg->FillPixel(thisX,thisY,thisZ);
		}//end loop pixels
	}//close else
	else if(mode==eSignificanceMap){
		for(unsigned int k=0;k<fPixelCollection.size();k++){
			Pixel thisPixel= fPixelCollection[k];
			double thisX= thisPixel.x;
			double thisY= thisPixel.y;
			double thisS= thisPixel.S;	
			double thisZ= thisPixel.Z;	
			sourceImg->FillPixel(thisX,thisY,thisZ);
		}//end loop pixels
	}//close else
	else if(mode==eSourceSignificanceMap){
		for(unsigned int k=0;k<fPixelCollection.size();k++){
			Pixel thisPixel= fPixelCollection[k];
			double thisX= thisPixel.x;
			double thisY= thisPixel.y;
			double thisS= thisPixel.S;	
			double sourceZ= (thisS-fMedian)/fMedianRMS;	
			sourceImg->FillPixel(thisX,thisY,sourceZ);
		}//end loop pixels
	}//close else if
	else if(mode==eFluxCurvMap){
		for(unsigned int k=0;k<fPixelCollection.size();k++){
			Pixel thisPixel= fPixelCollection[k];
			double thisX= thisPixel.x;
			double thisY= thisPixel.y;
			double thisCurv= thisPixel.Curv;	
			sourceImg->FillPixel(thisX,thisY,thisCurv);
		}//end loop pixels
	}//close else if
	else if(mode==eFluxCurvMixtureMap || mode==eFluxCurvMixtureSignificanceMap){//take a mixture 
		double normmin= 1;
		double normmax= 256;
		for(unsigned int k=0;k<fPixelCollection.size();k++){
			Pixel thisPixel= fPixelCollection[k];
			double thisX= thisPixel.x;
			double thisY= thisPixel.y;
			double thisS= thisPixel.S;	
			double thisCurv= thisPixel.Curv;	
			double thisS_norm= normmin + (normmax-normmin)*(thisS-fSmin)/(fSmax-fSmin);
			double thisCurv_norm= normmin + (normmax-normmin)*(thisS-fCurvMin)/(fCurvMax-fCurvMin);
			double weightedSignal= curvMixtureFract*thisCurv_norm + (1.-curvMixtureFract)*thisS_norm;
			sourceImg->FillPixel(thisX,thisY,weightedSignal);
		}//end loop pixels

		
		if(mode==eFluxCurvMixtureSignificanceMap){
			sourceImg->ComputeStats(true,false,true);//recomputing stats		
			sourceImg->ComputeBkg(Img::eMedianBkg); 
			sourceSignificanceMap= sourceImg->GetSignificanceMap(false);	
			return sourceSignificanceMap;
		} 

	}//close else if
	else if(mode==eMeanFluxMap){
		for(unsigned int k=0;k<fPixelCollection.size();k++){
			Pixel thisPixel= fPixelCollection[k];
			double thisX= thisPixel.x;
			double thisY= thisPixel.y;
			double thisCurv= thisPixel.Curv;	
			sourceImg->FillPixel(thisX,thisY,fMean);
		}//end loop pixels
	}

	return sourceImg;

}//close Source::GetSourceImage()

bool Source::IsOverlapping(Source* aSource,double& overlapArea,double overlapThreshold){

	bool hasOverlap= false;
	overlapArea= -1;

	if(!aSource || !aSource->HasParameters() || this->fHasParameters){
		cerr<<"Source::IsOverlapping(): WARN: Invalid source to compare or sources have not computed parameters!"<<endl;
		overlapArea= -1;
		return false;
	}

	std::vector<Contour*> SourceContours;

	//Loop over contours to get ellipse
	double err;
	TGraph overlapAreaGraph;
	int algoChoice= 1;
	bool returnRelArea= true;
	double maxOverlapArea= -1.e+99;

	for(unsigned int i=0;i<fContourCollection.size();i++){
		TEllipse* thisEllipse= fContourCollection[i]->GetFittedEllipse();
		
		for(unsigned int j=0;j<SourceContours.size();j++){
			TEllipse* thisSourceEllipse= SourceContours[j]->GetFittedEllipse();

			if(!thisEllipse || !thisSourceEllipse) continue;
			
			double relOverlapArea = EllipseUtils::ComputeEllipseOverlap(thisEllipse,thisSourceEllipse,overlapAreaGraph,err,returnRelArea,algoChoice);
			if(relOverlapArea<0) continue;

			if(relOverlapArea>maxOverlapArea) maxOverlapArea= relOverlapArea;
		}//end source contour loop
	}//end this source contour loop

	overlapArea= maxOverlapArea;
	if(overlapArea>=overlapThreshold) hasOverlap= true;
	cout<<"Source::IsOverlapping(): INFO: overlapArea="<<overlapArea<<", hasOverlap? "<<hasOverlap<<" (thr="<<overlapThreshold<<")"<<endl;
	
	return hasOverlap;

}//close Source::IsOverlapping()


void Source::Draw(SourceImageMode mode,double curvMixtureFract,bool drawNested,bool drawOnlyNested,bool drawContour,bool drawEllipse,bool drawBoundingBox,bool drawLegends){

	cout<<"Source::Draw(): INFO: Filling source image..."<<fName<<endl;
	Img* sourceImg= GetImage(mode,curvMixtureFract);

	//## Draw stuff for check
	cout<<"Source::Draw(): INFO: Drawing source "<<fName<<endl;
	TString canvasName= Form("SourcePlot_%s",fName.c_str());
	TCanvas* SourcePlot= new TCanvas(canvasName,canvasName);
	SourcePlot->cd();

	if(sourceImg && !drawOnlyNested) {
		sourceImg->SetStats(0);
		sourceImg->Draw("COLZ");
	}
	
	cout<<"Source::Draw(): INFO: Drawing source "<<fName<<" contours..."<<endl;
	
	//Drawing contours?
	for(unsigned int i=0;i<fContourCollection.size();i++){
		cout<<"Source::Draw(): INFO: Drawing source "<<fName<<" contour no. "<<i<<"..."<<endl;

		if(drawContour){			
			TGraph* thisContourGraph= fContourCollection[i]->GetGraph();
			if(thisContourGraph) {
				thisContourGraph->SetMarkerSize(8);
				thisContourGraph->SetMarkerSize(0.3);
				thisContourGraph->SetMarkerColor(kBlack);
				thisContourGraph->SetLineColor(kBlack);
				thisContourGraph->SetLineStyle(kSolid);
				thisContourGraph->SetLineWidth(2);
				thisContourGraph->Draw("PLsame");
			}
		}//close if draw contour

		//Get bounding box
		if(drawBoundingBox){
			cout<<"Source::Draw(): INFO: Drawing source "<<fName<<" bounding box..."<<endl;
			TPolyLine* thisBoundingBox= fContourCollection[i]->GetBoundingBoxLine();
			if(thisBoundingBox){
				thisBoundingBox->SetLineColor(kGray+2);
				thisBoundingBox->SetLineStyle(kDashed);
				thisBoundingBox->Draw("lsame");
			}
		}

		//Get fitted ellipse
		if(drawEllipse){
			cout<<"Source::Draw(): INFO: Drawing source "<<fName<<" fitted ellipse..."<<endl;
			TEllipse* thisFittedEllipse= fContourCollection[i]->GetFittedEllipse();
			if(thisFittedEllipse){
				thisFittedEllipse->SetLineColor(kGray+2);
				thisFittedEllipse->SetLineStyle(kDotted);
				thisFittedEllipse->SetLineWidth(2);
				thisFittedEllipse->SetFillColor(0);
				thisFittedEllipse->SetFillStyle(0);
				thisFittedEllipse->Draw("lsame");
			}
		}//close if draw ellipse

			//Get parameter info box
			if(drawLegends){
				cout<<"Source::Draw(): INFO: Drawing source "<<fName<<" parameter info box..."<<endl;
				TPaveText* thisParamInfoBox= fContourCollection[i]->GetParamInfoBox();
				if(thisParamInfoBox){
					thisParamInfoBox->Draw("same");
				}
			}//close if

		}//end loop contours

		if(fHasParameters && drawLegends){	
			cout<<"Source::Draw(): INFO: Drawing source "<<fName<<" source parameters..."<<endl;
			TPaveText* sourceParamText = new TPaveText(0.4,0.15,0.8,0.3,"NDC");
			sourceParamText->AddText(Form("Id: %d, Type: %d, nNestedSource: %d, AtEdge?%d",fId,fType,fNestedSourceCollection.size(),fHasPixelsAtEdge));
			sourceParamText->AddText(Form("Xmin/Xmax: %1.2f/%1.2f, Ymin/Ymax: %1.2f/%1.2f",fXmin,fXmax,fYmin,fYmax));
			sourceParamText->AddText(Form("C(%1.2f,%1.2f), Wc(%1.2f,%1.2f)",fX0,fY0,fWx0,fWy0));
			sourceParamText->AddText(Form("Stot(#muJy): %1.1f, Smin/Smax(#muJy): %1.2f/%1.2f",fS*1.e+6,fSmin*1.e+6,fSmax*1.e+6));	
			sourceParamText->AddText(Form("Moments(#muJy): (%1.1f,%1.1f,%1.1f,%1.1f)",fMean*1.e+6,fRMS*1.e+6,fSkewness*1.e+6,fKurtosis*1.e+6));
			sourceParamText->AddText(Form("Median(#muJy): (%1.1f,%1.1f)",fMedian*1.e+6,fMedianRMS*1.e+6));
			sourceParamText->SetTextAlign(12);
			sourceParamText->SetTextSize(0.02);
			sourceParamText->SetTextFont(52);
			sourceParamText->SetFillColor(0);
			sourceParamText->SetBorderSize(1);	
			sourceParamText->Draw("same");
		}

		if(fHasFitInfo && drawLegends){
			cout<<"Source::Draw(): INFO: Drawing source "<<fName<<" fit info ..."<<endl;
			TEllipse* fitComponentEllipse= 0;
		
			for(unsigned int k=0;k<fFitInfo.size();k++){
				double Cx= fFitInfo[k]->Cx;	
				double Cy= fFitInfo[k]->Cy;	
				double theta= fFitInfo[k]->theta;
				double MajAxis= fFitInfo[k]->ellMaj;
				double MinAxis= fFitInfo[k]->ellMin;
				fitComponentEllipse= new TEllipse(Cx,Cy,MajAxis,MinAxis,0,360,theta);
				fitComponentEllipse->SetLineWidth(2);
				fitComponentEllipse->SetFillColor(0);
				fitComponentEllipse->SetFillStyle(0);
				fitComponentEllipse->Draw("lsame");
			}//end loop components
		}//close if
	

	//## Draw nested source if any	
	if(drawNested){
		cout<<"Source::Draw(): INFO: Drawing nested contours ..."<<endl;
	
		if(drawOnlyNested){
			Img* nestedSourceImg= 0;
			for(unsigned int k=0;k<fNestedSourceCollection.size();k++){
				if(k==0) nestedSourceImg= fNestedSourceCollection[k]->GetImage(mode,curvMixtureFract);
				else nestedSourceImg->Add(fNestedSourceCollection[k]->GetImage(mode,curvMixtureFract));
			}//end loop nested sources
			if(nestedSourceImg) {
				nestedSourceImg->SetStats(0);	
				nestedSourceImg->Draw("COLZ");
			}
		}//close if

		for(unsigned int k=0;k<fNestedSourceCollection.size();k++){
			std::vector<Contour*> NestedSourceContours= fNestedSourceCollection[k]->fContourCollection;
			if(!fNestedSourceCollection[k]->fHasParameters) cerr<<"No parameters for this source!"<<endl;
		
			for(unsigned int i=0;i<NestedSourceContours.size();i++){
				cout<<"Source::Draw(): INFO: Drawing nested contours no. "<<k<<"..."<<endl;
	
				if(drawContour && !drawOnlyNested){
					TGraph* thisNestedContourGraph= NestedSourceContours[i]->GetGraph();
					if(thisNestedContourGraph) {
						thisNestedContourGraph->SetMarkerSize(8);
						thisNestedContourGraph->SetMarkerSize(0.3);
						thisNestedContourGraph->SetMarkerColor(kRed);
						thisNestedContourGraph->SetLineColor(kRed);
						thisNestedContourGraph->SetLineStyle(kSolid);
						thisNestedContourGraph->SetLineWidth(2);
						thisNestedContourGraph->Draw("PLsame");
					}
				}//close if draw contours

				if(drawBoundingBox && !drawOnlyNested){
					//Get bounding box
					cout<<"Source::Draw(): INFO: Drawing nested bounding box no. "<<k<<"..."<<endl;
	
					TPolyLine* thisNestedBoundingBox= NestedSourceContours[i]->GetBoundingBoxLine();
					if(thisNestedBoundingBox){
						thisNestedBoundingBox->SetLineColor(kGray+2);
						thisNestedBoundingBox->SetLineStyle(kDashed);
						thisNestedBoundingBox->Draw("lsame");
					}
				}//close if draw bounding box

				if(drawEllipse){
					//Get fitted ellipse
					cout<<"Source::Draw(): INFO: Drawing nested fitted ellipse no. "<<k<<"..."<<endl;
					TEllipse* thisNestedFittedEllipse= NestedSourceContours[i]->GetFittedEllipse();
					if(thisNestedFittedEllipse){
						thisNestedFittedEllipse->SetLineColor(kGray+2);
						thisNestedFittedEllipse->SetLineStyle(kDotted);
						thisNestedFittedEllipse->SetLineWidth(2);
						thisNestedFittedEllipse->SetFillColor(0);
						thisNestedFittedEllipse->SetFillStyle(0);
						thisNestedFittedEllipse->Draw("lsame");
					}
				}//close if draw ellipse
			}//end loop nested source contours
		}//end loop nested sources
	}//close if draw nested

}//close Draw()


