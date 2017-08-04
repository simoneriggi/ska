/**
* @file Region.cc
* @class Region
* @brief Region
*
* Image region
* @author S. Riggi
* @date 22/06/2015
*/

#include <Region.h>
#include <Img.h>
#include <Source.h>

#include <Rcpp.h>
#include <RInside.h>

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
//using namespace cv;

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <chrono>
using namespace std;
using namespace std::chrono;

ClassImp(Region)

Region::Region(){

	fId= -999;
	m_HasStats= false;

	ResetMoments();

	fMean= 0;
	fRMS= 0;
	fKurtosis= 0;
	fSkewness= 0;
	fMedian= 0;
	fMedianRMS= 0;
	fMedian_curv= 0;
	fMedianRMS_curv= 0;

	fMean_curv= 0;
	fRMS_curv= 0;
	fH= 0;
	fParamH= 0;

	fIsSignificative= true;
	fMahalanobisDistance= 0;
	fIsSalient= true;
	fSaliency= 0;
	fTag= eUntagged;
	
	fImageMinX= 0;
	fImageMaxX= 0;
	fImageMinY= 0;
	fImageMaxY= 0;
	fImageSizeX= 0;
	fImageSizeY= 0;
	fImageMinS= 0;
	fImageMaxS= 0;
	fImageMinScurv= 0;
	fImageMaxScurv= 0;
	fImageMinSedge= 0;
	fImageMaxSedge= 0;
	fImageRMS= 0;

	fPixelCollection.clear();
	fPixelCollection.resize(0);

	fNeighbourRegions.clear();
	fNeighbourRegions.resize(0);
	fNeighbourRegionInfo.clear();
	fNeighbourRegionInfo.resize(0);

	fHasParameters= false;
	fContour= 0;
	fContourCollection.clear();
	fContourCollection.resize(0);	

	fContourPointIndex.clear();
	fAroundContourPoints.clear();

	fSubRegionIds.clear();

	fAbsDissMin= 1.e+99;
	fAbsDissMax= -1.e+99;
	fIsOccluded= false;

	//Added (CHECK FOR BUG!!!)
  fNPix= 0;
	fX0= 0;
	fY0= 0;

}//close costructor



Region::~Region(){
	
	//cout<<"Region::~Region(): INFO: Deleting region "<<fId<<" ..."<<endl;
  	
	/*
	if(fContour) {
		cout<<"delete contour: pto 1"<<endl;
	
		//fContour->Delete();
		delete fContour;
		cout<<"delete contour: pto 2"<<endl;
	
		fContour= 0;
		cout<<"delete contour: pto 3"<<endl;
	}
	
	if(fRegionImg){
		cout<<"delete image: pto 1"<<endl;
		//fRegionImg->Delete();
		delete fRegionImg;
		cout<<"delete img: pto 2"<<endl;
		fRegionImg= 0;
	}
	*/

}//close destructor




void Region::ResetMoments(){
	fNPix= 0;//npixels	
  fM1= 0;//1st moments
  fM2= 0;//2nd moment
	fM3= 0;//3rd moment
	fM4= 0;//4th moment
	fS= 0;
	fS_curv= 0.;
	fM1_curv= 0;//1st moments
  fM2_curv= 0;//2nd moment
	fS_edge= 0.;
	fSmax= -1.e+99;
	fSmin= 1.e+99;
	fSxx= 0;
	fSyy= 0;
	fSxy= 0;
	fPixIdmax= -1;
	fPixIdmin= -1;
	fWx0= 0;
	fWy0= 0;
	fX0= 0;
	fY0= 0;
	fXmin= 1.e+99;
	fXmax= -1.e+99;
	fYmin= 1.e+99;
	fYmax= -1.e+99;
	fIx_min= 1.e+99;
	fIx_max= -1.e+99;
	fIy_min= 1.e+99;
	fIy_max= -1.e+99;
	fColorSum.SetXYZ(0,0,0);
	fColor.SetXYZ(0,0,0);
	m_HasStats= false;
}


void Region::UpdateMoments(Pixel pixel){

	//Update moments
	double w= pixel.S;
	double x= pixel.x;	
	double y= pixel.y;
	int ix= pixel.ix;
	int iy= pixel.iy;
	int id= pixel.id;
	double S_edge= pixel.S_edge;
	double S_curv= pixel.S_curv;
	TVector3 col= pixel.color;

	if(w<fSmin) {
		fSmin= w;
		fPixIdmin= id;
	}
	if(w>fSmax) {
		fSmax= w;
		fPixIdmax= id;
	}
	fS+= w;
	fSxx+= w*x*x;
	fSyy+= w*y*y;
	fSxy+= w*x*y;
	fS_curv+= S_curv;
	fS_edge+= S_edge;

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

	fColorSum+= col;

	double delta_curv = S_curv - fM1_curv;
  double delta_curv_n = delta_curv/fNPix;
  double delta_curv_n2 = delta_curv_n * delta_curv_n;
  double f_curv = delta_curv * delta_curv_n * (fNPix-1);
  fM1_curv+= delta_curv_n;
  fM2_curv+= f_curv;	

}//close UpdateMoments()


void Region::AddPixel(Pixel pixel){

	//Append pixel to list
	fPixelCollection.push_back(pixel);

	//Update moment counts
	UpdateMoments(pixel);

}//close AddPixel()


void Region::AddRegion(Region* aRegion,bool addPixels){

	if(!aRegion) {
		cerr<<"Region::AddRegion(): ERROR: Null prt to given region, nothing will be added!"<<endl;
		return;
	}
	//cout<<"Region::AddRegion(): INFO: Adding current region (id="<<fId<<") with region "<<aRegion->fId<<" ..."<<endl;
		

	//Update moments for the region
	double S_B= aRegion->fS;
	double Scurv_B= aRegion->fS_curv;
	double Sedge_B= aRegion->fS_edge;
	double Smin_B= aRegion->fSmin;
	double Smax_B= aRegion->fSmax;
	int PixIdmin_B= aRegion->fPixIdmin;
	int PixIdmax_B= aRegion->fPixIdmax;
	double Sxx_B= aRegion->fSxx;
	double Syy_B= aRegion->fSyy;
	double Sxy_B= aRegion->fSxy;
	double Wx0_B= aRegion->fWx0;
	double Wy0_B= aRegion->fWy0;
	double X0_B= aRegion->fX0;
	double Y0_B= aRegion->fY0;
	int Xmin_B= aRegion->fXmin;
	int Xmax_B= aRegion->fXmax;
	int Ymin_B= aRegion->fYmin;
	int Ymax_B= aRegion->fYmax;	
	int Ixmin_B= aRegion->fIx_min;
	int Ixmax_B= aRegion->fIx_max;
	int Iymin_B= aRegion->fIy_min;
	int Iymax_B= aRegion->fIy_max;
	double N_B= aRegion->fNPix;
	double M1_B= aRegion->fM1;
	double M2_B= aRegion->fM2;
	double M3_B= aRegion->fM3;
	double M4_B= aRegion->fM4;
	double H_B= aRegion->fH;
	double M1Curv_B= aRegion->fM1_curv;
	double M2Curv_B= aRegion->fM2_curv;
	
	TVector3 ColorSum_B= aRegion->fColorSum;
	
	double S_A= this->fS;
	double Scurv_A= this->fS_curv;
	double Sedge_A= this->fS_edge;
	double Smin_A= this->fSmin;
	double Smax_A= this->fSmax;
	int PixIdmin_A= this->fPixIdmin;
	int PixIdmax_A= this->fPixIdmax;
	double Sxx_A= this->fSxx;
	double Syy_A= this->fSyy;
	double Sxy_A= this->fSxy;
	double Wx0_A= this->fWx0;
	double Wy0_A= this->fWy0;
	double X0_A= this->fX0;
	double Y0_A= this->fY0;
	double Xmin_A= this->fXmin;
	double Xmax_A= this->fXmax;
	double Ymin_A= this->fYmin;
	double Ymax_A= this->fYmax;
	int Ixmin_A= this->fIx_min;
	int Ixmax_A= this->fIx_max;
	int Iymin_A= this->fIy_min;
	int Iymax_A= this->fIy_max;
	double N_A= this->fNPix;
	double M1_A= this->fM1;
	double M2_A= this->fM2;
	double M3_A= this->fM3;
	double M4_A= this->fM4;
	double H_A= this->fH;
	double M1Curv_A= this->fM1_curv;
	double M2Curv_A= this->fM2_curv;
	
	TVector3 ColorSum_A= this->fColorSum;
	
	double delta= (N_A*M1_B-N_B*M1_A)/(N_A*N_B);

	fNPix= N_A + N_B;
	fM1= (N_A*M1_A + N_B*M1_B)/fNPix;
	fM2= M2_A + M2_B + pow(N_B*M1_A-N_A*M1_B,2)/(N_A*N_B*fNPix);
	fM3= M3_A + M3_B + pow(delta,3)*N_A*N_B*(N_A-N_B)/(fNPix*fNPix) + 3*delta*(N_A*M2_B-N_B*M2_A)/fNPix;
	fM4= M4_A + M4_B + pow(delta,4)*N_A*N_B*(N_A*N_A-N_A*N_B+N_B*N_B)/pow(fNPix,3) + 6*pow(delta,2)*(N_A*N_A*M2_B+N_B*N_B*M2_A)/pow(fNPix,2) + 4*delta*(N_A*M3_B-N_B*M3_A)/fNPix;

	fM1_curv= (N_A*M1Curv_A + N_B*M1Curv_B)/fNPix;
	fM2_curv= M2Curv_A + M2Curv_B + pow(N_B*M1Curv_A-N_A*M1Curv_B,2)/(N_A*N_B*fNPix);

	if(Smin_A<Smin_A) {
		fSmin= Smin_A;
		fPixIdmin= PixIdmin_A;
	}
	else{
		fSmin= Smin_B;
		fPixIdmin= PixIdmin_B;
	}

	if(Smax_A>Smax_A) {
		fSmax= Smax_A;
		fPixIdmax= PixIdmax_A;
	}
	else{
		fSmax= Smax_B;
		fPixIdmax= PixIdmax_B;
	}

	fS= S_A + S_B;
	fS_curv= Scurv_A + Scurv_B;
	fS_edge= Sedge_A + Sedge_B;
	fSxx= Sxx_A + Sxx_B;
	fSyy= Syy_A + Syy_B;
	fSxy= Sxy_A + Sxy_B;
	fWx0= Wx0_A + Wx0_B;
	fWy0= Wy0_A + Wy0_B;
	fX0= X0_A + X0_B;
	fY0= Y0_A + Y0_B;

	if(Xmin_A<Xmin_B) fXmin= Xmin_A;
	else fXmin= Xmin_B;

	if(Xmax_A>Xmax_B) fXmax= Xmax_A;
	else fXmax= Xmax_B;

	if(Ymin_A<Ymin_B) fYmin= Ymin_A;
	else fYmin= Ymin_B;

	if(Ymax_A>Ymax_B) fYmax= Ymax_A;
	else fYmax= Ymax_B;

	if(Ixmin_A<Ixmin_B) fIx_min= Ixmin_A;
	else fIx_min= Ixmin_B;

	if(Ixmax_A>Ixmax_B) fIx_max= Ixmax_A;
	else fIx_max= Ixmax_B;

	if(Iymin_A<Iymin_B) fIy_min= Iymin_A;
	else fIy_min= Iymin_B;

	if(Iymax_A>Iymax_B) fIy_max= Iymax_A;
	else fIy_max= Iymax_B;

	fColorSum= ColorSum_A + ColorSum_B; 
	fH= H_A + H_B;
	
	//Add new pixels at list end
	if(addPixels){
		std::vector<Pixel>::iterator it= fPixelCollection.end();
		fPixelCollection.insert(it,(aRegion->fPixelCollection).begin(),(aRegion->fPixelCollection).end());
	}

	//Append subregion ids 
	int Id_B= aRegion->fId;
	std::vector<int> SubRegionIds_B= aRegion->fSubRegionIds;
	fSubRegionIds.push_back(Id_B);
	std::vector<int>::iterator it2= fSubRegionIds.end();
	fSubRegionIds.insert(it2,SubRegionIds_B.begin(),SubRegionIds_B.end());

}//close Region::AddRegion()


void Region::SortNeighborInfo(bool isAscending,bool useTotDistance,double beta){

	if(fNeighbourRegionInfo.size()<=0) return;

	if(useTotDistance){
		std::sort(fNeighbourRegionInfo.begin(), fNeighbourRegionInfo.end(), sorter(isAscending,beta) );
	}
	else{
		if(isAscending) std::sort(fNeighbourRegionInfo.begin(), fNeighbourRegionInfo.end(), compareNeighborByAscendingDist);
		else std::sort(fNeighbourRegionInfo.begin(), fNeighbourRegionInfo.end(), compareNeighborByDescendingDist);
	}

}//close Region::SortNeighborInfo()





double Region::GetDissimilarity(Region* aRegion,bool useRobustParams,bool normalizeParams,bool addCurvDist){

	double diss= 1.e+99;
	std::pair<double,double> dists= GetDistance(aRegion,useRobustParams,normalizeParams,addCurvDist);
	double dist_color= dists.first;
	double dist_space= dists.second;
	diss= dist_color/(1.+dist_space);

	return diss;

}//close GetDissimilarity()


std::pair<double,double> Region::GetDistance(Region* aRegion,bool useRobustParams,bool normalizeParams,bool addCurvDist){

	double dist_color= 1.e+99;
	double dist_space= 1.e+99;
	std::pair<double,double> dists= std::make_pair(dist_color,dist_space);
	if(!aRegion){
		cerr<<"Region::GetDistance(): ERROR: Null ptr to given region...returning inf dists!"<<endl;
		return dists;
	}

	/*
	double w_this= this->fMean;		
	double wmin= fImageMinS;
	double wmax= fImageMaxS;	
	double NormMean= NormMin + (NormMax-NormMin)*(this->fMean-fImageMinS)/(fImageMaxS-fImageMinS);
	double NormLogMean= NormMin + log10(NormMean/NormMin)/log10(NormMax/NormMin) * (NormMax-NormMin);
	double NormMean_region= NormMin + (NormMax-NormMin)*(aRegion->fMean-fImageMinS)/(fImageMaxS-fImageMinS);
	double NormLogMean_region= NormMin + log10(NormMean_region/NormMin)/log10(NormMax/NormMin) * (NormMax-NormMin);
	double MeanDiff= NormMean-NormMean_region;
	double MeanLogDiff= NormLogMean-NormLogMean_region;
	*/

	//## Normalization of parameters for euclidean distance computation
	//## Under linear transformation A+Bx --> mean'=A+B*mean, rms'=B*rms
	double NormMin= 0;//1;
	double NormMax= 1;//256;
	
	double A= NormMin - (NormMax-NormMin)*fImageMinS/(fImageMaxS-fImageMinS);
	double B= (NormMax-NormMin)/(fImageMaxS-fImageMinS);
	double Acurv= NormMin - (NormMax-NormMin)*fImageMinScurv/(fImageMaxScurv-fImageMinScurv);	
	double Bcurv= (NormMax-NormMin)/(fImageMaxScurv-fImageMinScurv);

	double Ax0= NormMin - (NormMax-NormMin)*fImageMinX/(fImageMaxX-fImageMinX);
	double Bx0= (NormMax-NormMin)/(fImageMaxX-fImageMinX);
	double Ay0= NormMin - (NormMax-NormMin)*fImageMinY/(fImageMaxY-fImageMinY);
	double By0= (NormMax-NormMin)/(fImageMaxY-fImageMinY);

	//Color difference
	double colorDiff2= (this->fColor - aRegion->fColor).Mag2();
	
	//Mean difference
	double mean= this->fMean;
	double mean_norm= A+B*mean;//mean under linear transformation
	double meanN= aRegion->fMean;
	double meanN_norm= A+B*meanN;
	double meanDiff= mean-meanN;
	double meanDiff_norm= mean_norm-meanN_norm;
	
	//RMS difference
	double rms= this->fRMS;
	double rms_norm= B*rms;
	double rmsN= aRegion->fRMS;
	double rmsN_norm= B*rmsN;
	double rmsDiff= rms-rmsN;
	double rmsDiff_norm= rms_norm-rmsN_norm;
	
	//Median difference
	double median= this->fMedian;
	double median_norm= A+B*median;
	double medianN= aRegion->fMedian;
	double medianN_norm= A+B*medianN;
	double medianDiff= median-medianN;
	double medianDiff_norm= median_norm-medianN_norm;

	//Median MAD difference
	double mad= this->fMedianRMS;
	double mad_norm= B*mad; 
	double madN= aRegion->fMedianRMS;
	double madN_norm= B*madN;
	double madDiff= mad-madN;
	double madDiff_norm= mad_norm-madN_norm;
	
	//Curvature mean
	double curvMean= this->fMean_curv;
	double curvMean_norm= Acurv+Bcurv*curvMean;
	double curvMeanN= aRegion->fMean_curv;
	double curvMeanN_norm= Acurv+Bcurv*curvMeanN;
	double curvMeanDiff= curvMean-curvMeanN;
	double curvMeanDiff_norm= curvMean_norm-curvMeanN_norm;

	//Curvature rms
	double curvRMS= this->fRMS_curv;
	double curvRMS_norm= Bcurv*curvRMS;
	double curvRMSN= aRegion->fRMS_curv;
	double curvRMSN_norm= Bcurv*curvRMSN;
	double curvRMSDiff= curvRMS-curvRMSN;
	double curvRMSDiff_norm= curvRMS_norm-curvRMSN_norm;

	//Centroid distance
	double X0= this->fX0;
 	double Y0= this->fY0;
 	double X0_norm= Ax0+Bx0*X0; 
	double Y0_norm= Ay0+By0*Y0; 
	double X0N= aRegion->fX0;
 	double Y0N= aRegion->fY0;
	double X0N_norm= Ax0+Bx0*X0N;
	double Y0N_norm= Ay0+By0*Y0N;
	double centroidDiffX= X0-X0N;
	double centroidDiffY= Y0-Y0N;
	double centroidDiffX_norm= X0_norm-X0N_norm;
	double centroidDiffY_norm= Y0_norm-Y0N_norm;

	double centroidDistMin= 0;
	double centroidDistMax= sqrt(pow(fImageMaxX-fImageMinX,2) + pow(fImageMaxY-fImageMinY,2));
	double centroidDist= sqrt(pow(X0-X0N,2) + pow(Y0-Y0N,2));
	double centroidDistNorm= NormMin + (NormMax-NormMin)*(centroidDist-centroidDistMin)/(centroidDistMax-centroidDistMin);
	

	//## Compute Euclidean distance
	double dist2= 0;
	double dist2_space= 0;
	if(normalizeParams){
		if(useRobustParams){
			dist2+= medianDiff_norm*medianDiff_norm;//MEDIAN
			dist2+= madDiff_norm*madDiff_norm;//MEDIAN MAD
		}
		else{
			dist2+= meanDiff_norm*meanDiff_norm;//MEAN
			dist2+= rmsDiff_norm*rmsDiff_norm;//RMS
		}

		if(addCurvDist){
			dist2+= curvMeanDiff_norm*curvMeanDiff_norm;//MEAN Curvature
			dist2+= curvRMSDiff_norm*curvRMSDiff_norm;//RMS Curvature
		}

		//dist2_space+= centroidDiffX_norm*centroidDiffX_norm + centroidDiffY_norm*centroidDiffY_norm;//CENTROIDS
		dist2_space+= centroidDistNorm*centroidDistNorm;
		
	}//close if normalize pars
	else{
		if(useRobustParams){
			dist2+= medianDiff*medianDiff;//MEDIAN
			dist2+= madDiff*madDiff;//MEDIAN MAD
		}
		else{
			dist2+= meanDiff*meanDiff;//MEAN
			dist2+= rmsDiff*rmsDiff;//RMS
		}

		if(addCurvDist){
			dist2+= curvMeanDiff*curvMeanDiff;//MEAN Curvature
			dist2+= curvRMSDiff*curvRMSDiff;//RMS Curvature
		}

		//dist2_space+= centroidDiffX*centroidDiffX + centroidDiffY*centroidDiffY;//CENTROIDS
		dist2_space+= centroidDist*centroidDist;
	}

	dist_color= sqrt(dist2);		
	dist_space= sqrt(dist2_space);
	
	return std::make_pair(dist_color,dist_space);	
	
}//close Region::GetDistance()


Region::Parameters* Region::GetParams(bool includeCurvPar){

	int nPars= 2;
	int nPars_robust= 2;//3
	int nPars_spatial= 2;
	if(includeCurvPar){
		nPars+= 2;
		nPars_robust+= 2;
	}
	
	Region::Parameters* params= new Region::Parameters;
	params->pars= new TVectorD(nPars);
	params->robustPars= new TVectorD(nPars_robust);
	params->spatialPars= new TVectorD(nPars_spatial);

	//== Appearance pars ==
	(*(params->pars))(0)= fMean;
	(*(params->pars))(1)= fRMS;
	
	//== Robust appearance pars ==
	(*(params->robustPars))(0)= fMedian;
	(*(params->robustPars))(1)= fMedianRMS;
	//(*(params->robustPars))(2)= fH;

	if(includeCurvPar){
		(*(params->pars))(2)= fMean_curv;
		(*(params->pars))(3)= fRMS_curv;

		(*(params->robustPars))(2)= fMedian_curv;
		(*(params->robustPars))(3)= fMedianRMS_curv;
	}

	//== Spatial pars ==
	(*(params->spatialPars))(0)= fX0;
	(*(params->spatialPars))(1)= fY0;

	return params;

}//close GetParams()


TVectorD Region::GetParamVector(bool includeSpatialPar,bool includeCurvPar, bool useRobustPars){

	//Compute norm pars
	double NormMin= 0;
	double NormMax= 1;
	double A= NormMin - (NormMax-NormMin)*fImageMinS/(fImageMaxS-fImageMinS);
	double B= (NormMax-NormMin)/(fImageMaxS-fImageMinS);
	double Acurv= NormMin - (NormMax-NormMin)*fImageMinScurv/(fImageMaxScurv-fImageMinScurv);	
	double Bcurv= (NormMax-NormMin)/(fImageMaxScurv-fImageMinScurv);
	double Ax0= NormMin - (NormMax-NormMin)*fImageMinX/(fImageMaxX-fImageMinX);
	double Bx0= (NormMax-NormMin)/(fImageMaxX-fImageMinX);
	double Ay0= NormMin - (NormMax-NormMin)*fImageMinY/(fImageMaxY-fImageMinY);
	double By0= (NormMax-NormMin)/(fImageMaxY-fImageMinY);
	
	double MeanNorm= A+B*fMean;
	double RMSNorm= B*fRMS;
	double MeanCurvNorm= A+B*fMean_curv;
	double RMSCurvNorm= B*fRMS_curv;

	double MedianNorm= A+B*fMedian;
	double MedianRMSNorm= B*fMedianRMS;
	double MedianCurvNorm= Acurv+Bcurv*fMedian_curv;
	double MedianRMSCurvNorm= Bcurv*fMedianRMS_curv;

	double X0= this->fX0;
 	double Y0= this->fY0;
 	double X0Norm= Ax0+Bx0*X0; 
	double Y0Norm= Ay0+By0*Y0; 

	int nPars= 2;
	if(includeCurvPar) nPars+= 2;
	if(includeSpatialPar) nPars+= 2;
	
	TVectorD v(nPars);
	if(useRobustPars){
		v(0)= MedianNorm;
		v(1)= MedianRMSNorm;
	}
	else{
		v(0)= MeanNorm;
		v(1)= RMSNorm;
	}

	if(includeCurvPar){
		if(useRobustPars){
			v(2)= MedianCurvNorm;
			v(3)= MedianRMSCurvNorm;
		}
		else{
			v(2)= MeanCurvNorm;
			v(3)= RMSCurvNorm;
		}
		if(includeSpatialPar){
			v(4)= X0Norm;
			v(5)= Y0Norm;
		}
	}
	else{
		if(includeSpatialPar){
			v(2)= X0Norm;
			v(3)= Y0Norm;
		}
	}
 
	return v;

}//close Region::GetParamVector()

std::pair<double,double> Region::GetAsymmDistance(Region* aRegion,bool useRobustParams,bool normalizeParams,bool addSpatialDist,bool addCurvDist){

	std::chrono::high_resolution_clock::time_point t1, t2;

	t1 = high_resolution_clock::now();
	
	double dist= 1.e+99;
	double dist_neighbor= 1.e+99;
	std::pair<double,double> dists= std::make_pair(dist,dist_neighbor);
	if(!aRegion){
		cerr<<"Region::GetAsymmDistance(): ERROR: Null ptr to given region...returning inf dists!"<<endl;
		return dists;
	}
	
	int regionId= fId;
	int regionId_neighbor= aRegion->fId;

	if(useRobustParams){
		//Create a merged region between current and neighbor	
		Region mergedRegion= *(this);
		mergedRegion.AddRegion(aRegion,useRobustParams);
		mergedRegion.ComputeParameters(false,true,false);
		
		//Find distance between this region and merged region
		//dist= this->GetDistance(&mergedRegion,useRobustParams,normalizeParams);
		std::pair<double,double> dists= this->GetDistance(&mergedRegion,useRobustParams,normalizeParams,addCurvDist);
		dist= dists.first;

		//Find distance between neighbor region and merged region
		//dist_neighbor= aRegion->GetDistance(&mergedRegion,useRobustParams,normalizeParams);
		std::pair<double,double> dists_neighbor= aRegion->GetDistance(&mergedRegion,useRobustParams,normalizeParams,addCurvDist);
		dist_neighbor= dists_neighbor.first;

	}//close if
	else{//fast computing for mean & rms params
		double N_A= (double)fNPix;
		double N_B= (double)aRegion->fNPix; 
		double N= (double)N_A + (double)N_B;
		double M1_A= fM1;
		double M1_B= aRegion->fM1;
		double M2_A= fM2;
		double M2_B= aRegion->fM2;
		double M1Curv_A= fM1_curv;
		double M1Curv_B= aRegion->fM1_curv;
		double M2Curv_A= fM2_curv;
		double M2Curv_B= aRegion->fM2_curv;
		double Mean_A= fMean;
		double Mean_B= aRegion->fMean;
		double RMS_A= fRMS;
		double RMS_B= aRegion->fRMS;
		double MeanCurv_A= fMean_curv;
		double MeanCurv_B= aRegion->fMean_curv;
		double RMSCurv_A= fRMS_curv;
		double RMSCurv_B= aRegion->fRMS_curv;

		double S_A= fS;
		double Cx_A= fWx0;
		double Cy_A= fWy0;
		double S_B= aRegion->fS;
		double Cx_B= aRegion->fWx0;
		double Cy_B= aRegion->fWy0;
		double Cx= Cx_A*S_A/(S_A+S_B) + Cx_B*S_B/(S_A+S_B);
		double Cy= Cy_A*S_A/(S_A+S_B) + Cy_B*S_B/(S_A+S_B);

		double M1= (N_A*M1_A + N_B*M1_B)/N;
		double M2= M2_A + M2_B + pow(N_B*M1_A-N_A*M1_B,2)/(N_A*N_B*N);	
		double M1Curv= (N_A*M1Curv_A + N_B*M1Curv_B)/N;
		double M2Curv= M2Curv_A + M2Curv_B + pow(N_B*M1Curv_A-N_A*M1Curv_B,2)/(N_A*N_B*N);	
		double Mean= M1; 
		double MeanCurv= M1Curv; 
		double RMS= 0;
		double RMSCurv= 0;
		if(N>1) {
			RMS= sqrt(M2/(N-1));
			RMSCurv= sqrt(M2Curv/(N-1));
		}

		//Compute normalized params
		double NormMin= 0;
		double NormMax= 1;
		double A= NormMin - (NormMax-NormMin)*fImageMinS/(fImageMaxS-fImageMinS);
		double B= (NormMax-NormMin)/(fImageMaxS-fImageMinS);
		double Acurv= NormMin - (NormMax-NormMin)*fImageMinScurv/(fImageMaxScurv-fImageMinScurv);	
		double Bcurv= (NormMax-NormMin)/(fImageMaxScurv-fImageMinScurv);
		double Ax0= NormMin - (NormMax-NormMin)*fImageMinX/(fImageMaxX-fImageMinX);
		double Bx0= (NormMax-NormMin)/(fImageMaxX-fImageMinX);
		double Ay0= NormMin - (NormMax-NormMin)*fImageMinY/(fImageMaxY-fImageMinY);
		double By0= (NormMax-NormMin)/(fImageMaxY-fImageMinY);
		double MeanNorm_A= A+B*Mean_A;
		double MeanNorm_B= A+B*Mean_B;	
		double MeanNorm= A+B*Mean;
		double RMSNorm_A= B*RMS_A;
		double RMSNorm_B= B*RMS_B;	
		double RMSNorm= B*RMS;
		double MeanCurvNorm_A= A+B*MeanCurv_A;
		double MeanCurvNorm_B= A+B*MeanCurv_B;	
		double MeanCurvNorm= A+B*MeanCurv;
		double RMSCurvNorm_A= B*RMSCurv_A;
		double RMSCurvNorm_B= B*RMSCurv_B;	
		double RMSCurvNorm= B*RMSCurv;

		double CxNorm_A= Ax0+Bx0*Cx_A;
		double CyNorm_A= Ay0+By0*Cy_A;
		double CxNorm_B= Ax0+Bx0*Cx_B;
		double CyNorm_B= Ay0+By0*Cy_B;
		double CxNorm= Ax0+Bx0*Cx;
		double CyNorm= Ay0+By0*Cy;

		double dist2= 0;
		double dist2_neighbor= 0;

		
		if(normalizeParams){
			dist2+= pow(MeanNorm_A-MeanNorm,2);//mean diff
			dist2+= pow(RMSNorm_A-RMSNorm,2);//rms diff
			
			dist2_neighbor+= pow(MeanNorm_B-MeanNorm,2);//mean diff
			dist2_neighbor+= pow(RMSNorm_B-RMSNorm,2);//rms diff

			if(addCurvDist){
				dist2+= pow(MeanCurvNorm_A-MeanCurvNorm,2);//mean curvature diff		
				dist2+= pow(RMSCurvNorm_A-RMSCurvNorm,2);//rms curvature diff		

				dist2_neighbor+= pow(MeanCurvNorm_B-MeanCurvNorm,2);//mean curvature diff		
				dist2_neighbor+= pow(RMSCurvNorm_B-RMSCurvNorm,2);//rms curvature diff
			}
			

			if(addSpatialDist){
				dist2+= pow(CxNorm_A-CxNorm,2);//Cx diff
				dist2+= pow(CyNorm_A-CyNorm,2);//Cy diff

				dist2_neighbor+= pow(CxNorm_B-CxNorm,2);//Cx diff
				dist2_neighbor+= pow(CyNorm_B-CyNorm,2);//Cx diff
			}			

		}
		else{
			dist2+= pow(Mean_A-Mean,2);//mean diff
			dist2+= pow(RMS_A-RMS,2);//rms diff
				
			dist2_neighbor+= pow(Mean_B-Mean,2);//mean diff
			dist2_neighbor+= pow(RMS_B-RMS,2);//rms diff
			
			if(addCurvDist){
				dist2+= pow(MeanCurv_A-MeanCurv,2);//mean curvature diff		
				dist2+= pow(RMSCurv_A-RMSCurv,2);//rms curvature diff		
						
				dist2_neighbor+= pow(MeanCurv_B-MeanCurv,2);//mean curvature diff		
				dist2_neighbor+= pow(RMSCurv_B-RMSCurv,2);//rms curvature diff			
			}

			if(addSpatialDist){
				dist2+= pow(Cx_A-Cx,2);//Cx diff
				dist2+= pow(Cy_A-Cy,2);//Cy diff

				dist2_neighbor+= pow(Cx_B-Cx,2);//Cx diff
				dist2_neighbor+= pow(Cy_B-Cy,2);//Cx diff
			}	
		}

		dist= sqrt(dist2);		
		dist_neighbor= sqrt(dist2_neighbor);

		//cout<<"Region link "<<regionId<<"-"<<regionId_neighbor<<": A="<<A<<" B="<<B<<" Acurv="<<Acurv<<" Bcurv="<<Bcurv<<" N_A="<<N_A<<" N_B="<<N_B<<" M1_A="<<M1_A<<" M1_B="<<M1_B<<" M2_A="<<M2_A<<" M2_B="<<M2_B<<" M1="<<M1<<" M2="<<M2<<" Mean="<<Mean<<" RMS="<<RMS<<" MeanCurv="<<MeanCurv<<" RMSCurv="<<RMSCurv<<" Mean_A="<<Mean_A<<" RMS_A="<<RMS_A<<" MeanCurv_A="<<MeanCurv_A<<" RMSCurv_A="<<RMSCurv_A<<" Mean_B="<<Mean_B<<"  RMS_B="<<RMS_B<<" MeanCurv_B="<<MeanCurv_B<<" RMSCurv_B="<<RMSCurv_B<<" dist2="<<dist2<<" dist2_neighbor="<<dist2_neighbor<<endl;
	
	}//close else
	t2 = high_resolution_clock::now();
	auto dt = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
	
	return std::make_pair(dist,dist_neighbor);

}//close GetAsymmDistance()


void Region::AddNeighborInfo(NeighborInfo aNeighborInfo){
	
	//Check if neighbor id is already present or coincides with this region id
	long int neighborId= aNeighborInfo.id;
	std::vector<NeighborInfo>::iterator it = std::find_if(fNeighbourRegionInfo.begin(), fNeighbourRegionInfo.end(), FindNeighborId(neighborId));
	if ( fId==neighborId || it!= fNeighbourRegionInfo.end() ) {//already present in the list
		return;
	}

	//Add neighbor info to the list
	fNeighbourRegionInfo.push_back(aNeighborInfo);

}//close Region::AddNeighborInfo()


int Region::AddNeighborInfo(Region* aNeighborRegion,int order,bool useRobustParams,bool normalizeParams,bool addCurvDist){

	//##############################################################
	//##   ADD REGION TO THE NEIGHBOR LIST
	//##   To speed-up add also this region to the passed region
	//##############################################################
	std::chrono::high_resolution_clock::time_point tstart, tstop, t1, t2;
	tstart= high_resolution_clock::now();
	
	//Check given neighbor region
	if(!aNeighborRegion) {
		cerr<<"Region::AddNeighborInfo(): ERROR: Null ptr to given neighbor region!"<<endl;
		return -1;
	}

	//Check if neighbor id is already present or coincides with this region id
	t1 = high_resolution_clock::now();
	long int neighborId= aNeighborRegion->fId;
	std::vector<NeighborInfo>::iterator it = std::find_if(fNeighbourRegionInfo.begin(), fNeighbourRegionInfo.end(), FindNeighborId(neighborId));
	if ( fId==neighborId || it!= fNeighbourRegionInfo.end() ) {
		//cout<<"Region::AddNeighborInfo(): INFO: Region "<<fId<<": Neighbor region "<<neighborId<<" already present in the list!"<<endl;
		return 0;//already present in the list
	}
	t2 = high_resolution_clock::now();
	auto dt_search = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();

	double dist= 0;
	double dist_neighbor= 0;
	t1 = high_resolution_clock::now();

	if(useRobustParams){
		//Create a merged region between current and neighbor	
		Region mergedRegion= *(this);
		
		mergedRegion.AddRegion(aNeighborRegion,useRobustParams);
		mergedRegion.ComputeParameters(false,true,false);
		
		//Find distance between this region and merged region
		std::pair<double,double> dists= this->GetDistance(&mergedRegion,useRobustParams,normalizeParams,addCurvDist);
		dist= dists.first;
		
		//Find distance between neighbor region and merged region
		std::pair<double,double> dists_neighbor= aNeighborRegion->GetDistance(&mergedRegion,useRobustParams,normalizeParams,addCurvDist);
		dist_neighbor= dists_neighbor.first;
		
	}//close if
	else{//fast computing for mean & rms params
		double N_A= fNPix;
		double N_B= aNeighborRegion->fNPix; 
		double N= (double)N_A + (double)N_B;
		double M1_A= fM1;
		double M1_B= aNeighborRegion->fM1;
		double M2_A= fM2;
		double M2_B= aNeighborRegion->fM2;
		double M1Curv_A= fM1_curv;
		double M1Curv_B= aNeighborRegion->fM1_curv;
		double M2Curv_A= fM2_curv;
		double M2Curv_B= aNeighborRegion->fM2_curv;
		double Mean_A= fMean;
		double Mean_B= aNeighborRegion->fMean;
		double RMS_A= fRMS;
		double RMS_B= aNeighborRegion->fRMS;
		double MeanCurv_A= fMean_curv;
		double MeanCurv_B= aNeighborRegion->fMean_curv;
		double RMSCurv_A= fRMS_curv;
		double RMSCurv_B= aNeighborRegion->fRMS_curv;
		double M1= (N_A*M1_A + N_B*M1_B)/N;
		double M2= M2_A + M2_B + pow(N_B*M1_A-N_A*M1_B,2)/(N_A*N_B*N);	
		double M1Curv= (N_A*M1Curv_A + N_B*M1Curv_B)/N;
		double M2Curv= M2Curv_A + M2Curv_B + pow(N_B*M1Curv_A-N_A*M1Curv_B,2)/(N_A*N_B*N);	
		double Mean= M1; 
		double MeanCurv= M1Curv; 
		double RMS= 0;
		double RMSCurv= 0;
		if(N>1) {
			RMS= sqrt(M2/(N-1));
			RMSCurv= sqrt(M2Curv/(N-1));
		}
		double dist2= 0;
		dist2+= pow(Mean_A-Mean,2);//mean diff
		dist2+= pow(RMS_A-RMS,2);//rms diff

		double dist2_neighbor= 0;
		dist2_neighbor+= pow(Mean_B-Mean,2);//mean diff
		dist2_neighbor+= pow(RMS_B-RMS,2);//rms diff

		if(addCurvDist){
			dist2+= pow(MeanCurv_A-MeanCurv,2);//mean curvature diff		
			dist2+= pow(RMSCurv_A-RMSCurv,2);//rms curvature diff		

			dist2_neighbor+= pow(MeanCurv_B-MeanCurv,2);//mean curvature diff		
			dist2_neighbor+= pow(RMSCurv_B-RMSCurv,2);//rms curvature diff		
		}
		
		dist= sqrt(dist2);
		dist_neighbor= sqrt(dist2_neighbor);
	}//close else
	t2 = high_resolution_clock::now();
	auto dt_params = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();

	//Compute edgeness on border pixels only for 1st order neighbors
	t1 = high_resolution_clock::now();
	double Edgeness= 1;//fixed penalty for 2-nd order neighbors
	double Edgeness_neighbor= 1;
	
	t2 = high_resolution_clock::now();
	auto dt_computeEdgeness = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
	

	t1 = high_resolution_clock::now();
	//Append data to neighbor info list
	NeighborInfo info;
	info.id= neighborId;
	info.D= dist;
	info.order= order;
	//info.H= mergedRegion.fH;
	info.Edgeness= Edgeness;//TO BE ESTIMATED BETTER!
	fNeighbourRegionInfo.push_back(info);

	//Append data to region neighbors
	NeighborInfo info_neighbor;
	info_neighbor.id= fId;
	info_neighbor.D= dist_neighbor;
	info_neighbor.order= order;
	//info_neighbor.H= mergedRegion.fH;
	info_neighbor.Edgeness= Edgeness_neighbor;
	(aNeighborRegion->fNeighbourRegionInfo).push_back(info_neighbor);
	t2 = high_resolution_clock::now();
	auto dt_addInfo = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();


	tstop= high_resolution_clock::now();
	double dt_tot= std::chrono::duration_cast<std::chrono::microseconds>(tstop-tstart).count();

	//cout<<"Region::AddNeighborInfo(): INFO: dt_tot="<<dt_tot<<" dt_search="<<dt_search/dt_tot*100<<" dt_copy="<<dt_copyMergedRegion/dt_tot*100<<" dt_params="<<dt_params/dt_tot*100<<" dt_dist="<<dt_computeDistance/dt_tot*100<<" dt_edgeness="<<dt_computeEdgeness/dt_tot*100<<" dt_addInfo="<<dt_addInfo/dt_tot*100<<endl;
	//cout<<"Region::AddNeighborInfo(): INFO: dt_tot="<<dt_tot<<" dt_search="<<dt_search/dt_tot*100<<" dt_params="<<dt_params/dt_tot*100<<" dt_edgeness="<<dt_computeEdgeness/dt_tot*100<<" dt_addInfo="<<dt_addInfo/dt_tot*100<<endl;
	
	return 0;

}//close AddNeighborInfo()


Source* Region::GetSource(){

	Source* aSource= new Source;
	aSource->SetId(fId);	
	aSource->SetName(Form("S%d",fId));
	aSource->SetType(Source::eExtended);

	//bool fIsSignificant;
	//bool fIsSalient;

	//Fill pixels
	for(unsigned int l=0;l<fPixelCollection.size();l++){
		Region::Pixel pixel= fPixelCollection[l];
		
		Source::Pixel aPixel;
		aPixel.S= pixel.S;
		aPixel.Type= Source::eNormal;
		aPixel.Z= 0;
		aPixel.Curv= pixel.S_curv;
		aPixel.BkgLevel= 0;
		aPixel.id= pixel.id;
		aPixel.ix= pixel.ix;
		aPixel.iy= pixel.iy;
		aPixel.x= pixel.x;
		aPixel.y= pixel.y;
		if( pixel.x<=fImageMinX || pixel.y<=fImageMinY || pixel.x>=fImageMaxX || pixel.y>=fImageMaxY) 
			aSource->SetEdgeFlag(true);

		aSource->AddPixel(aPixel);
	}//end loop cluster pixels
	
	//## Compute stats parameters
	cout<<"Region::GetSource(): INFO: Computing source stats..."<<endl;
	aSource->ComputeStats();
	
	//## Compute morphology parameters
	cout<<"Region::GetSource(): INFO: Computing morphology params..."<<endl;
	aSource->ComputeMorphologyParams();

	return aSource;

}//close GetSource()


int Region::ComputeParameters(bool computeContours,bool computeRobustStats,bool forceRecomputing){

	//## Compute region stats & shape parameters
	if(fNPix<=0) return -1;

	//## Recomputing moments?
	//cout<<"Npix (before force recomputing)="<<fNPix<<" X0="<<fX0<<", Y0="<<fY0<<endl;
	if(forceRecomputing){
		ResetMoments();//reset moments
		for(unsigned int k=0;k<fPixelCollection.size();k++) UpdateMoments(fPixelCollection[k]);//update moments
		//cout<<"Npix (after moment reset & update)="<<fNPix<<" X0="<<fX0<<", Y0="<<fY0<<endl;
	}
	
	//## NB: Compute stats parameters only if not already done or if forced
	//##     This was a bug in older version affecting saliency calculation (the spatial component in region distance) and eventually the extended source extraction!!!
	//##     Need to introduce a rescaling of spatial vs color distance (a factor 100 more or less?) 
	if(!m_HasStats || forceRecomputing){
		fX0/= (double)(fNPix);
		fY0/= (double)(fNPix);
		fWx0/= fS;
		fWy0/= fS;
		fSxx/= fS;
		fSyy/= fS;
		fSxy/= fS;
		fColor*= 1./(double)(fNPix);
		m_HasStats= true;
	}

	//cout<<"Npix (after force recomputing)="<<fNPix<<" X0="<<fX0<<", Y0="<<fY0<<endl;
	

	fMean= fM1;
	fMean_curv= fM1_curv;
	fRMS= 0;
	fRMS_curv= 0;	
	if(fNPix>1) {
		fRMS= sqrt(fM2/(fNPix-1));
		fRMS_curv= sqrt(fM2_curv/(fNPix-1));
	}

	fSkewness= 0;
	fKurtosis= 0;
	if(fM2!=0) {
  	fSkewness= sqrt(fNPix)*fM3/pow(fM2,1.5);//need to adjust for finite population?
  	fKurtosis= fNPix*fM4/(fM2*fM2);//also given without -3
	}
	
	//## Compute robust stats (median, MAD, Entropy, ...)
	if(computeRobustStats){	
		std::chrono::high_resolution_clock::time_point t1, t2;
		t1 = high_resolution_clock::now();
		
		/*
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

		t2 = high_resolution_clock::now();
		auto dt_median = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();

		//## Compute pixel sample entropy	
		t1 = high_resolution_clock::now();
		double Entropy= -1;
		int status= this->ComputeEntropy(Entropy);
		if(status<0){
			cerr<<"Region::ComputeParameters(): ERROR: Computation of Entropy failed!"<<endl;
		}
		fH= Entropy;
		
		t2 = high_resolution_clock::now();
		auto dt_entropy = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
		cout<<"Region::ComputeParameters(): INFO: dt_median="<<dt_median<<" dt_entropy="<<dt_entropy<<endl;
		*/

		RobustStats rstats;
		int status= ComputeRobustStats(rstats);
		if(status<0) {
			cerr<<"Region::ComputeParameters(): ERROR: Robust stats params computing failed!"<<endl;
			return -1;
		}		
		fH= rstats.entropy;
		fMedian= rstats.median;
		fMedianRMS= rstats.mad;
		fMedian_curv= rstats.median_curv;
		fMedianRMS_curv= rstats.mad_curv;	
		t2 = high_resolution_clock::now();
		auto dt_rstats = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
		//cout<<"Region::ComputeParameters(): INFO: dt_rstats="<<dt_rstats<<endl;
	}//close if computeRobustStats

	if(computeContours){
		//######################################
		//## Find the region bounding box
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

		cv::Mat binarizedImg = cv::Mat::zeros(nBoxIY, nBoxIX, CV_8UC1);

		for(unsigned int k=0;k<fPixelCollection.size();k++){
			Pixel thisPixel= fPixelCollection[k];
			double thisS= thisPixel.S;	
			int thisId= thisPixel.id;
			double thisX= thisPixel.x;
			double thisY= thisPixel.y;
			int ix= thisPixel.ix;
			int iy= thisPixel.iy;
			ix-= boundingBoxIX[0];
			iy-= boundingBoxIY[0];
			int rowId= nBoxIY-1-iy;
			int colId= nBoxIX-1-ix;
			binarizedImg.at<uchar>(rowId, colId, 0) = 1;		
		}//end loop pixels


		//## Compute contour
		cv::Mat binarizedImg_clone= binarizedImg.clone();//image will be modified by findContours! 
		std::vector<std::vector<cv::Point>> contours; // Vector for storing contour
  	std::vector<cv::Vec4i> hierarchy;
		cv::findContours( binarizedImg_clone, contours, hierarchy,CV_RETR_EXTERNAL,CV_CHAIN_APPROX_NONE, cv::Point(0,0) ); // Find only the external contours in the image
  
		for(unsigned int i=0; i<contours.size(); i++){ // iterate through each contour
			int nContourPts= (int)contours[i].size();
			if(nContourPts<=0) continue;

			//Create and fill contour
			fContour= new Contour;

			//cout<<"Region::ComputeParameters(): INFO: Contour no. "<<i+1<<": (";
		
			for(int j=0;j<nContourPts;j++){
				int contx= contours[i][j].x;
				int conty= contours[i][j].y;
				int rowId= nBoxIY-1-conty;
				int colId= nBoxIX-1-contx;
				int contx_transf= colId + boundingBoxX[0];
				int conty_transf= rowId + boundingBoxY[0];
			
				fContour->AddPoint(cv::Point2f(contx_transf,conty_transf));
				//cout<<"("<<contx_transf<<","<<conty_transf<<"), ";
			}//end loop points in contour
			//cout<<")"<<endl;
		
			//Compute contour parameters
			/*
			if(fContour->ComputeParameters()<0){
				cerr<<"Region::ComputeParameters(): WARN: Failed to compute parameters for contour no. "<<i<<"!"<<endl;
				delete fContour;
				continue;
			}
			*/
		
			//Add contour to the list
			fContourCollection.push_back(fContour);	

		}//end loop contours


		//## Tag edge points
		fContourPointIndex.clear();
		fAroundContourPoints.clear();
		//fAroundContourPointsMap.clear();
 
		if(fContourCollection.size()>0){
			bool measureDist= true;
			int dx8[8] = {-1, -1, 0, 1, 1, 1, 0,-1};
			int dy8[8] = { 0, -1,-1,-1, 0, 1, 1, 1};

			TMatrixD isVisited(nBoxIY, nBoxIX);
			isVisited.Zero();

			for(unsigned int k=0;k<fPixelCollection.size();k++){
				double x= fPixelCollection[k].x;
				double y= fPixelCollection[k].y;
				int ix= fPixelCollection[k].ix;
				int iy= fPixelCollection[k].iy;
				double S_curv= fPixelCollection[k].S_curv; 
				double S_edge= fPixelCollection[k].S_edge; 

				//Loop over external contours (for disjoint regions, i.e. those merged, there are more than 1 contour)
				for(unsigned int s=0;s<fContourCollection.size();s++){

					std::vector<cv::Point2f> contour_points= fContourCollection[s]->GetPoints();
					double distanceToEdge= cv::pointPolygonTest(contour_points, cv::Point2f(x,y), measureDist);

					//If outside this contour, skip and test next
					if(distanceToEdge<0){
						continue;
					}

					if(distanceToEdge>0) {//inside contour
						fPixelCollection[k].isOnEdge= false; 
						if(distanceToEdge<fPixelCollection[k].distanceToEdge) 
						fPixelCollection[k].distanceToEdge= distanceToEdge;	
					}//close if inside contour

					else if(distanceToEdge==0) {//on contour edge
						fPixelCollection[k].isOnEdge= true;
						fPixelCollection[k].distanceToEdge= 0;
						fContourPointIndex.push_back(k);
	
						//Fill list of pixels surrounding border
						for (int l=0; l<8; l++) {
      				double x_neighbor = x + dx8[l];
							double y_neighbor = y + dy8[l];
							int ix_neighbor = ix + dx8[l];
							int iy_neighbor = iy + dy8[l];
							int rowId= nBoxIY-1-(iy_neighbor-boundingBoxIY[0]);
							int colId= nBoxIX-1-(ix_neighbor-boundingBoxIX[0]);
							double distanceToNeighbor= cv::pointPolygonTest(contour_points, cv::Point2f(x_neighbor,y_neighbor), false);
		
							EdgeInfo edge_info;
							edge_info.x= x_neighbor;
							edge_info.y= y_neighbor;
							edge_info.labelId= -1;	
							edge_info.S_curv= 0;
							edge_info.S_edge= 0;

							/*
							if(distanceToNeighbor<0 ){//outside contour
								fAroundContourPointsMap[k].push_back( cv::Point2f(x_neighbor,y_neighbor) );
							}
							*/
							if(distanceToNeighbor<0 && isVisited(rowId,colId)==0 ){//outside contour
								//fAroundContourPoints.push_back( cv::Point2f(x_neighbor,y_neighbor) );
								fAroundContourPoints.push_back(edge_info);
								isVisited(rowId,colId)= 1;
							}
						
						}//end loop neighbors
					}//close else if on edge
				
					//Break contour loop for this pixel because we already found the contour edge, so no need to test the other
					break;

				}//end loop contours
			}//end loop pixels
		}//close if
		else{
			cerr<<"Region::ComputeParameters(): WARN: Contour computation for this region failed!"<<endl;
		}
	}//close if compute contours

	return 0;

}//close Region::ComputeParameters()


int Region::ComputeRobustStats(RobustStats& rstats){

	//## Compute Kullback-Leibler divergence & Entropy of region
	int nNearestNeighbors= 5;//used for density estimation
	rstats.entropy= TMath::Infinity();//init to inf
	rstats.median= TMath::Infinity();
	rstats.mad= TMath::Infinity();

	try{
		RInside* fR= RInside::instancePtr();
		if(!fR){
			cerr<<"Region::ComputeRobustStats(): ERROR: Cannot retrieve RInside instance (did you create one in your main application?)!"<<endl;
			fR= new RInside;
		} 

		//## Load R library for outlier detection
		//cout<<"Region::ComputeEntropy(): INFO: Loading needed R packages..."<<endl;
		std::string RCmd= std::string("library(\"FNN\");");
		fR->parseEvalQ(RCmd);

		//## Import data matrix in R environment
		//cout<<"Region::ComputeEntropy(): INFO: Importing data matrix in R environment..."<<endl;
		int N= (int)fPixelCollection.size();
		Rcpp::NumericMatrix X(N,1);	
		Rcpp::NumericMatrix X_curv(N,1);	
		nNearestNeighbors= sqrt(N);
		
		double NormMin= 1;
		double NormMax= 256;

		for(int i=0;i<N;i++){
			double S= fPixelCollection[i].S;		
			double S_curv= fPixelCollection[i].S_curv;
			double SNorm= NormMin + (NormMax-NormMin)*(S-fImageMinS)/(fImageMaxS-fImageMinS);
			//X(i,0)= SNorm;
			X(i,0)= S;
			X_curv(i,0)= S_curv;	
		}

		(*fR)["X"]= X;
		(*fR)["X_curv"]= X_curv;
		(*fR)["nNearestNeighbors"]= nNearestNeighbors;
		
		//## Compute median & mad
		RCmd= std::string("median<-median(X); mad<-mad(X); median_curv<-median(X_curv); mad_curv<-mad(X_curv);");
		fR->parseEval(RCmd);
		rstats.median= fR->parseEval(std::string("median"));
		rstats.mad= fR->parseEval(std::string("mad"));
		rstats.median_curv= fR->parseEval(std::string("median_curv"));
		rstats.mad_curv= fR->parseEval(std::string("mad_curv"));

		//### Compute entropy
		if(nNearestNeighbors<N){
			RCmd= std::string("E<-entropy(X=X, k=nNearestNeighbors, \"kd_tree\");");
			fR->parseEval(RCmd);
			Rcpp::NumericVector Entr = fR->parseEval(std::string("E"));
			rstats.entropy= Entr(nNearestNeighbors-1);
		}

	}//close try block
	catch( std::exception &ex ) {
		cerr << "Region::ComputeRobustStats(): ERROR: Exception catched: " << ex.what() << endl;
		return -1;
  } 
	catch(...) { 
		cerr << "Region::ComputeRobustStats(): ERROR: C++ exception (unknown reason)" << endl;
		return -1;
  }	

	return 0;

}//close ComputeRobustStats()


std::vector<double> Region::GetEntropyDistance(Region* aRegion){

	std::vector<double> entropyPars;
	entropyPars.clear();
	if(!aRegion) return entropyPars;

	int nNearestNeighbors= 5;//used for density estimation
	try{
		RInside* fR= RInside::instancePtr();
		if(!fR){
			cerr<<"Region::GetEntropyDistance(): ERROR: Cannot retrieve RInside instance (did you create one in your main application?)!"<<endl;
			fR= new RInside;
		} 

		//## Load R library for outlier detection
		std::string RCmd= std::string("library(\"FNN\");");
		fR->parseEvalQ(RCmd);

		//## Import data matrix in R environment
		int N1= (int)fPixelCollection.size();
		Rcpp::NumericMatrix X1(N1,1);	
		int N2= (int)aRegion->fPixelCollection.size();
		Rcpp::NumericMatrix X2(N2,1);
		nNearestNeighbors= sqrt(std::min(N1,N2));	
	
		double NormMin= 1;
		double NormMax= 256;

		//cout<<"N1="<<N1<<" N2="<<N2<<endl;
		//cout<<"X1<-rbind(";
		for(int i=0;i<N1;i++){
			double S= fPixelCollection[i].S;
			double SNorm= NormMin + (NormMax-NormMin)*(S-fImageMinS)/(fImageMaxS-fImageMinS);
			//cout<<S<<",";
			X1(i,0)= SNorm;
		}
		//cout<<")"<<endl;
		//cout<<endl;
		//cout<<endl;
		//cout<<endl;
	
		//cout<<"X2<-rbind(";
		for(int i=0;i<N2;i++){
			double S= aRegion->fPixelCollection[i].S;
			double SNorm= NormMin + (NormMax-NormMin)*(S-fImageMinS)/(fImageMaxS-fImageMinS);
			//cout<<S<<",";
			X2(i,0)= SNorm;
		}	
		//cout<<")"<<endl;
		//cout<<endl;
		//cout<<endl;
		//cout<<endl;

		(*fR)["X1"]= X1;
		(*fR)["X2"]= X2;	
		(*fR)["k"]= nNearestNeighbors;

		//### Compute Entropy params
		RCmd= std::string("H1List<-entropy(X1, k, \"kd_tree\"); H2List<-entropy(X2, k, \"kd_tree\"); CE12List<-crossentropy(X1,X2,k, \"kd_tree\");CE21List<-crossentropy(X2,X1,k, \"kd_tree\");");
		fR->parseEval(RCmd);

		RCmd= std::string("H1<-H1List[k]; H2<-H2List[k]; CE_12<-CE12List[k]; CE_21<-CE21List[k]");
		fR->parseEval(RCmd);

		//## Importing results
		double H1= Rcpp::as< double >(fR->parseEval(std::string("H1")));
		double H2= Rcpp::as< double >(fR->parseEval(std::string("H2")));
		double CE_12= Rcpp::as< double >(fR->parseEval(std::string("CE_12")));
		double CE_21= Rcpp::as< double >(fR->parseEval(std::string("CE_21")));
		cout<<"H1="<<H1<<" H2="<<H2<<" CE_12="<<CE_12<<" CE_21="<<CE_21<<endl;
		entropyPars.push_back(H1);
		entropyPars.push_back(H2);
		entropyPars.push_back(CE_12);
		entropyPars.push_back(CE_21);
	}//close try block
	catch( std::exception &ex ) {
		cerr << "Region::GetEntropyDistance(): ERROR: Exception catched: " << ex.what() << endl;
		return entropyPars;
  } 
	catch(...) { 
		cerr << "Region::GetEntropyDistance(): ERROR: C++ exception (unknown reason)" << endl;
		return entropyPars;
  }	

	return entropyPars;

}//close Region::GetEntropyDistance()


int Region::ComputeKLD(Region* aRegion,double& KLD){

	if(!aRegion) return -1;
	
	//## Compute Kullback-Leibler divergence & Entropy of region
	int nNearestNeighbors= 5;//used for density estimation
	KLD= TMath::Infinity();//init to inf

	try{
		RInside* fR= RInside::instancePtr();
		if(!fR){
			cerr<<"Region::ComputeKLD(): ERROR: Cannot retrieve RInside instance (did you create one in your main application?)!"<<endl;
			fR= new RInside;
		} 

		//## Load R library for outlier detection
		//cout<<"Region::ComputeKLD(): INFO: Loading needed R packages..."<<endl;
		std::string RCmd= std::string("library(\"FNN\");");
		fR->parseEvalQ(RCmd);

		//## Import data matrix in R environment
		//cout<<"Region::ComputeKLD(): INFO: Importing data matrix in R environment..."<<endl;
		int N1= (int)fPixelCollection.size();
		Rcpp::NumericMatrix X1(N1,1);	
		int N2= (int)aRegion->fPixelCollection.size();
		Rcpp::NumericMatrix X2(N2,1);
		nNearestNeighbors= sqrt(std::min(N1,N2));	
		
		double NormMin= 1;
		double NormMax= 256;

		for(int i=0;i<N1;i++){
			double S= fPixelCollection[i].S;
			double SNorm= NormMin + (NormMax-NormMin)*(S-fImageMinS)/(fImageMaxS-fImageMinS);
			X1(i,0)= SNorm;
		}
		for(int i=0;i<N2;i++){
			double S= aRegion->fPixelCollection[i].S;
			double SNorm= NormMin + (NormMax-NormMin)*(S-fImageMinS)/(fImageMaxS-fImageMinS);
			X2(i,0)= SNorm;
		}

		(*fR)["X1"]= X1;
		(*fR)["X2"]= X2;	
		(*fR)["nNearestNeighbors"]= nNearestNeighbors;

		//### Compute Kullback-Leibler divergence & distance
		//cout<<"Region::ComputeKLD(): INFO: Running KLD algo (N1="<<N1<<", N2="<<N2<<"k="<<nNearestNeighbors<<")..."<<endl;
		//RCmd= std::string("kld<-KL.divergence(X1, X2, nNearestNeighbors, \"kd_tree\");");
		RCmd= std::string("kld<-KL.dist(X1, X2, nNearestNeighbors, \"kd_tree\"); E<-entropy(X1, nNearestNeighbors, \"kd_tree\"); CE<-crossentropy(X1,X2,nNearestNeighbors, \"kd_tree\");");
		fR->parseEval(RCmd);
		
		cout<<"== X1 =="<<endl;
		fR->parseEval("print(X1)");
		cout<<"== X2 =="<<endl;
		fR->parseEval("print(X2)");
		cout<<"== KLD =="<<endl;
		fR->parseEval("print(kld)");
		cout<<"== E =="<<endl;
		fR->parseEval("print(E)");
		cout<<"== CE =="<<endl;
		fR->parseEval("print(CE)");

		//## Importing results
		Rcpp::NumericVector kld = fR->parseEval(std::string("kld"));
		Rcpp::NumericVector Entr = fR->parseEval(std::string("E"));
		Rcpp::NumericVector CrossEntr = fR->parseEval(std::string("CE"));
		KLD= kld(nNearestNeighbors-1);
		double Entropy= Entr(nNearestNeighbors-1);
		double CrossEntropy=  CrossEntr(nNearestNeighbors-1);
		cout<<"KLD="<<KLD<<" CE-E="<<CrossEntropy-Entropy<<endl;

	}//close try block
	catch( std::exception &ex ) {
		cerr << "Region::ComputeKLD(): ERROR: Exception catched: " << ex.what() << endl;
		return -1;
  } 
	catch(...) { 
		cerr << "Region::ComputeKLD(): ERROR: C++ exception (unknown reason)" << endl;
		return -1;
  }	

	return 0;
	
}//close Region::ComputeKLD()




Img* Region::GetImage(RegionImageMode mode){

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
	TString imgName= Form("RegionImg_%d",fId);
	Img* regionImg= new Img(imgName,imgName,nBoxX,boundingBoxX[0]-0.5,boundingBoxX[1]+0.5,nBoxY,boundingBoxY[0]-0.5,boundingBoxY[1]+0.5);
	Img* regionSignificanceMap= 0;

	if(mode==eBinarizedMap){
		for(unsigned int k=0;k<fPixelCollection.size();k++){
			Pixel thisPixel= fPixelCollection[k];
			double thisX= thisPixel.x;
			double thisY= thisPixel.y;
			regionImg->FillPixel(thisX,thisY,1);
		}//end loop pixels
	}//close if
	else if(mode==eFluxMap){
		for(unsigned int k=0;k<fPixelCollection.size();k++){
			Pixel thisPixel= fPixelCollection[k];
			double thisX= thisPixel.x;
			double thisY= thisPixel.y;
			double thisS= thisPixel.S;	
			regionImg->FillPixel(thisX,thisY,thisS);
		}//end loop pixels
	}//close else
	else if(mode==eMeanMap){
		for(unsigned int k=0;k<fPixelCollection.size();k++){
			Pixel thisPixel= fPixelCollection[k];
			double thisX= thisPixel.x;
			double thisY= thisPixel.y;
			double thisS= thisPixel.S;	
			regionImg->FillPixel(thisX,thisY,fMean);
		}//end loop pixels
	}//close else
	else if(mode==eRMSMap){
		for(unsigned int k=0;k<fPixelCollection.size();k++){
			Pixel thisPixel= fPixelCollection[k];
			double thisX= thisPixel.x;
			double thisY= thisPixel.y;
			double thisS= thisPixel.S;	
			regionImg->FillPixel(thisX,thisY,fRMS);
		}//end loop pixels
	}//close else

	return regionImg;

}//close Region::GetImage()


void Region::Draw(RegionImageMode mode){

	cout<<"Region::Draw(): INFO: Filling region image..."<<fId<<endl;
	Img* regionImg= GetImage(mode);

	//## Draw stuff for check
	cout<<"Region::Draw(): INFO: Drawing region id "<<fId<<endl;
	TString canvasName= Form("RegionPlot_%d",fId);
	TCanvas* RegionPlot= new TCanvas(canvasName,canvasName);
	RegionPlot->cd();
	
	if(regionImg) regionImg->Draw("COLZ");
	
	cout<<"Region::Draw(): INFO: Drawing region "<<fId<<" contours..."<<endl;
	
	for(unsigned int i=0;i<fContourCollection.size();i++){
		//Get graph 
		cout<<"Region::Draw(): INFO: Drawing region "<<fId<<" contour no. "<<i<<"..."<<endl;
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
	}//end loop contours


}//close Draw()


