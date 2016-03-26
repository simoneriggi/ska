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
* @file Region.cc
* @class Region
* @brief Region data class
*
* Superpixel Region data
* @author S. Riggi
* @date 20/01/2015
*/

#include <Region.h>
#include <Pixel.h>
#include <Blob.h>
#include <Img.h>
#include <StatsUtils.h>

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


ClassImp(Caesar::Region)


namespace Caesar {

Region::Region() {
	Tag= eUntagged;
	m_SubRegionIds.clear();
}

Region::~Region(){
	
}//close destructor


Region::RegionPars* Region::GetParams(bool includeCurvPar){

	if(!m_HasStats){
		cerr<<"Region::GetParams(): ERROR: No stats computed for this region!"<<endl;
		return 0;
	}

	int nPars= 2;
	int nPars_robust= 2;
	int nPars_spatial= 2;
	if(includeCurvPar){
		nPars+= 2;
		nPars_robust+= 2;
	}
	
	Region::RegionPars* params= new Region::RegionPars;
	params->pars= new TVectorD(nPars);
	params->robustPars= new TVectorD(nPars_robust);
	params->spatialPars= new TVectorD(nPars_spatial);

	//== Appearance pars ==
	(*(params->pars))(0)= Mean;
	(*(params->pars))(1)= RMS;
	
	//== Robust appearance pars ==
	(*(params->robustPars))(0)= Median;
	(*(params->robustPars))(1)= MedianRMS;
	
	if(includeCurvPar){
		(*(params->pars))(2)= Mean_curv;
		(*(params->pars))(3)= RMS_curv;

		(*(params->robustPars))(2)= Median_curv;
		(*(params->robustPars))(3)= MedianRMS_curv;
	}

	//== Spatial pars ==
	(*(params->spatialPars))(0)= X0;
	(*(params->spatialPars))(1)= Y0;

	return params;

}//close GetParams()


int Region::AddRegion(Region* aRegion,bool addPixels,bool copyPixels){

	if(!aRegion) {
		cerr<<"Region::AddRegion(): ERROR: Null prt to given region, nothing will be added!"<<endl;
		return -1;
	}
	
	//Update moments for the region
	double S_B= aRegion->m_S;
	double Scurv_B= aRegion->m_S_curv;
	double Sedge_B= aRegion->m_S_edge;
	double Smin_B= aRegion->m_Smin;
	double Smax_B= aRegion->m_Smax;
	int PixIdmin_B= aRegion->m_PixIdmin;
	int PixIdmax_B= aRegion->m_PixIdmax;
	double Sxx_B= aRegion->m_Sxx;
	double Syy_B= aRegion->m_Syy;
	double Sxy_B= aRegion->m_Sxy;
	double Sx_B= aRegion->m_Sx;
	double Sy_B= aRegion->m_Sy;
	double X0_B= aRegion->X0;
	double Y0_B= aRegion->Y0;
	int Xmin_B= aRegion->m_Xmin;
	int Xmax_B= aRegion->m_Xmax;
	int Ymin_B= aRegion->m_Ymin;
	int Ymax_B= aRegion->m_Ymax;	
	int Ixmin_B= aRegion->m_Ix_min;
	int Ixmax_B= aRegion->m_Ix_max;
	int Iymin_B= aRegion->m_Iy_min;
	int Iymax_B= aRegion->m_Iy_max;
	double N_B= aRegion->NPix;
	double M1_B= aRegion->m_M1;
	double M2_B= aRegion->m_M2;
	double M3_B= aRegion->m_M3;
	double M4_B= aRegion->m_M4;
	double M1Curv_B= aRegion->m_M1_curv;
	double M2Curv_B= aRegion->m_M2_curv;
	
	
	double S_A= this->m_S;
	double Scurv_A= this->m_S_curv;
	double Sedge_A= this->m_S_edge;
	double Smin_A= this->m_Smin;
	double Smax_A= this->m_Smax;
	int PixIdmin_A= this->m_PixIdmin;
	int PixIdmax_A= this->m_PixIdmax;
	double Sxx_A= this->m_Sxx;
	double Syy_A= this->m_Syy;
	double Sxy_A= this->m_Sxy;
	double Sx_A= this->m_Sx;
	double Sy_A= this->m_Sy;
	double X0_A= this->X0;
	double Y0_A= this->Y0;
	double Xmin_A= this->m_Xmin;
	double Xmax_A= this->m_Xmax;
	double Ymin_A= this->m_Ymin;
	double Ymax_A= this->m_Ymax;
	int Ixmin_A= this->m_Ix_min;
	int Ixmax_A= this->m_Ix_max;
	int Iymin_A= this->m_Iy_min;
	int Iymax_A= this->m_Iy_max;
	double N_A= this->NPix;
	double M1_A= this->m_M1;
	double M2_A= this->m_M2;
	double M3_A= this->m_M3;
	double M4_A= this->m_M4;
	double M1Curv_A= this->m_M1_curv;
	double M2Curv_A= this->m_M2_curv;
	
	
	double delta= (N_A*M1_B-N_B*M1_A)/(N_A*N_B);

	NPix= N_A + N_B;
	m_M1= (N_A*M1_A + N_B*M1_B)/NPix;
	m_M2= M2_A + M2_B + pow(N_B*M1_A-N_A*M1_B,2)/(N_A*N_B*NPix);
	m_M3= M3_A + M3_B + pow(delta,3)*N_A*N_B*(N_A-N_B)/(NPix*NPix) + 3*delta*(N_A*M2_B-N_B*M2_A)/NPix;
	m_M4= M4_A + M4_B + pow(delta,4)*N_A*N_B*(N_A*N_A-N_A*N_B+N_B*N_B)/pow(NPix,3) + 6*pow(delta,2)*(N_A*N_A*M2_B+N_B*N_B*M2_A)/pow(NPix,2) + 4*delta*(N_A*M3_B-N_B*M3_A)/NPix;

	m_M1_curv= (N_A*M1Curv_A + N_B*M1Curv_B)/NPix;
	m_M2_curv= M2Curv_A + M2Curv_B + pow(N_B*M1Curv_A-N_A*M1Curv_B,2)/(N_A*N_B*NPix);

	if(Smin_A<Smin_A) {
		m_Smin= Smin_A;
		m_PixIdmin= PixIdmin_A;
	}
	else{
		m_Smin= Smin_B;
		m_PixIdmin= PixIdmin_B;
	}

	if(Smax_A>Smax_A) {
		m_Smax= Smax_A;
		m_PixIdmax= PixIdmax_A;
	}
	else{
		m_Smax= Smax_B;
		m_PixIdmax= PixIdmax_B;
	}

	m_S= S_A + S_B;
	m_S_curv= Scurv_A + Scurv_B;
	m_S_edge= Sedge_A + Sedge_B;
	m_Sxx= Sxx_A + Sxx_B;
	m_Syy= Syy_A + Syy_B;
	m_Sxy= Sxy_A + Sxy_B;
	m_Sx= Sx_A + Sx_B;
	m_Sy= Sy_A + Sy_B;
	X0= X0_A + X0_B;
	Y0= Y0_A + Y0_B;

	if(Xmin_A<Xmin_B) m_Xmin= Xmin_A;
	else m_Xmin= Xmin_B;

	if(Xmax_A>Xmax_B) m_Xmax= Xmax_A;
	else m_Xmax= Xmax_B;

	if(Ymin_A<Ymin_B) m_Ymin= Ymin_A;
	else m_Ymin= Ymin_B;

	if(Ymax_A>Ymax_B) m_Ymax= Ymax_A;
	else m_Ymax= Ymax_B;

	if(Ixmin_A<Ixmin_B) m_Ix_min= Ixmin_A;
	else m_Ix_min= Ixmin_B;

	if(Ixmax_A>Ixmax_B) m_Ix_max= Ixmax_A;
	else m_Ix_max= Ixmax_B;

	if(Iymin_A<Iymin_B) m_Iy_min= Iymin_A;
	else m_Iy_min= Iymin_B;

	if(Iymax_A>Iymax_B) m_Iy_max= Iymax_A;
	else m_Iy_max= Iymax_B;

	//Add new pixels at list end
	if(addPixels){
		if(copyPixels){//copying pixel to new pointers
			Pixel* aNewPixel= 0;
			for(int i=0;i<aRegion->GetNPixels();i++){
				Pixel* thisPixel= aRegion->GetPixel(i);
				aNewPixel= new Pixel;
				*aNewPixel= *thisPixel;
				m_Pixels.push_back(aNewPixel);
			}//end loop pixels
		}//clse if
		else{
			std::vector<Pixel*>::iterator it= m_Pixels.end();
			m_Pixels.insert(it,(aRegion->m_Pixels).begin(),(aRegion->m_Pixels).end());
		}
	}

	//Append subregion ids 
	int Id_B= aRegion->Id;
	//std::vector<int> SubRegionIds_B= aRegion->m_SubRegionIds;
	m_SubRegionIds.push_back(Id_B);

	std::vector<int>::iterator it2= m_SubRegionIds.end();
	m_SubRegionIds.insert(it2,(aRegion->m_SubRegionIds).begin(),(aRegion->m_SubRegionIds).end());

	return 0;

}//close Region::AddRegion()

int Region::GetDistance(double& dist_color,double& dist_space,Region* aRegion,bool useRobustParams,bool normalizeParams,bool addCurvDist){

	dist_color= 1.e+99;
	dist_space= 1.e+99;
	if(!aRegion){
		cerr<<"Region::GetDistance(): ERROR: Null ptr to given region...returning inf dists!"<<endl;
		return -1;
	}

	//## Normalization of parameters for euclidean distance computation
	//## Under linear transformation A+Bx --> mean'=A+B*mean, rms'=B*rms
	double NormMin= 0;//1;
	double NormMax= 1;//256;
	
	double A= NormMin - (NormMax-NormMin)*m_ImageMinS/(m_ImageMaxS-m_ImageMinS);
	double B= (NormMax-NormMin)/(m_ImageMaxS-m_ImageMinS);
	double Acurv= NormMin - (NormMax-NormMin)*m_ImageMinScurv/(m_ImageMaxScurv-m_ImageMinScurv);
	double Bcurv= (NormMax-NormMin)/(m_ImageMaxScurv-m_ImageMinScurv);

	//double Ax0= NormMin - (NormMax-NormMin)*m_ImageMinX/(m_ImageMaxX-m_ImageMinX);
	//double Bx0= (NormMax-NormMin)/(m_ImageMaxX-m_ImageMinX);
	//double Ay0= NormMin - (NormMax-NormMin)*m_ImageMinY/(m_ImageMaxY-m_ImageMinY);
	//double By0= (NormMax-NormMin)/(m_ImageMaxY-m_ImageMinY);

	//Mean difference
	double mean= this->Mean;
	double mean_norm= A+B*mean;//mean under linear transformation
	double meanN= aRegion->Mean;
	double meanN_norm= A+B*meanN;
	double meanDiff= mean-meanN;
	double meanDiff_norm= mean_norm-meanN_norm;
	
	//RMS difference
	double rms= this->RMS;
	double rms_norm= B*rms;
	double rmsN= aRegion->RMS;
	double rmsN_norm= B*rmsN;
	double rmsDiff= rms-rmsN;
	double rmsDiff_norm= rms_norm-rmsN_norm;
	
	//Median difference
	double median= this->Median;
	double median_norm= A+B*median;
	double medianN= aRegion->Median;
	double medianN_norm= A+B*medianN;
	double medianDiff= median-medianN;
	double medianDiff_norm= median_norm-medianN_norm;

	//Median MAD difference
	double mad= this->MedianRMS;
	double mad_norm= B*mad; 
	double madN= aRegion->MedianRMS;
	double madN_norm= B*madN;
	double madDiff= mad-madN;
	double madDiff_norm= mad_norm-madN_norm;
	
	//Curvature mean
	double curvMean= this->Mean_curv;
	double curvMean_norm= Acurv+Bcurv*curvMean;
	double curvMeanN= aRegion->Mean_curv;
	double curvMeanN_norm= Acurv+Bcurv*curvMeanN;
	double curvMeanDiff= curvMean-curvMeanN;
	double curvMeanDiff_norm= curvMean_norm-curvMeanN_norm;

	//Curvature rms
	double curvRMS= this->RMS_curv;
	double curvRMS_norm= Bcurv*curvRMS;
	double curvRMSN= aRegion->RMS_curv;
	double curvRMSN_norm= Bcurv*curvRMSN;
	double curvRMSDiff= curvRMS-curvRMSN;
	double curvRMSDiff_norm= curvRMS_norm-curvRMSN_norm;

	//Centroid distance
	double X0= this->X0;
 	double Y0= this->Y0;
 	//double X0_norm= Ax0+Bx0*X0; 
	//double Y0_norm= Ay0+By0*Y0; 
	double X0N= aRegion->X0;
 	double Y0N= aRegion->Y0;
	//double X0N_norm= Ax0+Bx0*X0N;
	//double Y0N_norm= Ay0+By0*Y0N;
	//double centroidDiffX= X0-X0N;
	//double centroidDiffY= Y0-Y0N;
	//double centroidDiffX_norm= X0_norm-X0N_norm;
	//double centroidDiffY_norm= Y0_norm-Y0N_norm;

	double centroidDistMin= 0;
	double centroidDistMax= sqrt(pow(m_ImageMaxX-m_ImageMinX,2) + pow(m_ImageMaxY-m_ImageMinY,2));
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
	
	return 0;
	
}//close Region::GetDistance()


int Region::GetAsymmDistance(double& dist,double& dist_neighbor,Region* aRegion,bool useRobustParams,bool normalizeParams,bool addSpatialDist,bool addCurvDist){

	dist= 1.e+99;
	dist_neighbor= 1.e+99;
	if(!aRegion){
		cerr<<"Region::GetAsymmDistance(): ERROR: Null ptr to given region...returning inf dists!"<<endl;
		return -1;
	}
	
	if(useRobustParams){
		//Create a merged region between current and neighbor	
		Region mergedRegion= *(this);
		mergedRegion.AddRegion(aRegion,useRobustParams);
		mergedRegion.ComputeStats(true,false);
		
		//Find distance between this region and merged region
		double dist_color, dist_space;
		int status= this->GetDistance(dist_color,dist_space,&mergedRegion,useRobustParams,normalizeParams,addCurvDist);
		if(status<0){
			cerr<<"Region::GetAsymmDistance(): ERROR: Failed to compute symm region distance!"<<endl;
			return -1;
		}
		dist= dist_color;

		//Find distance between neighbor region and merged region
		status= aRegion->GetDistance(dist_color, dist_space,&mergedRegion,useRobustParams,normalizeParams,addCurvDist);
		if(status<0){
			cerr<<"Region::GetAsymmDistance(): ERROR: Failed to compute symm neighbor region distance!"<<endl;
			return -1;
		}
		dist_neighbor= dist_color;

	}//close if
	else{//fast computing for mean & rms params
		double N_A= (double)NPix;
		double N_B= (double)aRegion->NPix; 
		double N= (double)N_A + (double)N_B;
		double M1_A= m_M1;
		double M1_B= aRegion->m_M1;
		double M2_A= m_M2;
		double M2_B= aRegion->m_M2;
		double M1Curv_A= m_M1_curv;
		double M1Curv_B= aRegion->m_M1_curv;
		double M2Curv_A= m_M2_curv;
		double M2Curv_B= aRegion->m_M2_curv;
		double Mean_A= Mean;
		double Mean_B= aRegion->Mean;
		double RMS_A= RMS;
		double RMS_B= aRegion->RMS;
		double MeanCurv_A= Mean_curv;
		double MeanCurv_B= aRegion->Mean_curv;
		double RMSCurv_A= RMS_curv;
		double RMSCurv_B= aRegion->RMS_curv;

		double S_A= m_S;
		double Cx_A= m_Sx;
		double Cy_A= m_Sy;
		double S_B= aRegion->m_S;
		double Cx_B= aRegion->m_Sx;
		double Cy_B= aRegion->m_Sy;
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
		double A= NormMin - (NormMax-NormMin)*m_ImageMinS/(m_ImageMaxS-m_ImageMinS);
		double B= (NormMax-NormMin)/(m_ImageMaxS-m_ImageMinS);
		double Acurv= NormMin - (NormMax-NormMin)*m_ImageMinScurv/(m_ImageMaxScurv-m_ImageMinScurv);	
		double Bcurv= (NormMax-NormMin)/(m_ImageMaxScurv-m_ImageMinScurv);
		double Ax0= NormMin - (NormMax-NormMin)*m_ImageMinX/(m_ImageMaxX-m_ImageMinX);
		double Bx0= (NormMax-NormMin)/(m_ImageMaxX-m_ImageMinX);
		double Ay0= NormMin - (NormMax-NormMin)*m_ImageMinY/(m_ImageMaxY-m_ImageMinY);
		double By0= (NormMax-NormMin)/(m_ImageMaxY-m_ImageMinY);
		double MeanNorm_A= A+B*Mean_A;
		double MeanNorm_B= A+B*Mean_B;	
		double MeanNorm= A+B*Mean;
		double RMSNorm_A= B*RMS_A;
		double RMSNorm_B= B*RMS_B;	
		double RMSNorm= B*RMS;
		double MeanCurvNorm_A= Acurv+Bcurv*MeanCurv_A;
		double MeanCurvNorm_B= Acurv+Bcurv*MeanCurv_B;	
		double MeanCurvNorm= Acurv+Bcurv*MeanCurv;
		double RMSCurvNorm_A= Bcurv*RMSCurv_A;
		double RMSCurvNorm_B= Bcurv*RMSCurv_B;	
		double RMSCurvNorm= Bcurv*RMSCurv;

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

		}//close if normalize pars
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
	}//close else
	
	return 0;

}//close GetAsymmDistance()

}//close namespace


