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
* @file MathUtils.h
* @class MathUtils
* @brief Utility functions for math tasks
*
* Utility functions for math tasks
* @author S. Riggi
* @date 15/01/2016
*/


#ifndef MathUtils_h
#define MathUtils_h 1

//OpenCV
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>


#include <TObject.h>
#include <TMath.h>
#include <TMatrixD.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

#include <vector>
#include <algorithm>
#include <map>
#include <time.h>
#include <ctime>

#include <complex>

using namespace std;

namespace Caesar {


class MathUtils : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    MathUtils();
		/**
		* \brief Class destructor: free allocated memory
		*/
   	virtual ~MathUtils();

		
	public:

		/**
		* \brief Return 2D grid partition given Nx x Ny pixels and box sizes
		*/
		static int Compute2DGrid(
			std::vector<long int>& ix_min,std::vector<long int>& ix_max,
			std::vector<long int>& iy_min,std::vector<long int>& iy_max,
			long int Nx,long int Ny,long int boxSizeX,long int boxSizeY,float gridStepSizeX,float gridStepSizeY
		);

		/**
		* \brief Perform bilinear interpolation on regular grid
		*/
		static int BiLinearInterpolation(
			std::vector<double>const& sampled_gridX,std::vector<double>const& sampled_gridY,
			std::vector<double>const& sampledZ,
			std::vector<double>const& interp_gridX,std::vector<double>const& interp_gridY,
			std::vector<double>& interpZ
		);

		static cv::Mat GetConvolution(cv::Mat I, cv::Mat kernel){
			// anchor Point
    	cv::Point anchor(kernel.cols - kernel.cols/2 - 1, kernel.rows - kernel.rows/2 - 1);

    	// flip the kernel
    	cv::Mat kernel_flipped;
    	int flip_code = -1; // -1-both axis, 0-x-axis, 1-y-axis
    	cv::flip(kernel, kernel_flipped, flip_code);

			//Apply convolution
			cv::Mat dst;
			filter2D(I, dst, I.depth(), kernel_flipped, anchor, cv::BORDER_CONSTANT);

			return dst;
		}//close GetConvolution()


		static cv::Mat GetConvolution2(cv::Mat I, cv::Mat kernel){
			int nrows= I.rows;
			int ncols= I.cols;
			int nrows_kernel= kernel.rows;
			int ncols_kernel= kernel.cols;
			cv::Mat dst= cv::Mat::zeros(nrows,ncols,CV_64FC1);

			for(int i=0;i<nrows;i++){//loop rows
				for(int j=0;j<ncols;j++){//loop cols

					for(int k=0;k<nrows_kernel;k++){//start loops on filter box rows
						int index_x= i-k+nrows_kernel/2;
						int mirror_index_x= GetMirrorIndex(index_x,nrows);

						for(int l=0;l<ncols_kernel;l++){//start loops on filter box cols
							int index_y= j-l+ncols_kernel/2;
							int mirror_index_y= GetMirrorIndex(index_y,ncols);	

							double K= kernel.at<double>(k,l);
							double f= I.at<double>(mirror_index_x,mirror_index_y);	
							dst.at<double>(i,j)= K*f;
						}//end loop kernel cols
					}//end loop kernel rows
				}//end loop y
			}//end loop x

			return dst;
		
		}//close GetConvolution()


		static cv::Mat GetATrousConvolution(cv::Mat I, cv::Mat kernel,int scale){
			int nrows= I.rows;
			int ncols= I.cols;
			int nrows_kernel= kernel.rows;
			int ncols_kernel= kernel.cols;
			int rowIndex= (nrows_kernel-1)/2;
			int colIndex= (ncols_kernel-1)/2;
			int filterGap= pow(2,scale-1); 
			cv::Mat dst= cv::Mat::zeros(nrows,ncols,CV_64FC1);
			
			//Compute convolution
			for(int i=0;i<nrows;i++){
				for(int j=0;j<ncols;j++){

					for(int k=-rowIndex;k<rowIndex;k++){//start loops on filter box rows
						int index_x= i + filterGap*k;
						int mirror_index_x= GetMirrorIndex(index_x,nrows);
		
						for(int l=-colIndex;l<colIndex;l++){//start loops on filter box cols
							int index_y= j + filterGap*l;
							int mirror_index_y= GetMirrorIndex(index_y,ncols);
		
							double K= kernel.at<double>(k+rowIndex,l+colIndex);
							double f= I.at<double>(mirror_index_x,mirror_index_y); 
							dst.at<double>(i,j)+= K*f;
						}//end loop filter cols
					}//end loop filter rows
				}//end loop cols
			}//end loop rows

			return dst;
		}//close GetATrousConvolution()


		static int GetMirrorIndex(int index,int N){
			int mirror_index= 0;
			if(index>=0 && index<=(N-1)) {
				mirror_index= index;
			}
			else if(index<0){
				mirror_index= -index;
			}
			else if(index>(N-1)){
				mirror_index= 2*(N-1)-index;
			}
			else{
				cerr<<"GetMirrorIndex(): ERROR: Invalid index of matrix size passed...exit!"<<endl;
				return -1;
			}	
			return mirror_index;
		}//close GetMirrorIndex()

		static double GetMatrixTrace(TMatrixD* T){			
			double trace= 0;
			for(int i=0;i<T->GetNrows();i++){
				trace+= (*T)(i,i);		
			}
  		return trace;
		}//close GetMatrixTrace()

		static std::vector< std::complex<double> > DFTShifted(std::vector< std::complex<double> > data, int n){
			int N= (int)data.size();
			if(n>N || n<0) n= N;//check size

			std::vector< std::complex<double> > Fn(n,0);
			for(int i=0; i<n; i++) {//loop over n
				Fn[i] = std::complex<double>(0.,0.);
				int s= -floor(N/2.) + i;
				//int s= i;
	
				for(int j=0;j<N;j++) {//loop over data size
					int k= j;
					//int k= -N/2 + j;
					double arg= 2.*TMath::Pi()*s*k/N;
					std::complex<double> prod= std::polar(1.,-arg);
					Fn[i]+= data[j]*prod;
				}//end loop data size
				Fn[i]/= (double)N;
			}//end loop n

			return Fn;
		}//close DFT()

		static std::vector< std::complex<double> > DFT(std::vector< std::complex<double> > data, int n){
			int N= (int)data.size();
			if(n>N || n<0) n= N;//check size

			std::vector< std::complex<double> > Fn(n,0);
			for(int i=0; i<n; i++) {//loop over n
				Fn[i] = std::complex<double>(0.,0.);
				//int s= -floor(N/2.) + i;
				int s= i;

				for(int j=0;j<N;j++) {//loop over data size
					int k= j;
					//int k= -floor(N/2.) + j;
					//double arg= 2.*TMath::Pi()*i*k/N;
					double arg= 2.*TMath::Pi()*s*k/N;
					std::complex<double> prod= std::polar(1.,-arg);
					Fn[i]+= data[j]*prod;
				}//end loop data size
			}//end loop n

			return Fn;
		}//close DFT()

		static std::vector< std::complex<double> > IDFT(std::vector< std::complex<double> > data, int n){
			int N= (int)data.size();
			if(n>N || n<0) n= N;//check size

			std::vector< std::complex<double> > fn(n,0);
			for(int i=0; i<n; i++) {//loop over n
				fn[i] = std::complex<double>(0.,0.);
				int s= i;

				for(int j=0;j<N;j++) {//loop over data size
					int k= j;
					//int k= -floor(N/2.) + j;
					
					double arg= 2.*TMath::Pi()*s*k/N;
					std::complex<double> prod= std::polar(1.,arg);
					fn[i]+= data[j]*prod;
				}//end loop data size
				fn[i]/= (double)N;	
			}//end loop n

			return fn;
		}//close DFT()

		static int EtaAuxiliaryFcn(int s,int N){
			int thr= -floor(N/2.) + N-1;
			int fval= 0;
			if(s<=thr) fval= s;
			else fval= N-s;
			return fval;
		}

		static std::vector<double> GetContourCurvature(std::vector< std::complex<double> > data,double sigma){
			int N= (int)data.size();
			std::vector< std::complex<double> > Us= DFT(data,N);

			//Compute Us smoothed with a gaussian
			std::vector< std::complex<double> > Us_firstDeriv; 
			std::vector< std::complex<double> > Us_secondDeriv; 
			std::vector< std::complex<double> > Us_firstDeriv_smoothed;
			std::vector< std::complex<double> > Us_secondDeriv_smoothed;
			for(int i=0;i<N;i++){
				//int s= -floor(N/2.) + i;
				//double Gs= sigma/sqrt(2*TMath::Pi())*exp(-pow(0.5*sigma*s,2));
				//double arg= 2*TMath::Pi()*s;
				int eta= EtaAuxiliaryFcn(i,N);
				double Gs= exp(-pow(sigma*eta,2));
				double arg= 2*TMath::Pi()*eta;
				std::complex<double> z(0,arg); 
				std::complex<double> thisUs_firstDeriv= z*Us[i];	
				std::complex<double> thisUs_secondDeriv= -pow(arg,2)*Us[i];
				std::complex<double> thisUs_firstDeriv_smoothed= thisUs_firstDeriv*Gs;
				std::complex<double> thisUs_secondDeriv_smoothed= thisUs_secondDeriv*Gs;
				Us_firstDeriv.push_back(thisUs_firstDeriv);
				Us_secondDeriv.push_back(thisUs_secondDeriv);
				Us_firstDeriv_smoothed.push_back(thisUs_firstDeriv_smoothed);
				Us_secondDeriv_smoothed.push_back(thisUs_secondDeriv_smoothed);
				
			}//end loop points

			//Compute ut'= IDST(Us')=IDST(i x 2pi x s x Us)	
			std::vector< std::complex<double> > ut_firstDeriv= IDFT(Us_firstDeriv,N);
			double L= 0;
			for(unsigned int i=0;i<ut_firstDeriv.size();i++) L+= std::abs(ut_firstDeriv[i]);	
			L*= 2.*TMath::Pi()/N;
		
			//Compute inverse transform of smoothed
			std::vector< std::complex<double> > ut_firstDeriv_smoothed= IDFT(Us_firstDeriv_smoothed, N);
			std::vector< std::complex<double> > ut_secondDeriv_smoothed= IDFT(Us_secondDeriv_smoothed, N);

			//Compute smoothed perymeter
			double L_smoothed= 0;
			for(int i=0;i<N;i++){
				std::complex<double> u1= ut_firstDeriv_smoothed[i];
				double u1_mod= std::abs(u1);  	
				L_smoothed+= u1_mod;
			}//end loop points
			L_smoothed*= 2.*TMath::Pi()/N;
		
			//Compute curvature
			double normFactor= L/L_smoothed;
			std::vector<double> curvature;
			for(int i=0;i<N;i++){
				//Normalize ut
				ut_firstDeriv_smoothed[i]*= normFactor;
				ut_secondDeriv_smoothed[i]*= normFactor;

				std::complex<double> u2_conj= std::conj(ut_secondDeriv_smoothed[i]); 
				std::complex<double> u1= ut_firstDeriv_smoothed[i];
				double u1_mod= std::abs(u1);  	
				double curv= -std::imag(u1*u2_conj)/pow(u1_mod,3);
				//curv*= L;//multiply by perymeter length
				//cout<<"Pnt no. "<<i<<" usmooth'="<<ut_firstDeriv_smoothed[i]<<" usmooth''="<<ut_secondDeriv_smoothed[i]<<" u1_mod="<<u1_mod<<" curv="<<curv<<endl;
				curvature.push_back(curv);
			}//end loop points

			return curvature;
		}//close GetCurveCurvature()

	private:
	
		ClassDef(MathUtils,1)
};

#ifdef __MAKECINT__
#pragma link C++ class MathUtils+;
#endif	

}//close namespace


#endif 
