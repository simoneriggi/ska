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
* @file StatsUtils.h
* @class StatsUtils
* @brief Utility functions for statistical tasks
*
* Utility functions for statistical tasks
* @author S. Riggi
* @date 15/01/2016
*/


#ifndef StatsUtils_h
#define StatsUtils_h 1

#include <TObject.h>
#include <TMatrixD.h>
#include <TMath.h>

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
#include <string>
#include <time.h>
#include <ctime>

#include <complex>

using namespace std;

namespace Caesar {


class StatsUtils : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    StatsUtils();
		/**
		* \brief Class destructor: free allocated memory
		*/
   	virtual ~StatsUtils();

		
	public:

		//Median Fast (should run in O(n))
		template < typename T >
			static T GetMedianFast( std::vector<T>&vec){
				size_t n = vec.size()/2;
				if(n<=0) return -999;
				nth_element(vec.begin(),vec.begin()+n, vec.end());
  			return vec[n];
			}//close GetMedianFast()

		//Get median (should run in O(nlog n))
		template < typename T >
			static T GetMedian( std::vector<T>&vec, bool isSorted=false){//this should run in O(nlog(n))
				size_t n = vec.size();
				if(n<=0) return -999;
  			if(!isSorted) std::sort(vec.begin(), vec.end());
				double median= 0;		
  			if(n%2==0) median = (vec[n/2-1] + vec[n/2])/2;	
				else median = vec[n/2];
  			return median;
			}//close GetMedian()

		//MAD
		template < typename T >
			static T GetMAD( std::vector<T>const &vec, T median){
				size_t n = vec.size();
				if(n<=0) return -999;
				std::vector<double> MADs;
  			for(int j=0;j<n;j++){
					double diff= fabs(vec[j]-median);
					MADs.push_back(diff);
				}
				std::sort(MADs.begin(),MADs.end());
				double MAD= StatsUtils::GetMedian(MADs,true);
  			return MAD;
			}//close GetMAD()

		template < typename T >
			static T GetMADFast( std::vector<T>const &vec, T median){
				size_t n = vec.size();
				if(n<=0) return -999;
				std::vector<double> MADs;
  			for(int j=0;j<n;j++){
					double diff= fabs(vec[j]-median);
					MADs.push_back(diff);
				}
				double MAD= StatsUtils::GetMedianFast(MADs);
  			return MAD;
			}//close GetMAD()

		//Clipped estimators
		template < typename T >
			static std::pair<T,T> GetClippedEstimators( std::vector<T>const &vec, double clipsig=3, int maxiter=5, double tol=0.1, bool isSorted=false){
				int n= (int)vec.size();
				if(n<=0) {
					return std::make_pair(-999,-999);
				}
				std::vector<T> data(vec.begin(),vec.end());
				if(!isSorted) std::sort(data.begin(), data.end());

				int iter = 0; 
				double median= -999;
				double lastmedian= -999;
				double rms= -999;
				double lastrms= -999;	
				double diff= 0;	
				double eps= 1.e+99;
  			
				while ( eps>tol && (iter < maxiter) ){
  				n = (int)data.size();
					median= GetMedian(data,true);	
					//rms= 1.4826*Utils::GetMAD(data,median);
					double sum = std::accumulate(data.begin(), data.end(), 0.0);
					double mean = sum/(double)n;		
					double accum = 0.0;
					std::for_each (std::begin(data), std::end(data), [&](const double d) {
    				accum += (d - mean) * (d - mean);
					});
					double stdev = sqrt(accum / (double)(n-1));
					rms= stdev;

					std::vector<T> data_tmp;
					for(int j=0;j<n;j++){ 
						diff= fabs(data[j]-median);
						if(diff<clipsig*rms) data_tmp.push_back(data[j]);
					}

					//Check data sizes
					n= (int)data_tmp.size();
					if(n<=0) {	
						median= lastmedian;
						rms= lastrms;
						break;
					}
	
					data.clear();
					data.assign(data_tmp.begin(),data_tmp.end());
						
					eps= fabs(lastrms-rms)/rms;
					lastmedian= median;
					lastrms= rms;	
					iter++;
				}// End of while loop

				std::pair<T,T> res= std::make_pair (median,rms);	
				return res;
			}//close GetClippedEstimators()


		//BiWeight estimators
		template < typename T >
			static std::pair<T,T> GetBiWeightEstimators(std::vector<T>const &vec,double locationPar,double scalePar,double C=6,double tol=0.0001,int nmaxIter=10){
				
				double N= (double)(vec.size());	
				double Tb= locationPar;	
			 	double S= scalePar;
				int niter= 0;	
				while(niter<nmaxIter){
					niter++;

					//Update location Tb		
					double Tb_old= Tb;
					double sumw= 0;
					double sum= 0;
					for(unsigned j=0;j<vec.size();j++){
						double x= vec[j];
						double u= (x-Tb)/(C*S);
						double u2= u*u;
						if(u2>1) continue;
						double w= pow(1.-u2,2);
						sumw+= w;
						sum+= w*x;
					}//end loop pixels

					Tb= Tb_old;
					if(sumw>0) Tb= sum/sumw;

					//Check termination condition
					double eps= fabs(Tb/Tb_old-1);
					if(eps<tol) break;

					//Update scale
					double S_old= S;
					double sumNum= 0;
					double sumDenom= 0;
					for(unsigned j=0;j<vec.size();j++){
						double x= vec[j];
						//double u= (x-Tb)/(C*medianRMS);
						double u= (x-Tb)/(C*S);
		  			double u2= u*u;
						if(u2>1) continue;
						double num= ( pow(x-Tb,2) * pow(1.-u2,4) );
						double denom= (1.-u2)*(1.-5.*u2);
						sumNum+= num;
						sumDenom+= denom;
					}//end loop pixels		

					S= S_old;
					double S2= S_old*S_old;
					if(sumDenom>1) {
						S2= N*sumNum/(sumDenom*(sumDenom-1.)); 
						S= sqrt(S2);
					}
				}//end while

				std::pair<T,T> res= std::make_pair (Tb,S);	
				return res;
			}//close GetBiWeightEstimators()

  	//Mahalanobis distance
		static double GetMahalanobisDistance(TMatrixD x, TMatrixD mean, TMatrixD Sigma,bool isInverted=false);

		//PageRank
		static int ComputePageRank(std::vector<double>& ranks,TMatrixD& M,double d=0.85,double tol=1.e-4);

	private:
	
		ClassDef(StatsUtils,1)
};

#ifdef __MAKECINT__
#pragma link C++ class StatsUtils+;

#pragma link C++ function StatsUtils::GetMedianFast<float>;
#pragma link C++ function StatsUtils::GetMedianFast<double>;
#pragma link C++ function StatsUtils::GetMedianFast<int>;
#pragma link C++ function StatsUtils::GetMedianFast<long int>;

#pragma link C++ function StatsUtils::GetMedian<float>;
#pragma link C++ function StatsUtils::GetMedian<double>;
#pragma link C++ function StatsUtils::GetMedian<int>;
#pragma link C++ function StatsUtils::GetMedian<long int>;

#pragma link C++ function StatsUtils::GetMAD<float>;
#pragma link C++ function StatsUtils::GetMAD<double>;
#pragma link C++ function StatsUtils::GetMAD<int>;
#pragma link C++ function StatsUtils::GetMAD<long int>;

#pragma link C++ function StatsUtils::GetMADFast<float>;
#pragma link C++ function StatsUtils::GetMADFast<double>;
#pragma link C++ function StatsUtils::GetMADFast<int>;
#pragma link C++ function StatsUtils::GetMADFast<long int>;

#pragma link C++ function StatsUtils::GetClippedEstimators<float>;
#pragma link C++ function StatsUtils::GetClippedEstimators<double>;
#pragma link C++ function StatsUtils::GetClippedEstimators<int>;
#pragma link C++ function StatsUtils::GetClippedEstimators<long int>;

#pragma link C++ function StatsUtils::GetBiWeightEstimators<float>;
#pragma link C++ function StatsUtils::GetBiWeightEstimators<double>;
#pragma link C++ function StatsUtils::GetBiWeightEstimators<int>;
#pragma link C++ function StatsUtils::GetBiWeightEstimators<long int>;

#endif	

}//close namespace


#endif 
