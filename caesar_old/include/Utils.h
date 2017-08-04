/**
* @file Utils.h
* @class Utils
* @brief Utility functions
*
* Utility functions
* @author S. Riggi
* @date 23/08/2010
*/


#ifndef Utils_h
#define Utils_h 1

#include "Img.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TF1.h>
#include <TGraph.h>
#include <TVector3.h>
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


class Img;

class Utils {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    Utils();
		/**
		* \brief Class destructor: free allocated memory
		*/
   ~Utils();

		
	public:

		/**
		* \brief Delete selected items from a vector
		*/
		template<class T>
			static void DeleteItems(std::vector<T>& data, const std::vector<int>& deleteIndices) {
    	std::vector<bool> markedElements(data.size(), false);
    	std::vector<T> tempBuffer;
    	tempBuffer.reserve(data.size()-deleteIndices.size());

   	 	for (std::vector<int>::const_iterator itDel = deleteIndices.begin(); itDel != deleteIndices.end(); itDel++)
      	markedElements[*itDel] = true;

    	for (size_t i=0; i<data.size(); i++) {
      	if (!markedElements[i]) {
					tempBuffer.push_back(data[i]);
				}
    	}
    	data = tempBuffer;
		}//close DeleteItems()


		/**
		* \brief Order vectors and get ordering index
		*/
		template<class T> struct index_cmp{

  		index_cmp(const T arr) : arr(arr) {}
  		bool operator()(const size_t a, const size_t b) const
 			{
    		return arr[a] < arr[b];
  		}
  		const T arr;
		};
		template<class T> struct descending_index_cmp{

  		descending_index_cmp(const T arr) : arr(arr) {}
  		bool operator()(const size_t a, const size_t b) const
 			{
    		return arr[a] > arr[b];
  		}
  		const T arr;
		};

		template< class T >
			static void reorder(std::vector<T> & unordered,std::vector<size_t> const & index_map,std::vector<T> & ordered){
  			// copy for the reorder according to index_map, because unsorted may also be
  			// sorted
  			std::vector<T> copy = unordered;
  			ordered.resize(index_map.size());
  			for(unsigned int i=0;i<index_map.size();i++)
					ordered[i] = copy[index_map[i]];
			}

		template <class T>
			static void sort(std::vector<T> & unsorted,std::vector<T> & sorted,std::vector<size_t> & index_map){
  			// Original unsorted index map
  			index_map.resize(unsorted.size());
 				for(size_t i=0;i<unsorted.size();i++)
					index_map[i] = i;
  
  			// Sort the index map, using unsorted for comparison
  			std::sort(index_map.begin(),index_map.end(),index_cmp<std::vector<T>& >(unsorted));
  			sorted.resize(unsorted.size());
  			reorder(unsorted,index_map,sorted);
			}
	
		template <class T>
			static void sort_descending(std::vector<T> & unsorted,std::vector<T> & sorted,std::vector<size_t> & index_map){
  			// Original unsorted index map
  			index_map.resize(unsorted.size());
 				for(size_t i=0;i<unsorted.size();i++)
					index_map[i] = i;
  
  			// Sort the index map, using unsorted for comparison
  			std::sort(index_map.begin(),index_map.end(),descending_index_cmp<std::vector<T>& >(unsorted));
  			sorted.resize(unsorted.size());
  			reorder(unsorted,index_map,sorted);
			}


		template < typename T >
			static T GetMedian( std::vector<T>& vec, bool isSorted=false){
				size_t n = vec.size();
				if(n<=0) return -999;

  			if(!isSorted) std::sort(vec.begin(), vec.end());

				double median= 0;		
  			if(n%2==0) median = (vec[n/2-1] + vec[n/2])/2;	
				else median = vec[n/2];

  			return median;
			}

		template < typename T >
			static T GetMAD( std::vector<T>& vec, double median){
  			size_t n = vec.size();
				if(n<=0) return -999;

				std::vector<double> MADs;
  			for(unsigned int j=0;j<n;j++){
					double diff= fabs(vec[j]-median);
					MADs.push_back(diff);
				}
				std::sort(MADs.begin(),MADs.end());
				double MAD= Utils::GetMedian(MADs,true);
  			return MAD;
			}
		
	

		template < typename T >
			static std::pair<T,T> GetClippedEstimators( std::vector<T>& vec, double clipsig=3, int maxiter=5, double tol=0.1, bool isSorted=false){
				int n= (int)vec.size();
				if(n<=0) {
					//cerr<<"Utils::GetClippedEstimators(): ERROR: Empty vector given!"<<endl;
					return std::make_pair(-999,-999);
				}
				if(!isSorted) std::sort(vec.begin(), vec.end());
				std::vector<T> data(vec.begin(),vec.end());

				int iter = 0; 
				double median= -999;
				double lastmedian= -999;
				double rms= -999;
				double lastrms= -999;	
				double diff= 0;	
				double eps= 1.e+99;
  			//while ( (c1 >= c2) && (iter < maxiter) ){
				while ( eps>tol && (iter < maxiter) ){
  				n = (int)data.size();
					median= Utils::GetMedian(data,true);	
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
		}//close function

		static double Trace(TMatrixD T){			
			double trace= 0;
			for(int i=0;i<T.GetNrows();i++){
				trace+= T(i,i);		
			}
  		return trace;
		}

		static double GetMahalanobisDistance(TMatrixD x, TMatrixD mean, TMatrixD Sigma,bool isInverted=false){
			double dist= 0;
			if(!isInverted) Sigma.Invert();
			TMatrixD diff= x-mean;
			TMatrixD diffT(TMatrixD::kTransposed,diff);
			TMatrixD M= diffT*Sigma*diff;
			dist= sqrt(M(0,0));
			return dist;
		}//close

		/*
		static timespec TimeDiff(timespec start, timespec end) {
			timespec temp;
			if ((end.tv_nsec-start.tv_nsec)<0) {
				temp.tv_sec = end.tv_sec-start.tv_sec-1;
				temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
			} 
			else {
				temp.tv_sec = end.tv_sec-start.tv_sec;
				temp.tv_nsec = end.tv_nsec-start.tv_nsec;
			}
			return temp;
		}

		static timespec TimeSum (timespec time1, timespec time2) {
			timespec  result ;
			result.tv_sec = time1.tv_sec + time2.tv_sec ;
    	result.tv_nsec = time1.tv_nsec + time2.tv_nsec ;
    	if (result.tv_nsec >= 1000000000L) {
      	result.tv_sec++ ;  result.tv_nsec = result.tv_nsec - 1000000000L ;
    	}
    	return (result) ;
		}
		
		static double TimeToSec (timespec time){
    	return ((double)time.tv_sec + (time.tv_nsec/1.e+09)) ;
		}

		static double TimeToNSec (timespec time){
			return (time.tv_sec*1.e+09 + (double)time.tv_nsec) ;
		}
		*/

		static void StringFindAndReplace(std::string& str, const std::string& oldstr, const std::string& newstr){
  		size_t pos = 0;
  		while((pos = str.find(oldstr, pos)) != std::string::npos){
    		str.replace(pos, oldstr.length(), newstr);
     		pos += newstr.length();
  		}	
		}//close function


		static bool PixelToWCSCoords(Img* image,int ix,int iy,double& xpos, double& ypos);


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

		
		static int GetEntropy(double& H, TMatrixD* data, int k=5);

		static int GetSaliency(TVectorD& saliencyVect,TMatrix* data,int k=5);

		static std::vector<double> ComputePageRank(TMatrixD M,double d=0.85,double tol=1.e-4);
	

	private:
	
	
};

#endif 
