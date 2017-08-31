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

#include <SysUtils.h>
#include <Logger.h>

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

//OpenMP headers
#ifdef OPENMP_ENABLED
  #include <omp.h>
	#include <parallel/partition.h>//C++ stdlib parallel algorithms
#endif



using namespace std;

namespace Caesar {

template <typename T>  
class ClippedStats : public TObject {

	public:
		T median;
		T mean;
		T stddev;

	public:	
		/** 
		\brief Constructor
 		*/
		ClippedStats(T _median,T _mean,T _stddev)
			: median(_median), mean(_mean), stddev(_stddev)
		{}

		/** 
		\brief Constructor
 		*/
		ClippedStats(){
			median= T(0);
			mean= T(0);
			stddev= T(0);
		}

		/** 
		\brief Destructor
 		*/
		virtual ~ClippedStats(){}

		/** 
		\brief Operator ==
 		*/
		ClippedStats& operator=(const ClippedStats &m) { 
  		if (this != &m) ((ClippedStats&)m).Copy(*this);
 		 	return *this;
		}

		/** 
		\brief Copy object
 		*/
		void Copy(TObject& obj) const {
			((ClippedStats&)obj).median= median;
			((ClippedStats&)obj).mean= mean;
			((ClippedStats&)obj).stddev= stddev;
		}

	ClassDef(ClippedStats,1)

};//close ClippedStats class

template <typename T>  
class StatMoments : public TObject {

	public:
		T N;
		T M1;
		T M2;
		T M3;
		T M4;
		T minVal;
		T maxVal;

	public:
		/** 
		\brief Constructor
 		*/
		StatMoments(){
			Reset();
		}
		/** 
		\brief Destructor
 		*/
		virtual ~StatMoments(){}

		/** 
		\brief Operator ==
 		*/
		StatMoments& operator=(const StatMoments &m) { 
  		if (this != &m) ((StatMoments&)m).Copy(*this);
 		 	return *this;
		}

		/** 
		\brief Copy object
 		*/
		void Copy(TObject& obj) const {
			((StatMoments&)obj).N= N;
			((StatMoments&)obj).M1= M1;
			((StatMoments&)obj).M2= M2;
			((StatMoments&)obj).M3= M3;
			((StatMoments&)obj).M4= M4;
			((StatMoments&)obj).minVal= minVal;
			((StatMoments&)obj).maxVal= maxVal;
		}


	public:
		/** 
		\brief Reset moments
 		*/
		void Reset(){
			N= T(0);
			M1= T(0);
			M2= T(0);
			M3= T(0);
			M4= T(0);
			minVal= std::numeric_limits<T>::max();
			maxVal= -std::numeric_limits<T>::max();
		}
		/** 
		\brief Print moments
 		*/
		void Print(){
			cout<<"*** STAT MOMENTS ***"<<endl;
			cout<<"N="<<N<<" min/max="<<minVal<<"/"<<maxVal<<endl;
			cout<<"M1: "<<M1<<", M2="<<M2<<endl;
			cout<<"M3: "<<M3<<", M4="<<M4<<endl;
			cout<<"*****************"<<endl;
		}

	ClassDef(StatMoments,1)

};//close StatMoments()


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
	
		/**
		* \brief Compute median using nth_element (should run in O(n))
		*/
		template < typename T >
		static T GetMedianFast(std::vector<T>& vec, bool useParallelVersion=false){

			//Check empty vector
			if(vec.size()<=0) {
				WARN_LOG("Empty data vector, returning zero!");
				return 0;
			}
			std::size_t n = vec.size()/2;
  		
			#ifdef OPENMP_ENABLED
				if(useParallelVersion) __gnu_parallel::__parallel_nth_element(vec.begin(),vec.begin()+n, vec.end(), std::less<T>());
				else nth_element(vec.begin(),vec.begin()+n, vec.end());
			#else 
				nth_element(vec.begin(),vec.begin()+n, vec.end());
			#endif
			
			double median= vec[n];//Odd number of elements
    	if(n%2==0){//Even numbers
				double a= *std::max_element(vec.begin(),vec.begin()+n);
				double median_even= (median+a)/2.;
				median= median_even;
			}

  		return median;
		}//close GetMedianFast()

		/**
		* \brief Compute MAD using nth_element (should run in O(n))
		*/
		template < typename T >
		static T GetMADFast( std::vector<T>const &vec, T median,bool useParallelVersion=false){
			size_t n = vec.size();
			if(n<=0) {
				WARN_LOG("Empty data vector, returning zero!");
				return 0;
			}
			std::vector<T> MADs;
			MADs.resize(n);
			
			#ifdef OPENMP_ENABLED
			#pragma omp parallel for
			#endif
			for(size_t j=0;j<n;j++) MADs[j]= fabs(vec[j]-median);
			/*
			#ifdef OPENMP_ENABLED
				#pragma omp declare reduction (merge : std::vector<T> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
				#pragma omp parallel for reduction(merge: MADs)
				for(size_t j=0;j<n;j++) MADs.push_back( fabs(vec[j]-median) );
			#else
  			for(size_t j=0;j<n;j++) MADs.push_back( fabs(vec[j]-median) );
			#endif
			*/

			T MAD= StatsUtils::GetMedianFast(MADs,useParallelVersion);
  		return MAD;
		}//close GetMADFast()


		/**
		* \brief Compute median using vector sorting (should run in O(nlog n))
		*/
		template < typename T >
		static T GetMedian( std::vector<T>&vec, bool isSorted=false){//this should run in O(nlog(n))
			size_t n = vec.size();
			if(n<=0) {
				WARN_LOG("Empty data vector, returning zero!");
				return 0;
			}
  		if(!isSorted) std::sort(vec.begin(), vec.end());
			double median= 0;		
  		if(n%2==0) median = (vec[n/2-1] + vec[n/2])/2;	
			else median = vec[n/2];
  		return median;
		}//close GetMedian()

		/**
		* \brief Compute MAD using vector sorting (should run in O(nlog n))
		*/
		template < typename T >
		static T GetMAD( std::vector<T>const &vec, T median){
			size_t n = vec.size();
			if(n<=0) {
				WARN_LOG("Empty data vector, returning zero!");
				return 0;
			}
			std::vector<double> MADs;
  		for(int j=0;j<n;j++){
				double diff= fabs(vec[j]-median);
				MADs.push_back(diff);
			}
			std::sort(MADs.begin(),MADs.end());
			double MAD= StatsUtils::GetMedian(MADs,true);
  		return MAD;
		}//close GetMAD()

		
		/**
		* \brief Compute mean & std dev
		*/
		template <typename T>
		static void ComputeMeanAndRMS(T& mean,T& stddev,std::vector<T>const &vec){

			mean= stddev= 0;
			if(vec.empty()) return;
			size_t n= vec.size();

			#ifdef OPENMP_ENABLED
				T sum= T(0);
				T s2= T(0);
				#pragma omp parallel for reduction(+: sum)
				for(size_t i=0;i<vec.size();i++) sum+= vec[i];
				mean= sum/n;

				#pragma omp parallel for reduction(+: s2)
				for(size_t i=0;i<vec.size();i++) s2 += (vec[i] - mean) * (vec[i] - mean);
    		stddev= sqrt(s2/(n-1));

			#else 
				mean= std::accumulate(vec.begin(), vec.end(), T(0))/n;
				T s2 = T(0);
    		for(auto x : vec) s2 += (x - mean) * (x - mean);
    		stddev= sqrt(s2/(n-1));
			#endif

		}//close ComputeMeanAndRMS()

		/**
		* \brief Compute clipped estimators
		*/
		template <typename T>
		static int UpdateClippedStats(ClippedStats<T>& stats_clipped,std::vector<T>& vec_clipped, std::vector<T>const& vec, ClippedStats<T>& stats, double clipsig=3, bool useParallelVersion=false){
			
			//Check if empty
			if(vec.empty()) return -1;
			
			//Compute mean & standard deviation
			//T mean= stats.mean;
			T stddev= stats.stddev;
			T median= stats.median;

			//Fill vector of truncated items	
			vec_clipped.clear();	
			for(auto x : vec) {
				T diff= fabs(x-median);
				if(diff<clipsig*stddev) vec_clipped.push_back(diff);
			}

			//Check if there are still data
			if(vec_clipped.empty()){
				INFO_LOG("Clipped data is empty, return.");
				return -1;
			}

			/*
			#ifdef OPENMP_ENABLED
				#pragma omp declare reduction (merge : std::vector<T> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
				#pragma omp parallel for reduction(merge: vec_clipped)
				for(size_t i=0;i<vec.size();i++) {
					T diff= fabs(vec[i]-median);
					if(diff<clipsig*stddev) vec_clipped.push_back(diff);
				}
				
			#else
				for(auto x : vec) {
					T diff= fabs(x-median);
					if(diff<clipsig*stddev) vec_clipped.push_back(diff);
				}
			#endif
			*/

			//Compute mean/stddev of clipped data
			T mean_clipped= T(0);
			T stddev_clipped= T(0);
			ComputeMeanAndRMS(mean_clipped,stddev_clipped,vec_clipped);
			
			//Compute median of clipped data
			T median_clipped= GetMedianFast(vec_clipped,useParallelVersion);	
			
			//Set stats				
			stats_clipped.median= median_clipped;
			stats_clipped.mean= mean_clipped;
			stats_clipped.stddev= stddev_clipped;
			
			return 0;

		}//close UpdateClippedStats()

		/**
		* \brief Compute clipped stat estimators
		*/
		template <typename T>
		static int GetClippedEstimators(ClippedStats<T>& stats, std::vector<T>const &vec, T median, T mean, T stddev, double clipsig=3, int maxiter=5, double tol=0.1, bool useParallelVersion=false){

			//Check if empty			
			if(vec.empty()) return -1;

			//Init stats & vector
			ClippedStats<T> stats_pre(median,mean,stddev);	
			std::vector<T> data_pre= vec;

			//Start iteration
			int iter = 0; 
			double eps= 1.e+99;
			

			while ( eps>tol && (iter < maxiter) ){

				//Update clipped stats
				std::vector<T> data;
				if(UpdateClippedStats(stats,data,data_pre,stats_pre,clipsig,useParallelVersion)<0){
					DEBUG_LOG("Stop clipping iteration as clipped vector has no more data!");
					stats= stats_pre;	
					data= data_pre;				
					break;
				}

				//Update tolerance
				eps= fabs(stats.stddev-stats_pre.stddev)/stats_pre.stddev;

				//Update stats pre and iter
				stats_pre= stats;
				data_pre= data;
				iter++;

			}//end iter loop	

			return 0;

		}//close GetClippedEstimators()

		/*
		template < typename T >
			static std::pair<T,T> GetClippedEstimators( std::vector<T>const &vec, double clipsig=3, int maxiter=5, double tol=0.1, bool isSorted=false){
				size_t n= vec.size();
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
  				n = data.size();
					median= GetMedian(data,true);
					//median= GetMedianFast(data,false);	
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
					n= data_tmp.size();
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
			*/

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

		template<typename T>
		static void UpdateMoments(StatMoments<T>& moments,T w){

			//Update moments
			if(w<moments.minVal) moments.minVal= w;
			if(w>moments.maxVal) moments.maxVal= w;

			moments.N++;
  		T delta = w - moments.M1;
  		T delta_n = delta/moments.N;
  		T delta_n2 = delta_n * delta_n;
  		T f = delta * delta_n * (moments.N-1);
  		moments.M1+= delta_n;
  		moments.M4+= f * delta_n2 * (moments.N*moments.N - 3*moments.N + 3) + 6 * delta_n2 * moments.M2 - 4 * delta_n * moments.M3;
  		moments.M3+= f * delta_n * (moments.N - 2) - 3 * delta_n * moments.M2;
  		moments.M2+= f;
	
		}//close UpdateMoments()


		template<typename T>
		static int ComputeMomentsFromParallel(StatMoments<T>& moments,std::vector<StatMoments<T>>const& parallel_moments)
		{
			//Check size	
			if(parallel_moments.empty()) return -1;

			//Compute moments
			if(parallel_moments.size()==1){
				moments.N= parallel_moments[0].N;
				moments.M1= parallel_moments[0].M1;
				moments.M2= parallel_moments[0].M2;
				moments.M3= parallel_moments[0].M3;
				moments.M4= parallel_moments[0].M4;
				moments.minVal= parallel_moments[0].minVal;
				moments.maxVal= parallel_moments[0].maxVal;
				return 0;
			}
		
			//Compute mean
			T N= T(0);
			T M1= T(0);
			T minVal= parallel_moments[0].minVal;
			T maxVal= parallel_moments[0].maxVal;
			for(size_t j=0;j<parallel_moments.size();j++){
				T minVal_j= parallel_moments[j].minVal;
				T maxVal_j= parallel_moments[j].maxVal;
				if(minVal_j<minVal) minVal= minVal_j;
				if(maxVal_j>maxVal) maxVal= maxVal_j;			
				T N_j= parallel_moments[j].N;
				T M1_j= parallel_moments[j].M1;
				N+= N_j;
				M1+= N_j*M1_j;
			}//end loop
			M1/= N;
		
			//Compute second moment: sum (M2_j + n_j*(mean-mean_j)^2 
			T M2= T(0);
			T mean= M1;
			for(size_t j=0;j<parallel_moments.size();j++){
				T N_j= parallel_moments[j].N;
				T mean_j= parallel_moments[j].M1;
				T M2_j= parallel_moments[j].M2;
				M2+= M2_j + N_j*(mean-mean_j)*(mean-mean_j);	
			}
	
			//Compute 3rd & 4th moments	
			T M1_A= parallel_moments[0].M1;
			T M2_A= parallel_moments[0].M2;
			T M3_A= parallel_moments[0].M3;
			T M4_A= parallel_moments[0].M4;
			T N_A= parallel_moments[0].N;
			T M1_AB= 0;
			T M2_AB= 0;
			T M3_AB= 0;
			T M4_AB= 0;
			T N_AB= 0;
			T M3= T(0);
			T M4= T(0);
			for(size_t j=1;j<parallel_moments.size();j++){
				T M1_B= parallel_moments[j].M1;
				T M2_B= parallel_moments[j].M2;
				T M3_B= parallel_moments[j].M3;
				T M4_B= parallel_moments[j].M4;
				T N_B= parallel_moments[j].N;
				T delta= M1_B-M1_A;
				N_AB= N_A+N_B;
				M1_AB= (N_A*M1_A+N_B*M1_B)/N_AB;
				//M2_AB= M2_A + M2_B + delta*delta*N_A*N_B/N_AB;
				M2_AB= M2_A + M2_B + N_A*(M1_AB-M1_A)*(M1_AB-M1_A) + N_B*(M1_AB-M1_B)*(M1_AB-M1_B);			
				M3_AB= M3_A + M3_B + pow(delta,3)*N_A*N_B*(N_A-N_B)/(N_AB*N_AB) + 3*delta*(N_A*M2_B-N_B*M2_A)/N_AB;
				M4_AB= M4_A + M4_B + pow(delta,4)*N_A*N_B*(N_A*N_A-N_A*N_B+N_B*N_B)/pow(N_AB,3) + 6*delta*delta*(N_A*N_A*M2_B+N_B*N_B*M2_A)/(N_AB*N_AB) + 4*delta*(N_A*M3_B-N_B*M3_A)/N_AB;

				//Update A partition to A+B
				N_A= N_AB;
				M1_A= M1_AB;
				M2_A= M2_AB;
				M3_A= M3_AB;
				M4_A= M4_AB;
			}//end loop partitions
		
			M3= M3_AB;
			M4= M4_AB;
					
			moments.N= N;
			moments.M1= M1;
			moments.M2= M2;
			moments.M3= M3;
			moments.M4= M4;
			moments.minVal= minVal;
			moments.maxVal= maxVal;
			
			return 0;

		}//close ComputeMomentsFromParallel()


		template<typename T,typename K> 
		static int ComputeStatsMoments(StatMoments<T>& moments,std::vector<K>const& data,bool skipNegativeValues=false)
		{
			if(data.empty()) return 0;
		
			#ifdef OPENMP_ENABLED
				std::vector<StatMoments<T>> moments_parallel;
				T N_t= T(0);
				T M1_t= T(0);
				T M2_t= T(0);
				T M3_t= T(0);
				T M4_t= T(0);
				T minVal_t= std::numeric_limits<T>::max();
				T maxVal_t= -std::numeric_limits<T>::max();

				#pragma omp parallel private(N_t,M1_t,M2_t,M3_t,M4_t,minVal_t,maxVal_t)
				{

					int thread_id= omp_get_thread_num();
					int nthreads= SysUtils::GetOMPThreads();
					DEBUG_LOG("Starting image moment computing in thread "<<thread_id<<" (nthreads="<<nthreads<<") ...");
			
					//Init moments
					minVal_t= std::numeric_limits<T>::max();
					maxVal_t= -std::numeric_limits<T>::max();
					N_t= T(0);
					M1_t= T(0);
					M2_t= T(0);
					M3_t= T(0);
					M4_t= T(0);

					#pragma omp single
   				{						
						moments_parallel.assign(nthreads,StatMoments<T>());
   				}

					#pragma omp for
					for(size_t i=0;i<data.size();i++){			
						T w= static_cast<T>(data[i]);
						if( w==0 || (skipNegativeValues && w<0) ) continue; 
						if(w<minVal_t) minVal_t= w;
						if(w>maxVal_t) maxVal_t= w;
						N_t++;
  					T delta = w - M1_t;
  					T delta_n = delta/N_t;
  					T delta_n2 = delta_n * delta_n;
  					T f = delta * delta_n * (N_t-1);
  					M1_t+= delta_n;
  					M4_t+= f * delta_n2 * (N_t*N_t - 3*N_t + 3) + 6 * delta_n2 * M2_t - 4 * delta_n * M3_t;
  					M3_t+= f * delta_n * (N_t - 2) - 3 * delta_n * M2_t;
  					M2_t+= f;
					}//end loop vector

					DEBUG_LOG("Thread id="<<omp_get_thread_num()<<": N="<<N_t<<", M1="<<M1_t<<", M2="<<M2_t<<", M3="<<M3_t<<", M4="<<M4_t);

					//Fill moments for this thread
					moments_parallel[thread_id].N= N_t;
					moments_parallel[thread_id].M1= M1_t;
					moments_parallel[thread_id].M2= M2_t;
					moments_parallel[thread_id].M3= M3_t;
					moments_parallel[thread_id].M4= M4_t;
					moments_parallel[thread_id].minVal= minVal_t;
					moments_parallel[thread_id].maxVal= maxVal_t;
					
				}//close parallel section
			
				//Compute aggregated moments
				if(ComputeMomentsFromParallel(moments,moments_parallel)<0){
					ERROR_LOG("Failed to aggregate moments from parallel estimates!");
					return -1;
				}

			#else
				T N= T(0);
				T M1= T(0);
				T M2= T(0);
				T M3= T(0);
				T M4= T(0);
				T minVal= std::numeric_limits<T>::max();
				T maxVal= -std::numeric_limits<T>::max();

				for(size_t i=0;i<data.size();i++){			
					T w= static_cast<T>(data[i]);
					if( w==0 || (skipNegativeValues && w<0) ) continue; 
					if(w<minVal) minVal= w;
					if(w>maxVal) maxVal= w;

					N++;
  				T delta = w - M1;
  				T delta_n = delta/N;
  				T delta_n2 = delta_n * delta_n;
  				T f = delta * delta_n * (N-1);
  				M1+= delta_n;
  				M4+= f * delta_n2 * (N*N - 3*N + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3;
  				M3+= f * delta_n * (N - 2) - 3 * delta_n * M2;
  				M2+= f;
				}//end loop pixels

				moments.N= N;
				moments.M1= M1;
				moments.M2= M2;
				moments.M3= M3;
				moments.M4= M4;
				moments.minVal= minVal;
				moments.maxVal= maxVal;
			#endif

			return 0;

		}//close ComputeStatsMoments()

		/** 
		\brief Compute Mahalanobis distance of input matrix
 		*/
		static double GetMahalanobisDistance(TMatrixD x, TMatrixD mean, TMatrixD Sigma,bool isInverted=false);

		
		/** 
		\brief Compute page rank of input matrix
 		*/
		static int ComputePageRank(std::vector<double>& ranks,TMatrixD& M,double d=0.85,double tol=1.e-4);

		/** 
		\brief Compute norm mean after linear transformation from [DataMin,DataMax] to [NormMin,NormMax]
 		*/
		template<typename T>
		static T ComputeNormMean(T mean,T DataMin,T DataMax,T NormMin,T NormMax){
			T A= NormMin - (NormMax-NormMin)*DataMin/(DataMax-DataMin);
			T B= (NormMax-NormMin)/(DataMax-DataMin);
			T mean_transf= A + B*mean;//mean under linear transformation
			return mean_transf;
		} 

		/** 
		\brief Compute norm difference between two values after linear transformation from [DataMin,DataMax] to [NormMin,NormMax]
 		*/
		template<typename T>
		static T ComputeNormDiff(T diff,T DataMin,T DataMax,T NormMin,T NormMax){
			T B= (NormMax-NormMin)/(DataMax-DataMin);
			T diff_transf= B*diff;//difference under linear transformation
			return diff_transf;
		}

		template<typename T>
		static T ComputeNormDiffSqr(T diff2,T DataMin,T DataMax,T NormMin,T NormMax){
			T B= (NormMax-NormMin)/(DataMax-DataMin);
			T diff2_transf= B*B*diff2;//difference under linear transformation
			return diff2_transf;
		} 

		/** 
		\brief Get mean/rms after linear transformation from [DataMin,DataMax] to [NormMin,NormMax]
 		*/
		template<typename T>
		static void ComputeNormMeanAndRMS(T& mean_transf,T& rms_transf,T mean,T rms,T DataMin,T DataMax,T NormMin,T NormMax){
			T A= NormMin - (NormMax-NormMin)*DataMin/(DataMax-DataMin);
			T B= (NormMax-NormMin)/(DataMax-DataMin);
			mean_transf= A + B*mean;//mean under linear transformation
			rms_transf= B*rms;
		} 

		/** 
		\brief Get value normalized in a range
 		*/
		template<typename T>
		static T GetNormValue(T x,T xmin,T xmax,T NormMin,T NormMax){
			T xnorm= NormMin + (NormMax-NormMin)*(x-xmin)/(xmax-xmin);
			return xnorm;
		}

		/** 
		\brief Normalize a vector	in a range
 		*/
		template<typename T>	
		static void NormalizeVector(std::vector<T>& data,T NormMin=0,T NormMax=1){
			//Find vector min & max
			T xmax = *max_element(data.begin(), data.end());
			T xmin = *min_element(data.begin(), data.end());

			//Transform vector
			std::transform(data.begin(), data.end(), data.begin(), 
   			[&NormMin,&NormMax,&xmin,&xmax](T x) -> T { return GetNormValue(x,xmin,xmax,NormMin,NormMax); }
			);
		}
		

	private:
	
		ClassDef(StatsUtils,1)
};

#ifdef __MAKECINT__
#pragma link C++ class StatsUtils+;

#pragma link C++ class Caesar::ClippedStats<int>+;
#pragma link C++ class Caesar::ClippedStats<long int>+;
#pragma link C++ class Caesar::ClippedStats<float>+;
#pragma link C++ class Caesar::ClippedStats<double>+;

#pragma link C++ class Caesar::StatMoments<int>+;
#pragma link C++ class Caesar::StatMoments<long int>+;
#pragma link C++ class Caesar::StatMoments<float>+;
#pragma link C++ class Caesar::StatMoments<double>+;

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
