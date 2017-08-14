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
* @file Image.cc
* @class Image
* @brief Image
*
* Image class
* @author S. Riggi
* @date 20/01/2015
*/

#include <Image.h>
#include <ImgStats.h>
#include <ImgMetaData.h>

//== Img headers
#include <ImgStats.h>

//== IO headers
#include <FITSReader.h>
#include <FITSWriter.h>

//== Utils headers
#include <StatsUtils.h>
#include <GraphicsUtils.h>
#include <CodeUtils.h>
#include <SysUtils.h>

//== Img processing headers
#include <BkgFinder.h>
#include <BkgData.h>
#include <Source.h>
#include <BlobFinder.h>
#include <ChanVeseSegmenter.h>
#include <Pixel.h>

//== Filter headers
#include <GuidedFilter.h>
#include <WTFilter.h>
#include <KirschFilter.h>
#include <LoGFilter.h>
#include <GradientFilter.h>
#include <MorphFilter.h>
#include <SaliencyFilter.h>

#include <Logger.h>

#include <TFile.h>
#include <TH2F.h>
#include <TMath.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TMath.h>
#include <TColor.h>
#include <TGaxis.h>
#include <TExec.h>
#include <TF1.h>

#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

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
#include <chrono>

using namespace std;


ClassImp(Caesar::Image)

namespace Caesar {

//Default constructor
Image::Image() 
{
	//Do not allocate memory in default constructor otherwise you will have a memory leak 
	//(see https://root.cern.ch/root/html534/guides/users-guide/AddingaClass.html#the-default-constructor)
  Init();
}//close costructor

Image::Image(long int nbinsx,long int nbinsy,float xlow,float ylow,std::string name)
{

	//Check mismatch between pixels size and dimx/dimy
	if(nbinsx<=0 || nbinsy<=0){
		std::stringstream ss;
		ss<<"Invalid image size (<=0) given!";
		ERROR_LOG(ss.str());
		throw std::out_of_range(ss.str().c_str());
	}

	//Init pars
  Init();

	//Set image name
	m_name= name;

	//Set image size
	SetSize(nbinsx,nbinsy,xlow,ylow);

}//close constructor


Image::Image(long int nbinsx,long int nbinsy,std::vector<float>const& pixels,float xlow,float ylow,std::string name)
{
	//Check mismatch between pixels size and dimx/dimy
	long int npixels= (long int)(pixels.size());
	if(nbinsx<=0 || nbinsy<=0 || npixels<=0 || npixels!=(nbinsx*nbinsy)){
		std::stringstream ss;
		ss<<"Invalid image size (<=0) or pixel size given (should be equal to Nx*Ny)!";
		ERROR_LOG(ss.str());
		throw std::out_of_range(ss.str().c_str());
	}

	//Init pars
  Init();

	//Set image name
	m_name= name;

	//Set image size
	SetSize(nbinsx,nbinsy,xlow,ylow);

	//Fill pixels	
	bool useNegativePixInStats= true;
	for(size_t i=0;i<pixels.size();i++){
		double w= pixels[i];
		if(FillPixel(i,w,useNegativePixInStats)<0) continue;
	}

}//close constructor


Image::Image(const Image &img)
{
	// Copy constructor
	DEBUG_LOG("Copy constuctor called...");
  ((Image&)img).Copy(*this);
}//close copy constructor



void Image::Copy(TObject& obj) const
{
	
	DEBUG_LOG("Copying main vars...");
	((Image&)obj).m_name= m_name;
	((Image&)obj).m_Nx= m_Nx;
	((Image&)obj).m_Ny= m_Ny;
	((Image&)obj).m_Xmin= m_Xmin;
	((Image&)obj).m_Ymin= m_Ymin;
	
	DEBUG_LOG("Copying pixel collection...");
	((Image&)obj).m_pixels= m_pixels;
		
	DEBUG_LOG("Copying stats vars...");
	((Image&)obj).m_HasMetaData= m_HasMetaData;
	((Image&)obj).m_HasStats= m_HasStats;
	((Image&)obj).m_Npix= m_Npix;
	((Image&)obj).m_M1= m_M1;
	((Image&)obj).m_M2= m_M2;
	((Image&)obj).m_M3= m_M3;
	((Image&)obj).m_M4= m_M4;
	((Image&)obj).m_PixelMin= m_PixelMin;
	((Image&)obj).m_PixelMax= m_PixelMax;
	
	DEBUG_LOG("Copying meta data...");
	if(m_MetaData){
		((Image&)obj).m_MetaData= new ImgMetaData;
		*((Image&)obj).m_MetaData = *m_MetaData;
	}
	
	DEBUG_LOG("Copying stats data...");
	if(m_Stats){
		((Image&)obj).m_Stats= new ImgStats;
		*((Image&)obj).m_Stats = *m_Stats;
	}
	
}//close Copy()


Image::~Image(){

	if(m_Stats){
		DEBUG_LOG("Deleting stats...");	
		delete m_Stats;
		m_Stats= 0;
	}

	if(m_MetaData){
		DEBUG_LOG("Deleting meta-data...");
		delete m_MetaData;
		m_MetaData= 0;	
	}	
 	
}//close destructor


Image& Image::operator=(const Image &img) { 
	// Operator =
  if (this != &img) ((Image&)img).Copy(*this);
  return *this;
}

//================================================================
//===    INITIALIZATION METHODS
//================================================================
void Image::Init(){
		
	//Set sizes
	m_Nx= m_Ny= 0;
	m_pixels.clear();
	m_Xmin= 0;
	m_Ymin= 0;

	//Meta-Data
	m_HasMetaData= false;
	m_MetaData= 0;
 
	//Stats
	m_HasStats= false;
	m_Stats= 0;
	ResetImgStats(true);//Reset stats & moments

}//close Init()


int Image::SetSize(long int Nx,long int Ny,float xlow,float ylow){
			
	//Check dimensions
	if(Nx<=0 || Ny<=0) {
		ERROR_LOG("Invalid (empty or negative) sizes given, nothing will be done!");
		return -1;
	}

	//Check range
	bool hasValidOffset= (xlow>=0 && ylow>=0); 
	if(!hasValidOffset){
		ERROR_LOG("Invalid xlow/ylow given (should be >=0)!");
		return -1;
	}

	//Reset & delete stats, reset moments
	ResetImgStats(true,true);

	//Clear existing vector (not needed)
	//m_pixels.clear();
		
	//Allocate new space
	long int N= Nx*Ny;
	try {
		m_pixels.resize(N);
	}
	catch(std::exception& e){
		ERROR_LOG("C++ exception occurred while allocating memory (N="<<N<<") for pixel vector!");
		return -1;
	}

	//Init sizes
	m_Nx= Nx;
	m_Ny= Ny;

	//Init min coordinates
	m_Xmin= xlow;
	m_Ymin= ylow;

	return 0;

}//close SetSize()

//================================================================
//===    READ METHODS
//================================================================
int Image::ReadFITS(std::string filename,int hdu_id,int ix_min,int ix_max,int iy_min,int iy_max){

	//Start timer		
	auto start = chrono::steady_clock::now();

	//Read image or tile
	bool check_file= true;
	Caesar::FITSFileInfo fits_info;
	int status= FITSReader::Read(*this,fits_info,filename,hdu_id,ix_min,ix_max,iy_min,iy_max,check_file);
	
	//Stop timer and print
	auto end = chrono::steady_clock::now();
	double dt= chrono::duration <double, milli> (end-start).count();
	INFO_LOG("Read FITS image "<<filename<<" in "<<dt<<" ms");	

	//Check error status
	if(status<0){
		ERROR_LOG("Failed to fill image from FITS file!");
		return -1;
	}

	return 0;

}//close ReadFITS()


Image* Image::GetTile(long int ix_min,long int ix_max,long int iy_min,long int iy_max,std::string imgname){

	//Extract pixel rectangular selection
	std::vector<float> tile_pixels;
	if(GetTilePixels(tile_pixels,ix_min,ix_max,iy_min,iy_max)<0){
		ERROR_LOG("Failed to extract pixel rectangular data, returning nullptr!");
		return nullptr;
	}

	long int TileSizeX= ix_max - ix_min + 1;
	long int TileSizeY= iy_max - iy_min + 1;
	float xlow= ix_min;
	float ylow= iy_min;
	std::string name= imgname;
	if(name=="") name= "tileimg";

	INFO_LOG("ix min/max="<<ix_min<<"/"<<ix_max<<", iy min/max="<<iy_min<<"/"<<iy_max<<" TileSizeX="<<TileSizeX<<", TileSizeY="<<TileSizeY<<", pixels size="<<tile_pixels.size());
	
	Image* tile= new Image(TileSizeX,TileSizeY,tile_pixels,xlow,ylow,name);
	tile->CopyMetaData(this->m_MetaData);
	
	return tile;

}//close GetTile()

//================================================================
//===    FILLING METHODS
//================================================================
int Image::FillPixel(long int gbin,double w,bool useNegativePixInStats){

	//Check bin & value
	if(gbin<0 || gbin>=(long int)(m_pixels.size())) {
		WARN_LOG("Invalid pixel bin ("<<gbin<<") given to be filled (hint: check if image size is not initialized or given bin exceed size)!");
		return -1;
	}
	if(TMath::IsNaN(w) || fabs(w)==TMath::Infinity()) {
		WARN_LOG("Given value is NaN or inf, skipping!");
		return -1;
	}

	//Compute binx & biny
	long int binx = GetBinX(gbin);
 	long int biny = GetBinY(gbin);

	//Check if pixel value has been already filled
	double w_old= m_pixels[gbin];
	if( w_old!=0 ){
		WARN_LOG("Pixel "<<gbin<<" ("<<binx<<","<<biny<<") has been already filled (w="<<w_old<<"), skipping...");
		return -1;
	}

	//Set pixel value
	m_pixels[gbin]= w;

	//Update moments
	if(w>=0 || (w<0 && useNegativePixInStats)) UpdateMoments(w);

	return 0;

}//close FillPixel()


int Image::FillPixel(long int ix,long int iy,double w,bool useNegativePixInStats){

	//Compute global bin
	long int gbin= GetBin(ix,iy);
	
	//Fill pixel
	return FillPixel(gbin,w,useNegativePixInStats);

}//close FillPixel()


//================================================================
//===    STATS METHODS
//================================================================
void Image::ResetImgStats(bool resetMoments,bool clearStats){

	//Reset stats
	if(m_Stats) m_Stats->Reset();

	//Reset moments
	if(resetMoments){
		m_PixelMin= +1.e+99;
		m_PixelMax= -1.e+99;
		m_Npix= 0;//npixels	
  	m_M1= 0;//1st moments
  	m_M2= 0;//2nd moment
		m_M3= 0;//3rd moment
		m_M4= 0;//4th moment
	}

	//Delete current stats data
	if(clearStats){
		ClearImgStats();
	}

}//close ResetImgStats()

void Image::ComputeMoments(bool skipNegativePixels){

	
	#ifdef OPENMP_ENABLED
		
		//Define variables for reduction (OMP does not allow class members in reduction operation)
		long long int nPix= 0;
		double pixelMin= +1.e+99; 
		double pixelMax= -1.e+99;
		double M1= 0;
		double M2= 0;
		double M3= 0;
		double M4= 0;
		double S= 0;

		long long int nPix_t= 0;
		double M1_t= 0;
		double M2_t= 0;
		double M3_t= 0;
		double M4_t= 0;
		double S_t= 0;

		std::vector<long long int> nPix_list;
		std::vector<double> M1_list;
		std::vector<double> M2_list;
		std::vector<double> M3_list;
		std::vector<double> M4_list;

		#pragma omp parallel private(nPix_t,M1_t,M2_t,M3_t,M4_t,S_t) reduction(max: pixelMax), reduction(min: pixelMin), reduction(+: nPix, S)
		{

			int thread_id= omp_get_thread_num();
			int nthreads= SysUtils::GetOMPThreads();
			INFO_LOG("Starting image moment computing in thread "<<thread_id<<" (nthreads="<<nthreads<<") ...");
			
			//Init moments
			nPix_t= 0;
			M1_t= 0;
			M2_t= 0;
			M3_t= 0;
			M4_t= 0;
			S_t= 0;

			#pragma omp single
   		{
     		nPix_list.assign(nthreads,0);
				M1_list.assign(nthreads,0);
				M2_list.assign(nthreads,0);	
				M3_list.assign(nthreads,0);
				M4_list.assign(nthreads,0);
   		}

			#pragma omp for
			for(size_t i=0;i<m_pixels.size();i++){			
				double w= m_pixels[i];
				if( w==0 || (skipNegativePixels && w<0) ) continue; 
				if(w<pixelMin) pixelMin= w;
				if(w>pixelMax) pixelMax= w;
				nPix_t++;
				S_t+= w;
  			double delta = w - M1_t;
  			double delta_n = delta/nPix_t;
  			double delta_n2 = delta_n * delta_n;
  			double f = delta * delta_n * (nPix_t-1);
  			M1_t+= delta_n;
  			M4_t+= f * delta_n2 * (nPix_t*nPix_t - 3*nPix_t + 3) + 6 * delta_n2 * M2_t - 4 * delta_n * M3_t;
  			M3_t+= f * delta_n * (nPix_t - 2) - 3 * delta_n * M2_t;
  			M2_t+= f;
			}//end loop pixels

			INFO_LOG("Thread id="<<omp_get_thread_num()<<": nPix_t="<<nPix_t<<", M1="<<M1_t<<", M2="<<M2_t<<", M3="<<M3_t<<", M4="<<M4_t);

			nPix+= nPix_t;
			//M1+= M1_t;
			//M2+= M2_t;
			//M3+= M3_t;
			//M4+= M4_t;
			S+= S_t;

			//Fill list
			nPix_list[thread_id]= nPix_t; 
			M1_list[thread_id]= M1_t; 
			M2_list[thread_id]= M2_t; 
			M3_list[thread_id]= M3_t; 
			M4_list[thread_id]= M4_t; 
	
		}//close parallel section

		
		//Compute moments
		if(M3_list.size()==1){
			M1= M1_list[0];
			M2= M2_list[0];
			M3= M3_list[0];
			M4= M4_list[0];
		}
		else{
			//Compute mean
			M1= S/(double)(nPix);
		
			//Compute second moment: sum (M2_j + n_j*(mean-mean_j)^2 
			M2= 0;
			double mean= M1;
			for(size_t j=0;j<M2_list.size();j++){
				double M2_j= M2_list[j];
				double mean_j= M1_list[j];
				double N_j= nPix_list[j];
				M2+= M2_j + N_j*(mean-mean_j)*(mean-mean_j);	
			}
	
			//Compute 3rd & 4th moments	
			double M1_A= M1_list[0];
			double M2_A= M2_list[0];
			double M3_A= M3_list[0];
			double M4_A= M4_list[0];
			double N_A= nPix_list[0];
			double M1_AB= 0;
			double M2_AB= 0;
			double M3_AB= 0;
			double M4_AB= 0;
			double N_AB= 0;
			for(size_t j=1;j<M3_list.size();j++){
				double M1_B= M1_list[j];
				double M2_B= M2_list[j];
				double M3_B= M3_list[j];
				double M4_B= M4_list[j];
				double N_B= nPix_list[j];
				double delta= M1_B-M1_A;
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
		}//close else


		//Set reduced vars to global class var
		m_Npix = nPix;
		m_PixelMin= pixelMin;
		m_PixelMax= pixelMax;
		m_M1= M1;
		m_M2= M2;
		m_M3= M3;
		m_M4= M4;

	#else
		for(size_t i=0;i<m_pixels.size();i++){			
			double w= m_pixels[i];
			if( w==0 || (skipNegativePixels && w<0) ) continue; 
			UpdateMoments(w);
		}//end loop pixels
	#endif

}//close ComputeMoments()


void Image::UpdateMoments(double w){

	//Update full image moments
	if(w<m_PixelMin) m_PixelMin= w;
	if(w> m_PixelMax) m_PixelMax= w;

	m_Npix++;
  double delta = w - m_M1;
  double delta_n = delta/m_Npix;
  double delta_n2 = delta_n * delta_n;
  double f = delta * delta_n * (m_Npix-1);
  m_M1+= delta_n;
  m_M4+= f * delta_n2 * (m_Npix*m_Npix - 3*m_Npix + 3) + 6 * delta_n2 * m_M2 - 4 * delta_n * m_M3;
  m_M3+= f * delta_n * (m_Npix - 2) - 3 * delta_n * m_M2;
  m_M2+= f;
	
}//close UpdateMoments()

void Image::ComputeStatsParams(bool computeRobustStats,bool skipNegativePixels){

	//## Reset previous stats params (Reset only Stats not moments!)
	ResetImgStats(false);


	//-- DEBUG
	DEBUG_LOG("Npix="<<m_Npix<<", M1="<<m_M1<<", m_M2="<<m_M2<<", m_M3="<<m_M3<<", m_M4="<<m_M4<<", pixel sizes="<<m_pixels.size());
	//----
	
	//## Compute Stats params
	m_Stats->n= m_Npix;
	m_Stats->min= m_PixelMin;
	m_Stats->max= m_PixelMax;
	m_Stats->mean= m_M1;
	m_Stats->rms= 0;
	if(m_Npix>2) m_Stats->rms= sqrt(m_M2/(m_Npix-1));
	m_Stats->skewness= 0;
	m_Stats->kurtosis= 0;
	if(m_M2!=0) {
  	m_Stats->skewness= sqrt(m_Npix)*m_M3/pow(m_M2,1.5);//need to adjust for finite population?
		m_Stats->kurtosis= m_Npix*m_M4/(m_M2*m_M2)-3;
	}
	

	//## Compute stats param errors
	m_Stats->meanErr= 0;
	if(m_Npix>0) m_Stats->meanErr= (m_Stats->rms)/sqrt(m_Npix);
	double varianceErr= 0;
	m_Stats->rmsErr= 0;
	if(m_Npix>1) {
		varianceErr= (m_M4-(m_Npix-3)/(m_Npix-1)*pow(m_Stats->rms,4))/m_Npix;
		m_Stats->rmsErr= varianceErr/(2*m_Stats->rms);
	}
		 
	m_Stats->skewnessErr= 0;
	m_Stats->kurtosisErr= 0;
	if(m_Npix>2) m_Stats->skewnessErr= sqrt(6.*m_Npix*(m_Npix-1)/((m_Npix-2)*(m_Npix+1)*(m_Npix+3)));//approximate for normal distribution
	
  
	//## End if no robust stats are to be computed
	if(!computeRobustStats) return;
	
	//## Remove negative values? 
	//## NB: Copy vector otherwise it is modified by sorting operation inside median and other robust estimators
	std::vector<float> pixels;
	pixels.clear();
	pixels.resize(0);
	if(skipNegativePixels){
		for(size_t i=0;i<m_pixels.size();i++){
			double pixelValue= m_pixels[i];
			if( pixelValue==0 || (skipNegativePixels && pixelValue<0) ) continue;
			pixels.push_back(pixelValue);
		}
	}

	//## Compute robust stats (median, MAD, ...)	
	//Sort and compute median for all image	
	float median= Caesar::StatsUtils::GetMedianFast<float>(pixels);
	m_Stats->median= median;

	//Compute MAD = median(|x_i-median|)
	float medianMAD= Caesar::StatsUtils::GetMADFast(pixels,median);	
	double medianRMS= medianMAD*1.4826;//0.6744888;
	m_Stats->medianRMS= medianRMS;

	
	//## Compute biweight robust estimators
	double C= 6.;
	double tol= 0.0001;
	double nmaxIter= 10;
	std::pair<float,float> biweightEstimators= Caesar::StatsUtils::GetBiWeightEstimators<float>(pixels,median,medianRMS,C,tol,nmaxIter);
	m_Stats->bwLocation= biweightEstimators.first;
	m_Stats->bwScale= biweightEstimators.second;

	//## Compute clipped estimators
	double clipSigma= 3;
	int clipMaxIter= 100;
	double clipTolerance= 0.1;
	bool useParallelVersion= false;
	float mean= m_Stats->mean;
	float rms= m_Stats->rms;
	ClippedStats<float> clipped_stats;
	Caesar::StatsUtils::GetClippedEstimators(clipped_stats,pixels,median,mean,rms,clipSigma,clipMaxIter,clipTolerance,useParallelVersion);
	m_Stats->clippedMedian= clipped_stats.median;
	m_Stats->clippedRMS= clipped_stats.stddev;
	
	
}//close ComputeStatsParams()



int Image::ComputeStats(bool computeRobustStats,bool skipNegativePixels,bool forceRecomputing){

	
	//## Start timer
	auto start = chrono::steady_clock::now();

	
	//## Check if image has already stats computed
	if(!HasStats()){
		m_Stats= new Caesar::ImgStats;
	}
	else{		
		WARN_LOG("Image has stats already computed...");
	}

	//## If recomputing is not requested (i.e. some pixels has been reset by the user, just set the stats params!
	if(!forceRecomputing){
		ComputeStatsParams(computeRobustStats,skipNegativePixels);
		m_HasStats= true;
		
		auto stop = chrono::steady_clock::now();
		double elapsed_time= chrono::duration <double, milli> (stop-start).count();
		INFO_LOG("Image stats computed in "<<elapsed_time<<" ms");
		return 0;
	}

	
	//## Recompute the moments and stats params
	DEBUG_LOG("Recomputing image stats...");

	//--> Reset stats
	ResetImgStats(true);

	//--> Recompute moments
	ComputeMoments(skipNegativePixels);

	//--> Recompute stats params
	ComputeStatsParams(true,skipNegativePixels);

	m_HasStats= true;

	//Stop timer
	auto end = chrono::steady_clock::now();
	double dt= chrono::duration <double, milli> (end-start).count();
	INFO_LOG("Image stats recomputed in "<<dt<<" ms");
	

	return 0;

}//close ComputeStats()


int Image::GetTilePixels(std::vector<float>& pixels,long int ix_min,long int ix_max,long int iy_min,long int iy_max,bool skipNegativePixels)
{

	//Check range given
	if(ix_min<0 || ix_max<0 || ix_max>=m_Nx || ix_min>=ix_max){
		ERROR_LOG("Invalid x range given!");
		return -1;
	}
	if(iy_min<0 || iy_max<0 || iy_max>=m_Ny || iy_min>=iy_max){
		ERROR_LOG("Invalid y range given!");
		return -1;
	}
	
	//Extract pixel sub vector
	pixels.clear();

	if(skipNegativePixels){

		for(long int j=iy_min;j<=iy_max;j++){
			//Get row start/end iterators
			long int gBin_min= GetBin(ix_min,j);		
			long int gBin_max= GetBin(ix_max,j);

			std::copy_if (
				m_pixels.begin() + gBin_min, 
				m_pixels.begin() + gBin_max, 
				std::back_inserter(pixels), 
				[](float w){
					return (w>=0);
				}
			);
		}//end loop rows
	}//close if

	else{

		for(long int j=iy_min;j<=iy_max;j++){
			//Get row start/end iterators
			long int gBin_min= GetBin(ix_min,j);		
			long int gBin_max= GetBin(ix_max+1,j);		

			//Extract row and append to vector
			pixels.insert(pixels.end(), m_pixels.begin() + gBin_min, m_pixels.begin() + gBin_max);

		}//end loop row
	}//close else
	
	return 0;

}//close GetTilePixels()


int Image::GetTileMeanStats(float& mean,float& stddev,long int& npix,long int ix_min,long int ix_max,long int iy_min,long int iy_max,bool skipNegativePixels)
{
	//Init
	mean= 0;
	stddev= 0;
	npix= 0;

	//Extract pixel rectangular selection
	std::vector<float> tile_pixels;
	if(GetTilePixels(tile_pixels,ix_min,ix_max,iy_min,iy_max,skipNegativePixels)<0){
		ERROR_LOG("Failed to extract pixel rectangular data!");
		return -1;
	}
	
	//Compute mean/rms
	Caesar::StatsUtils::ComputeMeanAndRMS(mean,stddev,tile_pixels);
	npix= static_cast<long int>(tile_pixels.size());

	return 0;

}//close GetTileMeanStats()
		
int Image::GetTileMedianStats(float& median,float& mad_rms,long int& npix,long int ix_min,long int ix_max,long int iy_min,long int iy_max,bool skipNegativePixels)
{

	//Init
	median= 0;
	mad_rms= 0;
	npix= 0;

	//Extract pixel rectangular selection
	std::vector<float> tile_pixels;
	if(GetTilePixels(tile_pixels,ix_min,ix_max,iy_min,iy_max,skipNegativePixels)<0){
		ERROR_LOG("Failed to extract pixel rectangular data!");
		return -1;
	}
	
	//Compute median/mad
	bool useParallelVersion= false;
	median= Caesar::StatsUtils::GetMedianFast(tile_pixels,useParallelVersion);
	mad_rms= Caesar::StatsUtils::GetMADFast(tile_pixels,median,useParallelVersion)*1.4826;
	npix= static_cast<long int>(tile_pixels.size());

	return 0;

}//close GetTileMedianStats()


int Image::GetTileClippedStats(ClippedStats<float>& clipped_stats,long int& npix,double clipSigma,long int ix_min,long int ix_max,long int iy_min,long int iy_max,bool skipNegativePixels)
{
	//Init
	npix= 0;
	
	//Extract pixel rectangular selection
	std::vector<float> tile_pixels;
	if(GetTilePixels(tile_pixels,ix_min,ix_max,iy_min,iy_max,skipNegativePixels)<0){
		ERROR_LOG("Failed to extract pixel rectangular data!");
		return -1;
	}
	
	//Compute mean, rms, median
	float mean= 0;
	float rms= 0;
	Caesar::StatsUtils::ComputeMeanAndRMS(mean,rms,tile_pixels);
	float median= Caesar::StatsUtils::GetMedianFast<float>(tile_pixels);
	
	//Compute clipped stats
	int clipMaxIter= 100;
	double clipTolerance= 0.1;
	bool useParallelVersion= false;
	Caesar::StatsUtils::GetClippedEstimators(clipped_stats,tile_pixels,median,mean,rms,clipSigma,clipMaxIter,clipTolerance,useParallelVersion);
	npix= static_cast<long int>(tile_pixels.size());

	return 0;

}//close GetTileClippedStats()
		

int Image::GetTileBiWeightStats(float& bwLocation,float& bwScale,long int& npix,double C,long int ix_min,long int ix_max,long int iy_min,long int iy_max,bool skipNegativePixels)
{
	//Init
	bwLocation= 0;
	bwScale= 0;
	npix= 0;

	//Extract pixel rectangular selection
	std::vector<float> tile_pixels;
	if(GetTilePixels(tile_pixels,ix_min,ix_max,iy_min,iy_max,skipNegativePixels)<0){
		ERROR_LOG("Failed to extract pixel rectangular data!");
		return -1;
	}
	
	//Compute median & MAD
	bool useParallelVersion= false;
	float median= Caesar::StatsUtils::GetMedianFast<float>(tile_pixels);
	float medianRMS= Caesar::StatsUtils::GetMADFast(tile_pixels,median,useParallelVersion)*1.4826;
	
	//## Compute biweight robust estimators
	double tol= 0.0001;
	double nmaxIter= 10;
	std::pair<float,float> biweightEstimators= Caesar::StatsUtils::GetBiWeightEstimators<float>(tile_pixels,median,medianRMS,C,tol,nmaxIter);
	bwLocation= biweightEstimators.first;
	bwScale= biweightEstimators.second;
	npix= static_cast<long int>(tile_pixels.size());

	return 0;

}//close GetTileBiWeightStats()

//================================================================
//===    BACKGROUND CALCULATION
//================================================================
ImgBkgData* Image::ComputeBkg(int estimator,bool computeLocalBkg,int boxSizeX,int boxSizeY, double gridStepSizeX,double gridStepSizeY,bool use2ndPass,bool skipOutliers,double seedThr,double mergeThr,int minPixels){
	
	//## Compute bkg data
	DEBUG_LOG("Using grid bkg method...");
	ImgBkgData* bkgData= BkgFinder::FindBkg(this,estimator,computeLocalBkg, boxSizeX, boxSizeY, gridStepSizeX,gridStepSizeY,use2ndPass,skipOutliers,seedThr,mergeThr,minPixels);	
	if(!bkgData){
		ERROR_LOG("Computation of local background failed for this image!");
		return 0;
	}
	return bkgData;
	
}//close ComputeBkg()


Image* Image::GetSignificanceMap(ImgBkgData* bkgData,bool useLocalBkg){

	//Check image
	if(!bkgData){
		ERROR_LOG("Null ptr to bkg data!");
		return 0;
	}

	//Integrity check for local bkg
	long int Nx= this->GetNx();
	long int Ny= this->GetNy();
	if(useLocalBkg){
		if(!bkgData->HasLocalBkg()){
			ERROR_LOG("Local bkg option requested but no local bkg data are available!");
			return 0;
		}
		if( Nx!=(bkgData->BkgMap)->GetNx() || Ny!=(bkgData->BkgMap)->GetNy() ||
				Nx!=(bkgData->NoiseMap)->GetNx() || Ny!=(bkgData->NoiseMap)->GetNy()
		){
			ERROR_LOG("Bkg/Noise maps have different size!");		
			return 0;
		}
	}//close if
		
	//Clone this image and reset content
	TString imgName= Form("%s_significance",this->GetName().c_str());
	Image* significanceMap= this->GetCloned(std::string(imgName),true,true);
	significanceMap->Reset();

	if(useLocalBkg){
		#ifdef OPENMP_ENABLED
		#pragma omp parallel for collapse(2)
		#endif
		for(long int i=0;i<Nx;i++){		
			for(long int j=0;j<Ny;j++){	
				long int gBin= GetBin(i,j);
				double w= m_pixels[gBin];											
				double bkgLevel= (bkgData->BkgMap)->GetBinContent(i,j);
				double bkgRMS= (bkgData->NoiseMap)->GetBinContent(i,j);
				if( w==0 || bkgRMS<=0) {
					DEBUG_LOG("Empty pixel or invalid bkg  ("<<i<<","<<j<<") skip it!");
					continue;
				}
				double Z= (w-bkgLevel)/bkgRMS;
				significanceMap->FillPixel(i,j,Z);
			}//end loop
		}//end loop 
	}//close if useLocalBkg
	else{
		double bkgLevel= bkgData->gBkg;
		double bkgRMS= bkgData->gNoise;
				
		#ifdef OPENMP_ENABLED
		#pragma omp parallel for collapse(2)
		#endif
		for(long int i=0;i<Nx;i++){		
			for(long int j=0;j<Ny;j++){	
				long int gBin= GetBin(i,j);
				double w= m_pixels[gBin];										
				if( w==0 || bkgRMS<=0) continue;
				
				double Z= (w-bkgLevel)/bkgRMS;
				significanceMap->FillPixel(i,j,Z);
			}//end loop
		}//end loop
	}//close else

	return significanceMap;

}//close GetSignificanceMap()

//================================================================
//===    CONVERSION METHODS
//================================================================
TH2D* Image::GetHisto2D(std::string histoname){

	//Check if image is not empty
	if(m_pixels.empty() || m_Nx<=0 || m_Ny<=0){
		WARN_LOG("Image is empty or has no data/size stored, returning nullptr!");
		return nullptr;
	}

	//Create 2D histo
	float xmin= m_Xmin - 0.5;
	float xmax= (m_Xmin+m_Nx-1) + 0.5;
	float ymin= m_Ymin - 0.5;
	float ymax= (m_Ymin+m_Ny-1) + 0.5;

	std::string hname= m_name;
	if(histoname!="") hname= histoname;
	TH2D* histo= new TH2D(hname.c_str(),hname.c_str(),m_Nx,xmin,xmax,m_Ny,ymin,ymax);
	histo->Sumw2();

	//Fill histo
	for(long int j=0;j<m_Ny;j++){
		for(long int i=0;i<m_Nx;i++){
			long int gBin= GetBin(i,j);	
			double w= m_pixels[gBin];
			histo->SetBinContent(i+1,j+1,w);
		}
	}

	return histo;

}//close GetHisto2D()



}//close namespace
