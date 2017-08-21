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
	: TNamed()
{
	//Do not allocate memory in default constructor otherwise you will have a memory leak 
	//(see https://root.cern.ch/root/html534/guides/users-guide/AddingaClass.html#the-default-constructor)
  Init();
}//close costructor

Image::Image(long int nbinsx,long int nbinsy,float xlow,float ylow,std::string name) 
	: TNamed(name.c_str(),name.c_str())
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
	: TNamed(name.c_str(),name.c_str())
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

Image::Image(long int nbinsx,long int nbinsy,float w,float xlow,float ylow,std::string name)
	: TNamed(name.c_str(),name.c_str())
{
	//Check mismatch between pixels size and dimx/dimy
	if(nbinsx<=0 || nbinsy<=0){
		std::stringstream ss;
		ss<<"Invalid image size (<=0) or pixel size given (should be equal to Nx*Ny)!";
		ERROR_LOG(ss.str());
		throw std::out_of_range(ss.str().c_str());
	}
	long int npixels= nbinsx*nbinsy;

	//Init pars
  Init();

	//Set image name
	m_name= name;

	//Set image size
	SetSize(nbinsx,nbinsy,xlow,ylow);

	//Fill pixels	
	bool useNegativePixInStats= true;
	for(long int i=0;i<npixels;i++){
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
	DEBUG_LOG("Copying parent TNamed...");
	TNamed::Copy((Image&)obj);

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
	((Image&)obj).m_StatMoments= m_StatMoments;
	/*
	((Image&)obj).m_Npix= m_Npix;
	((Image&)obj).m_M1= m_M1;
	((Image&)obj).m_M2= m_M2;
	((Image&)obj).m_M3= m_M3;
	((Image&)obj).m_M4= m_M4;
	((Image&)obj).m_PixelMin= m_PixelMin;
	((Image&)obj).m_PixelMax= m_PixelMax;
	*/

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
	m_StatMoments.Reset();
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
//===    READ/WRITE METHODS
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


int Image::WriteFITS(std::string outfilename){
		
	//## Write fits to file
	if(FITSWriter::WriteFITS(this,outfilename)<0){
		ERROR_LOG("Failed to write image to FITS file!");
		return -1;
	}
	
	return 0;

}//close WriteFITS()


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
int Image::CheckFillPixel(long int gbin,double w){

	//Check bin & value
	if(gbin<0 || gbin>=(long int)(m_pixels.size())) {
		WARN_LOG("Invalid pixel bin ("<<gbin<<") given to be filled (hint: check if image size is not initialized or given bin exceed size)!");
		return -1;
	}
	if(TMath::IsNaN(w) || fabs(w)==TMath::Infinity()) {
		WARN_LOG("Given value (w="<<w<<") for bin "<<gbin<<" is NaN or inf, skipping!");
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
	
	return 0;

}//close CheckFillPixel()

int Image::FillPixel(long int gbin,double w,bool useNegativePixInStats){

	//Check bin & value
	if(CheckFillPixel(gbin,w)<0) {
		WARN_LOG("Invalid bin given (gbin="<<gbin<<") or nan/inf value given (w="<<w<<") or pixel already filled!");
		return -1;
	}

	//Set pixel value
	m_pixels[gbin]= w;

	//Update moments
	if(w>=0 || (w<0 && useNegativePixInStats)) {
		Caesar::StatsUtils::UpdateMoments(m_StatMoments,w);
		//UpdateMoments(w);
	}

	return 0;

}//close FillPixel()


int Image::FillPixel(long int ix,long int iy,double w,bool useNegativePixInStats){

	//Compute global bin
	long int gbin= GetBin(ix,iy);
	
	//Fill pixel
	return FillPixel(gbin,w,useNegativePixInStats);

}//close FillPixel()

#ifdef OPENMP_ENABLED
int Image::FillPixelMT(Caesar::StatMoments<double>& moments,long int gbin,double w,bool useNegativePixInStats){

	//Check bin & value
	if(CheckFillPixel(gbin,w)<0) {
		WARN_LOG("Invalid bin given (gbin="<<gbin<<") or nan/inf value given (w="<<w<<") or pixel already filled!");
		return -1;
	}

	//Set pixel value
	m_pixels[gbin]= w;

	//Update moments
	if(w>=0 || (w<0 && useNegativePixInStats)) {
		Caesar::StatsUtils::UpdateMoments(moments,w);
	}

	return 0;

}//close FillPixelMT()
#endif

#ifdef OPENMP_ENABLED
int Image::FillPixelMT(Caesar::StatMoments<double>& moments,long int ix,long int iy,double w,bool useNegativePixInStats){

	//Compute global bin
	long int gbin= GetBin(ix,iy);
	
	//Fill pixel
	return FillPixelMT(moments,gbin,w,useNegativePixInStats);

}//close FillPixelMT()
#endif


void Image::FillFromMat(cv::Mat& mat,bool useNegativePixInStats){

	//Get image size
	long int nRows = mat.rows;
  long int nCols = mat.cols;
	long int Ny= this->GetNy();

	//Reset this image
	Reset();
		
	#pragma omp parallel for
	for(long int i=0;i<nRows;++i) {
		long int rowId= i;
		long int iy= Ny-1-rowId;
		double* p = mat.ptr<double>(i);
    for (long int j=0;j<nCols;++j){
			long int colId= j;
			long int ix= colId;
			double w= p[j];
			this->FillPixel(ix,iy,w,useNegativePixInStats);
    }//end loop cols
  }//end loop rows
	
}//close FillFromMat()


//================================================================
//===    STATS METHODS
//================================================================
void Image::ResetImgStats(bool resetMoments,bool clearStats){

	//Reset stats
	if(m_Stats) m_Stats->Reset();

	//Reset moments
	if(resetMoments){
		m_StatMoments.Reset();
		/*
		m_PixelMin= +1.e+99;
		m_PixelMax= -1.e+99;
		m_Npix= 0;//npixels	
  	m_M1= 0;//1st moments
  	m_M2= 0;//2nd moment
		m_M3= 0;//3rd moment
		m_M4= 0;//4th moment
		*/
	}

	//Delete current stats data
	if(clearStats){
		ClearImgStats();
	}

}//close ResetImgStats()

int Image::ComputeMoments(bool skipNegativePixels){

	//## Recompute stat moments
	//## NB: If OMP is enabled this is done in parallel and moments are aggregated to return the correct cumulative estimate 
	int status= Caesar::StatsUtils::ComputeStatsMoments(m_StatMoments,m_pixels,skipNegativePixels);
	if(status<0){
		ERROR_LOG("Failed to compute stat moments!");
	}

	/*
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
			Caesar::StatsUtils::UpdateMoments(m_StatMoments,w);
		}//end loop pixels
	#endif
	*/

	return status;

}//close ComputeMoments()

/*
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
*/

void Image::ComputeStatsParams(bool computeRobustStats,bool skipNegativePixels){

	//## Reset previous stats params (Reset only Stats not moments!)
	ResetImgStats(false);

	//-- DEBUG
	//DEBUG_LOG("Npix="<<m_Npix<<", M1="<<m_M1<<", m_M2="<<m_M2<<", m_M3="<<m_M3<<", m_M4="<<m_M4<<", pixel sizes="<<m_pixels.size());
	DEBUG_LOG("Npix="<<m_StatMoments.N<<", M1="<<m_StatMoments.M1<<", M2="<<m_StatMoments.M2<<", M3="<<m_StatMoments.M3<<", M4="<<m_StatMoments.M4<<", pixel sizes="<<m_pixels.size());
	//----
	
	//## Compute Stats params
	/*
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
	*/
	m_Stats->n= m_StatMoments.N;
	m_Stats->min= m_StatMoments.minVal;
	m_Stats->max= m_StatMoments.maxVal;
	m_Stats->mean= m_StatMoments.M1;
	m_Stats->rms= 0;
	if(m_StatMoments.N>2) m_Stats->rms= sqrt(m_StatMoments.M2/(m_StatMoments.N-1));
	m_Stats->skewness= 0;
	m_Stats->kurtosis= 0;
	if(m_StatMoments.M2!=0) {
  	m_Stats->skewness= sqrt(m_StatMoments.N)*m_StatMoments.M3/pow(m_StatMoments.M2,1.5);//need to adjust for finite population?
		m_Stats->kurtosis= m_StatMoments.N*m_StatMoments.M4/(m_StatMoments.M2*m_StatMoments.M2)-3;
	}


	//## Compute stats param errors
	/*
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
	*/
	m_Stats->meanErr= 0;
	if(m_StatMoments.N>0) m_Stats->meanErr= (m_Stats->rms)/sqrt(m_StatMoments.N);
	double varianceErr= 0;
	m_Stats->rmsErr= 0;
	if(m_StatMoments.N>1) {
		varianceErr= (m_StatMoments.M4-(m_StatMoments.N-3)/(m_StatMoments.N-1)*pow(m_Stats->rms,4))/m_StatMoments.N;
		m_Stats->rmsErr= varianceErr/(2*m_Stats->rms);
	}	 
	m_Stats->skewnessErr= 0;
	m_Stats->kurtosisErr= 0;
	if(m_StatMoments.N>2) m_Stats->skewnessErr= sqrt(6.*m_StatMoments.N*(m_StatMoments.N-1)/((m_StatMoments.N-2)*(m_StatMoments.N+1)*(m_StatMoments.N+3)));//approximate for normal distribution
	

  
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
	TString imgName= Form("%s_significance",m_name.c_str());
	Image* significanceMap= this->GetCloned(std::string(imgName),true,true);
	significanceMap->Reset();

	
	#ifdef OPENMP_ENABLED
		Caesar::StatMoments<double> moments_t;	
		std::vector<Caesar::StatMoments<double>> parallel_moments;
	
		//#pragma omp declare reduction (merge : std::vector<Caesar::StatMoments<double>> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
		#pragma omp parallel private(moments_t) //reduction(merge: parallel_moments)
		{
			int thread_id= omp_get_thread_num();
			int nthreads= SysUtils::GetOMPThreads();

			#pragma omp single
   		{
     		parallel_moments.assign(nthreads,Caesar::StatMoments<double>());
   		}

			if(useLocalBkg){
				#pragma omp for collapse(2)
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
						significanceMap->FillPixelMT(moments_t,i,j,Z);
					}//end loop
				}//end loop 

				//parallel_moments.push_back(moments_t);
				parallel_moments[thread_id]= moments_t;
			}//close if local bkg

			else{//global bkg
				double bkgLevel= bkgData->gBkg;
				double bkgRMS= bkgData->gNoise;
						
				#pragma omp for collapse(2)
				for(long int i=0;i<Nx;i++){		
					for(long int j=0;j<Ny;j++){	
						long int gBin= GetBin(i,j);
						double w= m_pixels[gBin];										
						if( w==0 || bkgRMS<=0) continue;
				
						double Z= (w-bkgLevel)/bkgRMS;
						significanceMap->FillPixelMT(moments_t,i,j,Z);
					}//end loop
				}//end loop

				//parallel_moments.push_back(moments_t);
				parallel_moments[thread_id]= moments_t;

			}//close else 
			
		}//close parallel section

		//Update moments from parallel estimates
		if(Caesar::StatsUtils::ComputeMomentsFromParallel(m_StatMoments,parallel_moments)<0){
			ERROR_LOG("Failed to compute cumulative moments from parallel estimates (NB: image will have wrong moments!)");
			if(significanceMap){
				delete significanceMap;
				significanceMap= 0;
			}
			return nullptr;
		}

	#else
			
		if(useLocalBkg){
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

	#endif

	return significanceMap;

}//close GetSignificanceMap()

//=======================================================
//==          SOURCE EXTRACTION
//=======================================================
int Image::FindCompactSource(std::vector<Source*>& sources,Image* floodImg,ImgBkgData* bkgData,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool mergeBelowSeed,bool findNestedSources,double nestedBlobThreshold)
{

	//Find sources
	int status= BlobFinder::FindBlobs(this,sources,floodImg,bkgData,seedThr,mergeThr,minPixels,findNegativeExcess,mergeBelowSeed);
	if(status<0){
		ERROR_LOG("Blob finder failed!");
		for(unsigned int k=0;k<sources.size();k++){
			if(sources[k]){
				delete sources[k];
				sources[k]= 0;
			}	
		}//end loop sources
		sources.clear();
		return -1;
	}

	//Find nested sources?
	if(findNestedSources && sources.size()>0){
		int status= FindNestedSource(sources,bkgData,minPixels,nestedBlobThreshold);
		if(status<0){
			WARN_LOG("Nested source search failed!");
		}
	}//close if

	return 0;

}//close FindCompactSource()


int Image::FindNestedSource(std::vector<Source*>& sources,ImgBkgData* bkgData,int minPixels,double nestedBlobThreshold){

	//Check if given mother source list is empty
	int nSources= static_cast<int>(sources.size());
	if(nSources<=0){
		WARN_LOG("Empty source list given!");
		return 0;
	}

	//Find image mask of found sources
	Image* sourceMask= this->GetSourceMask(sources,false);
	if(!sourceMask){
		ERROR_LOG("Null ptr to computed source mask!");
		return -1;
	}

	Image* curvMap= this->GetLaplacianImage(true);
	if(!curvMap){
		ERROR_LOG("Null ptr to computed curvature mask!");
		if(sourceMask) {
			delete sourceMask;
			sourceMask= 0;
		}
		return -1;
	}
	curvMap->ComputeStats(true,false,false);
	double curvMapRMS= curvMap->GetPixelStats()->medianRMS;
	double curvMapThr= curvMapRMS*nestedBlobThreshold;
	Image* blobMask= curvMap->GetBinarizedImage(curvMapThr);
	if(!blobMask){
		ERROR_LOG("Failed to compute blob mask!");
		if(sourceMask) {
			delete sourceMask;
			sourceMask= 0;
		}
		if(curvMap) {
			delete curvMap;
			curvMap= 0;
		}
		return -1;
	}

	//Find blob+source mask
	Image* sourcePlusBlobMask= sourceMask->GetMask(blobMask,true);
	if(!sourcePlusBlobMask){
		ERROR_LOG("Failed to compute (source+blob) mask!");
		if(sourceMask) {
			delete sourceMask;
			sourceMask= 0;
		}
		if(curvMap) {
			delete curvMap;
			curvMap= 0;
		}
		if(blobMask) {
			delete blobMask;
			blobMask= 0;
		}
		return -1;
	}

	//Find nested blobs 
	std::vector<Source*> NestedSources;
	double fgValue= 1;
	int status= BlobFinder::FindBlobs(this,NestedSources,sourcePlusBlobMask,bkgData,fgValue,fgValue,minPixels,false,false);
	if(status<0){
		ERROR_LOG("Nested blob finder failed!");
		if(sourceMask) {
			delete sourceMask;
			sourceMask= 0;
		}
		if(curvMap) {
			delete curvMap;
			curvMap= 0;
		}
		if(blobMask) {
			delete blobMask;
			blobMask= 0;
		}
		if(sourcePlusBlobMask) {
			delete sourcePlusBlobMask;
			sourcePlusBlobMask= 0;
		}
		for(unsigned int k=0;k<NestedSources.size();k++){
			if(NestedSources[k]){
				delete NestedSources[k];
				NestedSources[k]= 0;
			}	
		}//end loop sources
		NestedSources.clear();
		return -1;
	}//close if

	//Add nested sources to mother source
	int nNestedSources= (int)NestedSources.size();
	if(nNestedSources>=0){
		INFO_LOG("#"<<nNestedSources<<" nested sources found!");

		//## Find matching between mother and nested sources
		for(int j=0;j<nNestedSources;j++){
			bool isMotherFound= false;
			DEBUG_LOG("Finding matching for nested source no. "<<j);
			
			for(int i=0;i<nSources;i++){
				int sourceId= sources[i]->Id;
				bool isInside= NestedSources[j]->IsInsideSource(sources[i]);
				if(isInside){
					DEBUG_LOG("Nested source no. "<<j<<" added to source id="<<sourceId<<" ...");
					NestedSources[j]->ComputeStats();
					NestedSources[j]->ComputeMorphologyParams();
					sources[i]->AddNestedSource(NestedSources[j]);
					isMotherFound= true;
					break;
				}
			}//end loop mother sources
			if(!isMotherFound){
				WARN_LOG("Cannot find mother source for nested source no. "<<j<<"!");
				NestedSources[j]->Print();
			}			
		}//end loop nested sources							
	}//close nNestedBlobs>0

	//Clear
	if(sourceMask) {
		delete sourceMask;
		sourceMask= 0;
	}
	if(curvMap) {
		delete curvMap;
		curvMap= 0;
	}
	if(blobMask) {
		delete blobMask;
		blobMask= 0;
	}
	if(sourcePlusBlobMask) {
		delete sourcePlusBlobMask;
		sourcePlusBlobMask= 0;
	}

	return 0;

}//close FindNestedSources()

int Image::FindExtendedSource_CV(std::vector<Source*>& sources,ImgBkgData* bkgData,int minPixels,bool findNegativeExcess,double dt,double h,double lambda1,double lambda2,double mu,double nu,double p){

	//## Compute segmented image
	Image* segmentedImg= ChanVeseSegmenter::FindSegmentation(this,false,dt,h,lambda1,lambda2,mu,nu,p);
	if(!segmentedImg){
		ERROR_LOG("Failed to compute ChanVese image segmentation!");
		return -1;
	}
	
	//## Finding blobs in masked image
	double fgValue= 1;	
	int status= this->FindCompactSource(sources,segmentedImg,bkgData,fgValue,fgValue,minPixels,false,false,false);
	if(status<0){
		ERROR_LOG("Finding sources in Chan-Vese segmented mask failed!");
		return -1;
	}
	
	//## Clear segmented image
	if(segmentedImg){
		delete segmentedImg;
		segmentedImg= 0;
	}

	//## Remove sources of negative excess 
	if(!findNegativeExcess){

		//Compute image stats if not already present
		if(!this->HasStats()){
			this->ComputeStats(true,false,false);
		}	
		if(!m_Stats){
			ERROR_LOG("Failed to compute image stats, clearing all detected sources and returning empty list!");
			for(size_t k=0;k<sources.size();k++){
				if(sources[k]){
					delete sources[k];
					sources[k]= 0;
				}	
			}//end loop sources
			sources.clear();
			return -1;
		}//close if

		//Find and remove sources with flux below the input image median (THIS COULD BE IMPROVED!)
		std::vector<int> sourcesToBeRemoved;		
		int imgMedian= m_Stats->median;
		for(size_t k=0;k<sources.size();k++){
			//Tag sources as extended
			sources[k]->SetType(Source::eExtended);

			//Check if source is a "negative excess"
			double Smedian= sources[k]->Median;
			if(Smedian<imgMedian) sourcesToBeRemoved.push_back(k);
		}
		CodeUtils::DeleteItems(sources, sourcesToBeRemoved);
	}//close if

	return 0;
	
}//close FindExtendedSources_CV()

//==================================================
//==       FILTERING METHODS
//==================================================

Image* Image::GetMask(Image* mask,bool isBinary)
{

	//## Check input mask
	if(!mask) {
		ERROR_LOG("Null ptr to given image mask!");
		return 0;
	}
		
	//## Check mask bins
	long int Nx= mask->GetNx();
	long int Ny= mask->GetNy();
	if(Nx!=this->GetNx() || Ny!=this->GetNy()){
		ERROR_LOG("Mask binning is different than current image!");
		return 0;
	}	

	//## Clone map
	TString imgName= Form("%s_Mask",m_name.c_str());	
	Image* maskedImage= this->GetCloned(std::string(imgName),true,true);
	maskedImage->Reset();

	//## Loop over mask	
	if(isBinary){
		#ifdef OPENMP_ENABLED
		#pragma omp parallel for collapse(2)
		#endif
		for(long int i=0;i<Nx;i++){	
			for(long int j=0;j<Ny;j++){
				double binContent= this->GetPixelValue(i,j);
				if(binContent==0) continue;
				double maskContent= mask->GetPixelValue(i,j);
				if(maskContent!=0) maskedImage->SetPixelValue(i,j,1);
			}//end loop bins Y
		}//end loop bins X
	}//close if
	else{
		#ifdef OPENMP_ENABLED
		#pragma omp parallel for collapse(2)
		#endif
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				double binContent= this->GetPixelValue(i,j);	
				if(binContent==0) continue;
				double maskContent= mask->GetPixelValue(i,j);
				if(maskContent!=0) maskedImage->SetPixelValue(i,j,binContent);
			}//end loop bins Y
		}//end loop bins X
	}//close else

	//Force re-computation of stats (in parallel computation moments are wrong)
	maskedImage->ComputeStats(true,false,true);

	return maskedImage;

}//close GetMask()

Image* Image::GetSourceMask(std::vector<Source*>const& sources,bool isBinary,bool invert){

	//## Clone map
	long int Nx= this->GetNx();
	long int Ny= this->GetNy();
	bool copyMetaData= true;
	bool resetStats= true;
	TString imgName= Form("%s_SourceMask",m_name.c_str());	
	Image* maskedImage= this->GetCloned(std::string(imgName),copyMetaData,resetStats);
	
	//## Check source list
	int nSources= static_cast<int>(sources.size());
	if(nSources<=0) {
		WARN_LOG("Source list is empty, returning same image!");
		return maskedImage;	
	}

	if(invert){
		if(isBinary){
			#ifdef OPENMP_ENABLED
			#pragma omp parallel for
			#endif
			for(int k=0;k<nSources;k++){
				for(int l=0;l<sources[k]->GetNPixels();l++){
					long int id= (sources[k]->GetPixel(l))->id;
					maskedImage->SetPixelValue(id,0);
				}//end loop pixels
			}//end loop sources	

			#ifdef OPENMP_ENABLED
			#pragma omp parallel for
			#endif
			for(long int i=0;i<Nx;i++){
				for(long int j=0;j<Ny;j++){
					double w= maskedImage->GetPixelValue(i,j);		
					if(w==0) continue;
					maskedImage->SetPixelValue(i,j,1);				
				}
			}
		}//close if
		else{
			#ifdef OPENMP_ENABLED
			#pragma omp parallel for
			#endif
			for(int k=0;k<nSources;k++){
				for(int l=0;l<sources[k]->GetNPixels();l++){
					long int id= (sources[k]->GetPixel(l))->id;
					maskedImage->SetPixelValue(id,0);
				}//end loop pixels
			}//end loop sources		
		}//close else

		//Force re-computation of stats
		maskedImage->ComputeStats(true,false,true);

	}//close if invert
	else{
		//Reset map and loop over sources
		maskedImage->Reset();

		if(isBinary){	
			#ifdef OPENMP_ENABLED
			#pragma omp parallel for
			#endif		
			for(int k=0;k<nSources;k++){
				for(int l=0;l<sources[k]->GetNPixels();l++){
					long int id= (sources[k]->GetPixel(l))->id;
					//maskedImage->FillPixel(id,1);
					maskedImage->SetPixelValue(id,1);
				}//end loop pixels
			}//end loop sources		
		}//close if
		else{
			#ifdef OPENMP_ENABLED
			#pragma omp parallel for
			#endif
			for(int k=0;k<nSources;k++){
				for(int l=0;l<sources[k]->GetNPixels();l++){
					long int id= (sources[k]->GetPixel(l))->id;
					double w= this->GetPixelValue(id);
					//maskedImage->FillPixel(id,w);
					maskedImage->SetPixelValue(id,w);
				}//end loop pixels
			}//end loop sources		
		}//close else

		//Force re-computation of stats
		maskedImage->ComputeStats(true,false,true);

	}//close else

	return maskedImage;

}//close GetSourceMask()

Image* Image::GetSourceResidual(std::vector<Source*>const& sources,int KernSize,int dilateModel,int dilateSourceType,bool skipToNested,ImgBkgData* bkgData,bool useLocalBkg,bool randomize,double zThr){

	//Clone input image	
	TString imgName= Form("%s_Residual",m_name.c_str());	
	Image* residualImg= this->GetCloned(std::string(imgName),true,true);
	if(!residualImg){
		ERROR_LOG("Failed to clone input image, returning nullptr!");
		return nullptr;
	}

	//Dilate source pixels
	int status= MorphFilter::DilateAroundSources(residualImg,sources,KernSize,dilateModel,dilateSourceType,skipToNested,bkgData,useLocalBkg,randomize,zThr);

	if(status<0){
		ERROR_LOG("Failed to dilate sources!");
		if(residualImg) {
			delete residualImg;
			residualImg= 0;
		}
		return 0;		
	}

	//Re-Compute stats
	residualImg->ComputeStats(true,false,true);

	return residualImg;

}//close GetSourceResidual()

Image* Image::GetNormalizedImage(std::string normScale,int normmin,int normmax,bool skipEmptyBins){

	//Get image content min/max
	double wmin= (this->m_StatMoments).minVal;
	double wmax= (this->m_StatMoments).maxVal;
	
	//Create normalized image
	TString imgName= Form("%s_Normalized",m_name.c_str());	
	Image* norm_img= this->GetCloned(std::string(imgName),true,true);
	norm_img->Reset();
	
	//Fill norm image
	if(normScale=="LINEAR"){

		#ifdef OPENMP_ENABLED
		#pragma omp parallel for
		#endif
		for(size_t i=0;i<m_pixels.size();i++){
			double w= m_pixels[i];
			if(skipEmptyBins && w==0) continue;
			double w_norm= normmin + (normmax-normmin)*(w-wmin)/(wmax-wmin);
			norm_img->SetPixelValue(i,w_norm);
		}
	}//close if

	else if(normScale=="LOG"){
		double safemin= 1;
		double safemax= 256;

		#ifdef OPENMP_ENABLED
		#pragma omp parallel for
		#endif
		for(size_t i=0;i<m_pixels.size();i++){
			double w= m_pixels[i];
			if(skipEmptyBins && w==0) continue;
			double w_norm= safemin + (safemax-safemin)*(w-wmin)/(wmax-wmin);
			double w_log= normmin + log10(w_norm/safemin)/log10(safemax/safemin) * (normmax-normmin);
			norm_img->SetPixelValue(i,w_log);
		}

	}//close else if
	else{
		WARN_LOG("Invalid norm scale option selected ("<<normScale<<") no transform applied to original image!");
		return norm_img;
	}

	//Force recomputation of stats if present, otherwise recompute only moments
	bool skipNegativePixels= false;
	bool computeRobustStats= true;	
	bool forceRecomputing= true;
	int status= 0;
	if(this->HasStats()) status= norm_img->ComputeStats(computeRobustStats,skipNegativePixels,forceRecomputing);
	else status= norm_img->ComputeMoments(skipNegativePixels);
	if(status<0){
		WARN_LOG("Failed to re-compute moments/stats for normalized image!");
	}
	
	return norm_img;

}//close GetNormalizedImage()

Image* Image::GetGradientImage(bool invert)
{
	//Compute gradient filtered image
	Image* gradImg= GradientFilter::GetGradientFilter(this);
	if(invert) gradImg->Scale(-1);

	return gradImg;

}//close GetGradientImage()

Image* Image::GetKirschImage()
{
	
	//Compute Kirsh filtered image
	Image* kirschImg= KirschFilter::GetKirschFilter(this);

	return kirschImg;

}//close GetKirschImage()

Image* Image::GetLoGImage(bool invert)
{
	//Compute LoG filtered image
	Image* LoGImg= LoGFilter::GetLoGFilter(this);
	if(LoGImg && invert) LoGImg->Scale(-1);
	
	return LoGImg;

}//close GetLoGImage()

Image* Image::GetNormLoGImage(int size,double scale,bool invert)
{
	//Compute normalized LoG filtered image
	Image* normLoGImg= LoGFilter::GetNormLoGFilter(this,size,scale);
	if(normLoGImg && invert) normLoGImg->Scale(-1);

	return normLoGImg;

}//close GetNormLoGImage()

Image* Image::GetLaplacianImage(bool invert)
{

	//Compute laplacian filtered image
	Image* laplImg= GradientFilter::GetLaplaceFilter(this);
	if(invert) laplImg->Scale(-1);
	
	return laplImg;

}//close GetLaplacianImage()


Image* Image::GetGuidedFilterImage(int radius,double eps)
{

	//Normalize img
	Image* img_norm= this->GetNormalizedImage("LINEAR",1,256);
	img_norm->SetName("tmpImg");
	if(!img_norm) {
		ERROR_LOG("Failed to get normalized image!");
		return 0;
	}

	//## Convert image to OpenCV mat
	cv::Mat I= img_norm->GetOpenCVMat("64");
	cv::Mat p= I;
	eps *= 255*255;   // Because the intensity range of our images is [0, 255]

	//## Run guided filter
	cv::Mat dst = Caesar::guidedFilter(I, p, radius, eps);

	//## Fill filtered image
	TString imgName= Form("%s_GuidedFilter",m_name.c_str());
	Image* FilterImg= this->GetCloned(std::string(imgName),true,true);
	FilterImg->Reset();
	FilterImg->FillFromMat(dst);

	//## Clear allocated data
	if(img_norm) {
		delete img_norm;
		img_norm= 0;
	}

	return FilterImg;

}//close GetGuidedFilterImage()

Image* Image::GetSmoothedImage(int size_x,int size_y,double sigma_x,double sigma_y)
{

	//## Get OpenCV mat
	cv::Mat mat= this->GetOpenCVMat("64");
	
	//## Smooth matrix
	cv::Size smooth_size(size_x,size_y);
	cv::Mat smoothed_mat;
	cv::GaussianBlur(mat,smoothed_mat, smooth_size, sigma_x, sigma_y, cv::BORDER_DEFAULT);

	//## Fill smoothed image
	TString imgName= Form("%s_Smoothed",m_name.c_str());
	Image* SmoothedImg= this->GetCloned(std::string(imgName),true,true);
	SmoothedImg->Reset();
	SmoothedImg->FillFromMat(smoothed_mat);
	
	return SmoothedImg;

}//close GetSmoothedImage()

std::vector<Image*> Image::GetWaveletDecomposition(int nScales){

	DEBUG_LOG("Computing wavelet decomposition up to scale J="<<nScales<<" ...");
	std::vector<Image*> img_decomposition;
	img_decomposition= WTFilter::GetDecomposition(this,nScales);
	
	return img_decomposition;

}//close GetWaveletDecomposition()


Image* Image::GetSaliencyMap(int reso,double regFactor,int minRegionSize,double knnFactor,bool useRobust,double expFalloffPar,double distanceRegPar){

	//## Compute single-reso saliency map
	Image* saliencyMap= 0;
	saliencyMap= SaliencyFilter::ComputeSaliencyMap(this,reso,regFactor,minRegionSize,knnFactor,useRobust,expFalloffPar,distanceRegPar);
	if(!saliencyMap){
		ERROR_LOG("Saliency map estimation failed!");
		return nullptr;
	}

	return saliencyMap;

}//close GetSaliencyMap()


Image* Image::GetMultiResoSaliencyMap(int resoMin,int resoMax,int resoStep,double beta,int minRegionSize,double knnFactor,bool useRobustPars,double expFalloffPar,double distanceRegPar,double salientMultiplicityThrFactor,bool addBkgMap,bool addNoiseMap,ImgBkgData* bkgData,double saliencyThrFactor,double imgThrFactor){

	//## Compute multi-reso saliency map
	Image* saliencyMap= 0;
	saliencyMap= SaliencyFilter::ComputeMultiResoSaliencyMap(this,resoMin,resoMax,resoStep,beta,minRegionSize,knnFactor,useRobustPars,expFalloffPar,distanceRegPar, salientMultiplicityThrFactor,addBkgMap,addNoiseMap,bkgData,saliencyThrFactor,imgThrFactor);
	if(!saliencyMap){
		ERROR_LOG("Multi-resolution saliency map estimation failed!");
		return nullptr;
	}

	return saliencyMap;

}//close GetMultiResoSaliencyMap()


int Image::Add(Image* img,double c,bool computeStats)
{

	//Check input image
	if(!img){
		ERROR_LOG("Null ptr to given input image!");
		return -1;
	}
	long int Nx= img->GetNx();
	long int Ny= img->GetNy();
	long int PixDataSize= img->GetPixelDataSize();
	if(m_Nx!=Nx || m_Ny!=Ny){
		ERROR_LOG("Image to be added has different size ("<<Nx<<","<<Ny<<") wrt to this image ("<<m_Nx<<","<<m_Ny<<")!");
		return -1;
	}
	if(PixDataSize!=this->GetPixelDataSize()){
		ERROR_LOG("Image to be added has different pixel vector size ("<<PixDataSize<<") wrt to this image ("<<m_pixels.size()<<")!");
		return -1;
	}

	//Reset stats
	ResetImgStats(true,true);
	
	//Loop to sum vectors
	bool useNegativePixInStats= true;
	
	#ifdef OPENMP_ENABLED		
		Caesar::StatMoments<double> moments_t;	
		std::vector<Caesar::StatMoments<double>> parallel_moments;
	
		//#pragma omp declare reduction (merge : std::vector<Caesar::StatMoments<double>> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
		#pragma omp parallel private(moments_t) //reduction(merge: parallel_moments)
		{
			int thread_id= omp_get_thread_num();
			int nthreads= SysUtils::GetOMPThreads();
			DEBUG_LOG("Starting multithread image add (thread_id="<<thread_id<<", nthreads="<<nthreads<<")");

			#pragma omp single
   		{
     		parallel_moments.assign(nthreads,Caesar::StatMoments<double>());
   		}

			#pragma omp for 
			for(size_t i=0;i<m_pixels.size();i++){
				double w1= m_pixels[i];			
				double w2= img->GetPixelValue(i);	
				double w= w1 + c*w2;
				m_pixels[i]= 0;
				if(FillPixelMT(moments_t,i,w,useNegativePixInStats)<0) continue;
			}
			
			//parallel_moments.push_back(moments_t);
			parallel_moments[thread_id]= moments_t;
		}//close parallel section
		
		//Update moments from parallel estimates
		if(Caesar::StatsUtils::ComputeMomentsFromParallel(m_StatMoments,parallel_moments)<0){
			ERROR_LOG("Failed to compute cumulative moments from parallel estimates (NB: image will have wrong moments!)");
			return -1;
		}

	#else 
		for(size_t i=0;i<m_pixels.size();i++){
			double w1= m_pixels[i];			
			double w2= img->GetPixelValue(i);	
			double w= w1 + c*w2;
			m_pixels[i]= 0;
			if(FillPixel(i,w,useNegativePixInStats)<0) continue;
		}
	#endif
	

	if(computeStats){
		bool computeRobustStats= true;
		bool skipNegativePixels= false;
		bool forceRecomputing= false;
		if(ComputeStats(computeRobustStats,skipNegativePixels,forceRecomputing)<0){
			WARN_LOG("Failed to compute stats after adding the two images!");
			return -1;
		}	
	}//close if computeStats

	return 0;

}//close Add()

int Image::Scale(double c)
{
	//Multiply pixel values by a factor c
	#ifdef OPENMP_ENABLED
	#pragma omp parallel for
	#endif
	for(size_t i=0;i<m_pixels.size();i++){
		double w_old= m_pixels[i];
		double w= w_old*c;
		m_pixels[i]= w;
	}

	//Force recomputation of stats if present, otherwise recompute only moments
	bool skipNegativePixels= false;
	bool computeRobustStats= true;	
	bool forceRecomputing= true;
	int status= 0;
	if(this->HasStats()) status= this->ComputeStats(computeRobustStats,skipNegativePixels,forceRecomputing);
	else status= this->ComputeMoments(skipNegativePixels);
		
	return status;

}//close Scale()

//================================
//==    THRESHOLDING METHODS
//================================
TH1D* Image::GetPixelHisto(int nbins,bool normalize){

	//Check if image has stats computed
	if(!HasStats()){
		WARN_LOG("No stats computed, returning nullptr!");
		return nullptr;
	}

	double Smin= m_Stats->min;
	double Smax= m_Stats->max;
	double Srange= Smax-Smin;
	double tol= 0.0;
	double Smin_tol= Smin-tol*fabs(Srange);
	double Smax_tol= Smax+tol*fabs(Srange);

	TString histoName= Form("%s_histo",m_name.c_str());
	TH1D* histo= new TH1D(histoName,histoName,nbins,Smin_tol,Smax_tol);
		
	for(size_t i=0;i<m_pixels.size();i++){
		double w= m_pixels[i];
		if(w==0) continue;
		histo->Fill(w);
	}

	if(normalize) histo->Scale(1./histo->Integral());
	return histo;

}//close GetPixelHisto()


double Image::FindOtsuThreshold(int nbins){
	
	//## Get histo and normalize
	TH1D* hist= this->GetPixelHisto(nbins,true);
	if(!hist) {
		ERROR_LOG("Failed to compute pixel histo, return thr=0!");
		return 0;
	}

	int Nx= hist->GetNbinsX();
	double sum = 0;
	double Wxs[Nx];
	for (int t=0;t<Nx;t++) {
		double w= hist->GetBinContent(t+1);
		double x= hist->GetBinCenter(t+1);
		Wxs[t]= w*x;
		sum+= w*x;
	}

	double sumB = 0;
	double wB = 0;
	double wF = 0;

	double varMax = 0;
	double threshold = 0;
	
	for (int t=0;t<Nx;t++) {
		double binX= hist->GetBinCenter(t+1);
		double binContent= hist->GetBinContent(t+1);	
		double Wx= Wxs[t];
  	wB += binContent;               // Weight Background
   	if (wB == 0) continue;

   	wF = 1. - wB;                 // Weight Foreground
   	if (wF == 0) break;

		sumB+= Wx;
		double mB = sumB/wB;            // Mean Background
   	double mF = (sum - sumB)/wF;    // Mean Foreground

   	// Calculate Between Class Variance
   	double varBetween = wB*wF * (mB - mF) * (mB - mF);
		
   	// Check if new maximum found
   	if (varBetween > varMax) {
   		varMax = varBetween;
      threshold = binX;
   	}
	}//end loop bins

	if(hist) {
		hist->Delete();
		hist= 0;
	}

	return threshold;

}//close FindOtsuThreshold()


double Image::FindValleyThreshold(int nbins,bool smooth){

	//## Init dilate kernel size and histos
	const int nKernels= 3;
	int kernelSizes[]= {3,5,7};//{5,7,9}
	int maxStep= floor(kernelSizes[nKernels-1]/2.);

	//## Get pixel histo (invert to find peaks corresponding to valley in original histo)
	TH1D* histo= this->GetPixelHisto(nbins);
	if(smooth) histo->Smooth(1);
	histo->Scale(-1);
	double sMin= histo->GetMinimum();
	const double SMALL_NUMBER= 1.e-6;	
	for(int i=0;i<histo->GetNbinsX();i++){
		double w= histo->GetBinContent(i+1);
		double wnew= w-sMin+SMALL_NUMBER;
		histo->SetBinContent(i+1,wnew);
	}//end loop bins

	//## Get dilated histos
	TH1D* dilatedHisto= 0;
	std::vector<TH1D*> dilatedHistoList;

	for(int k=0;k<nKernels;k++){//loop over kernels
		TString histoName= Form("hdilate_%d",k+1);
		dilatedHisto= (TH1D*)histo->Clone(histoName);
		dilatedHisto->Reset();
		dilatedHistoList.push_back(dilatedHisto);
		
		int step= floor(kernelSizes[k]/2.);		

		for(int i=0;i<histo->GetNbinsX();i++){//loop over bins	
			double wmax= -1.e+99;
			for(int j=-step;j<step;j++){//Loop over kernel range
				int s= i+j;
				int binId= s+1;
				if(histo->IsBinUnderflow(binId) || histo->IsBinOverflow(binId)) continue;
				double w= histo->GetBinContent(binId);
				if(w>wmax) wmax= w;
			}//end loop kernel

			dilatedHistoList[k]->SetBinContent(i+1,wmax);
		}//end loop bins
	}//end loop kernels

	//## Find valleys
	TGraph* peaks= new TGraph;
	int npeaks= 0;
	double valleyThr= 0;
	for(int i=0;i<histo->GetNbinsX();i++){
		int binId= i+1;
		if(binId<maxStep+1) continue;
		double x= histo->GetXaxis()->GetBinCenter(i+1);
		double w= histo->GetBinContent(i+1);
		bool isPeak= true;
		for(int k=0;k<nKernels;k++) {	
			double wdilate= dilatedHistoList[k]->GetBinContent(binId); 
			if(wdilate!=w){
				isPeak= false;	
				break;
			}
		}//end loop kernels

		if(isPeak){
			peaks->SetPoint(npeaks,x,w);
			if(npeaks==0){
				valleyThr= x;
			}
			npeaks++;			
		}
	}//end loop bins

	DEBUG_LOG("#"<<npeaks<<" valleys detected!");
	
	//## Clear stuff
	if(histo){
		histo->Delete();
		histo= 0;
	}
	if(peaks) peaks->Delete();
	for(int k=0;k<nKernels;k++) {
		if(dilatedHistoList[k]) {
			dilatedHistoList[k]->Delete();
			dilatedHistoList[k]= 0;
		}
	}//end loop kernels
	dilatedHistoList.clear();

	return valleyThr;

}//close FindValleyThreshold()

Image* Image::GetBinarizedImage(double threshold,double fgValue,bool isLowerThreshold){

	TString imgName= Form("%s_Binarized",m_name.c_str());	
	Image* BinarizedImg= this->GetCloned(std::string(imgName),true,true);
	BinarizedImg->Reset();

	#ifdef OPENMP_ENABLED
	#pragma omp parallel for
	#endif
	for(size_t i=0;i<m_pixels.size();i++){
		double w= m_pixels[i];
		if(w==0) continue;
		if(w>=threshold && !isLowerThreshold) {
			BinarizedImg->SetPixelValue(i,fgValue);
		}
		else if(w<threshold && isLowerThreshold) {
			BinarizedImg->SetPixelValue(i,fgValue);
		}
	}//end loop pixels
	
	//Force recomputation of stats if present, otherwise recompute only moments
	bool skipNegativePixels= false;
	if(BinarizedImg->ComputeMoments(skipNegativePixels)<0){
		ERROR_LOG("Failed to re-compute moments of binarized image!");
		return nullptr;
	}	
	
	return BinarizedImg;

}//close GetBinarizedImage()

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

TMatrixD* Image::GetMatrix(){

	//Allocate matrix
	TMatrixD* M= new TMatrixD(m_Ny,m_Nx);
	M->Zero();

	//Fill matrix
	#ifdef OPENMP_ENABLED
	#pragma omp parallel for collapse(2)
	#endif
	for(long int i=0;i<m_Nx;i++){//rows
		for(long int j=0;j<m_Ny;j++){//columns
			double w= this->GetPixelValue(i,j);
			long int rowId= j;
			long int colId= i;
			(*M)(rowId,colId)= w;
		}//end loop cols
	}//end loop rows

	return M;

}//close Img::GetMatrix()

cv::Mat Image::GetOpenCVMat(std::string encoding){

	long int Nx= this->GetNx();
	long int Ny= this->GetNy();

	//## Fill OpenCV mat
	cv::Mat mat;
	if(encoding=="64") mat= cv::Mat::zeros(Ny,Nx,CV_64FC1);
	else if(encoding=="32") mat= cv::Mat::zeros(Ny,Nx,CV_32FC1);
	else{
		WARN_LOG("Invalid encoding selected, using default 64bit encoding");
		mat= cv::Mat::zeros(Ny,Nx,CV_64FC1);
	}

	//The fast way
	long int nRows = mat.rows;
  long int nCols = mat.cols;
	
	#ifdef OPENMP_ENABLED
	#pragma omp parallel for
	#endif
	for(long int i=0;i<nRows;++i) {
		long int rowId= i;
		long int iy= Ny-1-rowId;
  	double* p = mat.ptr<double>(i);
    for (long int j=0;j<nCols;++j){
			int colId= j;
			int ix= colId;
			double w= this->GetPixelValue(ix,iy);
    	p[j] = w;
    }
  }

	return mat;

}//close ImgToMat()

//================================================================
//===    DRAW  METHODS
//================================================================
int Image::Draw(int palette,bool drawFull,bool useCurrentCanvas,std::string units)
{
	//Set palette
	Caesar::GraphicsUtils::SetPalette(palette);

	//Get temp histogram
	TH2D* htemp= GetHisto2D("htemp");
	if(!htemp){
		ERROR_LOG("Failed to get histo from this image!");
		return -1;
	}
	
	//Set canvas
	int canvas_width= 720;
	int canvas_height= 700;
	TString canvasName= Form("%s_Plot",m_name.c_str());	
	TCanvas* canvas= 0;
	if(useCurrentCanvas && gPad) {
		canvas= gPad->GetCanvas();
		canvas->SetName(canvasName);
		canvas->SetTitle(canvasName);
	}
	else{
		canvas= new TCanvas(canvasName,canvasName,canvas_width,canvas_height);
	}

	if(!canvas){
		ERROR_LOG("Failed to retrieve or set canvas!");
		return -1;
	}

	//Draw full image (without borders)
	canvas->cd();
	htemp->SetStats(0);
	
	if(drawFull){
		canvas->ToggleEventStatus();
  	canvas->SetRightMargin(0.0);
  	canvas->SetLeftMargin(0.0);
  	canvas->SetTopMargin(0.0);
  	canvas->SetBottomMargin(0.0);
		htemp->Draw("COLA");
	}
	else{
		gStyle->SetPadTopMargin(0.1);
  	gStyle->SetPadBottomMargin(0.1);
  	gStyle->SetPadLeftMargin(0.15);
  	gStyle->SetPadRightMargin(0.15);

		gPad->SetTopMargin(0.1);
		gPad->SetBottomMargin(0.1);
		gPad->SetLeftMargin(0.15);
  	gPad->SetRightMargin(0.15);

		//Set palette axis title	
		htemp->GetZaxis()->SetTitle(units.c_str());
		htemp->GetZaxis()->SetTitleSize(0.05);
		htemp->GetZaxis()->SetTitleOffset(0.9);
		htemp->Draw();
		gPad->Update();
		htemp->Draw("COLZ");
		gPad->Update();
  }

	return 0;

}//close Draw()


int Image::Draw(std::vector<Source*>const& sources,int palette,bool drawFull,bool useCurrentCanvas,std::string units)
{
	
	//Draw image first
	Draw(palette,drawFull,useCurrentCanvas,units);

	//Retrieve canvas and draw sources	
	TCanvas* canvas= gPad->GetCanvas();
	if(!canvas){
		WARN_LOG("Failed to get access to current canvas!");
		return -1;
	}

	//Draw sources in current canvas
	for(unsigned int k=0;k<sources.size();k++){	
		int type= sources[k]->Type;
		int lineColor= kBlack;
		if(type==Source::eCompact)
			lineColor= kBlack;	
		else if(type==Source::ePointLike)
			lineColor= kRed;
		else if(type==Source::eExtended)
			lineColor= kGreen+1;	
		sources[k]->Draw(false,false,true,lineColor);
	}//end loop sources

	return 0;

}//close Draw()


}//close namespace