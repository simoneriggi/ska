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
* @file Img.cc
* @class Img
* @brief Img
*
* Image class
* @author S. Riggi
* @date 20/01/2015
*/

#include <Img.h>


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


ClassImp(Caesar::ImgMetaData)
ClassImp(Caesar::ImgStats)
ClassImp(Caesar::Img)

namespace Caesar {

//Default constructor
Img::Img() : TH2F()
{
	//Do not allocate memory in default constructor otherwise you will have a memory leak 
	//(see https://root.cern.ch/root/html534/guides/users-guide/AddingaClass.html#the-default-constructor)
  Init();
}//close costructor

Img::Img(const char *name,const char *title,int nbinsx,float xlow,float xup,int nbinsy,float ylow,float yup)
	:TH2F(name,title,nbinsx,xlow,xup,nbinsy,ylow,yup)
{
  Init();
}//close constructor

Img::Img(const char *name,const char *title,int nbinsx,const float *xbins,int nbinsy,const float *ybins)
	:TH2F(name,title,nbinsx,xbins,nbinsy,ybins)
{
  Init();
}//close constructor


Img::Img(const Img &img) : TH2F() {
	// Copy constructor
	DEBUG_LOG("Copy constuctor called...");
  ((Img&)img).Copy(*this);
}//close constructor


void Img::Copy(TObject& obj) const
{
	//cout<<"Img::Copy(): INFO: Copying image "<<((Img&)obj).GetName()<<"..."<<endl;
	TH2F::Copy((Img&)obj);
	//cout<<"done!"<<endl;

	//cout<<"Img::Copy(): INFO: Copying private vars..."<<endl;
	((Img&)obj).m_HasMetaData= m_HasMetaData;
	((Img&)obj).m_HasStats= m_HasStats;
	((Img&)obj).m_Npix= m_Npix;
	((Img&)obj).m_M1= m_M1;
	((Img&)obj).m_M2= m_M2;
	((Img&)obj).m_M3= m_M3;
	((Img&)obj).m_M4= m_M4;
	((Img&)obj).m_PixelMin= m_PixelMin;
	((Img&)obj).m_PixelMax= m_PixelMax;
	//cout<<"done!"<<endl;

	//cout<<"Img::Copy(): INFO: Copying meta data..."<<endl;
	if(m_MetaData){
		((Img&)obj).m_MetaData= new ImgMetaData;
		*((Img&)obj).m_MetaData = *m_MetaData;
	}
	//cout<<"done!"<<endl;
	//cout<<"Img::Copy(): INFO: Copying stats data..."<<endl;
	if(m_Stats){
		((Img&)obj).m_Stats= new ImgStats;
		*((Img&)obj).m_Stats = *m_Stats;
	}
	//cout<<"done!"<<endl;
	

	
  
}//close Copy()


Img::~Img(){

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


Img& Img::operator=(const Img &img) { 
	// Operator =
  if (this != &img)  ((Img&)img).Copy(*this);
  return *this;
}


void Img::Init(){
	
	//Meta-Data
	m_HasMetaData= false;
	m_MetaData= 0;
 
	//Stats
	m_HasStats= false;
	m_Stats= 0;
	ResetImgStats(true);//Reset stats & moments

}//close Init()



int Img::FillPixel(double x,double y,double w,bool useNegativePixInStats){

	//Check pixel value
	if( TMath::IsNaN(w) || fabs(w)==TMath::Infinity() ) return -1;

	//Check bin
	int binx, biny, bin;
	binx = fXaxis.FindBin(x);
  biny = fYaxis.FindBin(y);
	if (binx<0 || biny <0) return -1;
  bin  = biny*(fXaxis.GetNbins()+2) + binx;

	int ix= binx-1;
	int iy= biny-1;	
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();

	//Update pixel stats
	if( (ix>=0 && ix<Nx) && (iy>=0 || iy<Ny) &&
			( w>=0 || (w<0 && useNegativePixInStats) ) 
		) 
	{//update moments

		double w_old= GetBinContent(bin);

		if( w_old!=0 ){//bin already filled
			WARN_LOG("Pixel "<<bin<<" ("<<ix<<","<<iy<<") has been already filled (w="<<w_old<<"), skipping...");
			return -1;
		}//close if
		else{
			UpdateMoments(w);
		}
	}//close if
	
	//Copied from TH2::Fill(x,y,w)
	int filled_bin= Fill(x,y,w);
	
  return filled_bin;

}//close Img::FillPixel()

void Img::FillFromMat(cv::Mat mat,bool useNegativePixInStats){

	int nRows = mat.rows;
  int nCols = mat.cols;
	int Ny= this->GetNbinsY();
	
	for(int i=0;i<nRows;++i) {
		int rowId= i;
		int iy= Ny-1-rowId;
		double y= this->GetYaxis()->GetBinCenter(iy+1);
  	double* p = mat.ptr<double>(i);
    for (int j=0;j<nCols;++j){
			int colId= j;
			int ix= colId;
			double x= this->GetXaxis()->GetBinCenter(ix+1);
			double w= p[j];
			this->FillPixel(x,y,w,useNegativePixInStats);
    }//end loop cols
  }//end loop rows
	
}//close FillFromMat()

void Img::ResetImgStats(bool resetMoments,bool clearStats){

	if(m_Stats){
		m_Stats->Reset();
	}//close if

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

}//close Img::ResetImgStats()


void Img::ComputeMoments(bool skipNegativePixels){

	#ifdef OPENMP_ENABLED
		
		//Define variables for reduction (OMP does not allow class members in reduction operation)
		long long int nPix= m_Npix;
		double pixelMin= m_PixelMin; 
		double pixelMax= m_PixelMax;
		double M1= m_M1;
		double M2= m_M2;
		double M3= m_M3;
		double M4= m_M4;

		#pragma omp parallel
		{
			INFO_LOG("Starting image moment computing in thread "<<omp_get_thread_num()<<" (nthreads="<<SysUtils::GetOMPThreads()<<") ...");
			
			#pragma omp for reduction(max: pixelMax), reduction(min: pixelMin), reduction(+: nPix,M1,M2,M3,M4)
		
			for(int i=0;i<this->GetNbinsX();i++){
				for(int j=0;j<this->GetNbinsY();j++){
					double w= this->GetBinContent(i+1,j+1);
					if( w==0 || (skipNegativePixels && w<0) ) continue; 
					if(w<pixelMin) pixelMin= w;
					if(w>pixelMax) pixelMax= w;
					nPix++;
  				double delta = w - M1;
  				double delta_n = delta/nPix;
  				double delta_n2 = delta_n * delta_n;
  				double f = delta * delta_n * (nPix-1);
  				M1+= delta_n;
  				M4+= f * delta_n2 * (nPix*nPix - 3*nPix + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3;
  				M3+= f * delta_n * (nPix - 2) - 3 * delta_n * M2;
  				M2+= f;
				}
			}
		}//close parallel section

		//Set reduced vars to global class var
		m_Npix = nPix;
		m_PixelMin= pixelMin;
		m_PixelMax= pixelMax;
		m_M1= M1;
		m_M2= M2;
		m_M3= M3;
		m_M4= M4;

	#else
		for(int i=0;i<this->GetNbinsX();i++){
			for(int j=0;j<this->GetNbinsY();j++){
				double w= this->GetBinContent(i+1,j+1);
				if( w==0 || (skipNegativePixels && w<0) ) continue; 
				UpdateMoments(w);
			}
		}
	#endif

}//close ComputeMoments()


void Img::UpdateMoments(double w){

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
	
}//close Img::UpdateMoments()

void Img::ComputeStatsParams(bool computeRobustStats,bool skipNegativePixels){

	//## Reset previous stats params (Reset only Stats not moments!)
	ResetImgStats(false);
	
	//## Compute Stats params
	m_Stats->n= m_Npix;
	m_Stats->min= m_PixelMin;
	m_Stats->max= m_PixelMax;
	m_Stats->mean= m_M1;
	m_Stats->rms= 0;
	if(m_Npix>2) m_Stats->rms= sqrt(m_M2/(m_Npix-1));
	m_Stats->skewness= 0;
	if(m_M2!=0) {
  	m_Stats->skewness= sqrt(m_Npix)*m_M3/pow(m_M2,1.5);//need to adjust for finite population?
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
	if(m_Npix>2) m_Stats->skewnessErr= sqrt(6.*m_Npix*(m_Npix-1)/((m_Npix-2)*(m_Npix+1)*(m_Npix+3)));//approximate for normal distribution
	
  
	//## End if no robust stats are to be computed
	if(!computeRobustStats) return;
	
	//## Compute robust stats (median, MAD, ...)	
	std::vector<double> pixels;
	pixels.clear();
	pixels.resize(0);
	for(int i=0;i<this->GetNbinsX();i++){
		for(int j=0;j<this->GetNbinsY();j++){
			double pixelValue= this->GetBinContent(i+1,j+1);
			if( pixelValue==0 || (skipNegativePixels && pixelValue<0) ) continue;
			pixels.push_back(pixelValue);
		}
	}
	
	//Sort and compute median for all image	
	double median= Caesar::StatsUtils::GetMedianFast<double>(pixels);
	m_Stats->median= median;

	//Compute MAD = median(|x_i-median|)
	double medianMAD= Caesar::StatsUtils::GetMADFast(pixels,median);
	double medianRMS= medianMAD*1.4826;//0.6744888;
	m_Stats->medianRMS= medianRMS;

	//## Compute biweight robust estimators
	double C= 6.;
	double tol= 0.0001;
	double nmaxIter= 10;
	std::pair<double,double> biweightEstimators= Caesar::StatsUtils::GetBiWeightEstimators<double>(pixels,median,medianRMS,C,tol,nmaxIter);
	m_Stats->bwLocation= biweightEstimators.first;
	m_Stats->bwScale= biweightEstimators.second;

	//## Compute clipped estimators
	double clipSigma= 3;
	int clipMaxIter= 100;
	double clipTolerance= 0.1;
	bool useParallelVersion= false;
	//std::pair<double,double> clippedEstimators= Caesar::StatsUtils::GetClippedEstimators<double>(pixels,clipSigma,clipMaxIter,clipTolerance, false);
	//m_Stats->clippedMedian= clippedEstimators.first;
	//m_Stats->clippedRMS= clippedEstimators.second;
	
	ClippedStats<double> clipped_stats;
	Caesar::StatsUtils::GetClippedEstimators(clipped_stats,pixels,median,m_Stats->mean,m_Stats->rms,clipSigma,clipMaxIter,clipTolerance,useParallelVersion);
	m_Stats->clippedMedian= clipped_stats.median;
	m_Stats->clippedRMS= clipped_stats.stddev;
	
	
}//close Img::ComputeStatsParams()



int Img::ComputeStats(bool computeRobustStats,bool skipNegativePixels,bool forceRecomputing){

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



Img* Img::GetTile(int ix_min,int ix_max,int iy_min,int iy_max){

	//## Check args
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	
	//## Check tile sizes
	std::string errflag= "";
	if(ix_min<0 || ix_min>=Nx || ix_min>ix_max){
		ERROR_LOG("Invalid min tile X given ("<<ix_min<<")");
		return 0;
	}
	if(ix_max<0 || ix_max>=Nx || ix_max<ix_min){
		ERROR_LOG("Invalid max tile X given ("<<ix_max<<")");
		return 0;
	}
	if(iy_min<0 || iy_min>=Ny || iy_min>iy_max){
		ERROR_LOG("Invalid min tile Y given ("<<iy_min<<")");
		return 0;
	}
	if(iy_max<0 || iy_max>=Ny || iy_max<iy_min){
		ERROR_LOG("Invalid max tile Y given ("<<iy_max<<")");
		return 0;
	}

	int TileSizeX= ix_max-ix_min+1;
	int TileSizeY= iy_max-iy_min+1;

	double binX_min= this->GetXaxis()->GetBinCenter(ix_min+1);
	double binX_max= this->GetXaxis()->GetBinCenter(ix_max+1);
	double binY_min= this->GetYaxis()->GetBinCenter(iy_min+1);
	double binY_max= this->GetYaxis()->GetBinCenter(iy_max+1);
	
	//## Init tile
	Img* tile= new Img;
	TString tileName= Form("%s_x%d_%d_y%d_%d",this->GetName(),ix_min,ix_max,iy_min,iy_max);
	tile->SetNameTitle(tileName,tileName);
	tile->SetBins(TileSizeX,binX_min-0.5,binX_max+0.5,TileSizeY,binY_min-0.5,binY_max+0.5);
	tile->CopyMetaData(this->m_MetaData);
	tile->Reset();
	
	//## Read tile
	//cout<<"Img::GetTile(): INFO: Read tile ix("<<ix_min<<","<<ix_max<<") iy("<<iy_min<<","<<iy_max<<") binX("<<binX_min<<","<<binX_max<<")"<<" binY("<<binY_min<<","<<binY_max<<")"<<endl;
	for(int j=iy_min;j<=iy_max;j++){
		double binY= this->GetYaxis()->GetBinCenter(j+1);
		for(int i=ix_min;i<=ix_max;i++){		
			double binX= this->GetXaxis()->GetBinCenter(i+1);
			double pixValue= this->GetBinContent(i+1,j+1);
			if(pixValue==0 || TMath::IsNaN(pixValue) || fabs(pixValue)==TMath::Infinity()) continue; 
			tile->FillPixel(binX,binY,pixValue);	
		}//end loop columns
	}//end loop row
	
	return tile;

}//close Img::GetTile()



BkgData* Img::ComputeBkg(int estimator,bool computeLocalBkg,int boxSizeX,int boxSizeY, double gridStepSizeX,double gridStepSizeY,bool use2ndPass,bool skipOutliers,double seedThr,double mergeThr,int minPixels){
	
	//## Compute bkg data
	DEBUG_LOG("Using grid bkg method...");
	BkgData* bkgData= BkgFinder::FindBkg(this,estimator,computeLocalBkg, boxSizeX, boxSizeY, gridStepSizeX,gridStepSizeY,use2ndPass,skipOutliers,seedThr,mergeThr,minPixels);	
	if(!bkgData){
		ERROR_LOG("Computation of local background failed for this image!");
		return 0;
	}
	return bkgData;
	
}//close Img::ComputeBkg()

Img* Img::GetSignificanceMap(BkgData* bkgData,bool useLocalBkg){

	//Check image
	if(!bkgData){
		ERROR_LOG("Null ptr to bkg data!");
		return 0;
	}
	//Integrity check for local bkg
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	if(useLocalBkg){
		if(!bkgData->HasLocalBkg()){
			ERROR_LOG("Local bkg option requested but no local bkg data are available!");
			return 0;
		}
		if( Nx!=(bkgData->BkgMap)->GetNbinsX() || Ny!=(bkgData->BkgMap)->GetNbinsY() ||
				Nx!=(bkgData->NoiseMap)->GetNbinsX() || Ny!=(bkgData->NoiseMap)->GetNbinsY()
		){
			cerr<<"Img::GetSignificanceMap(): ERROR: Bkg/Noise maps have different size!"<<endl;			
			return 0;
		}
	}//close if
		
	//Clone this image and reset content
	TString imgName= Form("%s_significance",std::string(this->GetName()).c_str());
	Img* significanceMap= this->GetCloned(std::string(imgName),true,true);
	significanceMap->Reset();

	if(useLocalBkg){
		for(int i=0;i<Nx;i++){
			double binX= this->GetXaxis()->GetBinCenter(i+1);
			for(int j=0;j<Ny;j++){	
				double binY= this->GetYaxis()->GetBinCenter(j+1);
				double w= this->GetBinContent(i+1,j+1);	
											
				double bkgLevel= (bkgData->BkgMap)->GetBinContent(i+1,j+1);
				double bkgRMS= (bkgData->NoiseMap)->GetBinContent(i+1,j+1);
				if( w==0 || bkgRMS<=0) {
					//cerr<<"Img::GetSignificanceMap(): WARN: Empty pixel or invalid bkg  ("<<i<<","<<j<<") skip it!"<<endl;
					continue;
				}
				double Z= (w-bkgLevel)/bkgRMS;
				significanceMap->FillPixel(binX,binY,Z);
			}//end loop
		}//end loop 
	}//close if useLocalBkg
	else{
		double bkgLevel= bkgData->gBkg;
		double bkgRMS= bkgData->gNoise;
				
		for(int i=0;i<Nx;i++){
			double binX= this->GetXaxis()->GetBinCenter(i+1);
			for(int j=0;j<Ny;j++){	
				double binY= this->GetYaxis()->GetBinCenter(j+1);
				double w= this->GetBinContent(i+1,j+1);										
				if( w==0 || bkgRMS<=0) {
					continue;
				}
				double Z= (w-bkgLevel)/bkgRMS;
				significanceMap->FillPixel(binX,binY,Z);
			}//end loop
		}//end loop
	}//close else

	return significanceMap;

}//close GetSignificanceMap()


int Img::FindCompactSource(std::vector<Source*>& sources,Img* floodImg,BkgData* bkgData,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool mergeBelowSeed,bool findNestedSources,double nestedBlobThreshold){

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


int Img::FindNestedSource(std::vector<Source*>& sources,BkgData* bkgData,int minPixels,double nestedBlobThreshold){

	//Check if given mother source list is empty
	int nSources= (int)sources.size();
	if(nSources<=0){
		WARN_LOG("Empty source list given!");
		return 0;
	}

	//Find image mask of found sources
	Img* sourceMask= this->GetSourceMask(sources,false);
	if(!sourceMask){
		ERROR_LOG("Null ptr to computed source mask!");
		return -1;
	}

	Img* curvMap= this->GetLaplacianImage(true);
	if(!curvMap){
		ERROR_LOG("Null ptr to computed curvature mask!");
		if(sourceMask) sourceMask->Delete();
		return -1;
	}
	curvMap->ComputeStats(true,false,false);
	double curvMapRMS= curvMap->GetPixelStats()->medianRMS;
	double curvMapThr= curvMapRMS*nestedBlobThreshold;
	Img* blobMask= curvMap->GetBinarizedImage(curvMapThr);
	if(!blobMask){
		ERROR_LOG("Failed to compute blob mask!");
		if(sourceMask) sourceMask->Delete();
		if(curvMap) curvMap->Delete();
		return -1;
	}
	
	/*
	//Find blob mask
	int kernFWHM_safe= 5;
	int kernFWHM= kernFWHM_safe;
	if(this->HasMetaData()){
		kernFWHM= m_MetaData->GetBeamSizeInPixel();
		if(kernFWHM==0){
			cerr<<"Img::FindNestedSource(): WARN: No (or invalid) beam information present in the image, assuming a safe value!"<<endl;
			kernFWHM= kernFWHM_safe;
		}
	}
	double kernSigma= kernFWHM/(2*sqrt(2*log(2)));
	double kernSigmaStep= 1;
	int kernFactor= 6;
	int thrModel= 2;
	Img* blobMask= BlobFinder::GetMultiScaleBlobMask(this,kernFactor,kernSigma,kernSigma,kernSigmaStep,thrModel,nestedBlobThreshold);
	if(!blobMask){
		cerr<<"Img::FindNestedSource(): ERROR: Failed to compute blob mask!"<<endl;
		if(sourceMask) sourceMask->Delete();
		return -1;
	}
	*/

	//Find blob+source mask
	Img* sourcePlusBlobMask= sourceMask->GetMask(blobMask,true);
	if(!sourcePlusBlobMask){
		ERROR_LOG("Failed to compute (source+blob) mask!");
		if(sourceMask) sourceMask->Delete();
		if(curvMap) curvMap->Delete();
		if(blobMask) blobMask->Delete();
		return -1;
	}

	//Find nested blobs 
	std::vector<Source*> NestedSources;
	double fgValue= 1;
	int status= BlobFinder::FindBlobs(this,NestedSources,sourcePlusBlobMask,bkgData,fgValue,fgValue,minPixels,false,false);
	if(status<0){
		ERROR_LOG("Nested blob finder failed!");
		if(sourceMask) sourceMask->Delete();
		if(blobMask) blobMask->Delete();
		if(curvMap) curvMap->Delete();
		if(sourcePlusBlobMask) sourcePlusBlobMask->Delete();
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
	if(sourceMask) sourceMask->Delete();
	if(curvMap) curvMap->Delete();
	if(blobMask) blobMask->Delete();
	if(sourcePlusBlobMask) sourcePlusBlobMask->Delete();

	return 0;

}//close FindNestedSources()


int Img::FindExtendedSource_CV(std::vector<Source*>& sources,BkgData* bkgData,int minPixels,bool findNegativeExcess,double dt,double h,double lambda1,double lambda2,double mu,double nu,double p){

	//## Compute segmented image
	Img* segmentedImg= ChanVeseSegmenter::FindSegmentation(this,false,dt,h,lambda1,lambda2,mu,nu,p);
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

	//## Remove sources of negative excess 
	if(!findNegativeExcess){
		if(!this->HasStats()){
			this->ComputeStats(true,false,false);
		}	
		if(!m_Stats){
			ERROR_LOG("Failed to compute image stats!");
			segmentedImg->Delete();
			for(unsigned int k=0;k<sources.size();k++){
				if(sources[k]){
					delete sources[k];
					sources[k]= 0;
				}	
			}//end loop sources
			sources.clear();
			return -1;
		}//close if

		std::vector<int> sourcesToBeRemoved;		
		int imgMedian= m_Stats->median;
		for(unsigned int k=0;k<sources.size();k++){
			double Smedian= sources[k]->Median;
			if(Smedian<imgMedian) sourcesToBeRemoved.push_back(k);
		}
		CodeUtils::DeleteItems(sources, sourcesToBeRemoved);
	}//close if

	//## Tag sources as extended
	for(unsigned int i=0;i<sources.size();i++){
		sources[i]->SetType(Source::eExtended);
	}//end loop sources

	segmentedImg->Delete();
	
	return 0;
	
}//close FindExtendedSources_CV()


int Img::FindExtendedSource_HClust(std::vector<Source*>&,Img* saliencyImg,Img* edgeImg){

	//Find SLIC partition
	//SLICData* slicData= SLIC::SPGenerator(this,int regionSize,double regParam, int minRegionSize, bool useLogScaleMapping, Img* edgeImg);

	//Tag regions
	//int SLICUtils::TagRegions(std::vector<Region*>& regions,Img* binaryMap_bkg,Img* binaryMap_signal)

	//Run segmentation
	//FindSegmentation(SLICData const& slicData,SLICData& segmSlicData,double SPMergingRegularization,bool SPMergingIncludeSpatialPars,bool use2ndNeighborsInSPMerging,int minMergedSP,double SPMergingRatio, double SPMergingMaxDissRatio,double SPMergingMaxDissRatio_2ndNeighbor,double SPMergingDissThreshold);
	

	return 0;

}//close FindExtendedSource_HClust()

int Img::ReadImageFile(std::string filename){

	//## Detect file extension
	std::string extension= filename.substr(filename.find_last_of(".") + 1);
	if(extension!= "png" && extension!="jpg" && extension!="bmp" && extension!="gif" ) {
		ERROR_LOG("Unknown file extension detected: ext="<<extension<<" (valid formats are png/jpg/bmp/gif)!");
		return -1;
	}
	
	//## Load image from file and set a matrix
	cv::Mat mat = cv::imread(filename.c_str(), CV_LOAD_IMAGE_COLOR);

	//## Convert to gray scale
	cv::Mat mat_gray;
  cvtColor( mat, mat_gray, CV_RGB2GRAY );
	//mat.convertTo(mat_gray, CV_32FC1);

	//## Fill an image
	int Nx= mat.cols;
	int Ny= mat.rows;
	
	this->SetBins(Nx,-0.5,Nx-0.5,Ny,-0.5,Ny-0.5);
	this->Sumw2();
	this->SetContour(999);
	
	for(int j=0;j<mat_gray.rows ;j++){//nBoxY
		int rowId= Ny-1-j;
		for(int i=0;i<mat_gray.cols;i++){
			//int colId= Nx-1-i;
			int colId= i;
			unsigned int matrixElement= mat_gray.at<uchar>(j,i);	
			//float matrixElement= mat_gray.at<float>(j,i);		
			int ix= colId ;
			int iy= rowId ;
			this->FillPixel(ix,iy,matrixElement);
		}
	}
	
	return 0;

}//close Img::ReadFile()

int Img::ReadFITS(std::string filename,int ix_min,int ix_max,int iy_min,int iy_max){

	//Start timer		
	auto start = chrono::steady_clock::now();

	//Read image or tile
	int status= 0;
	Caesar::FITSFileInfo fits_info;
	if(ix_min==-1 && ix_max==-1 && iy_min==-1 && iy_max==-1) {
		status= FITSReader::Read(filename,*this,fits_info);
	}
	else{
		status= FITSReader::ReadTileFast(filename,*this,fits_info,ix_min,ix_max,iy_min,iy_max);
	}
	
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

int Img::WriteFITS(std::string outfilename){
		
	//## Write fits to file
	if(FITSWriter::WriteFITS(this,outfilename)<0){
		ERROR_LOG("Failed to write image to FITS file!");
		return -1;
	}
	
	return 0;

}//close WriteFITS()


TH1D* Img::GetPixelHisto(int nbins,bool normalize){

	if(!HasStats()){
		WARN_LOG("No stats computed!");
		return 0;
		//cerr<<"Img::GetPixelHisto(): WARN: No stats computed!"<<endl;
		//this->ComputeStats(true,false,false);
	}

	double Smin= m_Stats->min;
	double Smax= m_Stats->max;
	double Srange= Smax-Smin;
	double tol= 0.0;
	double Smin_tol= Smin-tol*fabs(Srange);
	double Smax_tol= Smax+tol*fabs(Srange);

	TString histoName= Form("%s_histo",this->GetName());
	TH1D* histo= new TH1D(histoName,histoName,nbins,Smin_tol,Smax_tol);
		
	for(int i=0;i<this->GetNbinsX();i++){
		for(int j=0;j<this->GetNbinsY();j++){
			double w= this->GetBinContent(i+1,j+1);
			if(w==0) continue;
			histo->Fill(w);
		}//end loop bins Y
	}//end loop bins X

	if(normalize) histo->Scale(1./histo->Integral());
	return histo;

}//close GetPixelHisto()

double Img::FindOtsuThreshold(int nbins){
	
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
		//cout<<"t="<<t<<": wB="<<wB<<" wF="<<wF<<" sumB="<<sumB<<" mB="<<mB<<" mF="<<mF<<" var="<<varBetween<<endl;

   	// Check if new maximum found
   	if (varBetween > varMax) {
   		varMax = varBetween;
      threshold = binX;
			//cout<<"--> thr="<<threshold<<endl;
   	}
	}//end loop bins

	if(hist) {
		hist->Delete();
		hist= 0;
	}

	return threshold;

}//close FindOtsuThreshold()

double Img::FindValleyThreshold(int nbins,bool smooth){

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

double Img::FindMedianThreshold(double thrFactor){

	if(!HasStats()){
		this->ComputeStats(true,false,false);	
	}
	double median= m_Stats->median;
	double medianThr= thrFactor*median;
	
	return medianThr;

}//close FindMedianThreshold()

double Img::FindOptimalGlobalThreshold(double thrFactor,int nbins,bool smooth){

	if(!HasStats()){
		this->ComputeStats(true,false,false);	
	}
	double otsuThr= this->FindOtsuThreshold(nbins);
	double valleyThr= this->FindValleyThreshold(nbins,smooth);
	double medianThr= this->FindMedianThreshold(thrFactor);
	double optimalThr= std::max(std::min(otsuThr,valleyThr),medianThr);

	return optimalThr;

}//close FindOptimalGlobalThreshold()

Img* Img::GetBinarizedImage(double threshold,double fgValue,bool isLowerThreshold){

	TString imgName= Form("%s_Binarized",std::string(this->GetName()).c_str());	
	Img* BinarizedImg= this->GetCloned(std::string(imgName),true,true);
	BinarizedImg->Reset();

	for(int i=0;i<this->GetNbinsX();i++){
		double x= this->GetXaxis()->GetBinCenter(i+1);
		for(int j=0;j<this->GetNbinsY();j++){
			double y= this->GetYaxis()->GetBinCenter(j+1);
			double w= this->GetBinContent(i+1,j+1);
			if(w!=0){
				if(w>=threshold && !isLowerThreshold) {
					BinarizedImg->FillPixel(x,y,fgValue);
				}
				else if(w<threshold && isLowerThreshold) {
					BinarizedImg->FillPixel(x,y,fgValue);
				}
			}
		}//end loop bins Y
	}//end loop bins X
	
	return BinarizedImg;

}//close Img::GetBinarizedImage()





Img* Img::GetMask(Img* mask,bool isBinary){

	//## Check input mask
	if(!mask) {
		ERROR_LOG("Null ptr to given image mask!");
		return 0;
	}
		
	//## Check mask bins
	int Nx= mask->GetNbinsX();
	int Ny= mask->GetNbinsY();
	if(Nx!=this->GetNbinsX() || Ny!=this->GetNbinsY()){
		ERROR_LOG("Mask binning is different than current image!");
		return 0;
	}	

	//## Clone map
	TString imgName= Form("%s_Mask",std::string(this->GetName()).c_str());	
	Img* maskedImage= this->GetCloned(std::string(imgName),true,true);
	maskedImage->Reset();

	//## Loop over mask	
	if(isBinary){
		for(int i=0;i<Nx;i++){
			double x= mask->GetXaxis()->GetBinCenter(i+1);
			for(int j=0;j<Ny;j++){
				double y= mask->GetYaxis()->GetBinCenter(j+1);
				double binContent= this->GetBinContent(i+1,j+1);
				if(binContent==0) continue;
				double maskContent= mask->GetBinContent(i+1,j+1);
				if(maskContent!=0) maskedImage->FillPixel(x,y,1);
			}//end loop bins Y
		}//end loop bins X
	}//close if
	else{
		for(int i=0;i<Nx;i++){
			double x= mask->GetXaxis()->GetBinCenter(i+1);
			for(int j=0;j<Ny;j++){
				double y= mask->GetYaxis()->GetBinCenter(j+1);
				double binContent= this->GetBinContent(i+1,j+1);	
				if(binContent==0) continue;
				double maskContent= mask->GetBinContent(i+1,j+1);
				if(maskContent!=0) maskedImage->FillPixel(x,y,binContent);
			}//end loop bins Y
		}//end loop bins X
	}//close else

	return maskedImage;

}//close Img::GetMask()


Img* Img::GetSourceMask(std::vector<Source*>const& sources,bool isBinary,bool invert){

	//## Clone map
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	TString imgName= Form("%s_SourceMask",std::string(this->GetName()).c_str());	
	Img* maskedImage= this->GetCloned(std::string(imgName),true,true);
	
	//## Check source list
	int nSources= (int)sources.size();
	if(nSources<=0) {
		WARN_LOG("Source list is empty, returning same image!");
		return maskedImage;	
	}

	if(invert){
		if(isBinary){
			for(int k=0;k<nSources;k++){
				//PixelCollection thisSourcePixels= sources[k]->GetPixels();
				//for(unsigned int l=0;l<thisSourcePixels.size();l++){
				for(int l=0;l<sources[k]->GetNPixels();l++){
					//int id= thisSourcePixels[l]->id;
					int id= (sources[k]->GetPixel(l))->id;
					maskedImage->SetBinContent(id,0);
				}//end loop pixels
			}//end loop sources	

			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					double w= maskedImage->GetBinContent(i+1,j+1);		
					if(w==0) continue;
					maskedImage->SetBinContent(i+1,j+1,1);				
				}
			}
		}//close if
		else{
			for(int k=0;k<nSources;k++){
				//PixelCollection thisSourcePixels= sources[k]->GetPixels();
				//for(unsigned int l=0;l<thisSourcePixels.size();l++){
				for(int l=0;l<sources[k]->GetNPixels();l++){
					//int id= thisSourcePixels[l]->id;
					int id= (sources[k]->GetPixel(l))->id;
					maskedImage->SetBinContent(id,0);
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
			for(int k=0;k<nSources;k++){
				//PixelCollection thisSourcePixels= sources[k]->GetPixels();
				//for(unsigned int l=0;l<thisSourcePixels.size();l++){
				for(int l=0;l<sources[k]->GetNPixels();l++){
					Pixel* thisSourcePixel= sources[k]->GetPixel(l); 
					//double x= thisSourcePixels[l]->x;
					//double y= thisSourcePixels[l]->y;
					double x= thisSourcePixel->x;
					double y= thisSourcePixel->y;
					maskedImage->FillPixel(x,y,1);
				}//end loop pixels
			}//end loop sources		
		}//close if
		else{
			for(int k=0;k<nSources;k++){
				//PixelCollection thisSourcePixels= sources[k]->GetPixels();
				//for(unsigned int l=0;l<thisSourcePixels.size();l++){
				for(int l=0;l<sources[k]->GetNPixels();l++){
					Pixel* thisSourcePixel= sources[k]->GetPixel(l); 
					//int id= thisSourcePixels[l]->id;
					//double x= thisSourcePixels[l]->x;
					//double y= thisSourcePixels[l]->y;
					int id= thisSourcePixel->id;
					double x= thisSourcePixel->x;
					double y= thisSourcePixel->y;
					double binContent= this->GetBinContent(id);				
					maskedImage->FillPixel(x,y,binContent);
				}//end loop pixels
			}//end loop sources		
		}//close else
	}//close else

	return maskedImage;

}//close Img::GetSourceMask()


Img* Img::GetSourceResidual(std::vector<Source*>const& sources,int KernSize,int dilateModel,int dilateSourceType,bool skipToNested,BkgData* bkgData,bool useLocalBkg,bool randomize,double zThr){

	//Clone input image	
	TString imgName= Form("%s_Residual",std::string(this->GetName()).c_str());	
	Img* residualImg= this->GetCloned(std::string(imgName),true,true);
	if(!residualImg){
		ERROR_LOG("Failed to clone input image!");
		return 0;
	}

	//Dilate source pixels
	int status= MorphFilter::DilateAroundSources(residualImg,sources,KernSize,dilateModel,dilateSourceType,skipToNested,bkgData,useLocalBkg,randomize,zThr);

	if(status<0){
		ERROR_LOG("Failed to dilate sources!");
		if(residualImg) residualImg->Delete();
		return 0;		
	}

	//Re-Compute stats
	residualImg->ComputeStats(true,false,true);

	return residualImg;

}//close GetSourceResidual()


TMatrixD* Img::GetMatrix(){

	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	
	TMatrixD* M= new TMatrixD(Ny,Nx);
	M->Zero();

	for(int i=0;i<Nx;i++){//rows
		for(int j=0;j<Ny;j++){//columns
			double w= this->GetBinContent(i+1,j+1);
			int rowId= j;
			int colId= i;
			(*M)(rowId,colId)= w;
		}//end loop cols
	}//end loop rows

	return M;

}//close Img::GetMatrix()

cv::Mat Img::ImgToMat(std::string encoding){

	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();

	//## Fill OpenCV mat
	cv::Mat mat;
	if(encoding=="64") mat= cv::Mat::zeros(Ny,Nx,CV_64FC1);
	else if(encoding=="32") mat= cv::Mat::zeros(Ny,Nx,CV_32FC1);
	else{
		WARN_LOG("Invalid encoding selected, using default 64bit encoding");
		mat= cv::Mat::zeros(Ny,Nx,CV_64FC1);
	}

	/*
	//The slow way
	for(int j=0;j<Ny;j++){
		int rowId= Ny-1-j;
		for(int i=0;i<Nx;i++){
			int colId= i;
			mat.at<double>(rowId,colId)= this->GetBinContent(i+1,j+1);
		}//end loop x
	}//end loop y
	*/

	//The fast way
	int nRows = mat.rows;
  int nCols = mat.cols;
	for(int i=0;i<nRows;++i) {
		int rowId= i;
		int iy= Ny-1-rowId;
  	double* p = mat.ptr<double>(i);
    for (int j=0;j<nCols;++j){
			int colId= j;
			int ix= colId;
			double w= this->GetBinContent(ix+1,iy+1);
    	p[j] = w;
    }
  }

	return mat;

}//close ImgToMat()


Img* Img::GetNormalizedImage(std::string normScale,int normmin,int normmax,bool skipEmptyBins){

	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	double wmin= this->GetMinimum();
	double wmax= this->GetMaximum();
	
	TString imgName= Form("%s_Normalized",std::string(this->GetName()).c_str());	
	Img* norm_img= this->GetCloned(std::string(imgName),true,true);
	norm_img->Reset();
	
	if(normScale=="LINEAR"){
		for (int i=0;i<Nx; i++) { 
			double x= this->GetXaxis()->GetBinCenter(i+1);
			for (int j=0;j<Ny;j++) {	
				double y= this->GetYaxis()->GetBinCenter(j+1);
				double w= this->GetBinContent(i+1,j+1);	
				if(skipEmptyBins && w==0) continue;			
				double w_norm= normmin + (normmax-normmin)*(w-wmin)/(wmax-wmin);
				norm_img->FillPixel(x,y,w_norm);
			}//end loop yy
		}//end loop x
	}//close if
	else if(normScale=="LOG"){
		double safemin= 1;
		double safemax= 256;
		for (int i=0;i<Nx; i++) { 
			double x= this->GetXaxis()->GetBinCenter(i+1);
			for (int j=0;j<Ny; j++) {	
				double y= this->GetYaxis()->GetBinCenter(j+1);
				double w= this->GetBinContent(i+1,j+1);	
				if(skipEmptyBins && w==0) continue;			
				double w_norm= safemin + (safemax-safemin)*(w-wmin)/(wmax-wmin);
				double w_log= normmin + log10(w_norm/safemin)/log10(safemax/safemin) * (normmax-normmin);
				norm_img->FillPixel(x,y,w_log);
			}//end loop yy
		}//end loop x
	}//close else if
	else{
		WARN_LOG("Invalid norm scale option selected ("<<normScale<<") no transform applied to original image!");
	}
	
	return norm_img;

}//close Img::GetNormalizedImage()


Img* Img::GetGuidedFilterImage(int radius,double eps){

	//Normalize img
	Img* img_norm= this->GetNormalizedImage("LINEAR",1,256);
	img_norm->SetNameTitle("tmpImg","tmpImg");
	if(!img_norm) {
		ERROR_LOG("Failed to get normalized image!");
		return 0;
	}

	//## Convert image to OpenCV mat
	cv::Mat I= img_norm->ImgToMat("64");
	cv::Mat p= I;
	eps *= 255*255;   // Because the intensity range of our images is [0, 255]

	
	//## Run guided filter
	cv::Mat dst = Caesar::guidedFilter(I, p, radius, eps);

	//## Fill filtered image
	TString imgName= Form("%s_GuidedFilter",std::string(this->GetName()).c_str());
	Img* FilterImg= this->GetCloned(std::string(imgName),true,true);
	FilterImg->Reset();
	FilterImg->FillFromMat(dst);

	if(img_norm) img_norm->Delete(); 

	return FilterImg;

}//close GetGuidedFilterImage()


std::vector<Img*> Img::GetWaveletDecomposition(int nScales){

	DEBUG_LOG("Computing wavelet decomposition up to scale J="<<nScales<<" ...");
	std::vector<Img*> img_decomposition;
	img_decomposition= WTFilter::GetDecomposition(this,nScales);
	
	return img_decomposition;

}//close GetWaveletDecomposition()

Img* Img::GetKirschImage(){
	Img* kirschImg= Caesar::KirschFilter::GetKirschFilter(this);
	return kirschImg;
}//close Img::GetKirschImage()

Img* Img::GetLoGImage(bool invert){
	Img* LoGImg= Caesar::LoGFilter::GetLoGFilter(this);
	if(LoGImg && invert) LoGImg->Scale(-1);
	return LoGImg;
}//close Img::GetLoGImage()

Img* Img::GetNormLoGImage(int size,double scale,bool invert){
	Img* normLoGImg= Caesar::LoGFilter::GetNormLoGFilter(this,size,scale);
	if(normLoGImg && invert) normLoGImg->Scale(-1);
	return normLoGImg;
}//close Img::GetNormLoGImage()

Img* Img::GetLaplacianImage(bool invert){
	Img* laplImg= GradientFilter::GetLaplaceFilter(this);
	if(invert) laplImg->Scale(-1);
	return laplImg;
}//close Img::GetLaplacianImage()

Img* Img::GetGradientImage(bool invert){
	Img* gradImg= GradientFilter::GetGradientFilter(this);
	if(invert) gradImg->Scale(-1);
	return gradImg;
}//close Img::GetGradientImage()

Img* Img::GetSmoothedImage(int size_x,int size_y,double sigma_x,double sigma_y){

	//## Get OpenCV mat
	cv::Mat mat= this->ImgToMat("64");
	
	//## Smooth matrix
	cv::Size smooth_size(size_x,size_y);
	cv::Mat smoothed_mat;
	cv::GaussianBlur(mat,smoothed_mat, smooth_size, sigma_x, sigma_y, cv::BORDER_DEFAULT);

	//## Fill smoothed image
	TString imgName= Form("%s_Smoothed",std::string(this->GetName()).c_str());
	Img* SmoothedImg= this->GetCloned(std::string(imgName),true,true);
	SmoothedImg->Reset();
	SmoothedImg->FillFromMat(smoothed_mat);
	
	return SmoothedImg;

}//close GetSmoothedImage()


Img* Img::GetSaliencyMap(int reso,double regFactor,int minRegionSize,double knnFactor,bool useRobust,double expFalloffPar,double distanceRegPar){

	//## Compute single-reso saliency map
	Img* saliencyMap= 0;
	saliencyMap= SaliencyFilter::ComputeSaliencyMap(this,reso,regFactor,minRegionSize,knnFactor,useRobust,expFalloffPar,distanceRegPar);
	if(!saliencyMap){
		ERROR_LOG("Saliency map estimation failed!");
		return 0;
	}

	return saliencyMap;

}//close GetSaliencyMap()


Img* Img::GetMultiResoSaliencyMap(int resoMin,int resoMax,int resoStep,double beta,int minRegionSize,double knnFactor,bool useRobustPars,double expFalloffPar,double distanceRegPar,double salientMultiplicityThrFactor,bool addBkgMap,bool addNoiseMap,BkgData* bkgData,double saliencyThrFactor,double imgThrFactor){

	//## Compute multi-reso saliency map
	Img* saliencyMap= 0;
	saliencyMap= SaliencyFilter::ComputeMultiResoSaliencyMap(this,resoMin,resoMax,resoStep,beta,minRegionSize,knnFactor,useRobustPars,expFalloffPar,distanceRegPar, salientMultiplicityThrFactor,addBkgMap,addNoiseMap,bkgData,saliencyThrFactor,imgThrFactor);
	if(!saliencyMap){
		ERROR_LOG("Multi-resolution saliency map estimation failed!");
		return 0;
	}

	return saliencyMap;

}//close GetSaliencyMap()




int Img::Plot(std::vector<Source*>const& sources,bool useCurrentCanvas,bool drawFull,int paletteStyle,bool drawColorPalette,bool putWCSAxis,int coordSystem,std::string units){
	
	//Set palette
	int ncolors= 999;
	gStyle->SetNumberContours(ncolors);

	DEBUG_LOG("paletteStyle="<<paletteStyle);

	switch(paletteStyle){
		case eRAINBOW :
			gStyle->SetPalette(55);
			break;
		case eBLACKWHITE :	
			Caesar::GraphicsUtils::SetBWPalette(ncolors);
			break;
		case eBLACKBODY :
			gStyle->SetPalette(53);
			break;
		case eHOT2COLD :
			Caesar::GraphicsUtils::SetHotColdPalette(ncolors);
			break;
		case eCOLD2HOT :
			Caesar::GraphicsUtils::SetColdHotPalette(ncolors);
			break;
		case eTHERMAL :
			Caesar::GraphicsUtils::SetThermalPalette(ncolors);
			break;
		default: 
			gStyle->SetPalette(55);
			break;
	}//close switch

	/*
	if(paletteStyle=="RAINBOW") gStyle->SetPalette(55);
	else if(paletteStyle=="BLACKWHITE") Caesar::GraphicsUtils::SetBWPalette(ncolors);
	else if(paletteStyle=="BLACKBODY") gStyle->SetPalette(53);
	else if(paletteStyle=="HOT2COLD") Caesar::GraphicsUtils::SetHotColdPalette(ncolors);
	else if(paletteStyle=="COLD2HOT") Caesar::GraphicsUtils::SetColdHotPalette(ncolors);
	else if(paletteStyle=="THERMAL") Caesar::GraphicsUtils::SetThermalPalette(ncolors);
	else gStyle->SetPalette(55);
	*/

	//Set canvas
	TString canvasName= Form("%s_Plot",std::string(this->GetName()).c_str());	
	TCanvas* canvas= 0;
	if(useCurrentCanvas && gPad) {
		canvas= gPad->GetCanvas();
		canvas->SetName(canvasName);
		canvas->SetTitle(canvasName);
	}
	else{
		canvas= new TCanvas(canvasName,canvasName,720,700);
	}

	if(!canvas){
		ERROR_LOG("Failed to retrieve or set canvas!");
		return -1;
	}

	//Draw full image (without borders)
	canvas->cd();
	this->SetStats(0);
	
	if(drawFull){
		canvas->ToggleEventStatus();
  	canvas->SetRightMargin(0.0);
  	canvas->SetLeftMargin(0.0);
  	canvas->SetTopMargin(0.0);
  	canvas->SetBottomMargin(0.0);
		this->Draw("COLA");
	}
	else{
		gStyle->SetPadTopMargin(0.1);
  	gStyle->SetPadBottomMargin(0.1);
  	gStyle->SetPadLeftMargin(0.15);
  	//gStyle->SetPadRightMargin(0.19);
		gStyle->SetPadRightMargin(0.15);

		
		gPad->SetTopMargin(0.1);
		gPad->SetBottomMargin(0.1);
		gPad->SetLeftMargin(0.15);
		//gPad->SetRightMargin(0.19);
  	gPad->SetRightMargin(0.15);
  		
		if(drawColorPalette) {
			//Set palette axis title	
			this->GetZaxis()->SetTitle(units.c_str());
			this->GetZaxis()->SetTitleSize(0.05);
			this->GetZaxis()->SetTitleOffset(0.9);
			this->Draw();
			gPad->Update();
			if(putWCSAxis) this->Draw("COLAZ");
			else this->Draw("COLZ");
			gPad->Update();
		}//close if
		else {
			if(putWCSAxis) this->Draw("COLA");
			else this->Draw("COL");
		}
	}

	

	//## Draw WCS axis?
	if(putWCSAxis){
		gPad->Update();
		TGaxis* xaxis_wcs= new TGaxis;
		TGaxis* yaxis_wcs= new TGaxis;
		int status= GraphicsUtils::SetWCSAxis(this,*xaxis_wcs,*yaxis_wcs,coordSystem);
		if(status>=0){
			TExec* ex = new TExec("ex","GraphicsUtils::PadUpdater()");
   		this->GetListOfFunctions()->Add(ex);

			xaxis_wcs->Draw("same");
			yaxis_wcs->Draw("same");
		}
		else{
			WARN_LOG("Failed to set gAxis!");
		}
	}//close if

	//## Draw sources
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

}//close Img::Plot()


void Img::SetDrawRange(double zmin,double zmax){

	if(zmin>=zmax) return;

	for(int i=0;i<this->GetNbinsX();i++){
		for(int j=0;j<this->GetNbinsY();j++){
			double w= this->GetBinContent(i+1,j+1);
			if(w<zmin) this->SetBinContent(i+1,j+1,zmin);
			if(w>zmax) this->SetBinContent(i+1,j+1,zmax);
		}//end loop bins y
	}//end loop bins x

}//close Img::SetDrawRange()



}//close namespace

