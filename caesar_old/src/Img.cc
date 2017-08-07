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
#include <Region.h>
#include <Source.h>
#include <VLSlicSegmentation.h>
#include <SLICSegmenter.h>
#include <ChanVeseSegmentation.h>
#include <WTFilter.h>
#include <BkgFinder.h>
#include <Contour.h>
#include <FITSReader.h>
#include <Interpolator.h>
#include <ZernikeMoments.h>
//#include <wcs.h>

#include <guidedfilter.h>

#include <Rcpp.h>
#include <RInside.h>
using namespace Rcpp;

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
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
#include <TFitResult.h>

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
//#include <opencv2/ximgproc.hpp>
//#include <opencv2/ximgproc/edge_filter.hpp>

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
//using namespace cv;

ClassImp(Img)


Img::Img():TH2F()
{
	//cout<<"Img::Img(): INFO: Calling constructor 1..."<<endl;
  Init();
}//close costructor

Img::Img(const char *name,const char *title,int nbinsx,float xlow,float xup,int nbinsy,float ylow,float yup)
	:TH2F(name,title,nbinsx,xlow,xup,nbinsy,ylow,yup)
{
	//cout<<"Img::Img(): INFO: Calling constructor 2..."<<endl;
  Init();
}//close constructor

Img::Img(const char *name,const char *title,int nbinsx,const float *xbins,int nbinsy,const float *ybins)
	:TH2F(name,title,nbinsx,xbins,nbinsy,ybins)
{
	//cout<<"Img::Img(): INFO: Calling constructor 3..."<<endl;
  Init();
}//close constructor


Img::Img(const Img &img) : TH2F() {
	// Copy constructor
	//cout<<"Img::Img(): INFO: Calling copy constructor ..."<<endl;  
  ((Img&)img).Copy(*this);
}//close constructor


Img::~Img(){

	//cout<<"Img::~Img(): INFO: Deleting image "<<this->GetName()<<" ..."<<endl;
	
	//cout<<"Img::~Img(): Deleting pixel stats..."<<endl;
	if(fPixelStats) {
		delete fPixelStats;
		fPixelStats= 0;
	}

	//cout<<"Img::~Img(): Deleting pixel histo..."<<endl;
	if(fPixelHisto) {
		delete fPixelHisto;
		fPixelHisto= 0;
	}

	//cout<<"Img::~Img(): Deleting bkg data..."<<endl;
	for(unsigned int j=0;j<fBkgData.size();j++) {
		if(fBkgData[j]) {
			delete fBkgData[j];
			fBkgData[j]= 0;
		}
	}
	fBkgData.clear();

	
	//cout<<"Img::~Img(): Deleting source collection..."<<endl;
	for(unsigned int j=0;j<fSourceCollection.size();j++) {
		if(fSourceCollection[j]) {
			delete fSourceCollection[j];
			fSourceCollection[j]= 0;
		}
	}
	fSourceCollection.clear();
	

	//cout<<"Img::~Img(): Deleting bkg map..."<<endl;
	if(fInterpolatedBackgroundLevelMap) {
		delete fInterpolatedBackgroundLevelMap;
		fInterpolatedBackgroundLevelMap= 0;
	}
	
	//cout<<"Img::~Img(): Deleting bkg noise map..."<<endl;
	if(fInterpolatedBackgroundRMSMap) {
		delete fInterpolatedBackgroundRMSMap;
		fInterpolatedBackgroundRMSMap= 0;
	}
	
	//cout<<"Img::~Img(): Deleting bkg map..."<<endl;
	if(fBackgroundLevelMap) {
		delete fBackgroundLevelMap;
		fBackgroundLevelMap= 0;	
	}
	//cout<<"Img::~Img(): Deleting rms map..."<<endl;
	if(fBackgroundRMSMap) {
		delete fBackgroundRMSMap;
		fBackgroundRMSMap= 0;
	}
	//cout<<"Img::~Img(): INFO: End image "<<this->GetName()<<" delete."<<endl;
 	

}//close destructor


int Img::ReadFITSFile(std::string filename,int ix_min,int ix_max,int iy_min,int iy_max){
		
	FITSReader reader;	
	bool status= true;
	if(ix_min==-1 && ix_max==-1 && iy_min==-1 && iy_max==-1) {
		status= reader.Read(*this,filename);
	}
	else{
		status= reader.ReadTile(*this,filename,ix_min,ix_max,iy_min,iy_max);
	}
	
	if(!status){
		cerr<<"Img::ReadFITSFile(): ERROR: Failed to fill image from FITS file!"<<endl;
		return -1;
	}

	return 0;

}//close ReadFITSFile()


int Img::ReadFile(std::string filename){

	//## Detect file extension
	std::string extension= filename.substr(filename.find_last_of(".") + 1);
	if(extension!= "png" && extension!="jpg" && extension!="bmp" && extension!="gif" ) {
		cerr<<"Img::ReadFile(): ERROR: Unknown file extension detected: ext="<<extension<<" (valid ones are png/jpg/bmp/gif)!"<<endl;
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


bool Img::GetWorldCoord(int ix,int iy,double& xpos, double& ypos) {
	return Utils::PixelToWCSCoords(this,ix,iy,xpos,ypos);
}


void Img::Init(){

	fHasMetaData= false; 
	fHasStats= false;
	fPixelStats= 0;
	fPixelHisto= 0;
	fBackgroundLevelMap= 0;
	fBackgroundRMSMap= 0;
	fInterpolatedBackgroundLevelMap= 0;
	fInterpolatedBackgroundRMSMap= 0;
	
	fHasBkgData= false;
	fBkgData.clear();
	fBkgData.resize(0);
	fHasSources= false;
	fSourceCollection.clear();
	fSourceCollection.resize(0);

	fTileSizeX= 10;//in number of pixels
	fTileSizeY= 10;
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	
	fNTilesX= ceil( (double)(Nx)/(double)(fTileSizeX) );
	if(fNTilesX<=0) fNTilesX= 1;
	fNTilesY= ceil( (double)(Ny)/(double)(fTileSizeY) );
	if(fNTilesY<=0) fNTilesY= 1;

	
	ResetStats(true);

	//cout<<"Init(): initialize tiles..."<<endl;
	//InitTiles(fTileSizeX,fTileSizeY);

}//close Init()


void Img::InitTiles(int TileSizeX,int TileSizeY){

	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	if(Nx%TileSizeX!=0){
		cerr<<"Img::InitTiles(): WARNING: Tile size X does not fully match image size, last tile will have different size!"<<endl;		
	}
	if(Ny%TileSizeY!=0){
		cerr<<"Img::InitTiles(): WARNING: Tile size X does not fully match image size, last tile will have different size!"<<endl;
	}

	fTileSizeX= std::min(TileSizeX,Nx);
	fTileSizeY= std::min(TileSizeY,Ny);
	fNTilesX= ceil( (double)(Nx)/(double)(fTileSizeX) );
	if(fNTilesX<=0) fNTilesX= 1;
	fNTilesY= ceil( (double)(Ny)/(double)(fTileSizeY) );
	if(fNTilesY<=0) fNTilesY= 1;

}//close InitTiles()




void Img::SetPixel(int binx,int biny,double w, bool includeNegativePixelsInStats){
	
	//Check pixel value
	if( TMath::IsNaN(w) || fabs(w)==TMath::Infinity() ) return;

	//Check bin
	if(binx<0 || binx>=fXaxis.GetNbins()) return;
	if(biny<0 || biny>=fYaxis.GetNbins()) return;
	int bin  = biny*(fXaxis.GetNbins()+2) + binx;

	//Update Stats
	if( w>=0 || (w<0 && includeNegativePixelsInStats) ) {
		UpdateMoments(binx,biny,w);
	}

	// Set bin content
	SetBinContent(bin,w);

}//close Img::SetPixel()


void Img::UpdateMoments(int ix,int iy,double w){

	//Update full image moments
	if(w<fPixelMin) fPixelMin= w;
	if(w>fPixelMax) fPixelMax= w;

	fNpix++;
  double delta = w - fM1;
  double delta_n = delta/fNpix;
  double delta_n2 = delta_n * delta_n;
  double f = delta * delta_n * (fNpix-1);
  fM1+= delta_n;
  fM4+= f * delta_n2 * (fNpix*fNpix - 3*fNpix + 3) + 6 * delta_n2 * fM2 - 4 * delta_n * fM3;
  fM3+= f * delta_n * (fNpix - 2) - 3 * delta_n * fM2;
  fM2+= f;
	
}//close Img::UpdateMoments()


int Img::FillPixel(double x,double y,double w,bool includeNegativePixelsInStats){

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
			( w>=0 || (w<0 && includeNegativePixelsInStats) ) 
		) 
	{//update moments

		double w_old= GetBinContent(bin);

		if( w_old!=0 ){//bin already filled
			cerr<<"Img::FillPixel(): WARN: Pixel "<<bin<<" ("<<ix<<","<<iy<<") has been already filled, skipping..."<<endl;
			return -1;
		}//close if
		else{
			UpdateMoments(ix,iy,w);
		}
	}//close if
	
	//Copied from TH2::Fill(x,y,w)
	int filled_bin= Fill(x,y,w);
	
  return filled_bin;

}//close Img::FillPixel()


void Img::ResetStats(bool resetMoments){

	if(fPixelStats){
		fPixelStats->n= 0;
		fPixelStats->min= 0;
		fPixelStats->max= 0;
		fPixelStats->mean= 0;
		fPixelStats->meanErr= 0;	
		fPixelStats->rms= 0;
		fPixelStats->rmsErr= 0;
		fPixelStats->skewness= 0;
		fPixelStats->skewnessErr= 0;
  	fPixelStats->kurtosis= 0;
		fPixelStats->kurtosisErr= 0;
		fPixelStats->median= 0;
		fPixelStats->medianRMS= 0;
		fPixelStats->bwLocation= 0;
		fPixelStats->bwScale= 0;
	}//close if

	if(resetMoments){
		fPixelMin= +1.e+99;
		fPixelMax= -1.e+99;
		fNpix= 0;//npixels	
  	fM1= 0;//1st moments
  	fM2= 0;//2nd moment
		fM3= 0;//3rd moment
		fM4= 0;//4th moment
	}

}//close Img::ResetStats()



void Img::ComputeStatsParams(bool computeRobustStats,bool skipNegativePixels){

	//## Reset previous stats params (JUST STATS NOT MOMENTS!)
	ResetStats(false);
	
	//## Compute Stats params
	fPixelStats->n= fNpix;
	fPixelStats->min= fPixelMin;
	fPixelStats->max= fPixelMax;
	
	//if(fNpix<=1) return;

	fPixelStats->mean= fM1;
	fPixelStats->rms= 0;
	if(fNpix>2) fPixelStats->rms= sqrt(fM2/(fNpix-1));
	fPixelStats->skewness= 0;
	fPixelStats->kurtosis= 0;
	if(fM2!=0) {
  	fPixelStats->skewness= sqrt(fNpix)*fM3/pow(fM2,1.5);//need to adjust for finite population?
  	//fPixelStats->kurtosis= fNpix*fM4/(fM2*fM2)-3;//also given without -3
		fPixelStats->kurtosis= fNpix*fM4/(fM2*fM2);//also given without -3
	}

	fPixelStats->meanErr= 0;
	if(fNpix>0) fPixelStats->meanErr= fPixelStats->rms/sqrt(fNpix);
	double varianceErr= 0;
	fPixelStats->rmsErr= 0;
	if(fNpix>1) {
		varianceErr= (fM4-(fNpix-3)/(fNpix-1)*pow(fPixelStats->rms,4))/fNpix;
		fPixelStats->rmsErr= varianceErr/(2*fPixelStats->rms);
	}
		 
	fPixelStats->skewnessErr= 0;
	if(fNpix>2) fPixelStats->skewnessErr= sqrt(6.*fNpix*(fNpix-1)/((fNpix-2)*(fNpix+1)*(fNpix+3)));//approximate for normal distribution
	
  fPixelStats->kurtosisErr= 0;
	if(fNpix>3) {
		double kurtosisVariance= 24.*fNpix*(fNpix-1)*(fNpix-1)/((fNpix-3.)*(fNpix-2.)*(fNpix+3.)*(fNpix+5.));
		fPixelStats->kurtosisErr= sqrt(kurtosisVariance);
 	}

	fPixelStats->median= 0;
	fPixelStats->medianRMS= 0;
	
	//## End if no robust stats are to be computed
	if(!computeRobustStats) return;
	
	//## Compute robust stats (median, MAD, ...)	
	std::vector<double> Pixels;
	Pixels.clear();
	Pixels.resize(0);
	
	//## Fill list & histo of pixel values
	//Init pixel histo (filled later)
	if(!fPixelHisto) {
		fPixelHisto= new TH1D();
		fPixelHisto->Sumw2();
		fPixelHisto->SetNameTitle("PixelHisto","PixelHisto");
	}
	
	double minPixelVal= fPixelMin;
	double maxPixelVal= fPixelMax;
	double pixelRange= fabs(maxPixelVal-minPixelVal);
	double histoMinVal= minPixelVal - 0.2*pixelRange;
	double histoMaxVal= maxPixelVal + 0.2*pixelRange;
	fPixelHisto->SetBins(100,histoMinVal,histoMaxVal);
	fPixelHisto->Reset();

	for(int i=0;i<this->GetNbinsX();i++){
		for(int j=0;j<this->GetNbinsY();j++){
			double pixelValue= this->GetBinContent(i+1,j+1);
			if( pixelValue==0 || (skipNegativePixels && pixelValue<0) ) continue;
			Pixels.push_back(pixelValue);
			fPixelHisto->Fill(pixelValue);
		}
	}

	//cout<<"INFO: min/max="<<fPixelHisto->GetXaxis()->GetXmin()<<"/"<<fPixelHisto->GetXaxis()->GetXmax()<<endl;

	//Sort and compute median for all image	
	std::sort(Pixels.begin(),Pixels.end());
	double median= Utils::GetMedian(Pixels,true);
	fPixelStats->median= median;

	//Compute MAD = median(|x_i-median|)
	//cout<<"Img::ComputeStatsParams(): INFO: Computing MAD for full image..."<<endl;
	std::vector<double> MADs;
	for(unsigned j=0;j<Pixels.size();j++){
		double mad= fabs(Pixels[j]-median);
		MADs.push_back(mad);
	}
	std::sort(MADs.begin(),MADs.end());
	double medianMAD= Utils::GetMedian(MADs,true);
	double medianRMS= medianMAD*1.4826;//0.6744888;
	fPixelStats->medianRMS= medianRMS;

	//## Compute biweight robust estimators
	double C= 6.;
	double tol= 0.0001;
	double nmaxIter= 10;
	double biweightLocation= 0;
	double biweightScale= 0;
	double sumNumLocation= 0;
	double sumDenomLocation= 0;
	double sumNum= 0;
	double sumDenom= 0;
	double N= (double)(Pixels.size());

	double eps= 1.e+99;
	int niter= 0;	
	double Tb= median;
	double S= medianRMS;

	while(niter<nmaxIter){
		niter++;

		//Update location Tb		
		double Tb_old= Tb;
		double sumw= 0;
		double sum= 0;
		for(unsigned j=0;j<Pixels.size();j++){
			double x= Pixels[j];
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
		if(eps<tol) {
			//cout<<"Tb_old="<<Tb<<" Tb="<<Tb<<" eps="<<eps<<endl;
			break;
		}

		//Update scale
		double S_old= S;
		double sumNum= 0;
		double sumDenom= 0;
		for(unsigned j=0;j<Pixels.size();j++){
			double x= Pixels[j];
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
		
		//cout<<"--> ITER no. "<<niter<<": Tb="<<Tb<<" S="<<S<<" Eps="<<eps<<endl;

	}//end while

	fPixelStats->bwLocation= Tb;
	fPixelStats->bwScale= S;

	/*
	for(unsigned j=0;j<Pixels.size();j++){
		double x= Pixels[j]
		double u= (x-median)/(C*medianRMS);
		double u2= u*u;
		if(u2>1) continue;
		double num= ( pow(x-median,2) * pow(1.-u2,4) );
		double denom= (1.-u2)*(1.-5.*u2);
		sumNum+= num;
		sumDenom+= denom;
		sumNumLocation+= (x-median)*pow(1.-u2,2);
		sumDenomLocation+= pow(1.-u2,2);
	}

	double biweightVar= -1;
	if(sumDenom>0) biweightVar= N*sumNum/(sumDenom*(sumDenom-1.));
	if(biweightVar>0) biweightScale= sqrt(biweightVar);
	else biweightScale= 0;

	if(sumDenomLocation>0) biweightLocation= median + sumNumLocation/sumDenomLocation;
	fPixelStats->bwScale= biweightLocation;
	fPixelStats->bwScale= biweightScale;
	*/

	//## Compute clipped estimators
	double clipSigma= 3;
	int clipMaxIter= 100;
	double clipTolerance= 0.1;
	std::pair<double,double> clippedEstimators= Utils::GetClippedEstimators(Pixels,clipSigma,clipMaxIter,clipTolerance, true);
	double clippedMedian= clippedEstimators.first;
	double clippedRMS= clippedEstimators.second;
	fPixelStats->clippedMedian= clippedMedian;
	fPixelStats->clippedRMS= clippedRMS;

}//close Img::ComputeStatsParams()



int Img::ComputeStats(bool computeRobustStats,bool skipNegativePixels,bool forceRecomputing){

	//## Check if image has already stats computed
	if(!this->HasStats()){
		fPixelStats= new StatsData;
	}
	else{		
		cerr<<"Img::ComputeStats(): WARN: Image has already stats computed, recomputing stats!"<<endl;
	}

	//## If recomputing is not requested (i.e. some pixels has been reset by the user, just set the stats params!
	if(!forceRecomputing){
		//cout<<"Img::ComputeStats(): INFO: Computing stats params"<<endl;
		ComputeStatsParams(computeRobustStats,skipNegativePixels);
		fHasStats= true;
		return 0;
	}

	//## Recompute the moments and stats params
	//Reset stats & moments
	ResetStats(true);

	//Recompute moments
	std::vector<double> fPixels;
	fPixels.clear();
	fPixels.resize(0);
	for(int i=0;i<this->GetNbinsX();i++){
		for(int j=0;j<this->GetNbinsY();j++){
			double pixelValue= this->GetBinContent(i+1,j+1);
			if( pixelValue==0 || (skipNegativePixels && pixelValue<0) ) continue; 
			fPixels.push_back(pixelValue);
			UpdateMoments(i,j,pixelValue);
		}
	}

	//Recompute stats params
	ComputeStatsParams(true,skipNegativePixels);

	fHasStats= true;

	return 0;

}//close ComputeStats()


int Img::ComputeLocalBkg(Img::LocalBkgMethod method, Img::BkgMethod estimator, int boxSizeX, int boxSizeY, double boxSlideOffsetX, double boxSlideOffsetY,
int SPSize,double SPRegularization,int SPMinArea,bool useTwoPass){

	if(fHasBkgData){
		cout<<"Img::ComputeLocalBkg(): INFO: Image has already bkg data stored, recomputing bkg..."<<endl;
	}
	fBkgData.clear();
	fBkgData.resize(0);

	int status = 0;
	BkgFinder bkgFinder;
	bkgFinder.SetSegmentationRegionSize(SPSize);
	bkgFinder.SetSegmentationRegularization(SPRegularization);
	bkgFinder.SetSegmentationMinRegionArea(SPMinArea);

	BkgFinder::BkgMapData* bkgMapData= 0;
	if(method==eSuperpixelBkg){
		cout<<"Img::ComputeLocalBkg(): INFO: Using superpixel method..."<<endl;
		bkgMapData= bkgFinder.FindSuperpixelBkg(this, estimator, boxSizeX, boxSizeY, boxSlideOffsetX, boxSlideOffsetY);
	}
	else if(method==eGridBkg){	
		cout<<"Img::ComputeLocalBkg(): INFO: Using grid method..."<<endl;
		bkgMapData= bkgFinder.FindGridBkg(this, estimator, boxSizeX, boxSizeY, boxSlideOffsetX, boxSlideOffsetY);
	}
	else{
		cerr<<"Img::ComputeLocalBkg(): ERROR: Invalid local background method selected!"<<endl;
		return -1;
	}

	if(!bkgMapData){
		cerr<<"Img::ComputeLocalBkg(): ERROR: Computation of local background failed for this image!"<<endl;
		return -1;
	}
	
	fBkgData= bkgMapData->bkgData;
	
	TString imgName= Form("%s_InterpBkgMap",std::string(this->GetName()).c_str());
	fInterpolatedBackgroundLevelMap= bkgMapData->BkgMap;
	fInterpolatedBackgroundLevelMap->SetNameTitle(imgName,imgName);

	imgName= Form("%s_InterpRMSMap",std::string(this->GetName()).c_str());
	fInterpolatedBackgroundRMSMap= bkgMapData->RMSMap;
	fInterpolatedBackgroundRMSMap->SetNameTitle(imgName,imgName);

	//## Improve rms by recomputing stuff from residual map 
	if(useTwoPass){
		cout<<"Img::ComputeLocalBkg(): INFO: Improving rms estimation with a 2nd pass..."<<endl;
		Img* residualMap= (Img*)this->Clone("residualMap");
		residualMap->ResetBkg();//reset existing bkg info
		residualMap->Add(fInterpolatedBackgroundLevelMap,-1);//subtract the bkg level model
		BkgFinder::BkgMapData* bkgResidualMapData= 0;
		if(method==eSuperpixelBkg){
			bkgResidualMapData= bkgFinder.FindSuperpixelBkg(residualMap, estimator, boxSizeX, boxSizeY, boxSlideOffsetX, boxSlideOffsetY);
		}
		else if(method==eGridBkg){
			bkgResidualMapData= bkgFinder.FindGridBkg(this, estimator, boxSizeX, boxSizeY, boxSlideOffsetX, boxSlideOffsetY);
		}

		if(!bkgResidualMapData){
			cerr<<"Img::ComputeLocalBkg(): ERROR: Computation of local background failed @ second pass for this image!"<<endl;
			return -1;
		}
	
		//Update rms data in image data
		for(int i=0;i<(bkgResidualMapData->bkgData).size();i++) {
			fBkgData[i]->bkgRMS= (bkgResidualMapData->bkgData)[i]->bkgRMS;
		}
		fInterpolatedBackgroundRMSMap= bkgResidualMapData->RMSMap;

		residualMap->Delete();
	}//close if


	fHasBkgData= true;

	return 0;

}//close ComputeLocalBkg()


int Img::ComputeBkg(Img::BkgMethod method){

	if(fHasBkgData){
		cout<<"Img::ComputeBkg(): INFO: Image has already bkg data stored, recomputing bkg..."<<endl;
	}
	fBkgData.clear();
	fBkgData.resize(0);

	Img::BkgData* bkgData= 0;

	BkgFinder bkgFinder;
	bkgData= bkgFinder.ComputeBkg(this,method);
	if(!bkgData){
		cerr<<"Img::ComputeBkg(): ERROR: Computation of background failed for this image!"<<endl;
		return -1;
	}
	fBkgData.push_back(bkgData);
	fHasBkgData= true;	
	
	return 0;

}//close ComputeBkg()



int Img::GetPixelBkg(int ix,int iy,BkgData& bkgInfo,bool useLocalBackground){	
		
	//## Check if img has bkg data			
	if(!this->HasBkgData()) {
		cerr<<"Img::GetPixelBkg(): ERROR: No background data are available!"<<endl;
		return -1;
	}
	
	//Use global background computed?
	if(!useLocalBackground){
		bkgInfo.bkgLevel= fBkgData[0]->bkgLevel;
		bkgInfo.bkgRMS= fBkgData[0]->bkgRMS;
		return 0;
	}

	//Local background calculation
	if (!fInterpolatedBackgroundLevelMap || !fInterpolatedBackgroundRMSMap) {
		cout<<"Img::GetPixelBkg(): ERROR: Background & RMS maps are not available (did you compute local background?)..."<<endl;
		return -1;	
	}

	bkgInfo.bkgLevel= fInterpolatedBackgroundLevelMap->GetBinContent(ix+1,iy+1);
	bkgInfo.bkgRMS= fInterpolatedBackgroundRMSMap->GetBinContent(ix+1,iy+1);
	
	return 0;

}//close GetPixelBkg()


Img* Img::GetSignificanceMap(bool useLocalBackground){

	//Check bkg data
	if(!this->HasBkgData()){
		cout<<"Img::GetSignificanceMap(): INFO: Image has no bkg information stored (hints: compute bkg first!), nothing to be done!"<<endl;
		return 0;
	}

	//Clone this image and reset content
	TString imgName= Form("%s_SignificanceMap",std::string(this->GetName()).c_str());
	Img* significanceMap= (Img*)this->Clone(imgName);
	significanceMap->Reset();

	//Fill significance map	
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	BkgData bkgInfo;
	bkgInfo.ix_min= -1;
	bkgInfo.iy_min= -1;
	bkgInfo.ix_max= -1;
	bkgInfo.iy_max= -1;
	bkgInfo.npix= 0;
	bkgInfo.bkgLevel= 0;
	bkgInfo.bkgRMS= 0;
	bkgInfo.isReliableBkg= true;

	for(int i=0;i<Nx;i++){
		double binX= this->GetXaxis()->GetBinCenter(i+1);
		
		for(int j=0;j<Ny;j++){	
			double binY= this->GetYaxis()->GetBinCenter(j+1);

			double w= this->GetBinContent(i+1,j+1);
			int status= GetPixelBkg(i,j,bkgInfo,useLocalBackground);
			if(status<0) {
				//cerr<<"Img::GetSignificanceMap(): WARN: Cannot get background for pixel ("<<i<<","<<j<<") skip it!"<<endl;
				continue;
			}
			
			double bkgLevel= bkgInfo.bkgLevel;
			double bkgRMS= bkgInfo.bkgRMS;
			bool isReliableBkg= bkgInfo.isReliableBkg;
			/*
			if( w==0 || (w<0 && skipNegativePixels) ||
					status<0 || !isReliableBkg || bkgRMS<=0
			) {
				continue;
			}
			*/
			if( w==0 || !isReliableBkg || bkgRMS<=0) {
				//cerr<<"Img::GetSignificanceMap(): WARN: Empty pixel or invalid bkg  ("<<i<<","<<j<<") skip it!"<<endl;
				continue;
			}
				
			double Z= (w-bkgLevel)/bkgRMS;
			significanceMap->FillPixel(binX,binY,Z);
		}//end loop
	}//end loop 

	return significanceMap;

}//close Img::GetSignificanceMap()

Img* Img::GetSignificanceMap(Img* BkgMap, Img* NoiseMap){

	//Check image
	if(!BkgMap || !NoiseMap){
		cerr<<"Null ptr to bkg/noise maps!"<<endl;
		return 0;
	}

	//Integrity checks
	long int Nx= this->GetNbinsX();
	long int Ny= this->GetNbinsY();
	long int Nx_bkg= BkgMap->GetNbinsX();
	long int Ny_bkg= BkgMap->GetNbinsY();
	long int Nx_noise= NoiseMap->GetNbinsX();
	long int Ny_noise= NoiseMap->GetNbinsY();
	if( Nx!=Nx_bkg || Ny!=Ny_bkg ||
			Nx!=Nx_noise || Ny!=Ny_noise ||
			Nx_bkg!=Nx_noise || Ny_bkg!=Ny_noise
	)
	{
		cerr<<"Img::GetSignificanceMap(): ERROR: Bkg/Noise/image maps have different size!"<<endl;			
		return 0;
	}
	
		
	//Clone this image and reset content
	TString imgName= Form("%s_significance",std::string(this->GetName()).c_str());
	Img* significanceMap= (Img*)this->Clone(imgName);
	significanceMap->Reset();

	for(int i=0;i<Nx;i++){
		double binX= this->GetXaxis()->GetBinCenter(i+1);
		for(int j=0;j<Ny;j++){	
			double binY= this->GetYaxis()->GetBinCenter(j+1);
			double w= this->GetBinContent(i+1,j+1);	
											
			double bkgLevel= BkgMap->GetBinContent(i+1,j+1);
			double bkgRMS= NoiseMap->GetBinContent(i+1,j+1);
			if( w==0 || bkgRMS<=0) continue;
			
			double Z= (w-bkgLevel)/bkgRMS;
			significanceMap->FillPixel(binX,binY,Z);
		}//end loop
	}//end loop 
	

	return significanceMap;

}//close GetSignificanceMap()


Source* Img::FindSeededBlob(int seedPixelId,double mergeThr,bool mergeBelowSeed){

	//## Init
	Source* aBlob= 0;
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	int Ntot= (Nx+2)*(Ny+2);

	//## Run flood-fill 	
	std::vector<int> clusterPixelIds= this->FloodFill(seedPixelId,mergeThr);
	int nClusterPixels= (int)clusterPixelIds.size();
	if(nClusterPixels<=0){
		cout<<"Img::FindSeededBlob(): INFO: No blob detected!"<<endl;
		return 0;
	}

	//## Create a source
	std::string sourceName= std::string(Form("Img%s_blobSeedId%d",std::string(this->GetName()).c_str(),seedPixelId));
		
	aBlob= new Source;
	aBlob->SetId(seedPixelId);	
	aBlob->SetName(sourceName);
		
	for(int l=0;l<nClusterPixels;l++){
		int clusterPixelId= clusterPixelIds[l];	
		int clusterPixelIdX, clusterPixelIdY, clusterPixelIdZ;
		this->GetBinXYZ(clusterPixelId,clusterPixelIdX,clusterPixelIdY,clusterPixelIdZ);
		double S= this->GetBinContent(clusterPixelId);			
		double Z= 0;
		double x= this->GetXaxis()->GetBinCenter(clusterPixelIdX);
		double y= this->GetYaxis()->GetBinCenter(clusterPixelIdY);
		double Curv= 0;
		int ix= clusterPixelIdX-1;
		int iy= clusterPixelIdY-1;

		double bkgLevel= 0;
		double noiseLevel= 0;	
	
		Source::Pixel clusterPixel;
		clusterPixel.S= S;
		clusterPixel.Type= Source::eNormal;
		clusterPixel.Z= Z;
		clusterPixel.Curv= Curv;
		clusterPixel.BkgLevel= bkgLevel;
		clusterPixel.NoiseLevel= noiseLevel;
		clusterPixel.id= clusterPixelId;
		clusterPixel.ix= clusterPixelIdX-1;
		clusterPixel.iy= clusterPixelIdY-1;
		clusterPixel.x= x;
		clusterPixel.y= y;
		if( clusterPixelIdX<=1 || clusterPixelIdY<=1 || clusterPixelIdX>=Nx || clusterPixelIdY>=Ny) 
			aBlob->SetEdgeFlag(true);

		aBlob->AddPixel(clusterPixel);
	}//end loop cluster pixels

	return aBlob;

}//close Img::FindSeededBlob()

std::vector<int> Img::FloodFill(int seedPixelId,double mergeThr,bool mergeBelowSeed){

	std::vector<int> clusterPixelIds;
	clusterPixelIds.clear();

	//Check given seed id
	if(this->IsBinOverflow(seedPixelId) || this->IsBinUnderflow(seedPixelId) ){
		cerr<<"Img::FloodFill(): WARN: Given seed id is outside image range!"<<endl;
		return clusterPixelIds;
	}
	double seedSignal= this->GetBinContent(seedPixelId);

	//Add seed to queue and loop over queue
	std::queue<int> pixelQueue;
	pixelQueue.push(seedPixelId);

	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	int Ntot= (Nx+2)*(Ny+2);
	std::vector<bool> isAddedInQueue(Ntot,false);	
	std::vector<bool> isAddedInCluster(Ntot,false);

	//cout<<"****** FLOOD-FILL ************"<<endl;
	//cout<<"seedPixelId="<<seedPixelId<<endl;

	while(!pixelQueue.empty()){

		//Take first pixel in queue, process it and then remove from the queue
		int gBinId= pixelQueue.front(); 
		int binIdX, binIdY, binIdZ; 
		this->GetBinXYZ(gBinId,binIdX,binIdY,binIdZ);
		pixelQueue.pop();

		if(mergeBelowSeed){
			//Loop on row pixels above threshold
    	while (binIdX-1>0 && fabs(this->GetBinContent(binIdX-1,binIdY))>=mergeThr && fabs(this->GetBinContent(binIdX-1,binIdY))<=seedSignal) {
    		binIdX--;
    	}//close while loop
    
			bool spanUp = false;
    	bool spanDown = false;
		 
    	while (binIdX<=Nx && fabs(this->GetBinContent(binIdX,binIdY))>=mergeThr && fabs(this->GetBinContent(binIdX,binIdY))<=seedSignal) {
				int gBinId_cluster= this->GetBin(binIdX,binIdY);
				if(!isAddedInCluster[gBinId_cluster]) {
					clusterPixelIds.push_back(gBinId_cluster);
					isAddedInCluster[gBinId_cluster]= true;
				}
			
				//search up pixel
				int gBinId_up= this->GetBin(binIdX,binIdY+1);

      	if (!spanUp && binIdY+1<=Ny && fabs(this->GetBinContent(binIdX,binIdY+1))>=mergeThr && fabs(this->GetBinContent(binIdX,binIdY+1))<=seedSignal) {
      		if(!isAddedInQueue[gBinId_up]) {
						pixelQueue.push(gBinId_up);
						isAddedInQueue[gBinId_up]= true;
						spanUp = true;
					} 
				} 
				else if (spanUp && binIdY+1<=Ny && (fabs(this->GetBinContent(binIdX,binIdY+1))<mergeThr || fabs(this->GetBinContent(binIdX,binIdY+1))>seedSignal) ) {
      		spanUp = false;
      	}

				//search down pixel
				int gBinId_down= this->GetBin(binIdX,binIdY-1);

     		if (!spanDown && binIdY-1>0 && fabs(this->GetBinContent(binIdX,binIdY-1))>=mergeThr && fabs(this->GetBinContent(binIdX,binIdY-1))<=seedSignal) {
      		if(!isAddedInQueue[gBinId_down]) {
						pixelQueue.push(gBinId_down);
						isAddedInQueue[gBinId_down]= true;
						spanDown = true;
					} 
      	} 
				else if (spanDown && binIdY-1>0 && (fabs(this->GetBinContent(binIdX,binIdY-1))<mergeThr || fabs(this->GetBinContent(binIdX,binIdY-1))>seedSignal) ) {
      		spanDown = false;
      	}
      	binIdX++;
			}//end while loop
		
		}//close if
		else{


			//Loop on row pixels above threshold
    	while (binIdX-1>0 && fabs(this->GetBinContent(binIdX-1,binIdY))>=mergeThr) {
    		binIdX--;
   	 	}//close while loop
    
			bool spanUp = false;
    	bool spanDown = false;
		 
    	while (binIdX<=Nx && fabs(this->GetBinContent(binIdX,binIdY))>=mergeThr) {
				int gBinId_cluster= this->GetBin(binIdX,binIdY);
				if(!isAddedInCluster[gBinId_cluster]) {
					clusterPixelIds.push_back(gBinId_cluster);
					isAddedInCluster[gBinId_cluster]= true;
				}
			
				//search up pixel
				int gBinId_up= this->GetBin(binIdX,binIdY+1);

    	  if (!spanUp && binIdY+1<=Ny && fabs(this->GetBinContent(binIdX,binIdY+1))>=mergeThr) {
    	  	if(!isAddedInQueue[gBinId_up]) {
						pixelQueue.push(gBinId_up);
						isAddedInQueue[gBinId_up]= true;
						spanUp = true;
					} 
				} 
				else if (spanUp && binIdY+1<=Ny && fabs(this->GetBinContent(binIdX,binIdY+1))<mergeThr) {
    	  	spanUp = false;
    	  }

				//search down pixel
				int gBinId_down= this->GetBin(binIdX,binIdY-1);

     		if (!spanDown && binIdY-1>0 && fabs(this->GetBinContent(binIdX,binIdY-1))>=mergeThr) {
      		if(!isAddedInQueue[gBinId_down]) {
						pixelQueue.push(gBinId_down);
						isAddedInQueue[gBinId_down]= true;
						spanDown = true;
					} 
      	} 
				else if (spanDown && binIdY-1>0 && fabs(this->GetBinContent(binIdX,binIdY-1))<mergeThr) {
      		spanDown = false;
      	}
      	binIdX++;		
  		}//end while loop 
		}//close else
	}//end queue loop

	//### END SEARCH ####
		
	//Append cluster pixels to a source object
	int nClusterPixels= clusterPixelIds.size();
	if(nClusterPixels==0) cout<<"Img::FloodFill(): INFO: No cluster pixels found for this seed!"<<endl;

	//## DEBUG
	//cout<<"("; 
	//for(int i=0;i<nClusterPixels;i++){
	//	cout<<clusterPixelIds[i]<<",";
	//}
	//cout<<")"<<endl;
	//cout<<"**************"<<endl;
	//########
	
	return clusterPixelIds;

}//close Img::FloodFill()



std::vector<Source*> Img::FindBlobs(Img* significanceMap,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool useLocalBackground,bool mergeBelowSeed){

	Source* aBlob= 0;
	std::vector<Source*> blobs;
	blobs.clear();

	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	int Ntot= (Nx+2)*(Ny+2);

	//## Check if background data are available for this image
	if(!fHasBkgData){
		cerr<<"Img::FindBlobs(): WARN: No bkg data available for this image, no bkg info will be stored to source!"<<endl;
	}
	BkgData bkgInfo;
	bkgInfo.ix_min= -1;
	bkgInfo.iy_min= -1;
	bkgInfo.ix_max= -1;
	bkgInfo.iy_max= -1;
	bkgInfo.npix= 0;
	bkgInfo.bkgLevel= 0;
	bkgInfo.bkgRMS= 0;
	bkgInfo.isReliableBkg= true;

	//## Check for given significance map 
	if(!significanceMap){
		cout<<"Img::FindBlobs(): ERROR: Null ptr to given significance map!"<<endl;
		return blobs;
	}
	
	//## Find curvature map (it is added to found sources)
	cout<<"Img::FindBlobs(): INFO: Finding image curvature..."<<endl;
	Img* curvatureMap= this->GetLoGImage(true);//invert map
	if(!curvatureMap) {
		cout<<"Img::FindBlobs(): ERROR: Null ptr to computed curvature map!"<<endl;
		return blobs;
	}
	
	//## Find seed pixels (above seed threshold)	
	std::vector<int> pixelSeeds;
	std::vector<bool> isNegativeExcessSeed;
	
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){	
			int gBin= significanceMap->GetBin(i+1,j+1);	
			double Z= significanceMap->GetBinContent(i+1,j+1);
			bool isNegative= false;
			if(fabs(Z)>=seedThr) {
				if(Z<0) isNegative= true;
				pixelSeeds.push_back(gBin);	
				isNegativeExcessSeed.push_back(isNegative);
			}
		}//end loop y
	}//end loop x
	
	cout<<"Img::FindBlobs(): INFO: "<<pixelSeeds.size()<<" seeds found ..."<<endl;

	//## Perform cluster finding starting from detected seeds
	int nSources= 0;
	std::vector<bool> isAddedInCluster(Ntot,false);

	for(unsigned int k=0;k<pixelSeeds.size();k++){
		int seedPixelId= pixelSeeds[k];
		int binX, binY, binZ;
		significanceMap->GetBinXYZ(seedPixelId,binX,binY,binZ);
		
		if(isAddedInCluster[seedPixelId]) continue;
		
		//Skip negative excess seed if not requested
		if(!findNegativeExcess && isNegativeExcessSeed[k]) continue;
		
		//Compute flooded pixels
		std::vector<int> clusterPixelIds= significanceMap->FloodFill(seedPixelId,mergeThr,mergeBelowSeed);

		//Append cluster pixels to a source object
		int nClusterPixels= (int)clusterPixelIds.size();
		if(nClusterPixels==0 || nClusterPixels<minPixels) {
			cout<<"Img::FindBlobs(): INFO: Blob pixels found @ (x,y)=("<<binX<<","<<binY<<") (N="<<nClusterPixels<<") below npix threshold (thr="<<minPixels<<"), skip blob!"<<endl; 
			continue;
		}
		cout<<"Img::FindBlobs(): INFO: Blob found @ (x,y)=("<<binX<<","<<binY<<") (N="<<nClusterPixels<<")"<<endl;

		nSources++;
		std::string sourceName= std::string(Form("Img%s_blobId%d",std::string(this->GetName()).c_str(),nSources));
		cout<<"Img::FindBlobs(): INFO: Adding new blob (# "<<nSources<<") to list "<<endl;
		
		aBlob= new Source;
		aBlob->SetId(nSources);	
		aBlob->SetName(sourceName);
		
		for(int l=0;l<nClusterPixels;l++){
			int clusterPixelId= clusterPixelIds[l];	
			if(isAddedInCluster[clusterPixelId]) continue;
			isAddedInCluster[clusterPixelId]= true;//do not forget to add to list of taken pixels!
			int clusterPixelIdX, clusterPixelIdY, clusterPixelIdZ;
			this->GetBinXYZ(clusterPixelId,clusterPixelIdX,clusterPixelIdY,clusterPixelIdZ);
			double S= this->GetBinContent(clusterPixelId);			
			double Z= significanceMap->GetBinContent(clusterPixelId);
			double x= this->GetXaxis()->GetBinCenter(clusterPixelIdX);
			double y= this->GetYaxis()->GetBinCenter(clusterPixelIdY);
			double Curv= curvatureMap->GetBinContent(clusterPixelId);
			int ix= clusterPixelIdX-1;
			int iy= clusterPixelIdY-1;

			double bkgLevel= 0;
			double noiseLevel= 0;	
	
			if(fHasBkgData && GetPixelBkg(ix,iy,bkgInfo,useLocalBackground)==0) {
				bkgLevel= bkgInfo.bkgLevel;
				noiseLevel= bkgInfo.bkgRMS;
			}

			Source::Pixel clusterPixel;
			clusterPixel.S= S;
			if(Z>=seedThr) clusterPixel.Type= Source::eSeed;
			else clusterPixel.Type= Source::eNormal;
			clusterPixel.Z= Z;
			clusterPixel.Curv= Curv;
			clusterPixel.BkgLevel= bkgLevel;
			clusterPixel.NoiseLevel= noiseLevel;
			clusterPixel.id= clusterPixelId;
			clusterPixel.ix= clusterPixelIdX-1;
			clusterPixel.iy= clusterPixelIdY-1;
			clusterPixel.x= x;
			clusterPixel.y= y;
			if( clusterPixelIdX<=1 || clusterPixelIdY<=1 || clusterPixelIdX>=Nx || clusterPixelIdY>=Ny) 
				aBlob->SetEdgeFlag(true);

			aBlob->AddPixel(clusterPixel);
		}//end loop cluster pixels

		blobs.push_back(aBlob);
		
	}//end loop seeds

	cout<<"Img::FindBlobs(): INFO: #"<<blobs.size()<<" blobs found!"<<endl;

	if(curvatureMap) curvatureMap->Delete(); 			

	return blobs;

}//close Img::FindBlobs()


int Img::FindCompactSource(Img* significanceMap,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool findNestedSources,bool useLocalBackground,bool mergeBelowSeed,double peakThreshold){
	
	//## Find significant blobs in image
	std::vector<Source*> Blobs= FindBlobs(significanceMap,seedThr,mergeThr,minPixels,findNegativeExcess,useLocalBackground,mergeBelowSeed);
	int nBlobs= (int)Blobs.size();
	if(nBlobs<=0){
		cout<<"Img::FindCompactSource(): INFO: No blobs found for image: "<<this->fName<<endl;
		return 0;
	}
	
	double fluxCorrection= 1;
	if(this->HasMetaData()){
		double fx= this->metadata.Bmaj;
		double fy= this->metadata.Bmin;
		fluxCorrection= TMath::Pi()*fx*fy/(4*log(2));
	}

	//## Loop over detected blobs and select valid source candidates
	int nSources= 0;
	cout<<"Img::FindCompactSource(): INFO: "<<nBlobs<<" blobs detected..."<<endl;

	for(int k=0;k<nBlobs;k++){
		if(!Blobs[k]) continue;

		//## Set source name & id & type
		std::string sourceName= std::string(Form("Img%s_sId%d",std::string(this->GetName()).c_str(),k));
		
		//Blobs[k]->SetId(nSources);	
		Blobs[k]->SetId(k);
		Blobs[k]->SetName(sourceName);
		Blobs[k]->SetType(Source::eCompact);

		//## Compute stats parameters
		cout<<"Img::FindCompactSource(): INFO: Computing source stats..."<<endl;
		Blobs[k]->ComputeStats();

		//## Compute morphology parameters
		cout<<"Img::FindCompactSource(): INFO: Computing morphology params..."<<endl;
		Blobs[k]->ComputeMorphologyParams();
	
		//## Set flux correction factor
		Blobs[k]->fFluxCorrection= fluxCorrection;

		cout<<"Img::FindCompactSource(): INFO: Adding new source (# "<<nSources<<") to list..."<<endl;
		fSourceCollection.push_back(Blobs[k]);
		nSources++;
	}//end loop blobs

	if(nSources>0) fHasSources= true;
	cout<<"Img::FindCompactSource(): INFO: "<<nSources<<" sources found!"<<endl;

	//##### FIND NESTED ########
	if(findNestedSources && nSources>0){
		cout<<"Img::FindCompactSource(): INFO: Finding nested sources..."<<endl;

		//## Find image mask of found sources
		Img* sourceMask= this->GetSourceMask(fSourceCollection,false);
		if(!sourceMask){
			cerr<<"Img::FindCompactSource(): ERROR: Null ptr to computed source mask!"<<endl;
			return -1;
		}
		Img* sourceMask_binary= this->GetSourceMask(fSourceCollection,true);
		if(!sourceMask_binary){
			cerr<<"Img::FindCompactSource(): ERROR: Null ptr to computed binary source mask!"<<endl;
			return -1;
		}

		//## Find curvature map (it is added to found sources)
		cout<<"Img::FindCompactSource(): INFO: Finding image curvature..."<<endl;
		//Img* curvatureMap= this->GetLoGImage(true);//invert map

		int kernFWHM= this->GetBeamSizeInPixel();
		double kernSigma= kernFWHM/(2*sqrt(2*log(2)));
		int kernFactor= 6;
		int kernSize= std::round(kernFactor*kernSigma);
		cout<<"Img::FindCompactSource(): INFO: kernFWHM="<<kernFWHM<<", kernSigma="<<kernSigma<<" kernSize="<<kernSize<<endl;
		Img* curvatureMap= this->GetNormLoGImage(kernSize,kernSigma,true);//invert map
		if(!curvatureMap) {
			cerr<<"Img::FindCompactSource(): ERROR: Null ptr to computed curvature map!"<<endl;
			return -1;
		}
		curvatureMap->ComputeStats(true,false,true);
		double curvRMS= (curvatureMap->GetPixelStats())->medianRMS;
	
		//## Apply threshold to curv map		
		//double curvThr= 0;
		double curvThr= 1.*fabs(curvRMS);//threshold is +1 sigma
		double bkgValue= 0; 
		double fgValue= 1;
		Img* curvatureMap_binary= curvatureMap->GetBinarized(curvThr,bkgValue,fgValue);

		//## Get source+curvature mask		
		//Img* sourcePlusCurvatureMask= curvatureMap->GetMask(sourceMask_binary,false);
		Img* sourcePlusCurvatureMask= sourceMask->GetMask(curvatureMap_binary,false);
		if(!sourcePlusCurvatureMask){
			cerr<<"Img::FindCompactSource(): ERROR: Null ptr to computed (source+curvature) mask!"<<endl;
			return -1;
		}

		//## Find peaks in curvature map
		int tol= 1;//pixel-tolerance
		std::vector<int> peakIds= sourcePlusCurvatureMask->FindPeakIds(tol);
		int nPeaks= (int)peakIds.size();

		int nGoodPeaks= 0;
		for(int k=0;k<nPeaks;k++){
			int peakBinId= peakIds[k];
			
			double peakZ= significanceMap->GetBinContent(peakBinId);
			if(peakZ<peakThreshold) continue;//skip peaks below threshold
			int peakBinX, peakBinY, peakBinZ;
			curvatureMap_binary->GetBinXYZ(peakBinId,peakBinX,peakBinY,peakBinZ);
			cout<<"--> Good peak @ ("<<curvatureMap_binary->GetXaxis()->GetBinCenter(peakBinX)<<","<<curvatureMap_binary->GetYaxis()->GetBinCenter(peakBinY)<<")"<<endl;
			curvatureMap_binary->AddBinContent(peakBinId,1);
			nGoodPeaks++;
		}	
		cout<<"Source::FindNestedSource(): INFO: #"<<nGoodPeaks<<"/"<<nPeaks<< " good peaks detected!"<<endl;

		std::vector<Source*> NestedSources;
		bool mergeBelowSeed_nested= false;
		if(nGoodPeaks>0) NestedSources= sourceMask->FindBlobs(curvatureMap_binary,fgValue+1,fgValue,minPixels,findNegativeExcess,useLocalBackground,mergeBelowSeed_nested);

		int nNestedSources= (int)NestedSources.size();
		if(nNestedSources>=0){
			cout<<"Img::FindCompactSource(): INFO: #"<<nNestedSources<<" nested sources found!"<<endl;

			//## Find matching between mother and nested sources
			for(int j=0;j<nNestedSources;j++){
				bool isMotherFound= false;
				cout<<"Img::FindCompactSource(): INFO: Finding matching for nested source no. "<<j<<endl;
				NestedSources[j]->SetType(Source::eCompact);
				NestedSources[j]->ComputeStats();
				NestedSources[j]->ComputeMorphologyParams();
				for(int i=0;i<nSources;i++){
					int sourceId= fSourceCollection[i]->fId;
					bool isInside= NestedSources[j]->IsInsideSource(fSourceCollection[i]);
					if(isInside){
						cout<<"Img::FindCompactSource(): INFO: Nested source no. "<<j<<" added to source id="<<sourceId<<" ..."<<endl;
						fSourceCollection[i]->AddNestedSource(NestedSources[j]);
						isMotherFound= true;
						break;
					}
				}//end loop mother sources
				if(!isMotherFound){
					cerr<<"Img::FindCompactSource(): WARN: Cannot find mother source for nested source no. "<<j<<"!"<<endl;
					NestedSources[j]->Dump();
				}			
			}//end loop nested sources
								
		}//close nNestedBlobs>0
		
		//Clear up
		if(curvatureMap) {
			delete curvatureMap;
			curvatureMap= 0;
		}
		if(curvatureMap_binary) {
			delete curvatureMap_binary;
			curvatureMap_binary= 0;
		}
		if(sourceMask) {
			delete sourceMask;
			sourceMask= 0;
		}
		if(sourceMask_binary){
			delete sourceMask_binary;
			sourceMask_binary= 0;
		}
		if(sourcePlusCurvatureMask){
			delete sourcePlusCurvatureMask;
			sourcePlusCurvatureMask= 0;
		}	
	}//close if find nested sources
	


	return 0;

}//close Img::FindCompactSource()


TH1D* Img::GetPixelHisto(int nbins,bool normalize){

	if(!this->HasStats()){
		cerr<<"Img::GetPixelHisto(): WARN: No stats computed!"<<endl;
		return 0;
	}

	double Smin= fPixelStats->min;
	double Smax= fPixelStats->max;
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
	double sMax= histo->GetMaximum();
	double SMALL_NUMBER= 1.e-6;	
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

	cout<<"Img::FindValleyThreshold(): INFO: "<<npeaks<<" valleys detected!"<<endl;
	
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


double Img::FindOtsuThreshold(int nbins){
	
	//## Get histo and normalize
	TH1D* hist= this->GetPixelHisto(nbins,true);
	if(!hist) {
		cerr<<"Img::FindOtsuThreshold(): ERROR: Failed to compute pixel histo, return thr=0!"<<endl;
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


Img* Img::FindVLSLICSegmentation(std::vector<Region*>& regions,int RegionSize,double regularization,int minRegionSize,bool useLogContrast,bool mergeRegions,int algoType){

	//## Init stuff
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	
	cout<<"Img::FindSLICSegmentation(): INFO: RegionSize="<<RegionSize<<" regularization="<<regularization<<" minRegionSize="<<minRegionSize<<endl;

	//## Perform segmentation
	VLSlicSegmentation slic;
	slic.SetLogContrastMapping(useLogContrast);
  int status= slic.RunSegmentation(this,RegionSize,regularization,minRegionSize,mergeRegions,algoType);
	if(status<0){
		cerr<<"Img::FindSLICSegmentation(): ERROR: Segmentation failed!"<<endl;
		return 0;
	}
	
	//## Getting results
	regions= slic.GetRegions();//list of segmented regions
	Img* coloredImg= slic.GetClusterColoredImage(this,regions);
	
	return coloredImg;

}//close FindVLSLICSegmentation()



std::vector<Region*> Img::FindSegmentation(int RegionSize,double regularization, int minRegionSize,double eps, TGraph& contours,Img& coloredImg){

	//## Init stuff
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	int Npix= Nx*Ny;
	
	cout<<"Img::FindSLICSegmentation(): INFO: RegionSize="<<RegionSize<<" regularization="<<regularization<<" minRegionSize="<<minRegionSize<<endl;

	//## Perform segmentation
	VLSlicSegmentation slic;
	slic.SetLogContrastMapping(false);
  slic.RunSegmentation(this, RegionSize, regularization,minRegionSize,true,eps);

	//## Getting results
	std::vector< std::vector<long int> > pixelLabels= slic.GetPixelClusterIds();
	std::vector<Region*> clusters= slic.GetRegions();

	contours= *(slic.GetClusterContours());
	cout<<"Img::FindSLICSegmentation(): INFO: nClusterContours="<<contours.GetN()<<endl;
	
	coloredImg= *(slic.GetClusterColoredImage(this,clusters));
	
	return clusters;

}//close FindSegmentation()


Img* Img::FindCVSegmentation(double dt,double h,double lambda1,double lambda2,double mu,double nu,double p){

	cout<<"Img::FindCVSegmentation(): INFO: dt="<<dt<<" h="<<h<<" lambda1="<<lambda1<<" lambda2="<<lambda2<<" mu="<<mu<<" nu="<<nu<<" p="<<p<<endl;

	//## Perform segmentation
	ChanVeseSegmentation chanvese;
  int status= chanvese.RunSegmentation(this,dt,h,lambda1,lambda2,mu,nu,p);
	if(status<0){
		cerr<<"Img::FindCVSegmentation(): ERROR: ChanVese Segmentation failed!"<<endl;
		return 0;
	}

	//## Getting results
	Contour* cvContours= chanvese.GetContour();
	Img* segmentedImg= chanvese.GetSegmentedImage();	

	return segmentedImg;

}//close FindCVSegmentation()



int Img::FindExtendedSource(int RegionSize,double regularization,int minRegionSize,double mergeEps,double threshold){

	//## Check if background data are available
	if(!this->HasBkgData()){
		cout<<"Img::FindExtendedSource(): INFO: Image has no bkg information stored (hints: compute bkg first!), nothing to be done!"<<endl;
		return -1;
	}

	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	
	//## Get significance map
	Img* significanceMap= this->GetSignificanceMap(false);
	if(!significanceMap){
		cerr<<"Img::FindExtendedSource(): INFO: Cannot get significance map...nothing to be done!"<<endl;
		return -1;
	}


	//## Perform the image segmentation and retrieve segmented regions
	VLSlicSegmentation slic;
	slic.SetLogContrastMapping(false);
	int status= slic.RunSegmentation(this,RegionSize,regularization,minRegionSize,true,mergeEps);
	if(status<0){
		cerr<<"Img::FindExtendedSource(): ERROR: Image segmentation failed!"<<endl;
		return -1;
	}
	
	std::vector<Region*> regions= slic.GetRegions();
	cout<<"Img::FindExtendedSource(): INFO: "<<regions.size()<<" segmented regions found..."<<endl;

	
	//## Find background & significative regions
	bool includeSpatialPar= false;
	bool includeCurvPar= false;
	double CL= 0.975;
	Region* BkgRegion= slic.FindBackgroundRegion(CL,includeSpatialPar,includeCurvPar);
	if(!BkgRegion){
		cerr<<"Img::FindExtendedSource(): ERROR: Failed to compute bkg & significative regions!"<<endl;
		return -1;
	}
	

	
	//## Promote significative region to source if above threshold
	int nSources= 0;
	Source* aSource= 0;
	for(unsigned int k=0;k<regions.size();k++){
		Region* thisRegion= regions[k];
		bool isSignificative= thisRegion->fIsSignificative;
		if(!isSignificative) continue;

		//Create a new source from this region
		nSources++;
		std::string sourceName= std::string(Form("Img%s_sId%d",std::string(this->GetName()).c_str(),nSources));
		cout<<"Img::FindExtendedSource(): INFO: Adding new source (# "<<nSources<<") to list "<<endl;
		
		aSource= new Source;
		aSource->SetId(nSources);	
		aSource->SetName(sourceName);
		aSource->SetType(Source::eExtended);
		
		std::vector<Region::Pixel> thisPixelCollection= thisRegion->fPixelCollection;
		int nClusterPixels= (int)thisPixelCollection.size(); 
		for(int l=0;l<nClusterPixels;l++){
			int clusterPixelId= thisPixelCollection[l].id;
			int clusterPixelIdX, clusterPixelIdY, clusterPixelIdZ;
			this->GetBinXYZ(clusterPixelId,clusterPixelIdX,clusterPixelIdY,clusterPixelIdZ);
			double S= this->GetBinContent(clusterPixelId);			
			double Z= significanceMap->GetBinContent(clusterPixelId);
			double x= this->GetXaxis()->GetBinCenter(clusterPixelIdX);
			double y= this->GetYaxis()->GetBinCenter(clusterPixelIdY);
			
			Source::Pixel clusterPixel;
			clusterPixel.S= S;
			clusterPixel.Type= Source::eNormal;
			clusterPixel.Z= Z;
			clusterPixel.id= clusterPixelId;
			clusterPixel.ix= clusterPixelIdX-1;
			clusterPixel.iy= clusterPixelIdY-1;
			clusterPixel.x= x;//xwcs;
			clusterPixel.y= y;//ywcs;

			if( clusterPixelIdX<=1 || clusterPixelIdY<=1 || clusterPixelIdX>=Nx || clusterPixelIdY>=Ny) 
				aSource->SetEdgeFlag(true);

			aSource->AddPixel(clusterPixel);
		}//end loop cluster pixels

		//## Set flux correction factor
		double fluxCorrection= 1;
		if(this->HasMetaData()){
			double fx= this->metadata.Bmaj;
			double fy= this->metadata.Bmin;
			fluxCorrection= TMath::Pi()*fx*fy/(4*log(2));
			aSource->fFluxCorrection= fluxCorrection;
		}

		//## Compute stats parameters
		aSource->ComputeStats();

		//## Compute morphology parameters
		aSource->ComputeMorphologyParams();
	
		//## Do not consider 'line-like' sources (e.g. only 1 pixels wide in either x or y)
		if( !aSource->IsGoodSource() ) {
			cout<<"Img::FindExtendedSource(): INFO: Line-like source detected (only 1 pixels wide in x or y)...skip it!"<<endl;
			continue;
		}
		
		fHasSources= true;
		fSourceCollection.push_back(aSource);
	}//end loop regions

	return 0;

}//close FindExtendedSource()


Img* Img::GetLogNormalizedImage(int normmin,int normmax,bool skipEmptyBins){

	//Get first normalized image in range 1-256
	//Img* norm_img= GetNormalizedImage(1,256,skipEmptyBins);
	
	//int Nx= norm_img->GetNbinsX();
	//int Ny= norm_img->GetNbinsY();
	//double wmin= norm_img->GetMinimum();
	//double wmax= norm_img->GetMaximum();

	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	double wmin= this->GetMinimum();
	double wmax= this->GetMaximum();
	double safemin= 1;
	double safemax= 256;

	Img* norm_img= (Img*)this->Clone("norm_img");


	for (int i=0;i<Nx; i++) { 
		for (int j=0;j<Ny; j++) {	
			//double w= norm_img->GetBinContent(i+1,j+1);	
			double w= this->GetBinContent(i+1,j+1);	
			if(skipEmptyBins && w==0) continue;
				
			double w_norm= safemin + (safemax-safemin)*(w-wmin)/(wmax-wmin);

			//double w_log= normmin + log10(w/wmin)/log10(wmax/wmin) * (normmax-normmin);	
			double w_log= normmin + log10(w_norm/safemin)/log10(safemax/safemin) * (normmax-normmin);
			
			norm_img->SetBinContent(i+1,j+1,w_log);
		}
	}//end loop

	return norm_img;
	
}//close Img::GetLogNormalizedImage()



Img* Img::GetSigmoidNormalizedImage(int normmin,int normmax,double x0, double sigma, bool skipEmptyBins){

	//Get first normalized image in range 1-256
	
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	double wmin= this->GetMinimum();
	double wmax= this->GetMaximum();
	double safemin= 1;
	double safemax= 256;

	Img* norm_img= (Img*)this->Clone("normSigm_img");

	for (int i=0;i<Nx; i++) { 
		for (int j=0;j<Ny; j++) {	
			double w= this->GetBinContent(i+1,j+1);	
			if(skipEmptyBins && w==0) continue;
				
			double w_norm= safemin + (safemax-safemin)*(w-wmin)/(wmax-wmin);
			double w_sigm= normmin + (1./(1+exp(-(w_norm-x0)/sigma)))*(normmax-normmin);

			norm_img->SetBinContent(i+1,j+1,w_sigm);
		}
	}//end loop

	Img* norm2_img= norm_img->GetNormalizedImage(normmin,normmax,skipEmptyBins);

	return norm2_img;
	
}//close Img::GetSigmoidNormalizedImage()




Img* Img::GetNormalizedImage(int normmin,int normmax,bool skipEmptyBins){

	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	double zmin= this->GetMinimum();
	//double zmin= fPixelMin;
	double zmax= this->GetMaximum();
	//double zmax= fPixelMax;

	Img* norm_img= (Img*)this->Clone("norm_img");

	for (int i=0;i<Nx; i++) { 
		for (int j=0;j<Ny; j++) {	
			double w= this->GetBinContent(i+1,j+1);	
			if(skipEmptyBins && w==0) continue;
				
			double w_norm= normmin + (normmax-normmin)*(w-zmin)/(zmax-zmin);
			norm_img->SetBinContent(i+1,j+1,w_norm);
		}
	}//end loop

	return norm_img;

}//close Img::Normalize()

cv::Mat Img::ImgToMat(std::string encoding){

	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();

	/*
	double* pixels= 0;
	pixels= new double[Nx*Ny];
	int pixCounter= 0;

	for(int j=0;j<Ny;j++){
		for(int i=0;i<Nx;i++){
			float pixValue= this->GetBinContent(i+1,j+1);
			pixels[pixCounter]= pixValue;
			pixCounter++;
		}
	}

	//## Fill OpenCV mat image
	cv::Mat mat = cv::Mat(Nx, Ny, CV_32FC1, pixels);
	*/

	//## Fill OpenCV mat
	cv::Mat mat;
	if(encoding=="64") mat= cv::Mat::zeros(Ny,Nx,CV_64FC1);
	else if(encoding=="32") mat= cv::Mat::zeros(Ny,Nx,CV_32FC1);
	else{
		cerr<<"Img::ImgToMat(): WARN: Invalid encoding selected, using 64!"<<endl;
		mat= cv::Mat::zeros(Ny,Nx,CV_64FC1);
	}

	for(int j=0;j<Ny;j++){
		int rowId= Ny-1-j;
		for(int i=0;i<Nx;i++){
			int colId= i;
			mat.at<double>(rowId,colId)= this->GetBinContent(i+1,j+1);
		}//end loop x
	}//end loop y

	return mat;

}//close ImgToMat()

Img* Img::GetGuidedFilterImage(int radius,double eps){

	//## Init data
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	
	//Normalize img
	Img* img_norm= this->GetNormalizedImage(1,256);
	if(!img_norm) return 0;

	//## Convert image to OpenCV mat
	cv::Mat I= img_norm->ImgToMat("64");
	cv::Mat p= I;
	eps *= 255*255;   // Because the intensity range of our images is [0, 255]

	//## Run guided filter
	cout<<"Img::GetGuidedFilterImage(): Start filtering..."<<endl;
	//cv::Mat dst;	
	//cv::ximgproc::guidedFilter(I,p,dst,radius,eps,-1);
	cv::Mat dst = guidedFilter(I, p, radius, eps);

	//## Fill filtered image
	cout<<"Img::GetGuidedFilterImage(): Filter stage end..."<<endl;
	TString imgName= Form("%s_GuidedFilter",std::string(this->GetName()).c_str());
	Img* FilterImg= (Img*)this->Clone(imgName);
	FilterImg->Reset();
	
	for(int j=0;j<Ny;j++){
		int rowId= Ny-1-j;
		for(int i=0;i<Nx;i++){
			int colId= i;
			double matrixElement= dst.at<double>(rowId,colId);
			FilterImg->SetBinContent(i+1,j+1,matrixElement);
		}//end loop x
	}//end loop y

	if(img_norm) img_norm->Delete(); 

	return FilterImg;

}//close GetGuidedFilter()

Img* Img::GetBilateralFilterImage(int KernelSize, double sigmaColor, double sigmaSpace){

	//## Init data
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	
	//Normalize img
	Img* img_norm= this->GetNormalizedImage(0,255);

	//## Convert image to OpenCV mat
	cv::Mat src= img_norm->ImgToMat("32");

	//## Run bilateral filter
	cout<<"Img::GetBilateralFilterImage(): Start filtering..."<<endl;
	cv::Mat dst;
	cv::bilateralFilter(src, dst, KernelSize,sigmaColor,sigmaSpace,cv::BORDER_DEFAULT);
	
	//## Fill filtered image
	cout<<"Img::GetBilateralFilterImage(): Filter stage end..."<<endl;
	TString imgName= Form("%s_BilateralFilter",std::string(this->GetName()).c_str());
	Img* FilterImg= (Img*)this->Clone(imgName);
	FilterImg->Reset();
	
	for(int j=0;j<Ny;j++){
		int rowId= Ny-1-j;
		for(int i=0;i<Nx;i++){
			int colId= i;
			double matrixElement= (double)dst.at<float>(rowId,colId);
			FilterImg->SetBinContent(i+1,j+1,matrixElement);
		}//end loop x
	}//end loop y

	if(img_norm) img_norm->Delete(); 

	return FilterImg;

}//close GetBilateralFilterImage()

Img* Img::FindCannyEdges(std::vector<int>& edgePixelIds,int KernelSize,double lowThreshold,double thresholdRatio){

	//## Init data
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	
	//Normalize img
	Img* img_norm= this->GetNormalizedImage(0,255);

	//## Convert image to OpenCV mat
	cv::Mat mat= img_norm->ImgToMat();

	//## Convert it to gray
	mat.convertTo(mat, CV_8UC1);
  
	//## Apply the Canny edge detector
	//double lowThreshold= 200;//threshold for Canny edge
	//double thresholdRatio= 100;//threshold ratio
	cv::Canny( mat, mat, lowThreshold, lowThreshold*thresholdRatio, KernelSize, true);

  // Using Canny's output as a mask, we display our result
  //dst = Scalar::all(0);
  //src.copyTo( dst, detected_edges);

	Img* edges= (Img*)this->Clone("edges");
	edges->Reset();
	edgePixelIds.clear();

	for(int j=0;j<Ny;j++){
		int rowId= Ny-1-j;
		double y= this->GetYaxis()->GetBinCenter(j+1);
		for(int i=0;i<Nx;i++){
			double x= this->GetXaxis()->GetBinCenter(i+1);
			int colId= i;
			double mat_val= (double)mat.at<unsigned char>(rowId,colId);
			edges->SetBinContent(i+1,j+1,mat_val);			
		}//end loop x
	}//end loop y

	img_norm->Delete();
	return edges;

}//close FindCannyEdges()

Img* Img::FindEdges(std::vector<int>& edgePixelIds, std::string filter,int size,double scale){

	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	Img* img_filt= 0;
	if(filter=="LAPL") img_filt= this->GetLaplacianImage(false);
	else if(filter=="LoG") img_filt= this->GetLoGImage(false);
	else if(filter=="NormLoG") img_filt= this->GetNormLoGImage(size,scale,false);
	else if(filter=="NONE") img_filt= this;
	else{
		cerr<<"Img::FindEdges(): ERROR: Invalid filter specified!"<<endl;
		return 0;
	}
	Img* edges= (Img*)img_filt->Clone("edges");
	edges->Reset();
	edgePixelIds.clear();

	//Find zero crossings
	for (int s=0; s<Nx;++s) {
		int i= s+1;
		for (int t=0;t<Ny; ++t) {
			int j= t+1;
			int gBin= this->GetBin(i,j);
			double binContent= this->GetBinContent(i,j);
			if(binContent==0) continue;
			
			// Currently only checking interior pixels to avoid bounds checking
      if (s>0 && s<(Nx-1) && t>0 && t< (Ny-1)) {
				
				double w_ij= img_filt->GetBinContent(i,j);
				double w_iminus1_jminus1= img_filt->GetBinContent(i-1,j-1);
				double w_iplus1_jplus1= img_filt->GetBinContent(i+1,j+1);
				double w_iminus1_j= img_filt->GetBinContent(i-1,j);
				double w_iplus1_j= img_filt->GetBinContent(i+1,j);
				double w_iplus1_jminus1= img_filt->GetBinContent(i+1,j-1);
				double w_iminus1_jplus1= img_filt->GetBinContent(i-1,j+1);
				double w_i_jminus1= img_filt->GetBinContent(i,j-1);
				double w_i_jplus1= img_filt->GetBinContent(i,j+1);

				bool isEdge= false;

        if (w_ij==0) {
          if (0 != w_iminus1_jminus1
           || 0 != w_iminus1_j
           || 0 != w_iminus1_jplus1
           || 0 != w_i_jminus1
           || 0 != w_i_jplus1
           || 0 != w_iplus1_jminus1
           || 0 != w_iplus1_j
           || 0 != w_iplus1_jplus1)
           {
						 edges->SetBinContent(i,j,1);		
						 isEdge= true;
             //(**edges)(i,j) = fg;
           }
        }//close if
        else {
          if (abs(w_ij) < abs(w_iminus1_jminus1) && (w_ij>0) != (w_iminus1_jminus1>0)){
          	edges->SetBinContent(i,j,1);
						isEdge= true;
					}
     			else if (abs(w_ij) < abs(w_iminus1_j) && (w_ij>0) != (w_iminus1_j>0)){
						edges->SetBinContent(i,j,1);
						isEdge= true;		
					}
          else if (abs(w_ij) < abs(w_iminus1_jplus1) && (w_ij>0) != (w_iminus1_jplus1>0)){
						edges->SetBinContent(i,j,1);
						isEdge= true;
					}
          else if (abs(w_ij) < abs(w_i_jminus1) && (w_ij>0) != (w_i_jminus1>0)){
						edges->SetBinContent(i,j,1);		
						isEdge= true;
					}
          else if (abs(w_ij) < abs(w_i_jplus1) && (w_ij>0) != (w_i_jplus1>0)){
						 edges->SetBinContent(i,j,1);
             isEdge= true;
					}
     			else if (abs(w_ij) < abs(w_iplus1_jminus1) && (w_ij>0) != (w_iplus1_jminus1>0))	{
						 edges->SetBinContent(i,j,1);
						 isEdge= true;
          }
     			else if (abs(w_ij) < abs(w_iplus1_j) && (w_ij>0) != (w_iplus1_j>0)){
            edges->SetBinContent(i,j,1);
						isEdge= true;
					}
          else if (abs(w_ij) < abs(w_iplus1_jplus1) && (w_ij>0) != (w_iplus1_jplus1>0)){
				     edges->SetBinContent(i,j,1);  
						 isEdge= true;
					}
        }//close else
			
				if(isEdge) edgePixelIds.push_back(gBin);
      }//close if check boundary
    }//end loop bins y
  }//end loop bins x
	
	cout<<"Img::FindEdges(): INFO: "<<edgePixelIds.size()<<" edge pixels detected..."<<endl;

	return edges;

}//close FindEdges()


Img* Img::GetBinarized(double threshold,double bkgValue,double fgValue,bool isLowerThreshold){

	TString imgName= Form("%s_Binarized",std::string(this->GetName()).c_str());	
	Img* BinarizedImg= (Img*)this->Clone(imgName);
	BinarizedImg->Reset();

	for(int i=0;i<this->GetNbinsX();i++){
		for(int j=0;j<this->GetNbinsY();j++){
			double binContent= this->GetBinContent(i+1,j+1);
			if(binContent!=0){
				if(binContent>=threshold && !isLowerThreshold) BinarizedImg->SetBinContent(i+1,j+1,fgValue);
				else if(binContent<threshold && isLowerThreshold) BinarizedImg->SetBinContent(i+1,j+1,fgValue); 
			}
		}	
	}

	return BinarizedImg;

}//close Img::GetBinarized()

std::vector<Img*> Img::GetGradientImages(std::string KernelType){

	//## Compute gradients using GradientFilter helper class
	std::vector<Img*> gradImgList;
	Img* gradxImg= 0;
	Img* gradyImg= 0;

	if(KernelType=="SCHARR"){
		gradxImg= GradientFilter::GetFilteredImage(this,Img::eSCHARR_HORIZ);
		gradyImg= GradientFilter::GetFilteredImage(this,Img::eSCHARR_VERT);
	}
	else if(KernelType=="SOBEL"){
		gradxImg= GradientFilter::GetFilteredImage(this,Img::eSOBEL_HORIZ);
		gradyImg= GradientFilter::GetFilteredImage(this,Img::eSOBEL_VERT);
	}
	else {
		cerr<<"Img::GetGradientImages(): ERROR: Invalid kernel specified!"<<endl;
		return gradImgList;
	}

	gradImgList.push_back(gradxImg);
	gradImgList.push_back(gradyImg);

	return gradImgList;

}//close Img::GetGradientImages()


Img* Img::GetGradientImage(std::string KernelType){

	std::vector<Img*> gradImgs= this->GetGradientImages(KernelType);

	for(int i=0;i<this->GetNbinsX();i++){
		for(int j=0;j<this->GetNbinsY();j++){
			double gradx= gradImgs[0]->GetBinContent(i+1,j+1);
			double grady= gradImgs[1]->GetBinContent(i+1,j+1);
			double gradSum= sqrt(gradx*gradx+grady*grady);
			gradImgs[0]->SetBinContent(i+1,j+1,gradSum);
		}	
	}
	
	return gradImgs[0];

}//close Img::GetGradientImage()

Img* Img::GetKirschImage(){
	Img* kirschImg= GradientFilter::GetFilteredImage(this,Img::eKIRSCH);
	return kirschImg;
}//close Img::GetKirschImage()

Img* Img::GetLaplacianImage(bool invert){
	Img* laplImg= GradientFilter::GetFilteredImage(this,Img::eLAPLACE);
	if(invert) laplImg->Scale(-1);
	return laplImg;
}//close Img::GetLaplacianImage()


Img* Img::GetLaplacianWeightedImage(double f,std::string filtType,int size,double scale){

	Img* curvImg= 0;
	if(filtType=="LAPL") curvImg= this->GetLaplacianImage(true);
	else if(filtType=="LoG") curvImg= this->GetLoGImage(true);
	else if(filtType=="NormLoG") curvImg= this->GetNormLoGImage(size,scale,true);
	else{
		cerr<<"Img::GetLaplacianWeightedImage(): ERROR: Invalid filter selected!"<<endl;
		return 0;
	}

	if(!curvImg) return 0;
	
	Img* curvImg_norm= curvImg->GetNormalizedImage();
	Img* img_norm= this->GetNormalizedImage();


	//Apply weights
	for(int i=0;i<this->GetNbinsX();i++){
		for(int j=0;j<this->GetNbinsY();j++){
			double binContent= img_norm->GetBinContent(i+1,j+1);
			//double binContent= this->GetBinContent(i+1,j+1);
			double curv= curvImg_norm->GetBinContent(i+1,j+1);
			//double curv= curvImg->GetBinContent(i+1,j+1);
			double w= curv;
			//double binContent_weighted= binContent*w; 
			double binContent_weighted= (1-f)*binContent + f*curv; 
			curvImg->SetBinContent(i+1,j+1,binContent_weighted);
		}//end loop bins y
	}//end loop bins x

	curvImg_norm->Delete();
	img_norm->Delete();

	return curvImg;

	/*
	if(!this->HasStats()){
		cerr<<"Img::GetLaplacianWeightedImage(): ERROR: Image has no stats computed!"<<endl;
		return 0;
	}
	
	double I_min= this->GetMinimum();
	double I_mean= fPixelStats->mean;
	double I_mean_scaled= I_mean-I_min;
	Img* laplacianWeightedImg= this->GetLaplacianImage();
	if(!laplacianWeightedImg) return 0;

	//Apply weights
	for(int i=0;i<laplacianWeightedImg->GetNbinsX();i++){
		for(int j=0;j<laplacianWeightedImg->GetNbinsY();j++){
			double binContent= laplacianWeightedImg->GetBinContent(i+1,j+1);
			double Ixy= this->GetBinContent(i+1,j+1);
			double Ixy_scaled= binContent-I_min;
			double w= Ixy;
			if(Ixy_scaled>I_mean_scaled){
				w= I_mean_scaled + log(Ixy_scaled/I_mean_scaled);
			}
			double weightedBinContent= binContent*w; 
			laplacianWeightedImg->SetBinContent(i+1,j+1,weightedBinContent);
		}//end oop bins y
		
	}//end loop bins x
	
	return laplacianWeightedImg;
	*/
}//close GetLaplacianWeightedImage()

Img* Img::GetLoGImage(bool invert){
	Img* laplGausImg= GradientFilter::GetFilteredImage(this,Img::eLoG);
	if(laplGausImg && invert) laplGausImg->Scale(-1);
	return laplGausImg;
}//close Img::GetLoGImage()

Img* Img::GetNormLoGImage(int size,double scale,bool invert){
	Img* laplGausImg= GradientFilter::GetFilteredImage(this,Img::eNormLoG,size,scale);
	if(laplGausImg && invert) laplGausImg->Scale(-1);
	return laplGausImg;
}//close Img::GetNormLoGImage()


Img* Img::GetGaborFilteredImage(int KernelSize,int nThetaSteps, double sigma,double lambda,double gamma,double psi){

	//Get OpenCV image matrix
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	cv::Mat mat= this->ImgToMat();

	//Define kernels
	//sigma: standard deviation of the Gaussian function used in the filter
	//theta: orientation of the normal to the parallel stripes of the Gabor function
	//lambda: wavelength of the sinusoidal factor in the above equation
	//gamma: spatial aspect ratio
	//psi: phase offset
	//ktype: type and range of values that each pixel in the Gabor kernel can hold
	double theta= 0;
	double thetaStepSize= TMath::Pi()/nThetaSteps;
	int ktype= CV_64F;

	std::vector<cv::Mat> filteredImgs; 
	for(int k=0;k<nThetaSteps;k++){	
		//Define kernel given current orientation
		cv::Mat kernel = cv::getGaborKernel(cv::Size(KernelSize,KernelSize), sigma, theta, lambda, gamma, psi, ktype);
		//kernel/= 1.5*sum(kernel);
		
		//Get filtered image
		cv::Mat dest;
		cv::filter2D(mat, dest, CV_64F, kernel);
		filteredImgs.push_back(dest);
	
		theta+= thetaStepSize;
	}//end loop number of kernels

	//Compute filtered image (from the maximum of each pixels)
	TString imgName= Form("%s_GaborFilt",std::string(this->GetName()).c_str());	
	Img* GaborFilteredImg= (Img*)this->Clone(imgName);
	GaborFilteredImg->Reset();
	for(int j=0;j<Ny;j++){
		int rowId= Ny-1-j;
		for(int i=0;i<Nx;i++){
			int colId= i;
			
			double maxVal= -1.e+99;
			for(unsigned int k=0;k<filteredImgs.size();k++){
				double matrixElement= filteredImgs[k].at<double>(rowId,colId);
				if(matrixElement>maxVal) maxVal= matrixElement;
			}
			GaborFilteredImg->SetBinContent(i+1,j+1,maxVal);
		}//end loop x bins
	}//end loop y bins

	return GaborFilteredImg;

}//close GetGaborFilteredImage()


Img* Img::GetSkeleton(int KernelSize,double threshold,int maxIter){

	//## Init data
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	TEllipse* aCircle= 0;
	std::vector<TEllipse*> detectedCircles;
	detectedCircles.clear();

	//## Get normalized image
	Img* img_norm= this->GetNormalizedImage(0,255);

	//## Convert image to OpenCV mat
	cv::Mat img= img_norm->ImgToMat();

	//## Convert it to gray
	img.convertTo(img, CV_8UC1);

	//## Threshold image and define kernel
	cv::threshold(img, img, threshold, 255, cv::THRESH_BINARY); 
	cv::Mat skel(img.size(), CV_8UC1, cv::Scalar(0));
	cv::Mat temp;
	cv::Mat eroded;
	cv::Mat element = cv::getStructuringElement(cv::MORPH_CROSS, cv::Size(KernelSize,KernelSize));
 
	//## Compute image skeleton
	bool done;	
	int niter= 0;	
	do {
  	cv::erode(img, eroded, element);
  	cv::dilate(eroded, temp, element); // temp = open(img)
  	cv::subtract(img, temp, temp);
  	cv::bitwise_or(skel, temp, skel);
  	eroded.copyTo(img);
 
  	done = (cv::countNonZero(img) == 0);
		niter++;
		if(niter>=maxIter) break;
	} 	
	while (!done);
	cout<<"Img::GetSkeleton(): Completed in "<<niter<<" iterations..."<<endl;
	
	//Fill filtered image 
	TString imgName= Form("%s_Skeleton",std::string(this->GetName()).c_str());	
	Img* SkeletonImg= (Img*)this->Clone(imgName);
	SkeletonImg->Reset();
	for(int j=0;j<Ny;j++){
		int rowId= Ny-1-j;
		for(int i=0;i<Nx;i++){
			int colId= i;
			double binContent= this->GetBinContent(i+1,j+1);
			if(binContent==0) continue;
			double matrixElement= (double)(skel.at<unsigned char>(rowId,colId));
			SkeletonImg->SetBinContent(i+1,j+1,matrixElement);
		}//end loop x bins
	}//end loop y bins

	img_norm->Delete();

	return SkeletonImg;

}//close GetSkeleton()







Img* Img::GetSmoothedSaliencyMap(int reso,double regFactor,int minRegionSize,double sigmaS,double sigmaC,double saliencyKFactor,bool useRobust,bool addCurvDist){

	//## Compute segmentation in superpixels
	SLICSegmenter slic;
	int status= slic.RunSegmentation(this,reso,regFactor,minRegionSize,false);
	if(status<0){
		cerr<<"Img::GetSmoothedSaliencyMap(): ERROR: Superpixel segmentation failed!"<<endl;
		return 0;
	}
	std::vector<Region*> regions= slic.GetRegions(); 
	int nRegions= (int)regions.size();
	if(nRegions<=0) return 0;

	//## Compute saliency map
	Img* saliencyMap= GetSmoothedSaliencyMap(regions,sigmaS,sigmaC,saliencyKFactor,useRobust,addCurvDist);
	if(!saliencyMap){
		cerr<<"Img::GetSmoothedSaliencyMap(): ERROR: Failed to compute saliency map!"<<endl;
		return 0;
	}

	return saliencyMap;

}//close GetSmoothedSaliencyMap()

Img* Img::GetSmoothedSaliencyMap(std::vector<Region*> regions,double sigmaS,double sigmaC,double saliencyKFactor,bool useRobust,bool addCurvDist){

	//## Check regions
	int nRegions= (int)regions.size();
	if(nRegions<=0) return 0;
	
	//## Compute region parameters (including robust stats)
	int nPars= 2;
	int nPars_robust= 2;
	if(addCurvDist){
		nPars+= 2;
		nPars_robust+= 2;
	}
	int nPars_spatial= 2;
	TMatrixD params(nRegions,nPars);
	params.Zero();
	TMatrixD params_robust(nRegions,nPars_robust);
	params_robust.Zero();
	TMatrixD params_spatial(nRegions,nPars_spatial);
	params_spatial.Zero();

	cout<<"Img::GetSmoothedSaliencyMap(): INFO: Compute region pars..."<<endl;
	for(int i=0;i<nRegions;i++){
		//Compute pars
		//regions[i]->ComputeParameters(false,true,true);//SBAGLIATO
		regions[i]->ComputeParameters(false,true,false);//GIUSTO

		//Get params vectors
		Region::Parameters* regionPars= regions[i]->GetParams(addCurvDist);
		TVectorD* parVect= regionPars->pars;
		TVectorD* parVect_robust= regionPars->robustPars;
		TVectorD* parVect_spatial= regionPars->spatialPars;

		//Fill param matrix
		params[i]= *parVect;
		params_robust[i]= *parVect_robust;
		params_spatial[i]= *parVect_spatial;

		if(parVect) parVect->Delete();
		if(parVect_robust) parVect_robust->Delete();
		if(parVect_spatial) parVect_spatial->Delete();
		
	}//end loop regions

	//## Find min & max
	cout<<"Img::GetSmoothedSaliencyMap(): INFO: Compute param min/max..."<<endl;
	TVectorD params_min(nPars);
	TVectorD params_max(nPars);
	params_min.Zero();
	params_max.Zero();
	for(int j=0;j<nPars;j++){
		TVectorD v = TMatrixDColumn(params,j);
		params_min(j)= v.Min();
		params_max(j)= v.Max();
	}//end loop pars

	TVectorD params_robust_min(nPars_robust);
	TVectorD params_robust_max(nPars_robust);
	params_robust_min.Zero();
	params_robust_max.Zero();
	for(int j=0;j<nPars_robust;j++){
		TVectorD v = TMatrixDColumn(params_robust,j);
		params_robust_min(j)= v.Min();
		params_robust_max(j)= v.Max();
	}//end loop pars

	TVectorD params_spatial_min(nPars_spatial);
	TVectorD params_spatial_max(nPars_spatial);
	params_spatial_min.Zero();
	params_spatial_max.Zero();
	for(int j=0;j<nPars_spatial;j++){
		TVectorD v = TMatrixDColumn(params_spatial,j);
		params_spatial_min(j)= v.Min();
		params_spatial_max(j)= v.Max();
	}//end loop pars

	//## Normalize parameters to [0,1]	
	cout<<"Img::GetSmoothedSaliencyMap(): INFO: Normalize region pars.."<<endl;
	double NormMin= 0;
	double NormMax= 1;
	for(int i=0;i<nRegions;i++){
		for(int j=0;j<nPars;j++){
			double parValue= params(i,j); 
			double parValue_min= params_min(j);
			double parValue_max= params_max(j);
			double parValue_norm= NormMin + (NormMax-NormMin)*(parValue-parValue_min)/(parValue_max-parValue_min);
			params(i,j)= parValue;
		}//end loop pars
		for(int j=0;j<nPars_robust;j++){
			double parValue= params_robust(i,j); 
			double parValue_min= params_robust_min(j);
			double parValue_max= params_robust_max(j);
			double parValue_norm= NormMin + (NormMax-NormMin)*(parValue-parValue_min)/(parValue_max-parValue_min);
			params_robust(i,j)= parValue;
		}//end loop pars
		for(int j=0;j<nPars_spatial;j++){
			double parValue= params_spatial(i,j); 
			double parValue_min= params_spatial_min(j);
			double parValue_max= params_spatial_max(j);
			double parValue_norm= NormMin + (NormMax-NormMin)*(parValue-parValue_min)/(parValue_max-parValue_min);
			params_spatial(i,j)= parValue;
		}//end loop pars
	}//end loop regions


	//## Compute mutual region distances
	TMatrixD* DissimilarityMatrix= new TMatrixD(nRegions,nRegions);
	DissimilarityMatrix->Zero();
	TMatrixD* DissimilarityMatrix_robust= new TMatrixD(nRegions,nRegions);
	DissimilarityMatrix_robust->Zero();

	TMatrixD* SaliencySpatialWeightMatrix= new TMatrixD(nRegions,nRegions);
	SaliencySpatialWeightMatrix->Zero();//fill with 0s
	TMatrixD* SaliencyColorWeightMatrix= new TMatrixD(nRegions,nRegions);
	SaliencyColorWeightMatrix->Zero();//fill with 0s
	TMatrixD* SaliencyColorWeightMatrix_robust= new TMatrixD(nRegions,nRegions);
	SaliencyColorWeightMatrix_robust->Zero();//fill with 0s
	
	for(int i=0;i<nRegions;i++) {
		(*SaliencySpatialWeightMatrix)(i,i)= 1;
 		(*SaliencyColorWeightMatrix)(i,i)= 1;
		(*SaliencyColorWeightMatrix_robust)(i,i)= 1;
 	}
	
	for(int i=0;i<nRegions-1;i++){
		TVectorD x_i = TMatrixDRow(params,i);
		TVectorD xrobust_i = TMatrixDRow(params_robust,i);
		TVectorD xspatial_i = TMatrixDRow(params_spatial,i);
		for(int j=i+1;j<nRegions;j++){
			TVectorD x_j = TMatrixDRow(params,j);
			TVectorD xrobust_j = TMatrixDRow(params_robust,j);
			TVectorD xspatial_j = TMatrixDRow(params_spatial,j);
			double dist2= (x_i-x_j).Norm2Sqr();
			double dist2_robust= (xrobust_i-xrobust_j).Norm2Sqr();
			double dist2_spatial= (xspatial_i-xspatial_j).Norm2Sqr();
			double dissimilarity= sqrt(dist2)/(1+sqrt(dist2_spatial));
			double dissimilarity_robust= sqrt(dist2_robust)/(1+sqrt(dist2_spatial));

			double Dc= sqrt(dist2);
			double Dc_robust= sqrt(dist2_robust);
			double Ds= sqrt(dist2_spatial);

			double ws= exp(-Ds/(2.*sigmaS*sigmaS));
			double wc= exp(-Dc/(2.*sigmaC*sigmaC));
			double wc_robust= exp(-Dc_robust/(2.*sigmaC*sigmaC));

			(*SaliencySpatialWeightMatrix)(i,j)= ws;
			(*SaliencySpatialWeightMatrix)(j,i)= ws;
				
			(*SaliencyColorWeightMatrix)(i,j)= wc;
			(*SaliencyColorWeightMatrix)(j,i)= wc;

			(*SaliencyColorWeightMatrix_robust)(i,j)= wc_robust;
			(*SaliencyColorWeightMatrix_robust)(j,i)= wc_robust;

			(*DissimilarityMatrix)(i,j)= dissimilarity;
			(*DissimilarityMatrix)(j,i)= dissimilarity;
			(*DissimilarityMatrix_robust)(i,j)= dissimilarity_robust;
			(*DissimilarityMatrix_robust)(j,i)= dissimilarity_robust;
		}//end loop regions
	}//end loop regions


	//## Normalize similarity spatial and color Weights by rows so that sum_j(Wij)=1	
	TVectorD normWc(SaliencyColorWeightMatrix->GetNrows());
	TVectorD normWc_robust(SaliencyColorWeightMatrix_robust->GetNrows());
	TVectorD normWs(SaliencySpatialWeightMatrix->GetNrows());
	for(int i=0;i<SaliencySpatialWeightMatrix->GetNrows();i++){
		normWc(i)= ((TVectorD)TMatrixDRow(*SaliencyColorWeightMatrix,i)).Sum();
		normWc_robust(i)= ((TVectorD)TMatrixDRow(*SaliencyColorWeightMatrix_robust,i)).Sum();
		normWs(i)= ((TVectorD)TMatrixDRow(*SaliencySpatialWeightMatrix,i)).Sum();
	}//end loop rows
	SaliencySpatialWeightMatrix->NormByColumn(normWs);
	SaliencyColorWeightMatrix->NormByColumn(normWc);
	SaliencyColorWeightMatrix_robust->NormByColumn(normWc_robust);
	
	//## Compute mu (color-weighted mean position) matrix
	TMatrixD mu(nRegions,nPars_spatial);// regionsx2 (rows: regions, cols: x, y)
	mu= (*SaliencyColorWeightMatrix)*params_spatial;

	TMatrixD mu_robust(nRegions,nPars_spatial);// regionsx2 (rows: regions, cols: x, y)
	mu_robust= (*SaliencyColorWeightMatrix_robust)*params_spatial;
	
	//## Compute U & D terms
	std::vector<double> DVector;
	std::vector<double> UVector;
	std::vector<double> DVector_robust;
	std::vector<double> UVector_robust;
	
	for(int i=0;i<nRegions;i++){
		TVectorD x_i = TMatrixDRow(params,i);
		TVectorD xrobust_i = TMatrixDRow(params_robust,i);
		TVectorD mu_i= TMatrixDRow(mu,i);
		TVectorD mu_robust_i= TMatrixDRow(mu_robust,i);

		double salD= 0;
		double salU= 0;
		double saliency= 0;
		double salD_robust= 0;
		double salU_robust= 0;
		double saliency_robust= 0;
		for(int j=0;j<nRegions;j++) {		
			double ws= (*SaliencySpatialWeightMatrix)(i,j);//already normalized to unit sum
			double wc= (*SaliencyColorWeightMatrix)(i,j);//already normalized to unit sum
			double wc_robust= (*SaliencyColorWeightMatrix_robust)(i,j);//already normalized to unit sum

			//Compute D term
			TVectorD P_j= TMatrixDRow(params_spatial,j);
			double Pdist2= (P_j-mu_i).Norm2Sqr();	
			double Pdist= sqrt(Pdist2);
			double Pdist_weighted= Pdist*wc;
			salD+= Pdist_weighted;

			double Pdist2_robust= (P_j-mu_robust_i).Norm2Sqr();	
			double Pdist_robust= sqrt(Pdist2_robust);
			double Pdist_robust_weighted= Pdist_robust*wc_robust;
			salD_robust+= Pdist_robust_weighted;
				
			//Compute U term
			TVectorD x_j = TMatrixDRow(params,j);
			TVectorD xrobust_j = TMatrixDRow(params_robust,j);
			
			double dist2= (x_j-x_i).Norm2Sqr();
			double dist= sqrt(dist2);
			double u= dist*ws;
			salU+= u;

			double dist2_robust= (xrobust_j-xrobust_i).Norm2Sqr();
			double dist_robust= sqrt(dist2_robust);
			double u_robust= dist_robust*ws;
			salU_robust+= u_robust;
		}
		DVector.push_back(salD);	
		UVector.push_back(salU);
		DVector_robust.push_back(salD_robust);	
		UVector_robust.push_back(salU_robust);
	}//end loop rows

	//## Normalize D and U to [0,1]
	std::vector<double> saliencyU;
	saliencyU.assign(UVector.begin(),UVector.end());	
	std::vector<double> saliencyD;
	saliencyD.assign(DVector.begin(),DVector.end());
	std::vector<double> saliencyList;

	std::sort(DVector.begin(),DVector.end());
	std::sort(UVector.begin(),UVector.end());
	double minSalU= UVector[0];
	double maxSalU= UVector[UVector.size()-1];
	double minSalD= DVector[0];
	double maxSalD= DVector[DVector.size()-1];


	std::vector<double> saliencyU_robust;
	saliencyU_robust.assign(UVector_robust.begin(),UVector_robust.end());	
	std::vector<double> saliencyD_robust;
	saliencyD_robust.assign(DVector_robust.begin(),DVector_robust.end());
	std::vector<double> saliencyList_robust;

	std::sort(DVector_robust.begin(),DVector_robust.end());
	std::sort(UVector_robust.begin(),UVector_robust.end());
	double minSalU_robust= UVector_robust[0];
	double maxSalU_robust= UVector_robust[UVector_robust.size()-1];
	double minSalD_robust= DVector_robust[0];
	double maxSalD_robust= DVector_robust[DVector_robust.size()-1];

	
	TString imgName= Form("%s_saliency",this->GetName());
	Img* saliencyImg= (Img*)this->Clone(imgName);
	saliencyImg->SetNameTitle(imgName,imgName);
	saliencyImg->Reset();

	
	//## Compute saliency (METHOD2)
	double SMALL_NUMBER= 0.00000000001;
	for(unsigned int i=0;i<saliencyU.size();i++){
		double salU= saliencyU[i];
		double salU_norm= NormMin + (NormMax-NormMin)*(salU-minSalU)/(maxSalU-minSalU);
		double salD= saliencyD[i];
		double salD_norm= NormMin + (NormMax-NormMin)*(salD-minSalD)/(maxSalD-minSalD);
		salU_norm+= SMALL_NUMBER;
		salD_norm+= SMALL_NUMBER;
		saliencyU[i]= salU_norm;
		saliencyD[i]= salD_norm;
		double Saliency= salU_norm*exp(-saliencyKFactor*salD_norm);

		double salU_robust= saliencyU_robust[i];
		double salU_robust_norm= NormMin + (NormMax-NormMin)*(salU_robust-minSalU_robust)/(maxSalU_robust-minSalU_robust);
		double salD_robust= saliencyD_robust[i];
		double salD_robust_norm= NormMin + (NormMax-NormMin)*(salD_robust-minSalD_robust)/(maxSalD_robust-minSalD_robust);
		salU_robust_norm+= SMALL_NUMBER;
		salD_robust_norm+= SMALL_NUMBER;
		saliencyU_robust[i]= salU_robust_norm;
		saliencyD_robust[i]= salD_robust_norm;
		double Saliency_robust= salU_robust_norm*exp(-saliencyKFactor*salD_robust_norm);
		
		if(useRobust) regions[i]->fSaliency= Saliency_robust;
		else regions[i]->fSaliency= Saliency;
		
		//Fill image
		int nPixelsInRegion= (int)regions[i]->fNPix;
		for(int j=0;j<nPixelsInRegion;j++){//loop on pixels inside region
			int thisPixelId= (regions[i]->fPixelCollection)[j].id;
			saliencyImg->SetBinContent(thisPixelId,regions[i]->fSaliency);
		}//end loop pixels in region
	}//end loop regions
	

	if(DissimilarityMatrix) DissimilarityMatrix->Delete();
	if(DissimilarityMatrix_robust) DissimilarityMatrix_robust->Delete();

	
	return saliencyImg;

}//close GetSmoothedSaliencyMap()

Img* Img::GetSaliencyMap(int reso,double regFactor,int minRegionSize,double knnFactor,bool useRobustPars,double expFalloffPar,double distanceRegPar,bool interpolate){

	//## Compute segmentation in superpixels
	SLICSegmenter slic;
	int status= slic.RunSegmentation(this,reso,regFactor,minRegionSize,false);
	if(status<0){
		cerr<<"Img::GetSaliencyMap(): ERROR: Superpixel segmentation failed!"<<endl;
		return 0;
	}
	std::vector<Region*> regions= slic.GetRegions(); 
	int nRegions= (int)regions.size();
	if(nRegions<=0) return 0;

	//## Compute saliency map
	Img* saliencyMap= GetSaliencyMap(regions,knnFactor,useRobustPars,expFalloffPar,distanceRegPar,interpolate);
	if(!saliencyMap){
		cerr<<"Img::GetSaliencyMap(): ERROR: Failed to compute saliency map!"<<endl;
		return 0;
	}

	return saliencyMap;

}//close GetSaliencyMap()


Img* Img::GetSaliencyMap(std::vector<Region*>& regions,double knnFactor,bool useRobustPars,double expFalloffPar,double distanceRegPar,bool interpolate){

	//## Check regions
	int nRegions= static_cast<int>(regions.size());
	if(nRegions<=0) return 0;

	//## Compute number of neighbors to be used in saliency calculation
	int knn_min= 10;
	int knn_chosen= (int)(std::round(knnFactor*nRegions));//percentage of available regions
	int knn= std::max(knn_chosen,knn_min);

	int KNN= knn;
	if(knn>nRegions || knn<0) KNN= nRegions;

	//## Get image info
	double Xmin= this->GetXaxis()->GetXmin();
	double Xmax= this->GetXaxis()->GetXmax();
	double Ymin= this->GetYaxis()->GetXmin();
	double Ymax= this->GetYaxis()->GetXmax();
	double width= fabs(Xmax-Xmin);
	double height= fabs(Ymax-Ymin);
	double diagonal= sqrt(width*width + height*height);
	//double Xo= Xmin + width/2.;
 	//double Yo= Ymin + height/2.;

	//## Compute region pars
	cout<<"Img::GetSaliencyMap(): INFO: Compute region pars..."<<endl;
	for(int i=0;i<nRegions;i++) regions[i]->ComputeParameters(false,true,false);

	//## Compute dissimilarity matrix
	TMatrixD* ColorDistMatrix= new TMatrixD(nRegions,nRegions);
	ColorDistMatrix->Zero();

	TMatrixD* SpatialDistMatrix= new TMatrixD(nRegions,nRegions);
	SpatialDistMatrix->Zero();
	double dist_c_min= 1.e+99;
	double dist_c_max= -1.e+99;
	double dist_s_min= 1.e+99;
	double dist_s_max= -1.e+99;
	
	for(int i=0;i<nRegions-1;i++){

		//Region pars i-th
		double mu_i= regions[i]->fMean;
		double median_i= regions[i]->fMedian;
		double Xc_i= regions[i]->fX0;
		double Yc_i= regions[i]->fY0;
		
		for(int j=i+1;j<nRegions;j++){

			//Region pars j-th
			double n_j= (double)(regions[j]->fNPix);
			double mu_j= regions[j]->fMean;
			double median_j= regions[j]->fMedian;
			double Xc_j= regions[j]->fX0;
			double Yc_j= regions[j]->fY0;
					
			
			//Compute color & spatial distances
			double dist_c= fabs(mu_i-mu_j);
			if(useRobustPars) dist_c= fabs(median_i-median_j);
			double dist_s= sqrt( (Xc_i-Xc_j)*(Xc_i-Xc_j) + (Yc_i-Yc_j)*(Yc_i-Yc_j) ); 
			
			(*ColorDistMatrix)(i,j)= dist_c;
			(*ColorDistMatrix)(j,i)= dist_c;

			(*SpatialDistMatrix)(i,j)= dist_s;
			(*SpatialDistMatrix)(j,i)= dist_s;

			//Find min & max
			if(dist_c<dist_c_min) dist_c_min= dist_c;
			if(dist_c>dist_c_max) dist_c_max= dist_c;
			if(dist_s<dist_s_min) dist_s_min= dist_s;
			if(dist_s>dist_s_max) dist_s_max= dist_s;
	
		}//end loop regions
	}//end loop regions

	//## Normalize distances to [0,1]
	//## Color distances normalized to min & max
	//## Spatial distances normalized to image diagonal
	cout<<"Img::GetSaliencyMap(): INFO: Color dist min/max: "<<dist_c_min<<"/"<<dist_c_max<<", Spatial dist min/max: "<<dist_s_min<<"/"<<dist_s_max<<" img size("<<width<<" x "<<height<<" (diagonal="<<diagonal<<")"<<endl;
	
	double NormMin= 0;
	double NormMax= 1;

	for(int i=0;i<nRegions;i++){
		for(int j=i+1;j<nRegions;j++){
			
			double dist_c= (*ColorDistMatrix)(i,j);
			double dist_c_norm= NormMin + (NormMax-NormMin)*(dist_c-dist_c_min)/(dist_c_max-dist_c_min);
			double dist_s= (*SpatialDistMatrix)(i,j);
			double dist_s_norm= NormMin + (NormMax-NormMin)*(dist_s-dist_s_min)/(dist_s_max-dist_s_min);	
			//double dist_s_norm= dist_s/diagonal;

			(*ColorDistMatrix)(i,j)= dist_c_norm;
			(*ColorDistMatrix)(j,i)= dist_c_norm;

			(*SpatialDistMatrix)(i,j)= dist_s_norm;
			(*SpatialDistMatrix)(j,i)= dist_s_norm;
		}//end loop regions
	}//end loop regions
	
	cout<<"Img::GetSaliencyMap(): INFO: Color dist min/max: "<<ColorDistMatrix->Min()<<"/"<<ColorDistMatrix->Max()<<", Spatial dist min/max: "<<SpatialDistMatrix->Min()<<"/"<<SpatialDistMatrix->Max()<<endl;

	//## Create saliency image
	TString imgName= Form("%s_saliency",this->GetName());
	Img* saliencyImg= (Img*)this->Clone(imgName);
	saliencyImg->SetNameTitle(imgName,imgName);
	saliencyImg->Reset();

	//## Compute saliency 
	cout<<"Img::GetSaliencyMap(): INFO: Computing saliency ..."<<endl;	
	double Smin= 1.e+99;
	double Smax= -1.e+99;
	std::vector<double> SList;

	for(int i=0;i<nRegions;i++){

		std::vector<double> dissList;

		for(int j=0;j<nRegions;j++){
			//if(i==j) continue;
			double dist_c= (*ColorDistMatrix)(i,j);
			double dist_s= (*SpatialDistMatrix)(i,j);
			double dist= dist_c/(1 + distanceRegPar*dist_s);
			double dissimilarity= exp(-expFalloffPar*dist);

			dissList.push_back(dissimilarity);
		}//end loop regions

		//Sort color dissimilarities for region i-th to use only K-th neighbors in color
		std::vector<double> sorted;
		std::vector<size_t> sort_index;//sorting index
		Utils::sort(dissList,sorted,sort_index);	

		//Compute saliency over k-th neighbors
		double S= 0;
		for(int k=0;k<KNN;k++){
			size_t index= sort_index[k];
			double D= dissList[index];	
			S+= D;
		}
		S/= (double)(KNN);
		SList.push_back(S);

		
		if(S<Smin) Smin= S;
		if(S>Smax) Smax= S;

	}//end loop regions

	cout<<"Img::GetSaliencyMap(): INFO: Saliency min/max="<<Smin<<"/"<<Smax<<endl;	
	
	//## Delete matrix		
	if(ColorDistMatrix) ColorDistMatrix->Delete();
	if(SpatialDistMatrix) SpatialDistMatrix->Delete();

	
	//## Normalize saliency and fill maps
	TGraph2D* interpolationGraph= new TGraph2D;

	for(int i=0;i<nRegions;i++){
			
		//Normalize Saliency color
		double S= SList[i];
		double Snorm= NormMin + (NormMax-NormMin)*(S-Smin)/(Smax-Smin);
		double Saliency= 1.-Snorm;
	
		//Fill interpolation graph
		double Cx= regions[i]->fX0;	
		double Cy= regions[i]->fY0;	
		interpolationGraph->SetPoint(i,Cx,Cy,Saliency);

		//Fill image
		int nPixelsInRegion= (int)regions[i]->fNPix;
		for(int j=0;j<nPixelsInRegion;j++){//loop on pixels inside region
			int thisPixelId= (regions[i]->fPixelCollection)[j].id;
			saliencyImg->SetBinContent(thisPixelId,Saliency);
		}//end loop pixels in region
		
	}//end loop regions
	
	
	//## Interpolate?
	if(interpolate){
		cout<<"Img::GetSaliencyMap(): INFO: Interpolating saliency map..."<<endl;
		
		for(int i=0;i<saliencyImg->GetNbinsX();i++){	
			double x= saliencyImg->GetXaxis()->GetBinCenter(i+1);				
			for(int j=0;j<saliencyImg->GetNbinsY();j++){
				double y= saliencyImg->GetYaxis()->GetBinCenter(j+1);
				double interpSaliency= interpolationGraph->Interpolate(x,y);
				saliencyImg->SetBinContent(i+1,j+1,interpSaliency);
			}//end loop bins Y
		}//end loop binsX
	}//close if

	//## Delete allocated memory
	if(interpolationGraph) interpolationGraph->Delete();

	return saliencyImg;

}//close GetSaliencyMap()


Img* Img::GetSaliencyMap_MultiParVersion(int reso,double regFactor,int minRegionSize,double knnFactor,bool useRobustPars,bool addCurvDist,bool interpolate,double expFalloffPar,double distanceRegPar){

	//## Compute segmentation in superpixels
	SLICSegmenter slic;
	int status= slic.RunSegmentation(this,reso,regFactor,minRegionSize,false);
	if(status<0){
		cerr<<"Img::GetSaliencyMap(): ERROR: Superpixel segmentation failed!"<<endl;
		return 0;
	}
	std::vector<Region*> regions= slic.GetRegions(); 
	int nRegions= (int)regions.size();
	if(nRegions<=0) return 0;

	//## Compute saliency map
	Img* saliencyMap= GetSaliencyMap_MultiParVersion(regions,knnFactor,useRobustPars,addCurvDist,interpolate,expFalloffPar,distanceRegPar);
	if(!saliencyMap){
		cerr<<"Img::GetSaliencyMap(): ERROR: Failed to compute saliency map!"<<endl;
		return 0;
	}

	return saliencyMap;

}//close GetSaliencyMap()

Img* Img::GetSaliencyMap_MultiParVersion(std::vector<Region*>& regions,double knnFactor,bool useRobust,bool addCurvDist,bool interpolate,double expFalloffPar,double distanceRegPar){

	//## Check regions
	int nRegions= static_cast<int>(regions.size());
	if(nRegions<=0) return 0;

	//## Compute number of neighbors to be used in saliency calculation
	int knn_min= 10;
	int knn_chosen= (int)(std::round(knnFactor*nRegions));//percentage of available regions
	int knn= std::max(knn_chosen,knn_min);

	int KNN= knn;
	if(knn>nRegions || knn<0) KNN= nRegions;
	
	
	//## Compute region parameters (including robust stats)
	int nPars= 2;
	if(addCurvDist) nPars+= 2;
	int nPars_spatial= 2;

	TMatrixD params(nRegions,nPars);
	params.Zero();
	TMatrixD params_spatial(nRegions,nPars_spatial);
	params_spatial.Zero();

	cout<<"Img::GetSaliencyMap_MultiParVersion(): INFO: Compute region pars..."<<endl;
	for(int i=0;i<nRegions;i++){
		//Compute pars
		//regions[i]->ComputeParameters(false,true,true);//### SBAGLIATO
		regions[i]->ComputeParameters(false,true,false);//#### GIUSTO

		//regions[i]->Dump();

		//Get params vectors
		Region::Parameters* regionPars= regions[i]->GetParams(addCurvDist);
		TVectorD* parVect= regionPars->pars;
		TVectorD* parVect_robust= regionPars->robustPars;
		TVectorD* parVect_spatial= regionPars->spatialPars;

		//Fill param matrix
		if(useRobust) params[i]= *parVect_robust;	 
		else params[i]= *parVect;
		params_spatial[i]= *parVect_spatial;
		
		//Clear vectors
		if(parVect) parVect->Delete();
		if(parVect_robust) parVect_robust->Delete();
		if(parVect_spatial) parVect_spatial->Delete();
		
	}//end loop regions

	//## Find min & max
	cout<<"Img::GetSaliencyMap_MultiParVersion(): INFO: Compute parameter min/max..."<<endl;
	TVectorD params_min(nPars);
	TVectorD params_max(nPars);
	params_min.Zero();
	params_max.Zero();
	for(int j=0;j<nPars;j++){
		TVectorD v = TMatrixDColumn(params,j);
		params_min(j)= v.Min();
		params_max(j)= v.Max();
	}//end loop pars

	
	TVectorD params_spatial_min(nPars_spatial);
	TVectorD params_spatial_max(nPars_spatial);
	params_spatial_min.Zero();
	params_spatial_max.Zero();
	for(int j=0;j<nPars_spatial;j++){
		TVectorD v = TMatrixDColumn(params_spatial,j);
		params_spatial_min(j)= v.Min();
		params_spatial_max(j)= v.Max();
	}//end loop pars

	

	//## Normalize parameters to [0,1]	
	cout<<"Img::GetSaliencyMap(): INFO: Normalize region pars.."<<endl;
	double NormMin= 0;
	double NormMax= 1;
	for(int i=0;i<nRegions;i++){
		for(int j=0;j<nPars;j++){
			double parValue= params(i,j);
			double parValue_min= params_min(j);
			double parValue_max= params_max(j);
			double parValue_norm= NormMin + (NormMax-NormMin)*(parValue-parValue_min)/(parValue_max-parValue_min);
			//params(i,j)= parValue;//OLD VERSION
			params(i,j)= parValue_norm;

		}//end loop pars

		for(int j=0;j<nPars_spatial;j++){
			double parValue= params_spatial(i,j); 
			double parValue_min= params_spatial_min(j);
			double parValue_max= params_spatial_max(j);
			double parValue_norm= NormMin + (NormMax-NormMin)*(parValue-parValue_min)/(parValue_max-parValue_min);
			//params_spatial(i,j)= parValue;//OLD VERSION

			params_spatial(i,j)= parValue_norm;//DEBUG

		}//end loop pars

	}//end loop regions


	//## Compute mutual region distances
	TMatrixD* DissimilarityMatrix= new TMatrixD(nRegions,nRegions);
	DissimilarityMatrix->Zero();
	
	for(int i=0;i<nRegions-1;i++){
		TVectorD x_i = TMatrixDRow(params,i);
		TVectorD xspatial_i = TMatrixDRow(params_spatial,i);
		
		for(int j=i+1;j<nRegions;j++){
			TVectorD x_j = TMatrixDRow(params,j);
			TVectorD xspatial_j = TMatrixDRow(params_spatial,j);
			
			double dist2_c= (x_i-x_j).Norm2Sqr();
			double dist_c= sqrt(dist2_c);

			double dist2_spatial= (xspatial_i-xspatial_j).Norm2Sqr();
			double dist_spatial= sqrt(dist2_spatial);
			
			double dist= dist_c/(1.+ distanceRegPar*dist_spatial);
			double dissimilarity= exp(-expFalloffPar*dist);
			(*DissimilarityMatrix)(i,j)= dissimilarity;
			(*DissimilarityMatrix)(j,i)= dissimilarity;

		}//end loop regions
	}//end loop regions

	
	//## Compute saliency (METHOD1)
	double Smin= 1.e+99;
	double Smax= -1.e+99;
	std::vector<double> SList;
	
	for(int i=0;i<nRegions;i++){
		
		//Store dissimilarity in vector
		std::vector<double> dissList;
		
		for(int j=0;j<nRegions;j++){
			double D= (*DissimilarityMatrix)(i,j);
			dissList.push_back(D);
		}//end loop regions

		//Find the KNN neighbors
		std::vector<double> sorted;
		std::vector<size_t> sort_index;//sorting index
		Utils::sort( dissList,sorted,sort_index);	

		//Compute saliency over KNN neighbors
		double S= 0;
	
		for(int k=0;k<KNN;k++){
			size_t index= sort_index[k];
			double diss= dissList[index];
			S+= diss;
		}//end loop k-th neighbors 

		//Scale by KNN
		S/= (double)(KNN);

		if(S<Smin) Smin= S;
		if(S>Smax) Smax= S;
		SList.push_back(S);
		
	}//end loop regions
	
	//Delete matrix
	if(DissimilarityMatrix) DissimilarityMatrix->Delete();
	

	//## Create saliency map and interpolation graph
	TString imgName= Form("%s_saliency",this->GetName());
	Img* saliencyImg= (Img*)this->Clone(imgName);
	saliencyImg->SetNameTitle(imgName,imgName);
	saliencyImg->Reset();

	TGraph2D* interpolationGraph= new TGraph2D;
	

	//Normalize saliency term and fill image
	for(int i=0;i<nRegions;i++){
			
		//Normalize Saliency
		double S= SList[i];
		double Snorm= NormMin + (NormMax-NormMin)*(S-Smin)/(Smax-Smin);
		double Saliency= 1.-Snorm;
		
		//Fill interpolation graph
		double Cx= regions[i]->fX0;	
		double Cy= regions[i]->fY0;
		interpolationGraph->SetPoint(i,Cx,Cy,Saliency);

		//Fill image
		int nPixelsInRegion= (int)regions[i]->fNPix;
		for(int j=0;j<nPixelsInRegion;j++){//loop on pixels inside region
			int thisPixelId= (regions[i]->fPixelCollection)[j].id;
			saliencyImg->SetBinContent(thisPixelId,Saliency);
		}//end loop pixels in region
		
	}//end loop regions


	
	//## Interpolate?
	if(interpolate){
		cout<<"Img::GetSaliencyMap_MultiParVersion(): INFO: Interpolating saliency map..."<<endl;
		for(int i=0;i<saliencyImg->GetNbinsX();i++){	
			double x= saliencyImg->GetXaxis()->GetBinCenter(i+1);				
			for(int j=0;j<saliencyImg->GetNbinsY();j++){
				double y= saliencyImg->GetYaxis()->GetBinCenter(j+1);
				double interpSaliency= interpolationGraph->Interpolate(x,y);
				//if(interpSaliency!=0) 
					saliencyImg->SetBinContent(i+1,j+1,interpSaliency);
			}//end loop bins Y
		}//end loop binsX
	}//close if



	/*
	//## Return interpolated map
	Img* saliencyMap= 0;
	if(interpolate){
		cout<<"Img::GetSaliencyMap(): INFO: Interpolating saliency map..."<<endl;
		TMatrixD* saliencyMatrix= new TMatrixD(nRegions,3);
		for(int i=0;i<nRegions;i++){
			double Cx= regions[i]->fX0;	
			double Cy= regions[i]->fY0;	
			double Sal= regions[i]->fSaliency;
			(*saliencyMatrix)(i,0)= Cx;
			(*saliencyMatrix)(i,1)= Cy;
			(*saliencyMatrix)(i,2)= Sal;
		}//end loop regions
		Interpolator interpolator;
		interpolator.FindInterpolation(saliencyMatrix,this);
		saliencyMap= interpolator.GetInterpolatedMap();
		
		if(saliencyMatrix) saliencyMatrix->Delete();
		if(saliencyImg) saliencyImg->Delete();		
	}	//close if
	else{
		saliencyMap= saliencyImg;
	}
	*/

	//Clear memory
	if(interpolationGraph) interpolationGraph->Delete();

	return saliencyImg;

}//close GetSaliencyMap_MultiParVersion()



Img* Img::GetSaliencyMap_LuoMethod(int reso,double regFactor,int minRegionSize,double expFalloffPar,double distanceRegPar){

	//## Compute segmentation in superpixels
	SLICSegmenter slic;
	int status= slic.RunSegmentation(this,reso,regFactor,minRegionSize,false);
	if(status<0){
		cerr<<"Img::GetSaliencyMap_LuoMethod(): ERROR: Superpixel segmentation failed!"<<endl;
		return 0;
	}
	std::vector<Region*> regions= slic.GetRegions(); 
	int nRegions= (int)regions.size();
	if(nRegions<=0) return 0;

	//## Compute saliency map
	Img* saliencyMap= GetSaliencyMap_LuoMethod(regions,expFalloffPar,distanceRegPar);
	if(!saliencyMap){
		cerr<<"Img::GetSaliencyMap_LuoMethod(): ERROR: Failed to compute saliency map!"<<endl;
		return 0;
	}

	return saliencyMap;

}//close GetSaliencyMap_LuoMethod()

/*
Img* Img::GetSaliencyMap_LuoMethod(std::vector<Region*>& regions,double expFalloffPar){

	//## Check regions
	int nRegions= (int)regions.size();
	if(nRegions<=0) return 0;

	//## Compute image center
	double Xmin= this->GetXaxis()->GetXmin();
	double Xmax= this->GetXaxis()->GetXmax();
	double Ymin= this->GetYaxis()->GetXmin();
	double Ymax= this->GetYaxis()->GetXmax();
	double width= fabs(Xmax-Xmin);
	double height= fabs(Ymax-Ymin);
	double Xo= Xmin + width/2.;
 	double Yo= Ymin + height/2.;

	//## Compute region parameters (including robust stats)
	cout<<"Img::GetSaliencyMap_LuoMethod(): INFO: Compute region pars..."<<endl;

	std::vector<double> Ds_list;

	for(int i=0;i<nRegions;i++){
		//Compute pars
		regions[i]->ComputeParameters(false,true,false);

		//Get region pars
		int n_i= regions[i]->fNPix;
		double mu_i= regions[i]->fMean;
		double median_i= regions[i]->fMedian;
		double Xc_i= regions[i]->fX0;
		double Yc_i= regions[i]->fY0;

		//Compute spatial distance weight
		double Ds2_x= ((Xc_i-Xo)/width) * ((Xc_i-Xo)/width);
		double Ds2_y= ((Yc_i-Yo)/height) * ((Yc_i-Yo)/height);
		double Ds= sqrt(2.* (Ds2_x+Ds2_y) );//normalized to [0,1]

		Ds_list.push_back(Ds);

	}//end loop regions


	//## Compute color distances
	double EPS= 1.e-6;//a small number to prevent division by 0

	double S_color_min= 1.e+99;
	double S_color_max= -1.e+99;
	
	std::vector<double> S_color_list;
	std::vector<double> S_spatial_list;
	
	
	for(int i=0;i<nRegions;i++){
		int ix= i;
	
		//Region pars i-th
		int n_i= regions[i]->fNPix;
		double mu_i= regions[i]->fMean;
		double median_i= regions[i]->fMedian;
		double Xc_i= regions[i]->fX0;
		double Yc_i= regions[i]->fY0;

		//Accumulator
		double S_color= 0;
		double S_spatial= 0;
		double w_sum= 0;
		
		for(int j=0;j<nRegions;j++){
			if(i==j) continue;

			//Compute index
			int iy= j;
			if(j>i) {
				ix= j;
				iy= i;
			}
			int index= ix*(nRegions-1) - (ix-1)*ix/2 + iy - ix - 1;

			//Region pars j-th
			double n_j= (double)(regions[j]->fNPix);
			double mu_j= regions[j]->fMean;
			double median_j= regions[j]->fMedian;
			double Xc_j= regions[j]->fX0;
			double Yc_j= regions[j]->fY0;
					
			//Color distance
			double Dij_color= fabs(mu_i-mu_j);
			double Dij_color_scaled= Dij_color/n_j;
			
			//Compute color saliency term Sij= exp(-Dij/nj)
			double Sij_color= exp(-expFalloffPar*Dij_color);
			//double Sij_color= exp(-Dij_color_scaled);
			S_color+= Sij_color;
			
			//Compute spatial weights 
			double wij= 1./(Dij_color + EPS);
			w_sum+= wij;
			
			//Compute spatial saliency term Sij_s= wj x exp(1-Dj_s)
			double Dj_spatial= Ds_list[j];
			double Sij_spatial= wij * exp(1-Dj_spatial);
			S_spatial+= Sij_spatial;
			
		}//end loop regions

		//Fill list
		S_color_list.push_back(S_color);
		
		//Compute color saliency min & max
		if(S_color<S_color_min) S_color_min= S_color;
		if(S_color>S_color_max) S_color_max= S_color;
		
		//Normalize spatial saliency
		S_spatial/= w_sum;
		S_spatial_list.push_back(S_spatial);
		
	}//end loop regions

	cout<<"Img::GetSaliencyMap_LuoMethod(): INFO: Saliency min/max: "<<S_color_min<<"/"<<S_color_min<<endl;
	
	//Create saliency images
	TString imgName= Form("%s_saliency",this->GetName());
	Img* saliencyImg= (Img*)this->Clone(imgName);
	saliencyImg->SetNameTitle(imgName,imgName);
	saliencyImg->Reset();

	imgName= Form("%s_saliency_color",this->GetName());
	Img* saliencyImg_color= (Img*)this->Clone(imgName);
	saliencyImg_color->SetNameTitle(imgName,imgName);
	saliencyImg_color->Reset();

	imgName= Form("%s_saliency_spatial",this->GetName());
	Img* saliencyImg_spatial= (Img*)this->Clone(imgName);
	saliencyImg_spatial->SetNameTitle(imgName,imgName);
	saliencyImg_spatial->Reset();

	//Normalize color saliency and fill image
	double NormMin= 0;
	double NormMax= 1;
	for(int i=0;i<nRegions;i++){
		//Compute color saliency normalized to [0,1]
		double S_color= S_color_list[i];	
		double S_color_norm= NormMin + (NormMax-NormMin)*(S_color-S_color_min)/(S_color_max-S_color_min);
		double Saliency_color= 1-S_color_norm;

		//Get spatial saliency 
		double Saliency_spatial= S_spatial_list[i];
		
		//Compute global saliency
		double Saliency= Saliency_color*Saliency_spatial;
		
		//Fill image
		int nPixelsInRegion= (int)regions[i]->fNPix;

		
		for(int j=0;j<nPixelsInRegion;j++){//loop on pixels inside region
			int thisPixelId= (regions[i]->fPixelCollection)[j].id;
			saliencyImg_color->SetBinContent(thisPixelId,Saliency_color);
			saliencyImg_spatial->SetBinContent(thisPixelId,Saliency_spatial);
			saliencyImg->SetBinContent(thisPixelId,Saliency);
		}//end loop pixels in region
		
	}//end loop regions

	//---- DEBUG ----------------
	//Plot saliency maps
	TCanvas* ColorSaliencyPlot= new TCanvas("ColorSaliencyPlot","ColorSaliencyPlot");
	ColorSaliencyPlot->cd();
	saliencyImg_color->Draw("COLZ");

	TCanvas* SpatialSaliencyPlot= new TCanvas("SpatialSaliencyPlot","SpatialSaliencyPlot");
	SpatialSaliencyPlot->cd();
	saliencyImg_spatial->Draw("COLZ");
	//-----------------------------	

	return saliencyImg;

}//close GetSaliencyMap_LuoMethod()
*/

Img* Img::GetSaliencyMap_LuoMethod(std::vector<Region*>& regions,double expFalloffPar,double distanceRegPar){

	//## Check regions
	int nRegions= (int)regions.size();
	if(nRegions<=0) return 0;

	//## Compute image center
	double Xmin= this->GetXaxis()->GetXmin();
	double Xmax= this->GetXaxis()->GetXmax();
	double Ymin= this->GetYaxis()->GetXmin();
	double Ymax= this->GetYaxis()->GetXmax();
	double width= fabs(Xmax-Xmin);
	double height= fabs(Ymax-Ymin);
	double diagonal= sqrt(width*width + height*height);
	double Xo= Xmin + width/2.;
 	double Yo= Ymin + height/2.;

	//## Compute region parameters (including robust stats)
	cout<<"Img::GetSaliencyMap_LuoMethod(): INFO: Compute region pars..."<<endl;
	for(int i=0;i<nRegions;i++) regions[i]->ComputeParameters(false,true,false);

	//## Compute dissimilarity matrix
	TMatrixD* ColorDistMatrix= new TMatrixD(nRegions,nRegions);
	ColorDistMatrix->Zero();

	TMatrixD* SpatialDistMatrix= new TMatrixD(nRegions,nRegions);
	SpatialDistMatrix->Zero();

	double dist_c_min= 1.e+99;
	double dist_c_max= -1.e+99;
	double dist_s_min= 1.e+99;
	double dist_s_max= -1.e+99;
	
	
	for(int i=0;i<nRegions-1;i++){

		//Region pars i-th
		int n_i= regions[i]->fNPix;
		double mu_i= regions[i]->fMean;
		double Xc_i= regions[i]->fX0;
		double Yc_i= regions[i]->fY0;
		
		//Scale by number of pixels (DEBUG)
		//Xc_i/= (double)(n_i);
		//Yc_i/= (double)(n_i);

		for(int j=i+1;j<nRegions;j++){

			//Region pars j-th
			double n_j= (double)(regions[j]->fNPix);
			double mu_j= regions[j]->fMean;
			double median_j= regions[j]->fMedian;
			double Xc_j= regions[j]->fX0;
			double Yc_j= regions[j]->fY0;
					
			//Scale by number of pixels (DEBUG)
			//Xc_j/= (double)(n_j);
			//Yc_j/= (double)(n_j);

			//Compute color & spatial distances
			double dist_c= fabs(mu_i-mu_j);
			double dist_s= sqrt( (Xc_i-Xc_j)*(Xc_i-Xc_j) + (Yc_i-Yc_j)*(Yc_i-Yc_j) ); 
			
			(*ColorDistMatrix)(i,j)= dist_c;
			(*ColorDistMatrix)(j,i)= dist_c;

			(*SpatialDistMatrix)(i,j)= dist_s;
			(*SpatialDistMatrix)(j,i)= dist_s;

			//Find min & max
			if(dist_c<dist_c_min) dist_c_min= dist_c;
			if(dist_c>dist_c_max) dist_c_max= dist_c;
			if(dist_s<dist_s_min) dist_s_min= dist_s;
			if(dist_s>dist_s_max) dist_s_max= dist_s;
	
		}//end loop regions
	}//end loop regions

	//## Normalize distances to [0,1]
	//## Color distances normalized to min & max
	//## Spatial distances normalized to image diagonal
	cout<<"Img::GetSaliencyMap_LuoMethod(): INFO: Color dist min/max: "<<dist_c_min<<"/"<<dist_c_max<<", Spatial dist min/max: "<<dist_s_min<<"/"<<dist_s_max<<" img size("<<width<<" x "<<height<<" (diagonal="<<diagonal<<")"<<endl;
	
	
	double NormMin= 0;
	double NormMax= 1;

	
	for(int i=0;i<nRegions;i++){
		for(int j=i+1;j<nRegions;j++){
			
			double dist_c= (*ColorDistMatrix)(i,j);
			double dist_c_norm= NormMin + (NormMax-NormMin)*(dist_c-dist_c_min)/(dist_c_max-dist_c_min);
			double dist_s= (*SpatialDistMatrix)(i,j);
			double dist_s_norm= NormMin + (NormMax-NormMin)*(dist_s-dist_s_min)/(dist_s_max-dist_s_min);	
			//double dist_s_norm= dist_s/diagonal;

			(*ColorDistMatrix)(i,j)= dist_c_norm;
			(*ColorDistMatrix)(j,i)= dist_c_norm;

			(*SpatialDistMatrix)(i,j)= dist_s_norm;
			(*SpatialDistMatrix)(j,i)= dist_s_norm;
		}//end loop regions
	}//end loop regions
	

	cout<<"Img::GetSaliencyMap_LuoMethod(): INFO: Color dist min/max: "<<ColorDistMatrix->Min()<<"/"<<ColorDistMatrix->Max()<<", Spatial dist min/max: "<<SpatialDistMatrix->Min()<<"/"<<SpatialDistMatrix->Max()<<endl;

	//## Create saliency image
	TString imgName= Form("%s_saliency",this->GetName());
	Img* saliencyImg= (Img*)this->Clone(imgName);
	saliencyImg->SetNameTitle(imgName,imgName);
	saliencyImg->Reset();

	imgName= Form("%s_saliency_color",this->GetName());
	Img* saliencyImg_color= (Img*)this->Clone(imgName);
	saliencyImg_color->SetNameTitle(imgName,imgName);
	saliencyImg_color->Reset();

	imgName= Form("%s_saliency_spatial",this->GetName());
	Img* saliencyImg_spatial= (Img*)this->Clone(imgName);
	saliencyImg_spatial->SetNameTitle(imgName,imgName);
	saliencyImg_spatial->Reset();

	//## Compute saliency 
	cout<<"Img::GetSaliencyMap_LuoMethod(): INFO: Computing saliency ..."<<endl;	
	double Smin= 1.e+99;
	double Smax= -1.e+99;
	double Smin_spatial= 1.e+99;
	double Smax_spatial= -1.e+99;
	double Vmin= 1.e+99;
	double Vmax= -1.e+99;
	std::vector<double> SList;
	std::vector<double> SList_spatial;
	std::vector<double> VList;
	double EPS= 1.e-6;//a small number to prevent division by 0
	int KNN= 200;

	for(int i=0;i<nRegions;i++){

		int n_i= regions[i]->fNPix;
		double S= 0;
		double W= 0;
		double S_spatial= 0;
		double X0_s= 0;
		double Y0_s= 0;
		std::vector<double> dissColorList;
		std::vector<double> dissSpatialList;

		for(int j=0;j<nRegions;j++){
			if(i==j) continue;
			double Xc_j= regions[j]->fX0;
			double Yc_j= regions[j]->fY0;
			double dist_c= (*ColorDistMatrix)(i,j);
			double dist_s= (*SpatialDistMatrix)(i,j);
			double dist= dist_c/(1 + distanceRegPar*dist_s);
			//double Sij_c= exp(-expFalloffPar*dist_c);
			double Sij_c= exp(-expFalloffPar*dist);

			//double wij= 1./(dist_c+EPS);
			double wij= Sij_c;
			double Sij_s= wij*dist_s*dist_s;//variance

			if(i!=j){
				W+= wij; 
				S+= Sij_c;
				S_spatial+= Sij_s;

				X0_s+= Xc_j*wij/width;
				Y0_s+= Yc_j*wij/height;
			}

			dissColorList.push_back(dist_c);
 			dissSpatialList.push_back(dist_s);

		}//end loop regions

		S_spatial/= W;
		X0_s/= W;
		Y0_s/= W;

		//Sort color dissimilarities for region i-th to use only K-th neighbors in color
		std::vector<double> sorted;
		std::vector<size_t> sort_index;//sorting index
		Utils::sort(dissColorList,sorted,sort_index);	

		//Compute spatial saliency over k-th neighbors
		double V= 0;

		for(int k=0;k<KNN;k++){
			size_t index= sort_index[k];
			double dist_s= dissSpatialList[index];	
			int n_k= regions[index]->fNPix;	
			double Vik= (dist_s*dist_s);
			V+= Vik;
		}
		V/= (double)KNN;
		VList.push_back(V);

		if(V<Vmin) Vmin= V;
		if(V>Vmax) Vmax= V;

		/*
		//Calculate spatial variance (Perazzi)
		double V= 0;
		for(int j=0;j<nRegions;j++){
			double Xc_j= regions[j]->fX0;
			double Yc_j= regions[j]->fY0;
			double dist_c= (*ColorDistMatrix)(i,j);
			double dist_s= (*SpatialDistMatrix)(i,j);
			double Sij_c= exp(-expFalloffPar*dist_c);

			//double wij= 1./(dist_c+EPS);
			double wij= Sij_c;
			double Vij= (Xc_j/width-X0_s)*(Xc_j/width-X0_s) + (Yc_j/height-Y0_s)*(Yc_j/height-Y0_s);
			V+= wij* Vij;
		}//end loop regions
		
		if(V<Vmin) Vmin= V;
		if(V>Vmax) Vmax= V;
		*/

		if(S<Smin) Smin= S;
		if(S>Smax) Smax= S;
		if(S_spatial<Smin_spatial) Smin_spatial= S_spatial;
		if(S_spatial>Smax_spatial) Smax_spatial= S_spatial;
		SList.push_back(S);
		SList_spatial.push_back(S_spatial);
		//VList.push_back(V);

	}//end loop regions
		
	cout<<"Img::GetSaliencyMap_LuoMethod(): INFO: Saliency min/max="<<Smin<<"/"<<Smax<<", Saliency spatial min/max="<<Smin_spatial<<"/"<<Smax_spatial<<endl;	
	
	for(int i=0;i<nRegions;i++){
			
		//Normalize Saliency color
		double S_color= SList[i];
		double Snorm_color= NormMin + (NormMax-NormMin)*(S_color-Smin)/(Smax-Smin);
		double Saliency_color= 1.-Snorm_color;

		double S_spatial= SList_spatial[i];
		double Snorm_spatial= NormMin + (NormMax-NormMin)*(S_spatial-Smin_spatial)/(Smax_spatial-Smin_spatial);
		//double Saliency_spatial= 1.-Snorm_spatial;

		
		double V= VList[i];
		double Vnorm= NormMin + (NormMax-NormMin)*(V-Vmin)/(Vmax-Vmin);
		double Saliency_spatial= 1.-Vnorm;
		

		double Saliency= Saliency_color*Saliency_spatial;

		//Fill image
		int nPixelsInRegion= (int)regions[i]->fNPix;
		for(int j=0;j<nPixelsInRegion;j++){//loop on pixels inside region
			int thisPixelId= (regions[i]->fPixelCollection)[j].id;
			saliencyImg->SetBinContent(thisPixelId,Saliency);
			saliencyImg_color->SetBinContent(thisPixelId,Saliency_color);
			saliencyImg_spatial->SetBinContent(thisPixelId,Saliency_spatial);
		}//end loop pixels in region
		
	}//end loop regions
	
	//---- DEBUG ----------------
	//Plot saliency maps
	TCanvas* ColorSaliencyPlot= new TCanvas("ColorSaliencyPlot","ColorSaliencyPlot");
	ColorSaliencyPlot->cd();
	saliencyImg_color->Draw("COLZ");

	TCanvas* SpatialSaliencyPlot= new TCanvas("SpatialSaliencyPlot","SpatialSaliencyPlot");
	SpatialSaliencyPlot->cd();
	saliencyImg_spatial->Draw("COLZ");
	//-----------------------------	

	return saliencyImg;

}//close GetSaliencyMap_LuoMethod()



Img* Img::GetSaliencyMap_GofermanMethod(int reso,double regFactor,int minRegionSize,double knnFactor,double distanceRegPar,double expFalloffPar){

	//## Compute segmentation in superpixels
	SLICSegmenter slic;
	int status= slic.RunSegmentation(this,reso,regFactor,minRegionSize,false);
	if(status<0){
		cerr<<"Img::GetSaliencyMap_GofermanMethod(): ERROR: Superpixel segmentation failed!"<<endl;
		return 0;
	}
	std::vector<Region*> regions= slic.GetRegions(); 
	int nRegions= (int)regions.size();
	if(nRegions<=0) return 0;

	//## Compute saliency map
	Img* saliencyMap= GetSaliencyMap_GofermanMethod(regions,knnFactor,distanceRegPar,expFalloffPar);
	if(!saliencyMap){
		cerr<<"Img::GetSaliencyMap_GofermanMethod(): ERROR: Failed to compute saliency map!"<<endl;
		return 0;
	}

	return saliencyMap;

}//close GetSaliencyMap_GofermanMethod()


Img* Img::GetSaliencyMap_GofermanMethod(std::vector<Region*>& regions,double knnFactor,double distanceRegPar,double expFalloffPar){

	//## Check regions
	int nRegions= (int)regions.size();
	if(nRegions<=0) return 0;

	//## Compute number of neightbors to be used in saliency estimation
	int knn_min= 10;
	int knn_chosen= (int)(std::round(knnFactor*nRegions));//percentage of available regions
	int knn= std::max(knn_chosen,knn_min);
	int KNN= knn;
	if(knn>nRegions || knn<0) KNN= nRegions;

	//## Compute image center
	double Xmin= this->GetXaxis()->GetXmin();
	double Xmax= this->GetXaxis()->GetXmax();
	double Ymin= this->GetYaxis()->GetXmin();
	double Ymax= this->GetYaxis()->GetXmax();
	double width= fabs(Xmax-Xmin);
	double height= fabs(Ymax-Ymin);
	double diagonal= sqrt(width*width + height*height);
	double Xo= Xmin + width/2.;
 	double Yo= Ymin + height/2.;

	//## Compute region parameters (including robust stats)
	cout<<"Img::GetSaliencyMap_GofermanMethod(): INFO: Compute region pars..."<<endl;
	for(int i=0;i<nRegions;i++) regions[i]->ComputeParameters(false,true,false);

	//## Compute dissimilarity matrix
	TMatrixD* ColorDistMatrix= new TMatrixD(nRegions,nRegions);
	ColorDistMatrix->Zero();

	TMatrixD* SpatialDistMatrix= new TMatrixD(nRegions,nRegions);
	SpatialDistMatrix->Zero();

	TMatrixD* DissimilarityMatrix= new TMatrixD(nRegions,nRegions);
	DissimilarityMatrix->Zero();

	double dist_c_min= 1.e+99;
	double dist_c_max= -1.e+99;
	double dist_s_min= 1.e+99;
	double dist_s_max= -1.e+99;
	
	for(int i=0;i<nRegions-1;i++){

		//Region pars i-th
		int n_i= regions[i]->fNPix;
		double mu_i= regions[i]->fMean;
		double Xc_i= regions[i]->fX0;
		double Yc_i= regions[i]->fY0;
		
		//Scale by number of pixels
		//Xc_i/= (double)(n_i);
		//Yc_i/= (double)(n_i);

		for(int j=i+1;j<nRegions;j++){

			//Region pars j-th
			double n_j= (double)(regions[j]->fNPix);
			double mu_j= regions[j]->fMean;
			double median_j= regions[j]->fMedian;
			double Xc_j= regions[j]->fX0;
			double Yc_j= regions[j]->fY0;
					
			//Scale by number of pixels
			//Xc_j/= (double)(n_j);
			//Yc_j/= (double)(n_j);

			//Compute color & spatial distances
			double dist_c= fabs(mu_i-mu_j);
			double dist_s= sqrt( (Xc_i-Xc_j)*(Xc_i-Xc_j) + (Yc_i-Yc_j)*(Yc_i-Yc_j) ); 

			(*ColorDistMatrix)(i,j)= dist_c;
			(*ColorDistMatrix)(j,i)= dist_c;

			(*SpatialDistMatrix)(i,j)= dist_s;
			(*SpatialDistMatrix)(j,i)= dist_s;

			//Find min & max
			if(dist_c<dist_c_min) dist_c_min= dist_c;
			if(dist_c>dist_c_max) dist_c_max= dist_c;
			if(dist_s<dist_s_min) dist_s_min= dist_s;
			if(dist_s>dist_s_max) dist_s_max= dist_s;
	
		}//end loop regions
	}//end loop regions

	//## Normalize distances to [0,1]
	//## Color distances normalized to min & max
	//## Spatial distances normalized to image diagonal
	cout<<"Img::GetSaliencyMap_GofermanMethod(): INFO: Color dist min/max: "<<dist_c_min<<"/"<<dist_c_max<<", Spatial dist min/max: "<<dist_s_min<<"/"<<dist_s_max<<" img size("<<width<<" x "<<height<<" (diagonal="<<diagonal<<")"<<endl;
	
	
	double NormMin= 0;
	double NormMax= 1;

	
	for(int i=0;i<nRegions;i++){
		for(int j=i+1;j<nRegions;j++){
			
			double dist_c= (*ColorDistMatrix)(i,j);
			double dist_c_norm= NormMin + (NormMax-NormMin)*(dist_c-dist_c_min)/(dist_c_max-dist_c_min);
			double dist_s= (*SpatialDistMatrix)(i,j);
			double dist_s_norm= NormMin + (NormMax-NormMin)*(dist_s-dist_s_min)/(dist_s_max-dist_s_min);	
			//double dist_s_norm= dist_s/diagonal;

			(*ColorDistMatrix)(i,j)= dist_c_norm;
			(*ColorDistMatrix)(j,i)= dist_c_norm;

			(*SpatialDistMatrix)(i,j)= dist_s_norm;
			(*SpatialDistMatrix)(j,i)= dist_s_norm;
		}//end loop regions
	}//end loop regions
	

	cout<<"Img::GetSaliencyMap_GofermanMethod(): INFO: Color dist min/max: "<<ColorDistMatrix->Min()<<"/"<<ColorDistMatrix->Max()<<", Spatial dist min/max: "<<SpatialDistMatrix->Min()<<"/"<<SpatialDistMatrix->Max()<<endl;


	//## Create saliency image
	TString imgName= Form("%s_saliency",this->GetName());
	Img* saliencyImg= (Img*)this->Clone(imgName);
	saliencyImg->SetNameTitle(imgName,imgName);
	saliencyImg->Reset();

	//## Compute saliency 
	cout<<"Img::GetSaliencyMap_GofermanMethod(): INFO: Computing saliency (KNN="<<KNN<<" regions, distanceRegPar="<<distanceRegPar<<")"<<endl;	
	double Smin= 1.e+99;
	double Smax= -1.e+99;
	std::vector<double> SList;

	for(int i=0;i<nRegions;i++){
		
		std::vector<double> dissColorList;
		std::vector<double> dissList;
		for(int j=0;j<nRegions;j++){
			//if(i==j) continue;
			double dist_c= (*ColorDistMatrix)(i,j);
			double dist_s= (*SpatialDistMatrix)(i,j);
			double dist= dist_c/(1. + distanceRegPar*dist_s);
			
			dissColorList.push_back(dist_c);
			dissList.push_back(dist);
		}//end loop regions

		//Sort color dissimilarities for region i-th to use only K-th neighbors in color
		std::vector<double> sorted;
		std::vector<size_t> sort_index;//sorting index
		Utils::sort( dissColorList,sorted,sort_index);	

		//Compute saliency over k-th neighbors
		double sum= 0;

		for(int k=0;k<KNN;k++){
			size_t index= sort_index[k];
			double diss= dissList[index];
			sum+= diss;
		}
		double Sij= sum/(double)(KNN);
		double Si= exp(-expFalloffPar*Sij);
		double Saliency= 1.-Si;
		if(Si<Smin) Smin= Si;
		if(Si>Smax) Smax= Si;
		SList.push_back(Si);

		cout<<"Region "<<i<<": sum="<<sum<<", Sij="<<Sij<<", Si="<<Si<<", Saliency="<<Saliency<<endl;

			
		//Fill image
		int nPixelsInRegion= (int)regions[i]->fNPix;
		for(int j=0;j<nPixelsInRegion;j++){//loop on pixels inside region
			int thisPixelId= (regions[i]->fPixelCollection)[j].id;
			saliencyImg->SetBinContent(thisPixelId,Saliency);
		}//end loop pixels in region
		

	}//end loop regions
	
	cout<<"Img::GetSaliencyMap_GofermanMethod(): INFO: Saliency min/max="<<Smin<<"/"<<Smax<<endl;	
	
	/*
	//## Normalize saliency Sij and fill map
	for(int i=0;i<nRegions;i++){
		double S= SList[i];
		double Snorm= NormMin + (NormMax-NormMin)*(S-Smin)/(Smax-Smin);

		//double Saliency= 1.-exp(Snorm);
		double Saliency= 1.-Snorm;
			
		//Fill image
		int nPixelsInRegion= (int)regions[i]->fNPix;
		for(int j=0;j<nPixelsInRegion;j++){//loop on pixels inside region
			int thisPixelId= (regions[i]->fPixelCollection)[j].id;
			saliencyImg->SetBinContent(thisPixelId,Saliency);
		}//end loop pixels in region

	}//end loop regions
	*/

	if(ColorDistMatrix) ColorDistMatrix->Delete();
	if(SpatialDistMatrix) SpatialDistMatrix->Delete();
		
	return saliencyImg;

}//close GetSaliencyMap_GofermanMethod()


//Img* Img::GetMultiResoSaliencyMap(int resoMin,int resoMax,int resoStep,double beta,int minRegionSize,double knnFactor,bool useRobust,bool addCurvDist,double thr,bool addCurvMap,bool addBkgMap,bool addNoiseMap,int normalizationAcrossResoMode,double medianThrFactor,double medianImgThrFactor){
Img* Img::GetMultiResoSaliencyMap(int resoMin,int resoMax,int resoStep,double beta,int minRegionSize,double knnFactor,bool useRobustPars,double expFalloffPar, double distanceRegPar,double thr,bool addCurvMap,bool addBkgMap,bool addNoiseMap,int normalizationAcrossResoMode,double medianThrFactor,double medianImgThrFactor){

	//## Check if noise and bkg map are computed
	if(!fInterpolatedBackgroundLevelMap || !fInterpolatedBackgroundRMSMap) {
		cerr<<"Img::GetMultiResoSaliencyMap(): ERROR: Cannot get local bkg and noise maps (compute first local bkg!)"<<endl;
		return 0;
	}

	if(!this->HasStats()){
		cerr<<"Img::GetMultiResoSaliencyMap(): ERROR: No stats available for this image!"<<endl;
		return 0;
	}
	double imgMedian= fPixelStats->median;
	double imgMAD= fPixelStats->medianRMS;
	
	int nReso= (resoMax-resoMin)/resoStep + 1;
	
	TString imgName= Form("%s_saliencyMean",this->GetName());
	Img* saliencyImg_mean= (Img*)this->Clone(imgName);
	saliencyImg_mean->SetNameTitle(imgName,imgName);
	saliencyImg_mean->Reset();

	double NormMin= 0;
	double NormMax= 1;
	std::vector<Img*> salMaps;
	std::vector<double> salMapsThresholds;
	int nbins= 100;
	bool interpolateSaliency= false;
	
	for(int i=0;i<nReso;i++){
		int reso= resoMin + i*resoStep;
		cout<<"Img::GetMultiResoSaliencyMap(): INFO: Computing saliency map @ reso "<<reso<<" (step="<<resoStep<<")"<<endl;

		//Img* salMap= this->GetSaliencyMap_MultiParVersion(reso,beta,minRegionSize,knnFactor,useRobustPars,addCurvDist);
		Img* salMap= this->GetSaliencyMap(reso,beta,minRegionSize,knnFactor,useRobustPars,expFalloffPar,distanceRegPar,interpolateSaliency);
		Img* salMap_norm= salMap->GetNormalizedImage(NormMin,NormMax);
		saliencyImg_mean->Add(salMap_norm);
		
		if(salMap) salMap->Delete();
		salMaps.push_back(salMap_norm);

		//Compute stats		
		salMap_norm->ComputeStats(true,false,true);
		double salMedian= (salMap_norm->GetPixelStats())->median;
		double medianThr= medianThrFactor*salMedian;
		double otsuThr= salMap_norm->FindOtsuThreshold(nbins);
		double salThr= std::max(medianThr,otsuThr);
		salMapsThresholds.push_back(salThr);
	}//end loop reso
	
	//Normalize final saliency
	saliencyImg_mean->Scale(1./(double)nReso);
	
	cout<<"Img::GetMultiResoSaliencyMap(): INFO: Normalize saliency sum over reso..."<<endl;
	imgName= Form("%s_saliencyCombined",this->GetName());
	Img* saliencyImg= (Img*)saliencyImg_mean->Clone(imgName);
	saliencyImg->SetNameTitle(imgName,imgName);


	//Multi-resolution combination
	if(nReso>1){
		saliencyImg->Reset();

		double minSaliency= saliencyImg_mean->GetMinimum();
		double imgMedianThr= medianImgThrFactor*imgMedian;
	
		if(normalizationAcrossResoMode==1){//normalization with fixed threshold
			for(int i=0;i<saliencyImg->GetNbinsX();i++){
				for(int j=0;j<saliencyImg->GetNbinsY();j++){
					double imgBinContent= this->GetBinContent(i+1,j+1);
					//if(imgBinContent==0) continue;//GIUSTO
					//bool isPositiveExcess= (imgBinContent>imgMedianThr);//GIUSTO
					bool isPositiveExcess= (imgBinContent>imgMedianThr && imgBinContent!=0);
					double w= saliencyImg_mean->GetBinContent(i+1,j+1);
					double wmin= 1.e+99;
					double wmax= -1.e+99;
					for(int k=0;k<salMaps.size();k++){
						double thisw= salMaps[k]->GetBinContent(i+1,j+1);
						if(thisw<wmin) wmin= thisw;
						if(thisw>=wmax) wmax= thisw;
					}//end loop multi reso
					if(w>=thr && isPositiveExcess) saliencyImg->SetBinContent(i+1,j+1,wmax);
					else if(w>=thr && !isPositiveExcess) saliencyImg->SetBinContent(i+1,j+1,minSaliency);
					else saliencyImg->SetBinContent(i+1,j+1,wmin);
				}//end loop bins Y
			}//end loop bins X
		}//close if
	
		else if(normalizationAcrossResoMode==2){//normalization with adaptive threshold

			double salientMultiplicityFactorThr= thr;//this parameter is the fixed threshold in the previous mode
			int salientMultiplicityThr= std::round(salientMultiplicityFactorThr*nReso);

			for(int i=0;i<saliencyImg->GetNbinsX();i++){
				for(int j=0;j<saliencyImg->GetNbinsY();j++){
					double imgBinContent= this->GetBinContent(i+1,j+1);	
					//if(imgBinContent==0) continue;//GIUSTO
					//bool isPositiveExcess= (imgBinContent>imgMedianThr);//GIUSTO
					bool isPositiveExcess= (imgBinContent>imgMedianThr && imgBinContent!=0);
					double w= saliencyImg_mean->GetBinContent(i+1,j+1);
					double wmin= 1.e+99;
					double wmax= -1.e+99;
					int saliencyMultiplicity= 0;
					for(int k=0;k<salMaps.size();k++){
						double thisw= salMaps[k]->GetBinContent(i+1,j+1);
						double thisThreshold= salMapsThresholds[k];
						if(thisw>thisThreshold) saliencyMultiplicity++;
						if(thisw<wmin) wmin= thisw;
						if(thisw>=wmax) wmax= thisw;
					}//end loop multi reso
					if(saliencyMultiplicity>=salientMultiplicityThr && isPositiveExcess) saliencyImg->SetBinContent(i+1,j+1,wmax);
					else if(saliencyMultiplicity>=salientMultiplicityThr && !isPositiveExcess) saliencyImg->SetBinContent(i+1,j+1,minSaliency);
					else saliencyImg->SetBinContent(i+1,j+1,wmin);
				}//end loop bins Y
			}//end loop bins X
		}//close else if
	
		else if(normalizationAcrossResoMode==3){//normalization with max across scales
			for(int i=0;i<saliencyImg->GetNbinsX();i++){
				for(int j=0;j<saliencyImg->GetNbinsY();j++){
					double imgBinContent= this->GetBinContent(i+1,j+1);
					//if(imgBinContent==0) continue;		//GIUSTO		
					//bool isPositiveExcess= (imgBinContent>imgMedianThr);//GIUSTO
					bool isPositiveExcess= (imgBinContent>imgMedianThr && imgBinContent!=0);
					double w= saliencyImg_mean->GetBinContent(i+1,j+1);
					double wmin= 1.e+99;
					double wmax= -1.e+99;
					for(int k=0;k<salMaps.size();k++){
						double thisw= salMaps[k]->GetBinContent(i+1,j+1);
						if(thisw<wmin) wmin= thisw;
						if(thisw>=wmax) wmax= thisw;
					}//end loop multi reso
					if(isPositiveExcess) saliencyImg->SetBinContent(i+1,j+1,wmax);
					else saliencyImg->SetBinContent(i+1,j+1,wmin);
				}//end loop bins Y
			}//end loop bins X
		}//close else if
		else if(normalizationAcrossResoMode==4){//normalization with mean of scales
			for(int i=0;i<saliencyImg->GetNbinsX();i++){
				for(int j=0;j<saliencyImg->GetNbinsY();j++){
					double imgBinContent= this->GetBinContent(i+1,j+1);
					//if(imgBinContent==0) continue;//GIUSTO					
					double w= saliencyImg_mean->GetBinContent(i+1,j+1);
					saliencyImg->SetBinContent(i+1,j+1,w);
				}//end loop bins Y
			}//end loop bins X
		}//close else if
		else {//default normalization with mean of scales
			for(int i=0;i<saliencyImg->GetNbinsX();i++){
				for(int j=0;j<saliencyImg->GetNbinsY();j++){
					double imgBinContent= this->GetBinContent(i+1,j+1);
					//if(imgBinContent==0) continue;//GIUSTO				
					double w= saliencyImg_mean->GetBinContent(i+1,j+1);
					saliencyImg->SetBinContent(i+1,j+1,w);	
				}//end loop bins Y
			}//end loop bins X

		}//close else
	}//close if multi-resolution

	
	//Clear scale saliency maps
	for(unsigned int k=0;k<salMaps.size();k++){
		if(salMaps[k]) salMaps[k]->Delete();		
	}
	
	//Normalize bkg and noise maps
	if(addBkgMap){
		cout<<"Img::GetMultiResoSaliencyMap(): INFO: Normalize bkg map..."<<endl;
		Img* bkgImg= fInterpolatedBackgroundLevelMap->GetNormalizedImage(NormMin,NormMax);
		saliencyImg->Add(bkgImg);	
		if(bkgImg) bkgImg->Delete();
	}

	if(addNoiseMap){
		Img* noiseImg= fInterpolatedBackgroundRMSMap->GetNormalizedImage(NormMin,NormMax);
		saliencyImg->Add(noiseImg);
		if(noiseImg) noiseImg->Delete();
	}

	//Compute image curvature
	if(addCurvMap){
		Img* curvImg= this->GetLoGImage(true);
		Img* curvImg_norm= curvImg->GetNormalizedImage(NormMin,NormMax);
		saliencyImg->Add(curvImg_norm);
		if(curvImg) curvImg->Delete();
		if(curvImg_norm) curvImg_norm->Delete();
	}

	//Normalize final map
	cout<<"Img::GetMultiResoSaliencyMap(): INFO: Normalize final maps..."<<endl;
	imgName= Form("%s_saliencyMultiReso",this->GetName());
	Img* saliencyMap= saliencyImg->GetNormalizedImage(NormMin,NormMax);
	saliencyMap->SetNameTitle(imgName,imgName);
	if(saliencyImg) saliencyImg->Delete();

	for(int i=0;i<saliencyMap->GetNbinsX();i++){
		for(int j=0;j<saliencyMap->GetNbinsY();j++){
			double imgBinContent= this->GetBinContent(i+1,j+1);
			if(imgBinContent==0){
				saliencyMap->SetBinContent(i+1,j+1,0);
			}
		}
	}		

	return saliencyMap;
	
}//close GetMultiResoSaliencyMap()



Img* Img::GetMultiResoSmoothedSaliencyMap(int resoMin,int resoMax,int resoStep,double beta,int minRegionSize,double sigmaS,double sigmaC,double saliencyKFactor,bool useRobust,bool addCurvDist,double thr,bool addCurvMap,bool addBkgMap,bool addNoiseMap){

	//## Check if noise and bkg map are computed
	if(!fInterpolatedBackgroundLevelMap || !fInterpolatedBackgroundRMSMap) {
		cerr<<"Img::GetMultiResoSmoothedSaliencyMap(): ERROR: Cannot get local bkg and noise maps (compute first local bkg!)"<<endl;
		return 0;
	}

	int nReso= (resoMax-resoMin)/resoStep + 1;
	
	TString imgName= Form("%s_saliencyMean",this->GetName());
	Img* saliencyImg_mean= (Img*)this->Clone(imgName);
	saliencyImg_mean->SetNameTitle(imgName,imgName);
	saliencyImg_mean->Reset();

	double NormMin= 0;
	double NormMax= 1;
	std::vector<Img*> salMaps;

	for(int i=0;i<nReso;i++){
		int reso= resoMin + i*resoStep;
		cout<<"Img::GetMultiResoSmoothedSaliencyMap(): INFO: Computing saliency map @ reso "<<reso<<" (step="<<resoStep<<")"<<endl;
		
		Img* salMap= this->GetSmoothedSaliencyMap(reso,beta,minRegionSize,sigmaS,sigmaC,saliencyKFactor,useRobust,addCurvDist);
		Img* salMap_norm= salMap->GetNormalizedImage(NormMin,NormMax);
		saliencyImg_mean->Add(salMap_norm);
		//saliencyImg_mean->Add(salMap);

		if(salMap) salMap->Delete();
		//if(salMap_norm) salMap_norm->Delete();
		salMaps.push_back(salMap_norm);
		//salMaps.push_back(salMap);
	}//end loop reso
	
	//Combine multi-reso saliency maps
	saliencyImg_mean->Scale(1./(double)nReso);
	
	cout<<"Img::GetMultiResoSmoothedSaliencyMap(): INFO: Normalize saliency sum over reso..."<<endl;
	imgName= Form("%s_saliencyCombined",this->GetName());
	Img* saliencyImg= (Img*)saliencyImg_mean->Clone(imgName);
	saliencyImg->SetNameTitle(imgName,imgName);
	saliencyImg->Reset();
	
	for(int i=0;i<saliencyImg->GetNbinsX();i++){
		for(int j=0;j<saliencyImg->GetNbinsY();j++){
			double w= saliencyImg_mean->GetBinContent(i+1,j+1);
			double wmin= 1.e+99;
			double wmax= -1.e+99;
			for(int k=0;k<salMaps.size();k++){
				double thisw= salMaps[k]->GetBinContent(i+1,j+1);
				if(thisw<wmin) wmin= thisw;
				if(thisw>wmax) wmax= thisw;
			}//end loop multi reso
			if(w>=thr) saliencyImg->SetBinContent(i+1,j+1,wmax);
			else if(w>=0 && w<thr) saliencyImg->SetBinContent(i+1,j+1,wmin);
			else saliencyImg->SetBinContent(i+1,j+1,w);	
			
		}//end loop bins Y
	}//end loop bins X

	
	
	//if(saliencyImg_mean) saliencyImg_mean->Delete();
	for(unsigned int k=0;k<salMaps.size();k++){
		if(salMaps[k]) salMaps[k]->Delete();		
	}
	
	//Normalizee bkg and noise maps
	if(addBkgMap){
		cout<<"Img::GetMultiResoSmoothedSaliencyMap(): INFO: Normalize bkg map..."<<endl;
		Img* bkgImg= fInterpolatedBackgroundLevelMap->GetNormalizedImage(NormMin,NormMax);
		saliencyImg->Add(bkgImg);	
		if(bkgImg) bkgImg->Delete();
	}

	if(addNoiseMap){
		Img* noiseImg= fInterpolatedBackgroundRMSMap->GetNormalizedImage(NormMin,NormMax);
		saliencyImg->Add(noiseImg);
		if(noiseImg) noiseImg->Delete();
	}

	//Compute image curvature
	if(addCurvMap){
		Img* curvImg= this->GetLoGImage(true);
		Img* curvImg_norm= curvImg->GetNormalizedImage(NormMin,NormMax);
		saliencyImg->Add(curvImg_norm);
		if(curvImg) curvImg->Delete();
		if(curvImg_norm) curvImg_norm->Delete();
	}

	//Normalize final map
	cout<<"Img::GetMultiResoSmoothedSaliencyMap(): INFO: Normalize final maps..."<<endl;
	imgName= Form("%s_saliencyMultiReso",this->GetName());
	Img* saliencyMap= saliencyImg->GetNormalizedImage(NormMin,NormMax);
	saliencyMap->SetNameTitle(imgName,imgName);
	if(saliencyImg) saliencyImg->Delete();

	return saliencyMap;

}//close GetMultiResoSmoothedSaliencyMap()


Img* Img::GetSaliencyMap_SpectralRes(){
	
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();

	//## Convert image to OpenCV mat
	cv::Mat img= this->ImgToMat();

	//## Compute saliency map
	//instantiates the specific Saliency
  cv::Ptr<cv::saliency::Saliency> saliencyAlgorithm = cv::saliency::Saliency::create("SPECTRAL_RESIDUAL");

	cv::Mat salimg;
	//cv::saliency::Saliency::computeSaliency	(img,salimg);
	saliencyAlgorithm->computeSaliency( img,salimg );

	//## Fill image
	TString imgName= Form("%s_saliency",std::string(this->GetName()).c_str());	
	Img* SaliencyImg= (Img*)this->Clone(imgName);
	SaliencyImg->Reset();
	for(int j=0;j<Ny;j++){
		int rowId= Ny-1-j;
		for(int i=0;i<Nx;i++){
			int colId= i;
			double binContent= this->GetBinContent(i+1,j+1);
			if(binContent==0) continue;
			double matrixElement= (double)(salimg.at<double>(rowId,colId));
			SaliencyImg->SetBinContent(i+1,j+1,matrixElement);
		}//end loop x bins
	}//end loop y bins

	return SaliencyImg;

}//close GetSaliencyMap_SpectralRes()	


std::vector<Source*> Img::FindMultiScaleBlobs(Img* blobMask,double seedThr,double mergeThr,int minNPix){

	//## Init
	std::vector<Source*> Blobs;
	Blobs.clear();
	if(!blobMask) return Blobs;

	//## Find compact blobs
	bool findNegativeExcess= false;
	bool findNestedSources= false;
	bool useLocalBackground= false;
	bool mergeBelowSeed= false;
	double peakThreshold= seedThr;
	Blobs= this->FindBlobs(blobMask,seedThr,mergeThr,minNPix,findNegativeExcess,useLocalBackground,mergeBelowSeed);
	int nBlobs= (int)Blobs.size();
	cout<<"Img::FindMultiScaleBlobs(): INFO: "<<nBlobs<<" blobs detected..."<<endl;

	//## Loop over detected blobs and cmpute pars
	double fluxCorrection= 1;
	if(this->HasMetaData()){
		double fx= this->metadata.Bmaj;
		double fy= this->metadata.Bmin;
		fluxCorrection= TMath::Pi()*fx*fy/(4*log(2));
	}

	for(int k=0;k<nBlobs;k++){
		if(!Blobs[k]) continue;

		//## Set source name & id & type
		std::string sourceName= std::string(Form("Img%s_sId%d",std::string(this->GetName()).c_str(),k));
		Blobs[k]->SetId(k);
		Blobs[k]->SetName(sourceName);
		Blobs[k]->SetType(Source::eCompact);

		//## Compute stats parameters
		cout<<"Img::FindMultiScaleBlobs(): INFO: Computing source stats..."<<endl;
		Blobs[k]->ComputeStats();

		//## Compute morphology parameters
		cout<<"Img::FindMultiScaleBlobs(): INFO: Computing morphology params..."<<endl;
		Blobs[k]->ComputeMorphologyParams();
	
		//## Set flux correction factor
		Blobs[k]->fFluxCorrection= fluxCorrection;
	}//end loop blobs

	return Blobs;

}//close Img::FindMultiScaleBlobs()

std::vector<Source*> Img::FindMultiScaleBlobs(int kernelFactor,double sigmaMin,double sigmaMax,double sigmaStep,int thrModel,double thrFactor,double seedThrFactor,double mergeThrFactor,int minNPix){

	//## Init
	std::vector<Source*> Blobs;
	Blobs.clear();
	int nScales= (sigmaMax-sigmaMin)/sigmaStep + 1;
	double seedThr= std::round(seedThrFactor*nScales);
	double mergeThr= std::round(mergeThrFactor*nScales);

	//## Compute blobs
	Img* blobMask= this->GetMultiScaleBlobMask(kernelFactor,sigmaMin,sigmaMax,sigmaStep,thrModel,thrFactor);
	if(!blobMask) return Blobs;

	//## Find compact blobs
	Blobs= this->FindMultiScaleBlobs(blobMask,seedThr,mergeThr,minNPix);

	return Blobs;

}//close GetMultiScaleBlobs()

Img* Img::GetMultiScaleBlobMask(int kernelFactor,double sigmaMin,double sigmaMax,double sigmaStep,int thrModel,double thrFactor){

	//Get image stats and thresholds
	if(!this->HasStats()){
		cerr<<"Img::GetMultiScaleBlobMask(): ERROR: No stats available for this image!"<<endl;
		return 0;
	}
	double imgMedian= fPixelStats->median;
	double imgMAD= fPixelStats->medianRMS;
	
	//Init scales
	int nScales= (sigmaMax-sigmaMin)/sigmaStep + 1;
	
	TString imgName= Form("%s_blobMask",this->GetName());
	Img* blobMask= (Img*)this->Clone(imgName);
	blobMask->SetNameTitle(imgName,imgName);
	blobMask->Reset();

	double NormMin= 0;
	double NormMax= 1;
	int nbins= 100;
	std::vector<Img*> filterMaps;
	std::vector<double> thresholdLevels;
	
	for(int i=0;i<nScales;i++){
		double sigma= sigmaMin + i*sigmaStep;
		int kernelSize= kernelFactor*sigma;	
		if(kernelSize%2==0) kernelSize++;
		cout<<"Img::GetMultiScaleBlobMask(): INFO: Computing LoG map @ scale "<<sigma<<" (step="<<sigmaStep<<", kernsize="<<kernelSize<<")"<<endl;

		//Compute LoG filter
		Img* filterMap= this->GetNormLoGImage(kernelSize,sigma,true);
		filterMap->ComputeStats(true,false,true);
		filterMaps.push_back(filterMap);

		//Compute threshold levels	
		double median= (filterMap->GetPixelStats())->median;
		double medianRMS= (filterMap->GetPixelStats())->medianRMS;
		double medianThr= thrFactor*median;
		double medianRMSThr= thrFactor*medianRMS;
		double otsuThr= filterMap->FindOtsuThreshold(nbins);
		double optimalThr= std::max(medianThr,otsuThr);
		double thrLevel= medianRMSThr;
		if(thrModel==1) thrLevel= optimalThr;
		else if(thrModel==2) thrLevel= medianRMSThr;
		else thrLevel= medianRMSThr;
		thresholdLevels.push_back(thrLevel);	
	}//end loop reso
	
	//Find blobs across scales
	for(int i=0;i<blobMask->GetNbinsX();i++){
		for(int j=0;j<blobMask->GetNbinsY();j++){
			double binContent= this->GetBinContent(i+1,j+1);
			if(binContent==0) continue;

			double wsum= 0;
			int counter= 0;
			for(int k=0;k<filterMaps.size();k++){
				double w= filterMaps[k]->GetBinContent(i+1,j+1);
				if(w<thresholdLevels[k]) continue;
				wsum+= w;
				counter++;
			}//end loop scales

			blobMask->SetBinContent(i+1,j+1,counter);
		}//end loop y
	}//end loop x

	return blobMask;

}//close GetMultiScaleBlobMask()


Img* Img::GetMultiScaleBlobFilterMap(int kernelFactor,double sigmaMin,double sigmaMax,double sigmaStep,double significantScaleFractionThr,int normalizationAcrossScaleMode,double thrFactor,double imgThrFactor){

	//Get image stats and thresholds
	if(!this->HasStats()){
		cerr<<"Img::GetMultiScaleBlobFilterMap(): ERROR: No stats available for this image!"<<endl;
		return 0;
	}
	double imgMedian= fPixelStats->median;
	double imgMAD= fPixelStats->medianRMS;
	double imgMedianThr= imgThrFactor*imgMedian;
	
	//Init scales
	int nScales= (sigmaMax-sigmaMin)/sigmaStep + 1;
	
	TString imgName= Form("%s_blobFilterMean",this->GetName());
	Img* filterMap_mean= (Img*)this->Clone(imgName);
	filterMap_mean->SetNameTitle(imgName,imgName);
	filterMap_mean->Reset();

	double NormMin= 0;
	double NormMax= 1;
	std::vector<Img*> filterMaps;
	std::vector<Img*> filterMaps_norm;
	std::vector<double> optimalThresholds;
	std::vector<double> rmsThresholds;
	int nbins= 100;
	
	for(int i=0;i<nScales;i++){
		double sigma= sigmaMin + i*sigmaStep;
		int kernelSize= kernelFactor*sigma;	
		if(kernelSize%2==0) kernelSize++;
		cout<<"Img::GetMultiScaleBlobFilterMap(): INFO: Computing LoG map @ scale "<<sigma<<" (step="<<sigmaStep<<", kernsize="<<kernelSize<<")"<<endl;

		Img* filterMap= this->GetNormLoGImage(kernelSize,sigma,true);
		Img* filterMap_norm= filterMap->GetNormalizedImage(NormMin,NormMax);
		
		//if(filterMap) filterMap->Delete();
		filterMaps.push_back(filterMap_norm);
		filterMaps_norm.push_back(filterMap_norm);

		//Compute threshold levels for norm maps		
		filterMap_norm->ComputeStats(true,false,true);
		double median_norm= (filterMap_norm->GetPixelStats())->median;
		double medianRMS_norm= (filterMap_norm->GetPixelStats())->medianRMS;
		double medianThr_norm= thrFactor*median_norm;
		double medianRMSThr_norm= thrFactor*medianRMS_norm;
		double otsuThr_norm= filterMap_norm->FindOtsuThreshold(nbins);
		double optimalThr_norm= std::max(medianThr_norm,otsuThr_norm);
		optimalThresholds.push_back(optimalThr_norm);	

		//Compute threshold levels for maps		
		filterMap->ComputeStats(true,false,true);
		double median= (filterMap->GetPixelStats())->median;
		double medianRMS= (filterMap->GetPixelStats())->medianRMS;
		double medianThr= thrFactor*median;
		double medianRMSThr= thrFactor*medianRMS;
		double otsuThr= filterMap->FindOtsuThreshold(nbins);
		double optimalThr= std::max(medianThr,otsuThr);
		rmsThresholds.push_back(medianRMSThr);

		//Add this map to mean map	
		filterMap_mean->Add(filterMap_norm);
	}//end loop reso
	
	//Normalize final saliency
	filterMap_mean->Scale(1./(double)nScales);
	filterMap_mean->ComputeStats(true,false,true);
	double median_mean= (filterMap_mean->GetPixelStats())->median;
	double medianRMS_mean= (filterMap_mean->GetPixelStats())->medianRMS;
	double min_mean= (filterMap_mean->GetPixelStats())->min;
	double medianThr_mean= thrFactor*median_mean;
	
	cout<<"Img::GetMultiScaleBlobFilterMap(): INFO: Normalize saliency sum over reso..."<<endl;
	imgName= Form("%s_blobFilterCombined",this->GetName());
	Img* filterImg= (Img*)filterMap_mean->Clone(imgName);
	filterImg->SetNameTitle(imgName,imgName);
	filterImg->Reset();

	
	//## Normalization on mean filter map with fixed threshold
	if(normalizationAcrossScaleMode==1){
		for(int i=0;i<filterMap_mean->GetNbinsX();i++){
			for(int j=0;j<filterMap_mean->GetNbinsY();j++){
				double imgBinContent= this->GetBinContent(i+1,j+1);
				if(imgBinContent==0) continue;
				bool isPositiveExcess= (imgBinContent>imgMedianThr);
				double w= filterMap_mean->GetBinContent(i+1,j+1);

				//Find min/max across scales
				double wmin= 1.e+99;
				double wmax= -1.e+99;
				for(int k=0;k<filterMaps_norm.size();k++){
					double thisw= filterMaps_norm[k]->GetBinContent(i+1,j+1);
					if(thisw<wmin) wmin= thisw;
					if(thisw>=wmax) wmax= thisw;
				}//end loop multi reso

				//Apply normalization
				if(w>=medianThr_mean){
					if(isPositiveExcess) filterImg->SetBinContent(i+1,j+1,wmax);
					else filterImg->SetBinContent(i+1,j+1,min_mean);
				}
				else {
					filterImg->SetBinContent(i+1,j+1,wmin);
				}
			}//end loop bins Y
		}//end loop bins X
	}//close if
		
	//## Normalization with adaptive optimal threshold
	else if(normalizationAcrossScaleMode==2){

		double salientMultiplicityFactorThr= significantScaleFractionThr;
		int salientMultiplicityThr= std::round(salientMultiplicityFactorThr*nScales);

		for(int i=0;i<filterMap_mean->GetNbinsX();i++){
			for(int j=0;j<filterMap_mean->GetNbinsY();j++){
				double imgBinContent= this->GetBinContent(i+1,j+1);	
				if(imgBinContent==0) continue;
				bool isPositiveExcess= (imgBinContent>imgMedianThr);
				double w= filterMap_mean->GetBinContent(i+1,j+1);
				double wmin= 1.e+99;
				double wmax= -1.e+99;
				int saliencyMultiplicity= 0;
				for(int k=0;k<filterMaps_norm.size();k++){
					double thisw= filterMaps_norm[k]->GetBinContent(i+1,j+1);
					double thisThreshold= optimalThresholds[k];
					if(thisw>thisThreshold) saliencyMultiplicity++;
					if(thisw<wmin) wmin= thisw;
					if(thisw>=wmax) wmax= thisw;
				}//end loop multi reso

				//Apply normalization
				if(saliencyMultiplicity>=salientMultiplicityThr){ 
					if(isPositiveExcess) filterImg->SetBinContent(i+1,j+1,wmax);
					else filterImg->SetBinContent(i+1,j+1,min_mean);
				}
				else {
					filterImg->SetBinContent(i+1,j+1,wmin);
				}
			}//end loop bins Y
		}//end loop bins X
	}//close else if
	
	else if(normalizationAcrossScaleMode==3){//normalization with max across scales
		for(int i=0;i<filterMap_mean->GetNbinsX();i++){
			for(int j=0;j<filterMap_mean->GetNbinsY();j++){
				double imgBinContent= this->GetBinContent(i+1,j+1);
				if(imgBinContent==0) continue;				
				bool isPositiveExcess= (imgBinContent>imgMedianThr);
				double w= filterMap_mean->GetBinContent(i+1,j+1);
				double wmin= 1.e+99;
				double wmax= -1.e+99;
				for(int k=0;k<filterMaps_norm.size();k++){
					double thisw= filterMaps_norm[k]->GetBinContent(i+1,j+1);
					if(thisw<wmin) wmin= thisw;
					if(thisw>=wmax) wmax= thisw;
				}//end loop multi reso
				if(isPositiveExcess) filterImg->SetBinContent(i+1,j+1,wmax);
				else filterImg->SetBinContent(i+1,j+1,wmin);
			}//end loop bins Y
		}//end loop bins X

	}//close else if

	//## Normalization with mean of scales
	else if(normalizationAcrossScaleMode==4){
		for(int i=0;i<filterMap_mean->GetNbinsX();i++){
			for(int j=0;j<filterMap_mean->GetNbinsY();j++){
				double imgBinContent= this->GetBinContent(i+1,j+1);
				if(imgBinContent==0) continue;								
				double w= filterMap_mean->GetBinContent(i+1,j+1);
				filterImg->SetBinContent(i+1,j+1,w);
			}//end loop bins Y
		}//end loop bins X
	}//close else if

	//## Normalization with rms adaptive threshold
	else if(normalizationAcrossScaleMode==5){

		double salientMultiplicityFactorThr= significantScaleFractionThr;
		int salientMultiplicityThr= std::round(salientMultiplicityFactorThr*nScales);

		for(int i=0;i<filterMap_mean->GetNbinsX();i++){
			for(int j=0;j<filterMap_mean->GetNbinsY();j++){
				double imgBinContent= this->GetBinContent(i+1,j+1);	
				if(imgBinContent==0) continue;
				bool isPositiveExcess= (imgBinContent>imgMedianThr);
				double w= filterMap_mean->GetBinContent(i+1,j+1);
				double wmin= 1.e+99;
				double wmax= -1.e+99;
				int saliencyMultiplicity= 0;
				for(int k=0;k<filterMaps.size();k++){
					double thisw= filterMaps[k]->GetBinContent(i+1,j+1);
					double thisThreshold= rmsThresholds[k];
					if(thisw>thisThreshold) saliencyMultiplicity++;
					if(thisw<wmin) wmin= thisw;
					if(thisw>=wmax) wmax= thisw;
				}//end loop multi reso
			
				if(saliencyMultiplicity>=salientMultiplicityThr){
					if(isPositiveExcess) filterImg->SetBinContent(i+1,j+1,wmax);
				  else filterImg->SetBinContent(i+1,j+1,min_mean);
				}
				else {
					filterImg->SetBinContent(i+1,j+1,wmin);
				}
			}//end loop bins Y
		}//end loop bins X
	}//close else if

	else {//default normalization with mean of scales
		for(int i=0;i<filterMap_mean->GetNbinsX();i++){
			for(int j=0;j<filterMap_mean->GetNbinsY();j++){
				double imgBinContent= this->GetBinContent(i+1,j+1);
				if(imgBinContent==0) continue;				
				double w= filterMap_mean->GetBinContent(i+1,j+1);
				filterImg->SetBinContent(i+1,j+1,w);	
			}//end loop bins Y
		}//end loop bins X
	}//close else


	//## Clear maps
	for(unsigned int k=0;k<filterMaps.size();k++){
		if(filterMaps[k]) filterMaps[k]->Delete();	
		if(filterMaps_norm[k]) filterMaps_norm[k]->Delete();		
	}
	
	
	//## Normalize final map
	cout<<"Img::GetMultiScaleBlobFilterMap(): INFO: Normalize final maps..."<<endl;
	imgName= Form("%s_blobFilterMultiScale",this->GetName());
	Img* filterImg_final= filterImg->GetNormalizedImage(NormMin,NormMax);
	filterImg_final->SetNameTitle(imgName,imgName);
	if(filterImg) filterImg->Delete();

	for(int i=0;i<filterImg_final->GetNbinsX();i++){
		for(int j=0;j<filterImg_final->GetNbinsY();j++){
			double imgBinContent= this->GetBinContent(i+1,j+1);
			if(imgBinContent==0){
				filterImg_final->SetBinContent(i+1,j+1,0);
			}
		}
	}		

	return filterImg_final;
	
}//close GetMultiScaleBlobFilterMap()



std::vector<TLine*> Img::LinearHoughTransform(double edgeThreshold,double accumulatorReso,double thetaReso,int threshold, double minLineLength, double maxLineGap){

	//## Init data
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	TEllipse* aCircle= 0;
	std::vector<TEllipse*> detectedCircles;
	detectedCircles.clear();

	//## Get normalized image
	Img* img_norm= this->GetNormalizedImage(0,255);

	//## Convert image to OpenCV mat
	cv::Mat img= img_norm->ImgToMat();

	//## Convert it to gray
	img.convertTo(img, CV_8UC1);

	//## Threshold image and define kernel
	cv::threshold(img, img, threshold, 255, cv::THRESH_BINARY);

	//Mat dst, cdst;
 	//cv::Canny(mat, dst, 50, 200, 3);

	//## Compute Linear Hough transform
	double thetaReso_rad= thetaReso*TMath::DegToRad();
	std::vector<cv::Vec4i> lines;
  cv::HoughLinesP(img, lines, accumulatorReso,thetaReso_rad, threshold, minLineLength, maxLineGap);

	std::vector<TLine*> detected_lines;
	TLine* aLine= 0;

	int nLines= (int)lines.size();
	cout<<"Img::LinearHoughTransform(): INFO: "<<nLines<<" detected..."<<endl;

  for(int i = 0; i <nLines; i++ ) {
    cv::Vec4i l = lines[i];
		cv::Point startPt(l[0], l[1]);
		cv::Point endPt(l[2], l[3]);

		int ix1= startPt.x;
		int iy1= Ny-1-startPt.y; 
		double x1= this->GetXaxis()->GetBinCenter(ix1+1);
		double y1= this->GetYaxis()->GetBinCenter(iy1+1);
		int ix2= endPt.x;
		int iy2= Ny-1-endPt.y; 
		double x2= this->GetXaxis()->GetBinCenter(ix2+1);
		double y2= this->GetYaxis()->GetBinCenter(iy2+1);
	
		aLine= new TLine(x1,y1,x2,y2);
		aLine->SetLineWidth(2);
		aLine->SetLineStyle(kSolid);
		aLine->SetLineColor(kRed);
    detected_lines.push_back(aLine);
  }

	img_norm->Delete();
	return detected_lines;

}//close 


std::vector<TEllipse*> Img::HoughTransform(int accumulatorReso,int minCircleDistance,int minCircleRadius,int maxCircleRadius,double edgeThreshold,double accumulatorThreshold){

	//## Init data
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	TEllipse* aCircle= 0;
	std::vector<TEllipse*> detectedCircles;
	detectedCircles.clear();

	//## Get normalized image
	Img* img_norm= this->GetNormalizedImage(0,255);

	//## Convert image to OpenCV mat
	cv::Mat mat= img_norm->ImgToMat();

	//## Convert it to gray
	mat.convertTo(mat, CV_8UC1);
  
	//## Apply the Hough Transform to find the circles
	double param1= 200;//threshold for Canny edge
	double param2= 100;//threshold for accumulator
	std::vector<cv::Vec3f> circles;
  cv::HoughCircles( mat, circles, CV_HOUGH_GRADIENT, accumulatorReso, minCircleDistance, edgeThreshold, accumulatorThreshold, minCircleRadius, maxCircleRadius);

  //## Fill the detected circle list
	int nCircles= (int)circles.size();

	cout<<"Img::HoughTransform(): INFO: "<<nCircles<<" circles detected!"<<endl;

  for(int i=0;i<nCircles; i++) {
		cv::Point center(cvRound(circles[i][0]), cvRound(circles[i][1]));
    int radius = cvRound(circles[i][2]);
    int rowId= Ny-1-center.y;
		//int colId= Nx-1-center.x;
		int colId= center.x;
		//int centerx_transf= colId + boundingBoxX[0];
		//int centery_transf= rowId + boundingBoxY[0];
 
		int ix= center.x;
		int iy= Ny-1-center.y; 
		double x= this->GetXaxis()->GetBinCenter(ix+1);
		double y= this->GetYaxis()->GetBinCenter(iy+1);
		aCircle= new TEllipse(x,y,radius,radius);
		aCircle->SetLineColor(kRed);
		aCircle->SetLineWidth(2);
		aCircle->SetFillColor(0);
		aCircle->SetFillStyle(0);
		detectedCircles.push_back(aCircle);
  	cout<<"Img::HoughTransform(): INFO: Circle #"<<i<<": C("<<x<<","<<y<<") ixiy("<<ix<<","<<iy<<")"<<endl;
		
	}//end loop circles

	img_norm->Delete();

	return detectedCircles;

}//close HoughTransform()



TCanvas* Img::DrawPlain(bool useLogScale){

	Img* img_norm= 0;
	if(useLogScale) img_norm= this->GetLogNormalizedImage(1,256,false);
	else img_norm= this->GetNormalizedImage(1,256,false);

	TCanvas* canv = new TCanvas("image", "image");
	canv->cd();
  canv->ToggleEventStatus();
  canv->SetRightMargin(0.0);
  canv->SetLeftMargin(0.0);
  canv->SetTopMargin(0.0);
  canv->SetBottomMargin(0.0);

	img_norm->SetStats(0);
	img_norm->Draw("COLA");

	return canv;

}//close DrawPlain()

Img* Img::GetTile(int ix_min,int ix_max,int iy_min,int iy_max){

	//## Check args
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	
	//## Check tile sizes
	if(ix_min<0 || ix_min>=Nx || ix_min>ix_max){
		std::stringstream errMsg;
		errMsg<<"Invalid min tile X given ("<<ix_min<<")";
		errflag= errMsg.str();
		cerr<<"Img::GetTile(): ERROR: "<<errflag<<endl;
		return NULL;
	}
	if(ix_max<0 || ix_max>=Nx || ix_max<ix_min){
		std::stringstream errMsg;
		errMsg<<"Invalid max tile X given ("<<ix_max<<")";
		errflag= errMsg.str();
		cerr<<"Img::GetTile(): ERROR: "<<errflag<<endl;
		return NULL;
	}
	if(iy_min<0 || iy_min>=Ny || iy_min>iy_max){
		std::stringstream errMsg;
		errMsg<<"Invalid min tile Y given ("<<iy_min<<")";
		errflag= errMsg.str();
		cerr<<"Img::GetTile(): ERROR: "<<errflag<<endl;
		return NULL;
	}
	if(iy_max<0 || iy_max>=Ny || iy_max<iy_min){
		std::stringstream errMsg;
		errMsg<<"Invalid max tile Y given ("<<iy_max<<")";
		errflag= errMsg.str();
		cerr<<"Img::GetTile(): ERROR: "<<errflag<<endl;
		return NULL;
	}

	int TileSizeX= ix_max-ix_min+1;
	int TileSizeY= iy_max-iy_min+1;

	double binX_min= this->GetXaxis()->GetBinCenter(ix_min+1);
	double binX_max= this->GetXaxis()->GetBinCenter(ix_max+1);
	double binY_min= this->GetYaxis()->GetBinCenter(iy_min+1);
	double binY_max= this->GetYaxis()->GetBinCenter(iy_max+1);
	
	//## Init tile
	Img* tile= new Img;
	TString tileName= Form("%s_tile_x%d_%d_y%d_%d",this->GetName(),ix_min,ix_max,iy_min,iy_max);
	tile->SetNameTitle(tileName,tileName);
	tile->SetBins(TileSizeX,binX_min-0.5,binX_max+0.5,TileSizeY,binY_min-0.5,binY_max+0.5);
	tile->SetMetaData(this->metadata);
	tile->Reset();
	
	//## Read tile
	//cout<<"Img::GetTile(): INFO: Read tile ix("<<ix_min<<","<<ix_max<<") iy("<<iy_min<<","<<iy_max<<") binX("<<binX_min<<","<<binX_max<<")"<<" binY("<<binY_min<<","<<binY_max<<")"<<endl;
	for(int j=iy_min;j<=iy_max;j++){
		double binY= this->GetYaxis()->GetBinCenter(j+1);
		for(int i=ix_min;i<=ix_max;i++){		
			double binX= this->GetXaxis()->GetBinCenter(i+1);
			double pixValue= this->GetBinContent(i+1,j+1);
			if(pixValue==0) continue; 
			tile->FillPixel(binX,binY,pixValue);	
		}//end loop columns
	}//end loop row
	
	return tile;

}//close Img::GetTile()


Img* Img::GetTile(int tileId){

	//## Check tile id
	int nTiles= fNTilesX*fNTilesY;
  if(tileId<0 || tileId>=nTiles) {
		cerr<<"Img::GetTile(): ERROR: Cannot get tile id "<<tileId<<" (nTiles="<<nTiles<<")!"<<endl;
		return NULL;
	}

	//## Find tile id X Y ranges
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();

	int tileIdX= tileId%fNTilesX;
  int tileIdY = ((tileId-tileIdX)/fNTilesX)%fNTilesY;

	int ix_min= tileIdX*fTileSizeX;
	int ix_max= min(ix_min+fTileSizeX-1,Nx-1);
	
	int iy_min= tileIdY*fTileSizeY;
	int iy_max= min(iy_min+fTileSizeY-1,Ny-1);

	//## Check tile sizes
	if(ix_min<0 || ix_min>=Nx || ix_min>ix_max){
		std::stringstream errMsg;
		errMsg<<"Invalid min tile X given ("<<ix_min<<")";
		errflag= errMsg.str();
		cerr<<"Img::GetTile(): ERROR: "<<errflag<<endl;
		//throw std::invalid_argument(errMsg.str().c_str());
		return NULL;
	}
	if(ix_max<0 || ix_max>=Nx || ix_max<ix_min){
		std::stringstream errMsg;
		errMsg<<"Invalid max tile X given ("<<ix_max<<")";
		errflag= errMsg.str();
		cerr<<"Img::GetTile(): ERROR: "<<errflag<<endl;
		//throw std::invalid_argument(errMsg.str().c_str());
		return NULL;
	}
	if(iy_min<0 || iy_min>=Ny || iy_min>iy_max){
		std::stringstream errMsg;
		errMsg<<"Invalid min tile Y given ("<<iy_min<<")";
		errflag= errMsg.str();
		cerr<<"Img::GetTile(): ERROR: "<<errflag<<endl;
		//throw std::invalid_argument(errMsg.str().c_str());
		return NULL;
	}
	if(iy_max<0 || iy_max>=Ny || iy_max<iy_min){
		std::stringstream errMsg;
		errMsg<<"Invalid max tile Y given ("<<iy_max<<")";
		errflag= errMsg.str();
		cerr<<"Img::GetTile(): ERROR: "<<errflag<<endl;
		//throw std::invalid_argument(errMsg.str().c_str());
		return NULL;
	}
	
	int TileSizeX= ix_max-ix_min+1;
	int TileSizeY= iy_max-iy_min+1;

	double binX_min= this->GetXaxis()->GetBinCenter(ix_min+1);
	double binX_max= this->GetXaxis()->GetBinCenter(ix_max+1);
	double binY_min= this->GetYaxis()->GetBinCenter(iy_min+1);
	double binY_max= this->GetYaxis()->GetBinCenter(iy_max+1);
	
	//## Init tile
	Img* tile= new Img;
	//TString tileName= Form("%s-tile_x%d-%d_y%d-%d",this->GetName(),ix_min,ix_max,iy_min,iy_max);
	//TString tileName= Form("%s_tile%d",this->GetName(),tileId);
	TString tileName= Form("TileImg%d",tileId);
	tile->SetNameTitle(tileName,tileName);
	tile->SetBins(TileSizeX,binX_min-0.5,binX_max+0.5,TileSizeY,binY_min-0.5,binY_max+0.5);
	tile->SetMetaData(this->metadata);
	tile->Reset();
	
	//## Read tile
	cout<<"Img::GetTile(): INFO: Read tile ix("<<ix_min<<","<<ix_max<<") iy("<<iy_min<<","<<iy_max<<") binX("<<binX_min<<","<<binX_max<<")"<<" binY("<<binY_min<<","<<binY_max<<")"<<endl;
	for(int j=iy_min;j<=iy_max;j++){
		double binY= this->GetYaxis()->GetBinCenter(j+1);
		for(int i=ix_min;i<=ix_max;i++){		
			double binX= this->GetXaxis()->GetBinCenter(i+1);
			double pixValue= this->GetBinContent(i+1,j+1);
			if(pixValue==0) continue; 
			tile->FillPixel(binX,binY,pixValue);	
		}//end loop columns
	}//end loop row
	cout<<"Img::GetTile(): INFO: end tile read"<<endl;

	return tile;

}//close Img::GetTile()

Img* Img::MorphFilter(std::string operation,int KernSize,int structElementType) {

	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();

	//## Convert image to OpenCV format
	cv::Mat mat= this->ImgToMat();

	//## Init struct element
	int ptSize= (KernSize-1)/2;
	cv::Size kernel_size(KernSize,KernSize);
	cv::Mat element;
	if(structElementType==1) element= cv::getStructuringElement(cv::MORPH_RECT, kernel_size, cv::Point(ptSize,ptSize));
	else if(structElementType==2) element= cv::getStructuringElement(cv::MORPH_ELLIPSE, kernel_size, cv::Point(ptSize,ptSize));
	else if(structElementType==3) element= cv::getStructuringElement(cv::MORPH_CROSS, kernel_size, cv::Point(ptSize,ptSize));	
	else element= cv::getStructuringElement(cv::MORPH_RECT, kernel_size, cv::Point(ptSize,ptSize));

	//## Morphology operation
	//MORPH_OPEN - an opening operation
	//MORPH_CLOSE - a closing operation
	//MORPH_GRADIENT - a morphological gradient
	//MORPH_TOPHAT - “top hat”
	//MORPH_BLACKHAT - “black hat”
	cv::Mat mat_morph;
	if(operation=="OPENING") cv::morphologyEx(mat, mat_morph, cv::MORPH_OPEN, element);	
	else if(operation=="CLOSING") cv::morphologyEx(mat, mat_morph, cv::MORPH_CLOSE, element);	
	else if(operation=="GRADIENT") cv::morphologyEx(mat, mat_morph, cv::MORPH_GRADIENT, element);	
	else if(operation=="TOPHAT") cv::morphologyEx(mat, mat_morph, cv::MORPH_TOPHAT, element);	
	else if(operation=="BLACKHAT") cv::morphologyEx(mat, mat_morph, cv::MORPH_BLACKHAT, element);	
	else if(operation=="EROSION"){
		cv::erode(mat,mat_morph,element);
	}
	else if(operation=="DILATION"){
		cv::dilate(mat,mat_morph,element);
	}
	else{
		cerr<<"Img::MorphFilter(): ERROR: Invalid morph operation given ("<<operation<<")!"<<endl;
		return 0;
	}
	
	//## Convert back dilated image 
	TString imgName= Form("%s_Morph",std::string(this->GetName()).c_str());	
	Img* MorphImg= (Img*)this->Clone(imgName);
	MorphImg->Reset();

	for(int j=0;j<Ny;j++){
		int rowId= Ny-1-j;
		double y= this->GetYaxis()->GetBinCenter(j+1);
		for(int i=0;i<Nx;i++){
			double x= this->GetXaxis()->GetBinCenter(i+1);
			int colId= i;
			double binContent= this->GetBinContent(i+1,j+1);
			if(binContent==0) continue;
			double matrixElement= mat_morph.at<double>(rowId,colId);			
			MorphImg->SetBinContent(i+1,j+1,matrixElement);
		}//end loop x
	}//end loop y

	return MorphImg;

}//close Img::MorphFilter()


TGraph* Img::FindMorphPeakIds(int KernMin,int KernMax,int KernStep,int peakShiftTolerance,bool isValley){

	//## Init
	
	if(KernMax<KernMin || KernStep<0 || KernStep%2!=0 || KernMin%2==0 || KernMax%2==0) {
		return 0;
	}

	int structElementType= 1;
	std::string morphOp= "TOPHAT";
	if(isValley) morphOp= "BLACKHAT";
	std::vector<int> kernelSizes;
	for(int i=KernMin;i<KernMax;i+=KernStep){
		kernelSizes.push_back(i);
	}
	int nKernels= (int)kernelSizes.size();

	//## Find peaks with different kernel sizes
	bool hasPeaks= true;
	Img* peakImg= (Img*)this->Clone("peakImg");
	peakImg->Reset();

	for(int k=0;k<nKernels;k++){
		//Apply a morph operation
		Img* morphImg= this->MorphFilter(morphOp,kernelSizes[k],structElementType);

		//Peak/valley are pixels surviving the threshold
		Img* morphBinaryImg= morphImg->GetBinarized(0,0,1);

		//Sum peaks over multi-kernels
		peakImg->Add(morphBinaryImg);

		//Clear
		if(morphImg) morphImg->Delete();
		if(morphBinaryImg) morphBinaryImg->Delete();
	}//end loop kernels

	//Threshold peaks/valleys persisting in all kernels 
	double peakThr= nKernels;
	Img* peakBinaryImg= peakImg->GetBinarized(peakThr,0,1);
			
	int nPeaks= 0;
	TGraph* peaks= new TGraph;
	for(int i=0;i<peakBinaryImg->GetNbinsX();i++){
		double x= peakBinaryImg->GetXaxis()->GetBinCenter(i+1);
		for(int j=0;j<peakBinaryImg->GetNbinsY();j++){
			double y= peakBinaryImg->GetYaxis()->GetBinCenter(j+1);
			int gBin= peakBinaryImg->GetBin(i+1,j+1);
			double w= peakBinaryImg->GetBinContent(i+1,j+1);
			if(w>0){
				peaks->SetPoint(nPeaks,x,y);		
				nPeaks++;			
			}
		}//end loop 
	}//end loop

	peakImg->Delete();
	peakBinaryImg->Delete();

	return peaks;

}//close FindMorphPeakIds()


TGraph* Img::FindPeaks(int tol){

	int nKernels= 3;
	//int kernelSizes[]= {5,7,9};
	int kernelSizes[]= {3,5,7};
	std::vector< std::vector<TVector2> > points;  
	for(int k=0;k<nKernels;k++) points.push_back( std::vector<TVector2>() );

	//## Find peaks with different kernel sizes
	bool hasPeaks= true;
	for(int k=0;k<nKernels;k++){
		TGraph* thisPeaks= new TGraph;
		Img* dilatedImg= this->Dilate(*thisPeaks,kernelSizes[k]);
		int nPeaks= thisPeaks->GetN();
		//cout<<"--> nPeaks="<<nPeaks<<endl;
		if(nPeaks<=0){
			hasPeaks= false;
			if(dilatedImg) dilatedImg->Delete();
			if(thisPeaks) thisPeaks->Delete();
			break;
		}
		for(int j=0;j<nPeaks;j++){
			double x, y;
			thisPeaks->GetPoint(j,x,y);
			points[k].push_back(TVector2(x,y));
		}//end loop peak points
		if(dilatedImg) dilatedImg->Delete();
		if(thisPeaks) thisPeaks->Delete();
	}//end loop kernels

	if(!hasPeaks) {
		cout<<"Img::FindPeaks(): INFO: No peaks detected in one or all dilated kernel runs!"<<endl;
		return 0;
	}

	//## Match peaks found with different kernels (given a tolerance)
	TGraph* peaks= new TGraph;
	int npeaks= 0;
	//cout<<"--> npeaks("<<points[0].size()<<","<<points[1].size()<<","<<points[2].size()<<")"<<endl;
	for (int i=0;i<points[0].size();i++) {
  	for (int j=0; j<points[1].size(); j++) {
    	for (int k=0; k<points[2].size(); k++) {
				TVector2 P1= points[0][i];				
				TVector2 P2= points[1][j];     		
				TVector2 P3= points[2][k];
				double dist12_x= fabs(P1.X()-P2.X());
				double dist12_y= fabs(P1.Y()-P2.Y());
				double dist13_x= fabs(P1.X()-P3.X());
				double dist13_y= fabs(P1.Y()-P3.Y());
				double dist23_x= fabs(P2.X()-P3.X());
				double dist23_y= fabs(P2.Y()-P3.Y());

				if( dist12_x<=tol && dist13_x<=tol && dist23_x<=tol &&
						dist12_y<=tol && dist13_y<=tol && dist23_y<=tol
				){
					TVector2 PMean= (P1+P2+P3)*1/3.;
					peaks->SetPoint(npeaks,PMean.X(),PMean.Y());
					npeaks++;			
				}   	
      }//end loop k
    }//end loop j
  }//end loop i

	if(npeaks<=0){
		cout<<"Img::FindPeaks(): INFO: No matching peaks detected!"<<endl;
		if(peaks) peaks->Delete();
		return 0;
	}
	cout<<"Img::FindPeaks(): INFO: #"<<npeaks<<" detected!"<<endl;

	return peaks;

}//close Img::FindPeaks()

std::vector<int> Img::FindPeakIds(int peakShiftTolerance){

	//Init peak ids
	std::vector<int> peakIds;
	peakIds.clear();
	
	//Find peak graph
	TGraph* peakGraph= this->FindPeaks(peakShiftTolerance);
	if(!peakGraph) return peakIds;
	
	int nPeaks= peakGraph->GetN();
	for(int k=0;k<nPeaks;k++){
		double x, y;
		peakGraph->GetPoint(k,x,y);
		int peakBinId= this->FindBin(x,y);
		if(this->IsBinOverflow(peakBinId) || this->IsBinUnderflow(peakBinId) ) continue;
		peakIds.push_back(peakBinId);
	}	

	peakGraph->Delete();

	return peakIds;

}//close FindPeakIds()

Img* Img::Dilate(TGraph& peaks,int KernSize){
		
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();

	//## Convert image to OpenCV format
	cv::Mat mat= this->ImgToMat();

	//## Init dilation options
	cv::Size kernel_size(KernSize,KernSize);
	cv::Mat element= cv::getStructuringElement(cv::MORPH_RECT, kernel_size, cv::Point(-1,-1));
	
	//## Dilate image
	cv::Mat mat_dilated;
	int iterations= 1;
	cv::dilate(mat, mat_dilated, element, cv::Point(-1,-1),iterations,cv::BORDER_CONSTANT);

	//## Compare original and dilated image
	cv::Mat mat_cmp = cv::Mat::zeros(Ny,Nx,CV_8UC1);
	cv::compare(mat, mat_dilated, mat_cmp, cv::CMP_EQ);
	//cv::Mat mat_cmp = mat_dilated = mat;

	//## Convert back dilated image 
	TString imgName= Form("%s_Dilated",std::string(this->GetName()).c_str());	
	Img* DilatedImg= (Img*)this->Clone(imgName);
	DilatedImg->Reset();
	
	int npeaks= 0;
	peaks.Set(0);
	peaks.SetMarkerSize(1.3);
	peaks.SetMarkerColor(kRed);
	peaks.SetMarkerStyle(8);

	for(int j=0;j<Ny;j++){
		int rowId= Ny-1-j;
		double y= this->GetYaxis()->GetBinCenter(j+1);

		for(int i=0;i<Nx;i++){
			double x= this->GetXaxis()->GetBinCenter(i+1);

			int colId= i;
			double binContent= this->GetBinContent(i+1,j+1);
			if(binContent==0) continue;
			double matrixElement= mat_dilated.at<double>(rowId,colId);			
			DilatedImg->SetBinContent(i+1,j+1,matrixElement);
	
			float mat_comparison= (float)mat_cmp.at<unsigned char>(rowId,colId);
			//cout<<"(i,j)=("<<i<<","<<j<<") mat_comparison="<<mat_comparison<<endl;
			if(mat_comparison==255){

				//## Check surrounding pixels (do not tag as peak if the surrounding 3x3 is flat)
				bool isFlatArea= true;
				for(int ix=i-1;ix<i+1;ix++){
					for(int iy=j-1;iy<j+1;iy++){
						if(ix==i && iy==j) continue;
						double w= this->GetBinContent(ix+1,iy+1);
						if(w!=binContent) {
							isFlatArea= false;
							break;
						}
					}//end loop kernel y
				}//end loop kernel x

				if(!isFlatArea){
					peaks.SetPoint(npeaks,x,y);
					npeaks++;
					cout<<"Img::Dilate(): INFO: Peaks #"<<npeaks<<" detected @ ("<<i<<","<<j<<")"<<endl;
				}
			}
		}//end loop x
	}//end loop y

	return DilatedImg;

}//close Img::Dilate()



int Img::DilateSource(Source* source,bool useLocalBackground,int KernSize,bool skipToNested,int sourceType,int dilateModel,bool randomize,double randSigma){

	if(!source || !source->fIsGoodSource) return -1;

	double sigmaTrunc= randSigma;//trunc random gaussian to +-sigmaTrunc	
	RInside* fR= 0;
	std::string randomGenCmd= std::string("rtruncnorm(1, a=-sigmaTrunc, b=sigmaTrunc, mean = 0, sd = 1)");

	try{
		//Get instance to RInside
		fR= RInside::instancePtr();
		if(!fR){
			cerr<<"Img::DilateSource(): ERROR: Cannot retrieve RInside instance (did you create one in your main application?)!"<<endl;
			fR= new RInside;		
		} 
		(*fR)["sigmaTrunc"]= sigmaTrunc;	
		fR->parseEvalQ( std::string("library(\"truncnorm\");") );
	}
	catch( std::exception &ex ) {
		cerr << "Img::DilateSource(): ERROR: Exception catched: " << ex.what() << endl;
		return -1;
  } 
	catch(...) { 
		cerr << "Img::DilateSource(): ERROR: C++ exception (unknown reason)" << endl;
		return -1;
  }	

	BkgData bkgInfo;
	bkgInfo.ix_min= -1;
	bkgInfo.iy_min= -1;
	bkgInfo.ix_max= -1;
	bkgInfo.iy_max= -1;
	bkgInfo.npix= 0;
	bkgInfo.bkgLevel= 0;
	bkgInfo.bkgRMS= 0;
	bkgInfo.isReliableBkg= true;

	Source::PixelCollection sourcePixels= source->GetPixels();
	bool hasNestedSources= source->fHasNestedSources;
	double sourceMean= source->fMean;
	double sourceRMS= source->fRMS;
	double sourceMedian= source->fMedian;
	double sourceMedianRMS= source->fMedianRMS;
		
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	std::vector< std::vector<bool> > isVisited;
	for(int i=0;i<Nx;i++){
		std::vector<bool> v;
		for(int j=0;j<Ny;j++){
			v.push_back(false);			
		}//end loop bins y
		isVisited.push_back(v);
	}//end loop bins x

	std::vector<int> pixelsToBeDilated;
	int dilateSize= KernSize/2;
	cv::Mat element= cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(KernSize,KernSize));
	TMatrixD DilateKernel(KernSize,KernSize);	
	DilateKernel.Zero();
	for(int i=0;i<DilateKernel.GetNrows();i++){
		for(int j=0;j<DilateKernel.GetNcols();j++){
			DilateKernel(i,j)= (double)element.at<char>(i,j); 
		}
	}
	//cout<<"Img::DilateSource(): INFO: Dilate kernel"<<endl;
	//DilateKernel.Print();

	//## Mother source is dilated replacing image pixel values with random image background
	//if(sourceType==-1 || source->fType==sourceType){
	if(!hasNestedSources && (sourceType==-1 || source->fType==sourceType) ){
	 cout<<"Img::DilateSource(): INFO: Dilating mother source..."<<endl;

		for(unsigned int l=0;l<sourcePixels.size();l++){
			Source::Pixel thisPixel= sourcePixels[l];
			int id= thisPixel.id;
			int binx, biny, binz;
			this->GetBinXYZ(id,binx,biny,binz);
			int ix= binx-1;
			int iy= biny-1;			

			std::vector<int>::iterator it= std::find(pixelsToBeDilated.begin(),pixelsToBeDilated.end(),id);
			if( (it==pixelsToBeDilated.end() || pixelsToBeDilated.empty()) ) pixelsToBeDilated.push_back(id);
						
			for(int tx=-dilateSize;tx<=dilateSize;tx++){
				int binx_next= binx+tx;
				int colId= tx + dilateSize;
				for(int ty=-dilateSize;ty<=dilateSize;ty++){	
					int biny_next= biny+ty;
					int rowId= tx + dilateSize;
					if(ix+tx<Nx && ix+tx>=0 && iy+ty<Ny && iy+ty>=0){
						double kernValue= DilateKernel(rowId,colId);
						int gBinId= this->GetBin(binx_next,biny_next);
						std::vector<int>::iterator it= std::find(pixelsToBeDilated.begin(),pixelsToBeDilated.end(),gBinId);
						if( kernValue>0 && (it==pixelsToBeDilated.end() || pixelsToBeDilated.empty()) ) pixelsToBeDilated.push_back(gBinId);
					}
				}//end loop kernel
			}//end loop kernel
		}//end loop pixels

		/*
		for(unsigned int l=0;l<sourcePixels.size();l++){
			Source::Pixel thisPixel= sourcePixels[l];
			int id= thisPixel.id;
			double S= thisPixel.S;
			double X= thisPixel.x;
			double Y= thisPixel.y;
			int ix= thisPixel.ix;
			int iy= thisPixel.iy;
			if(ix>=Nx || iy>=Ny) {
				cout<<"Img::DilateSource(): INFO: Skipping source pixel ("<<ix<<","<<iy<<") as invalid (Nx,Ny)=("<<Nx<<","<<Ny<<")"<<endl;
				continue;
			}

			int binx, biny, binz;
			this->GetBinXYZ(id,binx,biny,binz);
				
			//int status= this->GetPixelBkg(ix,iy,bkgInfo,useLocalBackground);
			//if(status<0){
			//	cerr<<"Img::DilateSource(): WARN: Cannot get bkg for this pixel ("<<ix<<","<<iy<<"), skip it!"<<endl;
			//	continue;
			//}
			//double BkgRealization= gRandom->Gaus(bkgInfo.bkgLevel,bkgInfo.bkgRMS);
			//double nRandSigmas= Rcpp::as<double>( fR->parseEval(randomGenCmd) );
			//double BkgRealization= bkgInfo.bkgLevel + nRandSigmas*bkgInfo.bkgRMS;

			if(isVisited[ix][iy]) continue;
				
			pixelsToBeDilated.push_back(id);
			//this->SetBinContent(ix+1,iy+1,BkgRealization);
			isVisited[ix][iy]= true;	
		
			//Enlarge residual region (dilation) with a given kernel size
			for(int tx=-KernSize;tx<=KernSize;tx++){
				//if(tx==0) continue;
				int binx_next= binx+tx;
				for(int ty=-KernSize;ty<=KernSize;ty++){
					//if(ty==0) continue;
					int biny_next= biny+ty;
					if(ix+tx<Nx && ix+tx>=0 && iy+ty<Ny && iy+ty>=0 && !isVisited[ix+tx][iy+ty]){
						//status= GetPixelBkg(ix+tx,iy+ty,bkgInfo,useLocalBackground);
						//nRandSigmas= Rcpp::as<double>( fR->parseEval(randomGenCmd) );
						//BkgRealization= bkgInfo.bkgLevel + nRandSigmas*bkgInfo.bkgRMS;

						int gBinId= this->GetBin(binx_next,biny_next);
						std::vector<int>::iterator it= std::find(pixelsToBeDilated.begin(),pixelsToBeDilated.end(),gBinId);
						if(it==pixelsToBeDilated.end() || pixelsToBeDilated.empty()) pixelsToBeDilated.push_back(gBinId);
						//this->SetBinContent(binx_next,biny_next,BkgRealization);
						isVisited[ix+tx][iy+ty]= true;
					}//close if
				}//end dilation loop y
			}//end dilation loop x
		}//end loop pixels
		*/

	}//close if


	//## Loop over nested source pixels and dilate them
	

	if(skipToNested){
		//## Only nested source are dilated, all if sourceType==-1 or selected types according to sourceType flag
		//## Pixel values are replaced with mother source mean&rms NOT with image background
		cout<<"Img::DilateSource(): INFO: Dilating nested sources..."<<endl;

		//Get nested source list
		std::vector<Source*> nestedSources= source->fNestedSourceCollection;

		if(source->fHasNestedSources || nestedSources.size()>0) {
			
			for(unsigned int k=0;k<nestedSources.size();k++){
				int nestedSourceType= nestedSources[k]->fType;
				if(sourceType==-1 || nestedSourceType==sourceType){
					Source::PixelCollection nestedSourcePixels= nestedSources[k]->GetPixels();
					for(unsigned int l=0;l<nestedSourcePixels.size();l++){
						Source::Pixel thisPixel= nestedSourcePixels[l];
						int id= thisPixel.id;
						int binx, biny, binz;
						this->GetBinXYZ(id,binx,biny,binz);	
						int ix= binx-1;
						int iy= biny-1;			

						std::vector<int>::iterator it= std::find(pixelsToBeDilated.begin(),pixelsToBeDilated.end(),id);
						if( (it==pixelsToBeDilated.end() || pixelsToBeDilated.empty()) ) pixelsToBeDilated.push_back(id);
						
						for(int tx=-dilateSize;tx<=dilateSize;tx++){
							int binx_next= binx+tx;
							int colId= tx + dilateSize;
							for(int ty=-dilateSize;ty<=dilateSize;ty++){	
								int biny_next= biny+ty;
								int rowId= tx + dilateSize;
								if(ix+tx<Nx && ix+tx>=0 && iy+ty<Ny && iy+ty>=0){
									double kernValue= DilateKernel(rowId,colId);
									int gBinId= this->GetBin(binx_next,biny_next);
									std::vector<int>::iterator it= std::find(pixelsToBeDilated.begin(),pixelsToBeDilated.end(),gBinId);
									if( kernValue>0 && (it==pixelsToBeDilated.end() || pixelsToBeDilated.empty()) ) pixelsToBeDilated.push_back(gBinId);
								}
							}//end loop kernel
						}//end loop kernel

					}//end loop pixels
				}//close if source type
			}//end loop nested sources

			/*
			for(unsigned int k=0;k<nestedSources.size();k++){
				int nestedSourceType= nestedSources[k]->fType;
			
				if(sourceType==-1 || nestedSourceType==sourceType){
					Source::PixelCollection nestedSourcePixels= nestedSources[k]->GetPixels();
					for(unsigned int l=0;l<nestedSourcePixels.size();l++){
						Source::Pixel thisPixel= nestedSourcePixels[l];
						double S= thisPixel.S;
						double X= thisPixel.x;
						double Y= thisPixel.y; 
						int id= this->FindBin(X,Y);
						int binx, biny, binz;
						this->GetBinXYZ(id,binx,biny,binz);
						int ix= binx-1;
						int iy= biny-1;					
			 
						if(ix>=Nx || iy>=Ny) {
							cout<<"Img::DilateSource(): INFO: Skipping nested source pixel ("<<ix<<","<<iy<<") as invalid (Nx,Ny)=("<<Nx<<","<<Ny<<")"<<endl;
							continue;
						}

						if(isVisited[ix][iy]) continue;
		
						pixelsToBeDilated.push_back(id);

						//this->SetBinContent(ix+1,iy+1,BkgRealization);
						isVisited[ix][iy]= true;	
		
						//Enlarge residual region (dilation) with a given kernel size
						int dilateSize= KernSize/2;

						//for(int tx=-KernSize;tx<=KernSize;tx++){
						for(int tx=-dilateSize;tx<=dilateSize;tx++){
							int binx_next= binx+tx;
							int colId= tx + dilateSize;

						//for(int ty=-KernSize;ty<=KernSize;ty++){
						for(int ty=-dilateSize;ty<=dilateSize;ty++){	
							int biny_next= biny+ty;
							int rowId= tx + dilateSize;
							if(ix+tx<Nx && ix+tx>=0 && iy+ty<Ny && iy+ty>=0 && !isVisited[ix+tx][iy+ty]){
								//--> Random source mean
								//nRandSigmas= Rcpp::as<double>( fR->parseEval(randomGenCmd) );
								//BkgRealization= sourceMedian + nRandSigmas*sourceMedianRMS;
								//BkgRealization= sourceMedian;
								
								//--> Random bkg
								//status= GetPixelBkg(ix+tx,iy+ty,bkgInfo,useLocalBackground);
								//nRandSigmas= Rcpp::as<double>( fR->parseEval(randomGenCmd) );
								//BkgRealization= bkgInfo.bkgLevel + nRandSigmas*bkgInfo.bkgRMS;
								//BkgRealization= bkgInfo.bkgLevel + 1*bkgInfo.bkgRMS;
									
									double kernValue= DilateKernel(rowId,colId);
									int gBinId= this->GetBin(binx_next,biny_next);
									std::vector<int>::iterator it= std::find(pixelsToBeDilated.begin(),pixelsToBeDilated.end(),gBinId);
									if( kernValue>0 && (it==pixelsToBeDilated.end() || pixelsToBeDilated.empty()) ) pixelsToBeDilated.push_back(gBinId);
									//this->SetBinContent(binx_next,biny_next,BkgRealization);
									isVisited[ix+tx][iy+ty]= true;
								}//close if
							}//end dilation loop y
						}//end dilation loop x
					}//end loop pixels
				}//close if source type
			}//end loop nested sources
			*/
		}//close if
		else{
			cout<<"Img::DilateSource(): INFO: No nested source present, nothing will be dilated in this step"<<endl;
		}
	}//close if


	//## Replace selected pixels
	cout<<"Img::DilateSource(): INFO: #"<<pixelsToBeDilated.size()<<" pixels to be dilated..."<<endl;
	for(unsigned int l=0;l<pixelsToBeDilated.size();l++){
		int id= pixelsToBeDilated[l];			
		double BkgRealization= this->GetBinContent(id);
		double BkgRMS= 0;
		if(dilateModel==eSourceMedian){
			BkgRealization= sourceMedian;
			BkgRMS= sourceMedianRMS;
		}
		else if(dilateModel==eSourceMean){
			BkgRealization= sourceMean;
			BkgRMS= sourceRMS;
		}
		else if(dilateModel==eBkg){
			int binx, biny, binz;
			this->GetBinXYZ(id,binx,biny,binz);
			int ix= binx-1;
			int iy= biny-1;		
			int status= this->GetPixelBkg(ix,iy,bkgInfo,useLocalBackground);
			if(status<0){
				cerr<<"Img::DilateSource(): WARN: Cannot get bkg for this pixel ("<<ix<<","<<iy<<"), skip it!"<<endl;
				continue;
			}
			BkgRealization= bkgInfo.bkgLevel;
			BkgRMS= bkgInfo.bkgRMS;
		}
		
		if(randomize) {
			double nRandSigmas= Rcpp::as<double>( fR->parseEval(randomGenCmd) );
			BkgRealization+= nRandSigmas*BkgRMS;
		}

		this->SetBinContent(id,BkgRealization);
	}//end loop dilated pixels


	return 0;

}//close Img::DilateSource()


std::vector<double> Img::GetZernikeMoments(int order,double radius,int method){

	std::vector<double>	zm;
	zm.clear();
	if(method==1) zm= ZernikeMoments::GetZernike2D_Direct(this,order,radius);
	else if(method==2) zm= ZernikeMoments::GetZernike2D(this,order,radius);
	
	return zm;

}//close GetZernikeMoments()

std::vector<double> Img::GetHuMoments(){

	//## Convert to OpenCV
	cv::Mat mat= this->ImgToMat("64");

	//## Compute moments
	cv::Moments moments= cv::moments(mat, false);
		
	//## Compute HuMoments
	double HuMoments[7];
	cv::HuMoments(moments, HuMoments);

	std::vector<double>	hm;
	hm.clear();
	for(int i=0;i<7;i++) hm.push_back(HuMoments[i]);	

	return hm;

}//close GetHuMoments()


Img* Img::GetSourceResidual(bool useLocalBackground,int KernSize,bool skipToNested,int sourceType,int dilateModel,bool randomize,double randSigma){

	if(!fHasSources) {
		cerr<<"Img::GetSourceResidual(): INFO: No source available!"<<endl;
		return this;
	}
	if(!this->HasBkgData() || fBkgData.size()<=0) {
		cerr<<"Img::GetSourceResidual(): INFO: No background available!"<<endl;
		return this;
	}

	//## Clone map
	double pixelMin= this->GetMinimum();
	double pixelMax= this->GetMaximum();
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();

	TString imgName= Form("%s_SourceResidual",std::string(this->GetName()).c_str());

	//## Clone image
	Img* SourceResidual= (Img*)this->Clone(imgName);
	if(this->fHasMetaData) SourceResidual->SetMetaData(this->metadata);
	if(this->HasBkgData()) {
		SourceResidual->SetBkgData(this->fBkgData);
		SourceResidual->SetInterpolatedBkgLevelMap(this->fInterpolatedBackgroundLevelMap);
		SourceResidual->SetInterpolatedBkgRMSMap(this->fInterpolatedBackgroundRMSMap);
	}
		
	//cout<<"HasBkgData? "<<this->HasBkgData()<<", HasBkgDataSourceResidual? "<<SourceResidual->HasBkgData()<<endl;

	cout<<"Img::GetSourceResidual(): INFO: Dilating source pixels ("<<fSourceCollection.size()<<" sources available)"<<endl;
	for(unsigned int k=0;k<fSourceCollection.size();k++){
		Source* thisSource= fSourceCollection[k];
		int status= SourceResidual->DilateSource(thisSource,useLocalBackground,KernSize,skipToNested,sourceType,dilateModel,randomize,randSigma);
		if(status<0) {
			cerr<<"Img::GetSourceResidual(): WARN: Dilating source no. "<<k<<" failed!"<<endl;
			continue;
		}
	}//end loop sources

	//## Before returning residual image, clear background info and sources	
	SourceResidual->ResetBkg();
	SourceResidual->ResetSources();

	//SourceResidual->SetMinimum(pixelMin);
	//SourceResidual->SetMaximum(pixelMax);

	return SourceResidual;

}//close GetSourceResidual()

Img* Img::GetMask(Img* mask,bool isBinary){

	if(!mask) return 0;
		
	//## Check mask bins
	int Nx= mask->GetNbinsX();
	int Ny= mask->GetNbinsY();
	if(Nx!=this->GetNbinsX() || Ny!=this->GetNbinsY()){
		cerr<<"Img::GetMask(): ERROR: Mask binning is different than current image!"<<endl;
		return 0;
	}	

	//## Clone map
	TString imgName= Form("%s_Mask",std::string(this->GetName()).c_str());
	Img* maskedImage= new Img(imgName,imgName,this->GetNbinsX(),this->GetXaxis()->GetXmin(),this->GetXaxis()->GetXmax(),this->GetNbinsY(),this->GetYaxis()->GetXmin(),this->GetYaxis()->GetXmax());
	if(this->fHasMetaData) maskedImage->SetMetaData(this->metadata);
	if(this->HasBkgData()) {
		maskedImage->SetBkgData(this->fBkgData);
		maskedImage->SetInterpolatedBkgLevelMap(this->fInterpolatedBackgroundLevelMap);
		maskedImage->SetInterpolatedBkgRMSMap(this->fInterpolatedBackgroundRMSMap);
	}

	//## Loop over mask	
	if(isBinary){
		for(int i=0;i<mask->GetNbinsX();i++){
			for(int j=0;j<mask->GetNbinsY();j++){
				double maskContent= mask->GetBinContent(i+1,j+1);
				if(maskContent>0) {
					maskedImage->SetBinContent(i+1,j+1,1);
					continue;
				}
				maskedImage->SetBinContent(i+1,j+1,0);
			}//end loop bins Y
		}//end loop bins X
	}
	else{
		for(int i=0;i<mask->GetNbinsX();i++){
			for(int j=0;j<mask->GetNbinsY();j++){
				double binContent= this->GetBinContent(i+1,j+1);
				double maskContent= mask->GetBinContent(i+1,j+1);
				if(maskContent>0) {
					maskedImage->SetBinContent(i+1,j+1,binContent);
					continue;
				}
				maskedImage->SetBinContent(i+1,j+1,0);
			}//end loop bins Y
		}//end loop bins X
	}//close else

	return maskedImage;

}//close Img::GetMask()



Img* Img::GetSourceMask(std::vector<Source*> sources,bool isBinary){

	//## Clone map
	TString imgName= Form("%s_SourceMask",std::string(this->GetName()).c_str());
	Img* SourceMask= (Img*)this->Clone(imgName);
	//Img* SourceMask= new Img(imgName,imgName,this->GetNbinsX(),this->GetXaxis()->GetXmin(),this->GetXaxis()->GetXmax(),this->GetNbinsY(),this->GetYaxis()->GetXmin(),this->GetYaxis()->GetXmax());
	
	//## Copy metadata and background info
	if(this->fHasMetaData) SourceMask->SetMetaData(this->metadata);
	if(this->HasBkgData()) {
		SourceMask->SetBkgData(this->fBkgData);
		if(this->fInterpolatedBackgroundLevelMap) SourceMask->SetInterpolatedBkgLevelMap(this->fInterpolatedBackgroundLevelMap);
		if(this->fInterpolatedBackgroundRMSMap) SourceMask->SetInterpolatedBkgRMSMap(this->fInterpolatedBackgroundRMSMap);
	}

	int nSources= (int)sources.size();
	if(nSources<=0) {
		cerr<<"Img::GetSourceMask(): WARN: Source list is empty, returning empty image!"<<endl;
		return SourceMask;	
	}

	//## Loop over sources
	SourceMask->Reset();

	if(isBinary){
		for(int k=0;k<nSources;k++){
			Source* thisSource= sources[k];
			Source::PixelCollection thisSourcePixels= thisSource->GetPixels();
			for(unsigned int l=0;l<thisSourcePixels.size();l++){
				Source::Pixel thisPixel= thisSourcePixels[l];
				int id= thisPixel.id;
				SourceMask->SetBinContent(id,1);
			}//end loop pixels
		}//end loop sources		
	}
	else{
		for(int k=0;k<nSources;k++){
			Source* thisSource= sources[k];
			Source::PixelCollection thisSourcePixels= thisSource->GetPixels();
			for(unsigned int l=0;l<thisSourcePixels.size();l++){
				Source::Pixel thisPixel= thisSourcePixels[l];
				int id= thisPixel.id;
				double binContent= this->GetBinContent(id);				
				SourceMask->SetBinContent(id,binContent);
			}//end loop pixels
		}//end loop sources		
	}//close else
	
	return SourceMask;

}//close Img::GetSourceMask()

Img* Img::GetSourceMask(std::vector<Source*> sources,int mode){

	//## Clone map
	TString imgName= Form("%s_SourceMask",std::string(this->GetName()).c_str());
	Img* SourceMask= (Img*)this->Clone(imgName);
	//Img* SourceMask= new Img(imgName,imgName,this->GetNbinsX(),this->GetXaxis()->GetXmin(),this->GetXaxis()->GetXmax(),this->GetNbinsY(),this->GetYaxis()->GetXmin(),this->GetYaxis()->GetXmax());
	
	//## Copy metadata and background info
	if(this->fHasMetaData) SourceMask->SetMetaData(this->metadata);
	if(this->HasBkgData()) {
		SourceMask->SetBkgData(this->fBkgData);
		if(this->fInterpolatedBackgroundLevelMap) SourceMask->SetInterpolatedBkgLevelMap(this->fInterpolatedBackgroundLevelMap);
		if(this->fInterpolatedBackgroundRMSMap) SourceMask->SetInterpolatedBkgRMSMap(this->fInterpolatedBackgroundRMSMap);
	}

	int nSources= (int)sources.size();
	if(nSources<=0) {
		cerr<<"Img::GetSourceMask(): WARN: Source list is empty, returning empty image!"<<endl;
		return SourceMask;	
	}

	//## Loop over sources
	SourceMask->Reset();

	if(mode==1){
		for(int k=0;k<nSources;k++){
			Source* thisSource= sources[k];
			Source::PixelCollection thisSourcePixels= thisSource->GetPixels();
			for(unsigned int l=0;l<thisSourcePixels.size();l++){
				Source::Pixel thisPixel= thisSourcePixels[l];
				int id= thisPixel.id;
				SourceMask->SetBinContent(id,1);
			}//end loop pixels
		}//end loop sources		
	}
	else if(mode==2){
		for(int k=0;k<nSources;k++){
			Source* thisSource= sources[k];
			Source::PixelCollection thisSourcePixels= thisSource->GetPixels();
			for(unsigned int l=0;l<thisSourcePixels.size();l++){
				Source::Pixel thisPixel= thisSourcePixels[l];
				int id= thisPixel.id;
				double binContent= this->GetBinContent(id);				
				SourceMask->SetBinContent(id,binContent);
			}//end loop pixels
		}//end loop sources		
	}//close else
	else if(mode==3){
		for(int k=0;k<nSources;k++){
			Source* thisSource= sources[k];
			double sourceMean= thisSource->fMean;
			Source::PixelCollection thisSourcePixels= thisSource->GetPixels();
			for(unsigned int l=0;l<thisSourcePixels.size();l++){
				Source::Pixel thisPixel= thisSourcePixels[l];
				int id= thisPixel.id;
				SourceMask->SetBinContent(id,sourceMean);
			}//end loop pixels
		}//end loop sources		
	}//close else
	
	return SourceMask;

}//close Img::GetSourceMask()

Img* Img::GetSourceMap(bool isBinary){

	if(!fHasSources) return 0;

	//## Clone map
	double pixelMin= this->GetMinimum();
	double pixelMax= this->GetMaximum();

	//Img* SourceMask= (Img*)this->Clone("SourceMask");
	//SourceMask->Reset();

	TString imgName= Form("%s_SourceMask",std::string(this->GetName()).c_str());
	Img* SourceMask= new Img(imgName,imgName,this->GetNbinsX(),this->GetXaxis()->GetXmin(),this->GetXaxis()->GetXmax(),this->GetNbinsY(),this->GetYaxis()->GetXmin(),this->GetYaxis()->GetXmax());
	SourceMask->SetMetaData(this->metadata);

	//## Loop over sources
	if(isBinary){
		for(unsigned int k=0;k<fSourceCollection.size();k++){
			Source* thisSource= fSourceCollection[k];
			Source::PixelCollection thisSourcePixels= thisSource->GetPixels();
	
			for(unsigned int l=0;l<thisSourcePixels.size();l++){
				Source::Pixel thisPixel= thisSourcePixels[l];
				int id= thisPixel.id;
				double S= thisPixel.S;
				double X= thisPixel.x;
				double Y= thisPixel.y;
				int binx, biny, binz;
				this->GetBinXYZ(id,binx,biny,binz);
				SourceMask->SetBinContent(id,1);
			}//end loop pixels
		}//end loop sources		
	}
	else{
		for(unsigned int k=0;k<fSourceCollection.size();k++){
			Source* thisSource= fSourceCollection[k];
			Source::PixelCollection thisSourcePixels= thisSource->GetPixels();
	
			for(unsigned int l=0;l<thisSourcePixels.size();l++){
				Source::Pixel thisPixel= thisSourcePixels[l];
				int id= thisPixel.id;
				double S= thisPixel.S;				
				double X= thisPixel.x;
				double Y= thisPixel.y;
				int binx, biny, binz;
				this->GetBinXYZ(id,binx,biny,binz);
				double binContent= this->GetBinContent(id);
				//SourceMask->FillPixel(binx-1,biny-1,S);
				//SourceMask->FillPixel(X,Y,S);
				//SourceMask->Fill(X,Y,S);
				//SourceMask->Fill(X,Y,binContent);
				SourceMask->SetBinContent(id,binContent);
			}//end loop pixels
		}//end loop sources		
		//SourceMask->SetMinimum(pixelMin);
		//SourceMask->SetMaximum(pixelMax);

	}//close else
	
	return SourceMask;

}//close Img::GetSourceMap()


TMatrixD* Img::GetMatrix(){

	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
	
	//TMatrixD* M= new TMatrixD(Nx,Ny);
	TMatrixD* M= new TMatrixD(Ny,Nx);
	M->Zero();

	for(int i=0;i<Nx;i++){//rows
		for(int j=0;j<Ny;j++){//columns
			double w= this->GetBinContent(i+1,j+1);
			//int rowId= Ny-j-1;
			//int colId= i;

			int rowId= j;
			int colId= i;

			//M(i,j)= w;
			(*M)(rowId,colId)= w;
		}//end loop cols
	}//end loop rows

	return M;

}//close Img::GetMatrix()


ImgPtrCollection Img::GetWaveletDecomposition(int nScales){

	cout<<"Img::GetWaveletDecomposition(): INFO: Computing wavelet decomposition up to scale J="<<nScales<<" ..."<<endl;
	ImgPtrCollection img_decomposition;
	img_decomposition= WTFilter::GetDecomposition(this,nScales);
	
	return img_decomposition;

}//close GetWaveletDecomposition()


Img* Img::Smooth(int size_x,int size_y,double sigma_x,double sigma_y){

	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();
		
	//## Fill OpenCV mat
	cv::Mat mat = cv::Mat::zeros(Ny,Nx,CV_64FC1);

	for(int j=0;j<Ny;j++){
		int rowId= Ny-1-j;
		for(int i=0;i<Nx;i++){
			int colId= i;
			mat.at<double>(rowId,colId)= this->GetBinContent(i+1,j+1);
		}//end loop x
	}//end loop y
	
	//## Smooth matrix
	cv::Mat smoothed_mat;
	cv::Size smooth_size(size_x,size_y);
	cv::GaussianBlur(mat,smoothed_mat, smooth_size, sigma_x, sigma_y, cv::BORDER_DEFAULT);

	//## Fill smoothed image
	TString imgName= Form("%s_Smoothed",std::string(this->GetName()).c_str());
	//Img* SmoothedImg= new Img(imgName,imgName,Nx,this->GetXaxis()->GetXmin(),this->GetXaxis()->GetXmax(),Ny,this->GetYaxis()->GetXmin(),this->GetYaxis()->GetXmax());
	//SmoothedImg->Reset();	
	Img* SmoothedImg= (Img*)this->Clone(imgName);
	SmoothedImg->Reset();
	
	for(int j=0;j<Ny;j++){
		int rowId= Ny-1-j;
		for(int i=0;i<Nx;i++){
			int colId= i;
			double matrixElement= smoothed_mat.at<double>(rowId,colId);
			SmoothedImg->SetBinContent(i+1,j+1,matrixElement);
		}//end loop x
	}//end loop y


	return SmoothedImg;

}//close Smooth()


std::vector<Contour*> Img::FindContour(){
	
	Contour* aContour= 0;
	std::vector<Contour*> ContourCollection;
	ContourCollection.clear();
	int Nx= this->GetNbinsX();
	int Ny= this->GetNbinsY();

	//## Fill image and binarized image
	int nRows= Ny;
	int nCols= Nx;
	cv::Mat binarizedImg = cv::Mat::zeros(nRows, nCols, CV_8UC1);
	
	for(int j=0;j<Ny;j++){
		int rowId= Ny-1-j;
		for(int i=0;i<Nx;i++){
			int colId= i;
			double binContent= this->GetBinContent(i+1,j+1);
			if(binContent!=0) binarizedImg.at<uchar>(rowId, colId, 0) = 1;
		}//end loop x
	}//end loop y

	//## Compute contour
	std::vector<std::vector<cv::Point>> contours; // Vector for storing contour
  std::vector<cv::Vec4i> hierarchy;
	cv::findContours( binarizedImg, contours, hierarchy,CV_RETR_EXTERNAL,CV_CHAIN_APPROX_NONE, cv::Point(0,0) ); // Find only the external contours in the image

	int nContours= (int)contours.size();
	if(nContours<=0) {
		cerr<<"Img::FindContour(): WARN: Cannot find contour!"<<endl;
		return ContourCollection;
	}
	cout<<"Img::FindContour(): INFO: "<<nContours<<" contours found!"<<endl;

	//## Return external contours
	for(int k=0;k<nContours;k++){
		int nContourPts= (int)contours[k].size();
		if(nContourPts<=0) {
			cerr<<"Img::FindContour(): WARN: No points in contour no. "<<k<<"!"<<endl;
			continue;
		}

		aContour= new Contour;

		cout<<"Img::FindContour(): INFO: Contour (";
		for(int j=0;j<nContourPts;j++){
			int contx= contours[k][j].x;
			int conty= contours[k][j].y;
			int rowId= Ny-1-conty;
			int colId= contx;	
			double contx_transf= this->GetXaxis()->GetBinCenter(colId+1);
			double conty_transf= this->GetYaxis()->GetBinCenter(rowId+1);
	
			aContour->AddPoint( cv::Point2f(contx_transf,conty_transf) );
			cout<<"("<<contx_transf<<","<<conty_transf<<"), ";
		}//end loop points in contour
		cout<<")"<<endl;
	
		ContourCollection.push_back(aContour);
	}//end loop contours

	return ContourCollection;	

}//close FindContour()



