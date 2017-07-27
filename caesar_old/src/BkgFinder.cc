/**
* @file BkgFinder.cc
* @class BkgFinder
* @brief BkgFinder
*
* BkgFinder class
* @author S. Riggi
* @date 20/01/2015
*/

#include <BkgFinder.h>
#include <VLSlicSegmentation.h>

#include <linterp.h>

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

ClassImp(BkgFinder)


BkgFinder::BkgFinder(){

	fR= 0;
	//fOutputFile= 0;
	//fDataTree= 0;

	fRegionSize= 10;
	fRegularization= 100;
	fMinRegionArea= 5;
	fColorEps= 5;

	/*
	fBackgroundMap= 0;
	fRMSMap= 0;
	fInterpolatedBackgroundMap= 0;
	fInterpolatedRMSMap= 0;
	*/
}//close costructor


BkgFinder::~BkgFinder(){

	cout<<"BkgFinder::~BkgFinder(): INFO: Delete bkg finder..."<<endl;
	

	//if(fDataTree) fDataTree->Delete();
	//if(fOutputFile && fOutputFile->IsOpen()) fOutputFile->Close();
	
	//if(fBackgroundMap) fBackgroundMap->Delete();
	//if(fRMSMap) fRMSMap->Delete();
	
}//close destructor


void BkgFinder::Init(){

	//Init RInside 
	fR= RInside::instancePtr();
	if(!fR){
		cerr<<"BkgFinder::Init(): ERROR: Cannot retrieve RInside instance (did you create one in your main application?)!"<<endl;
		fR= new RInside;
	}	

	//Init data
	//fBkgData.clear();
	

}//close Init()


//int BkgFinder::FindGridBkg(Img* img, Img::BkgMethod method, int boxSizeX, int boxSizeY, double gridStepSizeX, double gridStepSizeY){
BkgFinder::BkgMapData* BkgFinder::FindGridBkg(Img* img, Img::BkgMethod method, int boxSizeX, int boxSizeY, double gridStepSizeX, double gridStepSizeY){

	//if(!img) return -1;
	if(!img) return 0;

	//## Check args
	int Nx= img->GetNbinsX();
	int Ny= img->GetNbinsY();
	if(boxSizeX<=0 || boxSizeX>=Nx || boxSizeY<=0 || boxSizeY>=Ny) {
		cerr<<"BkgFinder::FindGridBkg(): ERROR: Invalid box size given (too small or larger than image size)!"<<endl;
		//return -1;
		return 0;
	}
	if(gridStepSizeX<=0 || gridStepSizeY<=0 ){
		cerr<<"BkgFinder::FindGridBkg(): ERROR: Invalid grid step size given (null or negative)!"<<endl;
		//return -1;
		return 0;
	}
	
	int TileSizeX= boxSizeX;
	int TileSizeY= boxSizeY;
	
	cout<<"BkgFinder::FindGridBkg(): INFO: N("<<Nx<<","<<Ny<<") TileSize=("<<TileSizeX<<","<<TileSizeY<<") GridStepSize("<<gridStepSizeX<<","<<gridStepSizeY<<")"<<endl;

	
	BkgFinder::BkgMapData* bkgMapData= 0;

	try{

		//## Initialize data
		Init();
	
		bkgMapData= new BkgMapData;
		
		//## Count number of rows & cols
		int indexX_start= TileSizeX/2;
		int indexY_start= TileSizeY/2;
		int indexX_end= Nx-1-TileSizeX/2;
		int indexY_end= Ny-1-TileSizeY/2;
		int indexX= indexX_start;
		int indexY= indexY_start;
		int nTiles= 0;
		int nTilesX= 0;
		int nTilesY= 0;
		std::vector<int> ixList;
		std::vector<int> iyList;
		cout<<"BkgFinder::FindGridBkg(): INFO: Counting number of tiles from ("<<indexX<<","<<indexY<<") up to ("<<indexX_end<<","<<indexY_end<<")..."<<endl;

		while(indexY<=indexY_end){
			iyList.push_back(indexY);
			nTilesY++;
			indexY+= gridStepSizeY;
		}//end while loop

		while(indexX<=indexX_end){
			ixList.push_back(indexX);	
			nTilesX++;
			indexX+= gridStepSizeX;
		}
		nTiles= nTilesX*nTilesY;

		cout<<"BkgFinder::FindGridBkg(): INFO: nTilesX="<<nTilesX<<"("<<ixList.size()<<") nTilesY="<<nTilesY<<"("<<iyList.size()<<") nTiles="<<nTiles<<endl;


		//## Loop over all number of tiles and compute bkg info for each one
		Img* TileImg= 0;
		Img::StatsData* stats= 0;
		Img::BkgData* TileImgBkg= 0;
		Rcpp::NumericVector sampledX(nTilesX);
		Rcpp::NumericVector sampledY(nTilesY);
		Rcpp::NumericMatrix sampledBkgLevel(nTilesX,nTilesY);
		Rcpp::NumericMatrix sampledRMS(nTilesX,nTilesY);
		Rcpp::NumericVector xlim(2);//vector of length 2 giving lower and upper limit for range x coordinates used for output grid
		Rcpp::NumericVector ylim(2);//vector of length 2 giving lower and upper limit for range of y coordinates used for output grid
		xlim(0)= img->GetXaxis()->GetXmin();
		xlim(1)= img->GetXaxis()->GetXmax();
		ylim(0)= img->GetYaxis()->GetXmin();
		ylim(1)= img->GetYaxis()->GetXmax();
		int counter= 0;

		std::vector<double> xBins_min;
		std::vector<double> xBins_max;
		std::vector<double> yBins_min;
		std::vector<double> yBins_max;

		std::vector<double> sampledGridX(nTilesX,0);
		std::vector<double> sampledGridY(nTilesY,0);
	
		
		cout<<"BkgFinder::FindGridBkg(): INFO: Computing tile background..."<<endl;

		for(int i=0;i<nTilesX;i++){
			int ix= ixList[i];
			int ix_min= ix-TileSizeX/2;
			int ix_max= ix_min+TileSizeX-1;
			double x= img->GetXaxis()->GetBinCenter(ix+1);
			double x_min= img->GetXaxis()->GetBinCenter(ix_min+1);
			double x_max= img->GetXaxis()->GetBinCenter(ix_max+1);
			double dx_min= img->GetXaxis()->GetBinWidth(ix_min+1);
			double dx_max= img->GetXaxis()->GetBinWidth(ix_max+1);	
			xBins_min.push_back(x_min-0.5*dx_min);
			xBins_max.push_back(x_max+0.5*dx_max);
			
			sampledX(i)= x;	
			sampledGridX[i]= x;

			for(int j=0;j<nTilesY;j++){
				counter++;
				int iy= iyList[j];
				int iy_min= iy-TileSizeY/2;
				int iy_max= iy_min+TileSizeY-1;
				double y= img->GetYaxis()->GetBinCenter(iy+1);
				double y_min= img->GetYaxis()->GetBinCenter(iy_min+1);
				double y_max= img->GetYaxis()->GetBinCenter(iy_max+1);
				double dy_min= img->GetYaxis()->GetBinWidth(iy_min+1);
				double dy_max= img->GetYaxis()->GetBinWidth(iy_max+1);
				if(i==0){
					yBins_min.push_back(y_min-0.5*dy_min);
					yBins_max.push_back(y_max+0.5*dy_max);
				}

				sampledY(j)= y;
				sampledGridY[j]= y;

				//Init sampled values
				sampledBkgLevel(i,j)= 0;
				sampledRMS(i,j)= 0;	
				int nGoodPreviousSamples= 0;
				for(int s=1;s<=3;s++){
					if(j-s>=0 && sampledBkgLevel(i,j-s)!=0 && sampledRMS(i,j-s)!=0){
						sampledBkgLevel(i,j)+= sampledBkgLevel(i,j-s);
						sampledRMS(i,j)+= sampledRMS(i,j-s);
						nGoodPreviousSamples++;
					}
					if(i-s>=0 && sampledBkgLevel(i-s,j)!=0 && sampledRMS(i-s,j)!=0){
						sampledBkgLevel(i,j)+= sampledBkgLevel(i-s,j);
						sampledRMS(i,j)+= sampledRMS(i-s,j);
						nGoodPreviousSamples++;
					}
				}
				if(nGoodPreviousSamples>0) {
					sampledBkgLevel(i,j)/= (double)nGoodPreviousSamples;
					sampledRMS(i,j)/= (double)nGoodPreviousSamples;
				}

				//cout<<"BkgFinder::FindGridBkg(): INFO: Retrieving tile no. "<<counter<<" C("<<ix<<","<<iy<<" Cxy("<<x<<","<<y<<") ix("<<ix_min<<","<<ix_max<<") iy("<<iy_min<<","<<iy_max<<") x("<<x_min<<","<<x_max<<") y("<<y_min<<","<<y_max<<")"<<endl;	
				TileImg= img->GetTile(ix_min,ix_max,iy_min,iy_max);
			
				if(!TileImg) {
					cout<<"BkgFinder::FindGridBkg(): ERROR: Cannot get tile from image!"<<endl; 
					//return -1;
					throw std::runtime_error("Cannot get tile from image!");
				}
				TString tileName= Form("Tile%d",nTiles);
				TileImg->SetNameTitle(tileName,tileName);

				//## Compute and get stats for this tile
				if( TileImg->ComputeStats(true,false,false)<0 ){
					cerr<<"BkgFinder::FindGridBkg(): ERROR: Failed to compute stats for tile no. "<<nTiles<<"...skip tile!"<<endl;
					continue;		
				}
			
				//## Get access to tile stats
				stats= TileImg->GetPixelStats();
				if(!stats) {
					cerr<<"BkgFinder::FindGridBkg(): ERROR: Cannot get stats for tile no. "<<nTiles<<"...skip tile!"<<endl; 			
					continue;		
				}
				//TileImg->DumpStats();
	
				//## Check number of entries for this tile
				npix= stats->n;
				//cout<<"npix="<<npix<<endl;
			
				if(npix<10){
					//cout<<"BkgFinder::FindGridBkg(): WARN: Too few pixels ("<<npix<<") for tile no. "<<nTiles<<" to compute bkg...skip tile!"<<endl; 			
					continue;	
				}

				//## Compute bkg for this tile
				TileImgBkg= BkgFinder::ComputeBkg(TileImg,method);
				if(!TileImgBkg){
					cerr<<"BkgFinder::FindGridBkg(): ERROR: Background calculation failed for tile no. "<<nTiles<<"...skip tile!"<<endl; 			
					continue;
				}	
				sampledBkgLevel(i,j)= TileImgBkg->bkgLevel;
				sampledRMS(i,j)= TileImgBkg->bkgRMS;
				
				/*
				cout<<"== BKG DATA TILE NO. "<<counter<<" =="<<endl;
				cout<<"N="<<npix<<" xrange("<<TileImgBkg->ix_min<<","<<TileImgBkg->ix_max<<") xrange("<<TileImgBkg->iy_min<<","<<TileImgBkg->iy_max<<")"<<endl;
				cout<<"bkgLevel="<<TileImgBkg->bkgLevel<<" bkgRMS="<<TileImgBkg->bkgRMS<<endl;
				cout<<"sampledBkgLevel("<<i<<","<<j<<")="<<sampledBkgLevel(i,j)<<" sampledRMS("<<i<<","<<j<<")="<<sampledRMS(i,j)<<endl;
				cout<<"===================================="<<endl;
				*/
				//fBkgData.push_back(TileImgBkg);
				(bkgMapData->bkgData).push_back(TileImgBkg);

			}//end loop Y
		}//end loop X

		//## Fill sample bkg & rms maps	
		double xBinsLowEdge[xBins_min.size()+1];
		cout<<"BkgFinder::FindGridBkg(): INFO: xBinsLowEdge(";
		for(unsigned int i=0;i<xBins_min.size();i++){
			xBinsLowEdge[i]= xBins_min[i];
			cout<<xBinsLowEdge[i]<<", ";
		}
		xBinsLowEdge[xBins_min.size()]= xBins_max[xBins_max.size()-1];
		cout<<xBinsLowEdge[xBins_min.size()]<<")"<<endl;

		double yBinsLowEdge[yBins_min.size()+1];
		cout<<"BkgFinder::FindGridBkg(): INFO: yBinsLowEdge(";
	
		for(unsigned int i=0;i<yBins_min.size();i++){
			yBinsLowEdge[i]= yBins_min[i];
			cout<<yBinsLowEdge[i]<<", ";
		}
		yBinsLowEdge[yBins_min.size()]= yBins_max[yBins_max.size()-1];
		cout<<yBinsLowEdge[yBins_min.size()]<<")"<<endl;

		//fBackgroundMap= new TH2D("BackgroundMap","BackgroundMap",xBins_min.size(),xBinsLowEdge,yBins_min.size(),yBinsLowEdge);
		//fRMSMap= new TH2D("RMSMap","RMSMap",xBins_min.size(),xBinsLowEdge,yBins_min.size(),yBinsLowEdge);

		//int nBkgSamples= (int)fBkgData.size();
		int nBkgSamples= (int)(bkgMapData->bkgData).size();
		for(unsigned int i=0;i<nBkgSamples;i++){
			double bkgLevel= (bkgMapData->bkgData)[i]->bkgLevel;
			double bkgRMS= (bkgMapData->bkgData)[i]->bkgRMS;
			double sampleX= (bkgMapData->bkgData)[i]->ix_min;
			double sampleY= (bkgMapData->bkgData)[i]->iy_min;				
		}//end loop sampled tiles


		//## 2D interpolation
		// construct the grid in each dimension. 
  	// note that we will pass in a sequence of iterators pointing to the beginning of each grid
		cout<<"BkgFinder::FindGridBkg(): INFO: Build 2D grid for interpolation (nTilesX="<<nTilesX<<", nTilesY="<<nTilesY<<")..."<<endl;
  	std::vector< std::vector<double>::iterator > grid_iter_list;
  	grid_iter_list.push_back(sampledGridX.begin());
  	grid_iter_list.push_back(sampledGridY.begin());
  
  	// the size of the grid in each dimension
  	array<int,2> grid_sizes;
  	grid_sizes[0] = nTilesX;
  	grid_sizes[1] = nTilesY;
  
  	// total number of elements
  	long int num_elements = grid_sizes[0] * grid_sizes[1];

		// fill in the values of f(x,y) at the gridpoints. 
  	// we will pass in a contiguous sequence, values are assumed to be laid out C-style
		cout<<"BkgFinder::FindGridBkg(): INFO: Fill f(x,y) sampled values for interpolation..."<<endl;
  	
 	 	std::vector<double> fbkg_values(num_elements,0);
		std::vector<double> frms_values(num_elements,0);
		
  	for (int i=0; i<grid_sizes[0]; i++) {
    	for (int j=0; j<grid_sizes[1]; j++) {
				double x= sampledGridX[i];
				double y= sampledGridY[j];
				long long int gBin= i*grid_sizes[1] + j;
				//int gBin= biny*(fXaxis.GetNbins()+2) + binx;
	  		fbkg_values[gBin] = sampledBkgLevel(i,j);
				frms_values[gBin] = sampledRMS(i,j);
				//cout<<"gBin="<<gBin<<" (i,j)=("<<i<<","<<j<<"), (x,y)=("<<x<<","<<y<<") bkgInterp="<<sampledBkgLevel(i,j)<<" thisBkgRMSValue="<<sampledRMS(i,j)<<endl;
			}
  	}
  
  	// construct the interpolator. the last two arguments are pointers to the underlying data
		cout<<"BkgFinder::FindGridBkg(): INFO: Build the bkg interpolator..."<<endl;
  	InterpMultilinear<2, double> interpBkg_ML(grid_iter_list.begin(), grid_sizes.begin(), fbkg_values.data(), fbkg_values.data() + num_elements);
		cout<<"BkgFinder::FindGridBkg(): INFO: Build the RMS interpolator..."<<endl;
  	InterpMultilinear<2, double> interpRMS_ML(grid_iter_list.begin(), grid_sizes.begin(), frms_values.data(), frms_values.data() + num_elements);
  	
		//Construct interpolated grid
		// interpolate multiple values: create sequences for each coordinate
		cout<<"BkgFinder::FindGridBkg(): INFO: Build the interpolated grid Nx="<<Nx<<" Ny="<<Ny<<" x("<<xlim(0)<<","<<xlim(1)<<") y("<<ylim(0)<<","<<ylim(1)<<")..."<<endl;
  	
  	std::vector<double> interp_gridX = linspace(xlim(0),xlim(1), Nx);
		std::vector<double> interp_gridY = linspace(ylim(0),ylim(1), Ny);
  	int num_interp_elements = interp_gridX.size() * interp_gridY.size();
  	std::vector<double> interp_x(num_interp_elements);
  	std::vector<double> interp_y(num_interp_elements);
  	for (int i=0; i<interp_gridX.size(); i++) {
    	for (int j=0; j<interp_gridY.size(); j++) {
				//int gBin= i*interp_gridX.size() + j;
				long long int gBin= i*interp_gridY.size() + j;
	  		interp_x[gBin] = interp_gridX[i];
	  		interp_y[gBin] = interp_gridY[j];
			}
  	}
 	 	std::vector<double> interpolatedBkg(num_interp_elements);
  	std::vector<double> interpolatedRMS(num_interp_elements);
  	
		// pass in a sequence of iterators, one for each coordinate
  	std::vector< std::vector<double>::iterator > interp_list;
  	interp_list.push_back(interp_x.begin());
  	interp_list.push_back(interp_y.begin());
  
		//Interpolate sequence
		cout<<"BkgFinder::FindGridBkg(): INFO: Run the interpolation on grid..."<<endl;
		interpBkg_ML.interp_vec(num_interp_elements, interp_list.begin(), interp_list.end(), interpolatedBkg.begin());
  	interpRMS_ML.interp_vec(num_interp_elements, interp_list.begin(), interp_list.end(), interpolatedRMS.begin());
  	
		//## Fill images
		bkgMapData->BkgMap= (Img*)img->Clone("interpBkgMap");
		(bkgMapData->BkgMap)->Reset();

		bkgMapData->RMSMap= (Img*)img->Clone("interpRMSMap");
		(bkgMapData->RMSMap)->Reset();
		
		for (int i=0; i<interp_gridX.size(); i++) {
    	for (int j=0; j<interp_gridY.size(); j++) {
				//int gBin= i*interp_gridX.size() + j;
				int gBin= i*interp_gridY.size() + j;
				double x= interp_x[gBin];
				double y= interp_y[gBin];
				double thisBkgValue= interpolatedBkg[gBin];
	  		double thisBkgRMSValue= interpolatedRMS[gBin];
				if(thisBkgValue==0 || thisBkgRMSValue==0 ){
					//cout<<"BkgFinder::FindGridBkg(): WARN: INterpolated value is zero (bkg="<<thisBkgValue<<", rms="<<thisBkgRMSValue<<")"<<endl;
	
					/*
					int nGoodPreviousSamples= 0;
					for(int s=1;s<=3;s++){
						if(j-s>=0 && sampledBkgLevel(i,j-s)!=0 && sampledRMS(i,j-s)!=0){
							sampledBkgLevel(i,j)+= sampledBkgLevel(i,j-s);
							sampledRMS(i,j)+= sampledRMS(i,j-s);
							nGoodPreviousSamples++;
						}
						if(i-s>=0 && sampledBkgLevel(i-s,j)!=0 && sampledRMS(i-s,j)!=0){
							sampledBkgLevel(i,j)+= sampledBkgLevel(i-s,j);
							sampledRMS(i,j)+= sampledRMS(i-s,j);
							nGoodPreviousSamples++;
						}
					}
					if(nGoodPreviousSamples>0) {
						sampledBkgLevel(i,j)/= (double)nGoodPreviousSamples;
						sampledRMS(i,j)/= (double)nGoodPreviousSamples;
					}
					*/
				}

				(bkgMapData->BkgMap)->SetBinContent(i+1,j+1,thisBkgValue);
				(bkgMapData->RMSMap)->SetBinContent(i+1,j+1,thisBkgRMSValue);
			}//end loop bins Y
  	}//end loop bins X
		
	}//close try block
	catch( std::exception &ex ) {
		cerr << "BkgFinder::FindGridBkg(): ERROR: Exception catched: " << ex.what() << endl;
		//return -1;
		return 0;
  } 
	catch(...) { 
		cerr << "BkgFinder::FindGridBkg(): ERROR: C++ exception (unknown reason)" << endl;
		//return -1;
		return 0;
  }	
	
	return bkgMapData;
	//return 0;

}//end FindGridBkg()


// return an evenly spaced 1-d grid of doubles.
std::vector<double> BkgFinder::linspace(double first, double last, int len) {
  std::vector<double> result(len);
  double step = (last-first) / (len - 1);
  for (int i=0; i<len; i++) { result[i] = first + i*step; }
  return result;
}


//int BkgFinder::FindBkg(Img* img, Img::BkgMethod method, int boxSizeX, int boxSizeY, double boxSlideOffsetX, double boxSlideOffsetY){
BkgFinder::BkgMapData* BkgFinder::FindBkg(Img* img, Img::BkgMethod method, int boxSizeX, int boxSizeY, double boxSlideOffsetX, double boxSlideOffsetY){

	//if(!img) return -1;
	if(!img) return 0;
	
	//## Check args
	int Nx= img->GetNbinsX();
	int Ny= img->GetNbinsY();
	if(boxSizeX<=0 || boxSizeX>=Nx || boxSizeY<=0 || boxSizeY>=Ny) {
		cerr<<"BkgFinder::FindBkg(): ERROR: Invalid box size given (too small or larger than image size)!"<<endl;
		//return -1;
		return 0;
	}
	if(boxSlideOffsetX<=0 || boxSlideOffsetY<=0 ){
		cerr<<"BkgFinder::FindBkg(): ERROR: Invalid sliding box offset given (null or negative)!"<<endl;
		//return -1;
		return 0;
	}
	
	int TileSizeX= boxSizeX;
	int TileSizeY= boxSizeY;
	int TileOffsetX= boxSlideOffsetX*TileSizeX;
  int TileOffsetY= boxSlideOffsetY*TileSizeY;
	
	//## Check first tile offset compared to image size
	if( TileSizeX+TileOffsetX>=Nx || TileSizeY+TileOffsetY>=Ny ){
		cerr<<"BkgFinder::FindBkg(): WARNING: Tile offset or tile size are too large (only one tile can be accomodated in the image, are you sure?!)"<<endl;
	}
	cout<<"BkgFinder::FindBkg(): INFO: TileSize=("<<TileSizeX<<","<<TileSizeY<<") TileOffset("<<TileOffsetX<<","<<TileOffsetY<<")"<<endl;

	BkgFinder::BkgMapData* bkgMapData= 0;

	//Compute local bkg
	try {
	
		//## Initialize data
		Init();
		bkgMapData= new BkgFinder::BkgMapData;

		//## Loop over all number of tiles and compute bkg info for each one
		int nGoodTiles= 0;	
		int nGoodTilesX= 0;
		int nGoodTilesY= 0;

		nTiles= 0;
		nTilesX= 0;
		nTilesY= 0;
		int indexX= 0;
		int indexY= 0;
		Img* TileImg= 0;
		Img::StatsData* stats= 0;
		Img::BkgData* TileImgBkg= 0;
	
		std::vector<double> xBins_min;
		std::vector<double> xBins_max;
		std::vector<double> yBins_min;
		std::vector<double> yBins_max;
		
	
		while(indexY<=(Ny-1)){
			int iy_min= indexY;
			int iy_max= std::min(iy_min+TileSizeY-1,Ny-1);
			double y_min= img->GetYaxis()->GetBinLowEdge(iy_min+1);
			double y_max= img->GetYaxis()->GetBinLowEdge(iy_max+1)+img->GetYaxis()->GetBinWidth(iy_max+1);
			yBins_min.push_back(y_min);
			yBins_max.push_back(y_max);
			nTilesY++;
		
			indexX= 0;
			nTilesX= 0;
			xBins_min.clear();
			xBins_max.clear();

			while(indexX<=(Nx-1)){
				nTiles++;
				nTilesX++;
				int ix_min= indexX;
				int ix_max= std::min(ix_min+TileSizeX-1,Nx-1);
				double x_min= img->GetXaxis()->GetBinLowEdge(ix_min+1);
				double x_max= img->GetXaxis()->GetBinLowEdge(ix_max+1)+img->GetXaxis()->GetBinWidth(ix_max+1);
				xBins_min.push_back(x_min);
				xBins_max.push_back(x_max);
			
				indexX+= TileOffsetX;
			
				cout<<"BkgFinder::FindBkg(): INFO: Retrieving tile no. "<<nTiles<<" ix("<<ix_min<<","<<ix_max<<") iy("<<iy_min<<","<<iy_max<<") x("<<x_min<<","<<x_max<<") y("<<y_min<<","<<y_max<<")"<<endl;	
				TileImg= img->GetTile(ix_min,ix_max,iy_min,iy_max);
			
				if(!TileImg) {
					cout<<"BkgFinder::FindBkg(): ERROR: Cannot get tile from image!"<<endl; 
					throw std::runtime_error("Cannot get tile from image!");
				}
				TString tileName= Form("Tile%d",nTiles);
				TileImg->SetNameTitle(tileName,tileName);

				//## Compute and get stats for this tile
				if( TileImg->ComputeStats(true,false,false)<0 ){
					cerr<<"BkgFinder::FindBkg(): ERROR: Failed to compute stats for tile no. "<<nTiles<<"...skip tile!"<<endl;
					continue;		
				}
			
				//## Get access to tile stats
				stats= TileImg->GetPixelStats();
				if(!stats) {
					cerr<<"BkgFinder::FindBkg(): ERROR: Cannot get stats for tile no. "<<nTiles<<"...skip tile!"<<endl; 			
					continue;		
				}
				//TileImg->DumpStats();
	
				//## Check number of entries for this tile
				npix= stats->n;
				if(npix<10){
					cout<<"BkgFinder::FindBkg(): WARN: Too few pixels ("<<npix<<") for tile no. "<<nTiles<<" to compute bkg...skip tile!"<<endl; 			
					continue;	
				}

				//## Compute bkg for this tile
				TileImgBkg= BkgFinder::ComputeBkg(TileImg,method);
				if(!TileImgBkg){
					cerr<<"BkgFinder::FindBkg(): ERROR: Background calculation failed for tile no. "<<nTiles<<"...skip tile!"<<endl; 			
					continue;
				}	
				cout<<"== BKG DATA TILE NO. "<<nTiles<<" =="<<endl;
				cout<<"N="<<npix<<" xrange("<<TileImgBkg->ix_min<<","<<TileImgBkg->ix_max<<") xrange("<<TileImgBkg->iy_min<<","<<TileImgBkg->iy_max<<")"<<endl;
				cout<<"bkgLevel="<<TileImgBkg->bkgLevel<<" bkgRMS="<<TileImgBkg->bkgRMS<<endl;
				cout<<"===================================="<<endl;
				//fBkgData.push_back(TileImgBkg);
				(bkgMapData->bkgData).push_back(TileImgBkg);
				
				/*
				//## Dump bkg info to file
				tileId= nTiles;
				tileMinX= ix_min;
				tileMaxX= ix_max;
				tileMinY= iy_min;
				tileMaxY= iy_max;
				mean= stats->mean;
				rms= stats->rms;
				median= stats->median;
				medianRMS= stats->medianRMS;
				bkgLevel= TileImgBkg->bkgLevel;
				bkgRMS= TileImgBkg->bkgRMS;
				fDataTree->Fill();
				*/
			}//end loop X
			indexY+= TileOffsetY;
		}//end loop Y

		cout<<"BkgFinder::FindBkg(): INFO: nTiles="<<nTiles<<"("<<nTilesX<<" x "<<nTilesY<<")"<<endl;
		//fInfoTree->Fill();

	
		//## Interpolation 
		double xBinsLowEdge[xBins_min.size()+1];
		cout<<"BkgFinder::FindBkg(): INFO: xBinsLowEdge(";
		for(unsigned int i=0;i<xBins_min.size();i++){
			xBinsLowEdge[i]= xBins_min[i];
			cout<<xBinsLowEdge[i]<<", ";
		}
		xBinsLowEdge[xBins_min.size()]= xBins_max[xBins_max.size()-1];
		cout<<xBinsLowEdge[xBins_min.size()]<<")"<<endl;

		double yBinsLowEdge[yBins_min.size()+1];
		cout<<"BkgFinder::FindBkg(): INFO: yBinsLowEdge(";
	
		for(unsigned int i=0;i<yBins_min.size();i++){
			yBinsLowEdge[i]= yBins_min[i];
			cout<<yBinsLowEdge[i]<<", ";
		}
		yBinsLowEdge[yBins_min.size()]= yBins_max[yBins_max.size()-1];
		cout<<yBinsLowEdge[yBins_min.size()]<<")"<<endl;

		
		TH2D* BackgroundMap= new TH2D("BackgroundMap","BackgroundMap",xBins_min.size(),xBinsLowEdge,yBins_min.size(),yBinsLowEdge);
		TH2D* RMSMap= new TH2D("RMSMap","RMSMap",xBins_min.size(),xBinsLowEdge,yBins_min.size(),yBinsLowEdge);

		cout<<"BkgFinder::FindBkg(): INFO: BackgroundMap rangex("<<BackgroundMap->GetXaxis()->GetXmin()<<","<<BackgroundMap->GetXaxis()->GetXmax()<<") yrange("<<BackgroundMap->GetYaxis()->GetXmin()<<","<<BackgroundMap->GetYaxis()->GetXmax()<<")"<<endl;
	
		//int nBkgSamples= (int)fBkgData.size();
		int nBkgSamples= (int)(bkgMapData->bkgData).size();
	
		for(unsigned int i=0;i<nBkgSamples;i++){
			double bkgLevel= (bkgMapData->bkgData)[i]->bkgLevel;
			double bkgRMS= (bkgMapData->bkgData)[i]->bkgRMS;
			double sampleX= (bkgMapData->bkgData)[i]->ix_min;
			double sampleY= (bkgMapData->bkgData)[i]->iy_min;
			//cout<<"BkgFinder::FindBkg(): INFO: ix("<<fBkgData[i]->ix_min<<","<<fBkgData[i]->ix_max<<") iy("<<fBkgData[i]->iy_min<<","<<fBkgData[i]->iy_max<<") sampleXY("<<sampleX<<","<<sampleY<<") "<<endl;
		
			BackgroundMap->Fill(sampleX,sampleY,bkgLevel);
			RMSMap->Fill(sampleX,sampleY,bkgRMS);
		}//end loop sampled tiles


		//## Load R library for interpolation
		cout<<"BkgFinder::FindBkg(): INFO: Loading needed R packages..."<<endl;
		std::string RCmd= std::string("library(\"akima\");");
		fR->parseEvalQ(RCmd);

		int nSamplesX= BackgroundMap->GetNbinsX();
		int nSamplesY= BackgroundMap->GetNbinsY();
		int nTotSamples= nSamplesX*nSamplesY;
		Rcpp::NumericVector sampledX(nSamplesX);
		Rcpp::NumericVector sampledY(nSamplesY);
		Rcpp::NumericMatrix sampledBkgLevel(nSamplesX,nSamplesY);
		Rcpp::NumericMatrix sampledBkgRMS(nSamplesX,nSamplesY);
		double dx= 1;//output grid spacing in x direction.
		double dy= 1;//output grid spacing in y direction.
		Rcpp::NumericVector xlim(2);//vector of length 2 giving lower and upper limit for range x coordinates used for output grid
		Rcpp::NumericVector ylim(2);//vector of length 2 giving lower and upper limit for range of y coordinates used for output grid
		xlim(0)= img->GetXaxis()->GetXmin();
		xlim(1)= img->GetXaxis()->GetXmax();
		ylim(0)= img->GetYaxis()->GetXmin();
		ylim(1)= img->GetYaxis()->GetXmax();

		int index= 0;
		for(int i=0;i<nSamplesX;i++){
			double x= BackgroundMap->GetXaxis()->GetBinCenter(i+1);
			//sampledX(index)= x;
			sampledX(i)= x;
			for(int j=0;j<nSamplesY;j++){
				double y= BackgroundMap->GetYaxis()->GetBinCenter(j+1);				
				//sampledY(index)= y;	
				sampledY(j)= y;
				sampledBkgLevel(i,j)= BackgroundMap->GetBinContent(i+1,j+1);
				sampledBkgRMS(i,j)= RMSMap->GetBinContent(i+1,j+1);
				index++;
			}//end loop Y bins
		}//end loop X bins

		Rcpp::NumericVector interpX(Nx*Ny);
		Rcpp::NumericVector interpY(Nx*Ny);
		index= 0;
		for(int i=0;i<Nx;i++){	
			double x= img->GetXaxis()->GetBinCenter(i+1);
			for(int j=0;j<Ny;j++){
				double y= img->GetYaxis()->GetBinCenter(j+1);
				interpX(index)= x;
				interpY(index)= y;
				index++;
			}
		}

		//## Import data matrix in R environment
		(*fR)["sampledX"]= sampledX;
		(*fR)["sampledY"]= sampledY;
		(*fR)["sampledBkgLevel"]= sampledBkgLevel;
		(*fR)["sampledBkgRMS"]= sampledBkgRMS;
		(*fR)["interpX"]= interpX;
		(*fR)["interpY"]= interpY;
		(*fR)["xlim"]= xlim;
		(*fR)["ylim"]= ylim;
		(*fR)["dx"]= dx;
		(*fR)["dy"]= dy;
		

		//### Run interpolation algorithm
		cout<<"BkgFinder::FindBkg(): INFO: Run interpolation algorithm..."<<endl;
		RCmd= std::string("bkgLevelInterp <- bicubic.grid(sampledX,sampledY,sampledBkgLevel,xlim,ylim,dx,dy);");
		//RCmd= std::string("bkgLevelInterp <- bicubic(sampledX,sampledY,sampledBkgLevel,interpX,interpY);");
		fR->parseEval(RCmd);
		
		RCmd= std::string("bkgRMSInterp <- bicubic.grid(sampledX,sampledY,sampledBkgRMS,xlim,ylim,dx,dy);");
		//RCmd= std::string("bkgRMSInterp <- bicubic(sampledX,sampledY,sampledBkgRMS,interpX,interpY);");
		fR->parseEval(RCmd);
		
		//## Retrieve results
		cout<<"BkgFinder::FindBkg(): INFO: Retrieving results..."<<endl;
		Rcpp::NumericVector interpolatedGridX = fR->parseEval(std::string("bkgLevelInterp$x"));
		Rcpp::NumericVector interpolatedGridY = fR->parseEval(std::string("bkgLevelInterp$y"));
		Rcpp::NumericMatrix bkgLevelInterpolatedValues = fR->parseEval(std::string("bkgLevelInterp$z"));
		Rcpp::NumericMatrix bkgRMSInterpolatedValues = fR->parseEval(std::string("bkgRMSInterp$z"));
	
		//## Fill images
		cout<<"BkgFinder::FindBkg(): INFO: Filling background images..."<<endl;
		(bkgMapData->BkgMap)= (Img*)img->Clone("interpBkgMap");
		(bkgMapData->BkgMap)->Reset();

		(bkgMapData->RMSMap)= (Img*)img->Clone("interpRMSMap");
		(bkgMapData->RMSMap)->Reset();

		/*
		if(!fInterpolatedBackgroundMap){
			fInterpolatedBackgroundMap= (Img*)img->Clone("interpBkgMap");
		}
		fInterpolatedBackgroundMap->Reset();
		if(!fInterpolatedRMSMap){
			fInterpolatedRMSMap= (Img*)img->Clone("interpRMSMap");
		}
		fInterpolatedRMSMap->Reset();
		*/

		for(int i=0;i<bkgLevelInterpolatedValues.nrow();i++){
			for(int j=0;j<bkgLevelInterpolatedValues.ncol();j++){
				double thisBkgValue= bkgLevelInterpolatedValues(i,j);
				double thisBkgRMSValue= bkgRMSInterpolatedValues(i,j);
				//fInterpolatedBackgroundMap->SetBinContent(i+1,j+1,thisBkgValue);
				//fInterpolatedRMSMap->SetBinContent(i+1,j+1,thisBkgRMSValue);
				(bkgMapData->BkgMap)->SetBinContent(i+1,j+1,thisBkgValue);
				(bkgMapData->RMSMap)->SetBinContent(i+1,j+1,thisBkgRMSValue);
			}//end loop cols
		}//end loop nrows

	}//close try block
	catch( std::exception &ex ) {
		cerr << "BkgFinder::FindBkg(): ERROR: Exception catched: " << ex.what() << endl;
		//return -1;
		return 0;
  } 
	catch(...) { 
		cerr << "BkgFinder::FindBkg(): ERROR: C++ exception (unknown reason)" << endl;
		//return -1;
		return 0;
  }	
	
	//return 0;
	return bkgMapData;
	
}//close BkgFinder::FindBkg()



BkgFinder::BkgMapData* BkgFinder::FindSuperpixelBkg(Img* img,Img::BkgMethod method, int boxSizeX, int boxSizeY, double gridStepSizeX, double gridStepSizeY){
	
	//if(!img) return -1;
	if(!img) return 0;

	//## Check args
	int Nx= img->GetNbinsX();
	int Ny= img->GetNbinsY();
	if(fRegionSize<=0 || fRegionSize>=Nx || fRegionSize>=Ny) {
		cerr<<"BkgFinder::FindSuperpixelBkg(): ERROR: Invalid region size given (too small or larger than image size)!"<<endl;
		//return -1;
		return 0;
	}
	
	BkgFinder::BkgMapData* bkgMapData= 0;

	try{
		//## Initialize data	
		Init();

		RInside* fR= RInside::instancePtr();
		if(!fR){
			cerr<<"BkgFinder::FindSuperpixelBkg(): ERROR: Cannot retrieve RInside instance (did you create one in your main application?)!"<<endl;
			fR= new RInside;
		}	

		//## Segment the image
		VLSlicSegmentation slic;
		slic.SetLogContrastMapping(false);
  	int status= slic.RunSegmentation(img,fRegionSize,fRegularization,fMinRegionArea,false);
		if(status<0){
			cerr<<"BkgFinder::FindSuperpixelBkg(): ERROR: Segmentation failed!"<<endl;
			throw std::runtime_error("Image segmentation failed!");
		}
	
		
		//## Tag the regions as background or significative
		bool includeSpatialPar= false;
		bool includeCurvPar= false;
		double CL= 0.975;
		Region* bkgRegion= slic.FindBackgroundRegion(CL,includeSpatialPar,includeCurvPar);
		if(!bkgRegion){
			cerr<<"BkgFinder::FindSuperpixelBkg(): ERROR: Failed to select background region!"<<endl;
			throw std::runtime_error("Selection of background regions failed!");
		}
		std::vector<Region*> regions= slic.GetRegions();//get list of regions AFTER previous step
		cout<<"BkgFinder::FindSuperpixelBkg(): INFO: "<<regions.size()<<" regions found..."<<endl;	

		Img* bkgImage= slic.GetBackgroundImage();
		if(!bkgImage){
			cerr<<"BkgFinder::FindSuperpixelBkg(): ERROR: Failed to get background image!"<<endl;
			throw std::runtime_error("Cannot get bkgImage!");			
		}

		//## Find grid background estimate on the bkgImage
		//status= FindGridBkg(bkgImage,method,boxSizeX,boxSizeY,gridStepSizeX,gridStepSizeY);
		bkgMapData= FindGridBkg(bkgImage,method,boxSizeX,boxSizeY,gridStepSizeX,gridStepSizeY);
		
		//if(status<0){
		if(!bkgMapData){
			cerr<<"BkgFinder::FindSuperpixelBkg(): ERROR: Failed to compute grid background!"<<endl;
			throw std::runtime_error("Failed to compute grid background!");
		}

	}//close try block
	catch( std::exception &ex ) {
		cerr << "BkgFinder::FindSuperpixelBkg(): ERROR: Exception catched: " << ex.what() << endl;
		//return -1;
		return 0;
  } 
	catch(...) { 
		cerr << "BkgFinder::FindSuperpixelBkg(): ERROR: C++ exception (unknown reason)" << endl;
		//return -1;
		return 0;
  }	

	
	//return 0;
	return bkgMapData;

}//close BkgFinder::FindSuperpixelBkg()


Img::BkgData* BkgFinder::ComputeBkg(Img* img,Img::BkgMethod method){

	if(!img) return 0;


	//## Compute bkg
	Img::BkgData* bkg= 0;

	// == mean bkg ==
	if(method == Img::eMeanBkg){
		bkg= GetMeanBkg(img);
	}
		
	//== median bkg ==
	else if(method == Img::eMedianBkg){
		bkg= GetMedianBkg(img);
	}
		
	//== robust bkg ==
	else if(method == Img::eRobustBkg){
		bkg= GetRobustBkg(img,fRegionSize,fRegularization,fMinRegionArea, fColorEps);
	}
	//== simple robust bkg ==
	else if(method == Img::eSimpleRobustBkg){
		bkg= GetSimpleRobustBkg(img,fRegionSize,fRegularization,fMinRegionArea, fColorEps);
	}
	//== dummy robust bkg ==
	else if(method == Img::eDummyRobustBkg){
		bkg= GetDummyRobustBkg(img);
	}

	//== BIWEIGHT bkg ==
	else if(method == Img::eBiWeightBkg){
		bkg= GetBiWeightBkg(img);
	}
	//== CLIPPED MEDIAN bkg ==
	else if(method == Img::eMedianClippedBkg){
		bkg= GetMedianClippedBkg(img);
	}

	//invalid method
	else{
		cerr<<"BkgFinder::ComputeBkg(): ERROR: Invalid bkg method selected ("<<method<<")...exit!"<<endl;
		return 0;
	}

	//Set x, y range
	//int ix_min= img->GetXaxis()->GetXmin();	
	//int iy_min= img->GetYaxis()->GetXmin();	
	//int ix_max= img->GetXaxis()->GetXmax();	
	//int iy_max= img->GetYaxis()->GetXmax();	
	int ix_min= img->GetXaxis()->GetBinCenter(1);	
	int iy_min= img->GetYaxis()->GetBinCenter(1);	
	int ix_max= img->GetXaxis()->GetBinCenter(img->GetNbinsX());	
	int iy_max= img->GetYaxis()->GetBinCenter(img->GetNbinsY());	
	bkg->ix_min= ix_min;
	bkg->ix_max= ix_max;
	bkg->iy_min= iy_min;
	bkg->iy_max= iy_max;
	

	return bkg;

}//close BkgFinder::ComputeBkg()




Img::BkgData* BkgFinder::GetMeanBkg(Img* img){

	if(!img) return 0;

	//Get stats
	Img::StatsData* stats= img->GetPixelStats();
	if(!stats) return 0;

	Img::BkgData* bkg= new Img::BkgData;
	bkg->npix= stats->n;
	bkg->bkgLevel= stats->mean;
	bkg->bkgRMS= stats->rms;
	bkg->isReliableBkg= true;
	if(stats->n<3) bkg->isReliableBkg= false;

	return bkg;

}//close BkgFinder::GetMeanBkg()


Img::BkgData* BkgFinder::GetMedianBkg(Img* img){

	if(!img) return 0;

	//Get stats
	Img::StatsData* stats= img->GetPixelStats();
	if(!stats) return 0;

	Img::BkgData* bkg= new Img::BkgData;
	bkg->npix= stats->n;
	bkg->bkgLevel= stats->median;
	bkg->bkgRMS= stats->medianRMS;
	bkg->isReliableBkg= true;
	if(stats->n<3) bkg->isReliableBkg= false;

	return bkg;

}//close BkgFinder::GetMedianBkg()


Img::BkgData* BkgFinder::GetBiWeightBkg(Img* img){

	if(!img) return 0;

	//Get stats
	Img::StatsData* stats= img->GetPixelStats();
	if(!stats) return 0;

	Img::BkgData* bkg= new Img::BkgData;
	bkg->npix= stats->n;
	bkg->bkgLevel= stats->bwLocation;
	bkg->bkgRMS= stats->bwScale;
	bkg->isReliableBkg= true;
	if(stats->n<3) bkg->isReliableBkg= false;

	return bkg;

}//close GetBiWeightBkg()

Img::BkgData* BkgFinder::GetMedianClippedBkg(Img* img){

	if(!img) return 0;

	//Get stats
	Img::StatsData* stats= img->GetPixelStats();
	if(!stats) return 0;

	double bkgLevel= stats->clippedMedian;
	double bkgRMS= stats->clippedRMS;

	Img::BkgData* bkg= new Img::BkgData;
	bkg->npix= stats->n;
	bkg->bkgLevel= bkgLevel;
	bkg->bkgRMS= bkgRMS;
	bkg->isReliableBkg= true;
	if(stats->n<3) bkg->isReliableBkg= false;

	return bkg;

}//close BkgFinder::GetMedianClippedBkg()


Img::BkgData* BkgFinder::GetSimpleRobustBkg(Img* img, int regionSize,double regularization,int minRegionArea,double colorEps){

	if(!img) return 0;

	Img::BkgData* bkg= 0;

	//## Segment the image
	VLSlicSegmentation slic;
  int status= slic.RunSegmentation(img,regionSize,regularization,minRegionArea,true,colorEps);
	if(status<0){
		cerr<<"BkgFinder::GetRobustBkg(): ERROR: Segmentation failed!"<<endl;
		return 0;
	}

	TString imgName= img->GetName();
	TString imgColoredName= Form("Segmented%s",std::string(imgName).c_str());

	//std::vector< std::vector<long int> > pixelLabels= slic.GetPixelClusterIds();
	std::vector<Region*> regions= slic.GetRegions();
	//TGraph* contours= slic.ComputeClusterContours(img,pixelLabels,regions);
	Img* img_colored= slic.GetClusterColoredImage(img,regions);	
	img_colored->SetNameTitle(imgColoredName,imgColoredName);

	//## Write image to file
	//fOutputFile->cd();
	//if(img) img->Write();
	//if(img_colored) img_colored->Write(); 

	//Get background region
	Region* bkgRegion= slic.FindSimpleBackgroundRegion();
	if(!bkgRegion){
		cerr<<"BkgFinder::GetRobustBkg(): ERROR: Failed to identity bkg region, trying with median background method!"<<endl;

		bkg= GetMedianBkg(img);
		if(!bkg) {
			cerr<<"BkgFinder::GetRobustBkg(): ERROR: Also median background method failed...no background will be returned!"<<endl;
			return 0;
		}
		else {
			return bkg;
		}
	}//close if
	
	double median= bkgRegion->fMedian;
	double medianRMS= bkgRegion->fMedianRMS;
	cout<<"BkgFinder::GetRobustBkg(): INFO: median="<<median<<", medianRMS="<<medianRMS<<endl;

	bkg= new Img::BkgData;
	bkg->bkgLevel= median;
	bkg->bkgRMS= medianRMS;
	bkg->isReliableBkg= true;

	if(bkgRegion) {
		delete bkgRegion;
		bkgRegion= 0;
	}
	
	return bkg;

}//close BkgFinder::GetSimpleRobustBkg(Img* img)


Img::BkgData* BkgFinder::GetRobustBkg(Img* img, int regionSize,double regularization,int minRegionArea,double colorEps){
	
	if(!img) return 0;

	Img::BkgData* bkg= 0;

	//## Segment the image
	VLSlicSegmentation slic;
  int status= slic.RunSegmentation(img,regionSize,regularization,minRegionArea,false,colorEps);
	if(status<0){
		cerr<<"BkgFinder::GetRobustBkg(): ERROR: Segmentation failed!"<<endl;
		return 0;
	}

	TString imgName= img->GetName();
	TString imgColoredName= Form("Segmented%s",std::string(imgName).c_str());

	//std::vector< std::vector<long int> > pixelLabels= slic.GetPixelClusterIds();
	std::vector<Region*> regions= slic.GetRegions();
	//TGraph* contours= slic.ComputeClusterContours(img,pixelLabels,regions);
	Img* img_colored= slic.GetClusterColoredImage(img,regions);	
	img_colored->SetNameTitle(imgColoredName,imgColoredName);

	//## Write image to file
	//fOutputFile->cd();
	//if(img) img->Write();
	//if(img_colored) img_colored->Write(); 

	//Get background region
	bool includeSpatialPar= false;
	bool includeCurvPar= false;
	double CL= 0.975;
	Region* bkgRegion= slic.FindBackgroundRegion(CL,includeSpatialPar,includeCurvPar);
	if(!bkgRegion){
		cerr<<"BkgFinder::GetRobustBkg(): ERROR: Failed to identity bkg region, trying with median background method!"<<endl;

		bkg= GetMedianBkg(img);
		if(!bkg) {
			cerr<<"BkgFinder::GetRobustBkg(): ERROR: Also median background method failed...no background will be returned!"<<endl;
			return 0;
		}
		else {
			return bkg;
		}
	}//close if
	
	double median= bkgRegion->fMedian;
	double medianRMS= bkgRegion->fMedianRMS;
	int nPix= bkgRegion->fNPix;
	cout<<"BkgFinder::GetRobustBkg(): INFO: median="<<median<<", medianRMS="<<medianRMS<<endl;

	bkg= new Img::BkgData;
	bkg->npix= nPix;
	bkg->bkgLevel= median;
	bkg->bkgRMS= medianRMS;
	bkg->isReliableBkg= true;
	
	if(bkgRegion) {
		delete bkgRegion;
		bkgRegion= 0;
	}
	
	return bkg;

}//close BkgFinder::GetRobustBkg()


Img::BkgData* BkgFinder::GetDummyRobustBkg(Img* img){

	if(!img) return 0;

	//Get stats
	TH1D* pixelHisto= img->GetPixelHisto();
	if(!pixelHisto) return 0;

	TH1D* pixelTransformedHisto= (TH1D*)pixelHisto->Clone("pixelTransformedHisto");
	pixelTransformedHisto->Reset();
	int nPix= 0;

	for(int i=0;i<img->GetNbinsX();i++){
		for(int j=0;j<img->GetNbinsY();j++){
			double w= img->GetBinContent(i+1,j+1);
			if(w==0) continue;
			if(w<0) {
				nPix+= 2;
				pixelTransformedHisto->Fill(-w);
				pixelTransformedHisto->Fill(w);
			}
		}//end loop y
	}//end loop x

	Img::BkgData* bkg= new Img::BkgData;
	bkg->npix= nPix;
	bkg->bkgLevel= pixelTransformedHisto->GetMean();
	bkg->bkgRMS= pixelTransformedHisto->GetRMS();
	bkg->isReliableBkg= true;
	if(pixelHisto->GetEntries()<3) bkg->isReliableBkg= false;

	pixelTransformedHisto->Delete();
	pixelTransformedHisto= 0;

	return bkg;

}//close GetSimpleRobustBkg()
