/**
* @file WTFilter.cc
* @class WTFilter
* @brief WTFilter
*
* Stationary Wavelet Transform filter class
* @author S. Riggi
* @date 20/01/2015
*/

#include <WTFilter.h>

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

//#ifndef __CINT__
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>

using namespace std;
//#endif

//typedef std::vector<Img> ImgCollection;

TMatrixD* WTFilter::fH;


WTFilter::WTFilter() {

	fH= 0;
	Init();

}//close costructor


WTFilter::~WTFilter(){

	if(fH) fH->Delete();

}//close destructor


void WTFilter::Init(){

	//## Init filter matrix
	//## Filter is a B3-spline filter 
	if(!fH) {
		int HValues[5][5]= { {1,4,6,4,1}, {4,16,24,16,4}, {6,24,36,24,6}, {4,16,24,16,4}, {1,4,6,4,1} };
		double HScaleFactor= 256.;
	
		fH= new TMatrixD(5,5);
		for(int i=0;i<5;i++){
			for(int j=0;j<5;j++){
				(*fH)(i,j)= HValues[i][j]/HScaleFactor;
			}
		}
		fH->Print();
	}//close if

}//close WTFilter::Init()


//std::vector<Img*> WTFilter::GetDecomposition(Img image,int nScales){
std::vector<Img*> WTFilter::GetDecomposition(Img* image,int nScales){
	
	//## Init filters
	Init();

	//## Get image matrix from image class
	cout<<"WTFilter::GetDecomposition(): Getting image matrix ...."<<endl;
	//TMatrixD I= image.GetMatrix();
	TMatrixD I= *(image->GetMatrix());
	int nrows= I.GetNrows();
	int ncols= I.GetNcols();

	//## Start computation
	cout<<"WTFilter::GetDecomposition(): Start computation ...."<<endl;
	
	TMatrixD* helper;
	std::vector<TMatrixD> F;
	std::vector<TMatrixD> W;//W(0) is the smoothed array

	//TMatrixD F[nScales+1];
	//TMatrixD W[nScales+1];//W(0) is the smoothed array
	for(int k=0;k<nScales+1;k++){
		helper= new TMatrixD(nrows,ncols);
		helper->Zero();
		F.push_back(*helper);
		W.push_back(*helper);
	}

	cout<<"WTFilter::GetDecomposition(): Init F with full reso image ..."<<endl;
	F[0]= I;//init F with full reso image

	for(int n=1;n<nScales+1;n++){
		cout<<"WTFilter::GetDecomposition(): Compute convolution for scale "<<n<<" ..."<<endl;
	
		F[n].Zero();
		F[n]= GetConvolution(*fH,F[n-1],n);
		
		W[n-1]= F[n-1]-F[n];
	}//end loop scales

	W[nScales]= F[nScales];

	//## Convert TMatrixD collection in img objects
	std::vector<Img*> imgCollection;
	imgCollection.clear();
	imgCollection.resize(0);
	Img* image_helper= 0;

	//for(int n=1;n<nScales+1;n++){
	for(int n=0;n<nScales+1;n++){
		TString imageName= Form("%s-WT%d",image->GetName(),n);
		TString imageTitle= Form("%s-WT%d",image->GetTitle(),n);
		image_helper= new Img(imageName,imageName,image->GetNbinsX(),image->GetXaxis()->GetXmin(),image->GetXaxis()->GetXmax(),image->GetNbinsY(),image->GetYaxis()->GetXmin(),image->GetYaxis()->GetXmax());

		//image_helper= (Img*)image->Clone(imageName);
		//image_helper->SetNameTitle(imageName,imageTitle);
		image_helper->Reset();

		for(int i=0;i<nrows;i++){
			for(int j=0;j<ncols;j++){
				int ix= j; 
				int iy= i;
				double w= (W[n])(i,j);
				image_helper->SetBinContent(ix+1,iy+1,w);
			}//end loop matrix cols
		}//end loop matrix rows
		imgCollection.push_back(image_helper);
	}//end loop scales

	return imgCollection;

}//close WTFilter::GetDecomposition()


TMatrixD WTFilter::GetConvolution(TMatrixD H,TMatrixD F,int n){
	
	//#### EDGE TREATMENT  #############
	//mirroring: c(k + N) = c(N - k)
	//periodicity: c(k + N) = c(k)
	//continuity: c(k + N) = c(N)
	//#####################################


	//Init convolution matrix
	int nrows= F.GetNrows();
	int ncols= F.GetNcols();
	int nfilter_rows= H.GetNrows();
	int nfilter_cols= H.GetNcols();
	TMatrixD C(nrows,ncols);
	C.Zero();

  
	//Extend matrix at the borders with reflection conditions	
	int rowIndex= (nfilter_rows-1)/2;
	int colIndex= (nfilter_cols-1)/2;
	int filterGap= pow(2,n-1); 

	int maxoffset_x= filterGap*rowIndex;
	int maxoffset_y= filterGap*colIndex;
	//TMatrixD FExt(nrows+2*maxoffset_x,ncols+2*maxoffset_y);
	//FExt.Zero();
	
	/*
	for(int i=-maxoffset_x;i<nrows+maxoffset_x;i++){
		double Fx= 0;
		int index_x= GetMirrorIndex(i,nrows);
		
		for(int j=-maxoffset_y;j<ncols+maxoffset_y;j++){
			int index_y= GetMirrorIndex(j,ncols);

			FExt(i+maxoffset_x,j+maxoffset_y)= F(index_x, index_y);
		}//end loop cols
	}//end loop rows
	*/

	//Compute convolution
	for(int i=0;i<nrows;i++){
		for(int j=0;j<ncols;j++){

			for(int k=-rowIndex;k<rowIndex;k++){//start loops on filter box rows
				int index_x= i + filterGap*k;
				int mirror_index_x= GetMirrorIndex(index_x,nrows);
		
				for(int l=-colIndex;l<colIndex;l++){//start loops on filter box cols
					int index_y= j + filterGap*l;
					//int mirror_index_y= GetMirrorIndex(index_y,nrows);
					int mirror_index_y= GetMirrorIndex(index_y,ncols);
		
					double h= H(k+rowIndex,l+colIndex);
					//double f= FExt(index_x+rowIndex,index_y+colIndex);
					double f= F(mirror_index_x,mirror_index_y); 
					C(i,j)+= h*f;
				}//end loop filter cols
			}//end loop filter rows
		
		}//end loop cols
	}//end loop rows
	
	cout<<"WTFilter::GetConvolution(): INFO: Convolution matrix at stage "<<n<<" ..."<<endl;
	//C.Print();

	return C;

}//close WTFilter::GetConvolution()


int WTFilter::GetMirrorIndex(int index,int N){
	
	int mirror_index= 0;
	if(index>=0 && index<=(N-1)) {
		mirror_index= index;
	}
	else if(index<0){
		//mirror_index= GetMirrorIndex(-index,N);
		mirror_index= -index;
	}
	else if(index>(N-1)){
		//mirror_index= GetMirrorIndex(N-1-index,N);
		//mirror_index= GetMirrorIndex(2*(N-1)-index,N);
		mirror_index= 2*(N-1)-index;
	}
	else{
		cerr<<"WTFilter::GetMirrorIndex(): ERROR: Invalid index of matrix size passed...exit!"<<endl;
		exit(1);
	}	
	
	return mirror_index;

}//close GetMirrorIndex()

		

