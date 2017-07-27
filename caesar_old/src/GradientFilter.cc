/**
* @file GradientFilter.cc
* @class GradientFilter
* @brief GradientFilter
*
* Gradient Filter
* @author S. Riggi
* @date 20/01/2015
*/


#include <GradientFilter.h>

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


#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>

using namespace std;

TMatrixD* GradientFilter::fGx_Sobel;
TMatrixD* GradientFilter::fGy_Sobel;
TMatrixD* GradientFilter::fGx_Scharr;
TMatrixD* GradientFilter::fGy_Scharr;
TMatrixD* GradientFilter::fG_Laplace;
TMatrixD* GradientFilter::fG_GausLaplace;
TMatrixD* GradientFilter::fG_NormLoG;

TMatrixD* GradientFilter::fG_KirschN;
TMatrixD* GradientFilter::fG_KirschS;
TMatrixD* GradientFilter::fG_KirschW;
TMatrixD* GradientFilter::fG_KirschE;
TMatrixD* GradientFilter::fG_KirschNE;
TMatrixD* GradientFilter::fG_KirschSE;
TMatrixD* GradientFilter::fG_KirschNW;
TMatrixD* GradientFilter::fG_KirschSW;

	
GradientFilter::GradientFilter() {

	fGx_Sobel= 0;
	fGy_Sobel= 0;
	fGx_Scharr= 0;
	fGy_Scharr= 0;
	fG_Laplace= 0;
	fG_GausLaplace= 0;
	fG_NormLoG= 0;
	fG_KirschN= 0;
	fG_KirschS= 0;
	fG_KirschW= 0;
	fG_KirschE= 0;
	fG_KirschNE= 0;
	fG_KirschSE= 0;
	fG_KirschNW= 0;
	fG_KirschSW= 0;

	Init();

}//close costructor


GradientFilter::~GradientFilter(){

	if(fGx_Sobel) fGx_Sobel->Delete();
	if(fGy_Sobel) fGy_Sobel->Delete();
	if(fGx_Scharr) fGx_Scharr->Delete();
	if(fGy_Scharr) fGy_Scharr->Delete();
	if(fG_Laplace) fG_Laplace->Delete();
	if(fG_GausLaplace) fG_GausLaplace->Delete();
	if(fG_NormLoG) fG_NormLoG->Delete();
	if(fG_KirschN) fG_KirschN->Delete();
	if(fG_KirschS) fG_KirschS->Delete();
	if(fG_KirschW) fG_KirschW->Delete();
	if(fG_KirschE) fG_KirschE->Delete();
	if(fG_KirschNE) fG_KirschNE->Delete();
	if(fG_KirschSE) fG_KirschSE->Delete();
	if(fG_KirschNW) fG_KirschNW->Delete();
	if(fG_KirschSW) fG_KirschSW->Delete();

}//close destructor


void GradientFilter::Init(){

	//## Init kernel matrix
	if(!fGx_Sobel && !fGy_Sobel) {
		int GxKernel[3][3]= { {-1,0,+1}, {-2,0,+2}, {-1,0,+1} };
		int GyKernel[3][3]= { {-1,-2,-1}, {0,0,0}, {+1,+2,+1} };
	
		fGx_Sobel= new TMatrixD(3,3);
		fGy_Sobel= new TMatrixD(3,3);
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				(*fGx_Sobel)(i,j)= GxKernel[i][j];
				(*fGy_Sobel)(i,j)= GyKernel[i][j];
			}
		}
		//cout<<"== Sobel Kernels =="<<endl;
		//fGx_Sobel->Print();
		//fGy_Sobel->Print();
	}//close if

	if(!fGx_Scharr && !fGy_Scharr) {
		int GxKernel[3][3]= { {-3,0,+3}, {-10,0,+10}, {-3,0,+3} };
		int GyKernel[3][3]= { {-3,-10,-3}, {0,0,0}, {+3,+10,+3} };
	
		fGx_Scharr= new TMatrixD(3,3);
		fGy_Scharr= new TMatrixD(3,3);
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				(*fGx_Scharr)(i,j)= GxKernel[i][j];
				(*fGy_Scharr)(i,j)= GyKernel[i][j];
			}
		}
		//cout<<"== Scharr Kernels =="<<endl;
		//fGx_Scharr->Print();
		//fGy_Scharr->Print();
	}//close if

	if(!fG_Laplace) {
		int GKernel[3][3]= { {1, 1, 1}, {1,-8, 1}, {1, 1, 1} };
		
		fG_Laplace= new TMatrixD(3,3);
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				(*fG_Laplace)(i,j)= GKernel[i][j];
			}
		}
		//cout<<"== Laplace Kernels =="<<endl;
		//fG_Laplace->Print();
	}//close if

	if(!fG_GausLaplace) {
		int GKernel[7][7]= { {0,0,1,1,1,0,0}, {0,1,3,3,3,1,0}, {1,3,0,-7,0,3,1}, {1,3,-7,-24,-7,3,1}, {0,0,1,1,1,0,0}, {0,1,3,3,3,1,0}, {1,3,0,-7,0,3,1} };
		
		fG_GausLaplace= new TMatrixD(7,7);
		for(int i=0;i<7;i++){
			for(int j=0;j<7;j++){
				(*fG_GausLaplace)(i,j)= GKernel[i][j];
			}
		}
	}//close if
	
	//## Kirsch kernels
	if(!fG_KirschN){
		int GKernel[3][3]= { {-3,-3,-3}, {-3,0,-3}, {5,5,5} };
		
		fG_KirschN= new TMatrixD(3,3);
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				(*fG_KirschN)(i,j)= GKernel[i][j];
			}
		}
	}//close if
	if(!fG_KirschS){
		int GKernel[3][3]= { {5,5,5}, {-3,0,-3}, {-3,-3,-3} };
		
		fG_KirschS= new TMatrixD(3,3);
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				(*fG_KirschS)(i,j)= GKernel[i][j];
			}
		}
	}//close if
	if(!fG_KirschW){
		int GKernel[3][3]= { {-3,-3,5}, {-3,0,5}, {-3,-3,5} };
		
		fG_KirschW= new TMatrixD(3,3);
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				(*fG_KirschW)(i,j)= GKernel[i][j];
			}
		}
	}//close if
	if(!fG_KirschE){
		int GKernel[3][3]= { {5,-3,-3}, {5,0,-3}, {5,-3,-3} };
		
		fG_KirschE= new TMatrixD(3,3);
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				(*fG_KirschE)(i,j)= GKernel[i][j];
			}
		}
	}//close if
	if(!fG_KirschNW){
		int GKernel[3][3]= { {-3,-3,-3}, {-3,0,5}, {-3,5,5} };
		
		fG_KirschNW= new TMatrixD(3,3);
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				(*fG_KirschNW)(i,j)= GKernel[i][j];
			}
		}
	}//close if
	if(!fG_KirschNE){
		int GKernel[3][3]= { {-3,-3,-3}, {5,0,-3}, {5,5,-3} };
		
		fG_KirschNE= new TMatrixD(3,3);
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				(*fG_KirschNE)(i,j)= GKernel[i][j];
			}
		}
	}//close if

	if(!fG_KirschSW){
		int GKernel[3][3]= { {-3,5,5}, {-3,0,5}, {-3,-3,-3} };
		
		fG_KirschSW= new TMatrixD(3,3);
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				(*fG_KirschSW)(i,j)= GKernel[i][j];
			}
		}
	}//close if
	if(!fG_KirschSE){
		int GKernel[3][3]= { {5,5,-3}, {5,0,-3}, {-3,-3,-3} };
		
		fG_KirschSE= new TMatrixD(3,3);
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				(*fG_KirschSE)(i,j)= GKernel[i][j];
			}
		}
	}//close if

}//close GradientFilter::Init()


TMatrixD* GradientFilter::MakeKernel(int KernelType,int size,double scale){

	TMatrixD* Kernel= 0;

	if(KernelType==Img::eSOBEL_HORIZ){
		Kernel= fGx_Sobel;
	}
	else if(KernelType==Img::eSOBEL_VERT){
		Kernel= fGy_Sobel;
	}
	else if(KernelType==Img::eSCHARR_HORIZ){
		Kernel= fGx_Scharr;
	}
	else if(KernelType==Img::eSCHARR_VERT){
		Kernel= fGy_Scharr;
	}
	else if(KernelType==Img::eLAPLACE){
		Kernel= fG_Laplace;
	}
	else if(KernelType==Img::eLoG){
		Kernel= fG_GausLaplace;
	}
	else if(KernelType==Img::eNormLoG){
		Kernel= new TMatrixD(size,size);	
	
		int nrows= size;
		int ncols= size;
		int halfSize= (size-1)/2;
		for (int i=0;i<nrows;i++){
  		for (int j=0;j<ncols;j++){
  			double x = j - halfSize;
  	  	double y = i - halfSize;
				double kernValue= NormLoGKernel(x,y,scale);
				(*Kernel)(i,j)= kernValue;
			}//end loop j
		}//end loop i
	
		double sum= Kernel->Sum();
		double mean= sum/(size*size);
		(*Kernel)-= mean;

	}//close else if
	else {
		cerr<<"GradientFilter::MakeKernel(): WARNING: Invalid kernel specified!"<<endl;
		return 0;
	}

	return Kernel;

}//close MakeKernel()


Img* GradientFilter::GetFilteredImage(Img* image,int KernelType,int size, double scale){
	
	//## Init kernels
	Init();
	
	//## Get image matrix from image class
	cout<<"GradientFilter::GetFilteredImage(): Getting image matrix ...."<<endl;
	TMatrixD I= *(image->GetMatrix());
	int nrows= I.GetNrows();
	int ncols= I.GetNcols();

	//## Get kernel
	TMatrixD* Kernel= 0;
	TMatrixD* filteredMatrix= 0;
	if(KernelType==Img::eKIRSCH){
		TMatrixD* filteredEnsemble[8];
		for(int k=0;k<8;k++){
			if(k==0) filteredEnsemble[k]= GetConvolution(*fG_KirschN,I);
			else if(k==1) filteredEnsemble[k]= GetConvolution(*fG_KirschS,I);
			else if(k==2) filteredEnsemble[k]= GetConvolution(*fG_KirschW,I);
			else if(k==3) filteredEnsemble[k]= GetConvolution(*fG_KirschE,I);
			else if(k==4) filteredEnsemble[k]= GetConvolution(*fG_KirschNW,I);
			else if(k==5) filteredEnsemble[k]= GetConvolution(*fG_KirschNE,I);
			else if(k==6) filteredEnsemble[k]= GetConvolution(*fG_KirschSW,I);
			else if(k==7) filteredEnsemble[k]= GetConvolution(*fG_KirschSE,I);
		}
		filteredMatrix= new TMatrixD(filteredEnsemble[0]->GetNrows(),filteredEnsemble[0]->GetNcols());
		filteredMatrix->Zero();
		for(int i=0;i<filteredMatrix->GetNrows();i++){
			for(int j=0;j<filteredMatrix->GetNcols();j++){
				double maxVal= -1.e+99;
				for(int k=0;k<8;k++){
					double val= (*filteredEnsemble[k])(i,j); 
					if(val>maxVal) maxVal= val;
				}	
				(*filteredMatrix)(i,j)= maxVal;
			}//end loop cols
		}//end loop rows
	
		for(int k=0;k<8;k++) filteredEnsemble[k]->Delete(); 
	}//close if Kirsch kernel
	else {	
		Kernel= MakeKernel(KernelType,size,scale);
		if(!Kernel){
			cerr<<"GradientFilter::GetFilteredImage(): WARN: Failed to get/compute kernel!"<<endl;
			return 0;
		}	

		filteredMatrix= GetConvolution(*Kernel,I);
		if(!filteredMatrix){
			cerr<<"GradientFilter::GetFilteredImage(): WARN: Failed to compute convolution!"<<endl;
			return 0;
		}
	}//close else

	/*
	TMatrixD* filteredMatrix= 0;
	if(KernelType==Img::eSOBEL_HORIZ){
		filteredMatrix= GetConvolution(*fGx_Sobel,I);
	}
	else if(KernelType==Img::eSOBEL_VERT){
		filteredMatrix= GetConvolution(*fGy_Sobel,I);
	}
	else if(KernelType==Img::eSCHARR_HORIZ){
		filteredMatrix= GetConvolution(*fGx_Scharr,I);		
	}
	else if(KernelType==Img::eSCHARR_VERT){
		filteredMatrix= GetConvolution(*fGy_Scharr,I);
	}
	else if(KernelType==Img::eLAPLACE){
		filteredMatrix= GetConvolution(*fG_Laplace,I);
	}
	else if(KernelType==Img::eLoG){
		filteredMatrix= GetConvolution(*fG_GausLaplace,I);
	}
	else if(KernelType==Img::eKIRSCH){
		TMatrixD* filteredEnsemble[8];
		for(int k=0;k<8;k++){
			if(k==0) filteredEnsemble[k]= GetConvolution(*fG_KirschN,I);
			else if(k==1) filteredEnsemble[k]= GetConvolution(*fG_KirschS,I);
			else if(k==2) filteredEnsemble[k]= GetConvolution(*fG_KirschW,I);
			else if(k==3) filteredEnsemble[k]= GetConvolution(*fG_KirschE,I);
			else if(k==4) filteredEnsemble[k]= GetConvolution(*fG_KirschNW,I);
			else if(k==5) filteredEnsemble[k]= GetConvolution(*fG_KirschNE,I);
			else if(k==6) filteredEnsemble[k]= GetConvolution(*fG_KirschSW,I);
			else if(k==7) filteredEnsemble[k]= GetConvolution(*fG_KirschSE,I);
		}
		filteredMatrix= new TMatrixD(filteredEnsemble[0]->GetNrows(),filteredEnsemble[0]->GetNcols());
		filteredMatrix->Zero();
		for(int i=0;i<filteredMatrix->GetNrows();i++){
			for(int j=0;j<filteredMatrix->GetNcols();j++){
				double maxVal= -1.e+99;
				for(int k=0;k<8;k++){
					double val= (*filteredEnsemble[k])(i,j); 
					if(val>maxVal) maxVal= val;
				}	
				(*filteredMatrix)(i,j)= maxVal;
			}//end loop cols
		}//end loop rows
	
		for(int k=0;k<8;k++) filteredEnsemble[k]->Delete(); 
	}//close else if Kirsch
	else {
		cerr<<"GradientFilter::GetFilteredImage(): WARNING: Invalid kernel specified!"<<endl;
		return 0;
	}
	*/

	//## Convert TMatrixD in image
	TString imgName= Form("%s-Filter%d",image->GetName(),KernelType);
	Img* filteredImage= new Img(imgName,imgName,image->GetNbinsX(),image->GetXaxis()->GetXmin(),image->GetXaxis()->GetXmax(),image->GetNbinsY(),image->GetYaxis()->GetXmin(),image->GetYaxis()->GetXmax());
	
	for(int i=0;i<nrows;i++){
		for(int j=0;j<ncols;j++){
			int ix= j; 
			int iy= i;
			double binContent= image->GetBinContent(ix+1,iy+1);
			if(binContent==0) continue;
			double w= (*filteredMatrix)(i,j);
			filteredImage->SetBinContent(ix+1,iy+1,w);
			//filteredImage->FillPixel(ix,iy,w);
		}//end loop matrix cols
	}//end loop matrix rows

	filteredMatrix->Delete();

	return filteredImage;

}//close GradientFilter::GetFilteredImage()


TMatrixD* GradientFilter::GetConvolution(TMatrixD H,TMatrixD F){
	
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
	//TMatrixD C(nrows,ncols);
	//C.Zero();

	TMatrixD* C= new TMatrixD(nrows,ncols);
	C->Zero();

  
	//Extend matrix at the borders with reflection conditions	
	int rowIndex= (nfilter_rows-1)/2;
	int colIndex= (nfilter_cols-1)/2;
	//int filterGap= pow(2,n-1); 

	
	//Compute convolution
	for(int i=0;i<nrows;i++){		
		for(int j=0;j<ncols;j++){

			for(int k=0;k<nfilter_rows;k++){//start loops on filter box rows
				int index_x= i-k+nfilter_rows/2;
				int mirror_index_x= GetMirrorIndex(index_x,nrows);

				for(int l=0;l<nfilter_cols;l++){//start loops on filter box cols
					int index_y= j-l+nfilter_cols/2;
					int mirror_index_y= GetMirrorIndex(index_y,ncols);	

					double h= H(k,l);
					double f= F(mirror_index_x,mirror_index_y); 
					(*C)(i,j)+= h*f;
				}//end loop filter cols
	
			}//end loop filter rows

		}//end loop y
	}//end loop x

	/*
	for(int i=0;i<nrows;i++){
		for(int j=0;j<ncols;j++){

			for(int k=-rowIndex;k<rowIndex;k++){//start loops on filter box rows
				//int index_x= i + filterGap*k;
				int index_x= i + k;
				int mirror_index_x= GetMirrorIndex(index_x,nrows);
		
				for(int l=-colIndex;l<colIndex;l++){//start loops on filter box cols
					//int index_y= j + filterGap*l;
					int index_y= j + l;
					int mirror_index_y= GetMirrorIndex(index_y,ncols);
		
					double h= H(k+rowIndex,l+colIndex);
					double f= F(mirror_index_x,mirror_index_y); 
					//C(i,j)+= h*f;
					(*C)(i,j)+= h*f;
				}//end loop filter cols
			}//end loop filter rows
		
		}//end loop cols
	}//end loop rows
	*/
	

	return C;

}//close GradientFilter::GetConvolution()


int GradientFilter::GetMirrorIndex(int index,int N){
	
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
		mirror_index= 2*(N-1)-index;
	}
	else{
		cerr<<"GradientFilter::GetMirrorIndex(): ERROR: Invalid index of matrix size passed...exit!"<<endl;
		exit(1);
	}	
	
	return mirror_index;

}//close GetMirrorIndex()

		
double GradientFilter::NormLoGKernel(double x,double y,double sigma){

	double arg= (x*x+y*y)/(2.*sigma*sigma);
	double norm= 1./(TMath::Pi()*pow(sigma,4));
	double fcn= - pow(sigma,2)*norm*(1-arg)*exp(-arg);

	return fcn;

}//close GradientFilter::NormLoGKernel()


