
/**
* @file Interpolator.cc
* @class Interpolator
* @brief Interpolator
*
* Find data interpolation
* @author S. Riggi
* @date 26/06/2015
*/


#include <Interpolator.h>

#include <Rcpp.h>
using namespace Rcpp ;

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


#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
using namespace cv;

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>

using namespace std;

ClassImp(Interpolator)


Interpolator::Interpolator(){

	fInterpolatedImg= 0;

}//close costructor



Interpolator::~Interpolator(){
	
	//if(fInterpolatedImg) fInterpolatedImg->Delete();

}//close destructor



int Interpolator::FindInterpolation(TMatrixD* data,Img* inputImg){

	//## Check input data
	if(!data || !inputImg) {
		cerr<<"Interpolator::FindInterpolation(): ERROR: Null ptr to data matrix given!"<<endl;
		return -1;
	}


	//## Wrapper to call R commands from C++ code
	try{

		RInside* fR= RInside::instancePtr();
		if(!fR){
			cerr<<"Interpolator::FindInterpolation(): ERROR: Cannot retrieve RInside instance (did you create one in your main application?)!"<<endl;
			fR= new RInside;
			//return -1;
		} 

		//## Load R library for interpolation
		cout<<"Interpolator::FindInterpolation(): INFO: Loading needed R packages..."<<endl;
		std::string RCmd= std::string("library(\"akima\");");
		fR->parseEvalQ(RCmd);

		//## Import data matrix in R environment
		cout<<"Interpolator::FindInterpolation(): INFO: Importing data matrix in R environment..."<<endl;
		int nDim= data->GetNcols();
		int nEntries= data->GetNrows();
		Rcpp::NumericVector xCoords(nEntries);
		Rcpp::NumericVector yCoords(nEntries);
		Rcpp::NumericVector zCoords(nEntries);
		for(int i=0;i<nEntries;i++){
			xCoords(i)=  (*data)(i,0);
			yCoords(i)=  (*data)(i,1);
			zCoords(i)=  (*data)(i,2);
			cout<<"(x,y,z)= ("<<xCoords(i)<<","<<yCoords(i)<<","<<zCoords(i)<<")"<<endl;
		}

		int Nx= inputImg->GetNbinsX();
		int Ny= inputImg->GetNbinsY();
		double minX= inputImg->GetXaxis()->GetXmin();
		double maxX= inputImg->GetXaxis()->GetXmax();
		double minY= inputImg->GetYaxis()->GetXmin();
		double maxY= inputImg->GetYaxis()->GetXmax();
		if(!fInterpolatedImg){
			TString imgName= Form("%s_interp",inputImg->GetName());
			fInterpolatedImg= (Img*)inputImg->Clone(imgName);
			fInterpolatedImg->SetNameTitle(imgName,imgName);
		}
		fInterpolatedImg->Reset();

 		
		(*fR)["xCoords"]= xCoords;
		(*fR)["yCoords"]= yCoords;
		(*fR)["zCoords"]= zCoords;
		(*fR)["Nx"]= Nx;
		(*fR)["Ny"]= Ny;
		(*fR)["minX"]= minX;
		(*fR)["maxX"]= maxX;
		(*fR)["minY"]= minY;
		(*fR)["maxY"]= maxY;
		fR->parseEval(std::string("print(cbind(xCoords,yCoords,zCoords))"));
		
		//### Run Mahalanobis outlier algorithm
		cout<<"Interpolator::FindInterpolation(): INFO: Running interpolation algo..."<<endl;
		RCmd= std::string("interpRes<-interp(x=xCoords,y=yCoords,z=zCoords,xo=seq(minX,maxX,1),yo=seq(minY,maxY,1),linear=FALSE,extrap=TRUE);");
		fR->parseEval(RCmd);

		//interp(x, y=NULL, z, xo=seq(min(x), max(x), length = nx),yo=seq(min(y), max(y), length = ny),
		//		linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL,
		//		ncp = NULL, nx = 40, ny = 40)

		//## Getting results
		cout<<"Interpolator::FindInterpolation(): INFO: Retrieve results..."<<endl;
		Rcpp::NumericMatrix interpZ = fR->parseEval(std::string("interpRes$z;"));
		fR->parseEval(std::string("print(cbind(xCoords,yCoords,zCoords,interpRes$z));"));

		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				double w= interpZ(i,j);
				fInterpolatedImg->SetBinContent(i+1,j+1,w);
			}//end loop bins Y
		}//end loop bins X

		cout<<"minX="<<minX<<" maxX="<<maxX<<" minY="<<minY<<" maxY="<<maxY<<endl;
		
	}//close try block
	catch( std::exception &ex ) {
		cerr << "Interpolator::FindInterpolation(): ERROR: Exception catched: " << ex.what() << endl;
		return -1;
  } 
	catch(...) { 
		cerr << "Interpolator::FindInterpolation(): ERROR: C++ exception (unknown reason)" << endl;
		return -1;
  }	
	
	return 0;

}//close FindInterpolation()



