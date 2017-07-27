/**
* @file Utils.h
* @class Utils
* @brief Utility functions
*
* Utility functions
* @author S. Riggi
* @date 23/08/2010
*/



#include <Utils.h>

//#include <wcs.h>

#include <RInside.h>
#include <Rcpp.h>

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

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <ctime>

using namespace std;

Utils::Utils(){

}

Utils::~Utils(){

}

std::vector<double> Utils::ComputePageRank(TMatrixD M,double d,double tol){

	bool hasConverged= true;
	int maxIterToStop= 1000;
	int N= M.GetNcols();	

	//cout<<"== GOOGLE RANK ADJACENCY =="<<endl;
	//M.Print();
	//cout<<"===="<<endl;
	
	//Initialize rank with uniform prob (0.5)
	TVectorD v(N);
	TVectorD v_last(N);
	TMatrixD UnoMatrix(N,N);	
	UnoMatrix.Zero();
	for(int i=0;i<N;i++) {
		v(i)= 0.5;
		v_last(i)= 1.e+99;//TMath::Infinity();
		for(int j=0;j<N;j++) UnoMatrix(i,j)= 1;
	}
	
	double norm2= v.Norm2Sqr();
	v*= 1./sqrt(norm2);
	
	//Compute Mhat= (d x M) + (1-d)/N*UnoMatrix(N,N)
	TMatrixD M_hat(N,N);
	M_hat.Zero();
	M_hat= d*M + ((1.-d)/(double)(N))*UnoMatrix;
	
	int iter= 0;
	int nIterToConverge= 0;
	double diff= (v-v_last).Norm2Sqr();
	
	while( diff>tol ){
		v_last= v;
    v = M_hat*v;
		double normFactor= sqrt(v.Norm2Sqr());
    v*= 1./normFactor;
		
		diff= sqrt((v-v_last).Norm2Sqr());	
		iter++;
		nIterToConverge= iter;

    if(iter >= maxIterToStop){
    	hasConverged= false;
			break;
		}
	}//end loop 
	
	std::vector<double> ranks;
	ranks.clear();

	if(hasConverged){
		//cout<<"== GOOGLE RANKS =="<<endl;
		//cout<<"nIter="<<nIterToConverge<<" p(";
		for(int k=0;k<v.GetNoElements();k++) {
			ranks.push_back(v[k]);
			//cout<<v[k]<<",";
		}
		//cout<<")"<<endl;
		//cout<<"=================="<<endl;
	}
	else{
		cerr<<"Utils::ComputePageRank(): WARN: Page rank did not converge!"<<endl;
	}

	return ranks;

}//close Utils::ComputePageRank()


int Utils::GetSaliency(TVectorD& saliencyVect,TMatrix* data,int kn){

	if(!data) return -1;

	/*
	try{

		RInside* fR= RInside::instancePtr();
		if(!fR){
			cerr<<"Utils::GetSaliency(): ERROR: Cannot retrieve RInside instance (did you create one in your main application?)!"<<endl;
			fR= new RInside;
		} 

		//## Load R library for outlier detection
		cout<<"Utils::GetSaliency(): INFO: Loading needed R packages..."<<endl;
		std::string RCmd= std::string("library(\"rrcovHD\");library(\"FNN\")");
		fR->parseEvalQ(RCmd);

		//## Import data matrix in R environment
		cout<<"Utils::GetSaliency(): INFO: Importing data matrix in R environment..."<<endl;
		int nDim= data->GetNcols();
		int nEntries= data->GetNrows();
		Rcpp::NumericMatrix data_matrix(nEntries,nDim);
	
		for(int i=0;i<nEntries;i++){
			for(int j=0;j<nDim;j++){
				data_matrix(i,j)= (*data)(i,j);
			}
		}
		(*fR)["data"]= data_matrix;

		//## Run nearest neighbor and get indexes
		(*fR)["kn"]= kn;
		cout<<"Utils::GetSaliency(): INFO: Running nearest-neighbor search algo..."<<endl;
		RCmd= std::string("NNIndex<-knn.index(x,k=kn);");
		fR->parseEval(RCmd);
		Rcpp::NumericMatrix nnIndexes = fR->parseEval(std::string("NNIndex"));
		
		//## Compute saliency vector
		saliencyVect.Zero();
		for(int i=0;i<nnIndexes.nrow();i++){
			for(int j=0;j<nnIndexes.ncol();j++){
				int index= nnIndexes(i,j);
				double S= 
				saliencyVect(i)+= 
			}//end loop cols
		}//end loop rows
	
	}//close try block
	catch( std::exception &ex ) {
		cerr << "OutlierDetector::FindOutliers(): ERROR: Exception catched: " << ex.what() << endl;
		return -1;
  } 
	catch(...) { 
		cerr << "OutlierDetector::FindOutliers(): ERROR: C++ exception (unknown reason)" << endl;
		return -1;
  }	
	
	return 0;
	*/

}//close GetSaliency()


int Utils::GetEntropy(double& H,TMatrixD* data, int k){

	if(!data) return -1;
	H= 0;
	int N= data->GetNrows();
	int p= data->GetNcols();

	//## If less than 2 values are present in the data set entropy to E=0
	if(N<2){
		cerr<<"Utils::GetEntropy(): WARN: <2 values in data...returning zero entropy!"<<endl;
		H= 0;
		return 0;
	}

	//## Check k vs N
	if(k>=N){
		cerr<<"Utils::GetEntropy(): WARN: Selected k>N...setting k=min(5,N-1)!"<<endl;
		k= std::min(5,N-1);
	}
	
	try{
		RInside* fR= RInside::instancePtr();
		if(!fR){
			cerr<<"Utils::GetEntropy(): ERROR: Cannot retrieve RInside instance (did you create one in your main application?)!"<<endl;
			fR= new RInside;
		} 

		//## Load R library for outlier detection
		std::string RCmd= std::string("library(\"FNN\");");
		fR->parseEvalQ(RCmd);

		//## Import data matrix in R environment
		cout<<"Utils::GetEntropy(): INFO: Importing data in R"<<endl;
		data->Print();
		Rcpp::NumericMatrix X(N,p);
		for(int i=0;i<N;i++){
			for(int j=0;j<p;j++){
				double value= (*data)(i,j);
				X(i,j)= value;
			}
		}
		(*fR)["X"]= X;	
		(*fR)["k"]= k;
		fR->parseEval(std::string("print(X)"));

		//### Compute Entropy params
		cout<<"Utils::GetEntropy(): INFO: Computing entropy..."<<endl;
		RCmd= std::string("HList <- entropy(X,k,\"kd_tree\");");
		fR->parseEval(RCmd);
	
		cout<<"Utils::GetEntropy(): INFO: Retrieve results..."<<endl;
		//RCmd= std::string("H <- HList[k];");
		//fR->parseEval(RCmd);

		//## Importing results		
		H= Rcpp::as< double >(fR->parseEval(std::string("HList[k]")));
		cout<<"Utils::GetEntropy(): INFO: H="<<H<<endl;

	}//close try block
	catch( std::exception &ex ) {
		cerr << "Utils::GetEntropy(): ERROR: Exception catched: " << ex.what() << endl;
		return -1;
  } 
	catch(...) { 
		cerr << "Utils::GetEntropy(): ERROR: C++ exception (unknown reason)" << endl;
		return -1;
  }	

	return 0;

}//close Utils::GetEntropy()


bool Utils::PixelToWCSCoords(Img* image,int ix,int iy,double& xpos, double& ypos) {

	/*
	//Check pixel values in input
	if( !image || ix<0 || iy<0 || ix>=image->GetNbinsX() || iy>=image->GetNbinsY() )
		return false;	
	
	//Check image meta-data
	if(!image->HasMetaData() ){
    cerr<<"Utils::PixelToWCSCoords: WARNING: No metadata available in image!"<<endl;
		return false;
	}
	Img::MetaData metadata= image->GetMetaData();

	//Set a coord system
	WorldCoor* wcs= 0;
	try {
		wcs= wcskinit(	
			metadata.Nx,metadata.Ny,
			(char*)metadata.CoordTypeX.c_str(),(char*)metadata.CoordTypeY.c_str(),
			metadata.Cx,metadata.Cy,
			metadata.Xc,metadata.Yc,
			NULL,
			metadata.dX,metadata.dY,
			metadata.RotY,
			(int)(metadata.Epoch),
			metadata.Epoch
		);
	}//close try block
	catch(std::exception const & e) {
		cerr<<"Utils::PixelToWCSCoords(): ERROR: WCS coordinate system initialization failed (err="<<e.what()<<")";
    return false;
	}

	//pix2wcs (wcs,ix+1,iy+1,&xpos, &ypos);	
	pix2wcs (wcs,ix,iy,&xpos, &ypos);
	*/


	/*
	//Test conversion
	for(int j=0;j<image->GetNbinsY();j++){
		double iy= image->GetYaxis()->GetBinCenter(j+1);
		for(int i=0;i<image->GetNbinsX();i++){
			double ix= image->GetXaxis()->GetBinCenter(i+1);
			double xpos, ypos;
			pix2wcs (wcs,ix+1,iy+1,&xpos, &ypos);
			//cout<<"Pixel ("<<ix<<","<<iy<<") ";
			//cout<<std::setprecision(4)<<fixed<<"Pos("<<xpos<<","<<ypos<<")"<<endl;
		}//end loop bins X
	}//end loop bins Y
  */

	return true;
		
}//close PixelToWCSCoords()


