
/**
* @file OutlierDetector.cc
* @class OutlierDetector
* @brief OutlierDetector
*
* Detect outlier in a multidimensional data table. Based on robust Mahalanobis distance
* @author S. Riggi
* @date 26/06/2015
*/


#include <OutlierDetector.h>

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

//ClassImp(OutlierDetector)

//RInside OutlierDetector::fR;

OutlierDetector::OutlierDetector(){

	fOutliersIds.clear();
	fDataFlags.clear();
	fDataDistances.clear();
	fRobustMean= 0;
	fRobustCov= 0;
	//fCutoffValues= 0;

}//close costructor



OutlierDetector::~OutlierDetector(){
	

}//close destructor



int OutlierDetector::FindOutliers(TMatrixD* data){

	//## Check input data
	if(!data) {
		cerr<<"OutlierDetector::FindOutliers(): ERROR: Null ptr to data matrix given!"<<endl;
		return -1;
	}

	//## Wrapper to call R commands from C++ code
	try{

		RInside* fR= RInside::instancePtr();
		if(!fR){
			cerr<<"OutlierDetector::FindOutliers(): ERROR: Cannot retrieve RInside instance (did you create one in your main application?)!"<<endl;
			fR= new RInside;
			//return -1;
		} 

		//## Load R library for outlier detection
		cout<<"OutlierDetector::FindOutliers(): INFO: Loading needed R packages..."<<endl;
		std::string RCmd= std::string("library(\"rrcovHD\");");
		fR->parseEvalQ(RCmd);

		//## Import data matrix in R environment
		cout<<"OutlierDetector::FindOutliers(): INFO: Importing data matrix in R environment..."<<endl;
		int nDim= data->GetNcols();
		int nEntries= data->GetNrows();
		Rcpp::NumericMatrix data_matrix(nEntries,nDim);
	
		for(int i=0;i<nEntries;i++){
			for(int j=0;j<nDim;j++){
				data_matrix(i,j)= (*data)(i,j);
			}
		}
		(*fR)["data"]= data_matrix;
		//fR->parseEvalQ("print(data);");

		//### Run Mahalanobis outlier algorithm
		cout<<"OutlierDetector::FindOutliers(): INFO: Running Mahalanobis outlier algo..."<<endl;
		RCmd= std::string("outlierResults <- OutlierMahdist(data,trace=TRUE); outlierIds <- getOutliers(outlierResults); dataFlags<-getFlag(outlierResults); distances<-getDistance(outlierResults); cutoffValues<-getCutoff(outlierResults);");
		fR->parseEval(RCmd);

		cout<<"OutlierDetector::FindOutliers(): INFO: Get robust covariance and mean estimates..."<<endl;
		RCmd= std::string("covInfo<-slot(outlierResults,\"covobj\"); robustMean<-attributes(covInfo[[1]])$center; robustCov<-attributes(covInfo[[1]])$cov");
		fR->parseEval(RCmd);

		//fR->parseEvalQ("print(outlierIds)");
		//fR->parseEvalQ("print(dataFlags)");
		//fR->parseEvalQ("print(distances)");

		//## Importing results
		//## Outliers found by this package are based on a Chi2 test
		cout<<"OutlierDetector::FindOutliers(): INFO: Importing results from R environment..."<<endl;
		Rcpp::NumericVector outlierIds = fR->parseEval(std::string("outlierIds"));
		Rcpp::NumericVector dataFlags = fR->parseEval(std::string("dataFlags"));
		Rcpp::NumericVector dataDistances = fR->parseEval(std::string("distances"));
		double MD_mean= Rcpp::as<double>( fR->parseEval(std::string("mean(distances)")) );
		double MD_sigma= Rcpp::as<double>( fR->parseEval(std::string("sqrt(var(distances))")) );
		Rcpp::NumericVector robustMean = fR->parseEval(std::string("robustMean"));
		Rcpp::NumericMatrix robustCov = fR->parseEval(std::string("robustCov"));
		fCutoffValue= Rcpp::as<double>( fR->parseEval(std::string("cutoffValues")) );
		

		fOutliersIds.clear();
		fDataFlags.clear();
		fDataDistances.clear();
		
		//cout<<"== LIST OF OUTLIER IDS =="<<endl;
		//cout<<outlierIds.size()<<" outliers found (";
		for(int k=0;k<outlierIds.size();k++){
			fOutliersIds.push_back(outlierIds(k));
			//cout<<outlierIds(k)<<", ";
		}//end loop
		//cout<<")"<<endl;

		//cout<<"dataFlags (";
		for(int k=0;k<dataFlags.size();k++){
			fDataFlags.push_back(dataFlags(k));
			//cout<<dataFlags(k)<<", ";
		}//end loop
		//cout<<")"<<endl;

		//cout<<"== MAHALANOBIS DISTANCE INFO =="<<endl;
		//cout<<"dataDistances (";
		for(int k=0;k<dataDistances.size();k++){
			fDataDistances.push_back(dataDistances(k));
			//cout<<dataDistances(k)<<", ";
		}//end loop
		//cout<<")"<<endl;	
		//cout<<"MD_mean="<<MD_mean<<", MD_sigma="<<MD_sigma<<endl;
		//cout<<"=========================="<<endl;
			
		cout<<"OutlierDetector::FindOutliers(): INFO: Filling vector/matrix..."<<endl;
		
		if(fRobustMean) fRobustMean->Delete();
		fRobustMean= new TVectorD(nDim);
		if(fRobustCov) fRobustCov->Delete();
		fRobustCov= new TMatrixD(nDim,nDim);
		//if(fCutoffValues) fCutoffValues->Delete();
		//fCutoffValues= new TVectorD(nDim);
		
		for(int i=0;i<nDim;i++){
			(*fRobustMean)(i)= robustMean(i);
			//(*fCutoffValues)(i)= cutoffValues(i);
			for(int j=0;j<nDim;j++){
				(*fRobustCov)(i,j)= robustCov(i,j);
			}//end loop 
		}//end loop 

		/*
		//## Find the outliers with a Grubbs test
		//## Now search for excesses in the mahalanobis distance
		int N= (int)dataDistances.size();
		double alpha= 0.05;//significance level (conventional 5%)
		double ndf= N-2;//degrees of freedom to be used for the Grubb test
		double alpha_test= alpha/(2*N);//significance to be used for the Grubb test (two-sided)
		//double alpha_test= alpha/N;//one-sided 
		double tStudentCriticalValue= TMath::StudentQuantile(alpha,ndf,false);//taking the upper critical values at significance alpha
		double GrubbTestStatistics= (N-1)/sqrt(N)* sqrt( pow(tStudentCriticalValue,2)/(N-1+pow(tStudentCriticalValue,2)) );
		
		fOutliersIds.clear();
		fDataFlags.clear();
		fDataFlags.assign(N,1);
		for(int i=0;i<N;i++){
			double MahaDist= fDataDistances[i];
			double Z= fabs(MahaDist-MD_mean)/MD_sigma;
			if(Z>GrubbTestStatistics){//mark this data as outlier
				fOutliersIds.push_back(i);
				fDataFlags[i]= 0;
			}
		}//end loop data
		*/
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

}//close FindOutliers()



