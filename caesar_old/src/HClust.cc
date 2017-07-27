

/**
* @file HClust.cc
* @class HClust
* @brief HClust
*
* Find hierarchical clustering
* @author S. Riggi
* @date 26/06/2015
*/



#include <HClust.h>

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


HClust::HClust(){


}//close costructor



HClust::~HClust(){
	

}//close destructor



int HClust::FindClusters(TMatrixD* DissMatrix,int aggloMethod,int minClustSize,double maxHeightQ,int deepSplitLevel){

	//## Check input data
	if(!DissMatrix) {
		cerr<<"HClust::FindClusters(): ERROR: Null ptr to dissimilarity matrix given!"<<endl;
		return -1;
	}

	//## Wrapper to call R commands from C++ code
	try{

		RInside* fR= RInside::instancePtr();
		if(!fR){
			cerr<<"HClust::FindClusters(): ERROR: Cannot retrieve RInside instance (did you create one in your main application?)!"<<endl;
			fR= new RInside;
		} 

		//## Load R library for outlier detection
		cout<<"HClust::FindClusters(): INFO: Loading needed R packages..."<<endl;
		std::string RCmd= std::string("library(\"dynamicTreeCut\")");
		fR->parseEvalQ(RCmd);

		//## Import data matrix in R environment
		cout<<"HClust::FindClusters(): INFO: Importing data matrix in R environment..."<<endl;
		int nCols= DissMatrix->GetNcols();
		int nRows= DissMatrix->GetNrows();
		Rcpp::NumericMatrix diss_matrix(nRows,nCols);
	
		for(int i=0;i<nRows;i++){
			for(int j=0;j<nCols;j++){
				diss_matrix(i,j)= (*DissMatrix)(i,j);
			}
		}
		(*fR)["diss_matrix"]= diss_matrix;

		//### Convert diss matrix into a R diss object
		cout<<"HClust::FindClusters(): INFO: Importing data matrix in R dist object..."<<endl;
		RCmd= std::string("d<-as.dist(as.matrix(diss_matrix));");
		fR->parseEval(RCmd);

		//### Run hierarchical cluster algorithm
		cout<<"HClust::FindClusters(): INFO: Running hierarchical clustering algo..."<<endl;
		if(aggloMethod==eWard)
			RCmd= std::string("hc <- hclust(d=d, method=\"ward.D2\");");
		else if(aggloMethod==eSingleLinkage)
			RCmd= std::string("hc <- hclust(d=d, method=\"single\");");
		else if(aggloMethod==eCompleteLinkage)
			RCmd= std::string("hc <- hclust(d=d, method=\"complete\");");
		else if(aggloMethod==eAverageLinkage)
			RCmd= std::string("hc <- hclust(d=d, method=\"average\");");
		else{
			std::string errMsg= "Invalid agglomation method selected!";
			cerr<<"HClust::FindClusters(): ERROR: "<<errMsg<<endl;
			throw std::invalid_argument(errMsg.c_str());
		}
		fR->parseEval(RCmd);

		//## Get heights of dendrogram
		(*fR)["maxHeightQ"]= maxHeightQ;
		RCmd= std::string("heights <- hc$height; maxTreeHeight<-quantile(heights,maxHeightQ);");
		fR->parseEval(RCmd);
		
		//## Get and draw dendrogram
		RCmd= std::string("dendro <- as.dendrogram(hc); plot(dendro);");
		fR->parseEval(RCmd);

		//### Apply dynamic cut to hierarchy tree
		cout<<"HClust::FindClusters(): INFO: Applying dynamic cut to hierarchy tree..."<<endl;
		(*fR)["minClustSize"]= minClustSize;
		(*fR)["deepSplitLevel"]= deepSplitLevel;
		if(deepSplitLevel==0) RCmd= std::string("cutRes<-cutreeDynamicTree(hc,maxTreeHeight,deepSplit=FALSE,minModuleSize=minClustSize);");
		else RCmd= std::string("cutRes<-cutreeDynamicTree(hc,maxTreeHeight,deepSplit=TRUE,minModuleSize=minClustSize);");

		//RCmd= std::string("cutRes<-cutreeDynamic(dendro=hc,method=\"hybrid\",minClusterSize=minClustSize,distM=as.matrix(diss_matrix),deepSplit=deepSplitLevel);");
		//fR->parseEval(RCmd);


		//## Retrieve results
		cout<<"HClust::FindClusters(): INFO: Retrieving results..."<<endl;
		Rcpp::NumericVector clusterIds = fR->parseEval(std::string("cutRes"));
		fClusterIds.clear();
		for(unsigned int k=0;k<clusterIds.size();k++){
			fClusterIds.push_back(clusterIds(k));
		}//end loop
		
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



