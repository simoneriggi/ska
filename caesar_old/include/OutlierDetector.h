/**
* @file OutlierDetector.h
* @class OutlierDetector
* @brief OutlierDetector
*
* Detect outlier in a multidimensional data table. Based on robust Mahalanobis distance
* @author S. Riggi
* @date 26/06/2015
*/



#ifndef OutlierDetector_h
#define OutlierDetector_h 1

#include <RInside.h>  // for the embedded R via RInside

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TF1.h>
#include <TF12.h>
#include <TF2.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TStyle.h>
#include <TGraph.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

#include <vector>
#include <algorithm>
#include <map>
#include <string>

using namespace std;


class OutlierDetector {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    OutlierDetector();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~OutlierDetector();

	public:
		int FindOutliers(TMatrixD*);

		std::vector<int> GetOutliers(){return fOutliersIds;}
		std::vector<int> GetDataFlags(){return fDataFlags;}
		std::vector<double> GetDistances(){return fDataDistances;}
		TVectorD* GetRobustMean(){return fRobustMean;}
		TMatrixD* GetRobustCov(){return fRobustCov;}
		//TVectorD* GetCutoffValues(){return fCutoffValues;}
		double GetCutoff(){return fCutoffValue;}

	private:
	
		//static RInside fR;
		std::vector<int> fOutliersIds;
		std::vector<int> fDataFlags;
		std::vector<double> fDataDistances;
		TVectorD* fRobustMean;
		TMatrixD* fRobustCov;
		//TVectorD* fCutoffValues;
		double fCutoffValue;
		
};


#endif



