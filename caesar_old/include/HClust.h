/**
* @file HClust.h
* @class HClust
* @brief HClust
*
* Find hierarchical clustering
* @author S. Riggi
* @date 26/06/2015
*/



#ifndef HClust_h
#define HClust_h 1

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


class HClust {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    HClust();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~HClust();

		enum AggloMethod {eWard=1,eSingleLinkage=2,eCompleteLinkage=3,eAverageLinkage=4};

	public:
		int FindClusters(TMatrixD*,int aggloMethod=2,int minClustSize=2,double maxHeightQ=0.99,int deepSplitLevel=1);
		std::vector<int> GetClusterIds(){return fClusterIds;}

	private:
		std::vector<int> fClusterIds;
		
};


#endif



