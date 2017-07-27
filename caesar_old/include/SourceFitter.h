/**
* @file SourceFitter.cc
* @class SourceFitter
* @brief SourceFitter
*
* Class to fit a source image with a mixture of gaussian/skew normal/skew-t bivariate functions
* @author S. Riggi
* @date 01/09/2015
*/

#ifndef SourceFitter_h
#define SourceFitter_h 1

#include <Img.h>

#include <RInside.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TH3D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TF1.h>
#include <TF12.h>
#include <TF2.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TStyle.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TApplication.h>

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
#include <iostream>
#include <time.h>
#include <ctime>

using namespace std;


class SourceFitter {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    SourceFitter();
		/** 
		\brief Class destructor
 		*/
  	~SourceFitter();

	public:

		static int FitSource(Source* aSource,double curvThr,int componentMinNPix);
		
	private:
		static double Gaus2DGeomFcn(double* x, double* par);
		static double Gaus2DMixtureGeomFcn(double* x, double* p);

	private:
		
		
};

#ifdef __MAKECINT__
#pragma link C++ class SourceFitter+; 
#endif

#endif

