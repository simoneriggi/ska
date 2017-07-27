/**
* @file Interpolator.h
* @class Interpolator
* @brief Interpolator
*
* Perform data interpolation
* @author S. Riggi
* @date 26/06/2015
*/



#ifndef Interpolator_h
#define Interpolator_h 1

#include "Img.h"

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


class Interpolator {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    Interpolator();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~Interpolator();

	public:
		int FindInterpolation(TMatrixD*,Img*);
		
		Img* GetInterpolatedMap(){return fInterpolatedImg;}

	private:
	
		Img* fInterpolatedImg;
		
};


#endif



