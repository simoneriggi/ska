/**
* @file WTFilter.h
* @class WTFilter
* @brief WTFilter
*
* Wavelet Filter
* @author S. Riggi
* @date 20/01/2015
*/



#ifndef WTFilter_h
#define WTFilter_h 1

#include "Img.h"

#include <TH1D.h>
#include <TH2D.h>
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

class WTFilter {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    WTFilter();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~WTFilter();

		//typedef std::vector<Img> ImgCollection;

	public:
	
		static std::vector<Img*> GetDecomposition(Img* image,int nScales);
		
	private:

		static void Init();
		static TMatrixD GetConvolution(TMatrixD H,TMatrixD F,int scaleId);
		static int GetMirrorIndex(int index,int N);

	private:

		static TMatrixD* fH;
		bool isInitCalled;

};

#endif
