/**
* @file GradientFilter.h
* @class GradientFilter
* @brief GradientFilter
*
* Gradient Filter
* @author S. Riggi
* @date 20/01/2015
*/



#ifndef GradientFilter_h
#define GradientFilter_h 1

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

class Img;

class GradientFilter {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    GradientFilter();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~GradientFilter();

		
	public:
	
		static Img* GetFilteredImage(Img* img,int KernelType,int size=3, double scale=1);
		
	private:

		static void Init();
		static TMatrixD* GetConvolution(TMatrixD H,TMatrixD F);
		static int GetMirrorIndex(int index,int N);
		static TMatrixD* MakeKernel(int KernelType,int size=3,double scale=1);
		static double NormLoGKernel(double x,double y,double sigma);

	private:

		static TMatrixD* fGx_Sobel;
		static TMatrixD* fGy_Sobel;
		static TMatrixD* fGx_Scharr;
		static TMatrixD* fGy_Scharr;
		static TMatrixD* fG_Laplace;
		static TMatrixD* fG_GausLaplace;
		static TMatrixD* fG_NormLoG;
		
		static TMatrixD* fG_KirschN;
		static TMatrixD* fG_KirschS;
		static TMatrixD* fG_KirschW;
		static TMatrixD* fG_KirschE;
		static TMatrixD* fG_KirschNE;
		static TMatrixD* fG_KirschSE;
		static TMatrixD* fG_KirschNW;
		static TMatrixD* fG_KirschSW;
		

		bool isInitCalled;

};

#endif
