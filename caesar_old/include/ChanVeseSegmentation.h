/**
* @file ChanVeseSegmentation.h
* @class ChanVeseSegmentation
* @brief ChanVeseSegmentation
*
* @author S. Riggi
* @date 15/06/2015
*/

#ifndef CHANVESE_SEGMENTATION_H
#define CHANVESE_SEGMENTATION_H


#include <Img.h>
#include <Contour.h>

#include <cmath>


//#include <mathop.h>

#include <TVector2.h>
#include <TGraph.h>
#include <TText.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TCanvas.h>

#include <stdio.h>
#include <math.h>
#include <vector>
#include <float.h>
using namespace std;


class ChanVeseSegmentation {

	public:
  	/**
		\brief Class constructor
		*/
    ChanVeseSegmentation();
		/**
		\brief Class destructor
		*/
    ~ChanVeseSegmentation();
		
		// Define structure containing parameters of
		struct CVsetup {
  		double dt; // time step 
  		double h;  // pixel spacing
  		double lambda1;
  		double lambda2;
  		double mu; // contour length weighting parameter
  		double nu; // region area weighting parameter
  		unsigned int p; // length weight exponent
		};


	public:
		int RunSegmentation(Img* img,double dt=0.1,double h=1,double lambda1=1.0,double lambda2=2.0,double mu=0.5,double nu=0,double p=1,double initContourRadius=1);
			 
		Img* GetSegmentedImage(){return fSegmentedImg;}
		Contour* GetContour(){return fContour;}
		Img* GetContourImage(){return fContourImg;}

	private:
		// Main segmentation algorithm
		void CVSegmentation(TMatrixD* img,TMatrixD* phi0,TMatrixD** phi,struct CVsetup* pCVinputs);
    
		// Compute gray level averages in foreground and background regions defined by level set function phi
		void GetRegionAverages(TMatrixD* img, TMatrixD* phi,double epsilon,double &c1,double &c2);

		// Compute coefficients needed in Chan-Vese segmentation algorithm given current level set function
		void GetChanVeseCoefficients(TMatrixD* phi,struct CVsetup* pCVinputs,
														 unsigned int i,
														 unsigned int j,
                             double L,
                             double& F1,
                             double& F2,
                             double& F3,
                             double& F4,
                             double& F,
                             double& deltaPhi);
                             
                             
		// Reinitialize a function to the signed distance function to its zero contour
		void ReinitPhi(TMatrixD* phiIn,TMatrixD** psiOut,double dt,double h,unsigned int numIts);

		void ZeroCrossings(TMatrixD* imageIn,TMatrixD** edges,double fg,double bg);

	
	private:

		Contour* fContour;
		Img* fSegmentedImg;
		Img* fContourImg;
		TCanvas* fPlot;
		TCanvas* fInitPlot;

};//close class
		 

#ifdef __MAKECINT__
#pragma link C++ class ChanVeseSegmentation+;
#endif
                    
#endif


