/**
* @file ZernikeMoments.h
* @class ZernikeMoments
* @brief ZernikeMoments
*
* @author S. Riggi
* @date 15/06/2015
*/



#ifndef ZERNIKE_MOMENTS_H
#define ZERNIKE_MOMENTS_H

#include "Img.h"
#include "Region.h"
#include "Contour.h"

//#include <mathop.h>

#include <TVector2.h>
#include <TGraph.h>
#include <TText.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include <stdio.h>
#include <math.h>
#include <vector>
#include <float.h>
using namespace std;


class ZernikeMoments {

	public:
  	/**
		\brief Class constructor
		*/
    ZernikeMoments();
		/**
		\brief Class destructor
		*/
    ~ZernikeMoments();
		
		
	public:
    
		static std::vector<double> GetZernike2D_Direct(Img* img, double order, double radius);
		static std::vector<double> GetZernike2D (Img* img, double order, double rad);
		static std::vector<double> GetZernike2DOld(Img* img, double D, double R);
		static std::vector<double> mb_Znl(double *X, double *Y, double *P, int size, double D, double m10_m00, double m01_m00, double R, double psum);

	private:
		int factorial(int n) {
  		return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
		}
		
		
	private:
  	
		
};

#endif

