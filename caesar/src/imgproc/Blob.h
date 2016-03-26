// ***********************************************************************
// * License and Disclaimer                                              *
// *                                                                     *
// * Copyright 2016 Simone Riggi																			   *
// *																																	   *
// * This file is part of Caesar. 																		   *
// * Caesar is free software: you can redistribute it and/or modify it   *
// * under the terms of the GNU General Public License as published by   *
// * the Free Software Foundation, either * version 3 of the License,    *
// * or (at your option) any later version.                              *
// * Caesar is distributed in the hope that it will be useful, but 			 *
// * WITHOUT ANY WARRANTY; without even the implied warranty of          * 
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                *
// * See the GNU General Public License for more details. You should     * 
// * have received a copy of the GNU General Public License along with   * 
// * Caesar. If not, see http://www.gnu.org/licenses/.                   *
// ***********************************************************************
/**
* @file Blob.h
* @class Blob
* @brief Blob class
*
* Class representing an image blob
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef Blob_h
#define Blob_h 1

#include <Pixel.h>
#include <Img.h>

#include <TObject.h>
#include <TMatrixD.h>
#include <TColor.h>

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

namespace Caesar {

class Contour;

class Blob : public TObject {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		Blob();
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~Blob();

	public:
		void SetId(int id){Id=id;}
		void SetName(std::string name){Name=name;}
		bool IsAtEdge(){return HasPixelsAtEdge;}	
		void SetEdgeFlag(bool choice){HasPixelsAtEdge=choice;}

		//Pixel handlers
		//PixelCollection GetPixels(){return m_Pixels;}
		int GetNPixels(){return (int)(m_Pixels.size());}	
		const PixelCollection& GetPixels() const {return m_Pixels;}
		void AddPixel(Pixel* pixel);
		bool HasPixels(){return (GetNPixels()>0);}
		Pixel* GetPixel(int index) {
			if(!HasPixels() || index<0 || index>=GetNPixels() ) return 0;
			return m_Pixels[index];
		}
		
		//Parameters
		int ComputeStats(bool computeRobustStats=true,bool forceRecomputing=false);
		bool HasStats(){return m_HasStats;}
		int ComputeMorphologyParams();
		bool HasParameters(){return m_HasParameters;}

		//Image setters/getters
		void SetImageRange(double xmin,double xmax,double ymin,double ymax){
			m_ImageMinX= xmin;
			m_ImageMaxX= xmax;
			m_ImageMinY= ymin;
			m_ImageMaxY= ymax;
		}
		void GetImageRange(double& xmin,double& xmax,double& ymin,double& ymax){
			xmin= m_ImageMinX;
			xmax= m_ImageMaxX;
			ymin= m_ImageMinY;
			ymax= m_ImageMaxY;
		}

		void SetImageSRange(double Smin,double Smax){m_ImageMinS=Smin; m_ImageMaxS=Smax;}
		void GetImageSRange(double& Smin,double& Smax){Smin=m_ImageMinS; Smax=m_ImageMaxS;}

		void SetImageRMS(double rms){m_ImageRMS=rms;}
		double GetImageRMS(){return m_ImageRMS;}

		void SetImageScurvRange(double Smin,double Smax){m_ImageMinScurv=Smin; m_ImageMaxScurv=Smax;}
		void GetImageScurvRange(double& Smin,double& Smax){Smin=m_ImageMinScurv; Smax=m_ImageMaxScurv;}

		void SetImageSedgeRange(double Smin,double Smax){m_ImageMinSedge=Smin; m_ImageMaxSedge=Smax;}
		void GetImageSedgeRange(double& Smin,double& Smax){Smin=m_ImageMinSedge; Smax=m_ImageMaxSedge;}

		void SetImageSize(long int Nx,long int Ny){m_ImageSizeX=Nx; m_ImageSizeY=Ny;}
		void GetImageSize(long int& sizeX,long int& sizeY){sizeX=m_ImageSizeX;sizeY=m_ImageSizeY;}

		/**
		* \brief Dump region info
		*/
		void Print(){
			cout<<"*** REGION NO. "<<Id<<" ***"<<endl;
			cout<<"N= "<<NPix<<" Smin/Smax="<<m_Smin<<"/"<<m_Smax<<" Xmin/Xmax="<<m_Xmin<<"/"<<m_Xmax<<", Ymin/Ymax="<<m_Ymin<<"/"<<m_Ymax<<endl;
			cout<<"X0="<<X0<<" Y0="<<Y0<<" Mean="<<Mean<<" RMS="<<RMS<<" Median="<<Median<<" MedianRMS="<<MedianRMS<<endl;
			cout<<"****************************"<<endl;
		}
		/**
		* \brief Get image
		*/
		Img* GetImage(Img::ImgType mode);
		/**
		* \brief Return contours
		*/
		std::vector<Contour*> GetContours(){return m_Contours;}
		

	private:
		//Init functions
		void Init();
		void UpdateMoments(Pixel* pixel);
		void ResetMoments();
		void ResetStats();
		void ClearPixels();
		void ClearContours();

	public:
		
		bool HasPixelsAtEdge;

		//Main params
		long int Id;//Blob id
		std::string Name;//Blob name
			
		//Stats params
		long int NPix;//Number of pixels in blob
		double Mean;//mean = M1/N
		double RMS;
		double Skewness;
		double Median;
		double MedianRMS;
		double X0;//X position average
		double Y0;//Y position average
			
		//Curvature Moments
		double Mean_curv;
		double RMS_curv;
		double Median_curv;
		double MedianRMS_curv;

		//2D morphological pars
		std::vector<double> Moments;
		std::vector<double> HuMoments;
		std::vector<double> ZMMoments;	

	protected:
		bool m_HasStats;	
		bool m_HasParameters;	

		//Pixel intensity moments
		double m_M1;//1st moment
		double m_M2;//2nd moment
		double m_M3;//3rd moment
		double m_M4;//4th moment
				
		//Pixel curvature moments
		double m_M1_curv;//1st moment
		double m_M2_curv;//2nd moment

		//Moments accumulator
		double m_S;//sum of pixel signals
		double m_Smax;//max of pixel signals
		double m_Smin;//min of pixel signals
		double m_Sxx;
		double m_Syy;
		double m_Sxy;
		double m_Sx;//Signal-weighted X position average
		double m_Sy;//Signal-weighted Y position average
		int m_PixIdmax;//id of pixel with max signal
		int m_PixIdmin;//id of pixel with min signal
		
		double m_S_curv;//sum of pixel curvature
		double m_S_edge;//sum of edge estimator

		//Image ranges
		long int m_ImageSizeX;
		long int m_ImageSizeY;
		double m_ImageMinX;
		double m_ImageMaxX;
		double m_ImageMinY;
		double m_ImageMaxY;
		double m_ImageMinS;
		double m_ImageMaxS;
		double m_ImageMinScurv;
		double m_ImageMaxScurv;	
		double m_ImageMinSedge;
		double m_ImageMaxSedge;
		double m_ImageRMS;

		double m_Xmin;
		double m_Xmax;
		double m_Ymin;
		double m_Ymax;
		int m_Ix_min;
		int m_Ix_max;
		int m_Iy_min;
		int m_Iy_max;

		//Pixel collection
		PixelCollection m_Pixels;
			
		//Contour collection
		std::vector<Contour*> m_Contours;

	ClassDef(Blob,1)

};//close Blob()



#ifdef __MAKECINT__
#pragma link C++ class Blob+;
#pragma link C++ class vector<Blob>+;
#pragma link C++ class vector<Blob*>+;
#endif

}//close namespace

#endif

