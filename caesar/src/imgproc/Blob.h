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

#ifdef BUILD_CAESAR_SERVER
	#include <msgpack.hpp>
#endif

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


struct MatchPixelType {
	MatchPixelType(const int& type) : m_type(type) {}
 	bool operator()(const Pixel* obj) const {
  	return obj->type == m_type;
 	}
 	private:
  	const int& m_type;
};


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

		std::vector<int> GetSeedPixelIndexes(){
			std::vector<int> seedPixelIndexes;
			if(m_Pixels.empty()) return seedPixelIndexes;
			std::vector<Pixel*>::iterator it = m_Pixels.begin();
			while ((it = std::find_if(it, m_Pixels.end(), MatchPixelType(Pixel::eSeed))) != m_Pixels.end()) {
    		//int id= (*it)->id;
				int index = std::distance(m_Pixels.begin(), it);
				seedPixelIndexes.push_back(index);
    		it++;
			}//end loop
			return seedPixelIndexes;
		}//close GetSeedPixelIndexes()
		
		//Parameters
		int ComputeStats(bool computeRobustStats=true,bool forceRecomputing=false);
		bool HasStats(){return m_HasStats;}
		int ComputeMorphologyParams();
		bool HasParameters(){return m_HasParameters;}

		double GetM1(){return m_M1;}
		double GetM2(){return m_M2;}
		double GetM3(){return m_M3;}
		double GetM4(){return m_M4;}
		double GetM1Curv(){return m_M1_curv;}
		double GetM2Curv(){return m_M2_curv;}

		double GetS(){return m_S;}
		double GetSmax(){return m_Smax;}
		double GetSmin(){return m_Smin;}
		double GetSxx(){return m_Sxx;}
		double GetSyy(){return m_Syy;}
		double GetSxy(){return m_Sxy;}
		double GetSx(){return m_Sx;}
		double GetSy(){return m_Sy;}
		long int GetSmaxPixId() {return m_PixIdmax;}
		long int GetSminPixId() {return m_PixIdmin;}
		double GetScurv(){return m_S_curv;}
		double GetSedge(){return m_S_edge;}
		
		
		
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

		void GetSourceRange(double& xmin,double& xmax,double& ymin,double& ymax){
			xmin= m_Xmin; 
			xmax= m_Xmax;
			ymin= m_Ymin;
			ymax= m_Ymax;
		}
		void GetSourcePixelRange(long int& ixmin,long int& ixmax,long int& iymin,long int& iymax){
			ixmin= m_Ix_min; 
			ixmax= m_Ix_max;
			iymin= m_Iy_min;
			iymax= m_Iy_max;
		}

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
		long int m_PixIdmax;//id of pixel with max signal
		long int m_PixIdmin;//id of pixel with min signal
		
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
		long int m_Ix_min;
		long int m_Ix_max;
		long int m_Iy_min;
		long int m_Iy_max;

		//Pixel collection
		PixelCollection m_Pixels;
			
		//Contour collection
		std::vector<Contour*> m_Contours;
	
	ClassDef(Blob,1)

	public:
		#ifdef BUILD_CAESAR_SERVER
			MSGPACK_DEFINE(
				HasPixelsAtEdge,Id,Name,
				NPix,Mean,RMS,Skewness,Median,MedianRMS,X0,Y0,
				Mean_curv,RMS_curv,Median_curv,MedianRMS_curv,
				Moments,HuMoments,ZMMoments,
				m_HasStats,m_HasParameters,
				m_M1,m_M2,m_M3,m_M4,
				m_M1_curv,m_M2_curv,
				m_S,m_Smax,m_Smin,m_Sxx,m_Syy,m_Sxy,m_Sx,m_Sy,m_PixIdmax,m_PixIdmin,
				m_S_curv,m_S_edge,
				m_ImageSizeX,m_ImageSizeY,m_ImageMinX,m_ImageMaxX,m_ImageMinY,m_ImageMaxY,m_ImageMinS,m_ImageMaxS,
				m_ImageMinScurv,m_ImageMaxScurv,m_ImageMinSedge,m_ImageMaxSedge,m_ImageRMS,
				m_Xmin,m_Xmax,m_Ymin,m_Ymax,m_Ix_min,m_Ix_max,m_Iy_min,m_Iy_max,
				m_Pixels
			);
		#endif


};//close Blob()



#ifdef __MAKECINT__
#pragma link C++ class Blob+;
#pragma link C++ class vector<Blob>+;
#pragma link C++ class vector<Blob*>+;
#endif

}//close namespace

#endif

