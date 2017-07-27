/**
* @file Contour.h
* @class Contour
* @brief Contour
*
* Class representing image contour with methods for morphological parameter extraction 
* @author S. Riggi
* @date 11/07/2015
*/

#ifndef Contour_h
#define Contour_h 1

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
#include <TGraph.h>
#include <TPolyLine.h>
#include <TPaveText.h>
#include <TEllipse.h>
#include <TVectorD.h>
#include <TMath.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <complex>

#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>


using namespace std;
//using namespace cv;

class Contour {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    Contour();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~Contour();

		typedef std::vector<cv::Point2f> Points;
		

	public:
		/**
		* \brief Get number of points
		*/
		int GetN(){return (int)fPoints.size();}
		/**
		* \brief Get contour points
		*/
		void SetPoints(Points p){fPoints=p;}
		/**
		* \brief Get contour points
		*/
		Points GetPoints(){return fPoints;}
		/**
		* \brief Get contour point with given index
		*/
		cv::Point2f* GetPoint(int i){
			if(GetN()<=0 || i<0 || i>=GetN()) return 0;
			return (&fPoints[i]);
		}
		/**
		* \brief Add contour points
		*/
		void AddPoint(cv::Point2f p){	
			//cout<<"p("<<p.x<<","<<p.y<<")"<<endl;
			fPoints.push_back(p);
		}
		/**
		* \brief Reset contour
		*/
		void Reset(){	
			fPoints.clear();
			HasParameters= false;
		}
		/**
		* \brief Check if contour has points
		*/
		bool HasPoints(){return (fPoints.size()>0);}
		/**
		* \brief Return a graph object with contour points
		*/
		TGraph* GetGraph();
		/**
		* \brief Return a polyline object with bounding box
		*/	
		TPolyLine* GetBoundingBoxLine();
		/**
		* \brief Return an ellipse object fitted to contour
		*/
		TEllipse* GetFittedEllipse();
		/**
		* \brief Return a info box with parameter values
		*/
		TPaveText* GetParamInfoBox();
		/**
		* \brief Compute contour parameters
		*/
		int ComputeParameters();

		void Dump(){
			cout<<"== ContourINFO =="<<endl;
			cout<<"C("<<Centroid.x<<","<<Centroid.y<<") Area: "<<Area<<", Perimeter: "<<Perymeter<<" BoundingBox: ("<<BoundingBoxMin<<","<<BoundingBoxMaj<<","<<BoundingBoxAngle<<")"<<endl;
			cout<<"Elong: "<<Elongation<<" Rectangularity: "<<Rectangularity<<" Roundness="<<Roundness<<" Eccentricity="<<Eccentricity<<endl;
			cout<<"CircularityRatio: "<<CircularityRatio<<" EllipseAreaRatio="<<EllipseAreaRatio<<" Ellipse("<<EllipseCenter.x<<","<<EllipseCenter.y<<","<<EllipseMinAxis<<","<<EllipseMajAxis<<","<<EllipseRotAngle<<")"<<endl;
			cout<<"HuMoments=(";
			for(int k=0;k<7;k++) cout<<HuMoments[k]<<",";
			cout<<")"<<endl;
			cout<<"================="<<endl;
		}//close Dump()

	private:
		void ComputeArea();
		void ComputePerymeter();
		void ComputeCircularityRatio();
		void ComputeBoundingBox();
		void ComputeElongation();
		void ComputeRectangularity();
		void ComputeRoundness();
		void ComputeMoments();
		void ComputeHuMoments();
		void ComputeEccentricity();
		void ComputeFittedEllipse();

		std::vector< std::complex<double> > GetComplexPointRepresentation(bool translateToCentroid=false){
			std::vector< std::complex<double> > U;
			if(translateToCentroid) for(unsigned int i=0;i<fPoints.size();i++) U.push_back( std::complex<double>(fPoints[i].x-Centroid.x,fPoints[i].y-Centroid.y) );
			else for(unsigned int i=0;i<fPoints.size();i++) U.push_back( std::complex<double>(fPoints[i].x,fPoints[i].y) );	
			return U;
		}

		void ComputeFourierDescriptors();
		void ComputeCentroidDistanceFD();
		void ComputeBendingEnergy();

		TVectorD EllipseFitter(TGraph*);
		TVectorD ConicToParametric(const TVectorD &conic);
		double EllipseFcn(double x, double y, TVectorD params) {
  		double v = 9999.9;
			double x0= params[0];
			double y0= params[1];	
			double a= params[2];
			double b= params[3];
			double theta= params[4];
  		if ((a == 0.0) || (b == 0.0)) return v; // just a precaution
  		// shift the center
  		x-= x0;
  		y-= y0;
  		// un-rotate the axes
  		theta*= TMath::Pi()/180.0; // degrees -> radians
  		v = x;
  		x = x * std::cos(theta) + y * std::sin(theta);
  		y = y * std::cos(theta) - v * std::sin(theta);
  		// "scale" axes
  		x/= a;
  		y/= b;
  		// calculate the "normalized distance"
  		v = x * x + y * y;
  		v = sqrt(v);
  		return v;
		}
		double EllipseFitChi2(TGraph* contourGraph,TVectorD params)	{
  		if (!contourGraph) return 0; // just a precaution
			double v = 0.;
  		for (int i=0; i<contourGraph->GetN(); i++) {
				double x, y;
				contourGraph->GetPoint(i,x,y);
				double r = EllipseFcn(x,y,params);
    		r-= 1.0; // ellipse's "radius" in "normalized coordinates" is always 1
    		v += r*r;
  		}//end loop points
  		return v;
		}//close EllipseFitChi2()

		

	public:	
		bool HasParameters;
		double Area;
		double Perymeter;
		bool IsConvexContour;
		double CircularityRatio;
		cv::Point2f BoundingBoxCenter;
		double BoundingBoxMaj;
		double BoundingBoxMin;
		double BoundingBoxAngle;
		double Elongation;
		double Rectangularity;
		double Roundness;
		double Eccentricity;
		double TiltAngle;
		bool HasEllipseFit;
		cv::Point2f EllipseCenter;
		double EllipseMajAxis;
		double EllipseMinAxis;
		double EllipseRotAngle;
		double EllipseFitRedChi2;
		double EllipseAreaRatio;

		cv::Moments Moments;
		double HuMoments[7];
		cv::Point2f BoundingBoxVertex[4];	
		cv::Point2f Centroid;
		
		std::vector< std::complex<double> > FDs;
		std::vector<double> ModFDs;//module of complex Fourier descriptors
		std::vector<double> BendingEnergies;
		std::vector<double> CentroidDistanceModFDs;//module of complex Fourier descriptors

	private:
		Points fPoints;	

};//close class

#ifdef __MAKECINT__
#pragma link C++ class Contour+;
#pragma link C++ class vector<Contour>+;
#pragma link C++ class vector<Contour*>+;
#pragma link C++ class cv::Point2f+;
#pragma link C++ class vector<cv::Point2f>+;
#pragma link C++ class cv::Point+;
#pragma link C++ class vector<cv::Point>+;
#pragma link C++ class cv::Moments;
#endif

#endif


