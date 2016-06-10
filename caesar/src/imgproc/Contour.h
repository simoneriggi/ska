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



#include <Logger.h>

#include <TGraph.h>
#include <TPolyLine.h>
#include <TPaveText.h>
#include <TEllipse.h>
#include <TVectorD.h>
#include <TVector2.h>
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
#include <complex>


using namespace std;

#ifdef BUILD_CAESAR_SERVER
	#include <msgpack.hpp>

	//Serialization for TVector2
	namespace msgpack {
	MSGPACK_API_VERSION_NAMESPACE(MSGPACK_DEFAULT_API_NS) {
		namespace adaptor {

			template<>
			struct convert<TVector2> {
    		msgpack::object const& operator()(msgpack::object const& o, TVector2& v) const {
        	if (o.type != msgpack::type::ARRAY) throw msgpack::type_error();
        	if (o.via.array.size != 2) throw msgpack::type_error();
        	v = TVector2(
            o.via.array.ptr[0].as<double>(),
            o.via.array.ptr[1].as<double>());
        	return o;
    		}
			};//close struct

			template<>
			struct pack<TVector2> {
    		template <typename Stream>
    			packer<Stream>& operator()(msgpack::packer<Stream>& o, TVector2 const& v) const {
      	  	// packing member variables as an array.
      	  	o.pack_array(2);
      	  	o.pack(v.X());
      	  	o.pack(v.Y());
      	  	return o;
    			}
			};//close struct

			template <>
			struct object_with_zone<TVector2> {
    		void operator()(msgpack::object::with_zone& o, TVector2 const& v) const {
      		o.type = type::ARRAY;
        	o.via.array.size = 2;
        	o.via.array.ptr = static_cast<msgpack::object*>(
        	o.zone.allocate_align(sizeof(msgpack::object) * o.via.array.size));
        	o.via.array.ptr[0] = msgpack::object(v.X(), o.zone);
      		o.via.array.ptr[1] = msgpack::object(v.Y(), o.zone);
    		}
			};//close struct
		} // namespace adaptor
	} // MSGPACK_API_VERSION_NAMESPACE(MSGPACK_DEFAULT_API_NS)
} // namespace msgpack
#endif


namespace Caesar {

class Contour : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    Contour();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~Contour();

		//typedef std::vector<cv::Point2f> Points;
		typedef std::vector<TVector2> Points;
		

	public:
		/**
		* \brief Get number of points
		*/
		int GetN(){return (int)m_Points.size();}	
		/**
		* \brief Check if contour has points
		*/
		bool HasPoints(){return (m_Points.size()>0);}
		/**
		* \brief Set contour points
		*/
		void SetPoints(Points p){m_Points=p;}
		/**
		* \brief Get contour points
		*/
		Points GetPoints(){return m_Points;}
		/**
		* \brief Get contour point with given index
		*/
		TVector2* GetPoint(int i){
			if(GetN()<=0 || i<0 || i>=GetN()) return 0;
			return (&m_Points[i]);
		}
		/**
		* \brief Add contour points
		*/
		void AddPoint(TVector2 p){
			m_Points.push_back(p);
		}
		/**
		* \brief Reset contour
		*/
		void Reset(){	
			m_Points.clear();
			HasParameters= false;
		}
		
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
		/**
		* \brief Compute shape parameters
		*/
		void ComputeShapeParams();
		/**
		* \brief Dump
		*/
		void Dump(){
			cout<<"== ContourINFO =="<<endl;
			cout<<"C("<<Centroid.X()<<","<<Centroid.Y()<<") Area: "<<Area<<", Perimeter: "<<Perymeter<<" BoundingBox: ("<<BoundingBoxMin<<","<<BoundingBoxMaj<<","<<BoundingBoxAngle<<")"<<endl;
			cout<<"Elong: "<<Elongation<<" Rectangularity: "<<Rectangularity<<" Roundness="<<Roundness<<" Eccentricity="<<Eccentricity<<endl;
			cout<<"CircularityRatio: "<<CircularityRatio<<" EllipseAreaRatio="<<EllipseAreaRatio<<" Ellipse("<<EllipseCenter.X()<<","<<EllipseCenter.Y()<<","<<EllipseMinAxis<<","<<EllipseMajAxis<<","<<EllipseRotAngle<<")"<<endl;
			cout<<"HuMoments=(";
			for(int k=0;k<7;k++) cout<<HuMoments[k]<<",";
			cout<<")"<<endl;
			cout<<"================="<<endl;
		}//close Dump()

		/**
		* \brief Log 
		*/
		void Log(std::string level="INFO"){
			LOG(level,GetPrintable());
		}

		/**
		* \brief GetPrintable
		*/
		std::string GetPrintable(){
			std::stringstream ss;
			ss<<"ContourInfo: ";
			ss<<"C("<<Centroid.X()<<","<<Centroid.Y()<<") Area: "<<Area<<", Perimeter: "<<Perymeter<<" BoundingBox: ("<<BoundingBoxMin<<","<<BoundingBoxMaj<<","<<BoundingBoxAngle<<"), ";
			ss<<"Elong: "<<Elongation<<" Rectangularity: "<<Rectangularity<<" Roundness="<<Roundness<<" Eccentricity="<<Eccentricity<<", ";
			ss<<"CircularityRatio: "<<CircularityRatio<<" EllipseAreaRatio="<<EllipseAreaRatio<<" Ellipse("<<EllipseCenter.X()<<","<<EllipseCenter.Y()<<","<<EllipseMinAxis<<","<<EllipseMajAxis<<","<<EllipseRotAngle<<"), ";
			ss<<"HuMoments=(";
			for(int k=0;k<7;k++) ss<<HuMoments[k]<<",";
			ss<<")";
			return ss.str();
		}

	private:
		void ComputeArea();
		void ComputePerymeter();
		void ComputeCircularityRatio();
		void ComputeBoundingBox();
		void ComputeElongation();
		void ComputeRectangularity();
		void ComputeRoundness();
		void ComputeMoments();
		void ComputeEccentricity();
		void ComputeFittedEllipse();

		std::vector< std::complex<double> > GetComplexPointRepresentation(bool translateToCentroid=false){
			std::vector< std::complex<double> > U;
			if(translateToCentroid) {
				for(unsigned int i=0;i<m_Points.size();i++)
					U.push_back( std::complex<double>(m_Points[i].X()-Centroid.X(),m_Points[i].Y()-Centroid.Y()) );
			}
			else {
				for(unsigned int i=0;i<m_Points.size();i++) 
					U.push_back( std::complex<double>(m_Points[i].X(),m_Points[i].Y()) );	
			}
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
		TVector2 BoundingBoxCenter;
		double BoundingBoxMaj;
		double BoundingBoxMin;
		double BoundingBoxAngle;
		double Elongation;
		double Rectangularity;
		double Roundness;
		double Eccentricity;
		double TiltAngle;
		bool HasEllipseFit;
		TVector2 EllipseCenter;
		double EllipseMajAxis;
		double EllipseMinAxis;
		double EllipseRotAngle;
		double EllipseFitRedChi2;
		double EllipseAreaRatio;

		// spatial moments: m00, m10, m01, m20, m11, m02, m30, m21, m12, m03
		double m00;
		double m10;
		double m01;
		double m20;
		double m11;
		double m02;
		double m30;
		double m21;
		double m12;
		double m03;

		double mu20;
		double mu11;
		double mu02; 
		double mu30; 
		double mu21; 
		double mu12;
		double mu03;

 	 	double nu20; 
		double nu11; 
		double nu02; 
		double nu30; 
		double nu21;
		double nu12;
		double nu03;

		//double HuMoments[7];
		std::vector<double> HuMoments;
		//TVector2 BoundingBoxVertex[4];	
		std::vector<TVector2> BoundingBoxVertex;
		TVector2 Centroid;
		
		//std::vector< std::complex<double> > FDs;
		std::vector<double> RealFDs;
		std::vector<double> ImagFDs;
		std::vector<double> ModFDs;//module of complex Fourier descriptors
		std::vector<double> BendingEnergies;
		std::vector<double> CentroidDistanceModFDs;//module of complex Fourier descriptors

	private:
		Points m_Points;	
		
	ClassDef(Contour,1)

	public:	
		#ifdef BUILD_CAESAR_SERVER
			MSGPACK_DEFINE(
				HasParameters,Area,Perymeter,IsConvexContour,CircularityRatio,
				BoundingBoxMaj,BoundingBoxMin,BoundingBoxAngle,
				Elongation,Rectangularity,Roundness,Eccentricity,TiltAngle,
				HasEllipseFit,EllipseMajAxis,EllipseMinAxis,EllipseRotAngle,EllipseFitRedChi2,EllipseAreaRatio,
				BoundingBoxCenter,EllipseCenter,
				m00,m10,m01,m20,m11,m02,m30,m21,m12,m03,mu20,mu11,mu02,mu30,mu21,mu12,mu03,
		 		nu20,nu11,nu02,nu30,nu21,nu12,nu03,
				HuMoments,BoundingBoxVertex,Centroid,
				RealFDs,ImagFDs,ModFDs,BendingEnergies,CentroidDistanceModFDs,
				m_Points
			)
		#endif

};//close class

#ifdef __MAKECINT__
#pragma link C++ class Contour+;
#pragma link C++ class vector<Contour>+;
#pragma link C++ class vector<Contour*>+;
#endif

}//close namespace 

#endif


