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
* @file FITSReader.h
* @class FITSReader
* @brief FITSReader
*
* Image Reader class for FITS files
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef FITSReader_h
#define FITSReader_h 1

#include <SysUtils.h>

#include <TObject.h>
#include <TFITS.h>
#include <TMath.h>

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


namespace Caesar {

class Img;
class ImgMetaData;

class FITSHeader : public TObject {

	public:
		FITSHeader(){};
		virtual ~FITSHeader(){};
		
	public:
		void Print(){
			cout<<"*** HEADER INFO ***"<<endl;
			cout<<"Image Size: "<<Nx<<"x"<<Ny<<" pixels, nRec="<<nRec<<endl;
			cout<<"Obs Coords: ("<<ObsRA<<","<<ObsDEC<<")"<<endl;
			cout<<"BUnit: "<<BUnit<<endl;
			cout<<"Coords Type: ("<<CoordTypeX<<","<<CoordTypeY<<")"<<endl;
			cout<<"PixelCoordCenter: ("<<Cx<<","<<Cy<<")"<<endl;
			cout<<"CoordCenter: ("<<Xc<<","<<Yc<<")"<<endl;
			cout<<"PixelStep: ("<<dX<<","<<dY<<")"<<endl;
			cout<<"BeamSize: ("<<Bmaj<<","<<Bmin<<","<<Bpa<<")"<<endl;
			cout<<"Rot: ("<<RotX<<","<<RotY<<")"<<endl;
			cout<<"Epoch: "<<Epoch<<endl;
			cout<<"***********************"<<endl;	
		}

	public:
		int Nx;
		int Ny;
		int nRec;
		double ObsRA;
		double ObsDEC;
		std::string BUnit;
		std::string CoordTypeX;
		std::string CoordTypeY;
		int Cx;
		int Cy;
		double Xc;
		double Yc;
		double dX;
		double dY;
		double Bmaj;
		double Bmin;
		double Bpa;
		double Epoch;
		double RotX;
		double RotY;

	ClassDef(FITSHeader,1)

};//close class

#ifdef __MAKECINT__
#pragma link C++ class FITSHeader+;
#endif


class FITSFileInfo : public TObject {

	public:
		FITSFileInfo(){};
		virtual ~FITSFileInfo(){};
		
	public:
		void Print(bool printHeader=true){
			info.Print();	
			if(printHeader)	header.Print();
		}	
		void PrintHeader(){header.Print();}

	public:	
		FITSHeader header;
		Caesar::FileInfo info;

	ClassDef(FITSFileInfo,1)

};

#ifdef __MAKECINT__
#pragma link C++ class FITSFileInfo+;
#endif

class FITSReader : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    FITSReader();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~FITSReader();


	public:
	
		static int Read(std::string filename,Caesar::Img& img,Caesar::FITSFileInfo& info,bool checkFile=true);
		static int ReadTile(std::string filename,Caesar::Img& img,Caesar::FITSFileInfo& info,int xMin,int xMax,int yMin,int yMax,bool checkFile=true);	
		static int ReadTileFast(std::string filename,Caesar::Img& image,Caesar::FITSFileInfo& fits_info,int ix_min,int ix_max,int iy_min,int iy_max,bool checkFile=true);
	private:

		static bool ReadHeader(TFITSHDU* hdu,Caesar::FITSFileInfo& fits_info);
		static TFITSHDU* ReadFile(std::string filename,Caesar::FITSFileInfo& fits_info,bool checkFile=true);
		
	private:

		ClassDef(FITSReader,1)

};

#ifdef __MAKECINT__
#pragma link C++ class FITSReader+; 
#endif

}//close namespace 


#endif


