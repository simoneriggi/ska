/**
* @file ImgFITSReader.h
* @class ImgFITSReader
* @brief ImgFITSReader
*
* FITS Image Reader class
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef ImgFITSReader_h
#define ImgFITSReader_h 1

#include "Img.h"

#include <TFITS.h>
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

using namespace std;


class ImgFITSReader {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    ImgFITSReader();
		ImgFITSReader(std::string name);

		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~ImgFITSReader();

		
		struct FITSHeader {
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
		};

		struct FITSFileInfo {
			std::string path;
			std::string filename;
			std::string filename_wext;
			int size;
		};

	public:
	
		bool SetFile(std::string path);	
		void Read(Img&);
		void ReadTile(Img&,int xMin,int xMax,int yMin,int yMax);	
		void ReadHeader();
		void DumpFileInfo();
		void DumpHeaderInfo();

		FITSHeader GetHeaderInfo(){return fHeader;}
		FITSFileInfo GetFileInfo(){return fFileInfo;}

	private:

		
	private:

		std::string fFullFileName;
		bool fIsValidFile;
		FITSFileInfo fFileInfo;		
		
		bool fHasHeader;
		bool fIsValidHeader;
		TFITSHDU* hdu;
		FITSHeader fHeader;

		ClassDef(ImgFITSReader,1)

};

#ifdef __MAKECINT__
#pragma link C++ class ImgFITSReader+; 
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<ImgFITSReader*>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<ImgFITSReader>+;
#endif

#endif


