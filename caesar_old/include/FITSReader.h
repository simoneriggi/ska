/**
* @file FITSReader.h
* @class FITSReader
* @brief FITSReader
*
* FITS Image Reader class
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef FITSReader_h
#define FITSReader_h 1

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


class FITSReader {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    FITSReader();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~FITSReader();

		
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
	
		bool Read(Img& img,std::string filename);
		bool ReadTile(Img& img,std::string filename,int xMin,int xMax,int yMin,int yMax);	
		void DumpFileInfo();
		void DumpHeaderInfo();

		//FITSHeader GetHeaderInfo(){return fHeader;}
		//FITSFileInfo GetFileInfo(){return fFileInfo;}
		FITSHeader* GetHeaderInfo(){return fHeader;}
		FITSFileInfo* GetFileInfo(){return fFileInfo;}

	private:

		void Init();
		bool ReadHeader();
		bool SetFile(std::string path);	
		Img::MetaData* SetMetadata();

	private:

		bool fIsValidFile;
		//FITSFileInfo fFileInfo;		
		FITSFileInfo* fFileInfo;
		
		
		bool fHasHeader;
		bool fIsValidHeader;
		TFITSHDU* hdu;
		//FITSHeader fHeader;
		FITSHeader* fHeader;

		ClassDef(FITSReader,1)

};

#ifdef __MAKECINT__
#pragma link C++ class FITSReader+; 
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<FITSReader*>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<ImgFITSReader>+;
#endif

#endif


