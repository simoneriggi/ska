/**
* @file ImgFITSReader.cc
* @class ImgFITSReader
* @brief ImgFITSReader
*
* FITS Image Reader class
* @author S. Riggi
* @date 20/01/2015
*/

#include <ImgFITSReader.h>
#include "Img.h"

#include <boost/filesystem.hpp>

#include <TFITS.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TCut.h>
#include <TEventList.h>
#include <TMath.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TColor.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include <Math/WrappedTF1.h>
#include <Math/GSLIntegrator.h>
#include <Math/GSLMinimizer.h>
#include <Math/Functor.h>
#include <Math/WrappedFunction.h>
#include <Math/WrappedParamFunction.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>
#include <Math/SpecFunc.h>
#include <Math/DistFunc.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
using namespace std;


ClassImp(ImgFITSReader)


ImgFITSReader::ImgFITSReader() {

	fFileInfo.path= "";
	fIsValidFile= false;
	fHasHeader= false;
	fIsValidHeader= false;
	hdu= NULL;
	
}//close costructor

ImgFITSReader::ImgFITSReader(std::string filename) {

	fFileInfo.path= filename;
	fIsValidFile= false;
	fHasHeader= false;
	fIsValidHeader= false;
	hdu= NULL;

	fIsValidFile= SetFile(filename);
	if(!fIsValidFile) cerr<<"ImgFITSReader::ImgFITSReader(): ERROR: Given file is invalid!"<<endl;
	
}//close costructor


ImgFITSReader::~ImgFITSReader() {

}//close destructor


bool ImgFITSReader::SetFile(std::string path){

	fFileInfo.path= path;
	fIsValidFile= false;

	//## Check filename
	if(path==""){
		cerr<<"ImgReader::SetImageFile(): WARNING: Empty filename given!"<<endl;	
		return false;
	} 

	//## Check file given
	boost::filesystem::path filePath(path.c_str());

  try {
		//Check if file exists on filesystem
    if (boost::filesystem::exists(filePath)) {
      if (boost::filesystem::is_regular_file(filePath)){
				
				//Get filename and extension
				if(filePath.has_filename()){
					fFileInfo.filename= filePath.filename().string();
					fFileInfo.filename_wext= filePath.stem().string();
					fFileInfo.size= boost::filesystem::file_size(filePath);
        	
					if(filePath.has_extension()) {
						std::string fileExtension= filePath.extension().string();	
						fIsValidFile= true;
						
						if(fileExtension!=".fits"){
							cout<<"ImgReader::ImgReader(): ERROR: Unknown or unsupported file extension detected ("<<fileExtension<<")..."<<endl;
							return false;
						}
	
						DumpFileInfo();

					}//close if
					else{
						cerr << "ImgReader::SetImageFile(): ERROR: Given file without extension!"<<endl;
						return false;
					}
				}//close if has filename
				else{
					cerr<<"ImgReader::ImgReader(): ERROR: Given path has not a filename!"<<endl;
					return false;	
				}
			}//close if
      else if (boost::filesystem::is_directory(filePath)){
        cerr << "ImgReader::SetImageFile(): ERROR: Given file ("<<filePath<<") is a directory!"<<endl;
				return false;
      }
      else{
        cerr << "ImgReader::SetImageFile(): ERROR: Given file ("<<filePath<<") exists, but is neither a regular file nor a directory!"<<endl;
				return false;
			}
    }//close if exists
    else{
      cerr << "ImgReader::SetImageFile(): ERROR: Given file ("<<filePath<<") does not exist!"<<endl;
			return false;
		}
  }//close try block

  catch (const boost::filesystem::filesystem_error& ex) {
    cout << ex.what() << '\n';
		return false;
  }

	return true;

}//close ImgReader::SetImageFile()


void ImgFITSReader::Read(Img& image){

	//## Check file
	if(!fIsValidFile){
		cerr<<"ImgReader::Read(): WARNING: No file given or invalid file, no read operations will be performed!"<<endl;
		throw std::runtime_error("ImgFITSReader::Read(): ERROR: No file given or invalid file, no read operations will be performed!");
	}
	
	//## Read header
	try{
		if(!fHasHeader) ReadHeader();
	}
	catch(std::exception const & e) {
		std::stringstream errMsg;
   	errMsg << "ImgReader::Read(): ERROR: Image header read failed with status "<<e.what();
		throw std::runtime_error(errMsg.str().c_str());
	}

	//Img::MetaData metadata;
	Img::MetaData metadata;
	metadata.Nx= fHeader.Nx;
	metadata.Ny= fHeader.Ny;
	metadata.Cx= fHeader.Cx;
	metadata.Cy= fHeader.Cy;
	metadata.Xc= fHeader.Xc;
	metadata.Yc= fHeader.Yc;
	metadata.CoordTypeX= fHeader.CoordTypeX;
	metadata.CoordTypeY= fHeader.CoordTypeY;
	metadata.BUnit= fHeader.BUnit;
	metadata.Bmaj= fHeader.Bmaj;
	metadata.Bmin= fHeader.Bmin;
	metadata.Bpa= fHeader.Bpa;
	metadata.dX= fHeader.dX;
	metadata.dY= fHeader.dY;
	metadata.RotX= fHeader.RotX;
	metadata.RotY= fHeader.RotY;
	metadata.Epoch= fHeader.Epoch;

	//## Set image properties	
	image.Reset();
	image.SetBins(fHeader.Nx,-0.5,fHeader.Nx-0.5,fHeader.Ny,-0.5,fHeader.Ny-0.5);
	image.SetMetaData(metadata);
	//image.SetWCS();

	//## Read image
	for(int j=0;j<fHeader.Ny;j++){
		TVectorD* thisPixelRow= hdu->GetArrayRow(j);
		//cout<<"INFO: RowId="<<i<<" N="<<thisPixelRow->GetNoElements()<<" minY="<<binY_min<<" maxY="<<binY_max<<endl;
		for(int i=0;i<fHeader.Nx;i++){
			double pixValue= ((*thisPixelRow))(i);
			if(TMath::IsNaN(pixValue) || fabs(pixValue)==TMath::Infinity()) {
				//cerr<<"ImgFITSReader::Read(): WARNING: NaN or Inf pixel value..."<<endl;
			}
			else{
				//if(pixValue<=0) cerr<<"ImgFITSReader::Read(): WARNING: Negative or null pixel value..."<<endl;
				//cout<<"thisCell("<<i<<","<<j<<")= "<<pixValue<<endl;
				image.FillPixel(i,j,pixValue);
			}
		}//end loop columns
		if(thisPixelRow) thisPixelRow->Delete();
	}//end loop row

}//close ImgFITSReader::Read()


void ImgFITSReader::ReadTile(Img& image,int ix_min,int ix_max,int iy_min,int iy_max){

	//## Check file
	if(!fIsValidFile){
		cerr<<"ImgReader::Read(): WARNING: No file given or invalid file, no read operations will be performed!"<<endl;
		throw std::runtime_error("ImgFITSReader::Read(): ERROR: No file given or invalid file, no read operations will be performed!");
	}

	//## Read header
	try{
		if(!fHasHeader) ReadHeader();
	}
	catch(std::exception const & e) {
		std::stringstream errMsg;
   	errMsg << "ImgReader::Read(): ERROR: Image header read failed with status "<<e.what();
		throw std::runtime_error(errMsg.str().c_str());
	}

	//## Check tile sizes
	if(ix_min<0 || ix_min>=fHeader.Nx || ix_min>=ix_max){
		std::stringstream errMsg;
		errMsg<<"ImgReader::ReadTile(): ERROR: Invalid min tile X given ("<<ix_min<<")";
		throw std::invalid_argument(errMsg.str().c_str());
	}
	if(ix_max<0 || ix_max>=fHeader.Nx || ix_max<=ix_min){
		std::stringstream errMsg;
		errMsg<<"ImgReader::ReadTile(): ERROR: Invalid max tile Y given ("<<ix_max<<")";
		throw std::invalid_argument(errMsg.str().c_str());
	}
	if(iy_min<0 || iy_min>=fHeader.Ny || iy_min>=iy_max){
		std::stringstream errMsg;
		errMsg<<"ImgReader::ReadTile(): ERROR: Invalid min tile X given ("<<iy_min<<")";
		throw std::invalid_argument(errMsg.str().c_str());
	}
	if(iy_max<0 || iy_max>=fHeader.Ny || iy_max<=iy_min){
		std::stringstream errMsg;
		errMsg<<"ImgReader::ReadTile(): ERROR: Invalid max tile Y given ("<<iy_max<<")";
		throw std::invalid_argument(errMsg.str().c_str());
	}

	int TileSizeX= ix_max-ix_min+1;
	int TileSizeY= iy_max-iy_min+1;
	//cout<<"TileSizeX="<<TileSizeX<<" binX_min="<<binX_min<<" binX_max="<<binX_max<<endl;
	//cout<<"TileSizeY="<<TileSizeX<<" binY_min="<<binY_min<<" binY_max="<<binY_max<<endl;
		
	//Img::MetaData metadata;
	Img::MetaData metadata;
	metadata.Nx= fHeader.Nx;
	metadata.Ny= fHeader.Ny;
	metadata.Cx= fHeader.Cx;
	metadata.Cy= fHeader.Cy;
	metadata.Xc= fHeader.Xc;
	metadata.Yc= fHeader.Yc;
	metadata.CoordTypeX= fHeader.CoordTypeX;
	metadata.CoordTypeY= fHeader.CoordTypeY;
	metadata.BUnit= fHeader.BUnit;
	metadata.Bmaj= fHeader.Bmaj;
	metadata.Bmin= fHeader.Bmin;
	metadata.Bpa= fHeader.Bpa;
	metadata.dX= fHeader.dX;
	metadata.dY= fHeader.dY;
	metadata.RotX= fHeader.RotX;
	metadata.RotY= fHeader.RotY;
	metadata.Epoch= fHeader.Epoch;

	//## Set image properties
	cout<<"ImgReader::ReadTile(): INFO: Set image properties..."<<endl;
	image.Reset();
	image.SetBins(TileSizeX,ix_min-0.5,ix_max+0.5,TileSizeY,iy_min-0.5,iy_max+0.5);	
	image.SetMetaData(metadata);
	//image.SetWCS();

	//## Read tile
	cout<<"ImgReader::ReadTile(): INFO: Starting reading tile..."<<endl;
	int nFilledPixels= 0;
	int nZeroPixels= 0;

	for(int j=iy_min;j<=iy_max;j++){
		TVectorD* thisPixelRow= hdu->GetArrayRow(j);
		for(int i=ix_min;i<=ix_max;i++){		
			double pixValue= ((*thisPixelRow))(i);
			if(TMath::IsNaN(pixValue) || fabs(pixValue)==TMath::Infinity()) {
				//cerr<<"ImgFITSReader::ReadTile(): WARNING: NaN or Inf pixel value..."<<endl;
			}
			else if(pixValue==0) nZeroPixels++;
			else{
				image.FillPixel(i,j,pixValue);
				nFilledPixels++;
			}
		}//end loop columns
		if(thisPixelRow) thisPixelRow->Delete();
	}//end loop row

	cout<<"ImgReader::ReadTile(): INFO: nFilledPixels="<<nFilledPixels<<" nZeroPixels="<<nZeroPixels<<endl;

	
}//close ImgFITSReader::ReadTile()




void ImgFITSReader::ReadHeader(){

	//## Check file validity
	if(!fIsValidFile){
		throw std::runtime_error("ImgFITSReader::ReadHeader(): ERROR: File is invalid, no header reading will be performed!");
	}

	//## Open primary HDU from file if not already present
	if(hdu) {
		fHasHeader= false;
		delete hdu;
	}
	
	hdu = new TFITSHDU(fFileInfo.path.c_str());
  if (hdu == 0) throw std::runtime_error("ImgReader::ReadHeader(): ERROR: Cannot access the HDU!");
  	
	cout<<"ImgReader::ReadHeader(): INFO: Header read with success..."<<endl;
	 

	fIsValidHeader= true;

	//##### GET STANDARD & MANDATORY KEYWORDS (if not existing set an invalid header...)
	// Get image size field
	int Nx= 0;
	int Ny= 0;
	if(hdu->GetKeywordValue("NAXIS1")=="" || hdu->GetKeywordValue("NAXIS2")=="") {
		fIsValidHeader= false;
		throw std::runtime_error("ImgReader::ReadHeader(): ERROR: Invalid header detected (no NAXIS keywords found!)");
	}
	else{	
		Nx= hdu->GetKeywordValue("NAXIS1").Atoi();//image size X
		Ny= hdu->GetKeywordValue("NAXIS2").Atoi();//image size Y
	}

	// Get pixel content unit
	std::string BUnit= "";
	if(hdu->GetKeywordValue("BUNIT")=="") {
		fIsValidHeader= false;
		//throw std::runtime_error("ImgReader::ReadHeader(): ERROR: Invalid header detected (no BUNIT keyword found!)");
	}
	else{
		BUnit= std::string(hdu->GetKeywordValue("BUNIT").Data());	
		if(BUnit[0]=='\'') BUnit.erase(0,1);
		if(BUnit[BUnit.length()-1]=='\'') BUnit.erase(BUnit.length()-1,1);		
	}

	//Get Coordinate type and projection
	std::string CoordTypeX= "";
	std::string CoordTypeY= "";
	if(hdu->GetKeywordValue("CTYPE1")=="" || hdu->GetKeywordValue("CTYPE2")=="") {
		fIsValidHeader= false;
		//throw std::runtime_error("ImgReader::ReadHeader(): ERROR: Invalid header detected (no CTYPE keyword found!)");
	}
	else{
		CoordTypeX= std::string(hdu->GetKeywordValue("CTYPE1").Data());
		CoordTypeY= std::string(hdu->GetKeywordValue("CTYPE2").Data());

		if(CoordTypeX[0]=='\'') CoordTypeX.erase(0,1);
		if(CoordTypeX[CoordTypeX.length()-1]=='\'') CoordTypeX.erase(CoordTypeX.length()-1,1);
		if(CoordTypeY[0]=='\'') CoordTypeY.erase(0,1);
		if(CoordTypeY[CoordTypeY.length()-1]=='\'') CoordTypeY.erase(CoordTypeY.length()-1,1);
	}

	//Get reference pixel id
	double CenterPixIdX= 0;
	double CenterPixIdY= 0;
	if(hdu->GetKeywordValue("CRPIX1")=="" || hdu->GetKeywordValue("CRPIX2")==""){
		fIsValidHeader= false;
		//throw std::runtime_error("ImgReader::ReadHeader(): ERROR: Invalid header detected (no CRPIX keyword found!)");
	}
	else{
		CenterPixIdX= hdu->GetKeywordValue("CRPIX1").Atof();
		CenterPixIdY= hdu->GetKeywordValue("CRPIX2").Atof();
	}
	
	//Get reference pixel coordinates (RA/DEC)
	double CenterPixX= 0;
	double CenterPixY= 0;
	if(hdu->GetKeywordValue("CRVAL1")=="" || hdu->GetKeywordValue("CRVAL2")==""){
		fIsValidHeader= false;
		//throw std::runtime_error("ImgReader::ReadHeader(): ERROR: Invalid header detected (no CRVAL keyword found!)");
	}
	else{
		CenterPixX= hdu->GetKeywordValue("CRVAL1").Atof();
		CenterPixY= hdu->GetKeywordValue("CRVAL2").Atof();
	}
	
	double PixStepX= 0;
	double PixStepY= 0;
	if(hdu->GetKeywordValue("CDELT1")=="" || hdu->GetKeywordValue("CDELT2")==""){
		fIsValidHeader= false;
		//throw std::runtime_error("ImgReader::ReadHeader(): ERROR: Invalid header detected (no CDELT keyword found!)");
	}
	else{
		PixStepX= hdu->GetKeywordValue("CDELT1").Atof();
		PixStepY= hdu->GetKeywordValue("CDELT2").Atof();
	}

	//Beam information
	double Bmaj= 0;
	double Bmin= 0;
	if(hdu->GetKeywordValue("BMAJ")=="" || hdu->GetKeywordValue("BMIN")==""){
		fIsValidHeader= false;
		//throw std::runtime_error("ImgReader::ReadHeader(): ERROR: Invalid header detected (no BMAJ/BMIN/BPA keywords found!)");
	}
	else{
		Bmaj= hdu->GetKeywordValue("BMAJ").Atof();
		Bmin= hdu->GetKeywordValue("BMIN").Atof();
	}

	double Bpa= 0;
	if(hdu->GetKeywordValue("BPA")==""){
		cerr<<"ImgReader::ReadHeader(): WARN: Missing BPA keyword in header...setting it to 0!"<<endl;
		Bpa= 0;	
	}
	else{
		Bpa= hdu->GetKeywordValue("BPA").Atof();
	}
	
	
	//########  GET NON STANDARD/NON-MANDATORY KEYWORDS
	int nRec=	hdu->GetRecordNumber();

	// Get rotation matrix
	double RotX= 0;
	double RotY= 0;
	if(hdu->GetKeywordValue("CROTA1")=="" || hdu->GetKeywordValue("CROTA2")==""){
		cerr<<"ImgReader::ReadHeader(): WARNING: Cannot find CROTA keyword in FITS header, setting rotation to 0!"<<endl;
	}
	else{
		RotX= hdu->GetKeywordValue("CROTA1").Atof();
		RotY= hdu->GetKeywordValue("CROTA2").Atof();
	}

	// Get Observation Location
	double ObsRA= -999;
	double ObsDEC= -999;
	if(hdu->GetKeywordValue("OBSRA")=="" || hdu->GetKeywordValue("OBSDEC")==""){
		cerr<<"ImgReader::ReadHeader(): WARNING: Cannot find OBSRA/OBSDEC keywords in FITS header, setting them to fake values!"<<endl;
	}
	else{
		ObsRA= hdu->GetKeywordValue("OBSRA").Atof();
		ObsDEC= hdu->GetKeywordValue("OBSDEC").Atof();
	}

	// Get epoch information
	double Epoch= 2000;
	if(hdu->GetKeywordValue("EPOCH")==""){
		cerr<<"ImgReader::ReadHeader(): WARNING: Cannot find epoch keywords in FITS header, setting it to 2000!"<<endl;
	}
	else{
		Epoch= hdu->GetKeywordValue("EPOCH").Atof();
	}

	fHeader.Nx= Nx;
	fHeader.Ny= Ny;
	fHeader.nRec= nRec;
	fHeader.ObsRA= ObsRA;
	fHeader.ObsDEC= ObsDEC;
	fHeader.BUnit= BUnit;
	fHeader.CoordTypeX= CoordTypeX;
	fHeader.CoordTypeY= CoordTypeY;
	fHeader.Cx= CenterPixIdX;
	fHeader.Cy= CenterPixIdY;
	fHeader.Xc= CenterPixX;
	fHeader.Yc= CenterPixY;
	fHeader.dX= PixStepX;
	fHeader.dY= PixStepY;
	fHeader.Bmaj= Bmaj;
	fHeader.Bmin= Bmin;
	fHeader.Bpa= Bpa;
	fHeader.RotX= RotX;
	fHeader.RotY= RotY;
	fHeader.Epoch= Epoch;
	fHasHeader= true;
	
	DumpHeaderInfo();
	
}//close ImgFITSReader::ReadHeader()


void ImgFITSReader::DumpHeaderInfo(){

	// Dump the HDUs within the FITS file and also their metadata
  //if(hdu) hdu->Print("F+");

	cout<<"*** IMAGE META-DATA ***"<<endl;
	cout<<"Image Size: "<<fHeader.Nx<<"x"<<fHeader.Ny<<" pixels, nRec="<<fHeader.nRec<<endl;
	cout<<"Obs Coords: ("<<fHeader.ObsRA<<","<<fHeader.ObsDEC<<")"<<endl;
	cout<<"BUnit: "<<fHeader.BUnit<<endl;
	cout<<"Coords Type: ("<<fHeader.CoordTypeX<<","<<fHeader.CoordTypeY<<")"<<endl;
	cout<<"PixelCoordCenter: ("<<fHeader.Cx<<","<<fHeader.Cy<<")"<<endl;
	cout<<"CoordCenter: ("<<fHeader.Xc<<","<<fHeader.Yc<<")"<<endl;
	cout<<"PixelStep: ("<<fHeader.dX<<","<<fHeader.dY<<")"<<endl;
	cout<<"BeamSize: ("<<fHeader.Bmaj<<","<<fHeader.Bmin<<","<<fHeader.Bpa<<")"<<endl;
	cout<<"Rot: ("<<fHeader.RotX<<","<<fHeader.RotY<<")"<<endl;
	cout<<"Epoch: "<<fHeader.Epoch<<endl;
	cout<<"***********************"<<endl;	

}//close ImgFITSReader::DumpHeaderInfo()


void ImgFITSReader::DumpFileInfo(){			
	cout <<"*** IMAGE FILE INFO ***"<<endl;
	cout << "File Path: "<<fFileInfo.path<<endl;
	cout << "File Name: "<<fFileInfo.filename<<endl;
	cout << "File WExt: "<<fFileInfo.filename_wext<<endl;
	cout << "File size (kB): "<<fFileInfo.size/1024.<<endl;
	cout <<"***********************"<<endl;
}	


