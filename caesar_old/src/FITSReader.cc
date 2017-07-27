/**
* @file FITSReader.cc
* @class FITSReader
* @brief FITSReader
*
* FITS Image Reader class
* @author S. Riggi
* @date 20/01/2015
*/

#include <FITSReader.h>
//#include "Img.h"

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


ClassImp(FITSReader)


FITSReader::FITSReader() {

	fFileInfo= 0;
	hdu= 0;
	fHeader= 0;
	fIsValidFile= false;
	fHasHeader= false;
	fIsValidHeader= false;
	
}//close costructor


FITSReader::~FITSReader() {

	if(fHeader) {
		delete fHeader;
		fHeader= 0;
	}
	if(fFileInfo){
		delete fFileInfo;
		fFileInfo= 0;
	}
	if(hdu){
		delete hdu;
		hdu= 0;
	}
	
}//close destructor

void FITSReader::Init(){

	fIsValidFile= false;
	fHasHeader= false;
	fIsValidHeader= false;
	hdu= NULL;
	fFileInfo= 0;
	fHeader= 0;

}//close Init()

Img::MetaData* FITSReader::SetMetadata(){

	if(!fHeader){
		cerr<<"FITSReader::SetMetaData(): ERROR: Null ptr to FITS header, cannot set metadata!"<<endl;
		return 0;	
	}
	
	Img::MetaData* metadata= new Img::MetaData;
	metadata->Nx= fHeader->Nx;
	metadata->Ny= fHeader->Ny;
	metadata->Cx= fHeader->Cx;
	metadata->Cy= fHeader->Cy;
	metadata->Xc= fHeader->Xc;
	metadata->Yc= fHeader->Yc;
	metadata->CoordTypeX= fHeader->CoordTypeX;
	metadata->CoordTypeY= fHeader->CoordTypeY;
	metadata->BUnit= fHeader->BUnit;
	metadata->Bmaj= fHeader->Bmaj;
	metadata->Bmin= fHeader->Bmin;
	metadata->Bpa= fHeader->Bpa;
	metadata->dX= fHeader->dX;
	metadata->dY= fHeader->dY;
	metadata->RotX= fHeader->RotX;
	metadata->RotY= fHeader->RotY;
	metadata->Epoch= fHeader->Epoch;

	return metadata;

}//close SetMetaData()


bool FITSReader::Read(Img& image,std::string filename){

	//## Init data
	Init();

	//## Check file
	if(!SetFile(filename)){
		cerr<<"FITSReader::Read(): ERROR: No file given or invalid file, nothing will be read!"<<endl;
		return false;
	}
	
	
	//## Read header
	if(!ReadHeader()){
		cerr<<"FITSReader::Read(): ERROR: Failed to read FITS header, nothing will be read!"<<endl;
		return false;
	}

	Img::MetaData* metadata= SetMetadata();
	if(!metadata){
		cerr<<"FITSReader::Read(): ERROR: Cannot set metadata!"<<endl;
		return false;
	}
	
	//## Set image properties
	image.Reset();
	image.SetBins(fHeader->Nx,-0.5,fHeader->Nx-0.5,fHeader->Ny,-0.5,fHeader->Ny-0.5);
	image.SetMetaData(*metadata);
	//image.SetWCS();

	//## Read image
	for(int j=0;j<fHeader->Ny;j++){
		TVectorD* thisPixelRow= hdu->GetArrayRow(j);
		for(int i=0;i<fHeader->Nx;i++){
			double pixValue= ((*thisPixelRow))(i);
			if(TMath::IsNaN(pixValue) || fabs(pixValue)==TMath::Infinity()) continue;
			image.FillPixel(i,j,pixValue);
		}//end loop columns
		if(thisPixelRow) thisPixelRow->Delete();
	}//end loop row

	return true;
	
}//close FITSReader::Read()


bool FITSReader::ReadTile(Img& image,std::string filename,int ix_min,int ix_max,int iy_min,int iy_max){

	//## Init data
	Init();

	//## Check file
	if(!SetFile(filename)){
		cerr<<"FITSReader::ReadTile(): ERROR: No file given or invalid file, nothing will be read!"<<endl;
		return false;
	}

	//## Read header
	if(!ReadHeader()){
		cerr<<"FITSReader::ReadTile(): ERROR: Image header read failed!"<<endl;
		return false;
	}

	//## Check tile sizes
	if(ix_min<0 || ix_min>=fHeader->Nx || ix_min>=ix_max){
		cerr<<"FITSReader::ReadTile(): ERROR: Invalid min tile X given ("<<ix_min<<")"<<endl;	
		return false;
	}
	if(ix_max<0 || ix_max>=fHeader->Nx || ix_max<=ix_min){
		cerr<<"FITSReader::ReadTile(): ERROR: Invalid max tile Y given ("<<ix_max<<")"<<endl;
		return false;
	}
	if(iy_min<0 || iy_min>=fHeader->Ny || iy_min>=iy_max){
		cerr<<"FITSReader::ReadTile(): ERROR: Invalid min tile X given ("<<iy_min<<")"<<endl;
		return false;
	}
	if(iy_max<0 || iy_max>=fHeader->Ny || iy_max<=iy_min){
		cerr<<"FITSReader::ReadTile(): ERROR: Invalid max tile Y given ("<<iy_max<<")"<<endl;
		return false;
	}

	int TileSizeX= ix_max-ix_min+1;
	int TileSizeY= iy_max-iy_min+1;
		
	Img::MetaData* metadata= SetMetadata();
	if(!metadata){
		cerr<<"FITSReader::ReadTile(): ERROR: Cannot set metadata!"<<endl;
		return 0;
	}
	
	//## Set image properties
	cout<<"ImgReader::ReadTile(): INFO: Set image properties..."<<endl;
	image.Reset();
	image.SetBins(TileSizeX,ix_min-0.5,ix_max+0.5,TileSizeY,iy_min-0.5,iy_max+0.5);	
	image.SetMetaData(*metadata);
	//image.SetWCS();

	//## Read tile
	cout<<"ImgReader::ReadTile(): INFO: Starting reading tile..."<<endl;
	int nFilledPixels= 0;
	
	for(int j=iy_min;j<=iy_max;j++){
		TVectorD* thisPixelRow= hdu->GetArrayRow(j);
		for(int i=ix_min;i<=ix_max;i++){		
			double pixValue= ((*thisPixelRow))(i);
			if(TMath::IsNaN(pixValue) || fabs(pixValue)==TMath::Infinity()) continue;
			image.FillPixel(i,j,pixValue);
			nFilledPixels++;
		}//end loop columns
		if(thisPixelRow) thisPixelRow->Delete();
	}//end loop row

	return true;
	
}//close FITSReader::ReadTile()



bool FITSReader::ReadHeader(){

	//## Open primary HDU from file if not already present
	if(hdu) {
		fHasHeader= false;
		delete hdu;
	}
	
	hdu = new TFITSHDU(fFileInfo->path.c_str());
  if (hdu == 0) {
		cerr<<"FITSReader::ReadHeader(): ERROR: Cannot access the FITS HDU!"<<endl;
		return false;
	}
  	
	cout<<"ImgReader::ReadHeader(): INFO: Header read with success..."<<endl;
	 
	fIsValidHeader= true;

	//##### GET STANDARD & MANDATORY KEYWORDS (if not existing set an invalid header...)
	// Get image size field

	int Nchannels= 0;
	if(hdu->GetKeywordValue("NAXIS")=="") {
		fIsValidHeader= false;
		cerr<<"FITSReader::ReadHeader(): ERROR: Invalid header detected (no NAXIS keyword)!"<<endl;
		return false;
	}			

	Nchannels= hdu->GetKeywordValue("NAXIS").Atoi();//nchannels
	if(Nchannels!=2){
		fIsValidHeader= false;
		cerr<<"FITSReader::ReadHeader(): ERROR: Not supported number of channels ("<<Nchannels<<"), only 2 channel images are supported!"<<endl;
		return false;	
	}	

	int Nx= 0;
	int Ny= 0;
	if(hdu->GetKeywordValue("NAXIS1")=="" || hdu->GetKeywordValue("NAXIS2")=="") {
		fIsValidHeader= false;
		cerr<<"FITSReader::ReadHeader(): ERROR: Invalid header detected (no NAXIS keywords found!"<<endl;
		return false;
	}
	Nx= hdu->GetKeywordValue("NAXIS1").Atoi();//image size X
	Ny= hdu->GetKeywordValue("NAXIS2").Atoi();//image size Y

	// Get pixel content unit
	std::string BUnit= "";
	if(hdu->GetKeywordValue("BUNIT")=="") {
		fIsValidHeader= false;
		cerr<<"FITSReader::ReadHeader(): WARNING: BUNIT keyword not found in header!"<<endl;
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
		cerr<<"FITSReader::ReadHeader(): WARNING: CTYPE keyword not found in header!"<<endl;
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
		cerr<<"FITSReader::ReadHeader(): WARNING: CRPIX keyword not found in header!"<<endl;
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
		cerr<<"FITSReader::ReadHeader(): WARNING: CRVAL keyword not found in header!"<<endl;
	}
	else{
		CenterPixX= hdu->GetKeywordValue("CRVAL1").Atof();
		CenterPixY= hdu->GetKeywordValue("CRVAL2").Atof();
	}
	
	double PixStepX= 0;
	double PixStepY= 0;
	if(hdu->GetKeywordValue("CDELT1")=="" || hdu->GetKeywordValue("CDELT2")==""){
		fIsValidHeader= false;
		cerr<<"FITSReader::ReadHeader(): WARNING: CDELT keyword not found in header!"<<endl;
	}
	else{
		PixStepX= hdu->GetKeywordValue("CDELT1").Atof();
		PixStepY= hdu->GetKeywordValue("CDELT2").Atof();
	}

	//Beam information
	double Bmaj= 0;
	double Bmin= 0;
	double Bpa= 0;
	if(hdu->GetKeywordValue("BMAJ")=="" || hdu->GetKeywordValue("BMIN")==""){
		fIsValidHeader= false;
		cerr<<"FITSReader::ReadHeader(): WARNING: BMAJ/BMIN keywords not found in header!"<<endl;
	}
	else{
		Bmaj= hdu->GetKeywordValue("BMAJ").Atof();
		Bmin= hdu->GetKeywordValue("BMIN").Atof();
		
	}

	if(hdu->GetKeywordValue("BPA")==""){
		cerr<<"FITSReader::ReadHeader(): WARNING: BPA keyword not found in header, setting to zero!"<<endl;
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
		cerr<<"FITSReader::ReadHeader(): WARNING: CROTA keyword not found in header, setting rotation to 0!"<<endl;
	}
	else{
		RotX= hdu->GetKeywordValue("CROTA1").Atof();
		RotY= hdu->GetKeywordValue("CROTA2").Atof();
	}

	// Get Observation Location
	double ObsRA= -999;
	double ObsDEC= -999;
	if(hdu->GetKeywordValue("OBSRA")=="" || hdu->GetKeywordValue("OBSDEC")==""){
		cerr<<"FITSReader::ReadHeader(): WARNING: OBSRA/OBSDEC keywords not found in header, setting them to fake values!"<<endl;
	}
	else{
		ObsRA= hdu->GetKeywordValue("OBSRA").Atof();
		ObsDEC= hdu->GetKeywordValue("OBSDEC").Atof();
	}

	// Get epoch information
	double Epoch= 2000;
	if(hdu->GetKeywordValue("EPOCH")==""){
		cerr<<"FITSReader::ReadHeader(): WARNING: EPOCH keyword not found in header, setting it to 2000!"<<endl;
	}
	else{
		Epoch= hdu->GetKeywordValue("EPOCH").Atof();
	}

	
	if(!fHeader) fHeader= new FITSHeader;
	fHeader->Nx= Nx;
	fHeader->Ny= Ny;
	fHeader->nRec= nRec;
	fHeader->ObsRA= ObsRA;
	fHeader->ObsDEC= ObsDEC;
	fHeader->BUnit= BUnit;
	fHeader->CoordTypeX= CoordTypeX;
	fHeader->CoordTypeY= CoordTypeY;
	fHeader->Cx= CenterPixIdX;
	fHeader->Cy= CenterPixIdY;
	fHeader->Xc= CenterPixX;
	fHeader->Yc= CenterPixY;
	fHeader->dX= PixStepX;
	fHeader->dY= PixStepY;
	fHeader->Bmaj= Bmaj;
	fHeader->Bmin= Bmin;
	fHeader->Bpa= Bpa;
	fHeader->RotX= RotX;
	fHeader->RotY= RotY;
	fHeader->Epoch= Epoch;
	fHasHeader= true;
	
	DumpHeaderInfo();
	
	return true;

}//close FITSReader::ReadHeader()


bool FITSReader::SetFile(std::string path){

	if(!fFileInfo) fFileInfo= new FITSFileInfo;
	fFileInfo->path= path;
	fIsValidFile= false;

	//## Check filename
	if(path==""){
		cerr<<"FITSReader::SetFile(): WARNING: Empty filename given!"<<endl;	
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
					fFileInfo->filename= filePath.filename().string();
					fFileInfo->filename_wext= filePath.stem().string();
					fFileInfo->size= boost::filesystem::file_size(filePath);
        	
					if(filePath.has_extension()) {
						std::string fileExtension= filePath.extension().string();	
						fIsValidFile= true;
						
						if(fileExtension!=".fits"){
							cout<<"FITSReader::SetFile(): ERROR: Unknown or unsupported file extension detected ("<<fileExtension<<")..."<<endl;
							return false;
						}
	
						DumpFileInfo();

					}//close if
					else{
						cerr << "FITSReader::SetFile(): ERROR: Given file without extension!"<<endl;
						return false;
					}
				}//close if has filename
				else{
					cerr<<"FITSReader::SetFile(): ERROR: Given path has not a filename!"<<endl;
					return false;	
				}
			}//close if
      else if (boost::filesystem::is_directory(filePath)){
        cerr << "FITSReader::SetFile(): ERROR: Given file ("<<filePath<<") is a directory!"<<endl;
				return false;
      }
      else{
        cerr << "FITSReader::SetFile(): ERROR: Given file ("<<filePath<<") exists, but is neither a regular file nor a directory!"<<endl;
				return false;
			}
    }//close if exists
    else{
      cerr << "FITSReader::SetFile(): ERROR: Given file ("<<filePath<<") does not exist!"<<endl;
			return false;
		}
  }//close try block

  catch (const boost::filesystem::filesystem_error& ex) {
    cout << ex.what() << '\n';
		return false;
  }

	return true;

}//close FITSReader::SetFile()


void FITSReader::DumpHeaderInfo(){

	if(!fHeader) return;

	cout<<"*** IMAGE META-DATA ***"<<endl;
	cout<<"Image Size: "<<fHeader->Nx<<"x"<<fHeader->Ny<<" pixels, nRec="<<fHeader->nRec<<endl;
	cout<<"Obs Coords: ("<<fHeader->ObsRA<<","<<fHeader->ObsDEC<<")"<<endl;
	cout<<"BUnit: "<<fHeader->BUnit<<endl;
	cout<<"Coords Type: ("<<fHeader->CoordTypeX<<","<<fHeader->CoordTypeY<<")"<<endl;
	cout<<"PixelCoordCenter: ("<<fHeader->Cx<<","<<fHeader->Cy<<")"<<endl;
	cout<<"CoordCenter: ("<<fHeader->Xc<<","<<fHeader->Yc<<")"<<endl;
	cout<<"PixelStep: ("<<fHeader->dX<<","<<fHeader->dY<<")"<<endl;
	cout<<"BeamSize: ("<<fHeader->Bmaj<<","<<fHeader->Bmin<<","<<fHeader->Bpa<<")"<<endl;
	cout<<"Rot: ("<<fHeader->RotX<<","<<fHeader->RotY<<")"<<endl;
	cout<<"Epoch: "<<fHeader->Epoch<<endl;
	cout<<"***********************"<<endl;	

}//close FITSReader::DumpHeaderInfo()


void FITSReader::DumpFileInfo(){	
	if(!fFileInfo) return;		
	cout <<"*** IMAGE FILE INFO ***"<<endl;
	cout << "File Path: "<<fFileInfo->path<<endl;
	cout << "File Name: "<<fFileInfo->filename<<endl;
	cout << "File WExt: "<<fFileInfo->filename_wext<<endl;
	cout << "File size (kB): "<<fFileInfo->size/1024.<<endl;
	cout <<"***********************"<<endl;
}	

