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
* @file FITSReader.cc
* @class FITSReader
* @brief FITSReader
*
* FITS Image Reader class
* @author S. Riggi
* @date 20/01/2015
*/

#include <FITSReader.h>
#include <SysUtils.h>
#include <Img.h>
#include <Logger.h>
#include <CodeUtils.h>

#include <TObject.h>
#include <TFITS.h>
#include <TMath.h>
#include <TVectorD.h>


#ifdef OPENMP_ENABLED
  #include <omp.h>
#else
  #define omp_get_thread_num() 0
#endif

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
using namespace std;

ClassImp(Caesar::FITSHeader)
ClassImp(Caesar::FITSFileInfo)
ClassImp(Caesar::FITSReader)

namespace Caesar {


FITSReader::FITSReader() {

}//close costructor


FITSReader::~FITSReader() {

}//close destructor


bool FITSReader::ReadHeader(TFITSHDU* hdu,FITSFileInfo& fits_info){

	if(!hdu) {
		ERROR_LOG("Null ptr to given FITS HDU!");
		return false;
	}

	//##### GET STANDARD & MANDATORY KEYWORDS (if not existing set an invalid header...)
	// Get image size field
	int Nchannels= 0;
	if(hdu->GetKeywordValue("NAXIS")=="") {
		WARN_LOG("Invalid header detected (no NAXIS keyword)!");
		return false;
	}			

	Nchannels= hdu->GetKeywordValue("NAXIS").Atoi();//nchannels
	if(Nchannels!=2){
		ERROR_LOG("Unsupported number of channels ("<<Nchannels<<"), only 2 channel images are supported!");
		return false;	
	}	

	int Nx= 0;
	int Ny= 0;
	if(hdu->GetKeywordValue("NAXIS1")=="" || hdu->GetKeywordValue("NAXIS2")=="") {
		ERROR_LOG("Invalid header detected (no NAXIS keywords found!");
		return false;
	}
	Nx= hdu->GetKeywordValue("NAXIS1").Atoi();//image size X
	Ny= hdu->GetKeywordValue("NAXIS2").Atoi();//image size Y

	// Get pixel content unit
	std::string BUnit= "";
	if(hdu->GetKeywordValue("BUNIT")=="") {
		WARN_LOG("BUNIT keyword not found in header!");
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
		WARN_LOG("CTYPE keyword not found in header!");
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
		WARN_LOG("CRPIX keyword not found in header!");
	}
	else{
		CenterPixIdX= hdu->GetKeywordValue("CRPIX1").Atof();
		CenterPixIdY= hdu->GetKeywordValue("CRPIX2").Atof();
	}
	
	//Get reference pixel coordinates (RA/DEC)
	double CenterPixX= 0;
	double CenterPixY= 0;
	if(hdu->GetKeywordValue("CRVAL1")=="" || hdu->GetKeywordValue("CRVAL2")==""){
		WARN_LOG("CRVAL keyword not found in header!");
	}
	else{
		CenterPixX= hdu->GetKeywordValue("CRVAL1").Atof();
		CenterPixY= hdu->GetKeywordValue("CRVAL2").Atof();
	}
	
	double PixStepX= 0;
	double PixStepY= 0;
	if(hdu->GetKeywordValue("CDELT1")=="" || hdu->GetKeywordValue("CDELT2")==""){
		WARN_LOG("CDELT keyword not found in header!");
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
		WARN_LOG("BMAJ/BMIN keywords not found in header!");
	}
	else{
		Bmaj= hdu->GetKeywordValue("BMAJ").Atof();
		Bmin= hdu->GetKeywordValue("BMIN").Atof();
	}

	if(hdu->GetKeywordValue("BPA")==""){
		WARN_LOG("BPA keyword not found in header, setting to zero!");
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
		WARN_LOG("CROTA keyword not found in header, setting rotation to 0!");
	}
	else{
		RotX= hdu->GetKeywordValue("CROTA1").Atof();
		RotY= hdu->GetKeywordValue("CROTA2").Atof();
	}

	// Get Observation Location
	double ObsRA= -999;
	double ObsDEC= -999;
	if(hdu->GetKeywordValue("OBSRA")=="" || hdu->GetKeywordValue("OBSDEC")==""){
		WARN_LOG("OBSRA/OBSDEC keywords not found in header, setting them to fake values!");
	}
	else{
		ObsRA= hdu->GetKeywordValue("OBSRA").Atof();
		ObsDEC= hdu->GetKeywordValue("OBSDEC").Atof();
	}

	// Get epoch information
	double Epoch= 2000;
	if(hdu->GetKeywordValue("EPOCH")==""){
		WARN_LOG("EPOCH keyword not found in header, setting it to 2000!");
	}
	else{
		Epoch= hdu->GetKeywordValue("EPOCH").Atof();
	}

	
	//## Fill header info in FITS INFO struct
	(fits_info.header).Nx= Nx;
	(fits_info.header).Ny= Ny;
	(fits_info.header).nRec= nRec;
	(fits_info.header).ObsRA= ObsRA;
	(fits_info.header).ObsDEC= ObsDEC;
	(fits_info.header).BUnit= BUnit;
	(fits_info.header).CoordTypeX= CoordTypeX;
	(fits_info.header).CoordTypeY= CoordTypeY;
	(fits_info.header).Cx= CenterPixIdX;
	(fits_info.header).Cy= CenterPixIdY;
	(fits_info.header).Xc= CenterPixX;
	(fits_info.header).Yc= CenterPixY;
	(fits_info.header).dX= PixStepX;
	(fits_info.header).dY= PixStepY;
	(fits_info.header).Bmaj= Bmaj;
	(fits_info.header).Bmin= Bmin;
	(fits_info.header).Bpa= Bpa;
	(fits_info.header).RotX= RotX;
	(fits_info.header).RotY= RotY;
	(fits_info.header).Epoch= Epoch;
	
	return true;

}//close FITSReader::ReadHeader()


TFITSHDU* FITSReader::ReadFile(std::string filename,Caesar::FITSFileInfo& fits_info,bool checkFile){
	
	//## Check file
	if(checkFile && !SysUtils::CheckFile(filename,fits_info.info,true,".fits")){
		ERROR_LOG("Failed to read file "<<filename<<"!");
		return 0;
	}

	//## Open primary HDU from file if not already present
	TFITSHDU* hdu= 0;
	try {
		hdu = new TFITSHDU(filename.c_str());
	}
	catch(...){
		ERROR_LOG("Failed to read FITS HDU!");
		return 0;
	}
  if (!hdu) {
		ERROR_LOG("Cannot access the FITS HDU!");
		return 0;
	}
	
	//## Read file header	
	INFO_LOG("Reading FITS file header...");
	if(!ReadHeader(hdu,fits_info)){
		ERROR_LOG("Failed to read header from FITS file!");
		return 0;
	}

	return hdu;
	
}//close ReadFile()


int FITSReader::Read(std::string filename,Caesar::Img& image,Caesar::FITSFileInfo& fits_info,bool checkFile){

	//## Read file
	TFITSHDU* hdu= ReadFile(filename,fits_info,checkFile);
	if(!hdu){
		ERROR_LOG("Failed to read file "<<filename<<"!");
		return -1;
	}
	int Nx= (fits_info.header).Nx;
	int Ny= (fits_info.header).Ny;

	//Compute image meta-data	
	Caesar::ImgMetaData* metadata= new Caesar::ImgMetaData;
	metadata->SetFITSCards(fits_info);

	//## Set image properties
	image.Reset();
	image.SetBins(Nx,-0.5,Nx-0.5,Ny,-0.5,Ny-0.5);
	image.SetMetaData(metadata);

	//## Fill image
	#ifdef OPENMP_ENABLED
		#pragma omp parallel
		{
			INFO_LOG("Starting image filling in thread "<<omp_get_thread_num()<<" (nthreads="<<SysUtils::GetOMPThreads()<<") ...");
			
			#pragma omp for
			for(int j=0;j<Ny;j++){
				TVectorD* thisPixelRow= hdu->GetArrayRow(j);
				for(int i=0;i<Nx;i++){
					double pixValue= ((*thisPixelRow))(i);
					if(TMath::IsNaN(pixValue) || fabs(pixValue)==TMath::Infinity()) continue;
					image.FillPixel(i,j,pixValue);
				}//end loop columns
				if(thisPixelRow) thisPixelRow->Delete();
			}//end loop row

		}//close parallel section
	#else
		for(int j=0;j<Ny;j++){
			TVectorD* thisPixelRow= hdu->GetArrayRow(j);
			for(int i=0;i<Nx;i++){
				double pixValue= ((*thisPixelRow))(i);
				if(TMath::IsNaN(pixValue) || fabs(pixValue)==TMath::Infinity()) continue;
				image.FillPixel(i,j,pixValue);
			}//end loop columns
			if(thisPixelRow) thisPixelRow->Delete();
		}//end loop row
	#endif

	return 0;
	
}//close FITSReader::Read()


int FITSReader::ReadTileFast(std::string filename,Caesar::Img& image,Caesar::FITSFileInfo& fits_info,int ix_min,int ix_max,int iy_min,int iy_max,bool checkFile){

	//Check file
	if(checkFile && !SysUtils::CheckFile(filename,fits_info.info,true,".fits")){
		ERROR_LOG("Failed to read file "<<filename<<"!");
		return -1;
	}

	//Append coord filters to input file (see CFITSIO for details, used internally by ROOT)
	std::stringstream ss;
	ss<<filename<<"["<<ix_min<<":"<<ix_max<<","<<iy_min<<":"<<iy_max<<"]";
	std::string filename_wfilter= ss.str();
	
	//Read file
	TFITSHDU* hdu= ReadFile(filename_wfilter,fits_info,false);
	if(!hdu){
		ERROR_LOG("Failed to read file with filter "<<filename_wfilter<<"!");
		return -1;
	}
	int Nx= (fits_info.header).Nx;
	int Ny= (fits_info.header).Ny;
	int TileSizeX= ix_max-ix_min+1;
	int TileSizeY= iy_max-iy_min+1;
	if(Nx!=TileSizeX || Ny!=TileSizeY){
		ERROR_LOG("Failed to read file with filter "<<filename_wfilter<<"!");
		return -1;
	}
	INFO_LOG("(Nx,Ny)=("<<Nx<<","<<Ny<<")");

	//Compute image meta-data	
	Caesar::ImgMetaData* metadata= new Caesar::ImgMetaData;
	metadata->SetFITSCards(fits_info);

	
	//## Set image properties
	DEBUG_LOG("Set image properties...");
	image.Reset();
	image.SetBins(TileSizeX,ix_min-0.5,ix_max+0.5,TileSizeY,iy_min-0.5,iy_max+0.5);	
	image.SetMetaData(metadata);
	

	//## Fill tile image
	DEBUG_LOG("Start filling tile...");
	int nFilledPixels= 0;
	
	for(int j=0;j<Ny;j++){
		TVectorD* thisPixelRow= hdu->GetArrayRow(j);
		for(int i=0;i<Nx;i++){		
			double pixValue= ((*thisPixelRow))(i);
			if(TMath::IsNaN(pixValue) || fabs(pixValue)==TMath::Infinity()) continue;
			image.FillPixel(i+ix_min,j+iy_min,pixValue);
			nFilledPixels++;
		}//end loop columns
		if(thisPixelRow) thisPixelRow->Delete();
	}//end loop row

	return 0;

}//close ReadTileFast()


int FITSReader::ReadTile(std::string filename,Caesar::Img& image,Caesar::FITSFileInfo& fits_info,int ix_min,int ix_max,int iy_min,int iy_max,bool checkFile){

	//## Read file
	TFITSHDU* hdu= ReadFile(filename,fits_info,checkFile);
	if(!hdu){
		ERROR_LOG("Failed to read file "<<filename<<"!");
		return -1;
	}
	int Nx= (fits_info.header).Nx;
	int Ny= (fits_info.header).Ny;
	

	//## Check tile sizes
	if(ix_min<0 || ix_min>=Nx || ix_min>=ix_max){
		ERROR_LOG("Invalid min tile X given ("<<ix_min<<")");
		return -1;
	}
	if(ix_max<0 || ix_max>=Nx || ix_max<=ix_min){
		ERROR_LOG("Invalid max tile Y given ("<<ix_max<<")");
		return -1;
	}
	if(iy_min<0 || iy_min>=Ny || iy_min>=iy_max){
		ERROR_LOG("Invalid min tile X given ("<<iy_min<<")");
		return -1;
	}
	if(iy_max<0 || iy_max>=Ny || iy_max<=iy_min){
		ERROR_LOG("Invalid max tile X given ("<<iy_max<<")");
		return -1;
	}

	int TileSizeX= ix_max-ix_min+1;
	int TileSizeY= iy_max-iy_min+1;
		
	//Compute image meta-data	
	Caesar::ImgMetaData* metadata= new Caesar::ImgMetaData;
	metadata->SetFITSCards(fits_info);

	
	//## Set image properties
	DEBUG_LOG("Set image properties...");
	image.Reset();
	image.SetBins(TileSizeX,ix_min-0.5,ix_max+0.5,TileSizeY,iy_min-0.5,iy_max+0.5);	
	image.SetMetaData(metadata);
	//image.SetWCS();

	//## Read tile
	DEBUG_LOG("Starting reading tile...");
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

	return 0;
	
}//close FITSReader::ReadTile()


}//close namespace

