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
* @file SysUtils.cc
* @class SysUtils
* @brief Utility functions for system tasks
*
* Utility functions for system tasks
* @author S. Riggi
* @date 23/08/2010
*/


#include <SysUtils.h>
#include <CodeUtils.h>
#include <Logger.h>


#include <fitsio.h>

#include <boost/filesystem.hpp>

#include <TObject.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <ctime>

using namespace std;

ClassImp(Caesar::SysUtils)
ClassImp(Caesar::FileInfo)

namespace Caesar {

SysUtils::SysUtils(){

}

SysUtils::~SysUtils(){

}

bool SysUtils::CheckFile(std::string path,Caesar::FileInfo& info,bool match_extension,std::string extension){

	//Check input file path
	if(path==""){	
		WARN_LOG("Empty filename given!");
		return false;
	}
	if(!(&info)){
		ERROR_LOG("Null ptr to file info struct given!");
		return false;
	}

	//Check if file actually exists on filesystem
  try {
		//Check if file exists on filesystem
		boost::filesystem::path file_path(path.c_str());
		if (!boost::filesystem::exists(file_path)){
			ERROR_LOG("File "<<path<<" not found in local filesystem!");
			return false;
		}
		if (!boost::filesystem::is_regular_file(file_path)){
			ERROR_LOG("File "<<path<<" is not a regular file!");
			return false;
		}
		if (boost::filesystem::is_directory(file_path)){
    	ERROR_LOG("File ("<<file_path<<") is a directory!");
			return false;
    }
	
		//Get filename and extension
		if(!file_path.has_filename()){
			ERROR_LOG("File ("<<file_path<<") does not have a filename!");
			return false;
		}

		//Set file info
		info.filename= file_path.filename().string();
		info.filename_wext= file_path.stem().string();
		info.size= boost::filesystem::file_size(file_path);
        	
		//Check extension
		if(!file_path.has_extension()){
			ERROR_LOG("Given file without extension!");
			return false;
		}

		std::string file_extension= file_path.extension().string();	
		info.extension= file_extension;
						
		if(match_extension && file_extension!=extension){
			ERROR_LOG("Invalid file extension detected ("<<file_extension<<"!="<<extension<<")...");
			return false;
		}
	
		//Dump file info
		//info.Print();
		std::string info_printable= info.GetPrintable();
		INFO_LOG(info_printable);

  }//close try block

  catch (const boost::filesystem::filesystem_error& ex) {
    ERROR_LOG("Exception detected while checking file (err: "<<ex.what()<<")!");
		return false;
  }

	return true;

}//close CheckFile()

int SysUtils::GetFITSImageSize(const std::string& filename,long int& Nx,long int& Ny){

	//Init	
	Nx= 0;
	Ny= 0;

	//Open file
  int status = 0;
	fitsfile *fptr;//pointer to the FITS file, defined in fitsio.h
  fits_open_file(&fptr, filename.c_str(), READONLY, &status);
  if(status){
		ERROR_LOG("Failed to open given fits file "<<filename<<"!");	
		return -1;
	}

  //Read the NAXIS1 and NAXIS2 keyword to get image size
	int nfound;
	long naxes[2];
	fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status);
  if (status){
		ERROR_LOG("Failed to get NAXIS keyword from given fits file "<<filename<<"!");	
		return -1;
	}
        
	Nx= naxes[0];
	Ny= naxes[1];
	
	return 0;

}//close GetFITSImageSize()

}//close namespace



