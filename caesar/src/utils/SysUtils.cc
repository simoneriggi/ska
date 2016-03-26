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
		cerr<<"SysUtils::CheckFile(): WARNING: Empty filename given!"<<endl;	
		return false;
	}
	if(!(&info)){
		cerr<<"SysUtils::CheckFile(): WARNING: Null ptr to file info struct given!"<<endl;	
		return false;
	}

	//Check if file actually exists on filesystem
  try {
		//Check if file exists on filesystem
		boost::filesystem::path file_path(path.c_str());
		if (!boost::filesystem::exists(file_path)){
			cerr<<"SysUtils::CheckFile(): ERROR: File "<<path<<" not found in local filesystem!"<<endl;
			return false;
		}
		if (!boost::filesystem::is_regular_file(file_path)){
			cerr<<"SysUtils::CheckFile(): ERROR: File "<<path<<" is not a regular file!"<<endl;
			return false;
		}
		if (boost::filesystem::is_directory(file_path)){
    	cerr << "SysUtils::CheckFile(): ERROR: File ("<<file_path<<") is a directory!"<<endl;
			return false;
    }
	
		//Get filename and extension
		if(!file_path.has_filename()){
			cerr << "SysUtils::CheckFile(): ERROR: File ("<<file_path<<") does not have a filename!"<<endl;
			return false;
		}

		//Set file info
		info.filename= file_path.filename().string();
		info.filename_wext= file_path.stem().string();
		info.size= boost::filesystem::file_size(file_path);
        	
		//Check extension
		if(!file_path.has_extension()){
			cerr << "SysUtils::CheckFile(): ERROR: Given file without extension!"<<endl;
			return false;
		}

		std::string file_extension= file_path.extension().string();	
		info.extension= file_extension;
						
		if(match_extension && file_extension!=extension){
			cerr<<"SysUtils::CheckFile(): ERROR: Invalid file extension detected ("<<file_extension<<"!="<<extension<<")..."<<endl;
			return false;
		}
	
		//Dump file info
		info.Print();

  }//close try block

  catch (const boost::filesystem::filesystem_error& ex) {
    cerr<<"SysUtils::CheckFile(): ERROR: Exception detected ("<<ex.what()<<")!"<<endl;
		return false;
  }

	return true;

}//close CheckFile()

}//close namespace



