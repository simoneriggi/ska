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
* @file Serializer.h
* @class Serializer
* @brief Serializer class
*
* Class for serializing objects
* @author S. Riggi
* @date 20/01/2015
*/
#ifndef Serializer_H
#define Serializer_H

#include <Source.h>

#ifdef BUILD_CAESAR_SERVER
	#include <tango.h>
	#include <msgpack.hpp>
#endif

//# JSON CPP
//#include <json/json.h>

#include <string>
#include <ctime>
#include <iostream>
#include <stdio.h>
#include <iomanip>

#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <map>

namespace Caesar {

class Serializer : public TObject {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		Serializer();
		/**
		* \brief Class destructor: free allocated memory
		*/
		~Serializer();

	public: 
		
		#ifdef BUILD_CAESAR_SERVER
			static int SourceToString(Source* source,std::string& msg);
			static int SourceToDevString(Source* source,Tango::DevString& msg);
			static int SourceToBuffer(Source* source,msgpack::sbuffer& buffer);
		#endif

		
	public:
    
		ClassDef(Serializer,1)


};

#ifdef __MAKECINT__
#pragma link C++ class Serializer+;
#endif

}//close namespace

#endif



