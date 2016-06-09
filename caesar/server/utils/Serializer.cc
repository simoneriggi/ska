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
* @file Serializer.cc
* @class Serializer
* @brief Serializer class
*
* Class for serializing objects
* @author S. Riggi
* @date 20/01/2015
*/

#include <Serializer.h>
#include <Logger.h>

//# MSG PACK
#ifdef BUILD_CAESAR_SERVER
	#include <msgpack.hpp>
#endif


#include <string>
#include <ctime>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <sstream>

#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <map>
#include <exception>

#include <chrono>

using namespace std;

ClassImp(Caesar::Serializer)

namespace Caesar {

Serializer::Serializer() {	
	
}   


Serializer::~Serializer(){

}

#ifdef BUILD_CAESAR_SERVER
int Serializer::SourceToString(Source* source,std::string& msg){

	msgpack::sbuffer buffer;
	if(SourceToBuffer(source,buffer)<0) return -1;

	msg= std::string(buffer.data(),buffer.size());

	return 0;

}//close SourceToString()


int Serializer::SourceToDevString(Source* source,Tango::DevString& msg){

	msgpack::sbuffer buffer;
	if(SourceToBuffer(source,buffer)<0)	return -1;

	msg= CORBA::string_alloc(buffer.size());
	memcpy (msg, buffer.data(), buffer.size());	
	msg[buffer.size()]= 0;

	return 0;

}//close SourceToDevString()

int Serializer::SourceToBuffer(Source* source,msgpack::sbuffer& buffer){

	try{	
		msgpack::pack(&buffer,*source);
	}
	catch(std::exception const & e) {
		ERROR_LOG("Source encoding failed (err: "<<e.what()<<")");
		return -1;
	}

	return 0;

}//close SourceToBuffer()
#endif


/*
bool MessageUtils::encodeFromMsgPack(Buffer buffer,Message& msg){

	//## Check data integrity
	if(!buffer.msg || buffer.size<=0) return false;

	//## Encode to message
	try {
		msgpack::unpacked unpacked_msg;
  	msgpack::unpack(&unpacked_msg, buffer.msg, buffer.size);	
  	unpacked_msg.get().convert(&msg);
	
		//## Check message validity
		if(!msg.isValid()) return false;
	}
	catch(std::exception const & e) {
		cout << "ERROR: message encoding failed with status "<<e.what() <<endl;
		return false;
	}

	return true;

}//close MessageUtils::encodeFromMsgPack()


bool MessageUtils::encodeFromMsgPack(SBuffer buffer,Message& msg){

	//## Check data integrity
	if(buffer.msg.empty() || buffer.size<=0) return false;

	//## Encode to message
	try {
		msgpack::unpacked unpacked_msg;
  	msgpack::unpack(&unpacked_msg, (buffer.msg).c_str(), buffer.size);	
  	unpacked_msg.get().convert(&msg);
	
		//## Check message validity
		if(!msg.isValid()) return false;
	}
	catch(std::exception const & e) {
		cout << "ERROR: message encoding failed with status "<<e.what() <<endl;
		return false;
	}

	return true;

}//close MessageUtils::encodeFromMsgPack()


bool MessageUtils::encodeFromMsgPackCollection(BufferCollection buffers,MessageCollection& msgs){

	bool status= true;
	msgs.clear();
	msgs.resize(0);

	for(unsigned int i=0;i<buffers.size();i++){
		Message thisMsg;
		bool isValidDecoding= MessageUtils::encodeFromMsgPack(buffers[i],thisMsg);
		if(!isValidDecoding) status= false; 
		msgs.push_back(thisMsg);
	}//end loop buffers

	return status;

}//close MessageUtils::encodeFromMsgPackCollection()


bool MessageUtils::encodeFromMsgPackCollection(SBufferCollection buffers,MessageCollection& msgs){

	bool status= true;
	msgs.clear();
	msgs.resize(0);

	for(unsigned int i=0;i<buffers.size();i++){
		Message thisMsg;
		bool isValidDecoding= MessageUtils::encodeFromMsgPack(buffers[i],thisMsg);
		if(!isValidDecoding) status= false; 
		msgs.push_back(thisMsg);
	}//end loop buffers

	return status;

}//close MessageUtils::encodeFromMsgPackCollection()



bool MessageUtils::encodeFromMsgPack(char* buffer,int bufsize,Message& msg){

	//## Check data integrity
	if(!buffer || bufsize<=0) return false;

	//## Encode to message
	try {
		msgpack::unpacked unpacked_msg;
  	msgpack::unpack(&unpacked_msg, buffer, bufsize);	
  	unpacked_msg.get().convert(&msg);
	
		//## Check message validity
		if(!msg.isValid()) return false;
	}
	catch(std::exception const & e) {
		cout << "ERROR: message encoding failed with status "<<e.what() <<endl;
		return false;
	}

	return true;

}//close encodeFromMsgPack()

*/

}//close namespace
