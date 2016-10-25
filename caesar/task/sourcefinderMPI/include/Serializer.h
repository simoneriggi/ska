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

#include <Source.pb.h>
#include <TaskData.h>

#include <Source.h>
#include <Contour.h>

#include <TVector2.h> 

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

class SBuffer : public TObject {
	public:
		SBuffer(){};
		virtual ~SBuffer(){};
	public: 
		std::string data;
		long int size;

	ClassDef(SBuffer,1)
};

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
		
			//## SOURCE SERIALIZATION
			//Source --> Buffer
			static int EncodePointToProtobuf(SourcePB::Point& point_pb,TVector2& point);
			static int EncodeContourToProtobuf(SourcePB::Contour& contour_pb,Contour* contour);
			static int EncodePixelToProtobuf(SourcePB::Pixel& pixel_pb,Pixel* pixel);
			static int EncodeBlobToProtobuf(SourcePB::Blob& blob_pb,Source* source);
			static int EncodeSourceToProtobuf(SourcePB::Source& source_pb,Source* source);		
			static int SourceToBuffer(SBuffer& buffer,Source* source);
			
			//Buffer --> Source
			static int EncodeProtobufToPoint(TVector2& point,const SourcePB::Point& point_pb);
			static int EncodeProtobufToContour(Contour& contour,const SourcePB::Contour& contour_pb);
			static int EncodeProtobufToPixel(Pixel& pixel,const SourcePB::Pixel& pixel_pb);	
			static int EncodeProtobufToBlob(Source& source,const SourcePB::Blob& blob_pb);
			static int EncodeProtobufToSource(Source& source,const SourcePB::Source& source_pb);			
			static int BufferToSource(Source& source,SBuffer& buffer);

			//## WORKER DATA SERIALIZATION ###
			//TaskData --> Buffer
			static int EncodeTaskDataToProtobuf(SourcePB::TaskData& taskData_pb,TaskData* taskData);
			static int TaskDataToBuffer(SBuffer& buffer,TaskData* taskData);
			static char* TaskDataToCharArray(long int& buffer_size,TaskData* taskData);

			static int EncodeTaskDataCollectionToProtobuf(SourcePB::TaskDataCollection& taskDataCollection_pb,std::vector<TaskData*> taskDataCollection);
			static int TaskDataCollectionToBuffer(SBuffer& buffer,std::vector<TaskData*> taskDataCollection);
			static char* TaskDataCollectionToCharArray(long int& buffer_size,std::vector<TaskData*> taskDataCollection);

			//Buffer --> TaskData
			static int EncodeProtobufToTaskData(TaskData& taskData,const SourcePB::TaskData& taskData_pb);
			static int BufferToTaskData(TaskData& taskData,SBuffer& buffer);
			static int CharArrayToTaskData(TaskData& taskData,char* buffer,long int buffer_size);

			static int EncodeProtobufToTaskDataCollection(std::vector<TaskData*>& taskDataCollection,const SourcePB::TaskDataCollection& taskDataCollection_pb,bool isTaskCollectionPreAllocated=false);
			static int BufferToTaskDataCollection(std::vector<TaskData*>& taskDataCollection,SBuffer& buffer,bool isTaskCollectionPreAllocated=false);
			static int CharArrayToTaskDataCollection(std::vector<TaskData*>& taskDataCollection,char* buffer,long int buffer_size,bool isTaskCollectionPreAllocated=false);

	public:
    
		//ClassDef(Serializer,1)

};

/*
#ifdef __MAKECINT__
#pragma link C++ class Serializer+;
#endif
*/

}//close namespace

#endif



