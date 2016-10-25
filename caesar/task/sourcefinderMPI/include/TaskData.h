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
* @file WorkerData.h
* @class WorkerData
* @brief WorkerData class
*
* Class for serializing objects
* @author S. Riggi
* @date 20/01/2015
*/
#ifndef _TASK_DATA_H
#define _TASK_DATA_H


#include <Source.h>

#include <Source.pb.h>
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
#include <time.h>
#include <sys/time.h>

namespace Caesar {



class TaskData : public TObject {
	public:
		//Standard constructor
		TaskData();
				
		//Copy constructor
		TaskData(const TaskData& data);

		//Destructor
		virtual ~TaskData();

	public:
		// Operator =
		TaskData& operator=(const TaskData& data);

		//Copy method
		void Copy(TObject &obj) const;

		//Serialize methods
		int SerializeToProtobuf(SourcePB::TaskData& taskDataPB);
		int EncodeFromProtobuf(const SourcePB::TaskData& taskDataPB);
		
		/*
		int SerializeToJson(Json::Value& jsonObj);
		int SerializeToJsonString(std::string& jsonString,bool isMinified=true);		
		int EncodeFromJson(const Json::Value& jsonObj);
		int EncodeFromJsonString(std::string& jsonString);
		*/

	public: 

		//Task info		
		std::string filename; 
		std::string jobId;
		long int workerId;
		long int taskId;
		std::vector<long int> neighborTaskId;
		std::vector<long int> neighborWorkerId;
	
		//long int IdX;
		//long int IdY;
		long int ix_min;
		long int ix_max;
		long int iy_min;
		long int iy_max;

		Source* source;
		std::vector<Source*> sources;
		std::vector<Source*> ext_sources;
		
	//ClassDef(TaskData,1)
};

/*
#ifdef __MAKECINT__
#pragma link C++ class TaskData+;
#endif
*/

}//close namespace

#endif



