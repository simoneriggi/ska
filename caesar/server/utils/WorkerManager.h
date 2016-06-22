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
* @file WorkerManager.h
* @class WorkerManager
* @brief WorkerManager
*
* WorkerManager
* @author S. Riggi
* @date 16/06/2016
*/
#ifndef WorkerManager_H
#define WorkerManager_H

#include <Logger.h>

#include <tango.h>


#include <chrono>
#include <thread>
#include <mutex>

namespace Caesar {

class WorkerManager {

	public :
		//Constructor
		WorkerManager();
		
		//Destructor
		~WorkerManager();

	public:
		//Has worker
		bool HasWorker(const std::string& device_name);

		//Add worker to the list
		int AddWorker(const std::string& device_name);
		
		//Flush workers
		int FlushWorkers();

		//Get free worker group
		int GetFreeWorkers(Tango::Group& group,int nMaxWorkers=-1,bool requireExactly=false);		

		//Get worker names
		void GetWorkerNames(std::vector<std::string>& worker_names);

	private:
		void Init();

	private:
		mutable std::mutex m_mutex;
		Tango::Group* m_workers;
		Tango::Group* m_freeWorkers;
		Tango::Group* m_busyWorkers;

};//close WorkerManager

}//close namespace

#endif
