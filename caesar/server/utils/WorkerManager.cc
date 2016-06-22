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
* @file WorkerManager.cc
* @class WorkerManager
* @brief WorkerManager class
*
* WorkerManager class
* @author S. Riggi
* @date 20/01/2015
*/

#include <WorkerManager.h>
#include <Logger.h>

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


namespace Caesar {

//Standard Constructor
WorkerManager::WorkerManager() {	
	Init();
}


//Destructor
WorkerManager::~WorkerManager(){
	
	/*
	if(m_freeWorkers){
		delete m_freeWorkers;
		m_freeWorkers= 0;
	}
	if(m_busyWorkers){
		delete m_busyWorkers;
		m_busyWorkers= 0;
	}
	*/
	if(m_workers){
		delete m_workers;
		m_workers= 0;
	}

}//close destructor

void WorkerManager::Init(){

	m_workers= 0;
	m_workers= new Tango::Group("Workers");

	m_freeWorkers= 0;
	m_freeWorkers= new Tango::Group("FreeWorkers");

	m_busyWorkers= 0;
	m_busyWorkers= new Tango::Group("BusyWorkers");
	
	m_workers->add(m_freeWorkers);	
	m_workers->add(m_busyWorkers);	

}//close Init()


bool WorkerManager::HasWorker(const std::string& device_name){
	
	bool hasWorker= false;
	try{
		hasWorker= m_workers->contains(device_name);
	}
	catch(const Tango::DevFailed& e){
		ERROR_LOG("Exception occurred while checking if worker ("<<device_name<<") exists in Group returning false!");
		hasWorker= false;
	}
	return hasWorker;

}//close HasWorker()

int WorkerManager::AddWorker(const std::string& device_name){

	std::lock_guard<std::mutex> lock(m_mutex);

	//Check given device name
	if(device_name=="") return -1;
	
	//Check if device already exists
	INFO_LOG("Check if the device already exist in the list...");
	if(HasWorker(device_name)){
		WARN_LOG("Worker ("<<device_name<<") already added!");
		return 0;
	}

	//Add worker in Group (add in busy group initially)
	try{
		m_workers->get_group("BusyWorkers")->add(device_name);
	}
	catch(const Tango::DevFailed& e){
		ERROR_LOG("Exception occurred while adding worker "<<device_name<<" to Group...");
		return -1;
	}
	
	return 0;

}//close AddWorker()


int WorkerManager::FlushWorkers(){

	std::lock_guard<std::mutex> lock(m_mutex);

	try{
		m_workers->remove_all();
	}
	catch(const Tango::DevFailed& e){
		ERROR_LOG("Exception occurred while removing all workers in Group...");
		return -1;
	}

	return 0;

}//close FlushWorkers()


int WorkerManager::GetFreeWorkers(Tango::Group& group,int nMaxWorkers,bool requireExactly){
	
	std::lock_guard<std::mutex> lock(m_mutex);

	//Check for free workers
	int nFreeWorkers= static_cast<int>(m_workers->get_group("FreeWorkers")->get_size());
	if(nFreeWorkers<=0) {
			INFO_LOG("No free workers available!");
			return -1;
	}
	
	//Return all free worker group in this case
	std::vector<std::string> freeWorkersNames= m_workers->get_group("FreeWorkers")->get_device_list();
	if(nMaxWorkers==-1){
		try {		
			group.add(freeWorkersNames);
		}
		catch(Tango::DevFailed& e){
			ERROR_LOG("Failed to add all workers in group!");
			return -1;
		}
		return 0;
	}//close if
	
	
	if(nMaxWorkers>nFreeWorkers && requireExactly){
		WARN_LOG("Number of requested workers exceed the available free workers!");
		return -1;
	}

	//Return group
	try {	
		for(int i=0;i<nMaxWorkers;i++){
			group.add(freeWorkersNames[i]);
		}
	}
	catch(Tango::DevFailed& e){
		ERROR_LOG("Failed to add selected workers in group!");
		return -1;
	}

	return 0;

}//close GetFreeWorkers()

//Get worker names
void WorkerManager::GetWorkerNames(std::vector<std::string>& worker_names){
			
	std::lock_guard<std::mutex> lock(m_mutex);
	worker_names.clear();
	worker_names= m_workers->get_device_list();

}//close GetWOrkerNames()

}//close namespace
