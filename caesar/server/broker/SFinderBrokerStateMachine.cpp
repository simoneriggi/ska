/*----- PROTECTED REGION ID(SFinderBrokerStateMachine.cpp) ENABLED START -----*/
static const char *RcsId = "$Id:  $";
//=============================================================================
//
// file :        SFinderBrokerStateMachine.cpp
//
// description : State machine file for the SFinderBroker class
//
// project :     
//
// This file is part of Tango device class.
// 
// Tango is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// Tango is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Tango.  If not, see <http://www.gnu.org/licenses/>.
// 
// $Author:  $
//
// $Revision:  $
// $Date:  $
//
// $HeadURL:  $
//
//=============================================================================
//                This file is generated by POGO
//        (Program Obviously used to Generate tango Object)
//=============================================================================

#include <SFinderBroker.h>

/*----- PROTECTED REGION END -----*/	//	SFinderBroker::SFinderBrokerStateMachine.cpp

//================================================================
//  States  |  Description
//================================================================


namespace SFinderBroker_ns
{
//=================================================
//		Attributes Allowed Methods
//=================================================

//--------------------------------------------------------
/**
 *	Method      : SFinderBroker::is_dynStringAttr_allowed()
 *	Description : Execution allowed for dynStringAttr attribute
 */
//--------------------------------------------------------
bool SFinderBroker::is_dynStringAttr_allowed(TANGO_UNUSED(Tango::AttReqType type))
{

	//	Not any excluded states for dynStringAttr attribute in read access.
	/*----- PROTECTED REGION ID(SFinderBroker::dynStringAttrStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinderBroker::dynStringAttrStateAllowed_READ
	return true;
}


//=================================================
//		Commands Allowed Methods
//=================================================

//--------------------------------------------------------
/**
 *	Method      : SFinderBroker::is_RegisterWorker_allowed()
 *	Description : Execution allowed for RegisterWorker attribute
 */
//--------------------------------------------------------
bool SFinderBroker::is_RegisterWorker_allowed(TANGO_UNUSED(const CORBA::Any &any))
{
	//	Not any excluded states for RegisterWorker command.
	/*----- PROTECTED REGION ID(SFinderBroker::RegisterWorkerStateAllowed) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinderBroker::RegisterWorkerStateAllowed
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinderBroker::is_ListWorkers_allowed()
 *	Description : Execution allowed for ListWorkers attribute
 */
//--------------------------------------------------------
bool SFinderBroker::is_ListWorkers_allowed(TANGO_UNUSED(const CORBA::Any &any))
{
	//	Not any excluded states for ListWorkers command.
	/*----- PROTECTED REGION ID(SFinderBroker::ListWorkersStateAllowed) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinderBroker::ListWorkersStateAllowed
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinderBroker::is_ListFreeWorkers_allowed()
 *	Description : Execution allowed for ListFreeWorkers attribute
 */
//--------------------------------------------------------
bool SFinderBroker::is_ListFreeWorkers_allowed(TANGO_UNUSED(const CORBA::Any &any))
{
	//	Not any excluded states for ListFreeWorkers command.
	/*----- PROTECTED REGION ID(SFinderBroker::ListFreeWorkersStateAllowed) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinderBroker::ListFreeWorkersStateAllowed
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinderBroker::is_ListBusyWorkers_allowed()
 *	Description : Execution allowed for ListBusyWorkers attribute
 */
//--------------------------------------------------------
bool SFinderBroker::is_ListBusyWorkers_allowed(TANGO_UNUSED(const CORBA::Any &any))
{
	//	Not any excluded states for ListBusyWorkers command.
	/*----- PROTECTED REGION ID(SFinderBroker::ListBusyWorkersStateAllowed) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinderBroker::ListBusyWorkersStateAllowed
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinderBroker::is_PingWorker_allowed()
 *	Description : Execution allowed for PingWorker attribute
 */
//--------------------------------------------------------
bool SFinderBroker::is_PingWorker_allowed(TANGO_UNUSED(const CORBA::Any &any))
{
	//	Not any excluded states for PingWorker command.
	/*----- PROTECTED REGION ID(SFinderBroker::PingWorkerStateAllowed) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinderBroker::PingWorkerStateAllowed
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinderBroker::is_SubscribeWorkers_allowed()
 *	Description : Execution allowed for SubscribeWorkers attribute
 */
//--------------------------------------------------------
bool SFinderBroker::is_SubscribeWorkers_allowed(TANGO_UNUSED(const CORBA::Any &any))
{
	//	Not any excluded states for SubscribeWorkers command.
	/*----- PROTECTED REGION ID(SFinderBroker::SubscribeWorkersStateAllowed) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinderBroker::SubscribeWorkersStateAllowed
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinderBroker::is_SubmitSourceFinderJob_allowed()
 *	Description : Execution allowed for SubmitSourceFinderJob attribute
 */
//--------------------------------------------------------
bool SFinderBroker::is_SubmitSourceFinderJob_allowed(TANGO_UNUSED(const CORBA::Any &any))
{
	//	Not any excluded states for SubmitSourceFinderJob command.
	/*----- PROTECTED REGION ID(SFinderBroker::SubmitSourceFinderJobStateAllowed) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinderBroker::SubmitSourceFinderJobStateAllowed
	return true;
}


/*----- PROTECTED REGION ID(SFinderBroker::SFinderBrokerStateAllowed.AdditionalMethods) ENABLED START -----*/

//	Additional Methods

/*----- PROTECTED REGION END -----*/	//	SFinderBroker::SFinderBrokerStateAllowed.AdditionalMethods

}	//	End of namespace
