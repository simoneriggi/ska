/*----- PROTECTED REGION ID(SFinderBrokerClass.h) ENABLED START -----*/
//=============================================================================
//
// file :        SFinderBrokerClass.h
//
// description : Include for the SFinderBroker root class.
//               This class is the singleton class for
//                the SFinderBroker device class.
//               It contains all properties and methods which the 
//               SFinderBroker requires only once e.g. the commands.
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


#ifndef SFinderBrokerClass_H
#define SFinderBrokerClass_H

#include <tango.h>
#include <SFinderBroker.h>


/*----- PROTECTED REGION END -----*/	//	SFinderBrokerClass.h


namespace SFinderBroker_ns
{
/*----- PROTECTED REGION ID(SFinderBrokerClass::classes for dynamic creation) ENABLED START -----*/


/*----- PROTECTED REGION END -----*/	//	SFinderBrokerClass::classes for dynamic creation

//=========================================
//	Define classes for commands
//=========================================
//	Command RegisterWorker class definition
class RegisterWorkerClass : public Tango::Command
{
public:
	RegisterWorkerClass(const char   *name,
	               Tango::CmdArgType in,
				   Tango::CmdArgType out,
				   const char        *in_desc,
				   const char        *out_desc,
				   Tango::DispLevel  level)
	:Command(name,in,out,in_desc,out_desc, level)	{};

	RegisterWorkerClass(const char   *name,
	               Tango::CmdArgType in,
				   Tango::CmdArgType out)
	:Command(name,in,out)	{};
	~RegisterWorkerClass() {};
	
	virtual CORBA::Any *execute (Tango::DeviceImpl *dev, const CORBA::Any &any);
	virtual bool is_allowed (Tango::DeviceImpl *dev, const CORBA::Any &any)
	{return (static_cast<SFinderBroker *>(dev))->is_RegisterWorker_allowed(any);}
};

//	Command ListWorkers class definition
class ListWorkersClass : public Tango::Command
{
public:
	ListWorkersClass(const char   *name,
	               Tango::CmdArgType in,
				   Tango::CmdArgType out,
				   const char        *in_desc,
				   const char        *out_desc,
				   Tango::DispLevel  level)
	:Command(name,in,out,in_desc,out_desc, level)	{};

	ListWorkersClass(const char   *name,
	               Tango::CmdArgType in,
				   Tango::CmdArgType out)
	:Command(name,in,out)	{};
	~ListWorkersClass() {};
	
	virtual CORBA::Any *execute (Tango::DeviceImpl *dev, const CORBA::Any &any);
	virtual bool is_allowed (Tango::DeviceImpl *dev, const CORBA::Any &any)
	{return (static_cast<SFinderBroker *>(dev))->is_ListWorkers_allowed(any);}
};

//	Command ListFreeWorkers class definition
class ListFreeWorkersClass : public Tango::Command
{
public:
	ListFreeWorkersClass(const char   *name,
	               Tango::CmdArgType in,
				   Tango::CmdArgType out,
				   const char        *in_desc,
				   const char        *out_desc,
				   Tango::DispLevel  level)
	:Command(name,in,out,in_desc,out_desc, level)	{};

	ListFreeWorkersClass(const char   *name,
	               Tango::CmdArgType in,
				   Tango::CmdArgType out)
	:Command(name,in,out)	{};
	~ListFreeWorkersClass() {};
	
	virtual CORBA::Any *execute (Tango::DeviceImpl *dev, const CORBA::Any &any);
	virtual bool is_allowed (Tango::DeviceImpl *dev, const CORBA::Any &any)
	{return (static_cast<SFinderBroker *>(dev))->is_ListFreeWorkers_allowed(any);}
};

//	Command ListBusyWorkers class definition
class ListBusyWorkersClass : public Tango::Command
{
public:
	ListBusyWorkersClass(const char   *name,
	               Tango::CmdArgType in,
				   Tango::CmdArgType out,
				   const char        *in_desc,
				   const char        *out_desc,
				   Tango::DispLevel  level)
	:Command(name,in,out,in_desc,out_desc, level)	{};

	ListBusyWorkersClass(const char   *name,
	               Tango::CmdArgType in,
				   Tango::CmdArgType out)
	:Command(name,in,out)	{};
	~ListBusyWorkersClass() {};
	
	virtual CORBA::Any *execute (Tango::DeviceImpl *dev, const CORBA::Any &any);
	virtual bool is_allowed (Tango::DeviceImpl *dev, const CORBA::Any &any)
	{return (static_cast<SFinderBroker *>(dev))->is_ListBusyWorkers_allowed(any);}
};

//	Command PingWorker class definition
class PingWorkerClass : public Tango::Command
{
public:
	PingWorkerClass(const char   *name,
	               Tango::CmdArgType in,
				   Tango::CmdArgType out,
				   const char        *in_desc,
				   const char        *out_desc,
				   Tango::DispLevel  level)
	:Command(name,in,out,in_desc,out_desc, level)	{};

	PingWorkerClass(const char   *name,
	               Tango::CmdArgType in,
				   Tango::CmdArgType out)
	:Command(name,in,out)	{};
	~PingWorkerClass() {};
	
	virtual CORBA::Any *execute (Tango::DeviceImpl *dev, const CORBA::Any &any);
	virtual bool is_allowed (Tango::DeviceImpl *dev, const CORBA::Any &any)
	{return (static_cast<SFinderBroker *>(dev))->is_PingWorker_allowed(any);}
};


/**
 *	The SFinderBrokerClass singleton definition
 */

#ifdef _TG_WINDOWS_
class __declspec(dllexport)  SFinderBrokerClass : public Tango::DeviceClass
#else
class SFinderBrokerClass : public Tango::DeviceClass
#endif
{
	/*----- PROTECTED REGION ID(SFinderBrokerClass::Additionnal DServer data members) ENABLED START -----*/
	
	
	/*----- PROTECTED REGION END -----*/	//	SFinderBrokerClass::Additionnal DServer data members

	public:
		//	write class properties data members
		Tango::DbData	cl_prop;
		Tango::DbData	cl_def_prop;
		Tango::DbData	dev_def_prop;
	
		//	Method prototypes
		static SFinderBrokerClass *init(const char *);
		static SFinderBrokerClass *instance();
		~SFinderBrokerClass();
		Tango::DbDatum	get_class_property(string &);
		Tango::DbDatum	get_default_device_property(string &);
		Tango::DbDatum	get_default_class_property(string &);
	
	protected:
		SFinderBrokerClass(string &);
		static SFinderBrokerClass *_instance;
		void command_factory();
		void attribute_factory(vector<Tango::Attr *> &);
		void pipe_factory();
		void write_class_property();
		void set_default_property();
		void get_class_property();
		string get_cvstag();
		string get_cvsroot();
	
	private:
		void device_factory(const Tango::DevVarStringArray *);
		void create_static_attribute_list(vector<Tango::Attr *> &);
		void erase_dynamic_attributes(const Tango::DevVarStringArray *,vector<Tango::Attr *> &);
		vector<string>	defaultAttList;
		Tango::Attr *get_attr_object_by_name(vector<Tango::Attr *> &att_list, string attname);
};

}	//	End of namespace

#endif   //	SFinderBroker_H
