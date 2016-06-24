/*----- PROTECTED REGION ID(SFinderBrokerClass.cpp) ENABLED START -----*/
static const char *RcsId      = "$Id:  $";
static const char *TagName    = "$Name:  $";
static const char *CvsPath    = "$Source:  $";
static const char *SvnPath    = "$HeadURL:  $";
static const char *HttpServer = "http://www.esrf.eu/computing/cs/tango/tango_doc/ds_doc/";
//=============================================================================
//
// file :        SFinderBrokerClass.cpp
//
// description : C++ source for the SFinderBrokerClass.
//               A singleton class derived from DeviceClass.
//               It implements the command and attribute list
//               and all properties and methods required
//               by the SFinderBroker once per process.
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


#include <SFinderBrokerClass.h>

/*----- PROTECTED REGION END -----*/	//	SFinderBrokerClass.cpp

//-------------------------------------------------------------------
/**
 *	Create SFinderBrokerClass singleton and
 *	return it in a C function for Python usage
 */
//-------------------------------------------------------------------
extern "C" {
#ifdef _TG_WINDOWS_

__declspec(dllexport)

#endif

	Tango::DeviceClass *_create_SFinderBroker_class(const char *name) {
		return SFinderBroker_ns::SFinderBrokerClass::init(name);
	}
}

namespace SFinderBroker_ns
{
//===================================================================
//	Initialize pointer for singleton pattern
//===================================================================
SFinderBrokerClass *SFinderBrokerClass::_instance = NULL;

//--------------------------------------------------------
/**
 * method : 		SFinderBrokerClass::SFinderBrokerClass(string &s)
 * description : 	constructor for the SFinderBrokerClass
 *
 * @param s	The class name
 */
//--------------------------------------------------------
SFinderBrokerClass::SFinderBrokerClass(string &s):Tango::DeviceClass(s)
{
	cout2 << "Entering SFinderBrokerClass constructor" << endl;
	set_default_property();
	write_class_property();

	/*----- PROTECTED REGION ID(SFinderBrokerClass::constructor) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinderBrokerClass::constructor

	cout2 << "Leaving SFinderBrokerClass constructor" << endl;
}

//--------------------------------------------------------
/**
 * method : 		SFinderBrokerClass::~SFinderBrokerClass()
 * description : 	destructor for the SFinderBrokerClass
 */
//--------------------------------------------------------
SFinderBrokerClass::~SFinderBrokerClass()
{
	/*----- PROTECTED REGION ID(SFinderBrokerClass::destructor) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinderBrokerClass::destructor

	_instance = NULL;
}


//--------------------------------------------------------
/**
 * method : 		SFinderBrokerClass::init
 * description : 	Create the object if not already done.
 *                  Otherwise, just return a pointer to the object
 *
 * @param	name	The class name
 */
//--------------------------------------------------------
SFinderBrokerClass *SFinderBrokerClass::init(const char *name)
{
	if (_instance == NULL)
	{
		try
		{
			string s(name);
			_instance = new SFinderBrokerClass(s);
		}
		catch (bad_alloc &)
		{
			throw;
		}
	}
	return _instance;
}

//--------------------------------------------------------
/**
 * method : 		SFinderBrokerClass::instance
 * description : 	Check if object already created,
 *                  and return a pointer to the object
 */
//--------------------------------------------------------
SFinderBrokerClass *SFinderBrokerClass::instance()
{
	if (_instance == NULL)
	{
		cerr << "Class is not initialised !!" << endl;
		exit(-1);
	}
	return _instance;
}



//===================================================================
//	Command execution method calls
//===================================================================
//--------------------------------------------------------
/**
 * method : 		RegisterWorkerClass::execute()
 * description : 	method to trigger the execution of the command.
 *
 * @param	device	The device on which the command must be executed
 * @param	in_any	The command input data
 *
 *	returns The command output data (packed in the Any object)
 */
//--------------------------------------------------------
CORBA::Any *RegisterWorkerClass::execute(Tango::DeviceImpl *device, const CORBA::Any &in_any)
{
	cout2 << "RegisterWorkerClass::execute(): arrived" << endl;
	Tango::DevString argin;
	extract(in_any, argin);
	return insert((static_cast<SFinderBroker *>(device))->register_worker(argin));
}

//--------------------------------------------------------
/**
 * method : 		ListWorkersClass::execute()
 * description : 	method to trigger the execution of the command.
 *
 * @param	device	The device on which the command must be executed
 * @param	in_any	The command input data
 *
 *	returns The command output data (packed in the Any object)
 */
//--------------------------------------------------------
CORBA::Any *ListWorkersClass::execute(Tango::DeviceImpl *device, TANGO_UNUSED(const CORBA::Any &in_any))
{
	cout2 << "ListWorkersClass::execute(): arrived" << endl;
	return insert((static_cast<SFinderBroker *>(device))->list_workers());
}

//--------------------------------------------------------
/**
 * method : 		ListFreeWorkersClass::execute()
 * description : 	method to trigger the execution of the command.
 *
 * @param	device	The device on which the command must be executed
 * @param	in_any	The command input data
 *
 *	returns The command output data (packed in the Any object)
 */
//--------------------------------------------------------
CORBA::Any *ListFreeWorkersClass::execute(Tango::DeviceImpl *device, TANGO_UNUSED(const CORBA::Any &in_any))
{
	cout2 << "ListFreeWorkersClass::execute(): arrived" << endl;
	return insert((static_cast<SFinderBroker *>(device))->list_free_workers());
}

//--------------------------------------------------------
/**
 * method : 		ListBusyWorkersClass::execute()
 * description : 	method to trigger the execution of the command.
 *
 * @param	device	The device on which the command must be executed
 * @param	in_any	The command input data
 *
 *	returns The command output data (packed in the Any object)
 */
//--------------------------------------------------------
CORBA::Any *ListBusyWorkersClass::execute(Tango::DeviceImpl *device, TANGO_UNUSED(const CORBA::Any &in_any))
{
	cout2 << "ListBusyWorkersClass::execute(): arrived" << endl;
	return insert((static_cast<SFinderBroker *>(device))->list_busy_workers());
}

//--------------------------------------------------------
/**
 * method : 		PingWorkerClass::execute()
 * description : 	method to trigger the execution of the command.
 *
 * @param	device	The device on which the command must be executed
 * @param	in_any	The command input data
 *
 *	returns The command output data (packed in the Any object)
 */
//--------------------------------------------------------
CORBA::Any *PingWorkerClass::execute(Tango::DeviceImpl *device, const CORBA::Any &in_any)
{
	cout2 << "PingWorkerClass::execute(): arrived" << endl;
	Tango::DevString argin;
	extract(in_any, argin);
	((static_cast<SFinderBroker *>(device))->ping_worker(argin));
	return new CORBA::Any();
}

//--------------------------------------------------------
/**
 * method : 		SubscribeWorkersClass::execute()
 * description : 	method to trigger the execution of the command.
 *
 * @param	device	The device on which the command must be executed
 * @param	in_any	The command input data
 *
 *	returns The command output data (packed in the Any object)
 */
//--------------------------------------------------------
CORBA::Any *SubscribeWorkersClass::execute(Tango::DeviceImpl *device, TANGO_UNUSED(const CORBA::Any &in_any))
{
	cout2 << "SubscribeWorkersClass::execute(): arrived" << endl;
	((static_cast<SFinderBroker *>(device))->subscribe_workers());
	return new CORBA::Any();
}

//--------------------------------------------------------
/**
 * method : 		SubmitSourceFinderJobClass::execute()
 * description : 	method to trigger the execution of the command.
 *
 * @param	device	The device on which the command must be executed
 * @param	in_any	The command input data
 *
 *	returns The command output data (packed in the Any object)
 */
//--------------------------------------------------------
CORBA::Any *SubmitSourceFinderJobClass::execute(Tango::DeviceImpl *device, const CORBA::Any &in_any)
{
	cout2 << "SubmitSourceFinderJobClass::execute(): arrived" << endl;
	const Tango::DevVarLongStringArray *argin;
	extract(in_any, argin);
	return insert((static_cast<SFinderBroker *>(device))->submit_source_finder_job(argin));
}


//===================================================================
//	Properties management
//===================================================================
//--------------------------------------------------------
/**
 *	Method      : SFinderBrokerClass::get_class_property()
 *	Description : Get the class property for specified name.
 */
//--------------------------------------------------------
Tango::DbDatum SFinderBrokerClass::get_class_property(string &prop_name)
{
	for (unsigned int i=0 ; i<cl_prop.size() ; i++)
		if (cl_prop[i].name == prop_name)
			return cl_prop[i];
	//	if not found, returns  an empty DbDatum
	return Tango::DbDatum(prop_name);
}

//--------------------------------------------------------
/**
 *	Method      : SFinderBrokerClass::get_default_device_property()
 *	Description : Return the default value for device property.
 */
//--------------------------------------------------------
Tango::DbDatum SFinderBrokerClass::get_default_device_property(string &prop_name)
{
	for (unsigned int i=0 ; i<dev_def_prop.size() ; i++)
		if (dev_def_prop[i].name == prop_name)
			return dev_def_prop[i];
	//	if not found, return  an empty DbDatum
	return Tango::DbDatum(prop_name);
}

//--------------------------------------------------------
/**
 *	Method      : SFinderBrokerClass::get_default_class_property()
 *	Description : Return the default value for class property.
 */
//--------------------------------------------------------
Tango::DbDatum SFinderBrokerClass::get_default_class_property(string &prop_name)
{
	for (unsigned int i=0 ; i<cl_def_prop.size() ; i++)
		if (cl_def_prop[i].name == prop_name)
			return cl_def_prop[i];
	//	if not found, return  an empty DbDatum
	return Tango::DbDatum(prop_name);
}


//--------------------------------------------------------
/**
 *	Method      : SFinderBrokerClass::set_default_property()
 *	Description : Set default property (class and device) for wizard.
 *                For each property, add to wizard property name and description.
 *                If default value has been set, add it to wizard property and
 *                store it in a DbDatum.
 */
//--------------------------------------------------------
void SFinderBrokerClass::set_default_property()
{
	string	prop_name;
	string	prop_desc;
	string	prop_def;
	vector<string>	vect_data;

	//	Set Default Class Properties

	//	Set Default device Properties
	prop_name = "federatedBrokers";
	prop_desc = "List of federated broker devices";
	prop_def  = "";
	vect_data.clear();
	if (prop_def.length()>0)
	{
		Tango::DbDatum	data(prop_name);
		data << vect_data ;
		dev_def_prop.push_back(data);
		add_wiz_dev_prop(prop_name, prop_desc,  prop_def);
	}
	else
		add_wiz_dev_prop(prop_name, prop_desc);
	prop_name = "maxNTasksPerWorker_default";
	prop_desc = "Maximum number of tasks per worker allowed by default";
	prop_def  = "10";
	vect_data.clear();
	vect_data.push_back("10");
	if (prop_def.length()>0)
	{
		Tango::DbDatum	data(prop_name);
		data << vect_data ;
		dev_def_prop.push_back(data);
		add_wiz_dev_prop(prop_name, prop_desc,  prop_def);
	}
	else
		add_wiz_dev_prop(prop_name, prop_desc);
}

//--------------------------------------------------------
/**
 *	Method      : SFinderBrokerClass::write_class_property()
 *	Description : Set class description fields as property in database
 */
//--------------------------------------------------------
void SFinderBrokerClass::write_class_property()
{
	//	First time, check if database used
	if (Tango::Util::_UseDb == false)
		return;

	Tango::DbData	data;
	string	classname = get_name();
	string	header;
	string::size_type	start, end;

	//	Put title
	Tango::DbDatum	title("ProjectTitle");
	string	str_title("");
	title << str_title;
	data.push_back(title);

	//	Put Description
	Tango::DbDatum	description("Description");
	vector<string>	str_desc;
	str_desc.push_back("");
	description << str_desc;
	data.push_back(description);

	//	put cvs or svn location
	string	filename("SFinderBroker");
	filename += "Class.cpp";

	// check for cvs information
	string	src_path(CvsPath);
	start = src_path.find("/");
	if (start!=string::npos)
	{
		end   = src_path.find(filename);
		if (end>start)
		{
			string	strloc = src_path.substr(start, end-start);
			//	Check if specific repository
			start = strloc.find("/cvsroot/");
			if (start!=string::npos && start>0)
			{
				string	repository = strloc.substr(0, start);
				if (repository.find("/segfs/")!=string::npos)
					strloc = "ESRF:" + strloc.substr(start, strloc.length()-start);
			}
			Tango::DbDatum	cvs_loc("cvs_location");
			cvs_loc << strloc;
			data.push_back(cvs_loc);
		}
	}

	// check for svn information
	else
	{
		string	src_path(SvnPath);
		start = src_path.find("://");
		if (start!=string::npos)
		{
			end = src_path.find(filename);
			if (end>start)
			{
				header = "$HeadURL: ";
				start = header.length();
				string	strloc = src_path.substr(start, (end-start));
				
				Tango::DbDatum	svn_loc("svn_location");
				svn_loc << strloc;
				data.push_back(svn_loc);
			}
		}
	}

	//	Get CVS or SVN revision tag
	
	// CVS tag
	string	tagname(TagName);
	header = "$Name: ";
	start = header.length();
	string	endstr(" $");
	
	end   = tagname.find(endstr);
	if (end!=string::npos && end>start)
	{
		string	strtag = tagname.substr(start, end-start);
		Tango::DbDatum	cvs_tag("cvs_tag");
		cvs_tag << strtag;
		data.push_back(cvs_tag);
	}
	
	// SVN tag
	string	svnpath(SvnPath);
	header = "$HeadURL: ";
	start = header.length();
	
	end   = svnpath.find(endstr);
	if (end!=string::npos && end>start)
	{
		string	strloc = svnpath.substr(start, end-start);
		
		string tagstr ("/tags/");
		start = strloc.find(tagstr);
		if ( start!=string::npos )
		{
			start = start + tagstr.length();
			end   = strloc.find(filename);
			string	strtag = strloc.substr(start, end-start-1);
			
			Tango::DbDatum	svn_tag("svn_tag");
			svn_tag << strtag;
			data.push_back(svn_tag);
		}
	}

	//	Get URL location
	string	httpServ(HttpServer);
	if (httpServ.length()>0)
	{
		Tango::DbDatum	db_doc_url("doc_url");
		db_doc_url << httpServ;
		data.push_back(db_doc_url);
	}

	//  Put inheritance
	Tango::DbDatum	inher_datum("InheritedFrom");
	vector<string> inheritance;
	inheritance.push_back("TANGO_BASE_CLASS");
	inher_datum << inheritance;
	data.push_back(inher_datum);

	//	Call database and and values
	get_db_class()->put_property(data);
}

//===================================================================
//	Factory methods
//===================================================================

//--------------------------------------------------------
/**
 *	Method      : SFinderBrokerClass::device_factory()
 *	Description : Create the device object(s)
 *                and store them in the device list
 */
//--------------------------------------------------------
void SFinderBrokerClass::device_factory(const Tango::DevVarStringArray *devlist_ptr)
{
	/*----- PROTECTED REGION ID(SFinderBrokerClass::device_factory_before) ENABLED START -----*/
	
	//	Add your own code
	
	/*----- PROTECTED REGION END -----*/	//	SFinderBrokerClass::device_factory_before

	//	Create devices and add it into the device list
	for (unsigned long i=0 ; i<devlist_ptr->length() ; i++)
	{
		cout4 << "Device name : " << (*devlist_ptr)[i].in() << endl;
		device_list.push_back(new SFinderBroker(this, (*devlist_ptr)[i]));
	}

	//	Manage dynamic attributes if any
	erase_dynamic_attributes(devlist_ptr, get_class_attr()->get_attr_list());

	//	Export devices to the outside world
	for (unsigned long i=1 ; i<=devlist_ptr->length() ; i++)
	{
		//	Add dynamic attributes if any
		SFinderBroker *dev = static_cast<SFinderBroker *>(device_list[device_list.size()-i]);
		dev->add_dynamic_attributes();

		//	Check before if database used.
		if ((Tango::Util::_UseDb == true) && (Tango::Util::_FileDb == false))
			export_device(dev);
		else
			export_device(dev, dev->get_name().c_str());
	}

	/*----- PROTECTED REGION ID(SFinderBrokerClass::device_factory_after) ENABLED START -----*/
	
	//	Add your own code
	
	/*----- PROTECTED REGION END -----*/	//	SFinderBrokerClass::device_factory_after
}
//--------------------------------------------------------
/**
 *	Method      : SFinderBrokerClass::attribute_factory()
 *	Description : Create the attribute object(s)
 *                and store them in the attribute list
 */
//--------------------------------------------------------
void SFinderBrokerClass::attribute_factory(vector<Tango::Attr *> &att_list)
{
	/*----- PROTECTED REGION ID(SFinderBrokerClass::attribute_factory_before) ENABLED START -----*/
	
	//	Add your own code
	
	/*----- PROTECTED REGION END -----*/	//	SFinderBrokerClass::attribute_factory_before

	//	Create a list of static attributes
	create_static_attribute_list(get_class_attr()->get_attr_list());
	/*----- PROTECTED REGION ID(SFinderBrokerClass::attribute_factory_after) ENABLED START -----*/
	
	//	Add your own code
	
	/*----- PROTECTED REGION END -----*/	//	SFinderBrokerClass::attribute_factory_after
}
//--------------------------------------------------------
/**
 *	Method      : SFinderBrokerClass::pipe_factory()
 *	Description : Create the pipe object(s)
 *                and store them in the pipe list
 */
//--------------------------------------------------------
void SFinderBrokerClass::pipe_factory()
{
	/*----- PROTECTED REGION ID(SFinderBrokerClass::pipe_factory_before) ENABLED START -----*/
	
	//	Add your own code
	
	/*----- PROTECTED REGION END -----*/	//	SFinderBrokerClass::pipe_factory_before
	/*----- PROTECTED REGION ID(SFinderBrokerClass::pipe_factory_after) ENABLED START -----*/
	
	//	Add your own code
	
	/*----- PROTECTED REGION END -----*/	//	SFinderBrokerClass::pipe_factory_after
}
//--------------------------------------------------------
/**
 *	Method      : SFinderBrokerClass::command_factory()
 *	Description : Create the command object(s)
 *                and store them in the command list
 */
//--------------------------------------------------------
void SFinderBrokerClass::command_factory()
{
	/*----- PROTECTED REGION ID(SFinderBrokerClass::command_factory_before) ENABLED START -----*/
	
	//	Add your own code
	
	/*----- PROTECTED REGION END -----*/	//	SFinderBrokerClass::command_factory_before


	//	Command RegisterWorker
	RegisterWorkerClass	*pRegisterWorkerCmd =
		new RegisterWorkerClass("RegisterWorker",
			Tango::DEV_STRING, Tango::DEVVAR_LONGSTRINGARRAY,
			"Worker device name to be registered",
			"",
			Tango::OPERATOR);
	command_list.push_back(pRegisterWorkerCmd);

	//	Command ListWorkers
	ListWorkersClass	*pListWorkersCmd =
		new ListWorkersClass("ListWorkers",
			Tango::DEV_VOID, Tango::DEVVAR_STRINGARRAY,
			"",
			"List of worker devices registered",
			Tango::OPERATOR);
	command_list.push_back(pListWorkersCmd);

	//	Command ListFreeWorkers
	ListFreeWorkersClass	*pListFreeWorkersCmd =
		new ListFreeWorkersClass("ListFreeWorkers",
			Tango::DEV_VOID, Tango::DEVVAR_STRINGARRAY,
			"",
			"List of free workers",
			Tango::OPERATOR);
	command_list.push_back(pListFreeWorkersCmd);

	//	Command ListBusyWorkers
	ListBusyWorkersClass	*pListBusyWorkersCmd =
		new ListBusyWorkersClass("ListBusyWorkers",
			Tango::DEV_VOID, Tango::DEVVAR_STRINGARRAY,
			"",
			"List busy workers",
			Tango::OPERATOR);
	command_list.push_back(pListBusyWorkersCmd);

	//	Command PingWorker
	PingWorkerClass	*pPingWorkerCmd =
		new PingWorkerClass("PingWorker",
			Tango::DEV_STRING, Tango::DEV_VOID,
			"Worker name",
			"",
			Tango::OPERATOR);
	command_list.push_back(pPingWorkerCmd);

	//	Command SubscribeWorkers
	SubscribeWorkersClass	*pSubscribeWorkersCmd =
		new SubscribeWorkersClass("SubscribeWorkers",
			Tango::DEV_VOID, Tango::DEV_VOID,
			"",
			"",
			Tango::OPERATOR);
	command_list.push_back(pSubscribeWorkersCmd);

	//	Command SubmitSourceFinderJob
	SubmitSourceFinderJobClass	*pSubmitSourceFinderJobCmd =
		new SubmitSourceFinderJobClass("SubmitSourceFinderJob",
			Tango::DEVVAR_LONGSTRINGARRAY, Tango::DEVVAR_LONGSTRINGARRAY,
			"String arg\n[0]: input image filename\n[1]: config options\n\nLong arg\n[0]: Max number of workers to be allocated",
			"Long arg\n[0]: ack code\n\nString arg\n[0]: err description",
			Tango::OPERATOR);
	command_list.push_back(pSubmitSourceFinderJobCmd);

	/*----- PROTECTED REGION ID(SFinderBrokerClass::command_factory_after) ENABLED START -----*/
	
	//	Add your own code
	
	/*----- PROTECTED REGION END -----*/	//	SFinderBrokerClass::command_factory_after
}

//===================================================================
//	Dynamic attributes related methods
//===================================================================

//--------------------------------------------------------
/**
 * method : 		SFinderBrokerClass::create_static_attribute_list
 * description : 	Create the a list of static attributes
 *
 * @param	att_list	the ceated attribute list
 */
//--------------------------------------------------------
void SFinderBrokerClass::create_static_attribute_list(vector<Tango::Attr *> &att_list)
{
	for (unsigned long i=0 ; i<att_list.size() ; i++)
	{
		string att_name(att_list[i]->get_name());
		transform(att_name.begin(), att_name.end(), att_name.begin(), ::tolower);
		defaultAttList.push_back(att_name);
	}

	cout2 << defaultAttList.size() << " attributes in default list" << endl;

	/*----- PROTECTED REGION ID(SFinderBrokerClass::create_static_att_list) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinderBrokerClass::create_static_att_list
}


//--------------------------------------------------------
/**
 * method : 		SFinderBrokerClass::erase_dynamic_attributes
 * description : 	delete the dynamic attributes if any.
 *
 * @param	devlist_ptr	the device list pointer
 * @param	list of all attributes
 */
//--------------------------------------------------------
void SFinderBrokerClass::erase_dynamic_attributes(const Tango::DevVarStringArray *devlist_ptr, vector<Tango::Attr *> &att_list)
{
	Tango::Util *tg = Tango::Util::instance();

	for (unsigned long i=0 ; i<devlist_ptr->length() ; i++)
	{
		Tango::DeviceImpl *dev_impl = tg->get_device_by_name(((string)(*devlist_ptr)[i]).c_str());
		SFinderBroker *dev = static_cast<SFinderBroker *> (dev_impl);

		vector<Tango::Attribute *> &dev_att_list = dev->get_device_attr()->get_attribute_list();
		vector<Tango::Attribute *>::iterator ite_att;
		for (ite_att=dev_att_list.begin() ; ite_att != dev_att_list.end() ; ++ite_att)
		{
			string att_name((*ite_att)->get_name_lower());
			if ((att_name == "state") || (att_name == "status"))
				continue;
			vector<string>::iterator ite_str = find(defaultAttList.begin(), defaultAttList.end(), att_name);
			if (ite_str == defaultAttList.end())
			{
				cout2 << att_name << " is a UNWANTED dynamic attribute for device " << (*devlist_ptr)[i] << endl;
				Tango::Attribute &att = dev->get_device_attr()->get_attr_by_name(att_name.c_str());
				dev->remove_attribute(att_list[att.get_attr_idx()], true, false);
				--ite_att;
			}
		}
	}
	/*----- PROTECTED REGION ID(SFinderBrokerClass::erase_dynamic_attributes) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinderBrokerClass::erase_dynamic_attributes
}

//--------------------------------------------------------
/**
 *	Method      : SFinderBrokerClass::get_attr_by_name()
 *	Description : returns Tango::Attr * object found by name
 */
//--------------------------------------------------------
Tango::Attr *SFinderBrokerClass::get_attr_object_by_name(vector<Tango::Attr *> &att_list, string attname)
{
	vector<Tango::Attr *>::iterator it;
	for (it=att_list.begin() ; it<att_list.end() ; ++it)
		if ((*it)->get_name()==attname)
			return (*it);
	//	Attr does not exist
	return NULL;
}


/*----- PROTECTED REGION ID(SFinderBrokerClass::Additional Methods) ENABLED START -----*/

/*----- PROTECTED REGION END -----*/	//	SFinderBrokerClass::Additional Methods
} //	namespace
