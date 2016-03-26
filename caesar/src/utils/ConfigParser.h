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
* @file ConfigParser.h
* @class ConfigParser
* @brief Parse the configuration file containing program parameters
* 
* @author S. Riggi
* @date 25/04/2010
*/

#ifndef ConfigParser_h
#define ConfigParser_h 1

#include <TObject.h>
#include <TBranch.h>
#include <TTree.h>

#include <vector>
#include <string>

#include <Img.h>

namespace Caesar{



class OptionBase : public TObject {
	public:
    virtual ~OptionBase() {}
	public:
		virtual void Print() = 0;
		virtual int SetValueFromStream(std::stringstream& sstream) = 0;
		virtual int SetValueFromString(std::string const& s) = 0;
		virtual int AddBranch(TTree* tree) = 0;
	ClassDef(OptionBase,1)	
};
typedef OptionBase* OptionBasePtr;
//typedef std::shared_ptr<OptionBase> OptionBasePtr;
typedef std::map<std::string, OptionBasePtr> OptionMap;

template <typename T> 
class Option : public OptionBase {

	public:
		/** 
		\brief Default Class constructor (needed by ROOT IO)
 		*/
		Option(){};
		/** 
		\brief Class constructor
 		*/
  	Option(std::string name,T defaultValue,T minValue,T maxValue) : 
			m_name(name), m_defaultValue(defaultValue), m_minValue(minValue), m_maxValue(maxValue)
		{
			m_value= m_defaultValue;
			m_hasSetValue= false;
		};
		/** 
		\brief Class destructor
 		*/
  	virtual ~Option(){};

		typedef T OptionType;

	public:
		/** 
		\brief Set option value from stringstream
 		*/
		int SetValueFromStream(std::stringstream& sstream){
			T tmpValue;
			if(typeid(T) == typeid(bool)){
				sstream >> std::boolalpha >> tmpValue;	
			}
			else{
				sstream>>tmpValue;
			}
			
			if(CheckRange() && (tmpValue<m_minValue || tmpValue>m_maxValue) ) {
				return -1;
			}
			m_value= tmpValue;
			m_hasSetValue= true;
			return 0;
		}
		/** 
		\brief Set option value from string
 		*/
		int SetValueFromString(std::string const& s){
			std::stringstream sstream(s);
			return SetValueFromStream(sstream);
		}
		/** 
		\brief Set option value
 		*/
		int SetValue(T value){
			if(CheckRange() && (value<m_minValue || value>m_maxValue) ) return -1;
			m_value= value;
			m_hasSetValue= true;
			return 0;
		}
		/** 
		\brief Get option value
 		*/
		T GetValue(){	
			if(!m_hasSetValue) return GetDefaultValue();//return default if no value has been set
			return m_value;
		}
		/** 
		\brief Get default value
 		*/	
		T GetDefaultValue(){
			return m_defaultValue;
		}

		/** 
		\brief Check given value against registered range (non sense for bool, string/char)?
 		*/
		bool CheckRange(){
			if(typeid(T) == typeid(std::string)){
				return false;
			}
			if( (typeid(T)==typeid(char)) || (typeid(T)==typeid(char*)) ) {
				return false;
			}
			if(typeid(T) == typeid(bool)){
				return false;
			}
			return true;
		}//close CheckRange()


		/** 
		\brief Print option
 		*/
		void Print(){
			cout<<"Option: {\"name\": \""<<m_name<<"\", \"value\": "<<m_value<<", \"default\": "<<m_defaultValue<<", \"min\": "<<m_minValue<<" \"max\": "<<m_maxValue<<"}"<<endl;
		}

		/** 
		\brief Add a TBranch to input TTree
 		*/
		int AddBranch(TTree* tree){
			if(!tree) {
				cerr<<"Option::AddBranch(): WARN: Null tree ptr given!"<<endl;			
				return -1;
			}
			TBranch* branch= tree->Branch(m_name.c_str(),this);
			if(!branch) {
				cerr<<"Option::AddBranch(): WARN: Failed to create a branch for this option..."<<endl;
				return -1;
			}
			return 0;
		}

	private:
		
		std::string m_name;
		T m_value;
		T m_defaultValue;
		T m_minValue;
		T m_maxValue;
		bool m_hasSetValue;

	ClassDef(Option,1)	
};



class OptionFactory : public TObject {

	public:
		static OptionFactory& Instance() {
    	// Since it's a static variable, if the class has already been created,
      // It won't be created again.
      // And it is thread-safe in C++11.
      static OptionFactory myInstance;
 
      // Return a reference to our instance.
      return myInstance;
    }
 
    // delete copy and move constructors and assign operators
    OptionFactory(OptionFactory const&) = delete;             // Copy construct
    OptionFactory(OptionFactory&&) = delete;                  // Move construct
    OptionFactory& operator=(OptionFactory const&) = delete;  // Copy assign
    OptionFactory& operator=(OptionFactory &&) = delete;      // Move assign


	public:

		/** 
		\brief Create an option
 		*/
		template<typename T>
		static OptionBasePtr Create(std::string name,T defaultValue,T minValue,T maxValue) {	
			return OptionBasePtr(new Option<T>(name,defaultValue,minValue,maxValue));
		}//close Create()

		/** 
		\brief Register option
 		*/
		template <typename T>
  		int Register(std::string name,T defaultValue,T minValue,T maxValue) {
				//Check args
				if(name=="") {
					cerr<<"OptionFactory()::Register(): WARN: Invalid option name given!"<<endl;					
					return -1;
				}
				if(minValue>maxValue) {
					cerr<<"OptionFactory()::Register(): WARN: Invalid range parameter given!"<<endl;					
					return -1;
				}

				//Check if option with given name already exists
				OptionBasePtr thisOption= GetOptionBase(name);
				if(thisOption){
					cerr<<"OptionFactory()::Register(): WARN: Option "<<name<<" already registered, skip it!"<<endl;
					return 0;
				}
	
				//Option does not exist, create one and add to the map
				OptionBasePtr aNewOption= Create<T>(name,defaultValue,minValue,maxValue);
				m_RegisteredOptions.insert( std::pair<std::string,OptionBasePtr>(name,aNewOption) );
		
				return 0;
  		}//close RegisterOption()
		
		/** 
		\brief Has option?
 		*/
		bool HasOption(std::string name){
			OptionMap::iterator it= m_RegisteredOptions.find(name);
			if ( m_RegisteredOptions.empty() || it==m_RegisteredOptions.end() )
				return false;
			return true;
		}

		/** 
		\brief Get option base pointer
 		*/
		OptionBasePtr GetOptionBase(std::string name){
			OptionMap::iterator it= m_RegisteredOptions.find(name);
			if ( m_RegisteredOptions.empty() || it==m_RegisteredOptions.end() )
				return 0;
			return it->second;
		}

		/** 
		\brief Get option impl pointer
 		*/
		template <typename T>
		//std::shared_ptr<Option<T>> GetOption(std::string name){
		Option<T>* GetOption(std::string name){
			//Try to dynamical cast to template type
			OptionBasePtr thisOption= GetOptionBase(name);
			if(!thisOption){
				cerr<<"OptionFactory()::GetOption(): INFO: Option not found!"<<endl;
				return 0;
			}
			//std::shared_ptr<Option<T>> thisCastedOption= std::dynamic_pointer_cast<Option<T>>(thisOption);
			Option<T>* thisCastedOption= dynamic_cast<Option<T>*>(thisOption);
			
			if(!thisCastedOption){
				cerr<<"OptionFactory()::GetOption(): ERROR: Failed to cast option to given data type (option registered with different type?)!"<<endl;
				return 0;
			}
			return thisCastedOption;
		}

		/** 
		\brief Get option value
 		*/
		template <typename T>
			int GetOptionValue(std::string name,T& value){
				Option<T>* thisOption= GetOption<T>(name);
				if(!thisOption) return -1;
				value= thisOption->GetValue();
				return 0;
			}//close GetOptionValue()

		/** 
		\brief Set option value
 		*/
		template <typename T>
			int SetOption(std::string name,T value){
				//std::shared_ptr<Option<T>> thisOption= GetOption<T>(name);	
				Option<T>* thisOption= GetOption<T>(name);
				if(!thisOption){
					cerr<<"OptionFactory()::SetOption(): ERROR: Cannot set value ("<<value<<") in option "<<name<<" (option not registered)!"<<endl;	
					return -1;
				}
				thisOption->SetValue(value);
				return 0;
			}//close SetOption()

		/** 
		\brief Set option value from string
 		*/
		int SetOptionFromString(std::string name,std::string stringified_value){
				OptionBasePtr thisOption= GetOptionBase(name);
				if(!thisOption){
					cerr<<"OptionFactory()::SetOptionFromString(): ERROR: Cannot set value ("<<stringified_value<<") in option "<<name<<" (option not registered)!"<<endl;	
					return -1;
				}
				thisOption->SetValueFromString(stringified_value);
				return 0;
			}//close SetOption()
		

		/** 
		\brief Print registered options
 		*/
		void PrintOptions(){
			for (OptionMap::const_iterator it = m_RegisteredOptions.begin(); it!=m_RegisteredOptions.end(); it++){
    		(it->second)->Print();
			}
		}//close PrintOptions()

		/** 
		\brief Make a TTree with all options
 		*/	
		TTree* MakeTTree(std::string treeName="CfgInfo"){

			//Create tree
			TTree* tree= new TTree(treeName.c_str(),treeName.c_str());

			//Adding branches
			for (OptionMap::const_iterator it = m_RegisteredOptions.begin(); it!=m_RegisteredOptions.end(); it++){
    		if( (it->second)->AddBranch(tree)<0 ) {
					cerr<<"OptionFactory::MakeTTree(): WARN: Failed to add branch for current option... "<<endl;
					continue;
				}
			}

			//Fill tree
			tree->Fill();

			return tree;
		}//close MakeTTree()
		
	protected:
		OptionFactory(){};
		virtual ~OptionFactory(){
			//Deleting map with options
			cout<<"~OptionFactory(): INFO: Deleting registered options..."<<endl;
			for (OptionMap::const_iterator it = m_RegisteredOptions.begin(); it!=m_RegisteredOptions.end(); it++){
				if(it->second) delete it->second;
			}//end loop map
			m_RegisteredOptions.clear();
		};

	private:
		OptionMap m_RegisteredOptions;		
	
	ClassDef(OptionFactory,1)	
};



#define REGISTER_OPTION(name,type,default_value,min_value,max_value) (OptionFactory::Instance().Register<type>(#name,default_value,min_value,max_value))
#define GET_OPTION(name,type) (OptionFactory::Instance().GetOption<type>(#name))
#define SET_OPTION(name,type,value) (OptionFactory::Instance().SetOption<type>(#name,value))
#define PRINT_OPTIONS() (OptionFactory::Instance().PrintOptions())
//#define GET_OPTION_VALUE(name,value) (OptionFactory::Instance().GetOptionValue<typeid(value)>(#name,value))
#define GET_OPTION_VALUE(name,value) (OptionFactory::Instance().GetOptionValue(#name,value))



class ConfigParser : public TObject {
  
	public:
		static ConfigParser& Instance() {
    	// Since it's a static variable, if the class has already been created,
      // It won't be created again.
      // And it is thread-safe in C++11.
      static ConfigParser myInstance;
 
			if(!m_HasRegisteredOptions){
				RegisterPredefinedOptions();
			}

      // Return a reference to our instance.
      return myInstance;
    }
 
    // delete copy and move constructors and assign operators
    ConfigParser(ConfigParser const&) = delete;             // Copy construct
    ConfigParser(ConfigParser&&) = delete;                  // Move construct
    ConfigParser& operator=(ConfigParser const&) = delete;  // Copy assign
    ConfigParser& operator=(ConfigParser &&) = delete;      // Move assign

	protected:
		/** 
		\brief Class constructor
 		*/
  	ConfigParser();
		/** 
		\brief Class destructor
 		*/
  	virtual ~ConfigParser();

	public:
		/** 
		\brief Read the config file, parse and set info to be used by other classes
 		*/
		int Parse(std::string filename);
		
		/** 
		\brief Print registered options
 		*/
		void PrintOptions(){
			return OptionFactory::Instance().PrintOptions();
		}

		/** 
		\brief Register options
 		*/
		template <typename T>
			std::shared_ptr<Option<T>> GetOption(std::string name){
				return OptionFactory::Instance().GetOption<T>(name);
			}//close GetOption()
	
		/** 
		\brief Check registered options
 		*/
		bool HasOption(std::string name){//Find option name
			return OptionFactory::Instance().HasOption(name);
		}
		/** 
		\brief Set option from string
 		*/
		int SetOptionFromString(std::string name,std::string stringified_value){
			return OptionFactory::Instance().SetOptionFromString(name,stringified_value);
		}

		/** 
		\brief Get option value
 		*/
		template <typename T>
			int GetOptionValue(std::string name,T& value){
				return OptionFactory::Instance().GetOptionValue(name,value);
			}//close GetOptionValue()

		/** 
		\brief Get config option tree
 		*/
		TTree* GetConfigTree(std::string treeName="CfgInfo"){
			return OptionFactory::Instance().MakeTTree(treeName);
		}
	
	private:
		static int RegisterPredefinedOptions();
		
	private:

		std::string m_ConfigFile;
		static bool m_HasRegisteredOptions;
		
	ClassDef(ConfigParser,1)
		
};


#ifdef __MAKECINT__
#pragma link C++ class OptionBase+;
#pragma link C++ class OptionBase*+;
//#pragma link C++ class OptionBasePtr+;
#pragma link C++ class std::map<std::string,Caesar::OptionBase*>+;
//#pragma link C++ typedef OptionMap+;
#pragma link C++ class Option<int>+;
#pragma link C++ class Option<long int>+;
#pragma link C++ class Option<float>+;
#pragma link C++ class Option<double>+;
#pragma link C++ class Option<std::string>+;
#pragma link C++ class Option<bool>+;
#pragma link C++ class OptionFactory+;
#pragma link C++ class ConfigParser+;
#pragma link C++ MACRO REGISTER_OPTION;
#endif


}//close namespace 

#endif
 
