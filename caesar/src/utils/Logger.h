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
* @file Logger.h
* @class Logger
* @brief Logger
*
* Logger class for managing logging msg
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef Logger_h
#define Logger_h 1

#include <CodeUtils.h>

#include <TObject.h>
#include <TFITS.h>
#include <TMath.h>

#ifdef USE_TANGO
#include <tango.h>
#include <log4tango.h>
#endif

//== LOG4CXX HEADERS ==
#include <log4cxx/logger.h>
#include <log4cxx/xml/domconfigurator.h>
#include <log4cxx/simplelayout.h>
#include <log4cxx/patternlayout.h>
#include <log4cxx/consoleappender.h>
#include <log4cxx/fileappender.h>
#include <log4cxx/rollingfileappender.h>
#include <log4cxx/propertyconfigurator.h>
#include <log4cxx/net/syslogappender.h>
#include <log4cxx/helpers/exception.h>


/*
//Boost headers
#define BOOST_LOG_DYN_LINK 1 // necessary when linking the boost_log library dynamically
#include <boost/log/core/core.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/log/common.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/attributes.hpp>
#include <boost/log/sinks/sink.hpp>
#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sinks/syslog_backend.hpp>
#include <boost/log/trivial.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
namespace logging = boost::log;
namespace src = boost::log::sources;
namespace sinks = boost::log::sinks;
namespace keywords = boost::log::keywords;
*/


#include <syslog.h>



#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <stdlib.h>

#include <vector>
#include <algorithm>
#include <map>
#include <string>

using namespace std;


namespace Caesar {

/*
enum boost_severity_level {	
	boost_off= 0,
	boost_fatal= 1,
	boost_error= 2,
  boost_warn= 3,
	boost_info= 4,
  boost_debug= 5
};
BOOST_LOG_ATTRIBUTE_KEYWORD(severity, "Severity",boost_severity_level)
*/

class Logger : public TObject {

	public:
		Logger(std::string level="OFF",std::string tag="logger")
			: m_level(level), m_tag(tag)		
		{

		};
		virtual ~Logger(){};

	protected:
		#ifdef USE_TANGO
			virtual log4cxx::LevelPtr GetMappedLogLevel(int level) const {	
				if(level==log4tango::Level::DEBUG) return log4cxx::Level::getDebug();
				else if(level==log4tango::Level::INFO) return log4cxx::Level::getInfo();
				else if(level==log4tango::Level::WARN) return log4cxx::Level::getWarn();
				else if(level==log4tango::Level::ERROR) return log4cxx::Level::getError();
				else if(level==log4tango::Level::FATAL) return log4cxx::Level::getFatal();
				else return log4cxx::Level::getOff();
				return log4cxx::Level::getOff();
			}//close GetMappedLogLevel()
			
			/*
			virtual boost_severity_level GetMappedLogLevel(int level) const {	
				if(level==log4tango::Level::DEBUG) return boost_debug;
				else if(level==log4tango::Level::INFO) return boost_info;
				else if(level==log4tango::Level::WARN) return boost_warn;
				else if(level==log4tango::Level::ERROR) return boost_error;
				else if(level==log4tango::Level::FATAL) return boost_fatal;
				else return boost_off;
				return boost_off;
			}//close GetMappedLogLevel()
			*/
			virtual log4tango::Level::Value GetTangoLogLevelFromString(std::string sLevel) const {
  			if (sLevel == "DEBUG") return log4tango::Level::DEBUG;
				else if (sLevel == "INFO") return log4tango::Level::INFO; 
				else if (sLevel == "WARN") return log4tango::Level::WARN;
				else if (sLevel == "ERROR") return log4tango::Level::ERROR;
				else if (sLevel == "FATAL") return log4tango::Level::FATAL;
				else if (sLevel == "OFF") return log4tango::Level::OFF;
				else return log4tango::Level::OFF;
				return log4tango::Level::OFF;
			}	
			
		#endif

		/*
		boost_severity_level GetMappedLogLevelFromString(const std::string sLevel) const {
	 		if (sLevel == "DEBUG") return boost_debug;
			else if (sLevel == "INFO") return boost_info;
			else if (sLevel == "WARN") return boost_warn;
			else if (sLevel == "ERROR") return boost_error;
			else if (sLevel == "FATAL") return boost_fatal;
			else if (sLevel == "OFF") return boost_off;
			else return boost_off;
 			return boost_off;
		}//close GetMappedLogLevelFromString()
		*/

	public:
		virtual int Init() = 0;

		virtual std::string GetHost() const {
			char hostname[HOST_NAME_MAX];
			gethostname(hostname, HOST_NAME_MAX);
			std::string syslog_host(hostname); 
			return syslog_host;
		}

		
		virtual void Log(const std::string sLevel, const std::string msg, const std::string msg_prefix=""){
			log4cxx::LevelPtr logLevel= log4cxx::Level::toLevel(sLevel,log4cxx::Level::getOff());
			if(logger && logger->isEnabledFor(logLevel) ){
				log4cxx::helpers::MessageBuffer oss_; \
      	logger->forcedLog(logLevel, oss_.str(oss_ << msg_prefix << msg), log4cxx::spi::LocationInfo(m_tag.c_str(),__LOG4CXX_FUNC__,__LINE__) );
			}
		}	
		
		/*
		virtual void Log(const std::string sLevel, const std::string msg, const std::string msg_prefix=""){
			logging::core_ptr syslogger= logging::core::get();
			boost_severity_level syslevel= GetMappedLogLevelFromString(sLevel);
					
			if(syslogger && syslogger->get_logging_enabled() ){
				src::severity_logger_mt<boost_severity_level> lg(keywords::severity = boost_off);
				BOOST_LOG_SEV(lg, syslevel) << msg_prefix << msg;
			}
		}//close Log()
		*/

		virtual void SetLogLevel(const std::string sLevel){
			log4cxx::LevelPtr thisLogLevel= logger->getLevel();//retrieve the current log level
			if(!thisLogLevel) thisLogLevel= log4cxx::Level::getOff();
			logger->setLevel(log4cxx::Level::toLevel(sLevel,thisLogLevel));// In case of wrong string the Log level is set to the old or to OFF
		}//close SetLogLevel()

		#ifdef USE_TANGO
			virtual void Log(int level, const std::string msg, const std::string msg_prefix=""){
				log4cxx::LevelPtr logLevel= GetMappedLogLevel(level);					
				if(logger && logger->isEnabledFor(logLevel) ){
					log4cxx::helpers::MessageBuffer oss_; \
      		logger->forcedLog(logLevel, oss_.str(oss_ << msg_prefix << msg), log4cxx::spi::LocationInfo(m_tag.c_str(),__LOG4CXX_FUNC__,__LINE__) );
				}
			}//close Log()

			virtual void SetLogLevel(int level){
				log4cxx::LevelPtr logLevel= GetMappedLogLevel(level);
				logger->setLevel(logLevel);
			}//close SetLogLevel()
		#endif


	protected:
		std::string m_level;
		std::string m_tag;		
		log4cxx::LoggerPtr logger;
		log4cxx::LayoutPtr layout;
		log4cxx::AppenderPtr appender;		
		
	ClassDef(Logger,1)		

};//close Logger class


#ifdef __MAKECINT__
#pragma link C++ class Logger+;
#endif


class SysLogger : public Logger {

	public:
		SysLogger(const std::string level,const std::string tag,const std::string facility)
			: Logger(level,tag), m_facility(facility)
		{
			
		};

		virtual ~SysLogger(){};

	private:
		

	public:
		virtual int Init(){
			//Create logger
			logger= log4cxx::LoggerPtr(log4cxx::Logger::getLogger(m_tag.c_str()));
			if(!logger) return -1;
			return 0;

			//Define log layout
			layout= log4cxx::LayoutPtr( new log4cxx::PatternLayout("%c %m%n") );
			if(!layout) return -1;

			//Get host & facility code	
			std::string syslog_host= GetHost();
			int syslog_facility_code= GetFacilityCode(m_facility);

			//Create and add appender
			appender= log4cxx::AppenderPtr( new log4cxx::net::SyslogAppender(layout, syslog_host, syslog_facility_code) );
			appender->setOption("SYSLOGHOST",syslog_host);
			logger->addAppender(appender);
	
			//Set logging level (OFF is set if a wrong name is given!)
			log4cxx::LevelPtr level= log4cxx::Level::toLevel(m_level,log4cxx::Level::getOff());
			logger->setLevel(level);  
			
			return 0;
		}//close Init()

		
		virtual int GetFacilityCode(const std::string syslog_facility) const {
			int syslog_facility_code= LOG_LOCAL6;
			if(syslog_facility=="local0") syslog_facility_code= LOG_LOCAL0;
			else if(syslog_facility=="local1") syslog_facility_code= LOG_LOCAL1;
			else if(syslog_facility=="local2") syslog_facility_code= LOG_LOCAL2;
			else if(syslog_facility=="local3") syslog_facility_code= LOG_LOCAL3;
			else if(syslog_facility=="local4") syslog_facility_code= LOG_LOCAL4;
			else if(syslog_facility=="local5") syslog_facility_code= LOG_LOCAL5;
			else if(syslog_facility=="local6") syslog_facility_code= LOG_LOCAL6;
			else if(syslog_facility=="local7") syslog_facility_code= LOG_LOCAL7;
			else if(syslog_facility=="syslog") syslog_facility_code= LOG_SYSLOG;
			else if(syslog_facility=="user") syslog_facility_code= LOG_USER;
			else syslog_facility_code= LOG_LOCAL6;

			return syslog_facility_code;
		}//close GetFacilityCode()
		
	protected:
		std::string m_facility;
		
	ClassDef(SysLogger,1)		

};//close SysLogger class

#ifdef __MAKECINT__
#pragma link C++ class SysLogger+;
#endif


class FileLogger : public Logger {

	public:
		FileLogger(const std::string level,const std::string tag,const std::string filename,bool appendFlag,const std::string maxFileSize,int maxBackupFiles)
			: Logger(level,tag), m_filename(filename), m_appendFlag(appendFlag), m_maxFileSize(maxFileSize), m_maxBackupFiles(maxBackupFiles)
		{
			
		};

		virtual ~FileLogger(){};

	private:
		

	public:
		virtual int Init(){
			//Create logger
			logger= log4cxx::LoggerPtr(log4cxx::Logger::getLogger(m_tag));
			if(!logger) return -1;
			
			//Define log layout
			//layout= log4cxx::LayoutPtr( new log4cxx::PatternLayout("%d %-5p [%c] [%l] %m%n") );
			layout= log4cxx::LayoutPtr( new log4cxx::PatternLayout("%d %-5p [%c] %m%n") );
			if(!layout) return -1;

			//Create and add appender
			//appender= log4cxx::AppenderPtr( new log4cxx::FileAppender(layout,m_filename,m_appendFlag) );	
			appender= log4cxx::AppenderPtr( new log4cxx::RollingFileAppender(layout,m_filename,m_appendFlag) );	
			//appender->setMaxBackupIndex(m_maxBackupFiles);
			//appender->setMaxFileSize(m_maxFileSize);
			std::stringstream maxBackupFiles_stream;
			maxBackupFiles_stream << m_maxBackupFiles;
			appender->setOption("MAXBACKUPINDEX",maxBackupFiles_stream.str());//obsolete??
			appender->setOption("MAXFILESIZE",m_maxFileSize);
			appender->setOption("FILENAME",m_filename);
			logger->addAppender(appender);
	
			//Set logging level (OFF is set if a wrong name is given!)
			log4cxx::LevelPtr level= log4cxx::Level::toLevel(m_level,log4cxx::Level::getInfo());
			logger->setLevel(level);  

			return 0;
		}//close Init()

		
	protected:
		std::string m_filename;
		int m_maxBackupFiles; 
		bool m_appendFlag;
		std::string m_maxFileSize;

	ClassDef(FileLogger,1)		

};//close FileLogger class

#ifdef __MAKECINT__
#pragma link C++ class FileLogger+;
#endif


class ConsoleLogger : public Logger {

	public:
		ConsoleLogger(const std::string level,const std::string tag,const std::string target)
			: Logger(level,tag), m_target(target)
		{
			
		};

		virtual ~ConsoleLogger(){};

	private:
		

	public:
		virtual int Init(){
			//Create logger
			logger= log4cxx::LoggerPtr(log4cxx::Logger::getLogger(m_tag));
			if(!logger) return -1;
			
			//Define log layout
			//layout= log4cxx::LayoutPtr( new log4cxx::PatternLayout("%d %-5p [%c] [%l] %m%n") );
			layout= log4cxx::LayoutPtr( new log4cxx::PatternLayout("%d %-5p [%c] %m%n") );
			if(!layout) return -1;

			//Create and add appender
			appender= log4cxx::AppenderPtr( new log4cxx::ConsoleAppender(layout,m_target) );	
			appender->setOption("TARGET",m_target);
			logger->addAppender(appender);
	
			//Set logging level (OFF is set if a wrong name is given!)
			log4cxx::LevelPtr level= log4cxx::Level::toLevel(m_level,log4cxx::Level::getInfo());
			logger->setLevel(level);  

			return 0;
		}//close Init()

	
	protected:
		std::string m_target;
		

	ClassDef(ConsoleLogger,1)		

};//close ConsoleLogger class

#ifdef __MAKECINT__
#pragma link C++ class ConsoleLogger+;
#endif

enum LoggerTarget {
	eCONSOLE_TARGET= 1,
	eFILE_TARGET= 2,
	eSYSLOG_TARGET= 3
};

class LoggerManager : public TObject {
	
	public:
		static LoggerManager& Instance() {
    	// Since it's a static variable, if the class has already been created,
      // It won't be created again.
      // And it is thread-safe in C++11.
      static LoggerManager myInstance;
 
      // Return a reference to our instance.
      return myInstance;
    }
 
    // delete copy and move constructors and assign operators
    LoggerManager(LoggerManager const&) = delete;             // Copy construct
    LoggerManager(LoggerManager&&) = delete;                  // Move construct
    LoggerManager& operator=(LoggerManager const&) = delete;  // Copy assign
    LoggerManager& operator=(LoggerManager &&) = delete;      // Move assign

	public:
		enum LogTarget{eSysLog=1,eConsole=2,eFile=3};

  	static int CreateSysLogger(const std::string level,const std::string tag="syslogger",const std::string facility="local6") {
			//Skip if already created
			if(m_logger) {	
				cerr<<"CreateSysLogger(): WARN: Already created...skip!"<<endl;
				return -1;
			}

			//Create syslogger
			m_target= eSysLog;
			m_logger= new SysLogger(level,tag,facility);
			m_logger->Init();

			return 0;
		}//close CreateSysLogger()

		static int CreateFileLogger(const std::string level,const std::string tag="filelogger",const std::string filename="out.log",bool appendFlag=false,std::string maxFileSize="10MB",int maxBackupFiles=2) {
				
			//Skip if already created
			if(m_logger) {	
				return -1;
			}

			//Create syslogger
			m_target= eFile;
			m_logger= new FileLogger(level,tag,filename,appendFlag,maxFileSize,maxBackupFiles);
			m_logger->Init();
			
			return 0;
		}//close CreateFileLogger()

		static int CreateConsoleLogger(const std::string level,const std::string tag="consolelogger",const std::string target="System.out") {
				
			//Skip if already created
			if(m_logger) {	
				return -1;
			}

			//Create syslogger
			m_target= eFile;
			m_logger= new ConsoleLogger(level,tag,target);
			m_logger->Init();
			
			return 0;
		}//close CreateFileLogger()

		Logger* GetLogger(){return m_logger;}
		int GetLoggerTarget(){return m_target;}

		

	protected:
		LoggerManager(){
			
		};
		virtual ~LoggerManager(){
			if(m_logger){
				delete m_logger;
				m_logger= 0;
			}		
		};

	private:
		static int m_target;
		static Logger* m_logger;

	ClassDef(LoggerManager,1)	

};//close LoggerManager()

#ifdef __MAKECINT__
#pragma link C++ class LoggerManager+;
#endif

class ScopedLogger {
	public:
		//--> Constructor
		ScopedLogger(std::string level,std::string prefix="")
  		: m_level(level), m_msgprefix(prefix) 
		{}
			
		//--> Destructor
  	~ScopedLogger(){ 
			//logging command
			Logger* logger= LoggerManager::Instance().GetLogger();
			if(logger) logger->Log(m_level, m_sstream.str(), m_msgprefix);			
		}//close destructor

	public:
		//Get stream
  	std::stringstream& stream(){ 	
			return m_sstream; 
		}	
	
	private:
		std::stringstream m_sstream;
  	std::string m_level;
		std::string m_msgprefix;

};//close ScopedLogger


//== LOG MACROS ===
#define LOG_PREFIX \
	__CLASS_PREFIX__ + __FUNCTION__ + std::string("() - ")

#define LOG(Level, What) ScopedLogger(Level,LOG_PREFIX).stream() << What
#define INFO_LOG(What) LOG("INFO",What)
#define WARN_LOG(What) LOG("WARN",What)
#define DEBUG_LOG(What) LOG("DEBUG",What)
#define ERROR_LOG(What) LOG("ERROR",What)
#define FATAL_LOG(What) LOG("FATAL",What)

}//close namespace 

#endif

