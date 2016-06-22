#ifndef SFinderThread_H
#define SFinderThread_H

#include <SFinder.h>

#include <tango.h>

#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>

#include <Img.h>
#include <BkgData.h>
#include <Source.h>
#include <WorkerData.h>
namespace Caesar {
	class Img;
	class BkgData;
	class Source;
	class WorkerData;
}


namespace SFinder_ns {
  
class SFinder;

//class SFinderThread : public omni_thread, public Tango::LogAdapter {
class SFinderThread : public Tango::LogAdapter {
  
	private:
  	static SFinder* device;	
			
  public :
    SFinderThread(SFinder* dev);
    ~SFinderThread();

	public:
		/** 
		\brief Start the thread
 		*/		
		void Start(std::string filename,std::string runId,std::vector<long int> tileMinX,std::vector<long int> tileMaxX,std::vector<long int> tileMinY,std::vector<long int> tileMaxY){
			m_stopThread = false;
    	m_thread = std::thread(&SFinderThread::Run,this,filename,runId,tileMinX,tileMaxX,tileMinY,tileMaxY);
			//if(m_thread.joinable()) m_thread.join();
    }
		
		/** 
		\brief Stop
 		*/
		void Stop(){
			DEBUG_STREAM<<"SchedulerThread::Stop(): INFO: Called Stop()..."<<endl;
			std::lock_guard<std::mutex> lock( m_mutex );
			m_stopThread = true;
			if(m_thread.joinable()) m_thread.join();
			DEBUG_STREAM<<"SchedulerThread::Stop(): INFO: done!"<<endl;
		}
     		
	private:
		/** 
		\brief Main thread function 
 		*/
		void Run(std::string filename,std::string runId,std::vector<long int> tileMinX,std::vector<long int> tileMaxX,std::vector<long int> tileMinY,std::vector<long int> tileMaxY);
		/** 
		\brief Run source finder task for a single image
 		*/
		int RunTask(std::vector<Caesar::Source*>& sources,const std::string& filename,const std::string& runId,long int tileMinX=-1,long int tileMaxX=-1,long int tileMinY=-1,long int tileMaxY=-1);
		/** 
		\brief Read image
 		*/
    Caesar::Img* ReadImage(const std::string& filename,long int tileMinX=-1,long int tileMaxX=-1,long int tileMinY=-1,long int tileMaxY=-1);
 		/** 
		\brief Compute stats & bkg
 		*/
		Caesar::BkgData* ComputeStatsAndBkg(Caesar::Img* img);
		/** 
		\brief Find compact sources
 		*/
		//int FindCompactSources(std::vector<Caesar::Source*>& sources,Caesar::Img* inputImg,bool computeStatsAndBkg=true,Caesar::BkgData* inputBkgData=0);
		int FindCompactSources(Caesar::WorkerData& workerData,Caesar::Img* inputImg,bool computeStatsAndBkg=true,Caesar::BkgData* inputBkgData=0);
		/** 
		\brief Find sources
 		*/
		int FindSources(std::vector<Caesar::Source*>& sources,Caesar::Img* inputImg,bool computeStatsAndBkg=true,Caesar::BkgData* inputBkgData=0);
		
		int SelectSources(std::vector<Caesar::Source*>& sources);
		bool IsGoodSource(Caesar::Source* aSource);
		bool IsPointLikeSource(Caesar::Source* aSource);

		/** 
		\brief ThrowProgressEvent
 		*/
		int PushWorkerProgressEvent();

		/** 
		\brief PushWorkerDataEvent
 		*/
		int PushWorkerDataEvent(Caesar::WorkerData& workerData);
		

	private:	
  	static log4tango::Logger* m_logger;
    static SFinder* m_device;	
		
		mutable std::mutex m_mutex;	
		std::atomic<bool> m_stopThread;
		std::thread m_thread;

	friend class SFinder;

};//close class


}//close namespace

#endif
