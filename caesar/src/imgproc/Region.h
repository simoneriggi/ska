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
* @file Region.h
* @class Region
* @brief Region data class
*
* Superpixel Region data
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef Region_h
#define Region_h 1

#include <Blob.h>

#include <TObject.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <iostream>
#include <time.h>
#include <ctime>

using namespace std;


namespace Caesar {


class Region : public Blob {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		Region();
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~Region();

		enum RegionTag {eBkgTag=0,eSignalTag=1,eUntagged=2};
	
		struct RegionPars {
			TVectorD* pars;
			TVectorD* robustPars;
			TVectorD* spatialPars;
			RegionPars() {
				pars= 0;
				robustPars= 0;
				spatialPars= 0;
			}
			~RegionPars() { 
				if(pars) pars->Delete();
				if(robustPars) robustPars->Delete();
				if(spatialPars) spatialPars->Delete();	
			}
		};//close RegionPars()

	public:
		Region::RegionPars* GetParams(bool includeCurvPar=true);
		int GetDistance(double& dist_color,double& dist_space,Region* aRegion,bool useRobustParams=false,bool normalizeParams=true,bool addCurvDist=true);
		int GetAsymmDistance(double& dist,double& dist_neighbor,Region* aRegion,bool useRobustParams=false,bool normalizeParams=true,bool addSpatialDist=false,bool addCurvDist=true);
		int AddRegion(Region* aRegion,bool addPixels=true,bool copyPixels=false);
		void AddSubRegionId(int id){m_SubRegionIds.push_back(id);}
		int GetNSubRegions(){return (int)(m_SubRegionIds.size());}
		const std::vector<int>& GetSubRegionIds() const {return m_SubRegionIds;}
		int GetSubRegionId(int index){
			if(GetNSubRegions()<=0 || index<0 || index>=GetNSubRegions() )
				return -1;
			return m_SubRegionIds[index];
		}

	public:
		int Tag;

	protected:
		std::vector<int> m_SubRegionIds;

	ClassDef(Region,1)

};//close Region()


class RegionCollection : public TObject {
	
	public:
		RegionCollection(){
			regions.clear();
		};
		virtual ~RegionCollection(){
			for(unsigned int i=0;i<regions.size();i++){
				if(regions[i]) {
					delete regions[i];	
					regions[i]= 0;
				}
			}
			regions.clear();
		};

	private:
		struct MatchId {
 			MatchId(const int& id) : m_id(id) {}
 			bool operator()(const Region* obj) const {
   			return obj->Id == m_id;
 			}
 			private:
   			const int& m_id;
		};
		

	public: 
		void Add(Region* aRegion){regions.push_back(aRegion);}
		int GetN(){return (int)(regions.size());}

		/**
		* \brief Get region by id
		*/
		Region* FindRegionById(int id){
			int index = FindRegion(id);
			if(index<0) return 0;
			return regions[index];
		}
	
		/**
		* \brief Find region (return index with given id)
		*/
		int FindRegion(int id){
			if(GetN()<=0) return -1;
			std::vector<Region*>::iterator it = std::find_if(regions.begin(), regions.end(), MatchId(id));
			if (it==regions.end()) return 0;//not found in collection
			int index = it-regions.begin();
			return index;
		}
	
		/**
		* \brief Get region with given index
		*/
		Region* GetRegion(int index){
			if(index<0 || index>=GetN()) return 0;
			return regions[index];
		}

		/**
		* \brief Get region map regionId-->index
		*/
		std::map<int,int> GetRegionIdMap() const {
			std::map<int,int> regionIdMap;
			for(unsigned int k=0;k<regions.size();k++){
				int regionId= regions[k]->Id;
				regionIdMap.insert( std::pair<int,int>(regionId,k) );
			}
			return regionIdMap;
		}//close GetRegionMapping()

		/**
		* \brief Get region index map (index-->regionId)
		*/
		std::map<int,int> GetRegionIndexMap() const {
			std::map<int,int> regionIndexMap;
			for(unsigned int k=0;k<regions.size();k++){
				int regionId= regions[k]->Id;
				regionIndexMap.insert( std::pair<int,int>(k,regionId) );
			}
			return regionIndexMap;
		}//close GetRegionMapping()
	

	public:
		std::vector<Region*> regions;
		

	ClassDef(RegionCollection,1)

};//close RegionCollection


#ifdef __MAKECINT__
#pragma link C++ class Region+;
#pragma link C++ class RegionCollection+;
#endif

}//close namespace

#endif

