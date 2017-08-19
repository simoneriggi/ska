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
* @file SLICData.h
* @class SLICData
* @brief SLIC data class
*
* Superpixel data
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef SLICData_h
#define SLICData_h 1

#include <Region.h>
#include <Contour.h>
#include <Consts.h>

#include <TObject.h>
#include <TMatrixD.h>


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


namespace Caesar {

//enum SLICEdgeModel {eKirschEdge,eChanVeseEdge};

class Img;
class Image;

class SLICData : public TObject {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SLICData();
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~SLICData();

	public:	
		/**
		* \brief Clear data
		*/
		void Clear();
		/**
		* \brief Clear images
		*/
		void ClearImages();
		/**
		* \brief Clear regions
		*/
		void ClearRegions();

		/**
		* \brief Get number of regions
		*/
		int GetNRegions() const {return (int)(regions.size());}
		/**
		* \brief Get region
		*/
		Region* GetRegion(int index){
			if(GetNRegions()<=0 || index<0 || index>=GetNRegions()) return 0;
			return regions[index];
		}
		/**
		* \brief Add region
		*/
		void AddRegion(Region* aRegion){	
			if(!aRegion) return;
			regions.push_back(aRegion);
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

		Image* inputImage;//the image passed to SLIC generator (normalized)
		Image* edgeImage;//the image borders passed to SLIC generator
		Image* laplImage;//the image laplacian

		Img* inputImg;//the image passed to SLIC generator (normalized)
		Img* edgeImg;//the image borders passed to SLIC generator
		Img* laplImg;//the image laplacian
		//TMatrixD* pixelLabels;//matrix with pixel region labels
		std::vector< std::vector<long int> > labels;
		std::vector<Region*> regions;//list of generated superpixels

		ClassDef(SLICData,1)

};//close SLICData


typedef std::map<int,std::vector<int>> SLICBoundaryPixMap;//key: neighbour id, value: list of shared pix ids
typedef SLICBoundaryPixMap::iterator SLICBoundaryPixMapIterator;
typedef std::vector< std::vector<int> > SLICConnectedRegions;


class SLICContourData : public TObject {

	public:
	
		SLICContourData(){
			contour= 0;
			ResetList();
		}
		virtual ~SLICContourData() { 
			ResetContour();		
			ResetList();
		}

	public:
		void ResetContour(){
			if(!contour) return;
			delete contour;
			contour= 0;	
		}
		void ResetList(){
			connectedRegionIds.clear();
			boundaryData.clear();
		}
	
	public:
		Contour* contour;
		std::vector<SLICBoundaryPixMap> boundaryData;
		SLICConnectedRegions connectedRegionIds; 

	ClassDef(SLICContourData,1)

};//close SLICContourData


class SLICNeighborData : public TObject {
	
	public:
		SLICNeighborData(){
			Order= 0;
			Id= -1;
			Index= -1;
			Tag= -1;
			D= 0;
			D_n= 0;
			E= 0;
			E_n= 0;
			Dtot= 0;
			Dtot_n= 0;
			Dsym= 0;
		};
		virtual ~SLICNeighborData(){};

	public:
		bool operator==(const SLICNeighborData& aNeighborData) const {
    	return ( (aNeighborData.Id==Id) && (aNeighborData.Index==Index) );
    }
	
		void Print(){
			cout<<"** NEIGHBOR DATA #"<<Index<<" **"<<endl;
			cout<<"Id="<<Id<<", D="<<D<<" Dtot="<<Dtot<<" Dsym="<<Dsym<<endl;
			cout<<"*******************"<<endl;
		}

	public:
		int Order;
		int Id;
		int Index;
		int Tag;	
		double D;//dissimilarity Dij
		double D_n;//dissimilarty Dji
		double E;//edgeness
		double E_n;
		double Dtot;//total dissimilarity (normalized)
		double Dtot_n;
		double Dsym;//symmetric dissimilarity

	ClassDef(SLICNeighborData,1)

};//close SLICNeighborData
typedef std::vector< std::vector<SLICNeighborData> > SLICNeighbors;


class SLICNeighborCollection : public TObject {
	
	public:
		SLICNeighborCollection(){
			m_neighbors.clear();
		};
		virtual ~SLICNeighborCollection(){};

	private:
		struct MatchId {
 			MatchId(const int& id) : m_id(id) {}
 			bool operator()(const SLICNeighborData& obj) const {
   			return obj.Id == m_id;
 			}
 			private:
   			const int& m_id;
		};
		struct MatchIndex {
 			MatchIndex(const int& index) : m_index(index) {}
 			bool operator()(const SLICNeighborData& obj) const {
   			return obj.Index == m_index;
 			}
 			private:
   			const int& m_index;
		};

		static bool compareByDiss(const SLICNeighborData &a, const SLICNeighborData &b) {
   		return ( (a.D) < (b.D) );
		}
		static bool compareByEdgeness(const SLICNeighborData &a, const SLICNeighborData &b) {
   		return ( (a.E) < (b.E) );
		}
		static bool compareByTotDiss(const SLICNeighborData &a, const SLICNeighborData &b) {
   		return ( (a.Dtot) < (b.Dtot) );
		}

	public: 
		void Add(SLICNeighborData nn){m_neighbors.push_back(nn);}
		void SortByDiss(){
			std::sort(m_neighbors.begin(), m_neighbors.end(), compareByDiss);
		};
		void SortByEdgeness(){
			std::sort(m_neighbors.begin(), m_neighbors.end(), compareByEdgeness);
		};
		void SortByTotDiss(){
			std::sort(m_neighbors.begin(), m_neighbors.end(), compareByTotDiss);
		};
		int GetN(){return (int)(m_neighbors.size());}

		int FindById(int id){
			if(GetN()<=0) return -1;
			std::vector<SLICNeighborData>::iterator it = std::find_if(m_neighbors.begin(), m_neighbors.end(), MatchId(id));
			if (it==m_neighbors.end()) return -1;//not found in collection
			int pos = it-m_neighbors.begin();
			return pos;
		}
		int FindByIndex(int index){
			if(GetN()<=0) return -1;
			std::vector<SLICNeighborData>::iterator it = std::find_if(m_neighbors.begin(), m_neighbors.end(), MatchIndex(index));
			if (it== m_neighbors.end()) return -1;//not found in collection
			size_t pos = it-m_neighbors.begin();
			return pos;
		}

		int FindCloserByDiss(){
			if(GetN()<=0) return -1;
			std::vector<SLICNeighborData>::iterator it= std::min_element(m_neighbors.begin(),m_neighbors.end(),compareByDiss);
			size_t pos = it-m_neighbors.begin();
			return pos;
		}
		int FindCloserByDissTot(){
			if(GetN()<=0) return -1;
			std::vector<SLICNeighborData>::iterator it= std::min_element(m_neighbors.begin(),m_neighbors.end(),compareByTotDiss);
			size_t pos = it-m_neighbors.begin();
			return pos;
		}

		std::vector<SLICNeighborData> GetNSortedByDiss(int N){
			std::vector<SLICNeighborData> topSorted;
			topSorted.clear();
			int collSize= GetN();
			if(collSize<=0 || N>collSize) return topSorted;
			
			//Partially sort Nth elements
			std::vector<SLICNeighborData> tmp;
			tmp.assign(m_neighbors.begin(),m_neighbors.end());
			partial_sort(tmp.begin(), tmp.begin()+N, tmp.end(), compareByDiss);

			//Copy first N sorted to a new collection
			topSorted.assign(m_neighbors.begin(), m_neighbors.begin()+N);
			
			return topSorted;
		}

		std::vector<SLICNeighborData> GetNSortedByTotDiss(int N){
			std::vector<SLICNeighborData> topSorted;
			topSorted.clear();
			int collSize= GetN();
			if(collSize<=0 || N>collSize) return topSorted;
			
			//Partially sort Nth elements
			std::vector<SLICNeighborData> tmp;
			tmp.assign(m_neighbors.begin(),m_neighbors.end());
			partial_sort(tmp.begin(), tmp.begin()+N, tmp.end(), compareByTotDiss);

			//Copy first N sorted to a new collection
			topSorted.assign(tmp.begin(), tmp.begin()+N);
			
			return topSorted;
		}

		bool IsIdAmongNClosersByDiss(int id,int N){		
			//Get first N sorted in a new collection
			std::vector<SLICNeighborData> topSorted= GetNSortedByDiss(N);
			if(topSorted.size()<=0) return false;

			//Find if given id is among the top sorted list	
			std::vector<SLICNeighborData>::iterator it = std::find_if(topSorted.begin(),topSorted.end(), MatchId(id));
			if (it!=topSorted.end()) return true;//found in top sorted collection!

			return false;
		}

		bool IsIdAmongNClosersByTotDiss(int id,int N){		
			//Get first N sorted in a new collection
			std::vector<SLICNeighborData> topSorted= GetNSortedByTotDiss(N);
			if(topSorted.size()<=0) return false;

			//Find if given id is among the top sorted list	
			std::vector<SLICNeighborData>::iterator it = std::find_if(topSorted.begin(),topSorted.end(), MatchId(id));
			if (it!=topSorted.end()) return true;//found in top sorted collection!

			return false;
		}

		bool IsIndexAmongNClosersByDiss(int index,int N){
			//Get first N sorted in a new collection
			std::vector<SLICNeighborData> topSorted= GetNSortedByDiss(N);
			if(topSorted.size()<=0) return false;

			//Find if given id is among the top sorted list	
			std::vector<SLICNeighborData>::iterator it = std::find_if(topSorted.begin(),topSorted.end(), MatchIndex(index));
			if (it!=topSorted.end()) return true;//found in top sorted collection!

			return false;
		}		

		bool IsIndexAmongNClosersByTotDiss(int index,int N){
			//Get first N sorted in a new collection
			std::vector<SLICNeighborData> topSorted= GetNSortedByTotDiss(N);
			if(topSorted.size()<=0) return false;

			//Find if given id is among the top sorted list	
			std::vector<SLICNeighborData>::iterator it = std::find_if(topSorted.begin(),topSorted.end(), MatchIndex(index));
			if (it!=topSorted.end()) return true;//found in top sorted collection!

			return false;
		}		

		void Print(){
			for(int i=0;i<GetN();i++) m_neighbors[i].Print();
		}

		//std::vector<SLICNeighborData> GetNeighbors(){return m_neighbors;}
		const std::vector<SLICNeighborData>& GetNeighbors() const { return m_neighbors; }

		SLICNeighborData* GetNeighbor(int index){
			if(index<0 || index>=GetN()) return 0;
			return (m_neighbors.data()+index);
		}

	private:
		std::vector<SLICNeighborData> m_neighbors;
		
	ClassDef(SLICNeighborCollection,1)

};//close SLICNeighborCollection
typedef std::vector<SLICNeighborCollection> SLICNeighborCollections;




class SLICSimilarityData : public TObject {

	public:
	
		SLICSimilarityData(){
			//DissimilarityMatrix= 0;
			AdjacencyMatrix= 0;
			//AbsDissimilarityMatrix= 0;
			//SaliencyDissimilarityMatrix= 0;
			//NeighborMatrix= 0;
		}
		virtual ~SLICSimilarityData() { 
			//if(DissimilarityMatrix) DissimilarityMatrix->Delete();
			if(AdjacencyMatrix) AdjacencyMatrix->Delete();
			//if(AbsDissimilarityMatrix) AbsDissimilarityMatrix->Delete();
			//if(SaliencyDissimilarityMatrix) SaliencyDissimilarityMatrix->Delete();
			//if(NeighborMatrix) NeighborMatrix->Delete();
			//cout<<"SLICUtils::~SLICSimilarityData(): Clearing vector..."<<endl;
			//DissimilaritySortIndexMatrix.clear();
		}

	public:
		//TMatrixD* DissimilarityMatrix;
		TMatrixD* AdjacencyMatrix;
		//TMatrixD* AbsDissimilarityMatrix;
		//TMatrixD* SaliencyDissimilarityMatrix;
		//TMatrixD* NeighborMatrix;
		//std::vector< std::vector<int> > DissimilaritySortIndexMatrix;
		//std::vector<double> saliencyList;
		double Dmin;
		double Dmax;
		double Dmedian;
		double Dmedianrms;
		double Emin;
		double Emax;
			
	ClassDef(SLICSimilarityData,1)
	
};//close SLICSimilarityData


#ifdef __MAKECINT__
#pragma link C++ class SLICData+;
#pragma link C++ enum SLICEdgeModel+;
#pragma link C++ class SLICContourData+;
#pragma link C++ class SLICNeighborData+;
#pragma link C++ class SLICNeighborCollection+;
#pragma link C++ class SLICSimilarityData+;
#endif

}//close namespace

#endif

