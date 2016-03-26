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
* @file Source.h
* @class Source
* @brief Source class
*
* Class representing an image source
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef Source_h
#define Source_h 1

#include <Blob.h>
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

class Contour;
class Img;

class Source : public Blob {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		Source();
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~Source();

		enum SourceType {eUnknown=0,eCompact=1,ePointLike=2,eExtended=3};
		enum SourceFlag {eReal=1,eCandidate=2,eFake=3};

	public:
		void SetType(SourceType choice){Type=choice;}
		void SetFlag(SourceFlag choice){Flag=choice;}
		void SetBeamFluxIntegral(double val){m_BeamFluxIntegral= val;}
		bool IsGoodSource(){return m_IsGoodSource;}
		void SetGoodSourceFlag(bool flag){m_IsGoodSource=flag;}
		void SetDepthLevel(int level){m_DepthLevel=level;}		

		/**
		* \brief Is source inside given source
		*/
		bool IsInsideSource(Source* aSource){
			if(!aSource) return false;
			bool isInsideX= (m_Xmin>=aSource->m_Xmin && m_Xmax<=aSource->m_Xmax);
			bool isInsideY= (m_Ymin>=aSource->m_Ymin && m_Ymax<=aSource->m_Ymax);
			bool isInside= (isInsideX && isInsideY);			
			return isInside;
		}	
		/**
		* \brief Add nested sources
		*/
		void AddNestedSource(Source* aNestedSource){
			if(!aNestedSource) return;
			int nNestedSources= (int)m_NestedSources.size();
			int nestedId= nNestedSources+1;
			TString nestedName= Form("%s_N%d",Name.c_str(),nestedId);
			aNestedSource->Id= nestedId;
			aNestedSource->Type= aNestedSource->Type;
			aNestedSource->m_DepthLevel= this->m_DepthLevel+1;
			m_NestedSources.push_back(aNestedSource);
			m_HasNestedSources= true;
		}	
		/**
		* \brief Has nested sources?
		*/
		bool HasNestedSources(){return (m_HasNestedSources && m_NestedSources.size()>0);}
		/**
		* \brief Get nested sources
		*/
		std::vector<Source*>& GetNestedSources(){return m_NestedSources;}
		/**
		* \brief Draw contours
		*/
		void Draw(bool drawBoundingBox=false,bool drawFittedEllipse=false,bool drawNested=false,int lineColor=kBlack);

		/**
		* \brief Get DS9 region info
		*/
		const std::string GetDS9Region(bool dumpNestedSourceInfo=false);
		/**
		* \brief Get DS9 ellipse info
		*/
		const std::string GetDS9EllipseRegion(bool dumpNestedSourceInfo=false);

	public:
		int Type;
		int Flag;
	private:
		double m_BeamFluxIntegral;

		//Is good source?
		bool m_IsGoodSource;

		//Nested source info
		int m_DepthLevel;
		bool m_HasNestedSources;
		Source* m_NestedSource;
		std::vector<Source*> m_NestedSources;	

		ClassDef(Source,1)

};//close Source()

typedef std::vector<Source*> SourceCollection;

#ifdef __MAKECINT__
#pragma link C++ class Source+;
#pragma link C++ class vector<Source>+;
#pragma link C++ class vector<Source*>+;
#endif

}//close namespace

#endif
