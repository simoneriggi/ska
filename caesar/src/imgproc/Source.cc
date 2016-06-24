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
* @file Source.cc
* @class Source
* @brief Source class
*
* Class representing an image source
* @author S. Riggi
* @date 20/01/2015
*/


#include <Source.h>
#include <Blob.h>
#include <Img.h>
#include <Contour.h>

#include <TObject.h>
#include <TMatrixD.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <numeric>
#include <time.h>
#include <ctime>
#include <queue>

using namespace std;

ClassImp(Caesar::Source)

namespace Caesar {


Source::Source() : Blob() {

	Init();

}//close costructor

Source::~Source(){
	
}//close destructor


Source::Source(const Source& source) : Blob() {
  // Contour copy constructor
	DEBUG_LOG("Copy constuctor called...");
  Init();
  ((Source&)source).Copy(*this);
}

void Source::Copy(TObject &obj) const {

	// Copy this source to source obj
	Blob::Copy((Source&)obj);
  ((Source&)obj).Type = Type;
	((Source&)obj).Flag = Flag;	
	((Source&)obj).m_BeamFluxIntegral = m_BeamFluxIntegral;
	((Source&)obj).m_IsGoodSource = m_IsGoodSource;	
	((Source&)obj).m_DepthLevel = m_DepthLevel;
	((Source&)obj).m_HasNestedSources = m_HasNestedSources;
	
	//Delete first a previously existing vector
	for(unsigned int i=0;i<(((Source&)obj).m_NestedSources).size();i++){
		if( (((Source&)obj).m_NestedSources)[i] ){
			delete (((Source&)obj).m_NestedSources)[i];
			(((Source&)obj).m_NestedSources)[i]= 0;
		}
	}
	(((Source&)obj).m_NestedSources).clear();

	((Source&)obj).m_NestedSource= 0;
	for(unsigned int i=0;i<m_NestedSources.size();i++){
		((Source&)obj).m_NestedSource= new Source;
		*(((Source&)obj).m_NestedSource)= *(m_NestedSources[i]);
		(((Source&)obj).m_NestedSources).push_back( ((Source&)obj).m_NestedSource );
	}

}//close Copy()

Source& Source::operator=(const Source& source) { 
	// Operator =
  if (this != &source)  ((Source&)source).Copy(*this);
  return *this;
}

void Source::Init(){

	Type= eUnknown;
	Flag= eCandidate;
	m_BeamFluxIntegral= 0;
	
	m_IsGoodSource= true;

	m_DepthLevel= 0;
	m_HasNestedSources= false;
	m_NestedSource= 0;
	m_NestedSources.clear();

}//close Init()


void Source::Draw(bool drawBoundingBox,bool drawEllipse,bool drawNested,int lineColor){

	//Drawing contours?
	for(unsigned int i=0;i<m_Contours.size();i++){		
		TGraph* thisContourGraph= m_Contours[i]->GetGraph();
		if(thisContourGraph) {
			thisContourGraph->SetMarkerSize(8);
			thisContourGraph->SetMarkerSize(0.3);
			thisContourGraph->SetMarkerColor(lineColor);
			thisContourGraph->SetLineColor(lineColor);
			thisContourGraph->SetLineStyle(kSolid);
			thisContourGraph->SetLineWidth(2);
			thisContourGraph->Draw("Lsame");
		}//close if 
		
		//Get bounding box
		if(drawBoundingBox){
			TPolyLine* thisBoundingBox= m_Contours[i]->GetBoundingBoxLine();
			if(thisBoundingBox){
				thisBoundingBox->SetLineColor(lineColor);
				thisBoundingBox->SetLineStyle(kDashed);
				thisBoundingBox->Draw("lsame");
			}
		}//close if

		//Get fitted ellipse
		if(drawEllipse){
			TEllipse* thisFittedEllipse= m_Contours[i]->GetFittedEllipse();
			if(thisFittedEllipse){
				thisFittedEllipse->SetLineColor(lineColor);
				thisFittedEllipse->SetLineStyle(kDotted);
				thisFittedEllipse->SetLineWidth(2);
				thisFittedEllipse->SetFillColor(0);
				thisFittedEllipse->SetFillStyle(0);
				thisFittedEllipse->Draw("lsame");
			}
		}//close if draw ellipse
	}//end loop contours

	if(drawNested){
		for(unsigned int i=0;i<m_NestedSources.size();i++){
			if(m_NestedSources[i]) m_NestedSources[i]->Draw(drawBoundingBox,drawEllipse,drawNested,lineColor+2);
		}
	}//close if nested

}//close Draw()


const std::string Source::GetDS9Region(bool dumpNestedSourceInfo){

	//global color=red dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 image
	std::stringstream sstream;
	sstream<<"polygon ";
	for(unsigned int i=0; i<m_Contours.size(); i++){ 
		int nPoints= m_Contours[i]->GetN();
		for(int j=0;j<nPoints;j++){
			TVector2* contPnt= m_Contours[i]->GetPoint(j);
			if(!contPnt) continue;
			sstream<<(int)contPnt->X()+1<<" "<<(int)contPnt->Y()+1<<" ";
		}
	}
	sstream<<"# text={S"<<Id<<"}";

	if(dumpNestedSourceInfo && m_HasNestedSources){			
		sstream<<endl;
		for(unsigned int k=0;k<m_NestedSources.size();k++){
			sstream<<"polygon ";
			std::vector<Contour*> nestedContours= m_NestedSources[k]->m_Contours;
			for(unsigned int i=0; i<nestedContours.size(); i++){ 
				int nPoints= nestedContours[i]->GetN();
				for(int j=0;j<nPoints;j++){
					TVector2* contPnt= nestedContours[i]->GetPoint(j);
					if(!contPnt) continue;
					sstream<<(int)contPnt->X()+1<<" "<<(int)contPnt->Y()+1<<" ";
				}
			}//end loop contours
			sstream<<"# text={S"<<Id<<"_Nest"<<k<<"}";
			if(k!=m_NestedSources.size()-1) sstream<<endl;
		}//end loop nested sources
	}//close dumpNestedSourceInfo

	const std::string dsregions= sstream.str();
	return dsregions;

}//close GetDS9Region()


const std::string Source::GetDS9EllipseRegion(bool dumpNestedSourceInfo){
			
	//ellipse x y radius radius angle
	std::stringstream sstream;
	sstream<<"ellipse ";
	for(unsigned int i=0; i<m_Contours.size(); i++){ 
		if(!m_Contours[i]->HasEllipseFit) continue;
		double EllX= m_Contours[i]->EllipseCenter.X();
		double EllY= m_Contours[i]->EllipseCenter.Y();
		double EllMajAxis= m_Contours[i]->EllipseMajAxis;
		double EllMinAxis= m_Contours[i]->EllipseMinAxis;
		double EllRotAxis= m_Contours[i]->EllipseRotAngle-90;//-90 comes from DS9 strange format!
		//TVector2 BBoxCenter= m_Contours[i]->BoundingBoxCenter;
		//double BBoxMinAxis=  m_Contours[i]->BoundingBoxMin;	
		//double BBoxMajAxis= m_Contours[i]->BoundingBoxMaj;
		//double BBoxAngle= m_Contours[i]->BoundingBoxAngle;	
		sstream<<EllX+1<<" "<<EllY+1<<" "<<(EllMajAxis/2)<<" "<<(EllMinAxis/2)<<" "<<(EllRotAxis)<<" ";
	}
	sstream<<"# text={S"<<Id<<"}";

	if(dumpNestedSourceInfo && m_HasNestedSources){			
		sstream<<endl;
		for(unsigned int k=0;k<m_NestedSources.size();k++){
			sstream<<"ellipse ";
			std::vector<Contour*> nestedContours= m_NestedSources[k]->m_Contours;
			for(unsigned int i=0; i<nestedContours.size(); i++){ 
				if(!nestedContours[i]->HasEllipseFit) continue;
				double EllX= nestedContours[i]->EllipseCenter.X();
				double EllY= nestedContours[i]->EllipseCenter.Y();
				double EllMajAxis= nestedContours[i]->EllipseMajAxis;
				double EllMinAxis= nestedContours[i]->EllipseMinAxis;
				double EllRotAxis= nestedContours[i]->EllipseRotAngle-90;//-90 comes from DS9 strange format!
				//TVector2 BBoxCenter= nestedContours[i]->BoundingBoxCenter;
				//double BBoxMinAxis=  nestedContours[i]->BoundingBoxMin;	
				//double BBoxMajAxis= nestedContours[i]->BoundingBoxMaj;
				//double BBoxAngle= nestedContours[i]->BoundingBoxAngle;	
				sstream<<EllX+1<<" "<<EllY+1<<" "<<(EllMajAxis/2)<<" "<<(EllMinAxis/2)<<" "<<(EllRotAxis)<<" ";
			}//end loop contours
			sstream<<"# text={S"<<Id<<"_Nest"<<k<<"}";
			if(k!=m_NestedSources.size()-1) sstream<<endl;
		}//end loop nested sources
	}//close dumpNestedSourceInfo

	const std::string dsregions= sstream.str();
	return dsregions;

}//close GetDS9EllipseRegion()

}//close namespace


