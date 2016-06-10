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
* @file Serializer.cc
* @class Serializer
* @brief Serializer class
*
* Class for serializing objects
* @author S. Riggi
* @date 20/01/2015
*/

#include <Serializer.h>
#include <Logger.h>

#include <TVector2.h>

//# MSG PACK
//#ifdef BUILD_CAESAR_SERVER
	#include <tango.h>
	#include <msgpack.hpp>
//#endif

#include <Source.pb.h>

#include <string>
#include <ctime>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <sstream>

#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <map>
#include <exception>

#include <chrono>

using namespace std;

ClassImp(Caesar::SBuffer)
ClassImp(Caesar::Serializer)

namespace Caesar {

Serializer::Serializer() {	
	
}   


Serializer::~Serializer(){

}

//#ifdef BUILD_CAESAR_SERVER
int Serializer::EncodePointToProtobuf(SourcePB::Point& point_pb,TVector2& point){

	try {
		point_pb.set_x(point.X());
		point_pb.set_y(point.Y());
	}
	catch(std::exception const & e) {
		ERROR_LOG("Point encoding to protobuf failed with status "<<e.what());
		return -1;
	}

	return 0;

}//close EncodePointToProtobuf()

int Serializer::EncodeContourToProtobuf(SourcePB::Contour& contour_pb,Contour* contour){
	
	if(!contour) return -1;

	try {
		contour_pb.set_hasparameters(contour->HasParameters);
		contour_pb.set_area(contour->Area);
		contour_pb.set_perymeter(contour->Perymeter);
		contour_pb.set_isconvexcontour(contour->IsConvexContour);
		contour_pb.set_circularityratio(contour->CircularityRatio);
			
		SourcePB::Point* BoundingBoxCenter= new SourcePB::Point;
		if(EncodePointToProtobuf(*BoundingBoxCenter,contour->BoundingBoxCenter)<0){
			throw std::runtime_error("Failed to encode bounding box center field");
		}
		contour_pb.set_allocated_boundingboxcenter(BoundingBoxCenter);

		contour_pb.set_boundingboxmaj(contour->BoundingBoxMaj);
		contour_pb.set_boundingboxmin(contour->BoundingBoxMin);
		contour_pb.set_boundingboxangle(contour->BoundingBoxAngle);
		contour_pb.set_elongation(contour->Elongation);
		contour_pb.set_rectangularity(contour->Rectangularity);
		contour_pb.set_roundness(contour->Roundness);
		contour_pb.set_eccentricity(contour->Eccentricity);
		contour_pb.set_tiltangle(contour->TiltAngle);
		contour_pb.set_hasellipsefit(contour->HasEllipseFit);

		SourcePB::Point* EllipseCenter= new SourcePB::Point;
		if(EncodePointToProtobuf(*EllipseCenter,contour->EllipseCenter)<0){
			throw std::runtime_error("Failed to encode ellipse center field");
		}
		contour_pb.set_allocated_ellipsecenter(EllipseCenter);

		contour_pb.set_ellipsemajaxis(contour->EllipseMajAxis);
		contour_pb.set_ellipseminaxis(contour->EllipseMinAxis);
		contour_pb.set_ellipserotangle(contour->EllipseRotAngle);
		contour_pb.set_ellipsefitredchi2(contour->EllipseFitRedChi2);
		contour_pb.set_ellipsearearatio(contour->EllipseAreaRatio);

		contour_pb.set_m00(contour->m00);
		contour_pb.set_m10(contour->m10);
		contour_pb.set_m01(contour->m01);
		contour_pb.set_m20(contour->m20);
		contour_pb.set_m11(contour->m11);
		contour_pb.set_m02(contour->m02);
		contour_pb.set_m30(contour->m30);
		contour_pb.set_m21(contour->m21);
		contour_pb.set_m12(contour->m12);
		contour_pb.set_m03(contour->m03);

		contour_pb.set_mu20(contour->mu20);
		contour_pb.set_mu11(contour->mu11);
		contour_pb.set_mu02(contour->mu02);
		contour_pb.set_mu30(contour->mu30);
		contour_pb.set_mu21(contour->mu21);
		contour_pb.set_mu12(contour->mu12);
		contour_pb.set_mu03(contour->mu03);
		
		contour_pb.set_nu20(contour->nu20);
		contour_pb.set_nu11(contour->nu11);
		contour_pb.set_nu02(contour->nu02);
		contour_pb.set_nu30(contour->nu30);
		contour_pb.set_nu21(contour->nu21);
		contour_pb.set_nu12(contour->nu12);
		contour_pb.set_nu03(contour->nu03);

		for(unsigned int i=0;i<contour->HuMoments.size();i++){		
 			contour_pb.add_humoments(contour->HuMoments[i]);
		}
	
		for(unsigned int i=0;i<contour->BoundingBoxVertex.size();i++){		
			SourcePB::Point* thisBoundingBoxVertex = contour_pb.add_boundingboxvertex();
			if(EncodePointToProtobuf(*thisBoundingBoxVertex,(contour->BoundingBoxVertex)[i])<0){
				std::stringstream errMsg;
				errMsg<<"BoundingBoxVertex no. "<<i+1<<" encoding to protobuf failed!";
				throw std::runtime_error(errMsg.str().c_str());
			}
		}
	
		SourcePB::Point* Centroid= new SourcePB::Point;
		if(EncodePointToProtobuf(*Centroid,contour->Centroid)<0){
			throw std::runtime_error("Failed to encode ellipse center field");
		}
		contour_pb.set_allocated_centroid(EllipseCenter);

		for(unsigned int i=0;i<contour->RealFDs.size();i++){		
 			contour_pb.add_realfds(contour->RealFDs[i]);
		}
		for(unsigned int i=0;i<contour->ImagFDs.size();i++){		
 			contour_pb.add_imagfds(contour->ImagFDs[i]);
		}
		for(unsigned int i=0;i<contour->ModFDs.size();i++){		
 			contour_pb.add_modfds(contour->ModFDs[i]);
		}
		for(unsigned int i=0;i<contour->BendingEnergies.size();i++){		
 			contour_pb.add_bendingenergies(contour->BendingEnergies[i]);
		}	
		for(unsigned int i=0;i<contour->CentroidDistanceModFDs.size();i++){		
 			contour_pb.add_centroiddistancemodfds(contour->CentroidDistanceModFDs[i]);
		}

		for(unsigned int i=0;i<contour->GetN();i++){		
			SourcePB::Point* thisPoint = contour_pb.add_m_points();
			if(EncodePointToProtobuf(*thisPoint,*contour->GetPoint(i))<0){
				std::stringstream errMsg;
				errMsg<<"Point no. "<<i+1<<" encoding to protobuf failed!";
				throw std::runtime_error(errMsg.str().c_str());
			}
		}
	}//close try block
	catch(std::exception const & e) {
		ERROR_LOG("Point encoding to protobuf failed with status "<<e.what());
		return -1;
	}

	return 0;

}//close EncodeContourToProtobuf()


int Serializer::EncodePixelToProtobuf(SourcePB::Pixel& pixel_pb,Pixel* pixel){
	
	if(!pixel) return -1;

	try {
		//Set public source fields
		pixel_pb.set_id(pixel->id);
		if(pixel->type==Pixel::eNormal){
			pixel_pb.set_type(SourcePB::Pixel::eNormal);
		}
		else if(pixel->type==Pixel::eSeed){
			pixel_pb.set_type(SourcePB::Pixel::eSeed);
		}
		else if(pixel->type==Pixel::eHalo){
			pixel_pb.set_type(SourcePB::Pixel::eHalo);
		}
		else{
			throw std::runtime_error("Invalid pixel type, failed to encode!");
		}

		pixel_pb.set_s(pixel->S);
		pixel_pb.set_x(pixel->x);
		pixel_pb.set_y(pixel->y);
		pixel_pb.set_ix(pixel->ix);
		pixel_pb.set_iy(pixel->iy);
		pixel_pb.set_isonedge(pixel->isOnEdge);
		pixel_pb.set_distancetoedge(pixel->distanceToEdge);

		pixel_pb.set_s_curv(pixel->GetCurv());
		pixel_pb.set_s_edge(pixel->GetEdge());
		
		std::pair<double,double> bkginfo= pixel->GetBkg();
		pixel_pb.set_bkglevel(bkginfo.first);
		pixel_pb.set_noiselevel(bkginfo.second);

	}//close try block
	catch(std::exception const & e) {
		ERROR_LOG("Pixel encoding to protobuf failed with status "<<e.what());
		return -1;
	}


	return 0;
}

int Serializer::EncodeBlobToProtobuf(SourcePB::Blob& blob_pb,Source* source){
	
	if(!source) return -1;

	try{
		blob_pb.set_haspixelsatedge(source->HasPixelsAtEdge);
		blob_pb.set_id(source->Id);
		blob_pb.set_name(source->Name);
		blob_pb.set_npix(source->NPix);
		blob_pb.set_mean(source->Mean);
		blob_pb.set_rms(source->RMS);
		blob_pb.set_skewness(source->Skewness);	
		blob_pb.set_median(source->Median);
		blob_pb.set_medianrms(source->MedianRMS);
		blob_pb.set_x0(source->X0);
		blob_pb.set_y0(source->Y0);
		blob_pb.set_mean_curv(source->Mean_curv);
		blob_pb.set_rms_curv(source->RMS_curv);
		blob_pb.set_median_curv(source->Median_curv);
		blob_pb.set_medianrms_curv(source->MedianRMS_curv);
		
		for(unsigned int i=0;i<source->Moments.size();i++) blob_pb.add_moments((source->Moments)[i]);
		for(unsigned int i=0;i<source->HuMoments.size();i++) blob_pb.add_humoments((source->HuMoments)[i]);
		for(unsigned int i=0;i<source->ZMMoments.size();i++) blob_pb.add_zmmoments((source->ZMMoments)[i]);
		blob_pb.set_m_hasstats(source->HasStats());	
		blob_pb.set_m_hasparameters(source->HasParameters());
		blob_pb.set_m_m1(source->GetM1());
		blob_pb.set_m_m2(source->GetM2());
		blob_pb.set_m_m3(source->GetM3());
		blob_pb.set_m_m4(source->GetM4());
		blob_pb.set_m_m1_curv(source->GetM1Curv());
		blob_pb.set_m_m2_curv(source->GetM2Curv());

		blob_pb.set_m_s(source->GetS());
		blob_pb.set_m_smax(source->GetSmax());
		blob_pb.set_m_smin(source->GetSmin());
		blob_pb.set_m_sxx(source->GetSxx());
		blob_pb.set_m_syy(source->GetSyy());
		blob_pb.set_m_sxy(source->GetSxy());
		blob_pb.set_m_sx(source->GetSx());
		blob_pb.set_m_sy(source->GetSy());
		blob_pb.set_m_pixidmax(source->GetSmaxPixId());
		blob_pb.set_m_pixidmin(source->GetSminPixId());
		blob_pb.set_m_s_curv(source->GetScurv());
		blob_pb.set_m_s_edge(source->GetSedge());
	
		long int imgsizex, imgsizey;
		source->GetImageSize(imgsizex,imgsizey);
		blob_pb.set_m_imagesizex(imgsizex);
		blob_pb.set_m_imagesizey(imgsizey);

		double xmin, xmax, ymin, ymax;
		source->GetImageRange(xmin,xmax,ymin,ymax);
		blob_pb.set_m_imageminx(xmin);
		blob_pb.set_m_imagemaxx(xmax);
		blob_pb.set_m_imageminy(ymin);
		blob_pb.set_m_imagemaxy(ymax);

		double imgSmin, imgSmax;
		source->GetImageSRange(imgSmin, imgSmax);
		blob_pb.set_m_imagemins(imgSmin);
		blob_pb.set_m_imagemaxs(imgSmax);

		double imgSmin_curv, imgSmax_curv;
		source->GetImageScurvRange(imgSmin_curv, imgSmax_curv);
		blob_pb.set_m_imageminscurv(imgSmin_curv);
		blob_pb.set_m_imagemaxscurv(imgSmax_curv);

		double imgSmin_edge, imgSmax_edge;
		source->GetImageSedgeRange(imgSmin_edge, imgSmax_edge);
		blob_pb.set_m_imageminsedge(imgSmin_edge);
		blob_pb.set_m_imagemaxsedge(imgSmax_edge);

		blob_pb.set_m_imagerms(source->GetImageRMS());

		double sxmin, sxmax, symin, symax;
		source->GetSourceRange(sxmin,sxmax,symin,symax);
		blob_pb.set_m_xmin(sxmin);
		blob_pb.set_m_xmax(sxmax);
		blob_pb.set_m_ymin(symin);
		blob_pb.set_m_ymax(symax);

		long int ixmin, ixmax, iymin, iymax;
		source->GetSourcePixelRange(ixmin, ixmax, iymin, iymax);
		blob_pb.set_m_ix_min(ixmin);
		blob_pb.set_m_ix_max(ixmax);
		blob_pb.set_m_iy_min(iymin);
		blob_pb.set_m_iy_max(iymax);
		
		//Add pixel collection to blob
		for(int k=0;k<source->GetNPixels();k++){
			SourcePB::Pixel* thisPixel = blob_pb.add_m_pixels();
			if(EncodePixelToProtobuf(*thisPixel,source->GetPixel(k))<0){
				std::stringstream errMsg;
				errMsg<<"Pixel no. "<<k+1<<" encoding to protobuf failed!";
				throw std::runtime_error(errMsg.str().c_str());
			}
		}

		//Add contout to blob
		//...

	}//close try block
	catch(std::exception const & e) {
		ERROR_LOG("Blob encoding to protobuf failed with status "<<e.what());
		return -1;
	}

	return 0;

}//close EncodeBlobToProtobuf()


int Serializer::EncodeSourceToProtobuf(SourcePB::Source& source_pb,Source* source){

	if(!source) return -1;

	try {
		
		source_pb.set_type(source->Type);
		source_pb.set_flag(source->Flag);
	
		//Set private source fields
		source_pb.set_m_beamfluxintegral(source->GetBeamFluxIntegral());
		source_pb.set_m_isgoodsource(source->IsGoodSource());
		source_pb.set_m_depthlevel(source->GetDepthLevel());
		source_pb.set_m_hasnestedsources(source->HasNestedSources());
		
		//Add nested sources
		for(int i=0;i<source->GetNestedSourceNumber();i++){
			SourcePB::Source* thisNestedSource = source_pb.add_m_nestedsources();
			if(EncodeSourceToProtobuf(*thisNestedSource,source->GetNestedSource(i))<0){
				std::stringstream errMsg;
				errMsg<<"Nested source no. "<<i+1<<" encoding to protobuf failed!";
				throw std::runtime_error(errMsg.str().c_str());
			}
		}

		//Set blob field
		SourcePB::Blob* blob= new SourcePB::Blob;
		if(EncodeBlobToProtobuf(*blob,source)<0){
			throw std::runtime_error("Failed to encode blob!");
		}
		source_pb.set_allocated_blob(blob);

		/*
		blob->set_haspixelsatedge(source->HasPixelsAtEdge);
		blob->set_id(source->Id);
		blob->set_name(source->Name);
		blob->set_npix(source->NPix);
		blob->set_mean(source->Mean);
		blob->set_rms(source->RMS);
		blob->set_skewness(source->Skewness);	
		blob->set_median(source->Median);
		blob->set_medianrms(source->MedianRMS);
		blob->set_x0(source->X0);
		blob->set_y0(source->Y0);
		blob->set_mean_curv(source->Mean_curv);
		blob->set_rms_curv(source->RMS_curv);
		blob->set_median_curv(source->Median_curv);
		blob->set_medianrms_curv(source->MedianRMS_curv);
		
		for(unsigned int i=0;i<source->Moments.size();i++) blob->add_moments((source->Moments)[i]);
		for(unsigned int i=0;i<source->HuMoments.size();i++) blob->add_humoments((source->HuMoments)[i]);
		for(unsigned int i=0;i<source->ZMMoments.size();i++) blob->add_zmmoments((source->ZMMoments)[i]);
		blob->set_m_hasstats(source->HasStats());	
		blob->set_m_hasparameters(source->HasParameters());
		blob->set_m_m1(source->GetM1());
		blob->set_m_m2(source->GetM2());
		blob->set_m_m3(source->GetM3());
		blob->set_m_m4(source->GetM4());
		blob->set_m_m1_curv(source->GetM1Curv());
		blob->set_m_m2_curv(source->GetM2Curv());

		blob->set_m_s(source->GetS());
		blob->set_m_smax(source->GetSmax());
		blob->set_m_smin(source->GetSmin());
		blob->set_m_sxx(source->GetSxx());
		blob->set_m_syy(source->GetSyy());
		blob->set_m_sxy(source->GetSxy());
		blob->set_m_sx(source->GetSx());
		blob->set_m_sy(source->GetSy());
		blob->set_m_pixidmax(source->GetSmaxPixId());
		blob->set_m_pixidmin(source->GetSminPixId());
		blob->set_m_s_curv(source->GetScurv());
		blob->set_m_s_edge(source->GetSedge());
	
		long int imgsizex, imgsizey;
		source->GetImageSize(imgsizex,imgsizey);
		blob->set_m_imagesizex(imgsizex);
		blob->set_m_imagesizey(imgsizey);

		double xmin, xmax, ymin, ymax;
		source->GetImageRange(xmin,xmax,ymin,ymax);
		blob->set_m_imageminx(xmin);
		blob->set_m_imagemaxx(xmax);
		blob->set_m_imageminy(ymin);
		blob->set_m_imagemaxy(ymax);

		double imgSmin, imgSmax;
		source->GetImageSRange(imgSmin, imgSmax);
		blob->set_m_imagemins(imgSmin);
		blob->set_m_imagemaxs(imgSmax);

		double imgSmin_curv, imgSmax_curv;
		source->GetImageScurvRange(imgSmin_curv, imgSmax_curv);
		blob->set_m_imageminscurv(imgSmin_curv);
		blob->set_m_imagemaxscurv(imgSmax_curv);

		double imgSmin_edge, imgSmax_edge;
		source->GetImageSedgeRange(imgSmin_edge, imgSmax_edge);
		blob->set_m_imageminsedge(imgSmin_edge);
		blob->set_m_imagemaxsedge(imgSmax_edge);

		blob->set_m_imagerms(source->GetImageRMS());

		double sxmin, sxmax, symin, symax;
		source->GetSourceRange(sxmin,sxmax,symin,symax);
		blob->set_m_xmin(sxmin);
		blob->set_m_xmax(sxmax);
		blob->set_m_ymin(symin);
		blob->set_m_ymax(symax);

		long int ixmin, ixmax, iymin, iymax;
		source->GetSourcePixelRange(ixmin, ixmax, iymin, iymax);
		blob->set_m_ix_min(ixmin);
		blob->set_m_ix_max(ixmax);
		blob->set_m_iy_min(iymin);
		blob->set_m_iy_max(iymax);
		
		//Add pixel collection to blob
		for(unsigned int k=0;k<source->GetNPixels();k++){
			SourcePB::Pixel* thisPixel = blob.add_m_pixels();
			if(EncodePixelToProtobuf(*thisPixel,source->GetPixel(k))<0){
				std::stringstream errMsg;
				errMsg<<"Nested source no. "<<i+1<<" encoding to protobuf failed!";
				throw std::runtime_error(errMsg.str().c_str());
			}
		}
		source_pb.set_allocated_blob(blob);
		*/

		
	}//close try block
	catch(std::exception const & e) {
		ERROR_LOG("Source encoding to protobuf failed with status "<<e.what());
		return -1;
	}

	return 0;

}//close EncodeSourceToProtobuf()


int Serializer::SourceToBuffer(SBuffer& buffer,Source* source){

	//## Check input source
	if(!source) return -1;

	try {
		// Create google protobuf source message
		SourcePB::Source source_pb;
		if(EncodeSourceToProtobuf(source_pb,source)<0){
			throw std::runtime_error("Encoding failed!");
		}
		
		
		
		/*
		//Create header message
		MessagePB::Header* header= new MessagePB::Header;
  	header->set_name(msg.name);
		header->set_source(msg.source);
		header->set_id(msg.id);
		header->set_nmsg(msg.nmsg);
		header->set_timestamp(msg.time);

		short type= msg.type;
		if(type==Message::eREQUEST){
			header->set_type(MessagePB::Header::REQUEST);

			//Create and fill request message
			MessagePB::Request* request= new MessagePB::Request;
				
			//for (std::map<std::string,double>::iterator it=msg.args.begin();it!=msg.args.end();++it){
			for (std::map<std::string,float>::iterator it=msg.args.begin();it!=msg.args.end();++it){
      	std::string arg_name= it->first;
				//double arg_value= it->second;
				float arg_value= it->second;
				MessagePB::Argument* thisArg = request->add_args();
    		thisArg->set_name(arg_name);
    		thisArg->set_value(arg_value);
    	}//end loop args
			for (std::map<std::string,std::string>::iterator it=msg.options.begin();it!=msg.options.end();++it){
      	std::string option_name= it->first;
				std::string option_value= it->second;
      	MessagePB::Option* thisOption = request->add_options();
    		thisOption->set_name(option_name);
    		thisOption->set_value(option_value);
			}//end loop options
			for (std::map<std::string,Message::DataCollection>::iterator it=msg.data.begin();it!=msg.data.end();++it){
      	std::string data_name= it->first;
				Message::DataCollection data_values= it->second;
			
				MessagePB::Data* thisDataArray = request->add_data();
				thisDataArray->set_name(data_name);
				for(unsigned int j=0;j<data_values.size();j++){
					thisDataArray->add_value(data_values[j]);	
				}//end loop data values
			}//end loop data
	
			msg_pb->set_allocated_request(request);	

		}//close if REQUEST
		else if(type==Message::eREPLY){
			header->set_type(MessagePB::Header::REPLY);

			//Create and fill reply message
			MessagePB::Reply* reply= new MessagePB::Reply;
			reply->set_info(msg.info);

			short ack= msg.ack;	
			if(ack==Message::eOK) reply->set_ack(MessagePB::Reply::OK);
			else if(ack==Message::eFAIL) reply->set_ack(MessagePB::Reply::FAIL);
			else if(ack==Message::eINVALID) reply->set_ack(MessagePB::Reply::INVALID);
			else {
				std::stringstream errMsg;
   	 		errMsg << "ERROR: Invalid ack code (ack="<<ack<<")";
    		throw std::invalid_argument(errMsg.str().c_str());
			}
			
			for (std::map<std::string,Message::DataCollection>::iterator it=msg.data.begin();it!=msg.data.end();++it){
      	std::string data_name= it->first;
				Message::DataCollection data_values= it->second;
			
				MessagePB::Data* thisDataArray = reply->add_data();
				thisDataArray->set_name(data_name);
				for(unsigned int j=0;j<data_values.size();j++){
					thisDataArray->add_value(data_values[j]);	
				}//end loop data values
			}//end loop data

			msg_pb->set_allocated_reply(reply);	

		}//close if REPLY

		else if(type==Message::eEVENT){
			header->set_type(MessagePB::Header::EVENT);

			MessagePB::Event* event= new MessagePB::Event;
			event->set_topic(msg.topic);
			event->set_info(msg.info);
			event->set_action(msg.action);
		
			int severity_level= msg.severity_level;
			if(severity_level==Message::eLOG) event->set_level(MessagePB::Event::LOG);
			else if(severity_level==Message::eALARM) event->set_level(MessagePB::Event::ALARM);
			else if(severity_level==Message::eCRITICAL) event->set_level(MessagePB::Event::CRITICAL);
			else{
				std::stringstream errMsg;
   	 		errMsg << "ERROR: Invalid severity level (level="<<severity_level<<")";
    		throw std::invalid_argument(errMsg.str().c_str());
			}

			for (std::map<std::string,Message::DataCollection>::iterator it=msg.data.begin();it!=msg.data.end();++it){
      	std::string data_name= it->first;
				Message::DataCollection data_values= it->second;
			
				MessagePB::Data* thisDataArray = event->add_data();
				thisDataArray->set_name(data_name);
				for(unsigned int j=0;j<data_values.size();j++){
					thisDataArray->add_value(data_values[j]);	
				}//end loop data values
			}//end loop data

			msg_pb->set_allocated_event(event);	

		}//close if EVENT

		msg_pb->set_allocated_header(header);	
		*/

		//## Fill buffer 
		buffer.size = source_pb.ByteSize();
		buffer.data= source_pb.SerializeAsString();		
		
	}//close try block
	catch(std::exception const & e) {
		ERROR_LOG("Source encoding failed with status "<<e.what());
		return -1;
	}

	return 0;

}//close SourceToBuffer()


int Serializer::SourceToBuffer(Source* source,msgpack::sbuffer& buffer){

	try{	
		msgpack::pack(&buffer,*source);
	}
	catch(std::exception const & e) {
		ERROR_LOG("Source encoding failed (err: "<<e.what()<<")");
		return -1;
	}

	return 0;

}//close SourceToBuffer()

int Serializer::SourceToString(Source* source,std::string& msg){

	msgpack::sbuffer buffer;
	if(SourceToBuffer(source,buffer)<0) return -1;

	msg= std::string(buffer.data(),buffer.size());

	return 0;

}//close SourceToString()


int Serializer::SourceToDevString(Source* source,Tango::DevString& msg){

	msgpack::sbuffer buffer;
	if(SourceToBuffer(source,buffer)<0)	return -1;

	msg= CORBA::string_alloc(buffer.size());
	memcpy (msg, buffer.data(), buffer.size());	
	msg[buffer.size()]= 0;

	return 0;

}//close SourceToDevString()

int Serializer::SourceCollectionToBuffer(std::vector<Source*>& sources,msgpack::sbuffer& buffer){

	try{	
		msgpack::pack(&buffer,sources);
	}
	catch(std::exception const & e) {
		ERROR_LOG("Source encoding failed (err: "<<e.what()<<")");
		return -1;
	}

	return 0;

}//close SourceCollectionToBuffer()

int Serializer::SourceCollectionToString(std::vector<Source*>& sources,std::string& msg){

	msgpack::sbuffer buffer;
	if(SourceCollectionToBuffer(sources,buffer)<0) return -1;

	msg= std::string(buffer.data(),buffer.size());

	return 0;

}//close SourceCollectionToString()

int Serializer::SourceCollectionToDevString(std::vector<Source*>& sources,Tango::DevString& msg){

	msgpack::sbuffer buffer;
	if(SourceCollectionToBuffer(sources,buffer)<0)	return -1;

	msg= CORBA::string_alloc(buffer.size());
	memcpy (msg, buffer.data(), buffer.size());	
	msg[buffer.size()]= 0;

	return 0;

}//close SourceCollectionToDevString()

//#endif

/*
int Serializer::BufferToSource(msgpack::sbuffer& buffer){

	return 0;
}
*/

/*

bool MessageUtils::encodeFromMsgPack(SBuffer buffer,Message& msg){

	//## Check data integrity
	if(buffer.msg.empty() || buffer.size<=0) return false;

	//## Encode to message
	try {
		msgpack::unpacked unpacked_msg;
  	msgpack::unpack(&unpacked_msg, (buffer.msg).c_str(), buffer.size);	
  	unpacked_msg.get().convert(&msg);
	
		//## Check message validity
		if(!msg.isValid()) return false;
	}
	catch(std::exception const & e) {
		cout << "ERROR: message encoding failed with status "<<e.what() <<endl;
		return false;
	}

	return true;

}//close MessageUtils::encodeFromMsgPack()



bool MessageUtils::encodeFromMsgPackCollection(SBufferCollection buffers,MessageCollection& msgs){

	bool status= true;
	msgs.clear();
	msgs.resize(0);

	for(unsigned int i=0;i<buffers.size();i++){
		Message thisMsg;
		bool isValidDecoding= MessageUtils::encodeFromMsgPack(buffers[i],thisMsg);
		if(!isValidDecoding) status= false; 
		msgs.push_back(thisMsg);
	}//end loop buffers

	return status;

}//close MessageUtils::encodeFromMsgPackCollection()



bool MessageUtils::encodeFromMsgPack(char* buffer,int bufsize,Message& msg){

	//## Check data integrity
	if(!buffer || bufsize<=0) return false;

	//## Encode to message
	try {
		msgpack::unpacked unpacked_msg;
  	msgpack::unpack(&unpacked_msg, buffer, bufsize);	
  	unpacked_msg.get().convert(&msg);
	
		//## Check message validity
		if(!msg.isValid()) return false;
	}
	catch(std::exception const & e) {
		cout << "ERROR: message encoding failed with status "<<e.what() <<endl;
		return false;
	}

	return true;

}//close encodeFromMsgPack()

*/

}//close namespace
