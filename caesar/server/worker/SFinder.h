/*----- PROTECTED REGION ID(SFinder.h) ENABLED START -----*/
//=============================================================================
//
// file :        SFinder.h
//
// description : Include file for the SFinder class
//
// project :     Source finder worker
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


#ifndef SFinder_H
#define SFinder_H

#include <SFinderThread.h>

#include <tango.h>

//Jsoncpp headers
#include <json/json.h>
#include <json/reader.h>

//Caesar headers
#include <Img.h>
#include <BkgData.h>
#include <Source.h>
namespace Caesar {
	class Img;
	class BkgData;
	class Source;
}

/*----- PROTECTED REGION END -----*/	//	SFinder.h

/**
 *  SFinder class description:
 *    Device server for source finding.
 */

namespace SFinder_ns
{
/*----- PROTECTED REGION ID(SFinder::Additional Class Declarations) ENABLED START -----*/

//	Additional Class Declarations
	class SFinderThread;

/*----- PROTECTED REGION END -----*/	//	SFinder::Additional Class Declarations

class SFinder : public TANGO_BASE_CLASS
{

/*----- PROTECTED REGION ID(SFinder::Data Members) ENABLED START -----*/

//	Add your own data members
	protected:
		//Thread 
		bool m_StopThreadFlag;	
		omni_mutex* m_mutex;
		SFinderThread* m_WorkerThread;

		//Pipes/Blobs
		Tango::DevicePipeBlob m_compactSourcesPipeBlob;
		
		//Broker group
		Tango::Group* m_brokerGroup;

		
		//Bkg attributes
		Tango::DevBoolean	attr_useLocalBkg_write;	
		Tango::DevBoolean	attr_use2ndPassInLocalBkg_write;
		Tango::DevBoolean attr_skipOutliersInLocalBkg_write;
		Tango::DevBoolean attr_skipNegativePixels_write;
		Tango::DevShort attr_localBkgMethod_write;
		Tango::DevShort attr_bkgEstimator_write;
		Tango::DevBoolean attr_useBeamInfoInBkg_write;
		Tango::DevFloat attr_localBkgBoxSizeX_write;
		Tango::DevFloat attr_localBkgBoxSizeY_write;
		Tango::DevFloat attr_localBkgGridStepSizeX_write;
		Tango::DevFloat attr_localBkgGridStepSizeY_write;	

		//Source finding attributes
		Tango::DevBoolean attr_searchCompactSources_write;
		Tango::DevFloat	attr_seedThr_write;
		Tango::DevFloat	attr_mergeThr_write;
		Tango::DevLong attr_minNPix_write;
		Tango::DevBoolean	attr_mergeBelowSeed_write;
		Tango::DevBoolean	attr_searchNegativeExcess_write;

		//Nested source finding
		Tango::DevBoolean	attr_searchNestedSources_write;
		Tango::DevFloat	attr_nestedBlobThrFactor_write;

		//Source selection
		Tango::DevBoolean	attr_selectCompactSources_write;
		Tango::DevBoolean	attr_useCircRatioCut_write;
		Tango::DevFloat	attr_psCircRatioThr_write;
		Tango::DevBoolean	attr_useElongCut_write;
		Tango::DevFloat	attr_psElongThr_write;
		Tango::DevBoolean	attr_useEllipseAreaRatioCut_write;
		Tango::DevFloat	attr_psEllipseAreaRatioMinThr_write;
		Tango::DevFloat	attr_psEllipseAreaRatioMaxThr_write;
		Tango::DevBoolean	attr_useMaxNPixCut_write;
		Tango::DevLong attr_psMaxNPix_write;
		Tango::DevBoolean	attr_useBoundingBoxCut_write;
		Tango::DevFloat	attr_minBoundingBoxThr_write;

/*----- PROTECTED REGION END -----*/	//	SFinder::Data Members

//	Device property data members
public:
	//	useLocalBkg_default:	Default value for useLocalBkg property
	Tango::DevBoolean	useLocalBkg_default;
	//	use2ndPassInLocalBkg_default:	Use 2nd pass when computing the local bkg
	Tango::DevBoolean	use2ndPassInLocalBkg_default;
	//	skipNegativePixels_default:	Default flag to skip/keep negative pixels when computing 
	//  image stats
	Tango::DevBoolean	skipNegativePixels_default;
	//	skipOutliersInLocalBkg_default:	Default value of flag to enable/disable skipping outliers in local
	//  bkg computation
	Tango::DevBoolean	skipOutliersInLocalBkg_default;
	//	localBkgMethod_default:	Default value of method to be used to compute local bkg 
	//  (1=Grid, 2=Superpixel)
	Tango::DevShort	localBkgMethod_default;
	//	bkgEstimator_default:	Default bkg estimator 
	//  (1=Mean, 2=Median, 3=BiWeight, 4=MedianClipped)
	Tango::DevShort	bkgEstimator_default;
	//	useBeamInfoInBkg_default:	Default flag for using beam information in bkg computation
	Tango::DevBoolean	useBeamInfoInBkg_default;
	//	localBkgBoxSizeX_default:	Default box size X used for local bkg computation. 
	//  If beam information is used, this corresponds to a
	//  multiple of beam size (e.g. x10, 20, 30 typically), otherwise it 
	//  corresponds to a fraction of image size.
	Tango::DevFloat	localBkgBoxSizeX_default;
	//	localBkgBoxSizeY_default:	Default box size Y used for local bkg computation. 
	//  If beam information is used, this corresponds to a
	//  multiple of beam size (e.g. x10, 20, 30 typically), otherwise it 
	//  corresponds to a fraction of image size.
	Tango::DevFloat	localBkgBoxSizeY_default;
	//	localBkgGridStepSizeX_default:	Default grid step size X used in local bkg computation. 
	//  This corresponds to a fraction of the box size X (e.g. 0.2, 0.5 
	//  are typical values)
	Tango::DevFloat	localBkgGridStepSizeX_default;
	//	localBkgGridStepSizeY_default:	Default grid step size X used in local bkg computation. 
	//  This corresponds to a fraction of the box size X (e.g. 0.2, 0.5 
	//  are typical values)
	Tango::DevFloat	localBkgGridStepSizeY_default;
	//	seedThr_default:	Default seed threshold (in number of sigmas above significance 
	//  level) to be used in source finding (flood-fill).
	Tango::DevFloat	seedThr_default;
	//	mergeThr_default:	Default merge threshold (in number of sigmas above significance 
	//  level) to be used in source finding (flood-fill).
	Tango::DevFloat	mergeThr_default;
	//	minNPix_default:	Default minimum number of pixels in source fnding (blob size in
	//  flood-fill).
	Tango::DevLong	minNPix_default;
	//	mergeBelowSeed_default:	Default flag value to aggregate only pixels below seed 
	//  threshold in flood-fill
	Tango::DevBoolean	mergeBelowSeed_default;
	//	searchNegativeExcess_default:	Default flag to search negative excess together with positive 
	//  in compact source search
	Tango::DevBoolean	searchNegativeExcess_default;
	//	searchNestedSources_default:	Default flag to search for nested sources in compact source search.
	Tango::DevBoolean	searchNestedSources_default;
	//	nestedBlobThrFactor_default:	Default Threshold (in multiple of curvature RMS) to be used 
	//  to detect nested blobs
	Tango::DevFloat	nestedBlobThrFactor_default;
	//	searchCompactSources_default:	Default flag to enable/disable search of compact sources
	Tango::DevBoolean	searchCompactSources_default;
	//	selectCompactSources_default:	Default flag to enable/disable selection of compact source 
	//  on the basis of defined cuts
	Tango::DevBoolean	selectCompactSources_default;
	//	useCircRatioCut_default:	Use cut on circularity ratio in source selection
	Tango::DevBoolean	useCircRatioCut_default;
	//	psCircRatioThr_default:	Point-source circularity ratio default cut. 
	//  Source is not selected as point-like if circratio<cut (circ ratio
	//  is =1 for a circle).
	Tango::DevFloat	psCircRatioThr_default;
	//	useElongCut_default:	Default flag to enable/disable elongation point-like source cut.
	//  Source is not selected as point-source if elong>cut (Elong=0 
	//  for circle/square).
	Tango::DevBoolean	useElongCut_default;
	//	psElongThr_default:	Default elongation point-like source cut
	Tango::DevFloat	psElongThr_default;
	//	useEllipseAreaRatioCut_default:	Default flag to enable/disable ellipse area ratio cut
	Tango::DevBoolean	useEllipseAreaRatioCut_default;
	//	psEllipseAreaRatioMinThr_default:	Default min Ellipse Area Ratio cut
	Tango::DevFloat	psEllipseAreaRatioMinThr_default;
	//	psEllipseAreaRatioMaxThr_default:	Default max Ellipse Area Ratio cut
	Tango::DevFloat	psEllipseAreaRatioMaxThr_default;
	//	useMaxNPixCut_default:	Default flag to enable/disable max npixel point-like source cut
	Tango::DevBoolean	useMaxNPixCut_default;
	//	psMaxNPix_default:	Default maximum number of pixels to select point-like sources
	Tango::DevLong	psMaxNPix_default;
	//	useBoundingBoxCut_default:	Default flag to enable/disable bounding box cut
	Tango::DevBoolean	useBoundingBoxCut_default;
	//	minBoundingBoxThr_default:	Minimum default bounding box cut (source tagged as bad if below this threshold)
	Tango::DevFloat	minBoundingBoxThr_default;
	//	brokerList:	List of broker devices
	vector<string>	brokerList;

//	Attribute data members
public:
	Tango::DevBoolean	*attr_useLocalBkg_read;
	Tango::DevBoolean	*attr_use2ndPassInLocalBkg_read;
	Tango::DevBoolean	*attr_skipNegativePixels_read;
	Tango::DevBoolean	*attr_skipOutliersInLocalBkg_read;
	Tango::DevShort	*attr_localBkgMethod_read;
	Tango::DevShort	*attr_bkgEstimator_read;
	Tango::DevBoolean	*attr_useBeamInfoInBkg_read;
	Tango::DevFloat	*attr_localBkgBoxSizeX_read;
	Tango::DevFloat	*attr_localBkgBoxSizeY_read;
	Tango::DevFloat	*attr_localBkgGridStepSizeX_read;
	Tango::DevFloat	*attr_localBkgGridStepSizeY_read;
	Tango::DevFloat	*attr_seedThr_read;
	Tango::DevFloat	*attr_mergeThr_read;
	Tango::DevLong	*attr_minNPix_read;
	Tango::DevBoolean	*attr_mergeBelowSeed_read;
	Tango::DevBoolean	*attr_searchNegativeExcess_read;
	Tango::DevFloat	*attr_nestedBlobThrFactor_read;
	Tango::DevBoolean	*attr_searchCompactSources_read;
	Tango::DevBoolean	*attr_selectCompactSources_read;
	Tango::DevBoolean	*attr_useCircRatioCut_read;
	Tango::DevFloat	*attr_psCircRatioThr_read;
	Tango::DevBoolean	*attr_useElongCut_read;
	Tango::DevFloat	*attr_psElongThr_read;
	Tango::DevBoolean	*attr_useEllipseAreaRatioCut_read;
	Tango::DevFloat	*attr_psEllipseAreaRatioMinThr_read;
	Tango::DevFloat	*attr_psEllipseAreaRatioMaxThr_read;
	Tango::DevBoolean	*attr_useMaxNPixCut_read;
	Tango::DevLong	*attr_psMaxNPix_read;
	Tango::DevBoolean	*attr_useBoundingBoxCut_read;
	Tango::DevFloat	*attr_minBoundingBoxThr_read;
	Tango::DevString	*attr_compactSourceData_read;
	Tango::DevString	*attr_runProgress_read;

//	Constructors and destructors
public:
	/**
	 * Constructs a newly device object.
	 *
	 *	@param cl	Class.
	 *	@param s 	Device Name
	 */
	SFinder(Tango::DeviceClass *cl,string &s);
	/**
	 * Constructs a newly device object.
	 *
	 *	@param cl	Class.
	 *	@param s 	Device Name
	 */
	SFinder(Tango::DeviceClass *cl,const char *s);
	/**
	 * Constructs a newly device object.
	 *
	 *	@param cl	Class.
	 *	@param s 	Device name
	 *	@param d	Device description.
	 */
	SFinder(Tango::DeviceClass *cl,const char *s,const char *d);
	/**
	 * The device object destructor.
	 */
	~SFinder() {delete_device();};


//	Miscellaneous methods
public:
	/*
	 *	will be called at device destruction or at init command.
	 */
	void delete_device();
	/*
	 *	Initialize the device
	 */
	virtual void init_device();
	/*
	 *	Read the device properties from database
	 */
	void get_device_property();
	/*
	 *	Always executed method before execution command method.
	 */
	virtual void always_executed_hook();


//	Attribute methods
public:
	//--------------------------------------------------------
	/*
	 *	Method      : SFinder::read_attr_hardware()
	 *	Description : Hardware acquisition for attributes.
	 */
	//--------------------------------------------------------
	virtual void read_attr_hardware(vector<long> &attr_list);
	//--------------------------------------------------------
	/*
	 *	Method      : SFinder::write_attr_hardware()
	 *	Description : Hardware writing for attributes.
	 */
	//--------------------------------------------------------
	virtual void write_attr_hardware(vector<long> &attr_list);

/**
 *	Attribute useLocalBkg related methods
 *	Description: 
 *
 *	Data type:	Tango::DevBoolean
 *	Attr type:	Scalar
 */
	virtual void read_useLocalBkg(Tango::Attribute &attr);
	virtual void write_useLocalBkg(Tango::WAttribute &attr);
	virtual bool is_useLocalBkg_allowed(Tango::AttReqType type);
/**
 *	Attribute use2ndPassInLocalBkg related methods
 *	Description: 
 *
 *	Data type:	Tango::DevBoolean
 *	Attr type:	Scalar
 */
	virtual void read_use2ndPassInLocalBkg(Tango::Attribute &attr);
	virtual void write_use2ndPassInLocalBkg(Tango::WAttribute &attr);
	virtual bool is_use2ndPassInLocalBkg_allowed(Tango::AttReqType type);
/**
 *	Attribute skipNegativePixels related methods
 *	Description: Flag to skip/keep negative pixels when computing the image stats
 *
 *	Data type:	Tango::DevBoolean
 *	Attr type:	Scalar
 */
	virtual void read_skipNegativePixels(Tango::Attribute &attr);
	virtual void write_skipNegativePixels(Tango::WAttribute &attr);
	virtual bool is_skipNegativePixels_allowed(Tango::AttReqType type);
/**
 *	Attribute skipOutliersInLocalBkg related methods
 *	Description: Flag to enable/disable skipping outliers in local bkg computation
 *
 *	Data type:	Tango::DevBoolean
 *	Attr type:	Scalar
 */
	virtual void read_skipOutliersInLocalBkg(Tango::Attribute &attr);
	virtual void write_skipOutliersInLocalBkg(Tango::WAttribute &attr);
	virtual bool is_skipOutliersInLocalBkg_allowed(Tango::AttReqType type);
/**
 *	Attribute localBkgMethod related methods
 *	Description: Method used to compute local bkg (1=Grid, 2=Superpixel)
 *
 *	Data type:	Tango::DevShort
 *	Attr type:	Scalar
 */
	virtual void read_localBkgMethod(Tango::Attribute &attr);
	virtual void write_localBkgMethod(Tango::WAttribute &attr);
	virtual bool is_localBkgMethod_allowed(Tango::AttReqType type);
/**
 *	Attribute bkgEstimator related methods
 *	Description: Estimator to be used to compute bkg 
 *               (1=mean, 2=median, 3=biweight, 4=clipped median)
 *
 *	Data type:	Tango::DevShort
 *	Attr type:	Scalar
 */
	virtual void read_bkgEstimator(Tango::Attribute &attr);
	virtual void write_bkgEstimator(Tango::WAttribute &attr);
	virtual bool is_bkgEstimator_allowed(Tango::AttReqType type);
/**
 *	Attribute useBeamInfoInBkg related methods
 *	Description: Use beam information in bkg computation
 *
 *	Data type:	Tango::DevBoolean
 *	Attr type:	Scalar
 */
	virtual void read_useBeamInfoInBkg(Tango::Attribute &attr);
	virtual void write_useBeamInfoInBkg(Tango::WAttribute &attr);
	virtual bool is_useBeamInfoInBkg_allowed(Tango::AttReqType type);
/**
 *	Attribute localBkgBoxSizeX related methods
 *	Description: Box size X used for local bkg computation. If beam information is used, this corresponds to 
 *               a multiple of beam size (e.g. x10, 20, 30 typically), otherwise it corresponds to a fraction of 
 *               image size.
 *
 *	Data type:	Tango::DevFloat
 *	Attr type:	Scalar
 */
	virtual void read_localBkgBoxSizeX(Tango::Attribute &attr);
	virtual void write_localBkgBoxSizeX(Tango::WAttribute &attr);
	virtual bool is_localBkgBoxSizeX_allowed(Tango::AttReqType type);
/**
 *	Attribute localBkgBoxSizeY related methods
 *	Description: Box size Y used for local bkg computation. If beam information is used, this corresponds to a
 *               multiple of beam size (e.g. x10, 20, 30 typically), otherwise it corresponds to a fraction of 
 *               image size.
 *
 *	Data type:	Tango::DevFloat
 *	Attr type:	Scalar
 */
	virtual void read_localBkgBoxSizeY(Tango::Attribute &attr);
	virtual void write_localBkgBoxSizeY(Tango::WAttribute &attr);
	virtual bool is_localBkgBoxSizeY_allowed(Tango::AttReqType type);
/**
 *	Attribute localBkgGridStepSizeX related methods
 *	Description: Grid step size X used in local bkg computation. This corresponds to a fraction of the box 
 *               size X (e.g. 0.2, 0.5 are typical values)
 *
 *	Data type:	Tango::DevFloat
 *	Attr type:	Scalar
 */
	virtual void read_localBkgGridStepSizeX(Tango::Attribute &attr);
	virtual void write_localBkgGridStepSizeX(Tango::WAttribute &attr);
	virtual bool is_localBkgGridStepSizeX_allowed(Tango::AttReqType type);
/**
 *	Attribute localBkgGridStepSizeY related methods
 *	Description: Grid step size Y used in local bkg computation. This corresponds to a fraction of the box size Y 
 *               (e.g. 0.2, 0.5 are typical values)
 *
 *	Data type:	Tango::DevFloat
 *	Attr type:	Scalar
 */
	virtual void read_localBkgGridStepSizeY(Tango::Attribute &attr);
	virtual void write_localBkgGridStepSizeY(Tango::WAttribute &attr);
	virtual bool is_localBkgGridStepSizeY_allowed(Tango::AttReqType type);
/**
 *	Attribute seedThr related methods
 *	Description: Seed threshold (in number of sigmas above significance level) to be used in source 
 *               finding (flood-fill).
 *
 *	Data type:	Tango::DevFloat
 *	Attr type:	Scalar
 */
	virtual void read_seedThr(Tango::Attribute &attr);
	virtual void write_seedThr(Tango::WAttribute &attr);
	virtual bool is_seedThr_allowed(Tango::AttReqType type);
/**
 *	Attribute mergeThr related methods
 *	Description: Merge threshold (in number of sigmas above significance level) to be used in source 
 *               finding (flood-fill).
 *
 *	Data type:	Tango::DevFloat
 *	Attr type:	Scalar
 */
	virtual void read_mergeThr(Tango::Attribute &attr);
	virtual void write_mergeThr(Tango::WAttribute &attr);
	virtual bool is_mergeThr_allowed(Tango::AttReqType type);
/**
 *	Attribute minNPix related methods
 *	Description: Minimum number of pixels in source fnding (blob size in flood-fill).
 *
 *	Data type:	Tango::DevLong
 *	Attr type:	Scalar
 */
	virtual void read_minNPix(Tango::Attribute &attr);
	virtual void write_minNPix(Tango::WAttribute &attr);
	virtual bool is_minNPix_allowed(Tango::AttReqType type);
/**
 *	Attribute mergeBelowSeed related methods
 *	Description: Flag value to aggregate only pixels below seed threshold in flood-fill
 *
 *	Data type:	Tango::DevBoolean
 *	Attr type:	Scalar
 */
	virtual void read_mergeBelowSeed(Tango::Attribute &attr);
	virtual void write_mergeBelowSeed(Tango::WAttribute &attr);
	virtual bool is_mergeBelowSeed_allowed(Tango::AttReqType type);
/**
 *	Attribute searchNegativeExcess related methods
 *	Description: Flag to search negative excess together with positive in compact source search
 *
 *	Data type:	Tango::DevBoolean
 *	Attr type:	Scalar
 */
	virtual void read_searchNegativeExcess(Tango::Attribute &attr);
	virtual void write_searchNegativeExcess(Tango::WAttribute &attr);
	virtual bool is_searchNegativeExcess_allowed(Tango::AttReqType type);
/**
 *	Attribute nestedBlobThrFactor related methods
 *	Description: Threshold (in multiple of curvature RMS) to be used to detect nested blobs
 *
 *	Data type:	Tango::DevFloat
 *	Attr type:	Scalar
 */
	virtual void read_nestedBlobThrFactor(Tango::Attribute &attr);
	virtual void write_nestedBlobThrFactor(Tango::WAttribute &attr);
	virtual bool is_nestedBlobThrFactor_allowed(Tango::AttReqType type);
/**
 *	Attribute searchCompactSources related methods
 *	Description: Flag to enable/disable search of compact sources
 *
 *	Data type:	Tango::DevBoolean
 *	Attr type:	Scalar
 */
	virtual void read_searchCompactSources(Tango::Attribute &attr);
	virtual void write_searchCompactSources(Tango::WAttribute &attr);
	virtual bool is_searchCompactSources_allowed(Tango::AttReqType type);
/**
 *	Attribute selectCompactSources related methods
 *	Description: Flag to enable/disable selection of compact sources according to defined cuts
 *
 *	Data type:	Tango::DevBoolean
 *	Attr type:	Scalar
 */
	virtual void read_selectCompactSources(Tango::Attribute &attr);
	virtual void write_selectCompactSources(Tango::WAttribute &attr);
	virtual bool is_selectCompactSources_allowed(Tango::AttReqType type);
/**
 *	Attribute useCircRatioCut related methods
 *	Description: Use cut on circularity ratio in source selection
 *
 *	Data type:	Tango::DevBoolean
 *	Attr type:	Scalar
 */
	virtual void read_useCircRatioCut(Tango::Attribute &attr);
	virtual void write_useCircRatioCut(Tango::WAttribute &attr);
	virtual bool is_useCircRatioCut_allowed(Tango::AttReqType type);
/**
 *	Attribute psCircRatioThr related methods
 *	Description: Point-source circularity ratio cut. 
 *               Source is not selected as point-like if circratio<cut (circ ratio
 *               is =1 for a circle).
 *
 *	Data type:	Tango::DevFloat
 *	Attr type:	Scalar
 */
	virtual void read_psCircRatioThr(Tango::Attribute &attr);
	virtual void write_psCircRatioThr(Tango::WAttribute &attr);
	virtual bool is_psCircRatioThr_allowed(Tango::AttReqType type);
/**
 *	Attribute useElongCut related methods
 *	Description: Flag to enable/disable elongation point-like source cut.
 *               Source is not selected as point-source if elong>cut (Elong=0 
 *               for circle/square).
 *
 *	Data type:	Tango::DevBoolean
 *	Attr type:	Scalar
 */
	virtual void read_useElongCut(Tango::Attribute &attr);
	virtual void write_useElongCut(Tango::WAttribute &attr);
	virtual bool is_useElongCut_allowed(Tango::AttReqType type);
/**
 *	Attribute psElongThr related methods
 *	Description: Elongation point-like source cut
 *
 *	Data type:	Tango::DevFloat
 *	Attr type:	Scalar
 */
	virtual void read_psElongThr(Tango::Attribute &attr);
	virtual void write_psElongThr(Tango::WAttribute &attr);
	virtual bool is_psElongThr_allowed(Tango::AttReqType type);
/**
 *	Attribute useEllipseAreaRatioCut related methods
 *	Description:  Flag to enable/disable ellipse area ratio cut
 *
 *	Data type:	Tango::DevBoolean
 *	Attr type:	Scalar
 */
	virtual void read_useEllipseAreaRatioCut(Tango::Attribute &attr);
	virtual void write_useEllipseAreaRatioCut(Tango::WAttribute &attr);
	virtual bool is_useEllipseAreaRatioCut_allowed(Tango::AttReqType type);
/**
 *	Attribute psEllipseAreaRatioMinThr related methods
 *	Description: Min Ellipse Area Ratio cut
 *
 *	Data type:	Tango::DevFloat
 *	Attr type:	Scalar
 */
	virtual void read_psEllipseAreaRatioMinThr(Tango::Attribute &attr);
	virtual void write_psEllipseAreaRatioMinThr(Tango::WAttribute &attr);
	virtual bool is_psEllipseAreaRatioMinThr_allowed(Tango::AttReqType type);
/**
 *	Attribute psEllipseAreaRatioMaxThr related methods
 *	Description: Max Ellipse Area Ratio cut
 *
 *	Data type:	Tango::DevFloat
 *	Attr type:	Scalar
 */
	virtual void read_psEllipseAreaRatioMaxThr(Tango::Attribute &attr);
	virtual void write_psEllipseAreaRatioMaxThr(Tango::WAttribute &attr);
	virtual bool is_psEllipseAreaRatioMaxThr_allowed(Tango::AttReqType type);
/**
 *	Attribute useMaxNPixCut related methods
 *	Description: 
 *
 *	Data type:	Tango::DevBoolean
 *	Attr type:	Scalar
 */
	virtual void read_useMaxNPixCut(Tango::Attribute &attr);
	virtual void write_useMaxNPixCut(Tango::WAttribute &attr);
	virtual bool is_useMaxNPixCut_allowed(Tango::AttReqType type);
/**
 *	Attribute psMaxNPix related methods
 *	Description: Maximum number of pixels to select point-like sources
 *
 *	Data type:	Tango::DevLong
 *	Attr type:	Scalar
 */
	virtual void read_psMaxNPix(Tango::Attribute &attr);
	virtual void write_psMaxNPix(Tango::WAttribute &attr);
	virtual bool is_psMaxNPix_allowed(Tango::AttReqType type);
/**
 *	Attribute useBoundingBoxCut related methods
 *	Description: Flag to enable/disable bounding box cut
 *
 *	Data type:	Tango::DevBoolean
 *	Attr type:	Scalar
 */
	virtual void read_useBoundingBoxCut(Tango::Attribute &attr);
	virtual void write_useBoundingBoxCut(Tango::WAttribute &attr);
	virtual bool is_useBoundingBoxCut_allowed(Tango::AttReqType type);
/**
 *	Attribute minBoundingBoxThr related methods
 *	Description: Minimum default bounding box cut (source tagged as bad if below this threshold)
 *
 *	Data type:	Tango::DevFloat
 *	Attr type:	Scalar
 */
	virtual void read_minBoundingBoxThr(Tango::Attribute &attr);
	virtual void write_minBoundingBoxThr(Tango::WAttribute &attr);
	virtual bool is_minBoundingBoxThr_allowed(Tango::AttReqType type);
/**
 *	Attribute compactSourceData related methods
 *	Description: 
 *
 *	Data type:	Tango::DevString
 *	Attr type:	Scalar
 */
	virtual void read_compactSourceData(Tango::Attribute &attr);
	virtual bool is_compactSourceData_allowed(Tango::AttReqType type);
/**
 *	Attribute runProgress related methods
 *	Description: Run progress info
 *               [0]: run id
 *               [1]: status (RUNNING, COMPLETED, ABORTED, FAILED)
 *               [2]: progress fraction (0-100%)
 *               [3]: log info (INIT, READ_IMAGE, BKG, COMPACT_SOURCE, EXT_SOURCE)
 *               [4]: timestamp
 *
 *	Data type:	Tango::DevString
 *	Attr type:	Spectrum max = 10
 */
	virtual void read_runProgress(Tango::Attribute &attr);
	virtual bool is_runProgress_allowed(Tango::AttReqType type);


	//--------------------------------------------------------
	/**
	 *	Method      : SFinder::add_dynamic_attributes()
	 *	Description : Add dynamic attributes if any.
	 */
	//--------------------------------------------------------
	void add_dynamic_attributes();



//	pipe related methods
public:
	//	Pipe compactSourcesPipe
	bool is_compactSourcesPipe_allowed(Tango::PipeReqType);
	void read_compactSourcesPipe(Tango::Pipe &);

//	Command related methods
public:
	/**
	 *	Command ExtractSources related method
	 *	Description: Find sources in map using the configuration passed 
	 *               as argument.
	 *
	 *	@param argin String arg
	 *               [0]: filename
	 *               [1]: run guid (set by the broker)
	 *               [2]: configuration string  
	 *               
	 *               Long arg
	 *               [nmaps+0]: tile min x
	 *               [nmaps+1]: tile max x
	 *               [nmaps+2]: tile min y
	 *               [nmaps+3]: tile max y
	 *	@returns Long arg
	 *           [0]: ack code
	 *           
	 *           String arg
	 *           [0]: err description
	 */
	virtual Tango::DevVarLongStringArray *extract_sources(const Tango::DevVarLongStringArray *argin);
	virtual bool is_ExtractSources_allowed(const CORBA::Any &any);
	/**
	 *	Command Configure related method
	 *	Description: 
	 *
	 *	@param argin Configuration string
	 *	@returns 
	 */
	virtual Tango::DevVarLongStringArray *configure(Tango::DevString argin);
	virtual bool is_Configure_allowed(const CORBA::Any &any);
	/**
	 *	Command RegisterMe related method
	 *	Description: Register worker in brokers
	 *
	 *	@returns 
	 */
	virtual Tango::DevVarLongStringArray *register_me();
	virtual bool is_RegisterMe_allowed(const CORBA::Any &any);


	//--------------------------------------------------------
	/**
	 *	Method      : SFinder::add_dynamic_commands()
	 *	Description : Add dynamic commands if any.
	 */
	//--------------------------------------------------------
	void add_dynamic_commands();

/*----- PROTECTED REGION ID(SFinder::Additional Method prototypes) ENABLED START -----*/

//	Additional Method prototypes
	protected:
		int InitBrokerGroup();
		int LoadDefaultConfig();
		int ApplyConfig(std::string& config);
		int SetAttrFromConfig(Json::Value& optionObj);
		int SetScalarAttrValue(Tango::WAttribute& attr,Json::Value& optionObj);
		int SetSpectrumAttrValue(Tango::WAttribute& attr,Json::Value& optionObj);

		/*
		int RunSourceTask(std::vector<Caesar::Source*>& sources,const std::string& filename,long int tileMinX=-1,long int tileMaxX=-1,long int tileMinY=-1,long int tileMaxY=-1);
		Caesar::Img* ReadImage(const std::string& filename,long int tileMinX=-1,long int tileMaxX=-1,long int tileMinY=-1,long int tileMaxY=-1);
		Caesar::BkgData* ComputeStatsAndBkg(Caesar::Img* img);
		int FindSources(std::vector<Caesar::Source*>& sources,Caesar::Img* inputImg,bool computeStatsAndBkg=true,Caesar::BkgData* inputBkgData=0);
		int FindCompactSources(std::vector<Caesar::Source*>& sources,Caesar::Img* inputImg,bool computeStatsAndBkg=true,Caesar::BkgData* inputBkgData=0);
		int SelectSources(std::vector<Caesar::Source*>& sources);
		bool IsGoodSource(Caesar::Source* aSource);
		bool IsPointLikeSource(Caesar::Source* aSource);
		*/
	friend class SFinderThread;

/*----- PROTECTED REGION END -----*/	//	SFinder::Additional Method prototypes
};

/*----- PROTECTED REGION ID(SFinder::Additional Classes Definitions) ENABLED START -----*/

//	Additional Classes Definitions

/*----- PROTECTED REGION END -----*/	//	SFinder::Additional Classes Definitions

}	//	End of namespace

#endif   //	SFinder_H
