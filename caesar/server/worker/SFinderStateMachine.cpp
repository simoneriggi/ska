/*----- PROTECTED REGION ID(SFinderStateMachine.cpp) ENABLED START -----*/
static const char *RcsId = "$Id:  $";
//=============================================================================
//
// file :        SFinderStateMachine.cpp
//
// description : State machine file for the SFinder class
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

#include <SFinder.h>

/*----- PROTECTED REGION END -----*/	//	SFinder::SFinderStateMachine.cpp

//================================================================
//  States   |  Description
//================================================================
//  INIT     |  Initializing state: used when reserved of under configuration options
//  ON       |  Default state when server is free of tasks (e.g. it can be reserved)
//  FAULT    |  Fault state, e.g. when internal error occurred or when task failed
//  RUNNING  |  State activated when a task is started


namespace SFinder_ns
{
//=================================================
//		Attributes Allowed Methods
//=================================================

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_useLocalBkg_allowed()
 *	Description : Execution allowed for useLocalBkg attribute
 */
//--------------------------------------------------------
bool SFinder::is_useLocalBkg_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for useLocalBkg attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::useLocalBkgStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::useLocalBkgStateAllowed_WRITE

	//	Not any excluded states for useLocalBkg attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::useLocalBkgStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::useLocalBkgStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_use2ndPassInLocalBkg_allowed()
 *	Description : Execution allowed for use2ndPassInLocalBkg attribute
 */
//--------------------------------------------------------
bool SFinder::is_use2ndPassInLocalBkg_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for use2ndPassInLocalBkg attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::use2ndPassInLocalBkgStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::use2ndPassInLocalBkgStateAllowed_WRITE

	//	Not any excluded states for use2ndPassInLocalBkg attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::use2ndPassInLocalBkgStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::use2ndPassInLocalBkgStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_skipNegativePixels_allowed()
 *	Description : Execution allowed for skipNegativePixels attribute
 */
//--------------------------------------------------------
bool SFinder::is_skipNegativePixels_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for skipNegativePixels attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::skipNegativePixelsStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::skipNegativePixelsStateAllowed_WRITE

	//	Not any excluded states for skipNegativePixels attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::skipNegativePixelsStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::skipNegativePixelsStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_skipOutliersInLocalBkg_allowed()
 *	Description : Execution allowed for skipOutliersInLocalBkg attribute
 */
//--------------------------------------------------------
bool SFinder::is_skipOutliersInLocalBkg_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for skipOutliersInLocalBkg attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::skipOutliersInLocalBkgStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::skipOutliersInLocalBkgStateAllowed_WRITE

	//	Not any excluded states for skipOutliersInLocalBkg attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::skipOutliersInLocalBkgStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::skipOutliersInLocalBkgStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_localBkgMethod_allowed()
 *	Description : Execution allowed for localBkgMethod attribute
 */
//--------------------------------------------------------
bool SFinder::is_localBkgMethod_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for localBkgMethod attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::localBkgMethodStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::localBkgMethodStateAllowed_WRITE

	//	Not any excluded states for localBkgMethod attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::localBkgMethodStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::localBkgMethodStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_bkgEstimator_allowed()
 *	Description : Execution allowed for bkgEstimator attribute
 */
//--------------------------------------------------------
bool SFinder::is_bkgEstimator_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for bkgEstimator attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::bkgEstimatorStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::bkgEstimatorStateAllowed_WRITE

	//	Not any excluded states for bkgEstimator attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::bkgEstimatorStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::bkgEstimatorStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_useBeamInfoInBkg_allowed()
 *	Description : Execution allowed for useBeamInfoInBkg attribute
 */
//--------------------------------------------------------
bool SFinder::is_useBeamInfoInBkg_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for useBeamInfoInBkg attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::useBeamInfoInBkgStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::useBeamInfoInBkgStateAllowed_WRITE

	//	Not any excluded states for useBeamInfoInBkg attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::useBeamInfoInBkgStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::useBeamInfoInBkgStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_localBkgBoxSizeX_allowed()
 *	Description : Execution allowed for localBkgBoxSizeX attribute
 */
//--------------------------------------------------------
bool SFinder::is_localBkgBoxSizeX_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for localBkgBoxSizeX attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::localBkgBoxSizeXStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::localBkgBoxSizeXStateAllowed_WRITE

	//	Not any excluded states for localBkgBoxSizeX attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::localBkgBoxSizeXStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::localBkgBoxSizeXStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_localBkgBoxSizeY_allowed()
 *	Description : Execution allowed for localBkgBoxSizeY attribute
 */
//--------------------------------------------------------
bool SFinder::is_localBkgBoxSizeY_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for localBkgBoxSizeY attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::localBkgBoxSizeYStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::localBkgBoxSizeYStateAllowed_WRITE

	//	Not any excluded states for localBkgBoxSizeY attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::localBkgBoxSizeYStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::localBkgBoxSizeYStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_localBkgGridStepSizeX_allowed()
 *	Description : Execution allowed for localBkgGridStepSizeX attribute
 */
//--------------------------------------------------------
bool SFinder::is_localBkgGridStepSizeX_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for localBkgGridStepSizeX attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::localBkgGridStepSizeXStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::localBkgGridStepSizeXStateAllowed_WRITE

	//	Not any excluded states for localBkgGridStepSizeX attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::localBkgGridStepSizeXStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::localBkgGridStepSizeXStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_localBkgGridStepSizeY_allowed()
 *	Description : Execution allowed for localBkgGridStepSizeY attribute
 */
//--------------------------------------------------------
bool SFinder::is_localBkgGridStepSizeY_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for localBkgGridStepSizeY attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::localBkgGridStepSizeYStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::localBkgGridStepSizeYStateAllowed_WRITE

	//	Not any excluded states for localBkgGridStepSizeY attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::localBkgGridStepSizeYStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::localBkgGridStepSizeYStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_seedThr_allowed()
 *	Description : Execution allowed for seedThr attribute
 */
//--------------------------------------------------------
bool SFinder::is_seedThr_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for seedThr attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::seedThrStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::seedThrStateAllowed_WRITE

	//	Not any excluded states for seedThr attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::seedThrStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::seedThrStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_mergeThr_allowed()
 *	Description : Execution allowed for mergeThr attribute
 */
//--------------------------------------------------------
bool SFinder::is_mergeThr_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for mergeThr attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::mergeThrStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::mergeThrStateAllowed_WRITE

	//	Not any excluded states for mergeThr attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::mergeThrStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::mergeThrStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_minNPix_allowed()
 *	Description : Execution allowed for minNPix attribute
 */
//--------------------------------------------------------
bool SFinder::is_minNPix_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for minNPix attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::minNPixStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::minNPixStateAllowed_WRITE

	//	Not any excluded states for minNPix attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::minNPixStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::minNPixStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_mergeBelowSeed_allowed()
 *	Description : Execution allowed for mergeBelowSeed attribute
 */
//--------------------------------------------------------
bool SFinder::is_mergeBelowSeed_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for mergeBelowSeed attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::mergeBelowSeedStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::mergeBelowSeedStateAllowed_WRITE

	//	Not any excluded states for mergeBelowSeed attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::mergeBelowSeedStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::mergeBelowSeedStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_searchNegativeExcess_allowed()
 *	Description : Execution allowed for searchNegativeExcess attribute
 */
//--------------------------------------------------------
bool SFinder::is_searchNegativeExcess_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for searchNegativeExcess attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::searchNegativeExcessStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::searchNegativeExcessStateAllowed_WRITE

	//	Not any excluded states for searchNegativeExcess attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::searchNegativeExcessStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::searchNegativeExcessStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_nestedBlobThrFactor_allowed()
 *	Description : Execution allowed for nestedBlobThrFactor attribute
 */
//--------------------------------------------------------
bool SFinder::is_nestedBlobThrFactor_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for nestedBlobThrFactor attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::nestedBlobThrFactorStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::nestedBlobThrFactorStateAllowed_WRITE

	//	Not any excluded states for nestedBlobThrFactor attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::nestedBlobThrFactorStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::nestedBlobThrFactorStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_searchCompactSources_allowed()
 *	Description : Execution allowed for searchCompactSources attribute
 */
//--------------------------------------------------------
bool SFinder::is_searchCompactSources_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for searchCompactSources attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::searchCompactSourcesStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::searchCompactSourcesStateAllowed_WRITE

	//	Not any excluded states for searchCompactSources attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::searchCompactSourcesStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::searchCompactSourcesStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_selectCompactSources_allowed()
 *	Description : Execution allowed for selectCompactSources attribute
 */
//--------------------------------------------------------
bool SFinder::is_selectCompactSources_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for selectCompactSources attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::selectCompactSourcesStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::selectCompactSourcesStateAllowed_WRITE

	//	Not any excluded states for selectCompactSources attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::selectCompactSourcesStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::selectCompactSourcesStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_useCircRatioCut_allowed()
 *	Description : Execution allowed for useCircRatioCut attribute
 */
//--------------------------------------------------------
bool SFinder::is_useCircRatioCut_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for useCircRatioCut attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::useCircRatioCutStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::useCircRatioCutStateAllowed_WRITE

	//	Not any excluded states for useCircRatioCut attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::useCircRatioCutStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::useCircRatioCutStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_psCircRatioThr_allowed()
 *	Description : Execution allowed for psCircRatioThr attribute
 */
//--------------------------------------------------------
bool SFinder::is_psCircRatioThr_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for psCircRatioThr attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::psCircRatioThrStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::psCircRatioThrStateAllowed_WRITE

	//	Not any excluded states for psCircRatioThr attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::psCircRatioThrStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::psCircRatioThrStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_useElongCut_allowed()
 *	Description : Execution allowed for useElongCut attribute
 */
//--------------------------------------------------------
bool SFinder::is_useElongCut_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for useElongCut attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::useElongCutStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::useElongCutStateAllowed_WRITE

	//	Not any excluded states for useElongCut attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::useElongCutStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::useElongCutStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_psElongThr_allowed()
 *	Description : Execution allowed for psElongThr attribute
 */
//--------------------------------------------------------
bool SFinder::is_psElongThr_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for psElongThr attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::psElongThrStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::psElongThrStateAllowed_WRITE

	//	Not any excluded states for psElongThr attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::psElongThrStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::psElongThrStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_useEllipseAreaRatioCut_allowed()
 *	Description : Execution allowed for useEllipseAreaRatioCut attribute
 */
//--------------------------------------------------------
bool SFinder::is_useEllipseAreaRatioCut_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for useEllipseAreaRatioCut attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::useEllipseAreaRatioCutStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::useEllipseAreaRatioCutStateAllowed_WRITE

	//	Not any excluded states for useEllipseAreaRatioCut attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::useEllipseAreaRatioCutStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::useEllipseAreaRatioCutStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_psEllipseAreaRatioMinThr_allowed()
 *	Description : Execution allowed for psEllipseAreaRatioMinThr attribute
 */
//--------------------------------------------------------
bool SFinder::is_psEllipseAreaRatioMinThr_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for psEllipseAreaRatioMinThr attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::psEllipseAreaRatioMinThrStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::psEllipseAreaRatioMinThrStateAllowed_WRITE

	//	Not any excluded states for psEllipseAreaRatioMinThr attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::psEllipseAreaRatioMinThrStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::psEllipseAreaRatioMinThrStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_psEllipseAreaRatioMaxThr_allowed()
 *	Description : Execution allowed for psEllipseAreaRatioMaxThr attribute
 */
//--------------------------------------------------------
bool SFinder::is_psEllipseAreaRatioMaxThr_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for psEllipseAreaRatioMaxThr attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::psEllipseAreaRatioMaxThrStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::psEllipseAreaRatioMaxThrStateAllowed_WRITE

	//	Not any excluded states for psEllipseAreaRatioMaxThr attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::psEllipseAreaRatioMaxThrStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::psEllipseAreaRatioMaxThrStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_useMaxNPixCut_allowed()
 *	Description : Execution allowed for useMaxNPixCut attribute
 */
//--------------------------------------------------------
bool SFinder::is_useMaxNPixCut_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for useMaxNPixCut attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::useMaxNPixCutStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::useMaxNPixCutStateAllowed_WRITE

	//	Not any excluded states for useMaxNPixCut attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::useMaxNPixCutStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::useMaxNPixCutStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_psMaxNPix_allowed()
 *	Description : Execution allowed for psMaxNPix attribute
 */
//--------------------------------------------------------
bool SFinder::is_psMaxNPix_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for psMaxNPix attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::psMaxNPixStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::psMaxNPixStateAllowed_WRITE

	//	Not any excluded states for psMaxNPix attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::psMaxNPixStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::psMaxNPixStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_useBoundingBoxCut_allowed()
 *	Description : Execution allowed for useBoundingBoxCut attribute
 */
//--------------------------------------------------------
bool SFinder::is_useBoundingBoxCut_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for useBoundingBoxCut attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::useBoundingBoxCutStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::useBoundingBoxCutStateAllowed_WRITE

	//	Not any excluded states for useBoundingBoxCut attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::useBoundingBoxCutStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::useBoundingBoxCutStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_minBoundingBoxThr_allowed()
 *	Description : Execution allowed for minBoundingBoxThr attribute
 */
//--------------------------------------------------------
bool SFinder::is_minBoundingBoxThr_allowed(TANGO_UNUSED(Tango::AttReqType type))
{
	//	Not any excluded states for minBoundingBoxThr attribute in Write access.
	/*----- PROTECTED REGION ID(SFinder::minBoundingBoxThrStateAllowed_WRITE) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::minBoundingBoxThrStateAllowed_WRITE

	//	Not any excluded states for minBoundingBoxThr attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::minBoundingBoxThrStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::minBoundingBoxThrStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_sourceData_allowed()
 *	Description : Execution allowed for sourceData attribute
 */
//--------------------------------------------------------
bool SFinder::is_sourceData_allowed(TANGO_UNUSED(Tango::AttReqType type))
{

	//	Not any excluded states for sourceData attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::sourceDataStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::sourceDataStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_runProgress_allowed()
 *	Description : Execution allowed for runProgress attribute
 */
//--------------------------------------------------------
bool SFinder::is_runProgress_allowed(TANGO_UNUSED(Tango::AttReqType type))
{

	//	Not any excluded states for runProgress attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::runProgressStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::runProgressStateAllowed_READ
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_encodedSourceData_allowed()
 *	Description : Execution allowed for encodedSourceData attribute
 */
//--------------------------------------------------------
bool SFinder::is_encodedSourceData_allowed(TANGO_UNUSED(Tango::AttReqType type))
{

	//	Not any excluded states for encodedSourceData attribute in read access.
	/*----- PROTECTED REGION ID(SFinder::encodedSourceDataStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::encodedSourceDataStateAllowed_READ
	return true;
}

//=================================================
//		pipe Allowed Methods
//=================================================
//--------------------------------------------------------
/**
 *	Method      : SFinder::is_compactSourcesPipe_allowed()
 *	Description : Execution allowed for compactSourcesPipe pipe
 */
//--------------------------------------------------------
bool SFinder::is_compactSourcesPipe_allowed(TANGO_UNUSED(Tango::PipeReqType type))
{
	//	Not any excluded states for compactSourcesPipe pipe in read access.
	/*----- PROTECTED REGION ID(SFinder::compactSourcesPipeStateAllowed_READ) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::compactSourcesPipeStateAllowed_READ
	return true;
}

//=================================================
//		Commands Allowed Methods
//=================================================

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_ExtractSources_allowed()
 *	Description : Execution allowed for ExtractSources attribute
 */
//--------------------------------------------------------
bool SFinder::is_ExtractSources_allowed(TANGO_UNUSED(const CORBA::Any &any))
{
	//	Compare device state with not allowed states.
	if (get_state()==Tango::INIT ||
		get_state()==Tango::FAULT ||
		get_state()==Tango::RUNNING)
	{
	/*----- PROTECTED REGION ID(SFinder::ExtractSourcesStateAllowed) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::ExtractSourcesStateAllowed
		return false;
	}
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_Configure_allowed()
 *	Description : Execution allowed for Configure attribute
 */
//--------------------------------------------------------
bool SFinder::is_Configure_allowed(TANGO_UNUSED(const CORBA::Any &any))
{
	//	Not any excluded states for Configure command.
	/*----- PROTECTED REGION ID(SFinder::ConfigureStateAllowed) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::ConfigureStateAllowed
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_RegisterMe_allowed()
 *	Description : Execution allowed for RegisterMe attribute
 */
//--------------------------------------------------------
bool SFinder::is_RegisterMe_allowed(TANGO_UNUSED(const CORBA::Any &any))
{
	//	Not any excluded states for RegisterMe command.
	/*----- PROTECTED REGION ID(SFinder::RegisterMeStateAllowed) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::RegisterMeStateAllowed
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_Free_allowed()
 *	Description : Execution allowed for Free attribute
 */
//--------------------------------------------------------
bool SFinder::is_Free_allowed(TANGO_UNUSED(const CORBA::Any &any))
{
	//	Not any excluded states for Free command.
	/*----- PROTECTED REGION ID(SFinder::FreeStateAllowed) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::FreeStateAllowed
	return true;
}

//--------------------------------------------------------
/**
 *	Method      : SFinder::is_Reserve_allowed()
 *	Description : Execution allowed for Reserve attribute
 */
//--------------------------------------------------------
bool SFinder::is_Reserve_allowed(TANGO_UNUSED(const CORBA::Any &any))
{
	//	Not any excluded states for Reserve command.
	/*----- PROTECTED REGION ID(SFinder::ReserveStateAllowed) ENABLED START -----*/
	
	/*----- PROTECTED REGION END -----*/	//	SFinder::ReserveStateAllowed
	return true;
}


/*----- PROTECTED REGION ID(SFinder::SFinderStateAllowed.AdditionalMethods) ENABLED START -----*/

//	Additional Methods

/*----- PROTECTED REGION END -----*/	//	SFinder::SFinderStateAllowed.AdditionalMethods

}	//	End of namespace
