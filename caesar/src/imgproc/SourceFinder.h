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
* @file SourceFinder.h
* @class SourceFinder
* @brief Source finder class
*
* Class to perform source finding 
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef SourceFinder_h
#define SourceFinder_h 1

#include <TObject.h>
#include <TMatrixD.h>
#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>

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

class Img;
class Source;
class BkgData;

enum SegmAlgo {
	eWT= 1,
	eHClust= 2,
	eChanVese= 3,
};

class SourceFinder : public TObject {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SourceFinder();
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~SourceFinder();

		enum FileType{
			eROOT= 0,
			eFITS= 1
		};

	public:
		/**
		* \brief Run source finder
		*/
		int Run();

		/**
		* \brief Set options from ConfigParser singleton
		*/
		int Configure();

		//static int FindCompactSource(Img* img,std::vector<Source*>& sources,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool findNestedSources,bool mergeBelowSeed,double peakThreshold);

		/**
		* \brief Read image
		*/
		int ReadImage();

		/**
		* \brief Compute Stats and Bkg info
		*/
		BkgData* ComputeStatsAndBkg(Img*);

	public:
		/**
		* \brief Read only a tile from image
		*/
		int SetTileRead(double xmin,double xmax,double ymin,double ymax){
			if(xmin>=xmax || ymin>=ymax){
				m_ReadTile= false;
				return -1;
			}
			m_ReadTile= true;
			m_TileMinX= xmin;
			m_TileMaxX= xmax;
			m_TileMinY= ymin;
			m_TileMaxY= ymax;
			return 0;
		}	

	private:
		void InitOptions();
		int Init();
		int Save();	
		int FindSources(std::vector<Source*>& sources,Img* inputImg,double seedThr,double mergeThr);
		int FindCompactSources();
		int FindResidualMap();
		int FindExtendedSources();
		int FindExtendedSources_HClust(Img*);
		int FindExtendedSources_ChanVese(Img*);
		int FindExtendedSources_WT(Img*);
	

		int SelectSources(std::vector<Source*>& sources);
		bool IsGoodSource(Source* aSource);
		bool IsPointLikeSource(Source* aSource);
		int DrawSources(Img* image,std::vector<Source*>& sources);

	public:
		
		//Input data
		std::string m_InputFileName;
		std::string m_InputImgName;
		std::string m_InputFileExtension;
		int m_InputFileType;
		Img* m_InputImg;

		//Output data
		TApplication* m_Application;
		bool m_IsInteractiveRun;
		std::string m_OutputFileName;
		TFile* m_OutputFile;
		bool m_SaveToFile;
		bool m_SaveConfig;
		bool m_SaveDS9Region;
		FILE* m_DS9CatalogFilePtr;	
		std::string m_DS9CatalogFileName;
		int m_DS9RegionFormat;		
		TTree* m_SourceTree;
		bool m_SaveSources;
		bool m_SaveResidualMap;

		//Source
		Source* m_Source;
		std::vector<Source*> m_CompactSources;
		std::vector<Source*> m_ExtendedSources;
		std::vector<Source*> m_SourceCollection;

		//Read options
		bool m_ReadTile;
		double m_TileMinX;
		double m_TileMaxX;
		double m_TileMinY;
		double m_TileMaxY;

		//Bkg computation
		BkgData* m_BkgData;
		Img* m_SignificanceMap;
		bool m_UseLocalBkg;	
		bool m_Use2ndPassInLocalBkg;
		bool m_SkipOutliersInLocalBkg;
		int m_LocalBkgMethod;
		int m_BkgEstimator;
		bool m_UseBeamInfoInBkg;
		double m_BoxSizeX;
		double m_BoxSizeY;
		double m_GridSizeX;
		double m_GridSizeY;

		//Residual map
		Img* m_ResidualImg;
		BkgData* m_ResidualBkgData;
		bool m_DilateNestedSources;
		int m_DilateKernelSize;
		int m_DilatedSourceType;
		int m_DilateSourceModel;
		bool m_DilateRandomize;
		bool m_UseResidualInExtendedSearch;

		//Smoothing
		bool m_UsePreSmoothing;	
		int m_SmoothFilter;			
		int m_GausFilterKernSize;
		double m_GausFilterSigma;
		double m_GuidedFilterRadius;
		double m_GuidedFilterColorEps;


		//Compact source search
		bool m_SearchCompactSources;
		int m_NMinPix;
		double m_SeedBrightThr;
		double m_SeedThr;
		double m_MergeThr;
		bool m_MergeBelowSeed;
		bool m_SearchNegativeExcess;

		//Nested source search
		bool m_SearchNestedSources;
		double m_NestedBlobThrFactor;
		
		//Source selection
		bool m_ApplySourceSelection;
		double m_SourceMinBoundingBox;
		double m_psCircRatioThr;
		double m_psElongThr;
		double m_psEllipseAreaRatioMinThr;
		double m_psEllipseAreaRatioMaxThr;
		double m_psMaxNPix;

		//Saliency computation
		Img* m_SaliencyImg;
		double m_SaliencyThrFactor;
		double m_SaliencyBkgThrFactor;
		double m_SaliencyImgThrFactor;
		int m_SaliencyResoMin;
		int m_SaliencyResoMax;
		int m_SaliencyResoStep;
		bool m_SaliencyUseRobustPars;
		bool m_SaliencyUseBkgMap;
		bool m_SaliencyUseNoiseMap;
		bool m_SaliencyUseCurvInDiss;
		double m_SaliencyNNFactor;
		double m_SaliencySpatialRegFactor;
		double m_SaliencyMultiResoCombThrFactor;
		double m_SaliencyDissExpFalloffPar;
		double m_SaliencySpatialDistRegPar;

		//Extended sources
		bool m_SearchExtendedSources;
		int m_ExtendedSearchMethod;
		int m_wtScaleExtended;

		//Superpixel options
		int m_spSize;
		double m_spBeta;
		int m_spMinArea;
		bool m_spUseLogContrast;

		//CHan-Vese options
		double m_cvTimeStepPar;
		double m_cvWindowSizePar;
		double m_cvLambda1Par;
		double m_cvLambda2Par;
		double m_cvMuPar;
		double m_cvNuPar;
		double m_cvPPar;
		


	ClassDef(SourceFinder,1)

};//close SourceFinder


#ifdef __MAKECINT__
#pragma link C++ class SourceFinder+;
#endif

}//close namespace

#endif

