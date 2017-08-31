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
* @file SFinder.h
* @class SFinder
* @brief Source finder class
*
* Class to perform source finding 
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _SFINDER_h
#define _SFINDER_h 1

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

class Image;
class Source;
class ImgBkgData;



class SFinder : public TObject {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SFinder();
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~SFinder();


	public:
		/**
		* \brief Run source finder
		*/
		int Run();

		/**
		* \brief Set options from ConfigParser singleton
		*/
		int Configure();

		/**
		* \brief Read image
		*/
		int ReadImage();

		/**
		* \brief Compute Stats and Bkg info
		*/
		ImgBkgData* ComputeStatsAndBkg(Image*);

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

		/**
		* \brief Initialize class options from ConfigParser
		*/
		void InitOptions();
		/**
		* \brief Initialize class variables
		*/
		int Init();
		/**
		* \brief Save data to output file
		*/
		int Save();	

		/**
		* \brief Save source region file
		*/
		int SaveDS9RegionFile();

		/**
		* \brief Find sources from input image
		*/
		int FindSources(std::vector<Source*>& sources,Image* inputImg,double seedThr,double mergeThr);
		/**
		* \brief Find compact sources
		*/
		int FindCompactSources();
		/**
		* \brief Compute residual map
		*/
		int FindResidualMap();

		/**
		* \brief Find extended sources
		*/
		int FindExtendedSources(Image*);
		/**
		* \brief Find extended sources with hierarchical clustering method
		*/
		int FindExtendedSources_HClust(Image*);
		/**
		* \brief Find extended sources with Chan-Vese method
		*/
		int FindExtendedSources_ChanVese(Image*);
		/**
		* \brief Find extended sources with Wavelet Transform method
		*/
		int FindExtendedSources_WT(Image*);
		/**
		* \brief Find extended sources with Saliency Map thresholding method
		*/
		int FindExtendedSources_SalThr(Image*);

		
		/**
		* \brief Compute edge image
		*/
		Image* ComputeEdgeImage(Image* inputImg,int model);

		/**
		* \brief Compute laplacian image
		*/
		Image* ComputeLaplacianImage(Image* inputImg);

		/**
		* \brief Select sources according to quality cuts given in configuration
		*/
		int SelectSources(std::vector<Source*>& sources);
		/**
		* \brief Tag a source as good or bad
		*/
		bool IsGoodSource(Source* aSource);
		/**
		* \brief Tag a source as point-like or not
		*/
		bool IsPointLikeSource(Source* aSource);
		/**
		* \brief Draw image along with detected sources
		*/
		int DrawSources(Image* image,std::vector<Source*>& sources);

		/**
		* \brief Deblend sources
		*/
		int DeblendSources(std::vector<Source*>& sources);

		/**
		* \brief Print performance stats
		*/
		void PrintPerformanceStats();

	public:
		
		//Input data
		std::string m_InputFileName;
		std::string m_InputImgName;
		std::string m_InputFileExtension;
		int m_InputFileType;
		Image* m_InputImg;

		//Output data
		TApplication* m_Application;
		bool m_IsInteractiveRun;
		std::string m_OutputFileName;
		TFile* m_OutputFile;
		bool m_saveToFile;
		bool m_saveConfig;
		bool m_saveDS9Region;
		//FILE* m_DS9CatalogFilePtr;	
		std::string m_DS9CatalogFileName;
		int m_DS9RegionFormat;		
		TTree* m_SourceTree;
		bool m_saveSources;
		bool m_saveResidualMap;
		bool m_saveInputMap;
		bool m_saveSignificanceMap;
		bool m_saveBkgMap;
		bool m_saveNoiseMap;
		bool m_saveSaliencyMap;
		bool m_saveEdgenessMap;
		bool m_saveCurvatureMap;
		bool m_saveSegmentedMap;

		//Performance stats data
		TTree* m_PerfTree;
		double totTime;
		double initTime;
		double readImageTime;
		double compactSourceTime;
		double sourceSelectionTime;
		double imgResidualTime;
		double extendedSourceTime;
		double sourceDeblendTime;
		double saveTime;

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
		ImgBkgData* m_BkgData;
		Image* m_SignificanceMap;
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
		Image* m_ResidualImg;
		ImgBkgData* m_ResidualBkgData;
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

		//Source deblending
		bool m_deblendSources;
		double m_deblendCurvThr;
		double m_deblendComponentMinNPix;
		

		//Saliency computation
		Image* m_SaliencyImg;
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

		//Hierachical clustering data
		Image* m_LaplImg;
		Image* m_EdgeImg;
		Image* m_SegmImg;
		int m_spMergingNSegmentsToStop;
		double m_spMergingRatio;
		double m_spMergingRegPar;
		double m_spMergingMaxDissRatio;
		double m_spMergingMaxDissRatio2ndNeighbours;
		double m_spMergingDissThreshold;
		int m_spMergingEdgeModel;		
		bool m_spMergingUse2ndNeighbours;
		bool m_spMergingIncludeSpatialPars;
		bool m_spMergingAddCurvDist;
		bool m_spMergingUseRobustPars;
		
	ClassDef(SFinder,1)

};//close SourceFinder


#ifdef __MAKECINT__
#pragma link C++ class SFinder+;
#endif

}//close namespace

#endif

