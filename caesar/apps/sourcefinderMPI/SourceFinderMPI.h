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
* @file SourceFinderMPI.h
* @class SourceFinderMPI
* @brief Source finder MPI class
*
* Class to perform source finding in MPI
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef SourceFinderMPI_h
#define SourceFinderMPI_h 1


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

#include <mpi.h>


namespace Caesar {

class Image;
class Source;
class ImgBkgData;
class TaskData;


class SourceFinderMPI : public TObject {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SourceFinderMPI();
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~SourceFinderMPI();

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

		/**
		* \brief Read image
		*/
		Image* ReadImage(long int ix_min,long int ix_max,long int iy_min,long int iy_max);

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
		void InitOptions();
		int Init();
		int PrepareWorkerTasks();
		int Save();	
		int FindSources(std::vector<Source*>& sources,Image* inputImg,double seedThr,double mergeThr);
		int FindCompactSources(TaskData*,Image*);
		int FindResidualMap();
		int FindExtendedSources();
		int FindExtendedSources_HClust(Image*);
		int FindExtendedSources_ChanVese(Image*);
		int FindExtendedSources_WT(Image*);
	

		int SelectSources(std::vector<Source*>& sources);
		//int SelectSources(TaskData*);

		bool IsGoodSource(Source* aSource);
		bool IsPointLikeSource(Source* aSource);
		int DrawSources(Image* image,std::vector<Source*>& sources);

		int UpdateTaskDataFromWorkers();
		int FindSourcesAtEdge();
		int MergeSourcesAtEdge();

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
		bool m_SaveToFile;
		bool m_SaveConfig;
		bool m_SaveDS9Region;
		FILE* m_DS9CatalogFilePtr;	
		std::string m_DS9CatalogFileName;
		int m_DS9RegionFormat;		
		TTree* m_SourceTree;
		bool m_SaveSources;
		bool m_SaveResidualMap;
	
		TTree* m_TaskInfoTree;
		double m_xmin;
		double m_xmax;
		double m_ymin;
		double m_ymax;

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

		//Read distributed options
		long int m_TileSizeX;
		long int m_TileSizeY;
		bool m_UseTileOverlap;
		double m_TileStepSizeX;
		double m_TileStepSizeY;

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
		
		//Task data
		std::vector< std::vector<TaskData*> > m_taskDataPerWorkers;
	
		//MPI vars
		int nproc;
		int myid;
		MPI_Group m_WorldGroup;
		MPI_Group m_WorkerGroup;
		MPI_Comm m_WorkerComm;
		int worker_ranks;
		int nworkers;

	//ClassDef(SourceFinderMPI,1)

};//close SourceFinder

/*
#ifdef __MAKECINT__
#pragma link C++ class SourceFinderMPI+;
#endif
*/

}//close namespace

#endif

