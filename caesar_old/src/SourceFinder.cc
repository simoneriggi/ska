/**
* @file SourceFinder.cc
* @class SourceFinder
* @brief SourceFinder
*
* SourceFinder class
* @author S. Riggi
* @date 20/01/2015
*/

#include <SourceFinder.h>

#include <ConfigParser.h>
#include <ImgFITSReader.h>
#include <Img.h>
#include <Contour.h>
#include <Source.h>
//#include <VLSlicSegmentation.h>
#include <ChanVeseSegmentation.h>

#include <RInside.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TGraph2D.h>
#include <TF1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TCut.h>
#include <TEventList.h>
#include <TMath.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TColor.h>

#include <Math/WrappedTF1.h>
#include <Math/GSLIntegrator.h>
#include <Math/GSLMinimizer.h>
#include <Math/Functor.h>
#include <Math/WrappedFunction.h>
#include <Math/WrappedParamFunction.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>
#include <Math/SpecFunc.h>
#include <Math/DistFunc.h>


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

//ClassImp(SourceFinder)


SourceFinder::SourceFinder(){
	
	fR= 0;
	fInputImg= 0;
	fSource= 0;
	fDS9CatalogFilePtr= 0;
	fOutputFile= 0;
	fOutputTree= 0;
	fApplication= 0;
	fSourceCollection.clear();
	myStyle= 0;
	fSaveToFile= true;
	//fSLICSegmentation= 0;
	fSLICSegmenter= 0;
	fCVSegmentation= 0;
	fSegmentedImage= 0;
	fFinalSegmentedImage= 0;
	fFinalSignalSegmentedImage= 0;
	fSegmentationContourGraph= 0;
	fIsInteractive= true;
	fConfigInfo= 0;
	fSPMergingInfo= 0;

	fSaliencyMap= 0;
	fSumSaliencyMap= 0;
	fInitSaliencyMap= 0;
	
}//close costructor


SourceFinder::~SourceFinder(){
	
	/*
	cout<<"SourceFinder::~SourceFinder(): INFO: Clear SPSegmentation..."<<endl;
	if(fSLICSegmentation) {
		delete fSLICSegmentation;
		fSLICSegmentation= 0;
	}
	cout<<"SourceFinder::~SourceFinder(): INFO: Clear CVSegmentation..."<<endl;
	if(fCVSegmentation) {
		delete fCVSegmentation;
		fCVSegmentation= 0;
	}
	*/
}//close destructor


int SourceFinder::SetConfig(){

	//Get input file name
	fInputFileName= ConfigParser::fInputFileName;
	if (fInputFileName == "") {
		cerr<<"SourceFinder::SetConfig(): ERROR: Empty input file name!"<<endl;
  	return -1;
  }
	fInputFileExtension= fInputFileName.substr(fInputFileName.find_last_of(".") + 1);
	if(fInputFileExtension!= "fits" && fInputFileExtension!="root" &&  fInputFileExtension!="png" &&  fInputFileExtension!="jpg" &&  fInputFileExtension!="bmp") {
		cerr<<"SourceFinder::SetConfig(): ERROR: Invalid file extension ("<<fInputFileExtension<<")...nothing to be done!"<<endl;
		return -1;
	}

	//Get output files
	fOutputFileName= ConfigParser::fOutputFileName;
	fDS9CatalogFileName= ConfigParser::fDS9CatalogFileName;
	fSaveToFile= ConfigParser::fSaveToFile;
	fSaveImageToFile= ConfigParser::fSaveImageToFile;
	fSaveImageType= ConfigParser::fSaveImageType;
	fDS9RegionFormat= ConfigParser::fDS9RegionFormat;
	fDrawSources= ConfigParser::fDrawSources;
	fIsInteractive= ConfigParser::fIsInteractive;

	//Get bkg options
	fUseLocalBackground= ConfigParser::fUseLocalBkg;
	fLocalBkgModel= ConfigParser::fLocalBkgMethod;
	fBkgEstimator= ConfigParser::fBkgEstimator;
	fBkgGridSize= ConfigParser::fGridSize;
	fBkgBoxSize= ConfigParser::fBoxSize;
	fBkgGridSizeX= ConfigParser::fGridSizeX;
	fBkgGridSizeY= ConfigParser::fGridSizeY;
	fBkgBoxSizeX= ConfigParser::fBoxSizeX;
	fBkgBoxSizeY= ConfigParser::fBoxSizeY;

	//Get source finding options
	fDeblendSources= ConfigParser::fDeblendSources;	
	fCurvatureThreshold= ConfigParser::fCurvatureThreshold;
	fPeakThreshold= ConfigParser::fPeakThreshold;
	fSourceComponentMinNPix= ConfigParser::fSourceComponentMinNPix;
	fSearchNestedSources= ConfigParser::fSearchNestedSources;
	fNPixMin= ConfigParser::fNPixMin;
	fSeedBrightThreshold= ConfigParser::fSeedBrightThreshold;
	fSeedThreshold= ConfigParser::fSeedThreshold;
	fMergeThreshold= ConfigParser::fMergeThreshold;
	fSearchBrightSources= ConfigParser::fSearchBrightSources;
	fSearchExtendedSources= ConfigParser::fSearchExtendedSources;
	fSearchFaintSources= ConfigParser::fSearchFaintSources;
	fExtendedSearchMethod= ConfigParser::fExtendedSearchMethod;
	fSearchNegativeExcess= ConfigParser::fSearchNegativeExcess;
	fUseCurvatureMixture= ConfigParser::fUseCurvatureMixture;
	fCurvatureWeight= ConfigParser::fCurvatureWeight;
	fWTScaleForFaintSourceSearch= ConfigParser::fWTScaleForFaintSourceSearch;
	fWTScaleForExtendedSourceSearch= ConfigParser::fWTScaleForExtendedSourceSearch;
	fUseResidualImageInExtendedSearch= ConfigParser::fUseResidualImageInExtendedSearch;

	//Source selection	
	fApplySourceSelection= ConfigParser::fApplySourceSelection;
	fMinBoundingBox= ConfigParser::fMinBoundingBox;
	fTagPointSources= ConfigParser::fTagPointSources;
	fPointSourceCircRatioThr= ConfigParser::fPointSourceCircRatioThr;
	fPointSourceElongThr= ConfigParser::fPointSourceElongThr;
	fPointSourceMinEllipseAreaRatio= ConfigParser::fPointSourceMinEllipseAreaRatio;
	fPointSourceMaxEllipseAreaRatio= ConfigParser::fPointSourceMaxEllipseAreaRatio;
	fPointSourceMaxNPix= ConfigParser::fPointSourceMaxNPix;

	//Source residual options
	fSourceDilateKernelSize= ConfigParser::fSourceDilateKernelSize;
	if(fSourceDilateKernelSize%2==0){
		cerr<<"SourceFinder::SetConfig(): ERROR: Given dilate kernel size is not odd!"<<endl;
		return -1;
	}
	fDilateNestedSources= ConfigParser::fDilateNestedSources;
	fDilatedSourceType= ConfigParser::fDilatedSourceType;
	fDilateSourceModel= ConfigParser::fDilateSourceModel;
	fRandomizeInDilate= ConfigParser::fRandomizeInDilate;
	fRandSigmaInDilate= ConfigParser::fRandSigmaInDilate;

	//SP Segmentation options		
	fSPSize= ConfigParser::fSPSize;
	fSPRegularization= ConfigParser::fSPRegularization;
	fSPMinArea= ConfigParser::fSPMinArea;
	fSPMergingDistEps= ConfigParser::fSPMergingDistEps;
	fSPMergingAlgo= ConfigParser::fSPMergingAlgo;
	fSPMergingRatio= ConfigParser::fSPMergingRatio;
	fSPMergingRegularization= ConfigParser::fSPMergingRegularization;
	fUse2ndNeighborsInSPMerging= ConfigParser::fUse2ndNeighborsInSPMerging;
	fMinMergedSP= ConfigParser::fMinMergedSP;
	fSPMergingDistThreshold= ConfigParser::fSPMergingDistThreshold;
	fSPUseLogContrast= ConfigParser::fSPUseLogContrast;
	fUsePixelRatioCut= ConfigParser::fUsePixelRatioCut;
	fPixelRatioCut= ConfigParser::fPixelRatioCut;
	fTagSignificativeSP= ConfigParser::fTagSignificativeSP;
	fSPTaggingMethod= ConfigParser::fSPTaggingMethod;
	fSignificantSPRatio= ConfigParser::fSignificantSPRatio;
	fSPMergingEdgeModel= ConfigParser::fSPMergingEdgeModel;
	
	fSPMergingUseAdaptingDistThreshold= ConfigParser::fSPMergingUseAdaptingDistThreshold;
	fSPMergingAdaptingDistThresholdScale= ConfigParser::fSPMergingAdaptingDistThresholdScale;
	fSPMergingMaxDissRatio= ConfigParser::fSPMergingMaxDissRatio;
	fSPMergingMaxDissRatio2ndNeighbor= ConfigParser::fSPMergingMaxDissRatio2ndNeighbor;
	fSPMergingAggloMethod= ConfigParser::fSPMergingAggloMethod;
	fSPMergingMinClustSize= ConfigParser::fSPMergingMinClustSize;
	fSPMergingMaxHeightQ= ConfigParser::fSPMergingMaxHeightQ;
	fSPMergingDeepSplitLevel= ConfigParser::fSPMergingDeepSplitLevel;
	fUseCurvatureInSPMerging= ConfigParser::fUseCurvatureInSPMerging;

	//Saliency map
	fSaliencyThresholdFactor= ConfigParser::fSaliencyThresholdFactor;
	fBkgSaliencyThresholdFactor= ConfigParser::fBkgSaliencyThresholdFactor;
	fSaliencyImgThresholdFactor= ConfigParser::fSaliencyImgThresholdFactor;
	fSaliencyMinReso= ConfigParser::fSaliencyMinReso;
	fSaliencyMaxReso= ConfigParser::fSaliencyMaxReso;
	fSaliencyResoStepSize= ConfigParser::fSaliencyResoStepSize;
	fSaliencyUseRobustPars= ConfigParser::fSaliencyUseRobustPars;
	fSaliencyUseBkgMap= ConfigParser::fSaliencyUseBkgMap;
	fSaliencyUseNoiseMap= ConfigParser::fSaliencyUseNoiseMap;
	fSaliencyUseCurvatureMap= ConfigParser::fSaliencyUseCurvatureMap;
	fSaliencyNNFactor= ConfigParser::fSaliencyNNFactor;
	fSaliencyFilterThresholdFactor= ConfigParser::fSaliencyFilterThresholdFactor;
	fSaliencyNormalizationMode= ConfigParser::fSaliencyNormalizationMode;

	//Chan-Vese segmentation options
	fCVTimeStep= ConfigParser::fCVTimeStep;
	fCVWindowSize= ConfigParser::fCVWindowSize;
	fCVLambda1Par= ConfigParser::fCVLambda1Par;
	fCVLambda2Par= ConfigParser::fCVLambda2Par;
	fCVMuPar= ConfigParser::fCVMuPar;
	fCVNuPar= ConfigParser::fCVNuPar;
	fCVPPar= ConfigParser::fCVPPar;
		
	//Smoothing
	fUsePreSmoothing= ConfigParser::fUsePreSmoothing;
	fSmoothingAlgo= ConfigParser::fSmoothingAlgo;
	fSmoothKernelSize= ConfigParser::fSmoothKernelSize;
	fSmoothSigma= ConfigParser::fSmoothSigma;
	fGuidedSmoothRadius= ConfigParser::fGuidedSmoothRadius;
	fGuidedSmoothColorEps= ConfigParser::fGuidedSmoothColorEps;

	return 0;

}//close SetConfig()

int SourceFinder::Init(){

	//## Init RInside ptr
	fR= RInside::instancePtr();
	if(!fR){
		cerr<<"SourceFinder::Init(): ERROR: Cannot retrieve RInside instance (did you create one in your main application?)!"<<endl;
		fR= new RInside;
	}	

	//## Init ROOT app
	if(!fApplication && fIsInteractive){
		fApplication= new TApplication("Application", 0, 0);
	}

	//## Set config options
	if(SetConfig()<0) return -1;

	//## Set graphics style
	SetGraphicsStyle();
	
	//## Create region file
	fDS9CatalogFilePtr= fopen(fDS9CatalogFileName.c_str(),"w");
	fDS9CatalogFilePtr2= fopen("EllipseDS9Region.reg","w");

	//## Create output file
	fOutputFile= new TFile(fOutputFileName.c_str(),"RECREATE");	
	fOutputFile->cd();

	//## Init source tree
	fSource= 0;
	fOutputTree= new TTree("SourceInfo","SourceInfo");
	fOutputTree->Branch("Source",&fSource);

	fSourceCollection.clear();

	//## Init config tree
	if(!fConfigInfo){
		fConfigInfo= new TTree("ConfigInfo","ConfigInfo");
		fConfigInfo->Branch("l",&fSPMinArea,"l/I");
		fConfigInfo->Branch("beta",&fSPRegularization,"beta/D");
		fConfigInfo->Branch("lambda",&fSPMergingRegularization,"lambda/D");
		fConfigInfo->Branch("dissThr",&fSPMergingDistEps,"dissThr/D");
		fConfigInfo->Fill();
	}
	if(!fSPMergingInfo){
		fSPMergingInfo= new TTree("SPMergingInfo","SPMergingInfo");
		fSPMergingInfo->Branch("levelId",&fLevelId,"levelId/I");
		fSPMergingInfo->Branch("nDissEntries",&fNDissEntries,"nDissEntries/I");
		fSPMergingInfo->Branch("DissList",fDissList,"DissList[nDissEntries]/D");
		fSPMergingInfo->Branch("nDissEntriesMerged",&fNDissEntries_merged,"nDissEntriesMerged/I");
		fSPMergingInfo->Branch("DissListMerged",fDissList_merged,"DissListMerged[nDissEntriesMerged]/D");
		fSPMergingInfo->Branch("minDissMerged",&fDissMin_merged,"minDissMerged/D");
		fSPMergingInfo->Branch("maxDissMerged",&fDissMax_merged,"maxDissMerged/D");
		fSPMergingInfo->Branch("minDiss",&fDissMin,"minDiss/D");
		fSPMergingInfo->Branch("maxDiss",&fDissMax,"maxDiss/D");
		fSPMergingInfo->Branch("medianDiss",&fDissMedian,"medianDiss/D");	
		fSPMergingInfo->Branch("medianrmsDiss",&fDissMedianRMS,"medianrmsDiss/D");
		fSPMergingInfo->Branch("medianDiss0",&fDissMedian0,"medianDiss0/D");	
		fSPMergingInfo->Branch("medianrmsDiss0",&fDissMedianRMS0,"medianrmsDiss0/D");
		fSPMergingInfo->Branch("nr",&fNR,"nr/I");
		fSPMergingInfo->Branch("mse",&fMSE,"mse/D");
	}

	//## Init classes
	/*
	if(!fSLICSegmentation) fSLICSegmentation= new VLSlicSegmentation;
	fSLICSegmentation->SetLogContrastMapping(fSPUseLogContrast);
	fSLICSegmentation->SetSPMergingRatio(fSPMergingRatio);
	fSLICSegmentation->SetSPMergingRegularization(fSPMergingRegularization);
	fSLICSegmentation->Use2ndNeighborsInSPMerging(fUse2ndNeighborsInSPMerging);
	fSLICSegmentation->SetMinMergedSP(fMinMergedSP);
	fSLICSegmentation->SetSPMergingDistThreshold(fSPMergingDistThreshold);
	fSLICSegmentation->UseRobustParamsInSPMerging(false);
	fSLICSegmentation->UsePixelRatioCut(fUsePixelRatioCut);
	fSLICSegmentation->SetPixelRatioCut(fPixelRatioCut);
	fSLICSegmentation->TagSignificantSP(fTagSignificativeSP);
	fSLICSegmentation->SetSignificantSPTaggingMethod(fSPTaggingMethod);
	fSLICSegmentation->SetSignificantSPRatioCut(fSignificantSPRatio);
	fSLICSegmentation->SetSaliencyThresholdFactor(fSaliencyThresholdFactor);
	*/

	if(!fSLICSegmenter) fSLICSegmenter= new SLICSegmenter;
	fSLICSegmenter->SetLogContrastMapping(fSPUseLogContrast);
	fSLICSegmenter->SetSPMergingRatio(fSPMergingRatio);
	fSLICSegmenter->SetSPMergingRegularization(fSPMergingRegularization);
	fSLICSegmenter->SetMinMergedSP(fMinMergedSP);
	fSLICSegmenter->SetSPMergingDistThreshold(fSPMergingDistThreshold);
	fSLICSegmenter->UsePixelRatioCut(fUsePixelRatioCut);
	fSLICSegmenter->SetPixelRatioCut(fPixelRatioCut);
	fSLICSegmenter->SetSignificantSPRatioCut(fSignificantSPRatio);
	fSLICSegmenter->UseAdaptiveDistThreshold(fSPMergingUseAdaptingDistThreshold);
	fSLICSegmenter->SetAdaptiveThresholdScale(fSPMergingAdaptingDistThresholdScale);
	fSLICSegmenter->SetMaxDissRatio(fSPMergingMaxDissRatio);
	fSLICSegmenter->SetMaxDissRatioFor2ndNeighbors(fSPMergingMaxDissRatio2ndNeighbor);	
	fSLICSegmenter->SetSPMergingAggloMethod(fSPMergingAggloMethod);
	fSLICSegmenter->SetSPMergingMinClustSize(fSPMergingMinClustSize);
	fSLICSegmenter->SetSPMergingMaxHeightQ(fSPMergingMaxHeightQ);
	fSLICSegmenter->SetSPMergingDeepSplitLevel(fSPMergingDeepSplitLevel);
	fSLICSegmenter->UseCurvatureInSPMerging(fUseCurvatureInSPMerging);
	fSLICSegmenter->SetEdgeModel(fSPMergingEdgeModel);
	fSLICSegmenter->Use2ndNeighborsInSPMerging(fUse2ndNeighborsInSPMerging);

	fSLICSegmenter->SetSaliencyThresholdFactor(fSaliencyThresholdFactor);
	fSLICSegmenter->SetBkgSaliencyThresholdFactor(fBkgSaliencyThresholdFactor);
	fSLICSegmenter->SetSaliencyImgThresholdFactor(fSaliencyImgThresholdFactor);
	fSLICSegmenter->SetSaliencyResoPars(fSaliencyMinReso,fSaliencyMaxReso,fSaliencyResoStepSize);
	fSLICSegmenter->UseRobustParsInSaliency(fSaliencyUseRobustPars);
	fSLICSegmenter->UseBkgMapInSaliency(fSaliencyUseBkgMap);
	fSLICSegmenter->UseNoiseMapInSaliency(fSaliencyUseNoiseMap);
	fSLICSegmenter->UseCurvatureMapInSaliency(fSaliencyUseCurvatureMap);
	fSLICSegmenter->SetNNFactorInSaliency(fSaliencyNNFactor);
	fSLICSegmenter->SetFilterThresholdFactorInSaliency(fSaliencyFilterThresholdFactor);
	fSLICSegmenter->SetSaliencyNormalizationMode(fSaliencyNormalizationMode);
	
	fSLICSegmenter->SetBoxSizeInLocalBkgMap(fBkgBoxSizeX,fBkgBoxSizeX);
	fSLICSegmenter->SetGridSizeInLocalBkgMap(fBkgGridSizeX,fBkgGridSizeY);

	fSLICSegmenter->SetMinNPixInSaliencyThresholding(fNPixMin);

	fSLICSegmenter->SetCVTimeStep(fCVTimeStep);
	fSLICSegmenter->SetCVWindowSize(fCVWindowSize);
	fSLICSegmenter->SetCVLambda1Par(fCVLambda1Par);
	fSLICSegmenter->SetCVLambda2Par(fCVLambda2Par);
	fSLICSegmenter->SetCVMuPar(fCVMuPar);
	fSLICSegmenter->SetCVNuPar(fCVNuPar);
	fSLICSegmenter->SetCVPPar(fCVPPar);

	if(!fCVSegmentation) fCVSegmentation= new ChanVeseSegmentation;

	return 0;

}//close Init()


int SourceFinder::Run(){

	//## Init pipeline
	if(Init()<0){
		cerr<<"SourceFinder::Run(): ERROR: Initialization failed!"<<endl;
		return -1;
	}

	//## Read input image
	if(ReadImage()<0){
		cerr<<"SourceFinder::Run(): ERROR: Reading of input image failed!"<<endl;
		return -1;
	}

	//## Find bright sources
	if(fSearchBrightSources && FindBrightSource()<0){
		cerr<<"SourceFinder::Run(): ERROR: Bright source search failed!"<<endl;
		return -1;
	}

	//## Find faint sources
	if(fSearchFaintSources && FindFaintSource()<0){
		cerr<<"SourceFinder::Run(): ERROR: Faint source search failed!"<<endl;
		return -1;
	}

	//## Find extended sources
	if(fSearchExtendedSources && FindExtendedSource()<0){
		cerr<<"SourceFinder::Run(): ERROR: Extended source search failed!"<<endl;
		return -1;
	}

	//## Deblend sources
	if(fDeblendSources) DeblendSources(fInputImg);

	//## Draw final sources
	if(fDrawSources) DrawSources(fInputImg,false);

	//## Save to file
	if(fSaveToFile) Save();
	
	if(fApplication && fIsInteractive) fApplication->Run();

	return 0;

}//close Run()

int SourceFinder::DeblendSources(Img* img){

	if(!img) return -1;

	//## Loop over image sources and perform deblending stage for non-extended sources
	std::vector<Source*> sources= img->GetSources();
	if(sources.size()<=0) return 0;
	for(unsigned int i=0;i<sources.size();i++){
		int sourceType= sources[i]->fType;
		if(sourceType!=Source::ePointLike) continue;
		int status= sources[i]->Deblend(fCurvatureThreshold,fSourceComponentMinNPix);
		if(status<0) cerr<<"SourceFinder::DeblendSources(): WARN: Deblending source no. "<<i<<" failed!"<<endl;

		//Deblend nested sources
		if(sources[i]->fHasNestedSources){
			for(unsigned int j=0;j<sources[i]->fNestedSourceCollection.size();j++){
				if((sources[i]->fNestedSourceCollection)[j]) 
					(sources[i]->fNestedSourceCollection)[j]->Deblend(fCurvatureThreshold,fSourceComponentMinNPix);
			}//end loop nested sources
		}//close if nested sources
	}//end loop sources

	return 0;

}//close DeblendSources()

bool SourceFinder::IsGoodSource(Source* aSource){
	
	if(!aSource) return false;

	//## Check for pixels 	
	if(aSource->fNPix<=0 || aSource->fPixelCollection.size()<=0) return false;

	//## Check for line-like source
	if(aSource->fContourCollection.size()<=0) {
		cerr<<"SourceFinder::IsGoodSource(): WARN: No contour stored for this source, cannot perform check!"<<endl;
		return true;
	}

	double BoundingBoxMin= ((aSource->fContourCollection)[0])->BoundingBoxMin;
	if(BoundingBoxMin<fMinBoundingBox) {
		cerr<<"SourceFinder::IsGoodSource(): INFO: BoundingBox cut not passed (BoundingBoxMin="<<BoundingBoxMin<<"<2)"<<endl;
		return false;
	}

	//## Add other check here ...
	//...
	//...

	return true;

}//close SourceFinder::isGoodSource()


bool SourceFinder::IsCompactPointLike(Source* aSource){

	if(!aSource) return false;

	if(!aSource->fHasParameters) {
		cerr<<"SourceFinder::IsCompactSource(): WARN: No parameters are available for this source (did you compute them?)...test cannot be performed!"<<endl;
		return true;
	}

	std::string sourceName= aSource->fName;
	int sourceId= aSource->fId;

	//Loop over contours and check if all of them have circular features
	bool isPointLike= true;
	for(int i=0;i<(aSource->fContourCollection).size();i++){
		Contour* thisContour= aSource->fContourCollection[i];

		/*
		//Test circularity ratio: 1= circle
		if(thisContour->CircularityRatio<fPointSourceCircRatioThr) {
			cout<<"SourceFinder::IsCompactSource(): INFO: Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass CircularityRatio cut (CR="<<thisContour->CircularityRatio<<"<"<<fPointSourceCircRatioThr<<")"<<endl;
			isPointLike= false;
			break;
		}
		*/

		//Test elongation (how symmetrical is the shape): 0=circle,square
		if(thisContour->Elongation>fPointSourceElongThr) {
			cout<<"SourceFinder::IsCompactSource(): INFO: Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass Elongation cut (ELONG="<<thisContour->CircularityRatio<<">"<<fPointSourceElongThr<<")"<<endl;
			isPointLike= false;
			break;	
		}

		//Test ellipse fit
		if(thisContour->EllipseAreaRatio<fPointSourceMinEllipseAreaRatio || thisContour->EllipseAreaRatio>fPointSourceMaxEllipseAreaRatio) {
			cout<<"SourceFinder::IsCompactSource(): INFO: Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass EllipseAreaRatio cut (EAR="<<thisContour->EllipseAreaRatio<<" outside range ["<<fPointSourceMinEllipseAreaRatio<<","<<fPointSourceMaxEllipseAreaRatio<<"])"<<endl;
			isPointLike= false;
			break;	
		}

	}//end contour loop
	
	//Check number of pixels
	if(aSource->fNPix>fPointSourceMaxNPix){
		cout<<"SourceFinder::IsCompactSource(): INFO: Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass nMaxPix cut (NPix="<<aSource->fNPix<<">"<<fPointSourceMaxNPix<<")"<<endl;
		isPointLike= false;
	}

	if(!isPointLike) return false;

	return true;

}//close SourceFinder::IsCompactSource()


int SourceFinder::FindBrightSource(){
	
	//## Find bright sources in input image (or input+curvature mixture image if selected)
	int status= FindSource(fInputImg,fSeedBrightThreshold,fMergeThreshold);
	if(status<0){
		cerr<<"SourceFinder::FindBrightSource(): ERROR: Bright source search failed!"<<endl;
		return -1;
	}

	//## Retrieve found sources 
	std::vector<Source*> sources= fInputImg->GetSources();
	int nSources= (int)sources.size();
	cout<<"SourceFinder::FindBrightSource(): INFO: "<<nSources<<" bright sources detected in input image..."<<endl;
		
	if(nSources<=0) return 0;

	//## Apply source selection?
	int nSelSources= nSources;

	if(fApplySourceSelection){
		nSelSources= 0;

		for(int i=0;i<nSources;i++){	
			std::string sourceName= sources[i]->fName;
			int sourceId= sources[i]->fId;

			//Is bad source (i.e. line-like blob, etc...)?
			if(!IsGoodSource(sources[i])) {
				cout<<"SourceFinder::FindBrightSource(): INFO: Source no. "<<i<<" (name="<<sourceName<<",id="<<sourceId<<") tagged as bad source, skipped!"<<endl;
				sources[i]->fIsGoodSource= false;
				continue;
			}
			
			//Is point-like source?
			if( IsCompactPointLike(sources[i]) ){
				cout<<"SourceFinder::FindBrightSource(): INFO: Source no. "<<i<<" (name="<<sourceName<<",id="<<sourceId<<") tagged as a point-like source ..."<<endl;
				sources[i]->SetType(Source::ePointLike);
			}
			//else{//extended source, search for nested?
				if(fSearchNestedSources) {
					cout<<"SourceFinder::FindBrightSource(): INFO: Finding sources nested to source no. "<<i<<" (name="<<sourceName<<",id="<<sourceId<<")"<<endl;
					int nestedSourceSearchStatus= FindNestedSource(sources[i]);
					if(nestedSourceSearchStatus<0) cerr<<"SourceFinder::FindBrightSource(): WARN: Nested source search failed!"<<endl;
				}	
			//}

			//Add source to the list	
			fSourceCollection.push_back(sources[i]);
			nSelSources++;
		}//end loop sources

	}//close if source selection
	else{
		//Add all sources to the list	
		fSourceCollection.insert(fSourceCollection.end(),sources.begin(),sources.end());
	}

	cout<<"SourceFinder::FindBrightSource(): INFO: Added "<<nSelSources<<" bright sources to the list..."<<endl;	

	
	// Draw bright sources
	//cout<<"SourceFinder::FindBrightSource(): INFO: Drawing "<<sources.size()<<" sources in input image..."<<endl;
	//DrawSources(fInputImg);
	
	return 0;

}//close FindBrightSource()


int SourceFinder::FindFaintSource(){

	//## Find bright source residual map: all sources are dilated, no nested
	Img* sourceResidualImg= fInputImg->GetSourceResidual(fUseLocalBackground,fSourceDilateKernelSize,false,-1,fDilateSourceModel,fRandomizeInDilate,fRandSigmaInDilate);
	if(!sourceResidualImg){
		cerr<<"SourceFinder::FindFaintSource(): WARN: Failed to get bright source residual map...setting it to input image!"<<endl;
		sourceResidualImg= (Img*)fInputImg->Clone("img-SourceResidual");	
	}

	//## Find faint source in the W1 scale of the residual image
	cout<<"SourceFinder::FindFaintSource(): INFO: Find faint sources in the "<<fWTScaleForFaintSourceSearch<<"-WT scale of bright source residual image ..."<<endl;
	std::vector<Img*> wt_faint= sourceResidualImg->GetWaveletDecomposition(fWTScaleForFaintSourceSearch);
		
	int status= FindSource(wt_faint[fWTScaleForFaintSourceSearch],fSeedThreshold,fMergeThreshold);
	if(status<0){
		cerr<<"SourceFinder::FindFaintSource(): ERROR: Faint source finding failed!"<<endl;
		return -1;
	}
	std::vector<Source*> FaintSourceCollection= wt_faint[fWTScaleForFaintSourceSearch]->GetSources();
	int nSources= (int)FaintSourceCollection.size();				

	cout<<"SourceFinder::FindFaintSource(): INFO: "<<nSources<<" faint sources found!"<<endl;
	if(nSources>0) {
		//## Apply source selection?
		int nSelSources= nSources;

		if(fApplySourceSelection){
			nSelSources= 0;

			for(int i=0;i<nSources;i++){
				//Is bad source (i.e. line-like blob, etc...)?
				if(!IsGoodSource(FaintSourceCollection[i])) {
					cout<<"SourceFinder::FindFaintSource(): INFO: Source no. "<<i<<" tagged as bad source, skipped!"<<endl;
					FaintSourceCollection[i]->fIsGoodSource= false;
					continue;
				}
			
				//Is point-like source?
				if( IsCompactPointLike(FaintSourceCollection[i]) ){
					cout<<"SourceFinder::FindFaintSource(): INFO: Candidate bright source tagged as a point-like source ..."<<endl;
					FaintSourceCollection[i]->SetType(Source::ePointLike);
				}
				//else{//extended source, search for nested?
					if(fSearchNestedSources) {
						cout<<"SourceFinder::FindFaintSource(): INFO: Finding nested source in bright source: "<<FaintSourceCollection[i]->fName<<endl;
						int nestedSourceSearchStatus= FindNestedSource(FaintSourceCollection[i]);
						if(nestedSourceSearchStatus<0) cerr<<"SourceFinder::FindFaintSource(): WARN: Nested source search failed!"<<endl;
					}	
				//}

				//Add source to the list	
				fSourceCollection.push_back(FaintSourceCollection[i]);
				nSelSources++;
			}//end loop sources

		}//close if source selection
		else{
			//Add faint sources to the list
			fSourceCollection.insert(fSourceCollection.end(),FaintSourceCollection.begin(),FaintSourceCollection.end());
		}		
		cout<<"SourceFinder::FindFaintSource(): INFO: Added "<<nSelSources<<" faint sources to the list..."<<endl;

		//Draw faint sources
		//cout<<"SourceFinder::FindFaintSource(): INFO: Drawing "<<FaintSourceCollection.size()<<" faint sources in input image..."<<endl;
		//DrawSources(wt_faint[fWTScaleForFaintSourceSearch]);
	}
	
	//## Clear-up	
	cout<<"SourceFinder::FindFaintSource(): INFO: Clearing up..."<<endl;
	if(sourceResidualImg){
		delete sourceResidualImg;
		sourceResidualImg= 0;
	}
	
	return 0;

}//close FindFaintSource()


int SourceFinder::FindExtendedSource(){

	Img* inputImg= (Img*)fInputImg->Clone("inputImg");
	inputImg->SetNameTitle("inputImg","inputImg");
	if(fInputImg->HasMetaData()) inputImg->SetMetaData(fInputImg->GetMetaData());
	if(fInputImg->HasBkgData()) {
		inputImg->SetBkgData(fInputImg->GetBkgData());
		inputImg->SetInterpolatedBkgLevelMap(fInputImg->GetInterpolatedBkgLevelMap());
		inputImg->SetInterpolatedBkgRMSMap(fInputImg->GetInterpolatedBkgRMSMap());
	}
	inputImg->ResetSources();

	Img* imgToSegment= inputImg;

	//## Get source residual with POINTLIKE Sources removed?
	cout<<"SourceFinder::FindExtendedSource(): INFO: nSources="<<fSourceCollection.size()<<" fUseResidualImageInExtendedSearch="<<fUseResidualImageInExtendedSearch<<endl;
	Img* sourceResidualImg= 0;
	if(fUseResidualImageInExtendedSearch && fSourceCollection.size()>0){
		cout<<"SourceFinder::FindExtendedSource(): INFO: Get source residual image (only point-like sources removed)..."<<endl;
	
		Img* tmpImg= (Img*)fInputImg->Clone("tmpImg");
		tmpImg->SetNameTitle("tmpImg","tmpImg");
		if(fInputImg->HasMetaData()) tmpImg->SetMetaData(fInputImg->GetMetaData());
		if(fInputImg->HasBkgData()) {
			tmpImg->SetBkgData(fInputImg->GetBkgData());
			tmpImg->SetInterpolatedBkgLevelMap(fInputImg->GetInterpolatedBkgLevelMap());
			tmpImg->SetInterpolatedBkgRMSMap(fInputImg->GetInterpolatedBkgRMSMap());
		}
		tmpImg->ResetSources();
		tmpImg->SetSources(fSourceCollection);

		sourceResidualImg= tmpImg->GetSourceResidual(fUseLocalBackground,fSourceDilateKernelSize,fDilateNestedSources,fDilatedSourceType,fDilateSourceModel,fRandomizeInDilate,fRandSigmaInDilate);
		if(!sourceResidualImg) {
			cerr<<"SourceFinder::FindExtendedSource(): ERROR: Source residual image retrieval failed!"<<endl;
			if(tmpImg) tmpImg->Delete();
			if(inputImg) inputImg->Delete();
			return -1;
		}
		sourceResidualImg->SetNameTitle("sourceResidualImg","sourceResidualImg");
		if(tmpImg) tmpImg->Delete();

		imgToSegment= sourceResidualImg;
	}//close if use residual image
	

	//## Apply a smoothing stage?
	Img* smoothedImg= 0;
	if(fUsePreSmoothing){
		Img* smoothedInputImg= inputImg;
		if(fUseResidualImageInExtendedSearch && sourceResidualImg)  
			smoothedInputImg= sourceResidualImg;

		if(fSmoothingAlgo==eGaus){
			smoothedImg= smoothedInputImg->Smooth(fSmoothKernelSize,fSmoothKernelSize,fSmoothSigma,fSmoothSigma);
		}
		else if(fSmoothingAlgo==eGuided){
			cout<<"SourceFinder::FindExtendedSource(): INFO: Computing guided image filter..."<<endl;
			smoothedImg= smoothedInputImg->GetGuidedFilterImage(fGuidedSmoothRadius,fGuidedSmoothColorEps);
		}
		else{
			cerr<<"SourceFinder::FindExtendedSource(): ERROR: Invalid smoothing algo selected!"<<endl;
			if(inputImg) inputImg->Delete();
			if(sourceResidualImg) sourceResidualImg->Delete();
			return -1;
		}

		if(!smoothedImg){
			cerr<<"SourceFinder::FindExtendedSource(): ERROR: Source residual image retrieval failed!"<<endl;
			if(inputImg) inputImg->Delete();
			if(sourceResidualImg) sourceResidualImg->Delete();
			return -1;
		}
		smoothedImg->SetNameTitle("smoothedImg","smoothedImg");

		imgToSegment= smoothedImg;
	}//close if use smoothing
	

	int status= 0;
	if(fExtendedSearchMethod==eWT){
		status= FindExtendedSource_WT(imgToSegment);
	}
	else if(fExtendedSearchMethod==eSPSegmentation){
		//status= FindExtendedSource_SPSegmentation(imgToSegment);
		status= FindExtendedSource_SPSegmentation_new(imgToSegment);
	}
	else if(fExtendedSearchMethod==eCVSegmentation){
		status= FindExtendedSource_CVContour(imgToSegment);
	}
	else{
		cerr<<"SourceFinder::FindExtendedSource(): ERROR: Invalid extended source method ("<<fExtendedSearchMethod<<") selected!"<<endl;
		if(inputImg) inputImg->Delete();
		if(sourceResidualImg) sourceResidualImg->Delete();
		if(smoothedImg) smoothedImg->Delete();
		return -1;
	}

	if(status<0){
		cerr<<"SourceFinder::FindExtendedSource(): ERROR: Extended search method failed!"<<endl;
		if(inputImg) inputImg->Delete();
		if(sourceResidualImg) sourceResidualImg->Delete();
		if(smoothedImg) smoothedImg->Delete();
		return -1;
	}

	cout<<"SourceFinder::FindExtendedSource(): INFO: Clearing up..."<<endl;
	if(inputImg){
		delete inputImg;
		inputImg= 0;
	}
	if(sourceResidualImg){
		delete sourceResidualImg;
		sourceResidualImg= 0;
	}
	if(smoothedImg){
		delete smoothedImg;
		smoothedImg= 0;
	}

	return 0;

}//close FindExtendedSource()


int SourceFinder::FindExtendedSource_WT(Img* inputImg){

	if(!inputImg){
		cerr<<"SourceFinder::FindExtendedSource_WT(): ERROR: Null ptr to input image given!"<<endl;
		return -1;
	}
	
	//## Find extended sources in the W3, W5 scales of the residual image where ONLY POINT-LIKE SOURCES are removed
	cout<<"SourceFinder::FindExtendedSource_WT(): INFO: Find extended sources in the residual image WT-"<<fWTScaleForExtendedSourceSearch<<"  scale ..."<<endl;
	std::vector<Img*> wt_extended= inputImg->GetWaveletDecomposition(fWTScaleForExtendedSourceSearch);
	

	int status= FindSource(wt_extended[fWTScaleForExtendedSourceSearch],fSeedThreshold,fMergeThreshold);
	if(status<0){
		cerr<<"SourceFinder::FindExtendedSource_WT(): ERROR: Extended source finding failed!"<<endl;
		return -1;
	}
	std::vector<Source*> ExtendedSourceCollection= wt_extended[fWTScaleForExtendedSourceSearch]->GetSources();
	int nSources= (int)ExtendedSourceCollection.size();		

	//Get segmented image
	fSegmentedImage= inputImg->GetSourceMask(ExtendedSourceCollection,false);
	fSegmentedImage->SetNameTitle("segmentedImg","segmentedImg");

	fFinalSignalSegmentedImage= inputImg->GetSourceMask(ExtendedSourceCollection,true);
	fFinalSignalSegmentedImage->SetNameTitle("segmentedImg_final_binary","segmentedImg_final_binary");


	cout<<"SourceFinder::FindExtendedSource_WT(): INFO: "<<ExtendedSourceCollection.size()<<" extended sources found!"<<endl;
	if(nSources>0) {
		//## Apply source selection?
		int nSelSources= nSources;

		if(fApplySourceSelection){
			nSelSources= 0;

			for(int i=0;i<nSources;i++){
				//Is bad source (i.e. line-like blob, etc...)?
				if(!IsGoodSource(ExtendedSourceCollection[i])) {
					cout<<"SourceFinder::FindExtendedSource_WT(): INFO: Source no. "<<i<<" tagged as bad source, skipped!"<<endl;	
					ExtendedSourceCollection[i]->fIsGoodSource= false;
					continue;
				}
			
				//Is point-like source?
				if( IsCompactPointLike(ExtendedSourceCollection[i]) ){
					cout<<"SourceFinder::FindExtendedSource_WT(): INFO: Candidate bright source tagged as a point-like source ..."<<endl;
					ExtendedSourceCollection[i]->SetType(Source::ePointLike);
				}
				//else{//extended source, search for nested?
					ExtendedSourceCollection[i]->fType= Source::eExtended;
					if(fSearchNestedSources) {
						cout<<"SourceFinder::FindExtendedSource_WT(): INFO: Finding nested source in extended source: "<<ExtendedSourceCollection[i]->fName<<endl;
						int nestedSourceSearchStatus= FindNestedSource(ExtendedSourceCollection[i]);
						if(nestedSourceSearchStatus<0) cerr<<"SourceFinder::FindExtendedSource(): WARN: Nested source search failed!"<<endl;
					}	

					//Add source to the list	
					fSourceCollection.push_back(ExtendedSourceCollection[i]);
					nSelSources++;
				//}//close else
				
			}//end loop sources

		}//close if source selection
		else{
			//Add extended sources to the list
			fSourceCollection.insert(fSourceCollection.end(),ExtendedSourceCollection.begin(),ExtendedSourceCollection.end());
		}

		cout<<"SourceFinder::FindExtendedSource_WT(): INFO: Added "<<nSelSources<<" extended sources to the list..."<<endl;

		//Draw extended sources
		//cout<<"SourceFinder::FindExtendedSource(): INFO: Drawing "<<ExtendedSourceCollection.size()<<" extended sources in input image..."<<endl;
		//DrawSources(wt_extended[fWTScaleForExtendedSourceSearch]);

	}//close if nsources>0

	//## Clear-up	
	cout<<"SourceFinder::FindExtendedSource_WT(): INFO: Clearing up..."<<endl;
	/*
	if(inputImg){
		delete inputImg;
		inputImg= 0;
	}
	if(sourceResidualImg){
		delete sourceResidualImg;
		sourceResidualImg= 0;
	}
	*/

	return 0;
	
}//close SourceFinder::FindExtendedSource()

/*
int SourceFinder::FindExtendedSource_SPSegmentation(Img* inputImg){

	if(!inputImg){
		cerr<<"SourceFinder::FindExtendedSource_SPSegmentation(): ERROR: Null ptr to input image given!"<<endl;
		return -1;
	}
	
	//## Perform segmentation	
  int status= fSLICSegmentation->RunSegmentation(inputImg,fSPSize,fSPRegularization,fSPMinArea,true,fSPMergingAlgo);
	if(status<0){
		cerr<<"SourceFinder::FindExtendedSource_SPSegmentation(): ERROR: Segmentation failed!"<<endl;
		return -1;
	}

	//## Retrieve results	
	std::vector<Region*> regions= fSLICSegmentation->GetRegions();//list of segmented regions
	fSegmentedImage= fSLICSegmentation->GetClusterColoredImage(inputImg,regions);//Get segmented image
	if(!fSegmentedImage){
		cerr<<"SourceFinder::FindExtendedSource_SPSegmentation(): ERROR: Cannot retrieve segmented image!"<<endl;
		return -1;
	}
	fSegmentedImage->SetNameTitle("segmentedImg","segmentedImg");

	fFinalSegmentedImage= fSLICSegmentation->GetClusterColoredImage(inputImg,regions,true);//Get segmented image (only selected)
	if(!fFinalSegmentedImage){
		cerr<<"SourceFinder::FindExtendedSource_SPSegmentation(): ERROR: Cannot retrieve final segmented image!"<<endl;
		return -1;
	}
	fFinalSegmentedImage->SetNameTitle("segmentedImg_final","segmentedImg_final");

	//Saliency map
	fSaliencyMap= fSLICSegmentation->ComputeSaliencyMap(regions);//saliency @ final hierarchy level
	if(!fSaliencyMap){
		cerr<<"SourceFinder::FindExtendedSource_SPSegmentation(): ERROR: Cannot retrieve saliency map!"<<endl;
		return -1;
	}
	fSaliencyMap->SetNameTitle("saliencyMap","saliencyMap");
	
	fInitSaliencyMap= fSLICSegmentation->GetSaliencyMap();//init saliency 
	if(!fInitSaliencyMap){
		cerr<<"SourceFinder::FindExtendedSource_SPSegmentation(): ERROR: Cannot retrieve init saliency map!"<<endl;
		return -1;
	}
	fInitSaliencyMap->SetNameTitle("saliencyMap_init","saliencyMap_init");
	
	fSumSaliencyMap= fSLICSegmentation->GetSumSaliencyMap();
	if(!fSumSaliencyMap){
		cerr<<"SourceFinder::FindExtendedSource_SPSegmentation(): ERROR: Cannot retrieve sum saliency map!"<<endl;
		return -1;
	}
	fSumSaliencyMap->SetNameTitle("saliencyMap_sum","saliencyMap_sum");
	

	//Merging info
	std::vector<VLSlicSegmentation::SPMergingInfo> mergingInfo= fSLICSegmentation->GetMergingInfo();
	
	for(unsigned int k=0;k<mergingInfo.size();k++){
		fLevelId= mergingInfo[k].levelId;
		fNR= mergingInfo[k].NR;
		fDissMin= mergingInfo[k].DissMin;
		fDissMax= mergingInfo[k].DissMax;	
		fDissMin_merged= mergingInfo[k].MergedDissMin;
		fDissMax_merged= mergingInfo[k].MergedDissMax;			
		fDissMedian= mergingInfo[k].DissMedian;
		fDissMedianRMS= mergingInfo[k].DissMedianRMS;	
		fDissMedian0= mergingInfo[k].DissMedian0;
		fDissMedianRMS0= mergingInfo[k].DissMedianRMS0;	
		fMSE= mergingInfo[k].MSE;
		fNDissEntries= (mergingInfo[k].DissList).size();
		fNDissEntries_merged= (mergingInfo[k].MergedDissList).size();	
		for(int j=0;j<fNDissEntries;j++) fDissList[j]= (mergingInfo[k].DissList)[j];	
		for(int j=0;j<fNDissEntries_merged;j++) fDissList_merged[j]= (mergingInfo[k].MergedDissList)[j];
		fSPMergingInfo->Fill();
	}

	//## Create sources from segmented regions
	for(unsigned int i=0;i<regions.size();i++){
		bool isSignificant= regions[i]->fIsSignificative;
		if(!isSignificant) continue;
		Source* aSource= regions[i]->GetSource();
		if(!aSource) continue;
		fSourceCollection.push_back(aSource);
	}//end loop regions

	
	//## Draw results
	TCanvas* SPSegmentationPlot= new TCanvas("SPSegmentationPlot","SPSegmentationPlot");
	SPSegmentationPlot->cd();
	if(fSegmentedImage) fSegmentedImage->Draw("COLZ");
	//if(regionContours) regionContours->Draw("Psame");
	
	TCanvas* SPSaliencyPlot= new TCanvas("SPSaliencyPlot","SPSaliencyPlot");
	SPSaliencyPlot->cd();
	if(fSaliencyMap) fSaliencyMap->Draw("COLZ");
	

	//## Clear-up	
	cout<<"SourceFinder::FindExtendedSource(): INFO: Clearing up..."<<endl;
	//if(inputImg){
	//	delete inputImg;
	//	inputImg= 0;
	//}
	//if(sourceResidualImg){
	//	delete sourceResidualImg;
	//	sourceResidualImg= 0;
	//}
	
	return 0;

}//close FindExtendedSource_SPSegmentation()
*/



int SourceFinder::FindExtendedSource_SPSegmentation_new(Img* inputImg){

	if(!inputImg){
		cerr<<"SourceFinder::FindExtendedSource_SPSegmentation(): ERROR: Null ptr to input image given!"<<endl;
		return -1;
	}

	
	//## Compute stats and local bkg
	if(!inputImg->HasStats()){
		if(inputImg->ComputeStats(true,false,true)<0) {
			cerr<<"SourceFinder::FindExtendedSource_SPSegmentation_new(): ERROR: Cannot get input image stats!"<<endl;
			return -1;
		}
	}

	//## Compute local image bkg
	if(inputImg->ComputeLocalBkg((Img::LocalBkgMethod)fLocalBkgModel, (Img::BkgMethod)fBkgEstimator,fBkgBoxSizeX,fBkgBoxSizeY,fBkgGridSizeX,fBkgGridSizeY)<0) {
		cerr<<"SourceFinder::FindExtendedSource_SPSegmentation_new(): ERROR: Failed to compute local bkg!"<<endl;
		return -1;
	}
	

	//## Perform segmentation	
  int status= fSLICSegmenter->RunSegmentation(inputImg,fSPSize,fSPRegularization,fSPMinArea,true);
	if(status<0){
		cerr<<"SourceFinder::FindExtendedSource_SPSegmentation(): ERROR: Segmentation failed!"<<endl;
		return -1;
	}

	//## Retrieve results
	std::vector<Region*> regions= fSLICSegmenter->GetRegions();//list of segmented regions
	std::vector<Region*> regions_merged= fSLICSegmenter->GetMergedRegions();//list of merged regions
	std::vector< std::vector<long> > labels= fSLICSegmenter->GetPixelClusterIds();//pixel labels
	std::vector< std::vector<long> > labels_merged= fSLICSegmenter->GetMergedPixelClusterIds();//pixel merged labels

	//Superpixel segmented image	
	fSegmentedImage= SLICUtils::GetSegmentedImage(inputImg,regions,false,false,true);
	if(!fSegmentedImage){
		cerr<<"SourceFinder::FindExtendedSource_SPSegmentation(): ERROR: Cannot retrieve segmented image!"<<endl;
		return -1;
	}
	fSegmentedImage->SetNameTitle("segmentedImg","segmentedImg");

	
	//Final segmented image
	fFinalSegmentedImage= SLICUtils::GetSegmentedImage(inputImg,regions_merged,false,false,true,false);
	if(!fFinalSegmentedImage){
		cerr<<"SourceFinder::FindExtendedSource_SPSegmentation(): ERROR: Cannot retrieve final segmented image!"<<endl;
		return -1;
	}
	fFinalSegmentedImage->SetNameTitle("segmentedImg_final","segmentedImg_final");

	//Binary final segmented image
	fFinalSignalSegmentedImage= SLICUtils::GetSegmentedImage(inputImg,regions_merged,true,false,true,true);
	if(!fFinalSignalSegmentedImage){
		cerr<<"SourceFinder::FindExtendedSource_SPSegmentation(): ERROR: Cannot retrieve final segmented image!"<<endl;
		return -1;
	}
	fFinalSignalSegmentedImage->SetNameTitle("segmentedImg_final_binary","segmentedImg_final_binary");

	// Compute region contour info
	SLICUtils::SLICContourData* contourData= SLICUtils::ComputeBoundaryContours(inputImg,labels_merged,regions_merged);
	if(contourData) {
		Contour* contour= contourData->contour;
		if(contour) {
			fSegmentationContourGraph= contour->GetGraph();
			if(fSegmentationContourGraph) fSegmentationContourGraph->SetNameTitle("segmentationContour","segmentationContour");
		}
		delete contourData;
		contourData= 0;
	}
	

	//Saliency map
	fSaliencyMap= fSLICSegmenter->GetSaliencyMap();
	if(!fSaliencyMap){
		cerr<<"SourceFinder::FindExtendedSource_SPSegmentation(): ERROR: Cannot retrieve saliency map!"<<endl;
		return -1;
	}
	fSaliencyMap->SetNameTitle("saliencyMap","saliencyMap");
	
	/*
	//Merging info
	std::vector<SLICSegmenter::SPMergingInfo> mergingInfo= fSLICSegmenter->GetMergingInfo();
	
	for(unsigned int k=0;k<mergingInfo.size();k++){
		fLevelId= mergingInfo[k].levelId;
		fNR= mergingInfo[k].NR;
		fDissMin= mergingInfo[k].DissMin;
		fDissMax= mergingInfo[k].DissMax;		
		fDissMedian= mergingInfo[k].DissMedian;
		fDissMedianRMS= mergingInfo[k].DissMedianRMS;	
		fDissMedian0= mergingInfo[k].DissMedian0;
		fDissMedianRMS0= mergingInfo[k].DissMedianRMS0;	
		fMSE= mergingInfo[k].MSE;
		fSPMergingInfo->Fill();
	}
	*/

	//## Create sources from segmented regions
	std::vector<Source*> ExtSources;
	for(unsigned int i=0;i<regions_merged.size();i++){
		bool isSignificant= regions_merged[i]->fIsSignificative;
		if(!isSignificant) continue;
		Source* aSource= regions_merged[i]->GetSource();
		if(!aSource) continue;
		ExtSources.push_back(aSource);
		fSourceCollection.push_back(aSource);
	}//end loop regions

	bool findNestedSources= false;
	if(findNestedSources){
		for(unsigned int j=0;j<ExtSources.size();j++){
			Source* aSource= ExtSources[j];
		
			//Find nested sources
			int findNestedStatus= aSource->FindNestedClusterRegions(fSPSize,fSPRegularization,fSPMinArea,fSPMergingRatio,fSPMergingRegularization,fSPMergingMaxDissRatio2ndNeighbor);
			if(findNestedStatus<0){
				cerr<<"SourceFinder::FindExtendedSource_SPSegmentation(): WARN: Find nested source failed!"<<endl;
				continue;
			}
		}//end loop sources
	}//close if


	//## Draw results
	TCanvas* SPSegmentationPlot= new TCanvas("SPSegmentationPlot","SPSegmentationPlot");
	SPSegmentationPlot->cd();
	if(fSegmentedImage) fSegmentedImage->Draw("COLZ");
	if(fSegmentationContourGraph) {
		fSegmentationContourGraph->SetMarkerStyle(1);
		fSegmentationContourGraph->SetMarkerColor(kRed);
		fSegmentationContourGraph->Draw("P");
	}

	TCanvas* SPFinalSegmentationPlot= new TCanvas("SPFinalSegmentationPlot","SPFinalSegmentationPlot");
	SPFinalSegmentationPlot->cd();
	if(fFinalSegmentedImage) fFinalSegmentedImage->Draw("COLZ");
	
	TCanvas* SPSaliencyPlot= new TCanvas("SPSaliencyPlot","SPSaliencyPlot");
	SPSaliencyPlot->cd();
	if(fSaliencyMap) fSaliencyMap->Draw("COLZ");


	/*
	//## Get Source mask and clustered regions
	if(ExtSources.size()>0){
		Img* sourceMask= fInputImg->GetSourceMask(ExtSources,false);
		
		int status= fSLICSegmentation->RunSegmentation(sourceMask,fSPSize,fSPRegularization,fSPMinArea,true,fSPMergingAlgo);
		if(status<0){
			cerr<<"SourceFinder::FindExtendedSource_SPSegmentation_new(): ERROR: Signal segmentation failed!"<<endl;
			return -1;
		}

		// Retrieve results	
		std::vector<Region*> signal_regions= fSLICSegmentation->GetRegions();//list of segmented regions
		fFinalSignalSegmentedImage= fSLICSegmentation->GetClusterColoredImage(sourceMask,signal_regions);//Get segmented image
		if(!fFinalSignalSegmentedImage){
			cerr<<"SourceFinder::FindExtendedSource_SPSegmentation(): ERROR: Cannot retrieve segmented image!"<<endl;
			return -1;
		}
	
		TCanvas* SPFinalSignalSegmentationPlot= new TCanvas("SPFinalSignalSegmentationPlot","SPFinalSignalSegmentationPlot");
		SPFinalSignalSegmentationPlot->cd();
		if(fFinalSignalSegmentedImage) fFinalSignalSegmentedImage->Draw("COLZ");
			
	}//close if
	*/
	
	
	//## Clear-up	
	cout<<"SourceFinder::FindExtendedSource(): INFO: Clearing up..."<<endl;
	
	return 0;

}//close 



int SourceFinder::FindExtendedSource_CVContour(Img* inputImg){

	if(!inputImg){
		cerr<<"SourceFinder::FindExtendedSource_CVContour(): ERROR: Null ptr to input image given!"<<endl;
		return -1;
	}

	//## Perform segmentation
	//int status= fCVSegmentation->RunSegmentation(sourceResidualImg,fCVTimeStep,fCVWindowSize,fCVLambda1Par,fCVLambda2Par,fCVMuPar,fCVNuPar,fCVPPar);
	int status= fCVSegmentation->RunSegmentation(inputImg,fCVTimeStep,fCVWindowSize,fCVLambda1Par,fCVLambda2Par,fCVMuPar,fCVNuPar,fCVPPar);
	if(status<0){
		cerr<<"Img::FindExtendedSource_CVContour(): ERROR: ChanVese Segmentation failed!"<<endl;
		return -1;
	}

	
	//## Getting results
	Contour* cvContours= fCVSegmentation->GetContour();
	fSegmentedImage= fCVSegmentation->GetSegmentedImage();
	if(fSegmentedImage) fSegmentedImage->SetNameTitle("segmentedImg","segmentedImg");

	

	//## Finding blobs in masked image
	Img* maskedImg= inputImg->GetMask(fSegmentedImage,false);
	double fgValue= 1;	
	//int nPixMin= 2;
	status= maskedImg->FindCompactSource(fSegmentedImage,fgValue,fgValue,fNPixMin,false,false,false,false);
	if(status<0){
		cerr<<"SourceFinder::FindExtendedSource_CVContour(): ERROR: Finding blobs in Chan-Vese segmented mask failed!"<<endl;
		return -1;
	}
	std::vector<Source*> sources= maskedImg->GetSources();	
	Source* aSource= 0;
	for(unsigned int k=0;k<sources.size();k++){
		aSource= new Source;
		*aSource= *(sources[k]);	
		fSourceCollection.push_back(aSource);
	}
	if(maskedImg) maskedImg->Delete();

	//## Draw results
	TCanvas* CVSegmentationPlot= new TCanvas("CVSegmentationPlot","CVSegmentationPlot");
	CVSegmentationPlot->cd();
	if(fSegmentedImage) fSegmentedImage->Draw("COLZ");
	if(cvContours) {
		TGraph* contourGraph= cvContours->GetGraph();
		if(contourGraph) contourGraph->Draw("Psame");
	}

	//## Clear-up	
	cout<<"SourceFinder::FindExtendedSource_CVContour(): INFO: Clearing up..."<<endl;
	/*
	if(inputImg){
		delete inputImg;
		inputImg= 0;
	}
	if(sourceResidualImg){
		delete sourceResidualImg;
		sourceResidualImg= 0;
	}
	*/

	return 0;

}//close FindExtendedSource_CVContour()


int SourceFinder::FindSource(Img* image,double seedThreshold,double mergeThreshold){

	Img* zimg= 0;

	//## Use curvature mixed image?
	if(fUseCurvatureMixture){
		cout<<"SourceFinder::FindSource(): INFO: Use curvature mixed image..."<<endl;
		zimg= image->GetLaplacianWeightedImage(fCurvatureWeight,"LoG");
		zimg->SetNameTitle("zimg","zimg");
	}//close if
	else {
		zimg= (Img*)image->Clone("zimg");
	}

	//## Compute stats
	cout<<"SourceFinder::FindSource(): INFO: Computing stats..."<<endl;
	if(!zimg->ComputeStats(true,false,true)<0){
		cerr<<"SourceFinder::FindSource(): ERROR: Stats computing failed!"<<endl;
		if(zimg) zimg->Delete();
		return -1;
	}
	zimg->DumpStats();		
	
	//## Compute Bkg
	cout<<"SourceFinder::FindSource(): INFO: Computing bkg..."<<endl;
	int bkgStatus= 0;
	if(fUseLocalBackground) bkgStatus= zimg->ComputeLocalBkg((Img::LocalBkgMethod)fLocalBkgModel, (Img::BkgMethod)fBkgEstimator,fBkgBoxSizeX,fBkgBoxSizeY,fBkgGridSizeX,fBkgGridSizeY);
	else bkgStatus= zimg->ComputeBkg((Img::BkgMethod)fBkgEstimator);
	if(bkgStatus<0) {
		cerr<<"SourceFinder::FindSource(): ERROR: Bkg computing failed!"<<endl;
		if(zimg) zimg->Delete();
		return -1;
	}
	zimg->DumpBkg();
		
	//## Get significance map
	cout<<"SourceFinder::FindSource(): INFO: Computing significance map..."<<endl;
	Img* significanceMap= zimg->GetSignificanceMap(fUseLocalBackground);
	if(!significanceMap){
		cerr<<"SourceFinder::FindSource(): ERROR: Cannot get significance map!"<<endl;	
		if(zimg) zimg->Delete();
		return -1;
	}

	TString canvasName= Form("SignificancePlot_%s",image->GetName());
	TCanvas* SignificancePlot= new TCanvas(canvasName,canvasName);
	SignificancePlot->cd();
	significanceMap->Draw("COLZ");
	
	//## Find sources
	cout<<"SourceFinder::FindSource(): INFO: Finding sources..."<<endl;
	bool findNegativeExcess= false;
	bool findNestedSources= fSearchNestedSources;
	bool mergeBelowSeed= false;
	int status= image->FindCompactSource(significanceMap,seedThreshold,mergeThreshold,fNPixMin,fSearchNegativeExcess,fSearchNestedSources,fUseLocalBackground,mergeBelowSeed,fPeakThreshold);
	if(status<0) {
		cerr<<"SourceFinder::FindSource(): ERROR: Source finding failed!"<<endl;
		if(zimg) zimg->Delete();
		return -1;
	}

	//## Clearup
	cerr<<"SourceFinder::FindSource(): INFO: Clearing up zimg..."<<endl;	
	if(zimg) {
		delete zimg;
		zimg= 0;
	}
	cout<<"done!"<<endl;
	
	return 0;

}//close FindSource()


int SourceFinder::FindNestedSource(Source* aSource){

	if(!aSource) return -1;
	
	std::vector<Source*> NestedSources= aSource->fNestedSourceCollection;
	int nNestedSources= (int)(NestedSources.size());
	cout<<"SourceFinder::FindNestedSource(): INFO: Processing no. "<<nNestedSources<<" nested sources for source id="<<aSource->fId<<"..."<<endl;
	
	for(int i=0;i<NestedSources.size();i++) {
		//Apply source selection
		if(fApplySourceSelection){
			if(!IsGoodSource(NestedSources[i])){
				cout<<"SourceFinder::FindNestedSource(): INFO: Nested blob no. "<<i<<" tagged as bad source, skipped!"<<endl;	
				NestedSources[i]->fIsGoodSource= false;
				continue;
			}		
			//Point-like source tagger
			if(IsCompactPointLike(NestedSources[i]) ){
				cout<<"SourceFinder::FindNestedSource(): INFO: Nested blob no. "<<i<<" tagged as point-like!"<<endl;	
				NestedSources[i]->fType= Source::ePointLike;
			}
		}//close if apply source selection
	}//end loop sources

	
	cout<<"SourceFinder::FindNestedSource(): INFO: "<<nNestedSources<<" nested sources found after selection cuts!"<<endl;

	

	/*
	//Reset current nested source collection
	aSource->ResetNestedSources();
	
	//## Get source image
	Img* sourceImg= aSource->GetImage(Source::eFluxMap);
	Img* curvImg= aSource->GetImage(Source::eFluxCurvMap);

	//## Find peaks in curvature map
	int tol= 1;//1 pixel-tolerance
	TGraph* peaks= curvImg->FindPeaks(tol);
	if(!peaks) {
		cout<<"SourceFinder::FindNestedSource(): INFO: No peaks detected!"<<endl;
		sourceImg->Delete();
		curvImg->Delete();
		return 0;
	}
	int nPeaks= peaks->GetN();
	cout<<"SourceFinder::FindNestedSource(): INFO: #"<<nPeaks<<" peaks detected!"<<endl;

	//## Find blobs image by thresholding the curvature map and then by flood-fill
	Img* curvImg_binary= curvImg->GetBinarized(fCurvatureThreshold);
	Img* maskedImg= sourceImg->GetMask(curvImg_binary,false);

	for(int k=0;k<nPeaks;k++){
		double x, y;
		peaks->GetPoint(k,x,y);
		int peakBinId= sourceImg->FindBin(x,y);
		if(sourceImg->IsBinOverflow(peakBinId) || sourceImg->IsBinUnderflow(peakBinId) ) continue;
		curvImg_binary->AddBinContent(peakBinId,1);
	}	
	double fgValue= 1;	
	bool findNegativeExcess= false;
	bool findNestedSources= false;
	bool useLocalBackground= false;
	bool mergeBelowSeed= true;
	int status= maskedImg->FindCompactSource(curvImg_binary,fgValue+1,fgValue,fSourceComponentMinNPix,findNegativeExcess,findNestedSources,useLocalBackground,mergeBelowSeed);
	if(status<0){
		cerr<<"SourceFinder::FindNestedSource(): WARN: Blob flood-filling failed!"<<endl;
		sourceImg->Delete();
		curvImg->Delete();
		peaks->Delete();	
		curvImg_binary->Delete();
		maskedImg->Delete();
		return -1;
	}

	//## Set nested sources
	std::vector<Source*> NestedSources= maskedImg->GetSources();
	if(NestedSources.size()<=0){
		cout<<"Source::FindNestedSource(): INFO: No nested source found!"<<endl;
		sourceImg->Delete();
		curvImg->Delete();
		peaks->Delete();	
		curvImg_binary->Delete();
		maskedImg->Delete();
		return 0;
	}
	cout<<"Source::FindNestedSource(): INFO: #"<<NestedSources.size()<<" nested blobs found..."<<endl;	
	
	for(unsigned int i=0;i<NestedSources.size();i++) {

		//Apply source selection
		if(fApplySourceSelection){
			if(!IsGoodSource(NestedSources[i])){
				cout<<"SourceFinder::FindNestedSource(): INFO: Nested blob no. "<<i<<" tagged as bad source, skipped!"<<endl;	
				continue;
			}		
			
			//Point-like source tagger
			if(IsCompactPointLike(NestedSources[i]) ){
				cout<<"SourceFinder::FindNestedSource(): INFO: Nested blob no. "<<i<<" tagged as point-like!"<<endl;	
				NestedSources[i]->fType= Source::ePointLike;
			}
			//else{//Skip non point-like sources
			//	cout<<"SourceFinder::FindNestedSource(): INFO: Nested source no. "<<i<<" skipped as non point-like..."<<endl;
			//	continue;
			//}
		}//close if apply source selection

		//Add nested source to input source
		aSource->AddNestedSource(NestedSources[i]);
	}//end loop sources

	int nNestedSources= (int)( (aSource->fNestedSourceCollection).size() );
	cout<<"SourceFinder::FindNestedSource(): INFO: "<<nNestedSources<<" nested sources found after selection cuts!"<<endl;

	//## Clear-up
	if(sourceImg) {
		delete sourceImg;
		sourceImg= 0;
	}
	if(curvImg) {
		delete curvImg;
		curvImg= 0;
	}
	if(peaks) peaks->Delete();	
	if(curvImg_binary) {
		delete curvImg_binary;
		curvImg_binary= 0;
	}
	if(maskedImg) {
		delete maskedImg;
		maskedImg= 0;
	}
	*/	

	return 0;

}//close SourceFinder::FindNestedSource()


int SourceFinder::DrawSources(Img* image,bool drawSourcePlots){

	cout<<"SourceFinder::DrawSources(): INFO: Drawing sources..."<<endl;
	if(!image) return -1;

	Img* image_norm= image->GetNormalizedImage(1,256);
	TString imageName= image->GetName();
	
	//## Draw image
	TString canvasName= Form("Plot_%s",std::string(imageName).c_str());
	TCanvas* Plot= new TCanvas(canvasName,canvasName);
	Plot->cd();

	if(image_norm) image_norm->Draw("COLZ");

	//cout<<"SourceFinder::DrawSources(): INFO: Get source list..."<<endl;
	//std::vector<Source*> SourceCollection= image->GetSources();

	cout<<"== DUMP SOURCES =="<<endl;
	//for(unsigned int i=0;i<SourceCollection.size();i++) {
	for(unsigned int i=0;i<fSourceCollection.size();i++) {
		Source* thisSource= fSourceCollection[i];
		if(thisSource) {	
			thisSource->Dump(true,true);

			int thisSourceType= thisSource->fType;
			std::vector<Contour*> contours= thisSource->GetContours();

			for(unsigned int k=0;k<contours.size();k++) {
				TGraph* thisContourGraph= contours[k]->GetGraph();
				if(thisContourGraph) {
					if(thisSourceType==Source::eCompact) thisContourGraph->SetLineColor(kBlack);
					else if(thisSourceType==Source::ePointLike) thisContourGraph->SetLineColor(kRed);
					else if(thisSourceType==Source::eExtended) thisContourGraph->SetLineColor(kGreen);
					else if(thisSourceType==Source::eExtendedSegm) thisContourGraph->SetLineColor(kBlue-4);
					
					thisContourGraph->SetLineWidth(2);
					thisContourGraph->Draw("Lsame");
				}
			}//end loop contours
		}//close if source

		//Draw nested sources
		std::vector<Source*> NestedSources= thisSource->fNestedSourceCollection;
		for(unsigned int k=0;k<NestedSources.size();k++){
			int thisNestedSourceType= NestedSources[k]->fType;
			std::vector<Contour*> nestedContours= NestedSources[k]->GetContours();	
			for(unsigned int l=0;l<nestedContours.size();l++) {
				TGraph* thisContourGraph= nestedContours[l]->GetGraph();
				if(thisContourGraph) {
					if(thisNestedSourceType==Source::eCompact) thisContourGraph->SetLineColor(kBlack);
					else if(thisNestedSourceType==Source::ePointLike) thisContourGraph->SetLineColor(kRed);
					else if(thisNestedSourceType==Source::eExtended) thisContourGraph->SetLineColor(kGreen);
					else if(thisNestedSourceType==Source::eExtendedSegm) thisContourGraph->SetLineColor(kBlue-4);		
					thisContourGraph->SetLineWidth(2);
					thisContourGraph->Draw("Lsame");
				}
			}//end loop contours
		}//end loop nested sources
	}//end loop sources

	
	/*
	//## Draw sources
	if(drawSourcePlots){
		cout<<"== DUMP BRIGHT SOURCES =="<<endl;
		cout<<"INFO: "<< SourceCollection.size()<<" bright sources found!"<<endl;
		for(int i=0;i<SourceCollection.size();i++){
			if(SourceCollection[i]) {
				SourceCollection[i]->Dump(true,true);
				//SourceCollection[i]->Draw(Source::eSourceSignificanceMap);
				//SourceCollection[i]->Draw(Source::eFluxCurvMixtureSignificanceMap,0.7);
			}
		}
	}//close if
	*/

	//## Draw residual image
	Img* residualImg= image->GetSourceResidual(fUseLocalBackground,fSourceDilateKernelSize,fDilateNestedSources,fDilatedSourceType,fDilateSourceModel,fRandomizeInDilate,fRandSigmaInDilate);

	TCanvas* ResidualPlot= new TCanvas("ResidualPlot","ResidualPlot");
	ResidualPlot->cd();
	residualImg->Draw("COLZ");

	cout<<"DrawSources(): INFO: End source drawing..."<<endl;
	
	return 0;

}//close DrawSources()



int SourceFinder::ReadImage(){

	//--> ROOT reading
	if(fInputFileExtension=="root"){// Read image from ROOT file
		TFile* inputFile = new TFile(fInputFileName.c_str(),"READ");
		if(!inputFile || inputFile->IsZombie()){
			cerr<<"SourceFinder::ReadImage(): ERROR: Cannot open input file "<<fInputFileName<<"!"<<endl;
			return -1;
		}
		fInputImg=  (Img*)inputFile->Get("FullImage");
		if(!fInputImg){
			cerr<<"SourceFinder::ReadImage(): ERROR: Cannot get image from input file "<<fInputFileName<<"!"<<endl;
			return -1;
		}
	}//close if

	//--> FITS reading
	else if(fInputFileExtension=="fits"){// Read image from FITS file
		ImgFITSReader reader(fInputFileName.c_str());	
		reader.ReadHeader();

		ImgFITSReader::FITSHeader header= reader.GetHeaderInfo();
		int Nx= header.Nx;
		int Ny= header.Ny;
		double Bmaj= header.Bmaj;
		double Bmin= header.Bmin;
		double dX= header.dX;
		double dY= header.dY;
		int nBeamPix= fabs(Bmaj/dX);
		cout<<"SourceFinder::ReadImage(): INFO: Header info: Nx="<<Nx<<" Ny="<<Ny<<" nBeamPix="<<nBeamPix<<endl;
	
		if(fInputImg) fInputImg->Delete();
		fInputImg= new Img;
		reader.Read(*fInputImg);
	}//close else if

	//--> NATURAL image reading
	else if(fInputFileExtension=="png" || fInputFileExtension=="jpg" || fInputFileExtension=="bmp"){//Read natural image
		if(fInputImg) fInputImg->Delete();
		fInputImg= new Img;
		if(fInputImg->ReadFile(fInputFileName)<0){
			cerr<<"SourceFinder::ReadImage(): ERROR: Failed to read natural image!"<<endl;
			return -1;
		}
	}	

	//--> Invalid extension
	else{
		cerr<<"SourceFinder::ReadImage(): ERROR: Invalid file extension detected!"<<endl;
		return -1;
	}
	fInputImg->SetNameTitle("img","img");	

	//## Compute stats
	cout<<"SourceFinder::ReadImage(): INFO: Computing input image stats..."<<endl;
	if(!fInputImg->ComputeStats(true,false,true)<0){
		cerr<<"SourceFinder::ReadImage(): ERROR: Stats computing failed!"<<endl;
		return -1;
	}
	fInputImg->DumpStats();		

	//## Compute Bkg
	int bkgStatus= 0;
	if(fUseLocalBackground) bkgStatus= fInputImg->ComputeLocalBkg((Img::LocalBkgMethod)fLocalBkgModel, (Img::BkgMethod)fBkgEstimator, fBkgBoxSizeX,fBkgBoxSizeY, fBkgGridSizeX, fBkgGridSizeY);
	else bkgStatus= fInputImg->ComputeBkg((Img::BkgMethod)fBkgEstimator);
	if(bkgStatus<0) {
		cerr<<"SourceFinder::ReadImage(): ERROR: Bkg computing failed!"<<endl;
		return -1;
	}
	fInputImg->DumpBkg();

	return 0;

}//close ReadImage()


void SourceFinder::Save(){

	cout<<"SourceFinder::Save(): INFO: Storing results to file & catalog..."<<endl;

	if(!fOutputFile || !fOutputTree || !fDS9CatalogFilePtr) return;
	fOutputFile->cd();

	cout<<"SourceFinder::Save(): INFO: Saving DS9 region header..."<<endl;
	fprintf(fDS9CatalogFilePtr,"global color=red font=\"helvetica 12 normal\" edit=1 move=1 delete=1 include=1\n");
	fprintf(fDS9CatalogFilePtr,"image\n");

	fprintf(fDS9CatalogFilePtr2,"global color=blue font=\"helvetica 12 normal\" edit=1 move=1 delete=1 include=1\n");
	fprintf(fDS9CatalogFilePtr2,"image\n");

	cout<<"SourceFinder::Save(): INFO: Saving "<<fSourceCollection.size()<<" sources to file..."<<endl;
	for(unsigned int k=0;k<fSourceCollection.size();k++){
		cout<<"SourceFinder::Save(): INFO: Dumping DS9 region info for source no. "<<k<<" ..."<<endl;
		std::string regionInfo= "";
		if(fDS9RegionFormat==1) regionInfo= fSourceCollection[k]->DumpDS9RegionInfo(true);
		else if(fDS9RegionFormat==2) regionInfo= fSourceCollection[k]->DumpDS9EllipseRegionInfo(true);
		else continue;

		fprintf(fDS9CatalogFilePtr,"%s\n",regionInfo.c_str());

		std::string regionInfo2= fSourceCollection[k]->DumpDS9EllipseRegionInfo(true);
	  fprintf(fDS9CatalogFilePtr2,"%s\n",regionInfo2.c_str());
	  
		cout<<"SourceFinder::Save(): INFO: Set source ptr to source "<<k<<"..."<<endl;	
		fSource= fSourceCollection[k];
		cout<<"SourceFinder::Save(): INFO: Saving source no. "<<k<<" to tree..."<<endl;	
		fOutputTree->Fill();
	}//end loop sources

	cout<<"SourceFinder::Save(): INFO: Closing DS9 file region..."<<endl;	
	fclose(fDS9CatalogFilePtr);
	fclose(fDS9CatalogFilePtr2);

	cout<<"SourceFinder::Save(): INFO: Writing tree to file..."<<endl;	
	fOutputTree->Write();

	//Save image to file?
	if(fSaveImageToFile){
		if(fSaveImageType==eInputImage) fInputImg->Write();
		else if(fSaveImageType==eResidualImage){
			Img* residualImg= fInputImg->GetSourceResidual(fUseLocalBackground,fSourceDilateKernelSize,fDilateNestedSources,fDilatedSourceType,fDilateSourceModel,fRandomizeInDilate,fRandSigmaInDilate);
			if(residualImg) residualImg->Write();
		}
		else if(fSaveImageType==eSegmentedImage && fSearchExtendedSources){
			if(fSegmentedImage) fSegmentedImage->Write();
			if(fFinalSegmentedImage) fFinalSegmentedImage->Write();
			if(fFinalSignalSegmentedImage) fFinalSignalSegmentedImage->Write();
			if(fSaliencyMap) fSaliencyMap->Write();
			if(fSumSaliencyMap) fSumSaliencyMap->Write();
			if(fInitSaliencyMap) fInitSaliencyMap->Write();
			if(fSegmentationContourGraph) fSegmentationContourGraph->Write();
		}
		else{
			cerr<<"SourceFinder::Save(): WARN: Invalid save image flag specified!"<<endl;
		}

		if(fSPMergingInfo) fSPMergingInfo->Write();
		if(fConfigInfo) fConfigInfo->Write();
	}//close if save image to file

	cout<<"SourceFinder::Save(): INFO: Closing output file..."<<endl;	
	fOutputFile->Close();

	cout<<"SourceFinder::Save(): INFO: End save to file"<<endl;
		
}//close Save()

void SourceFinder::SetGraphicsStyle(){
	
	if(myStyle) myStyle->Delete();
	myStyle= new TStyle("myStyle","myStyle");
	myStyle->SetCanvasDefH(700); 
  myStyle->SetCanvasDefW(700); 

	myStyle->SetFrameBorderMode(0);
	myStyle->SetCanvasBorderMode(0);
  myStyle->SetPadBorderMode(0);
  myStyle->SetPadColor(0);
  myStyle->SetCanvasColor(0);
  myStyle->SetTitleFillColor(0);
  myStyle->SetTitleBorderSize(1);
  myStyle->SetStatColor(0);
  myStyle->SetStatBorderSize(1);
  myStyle->SetOptTitle(0);
  myStyle->SetOptStat(0);
  myStyle->SetOptFit(1);
	myStyle->SetOptLogx(0);
	myStyle->SetOptLogy(0);
  //myStyle->SetPalette(1,0);
  myStyle->SetTitleBorderSize(0);//border size of Title PavelLabel
  myStyle->SetTitleX(0.1f);
	myStyle->SetTitleW(0.8f);
  myStyle->SetStatY(0.975);                
  myStyle->SetStatX(0.95);                
  myStyle->SetStatW(0.2);                
  myStyle->SetStatH(0.15);                
  myStyle->SetTitleXOffset(0.8);
  myStyle->SetTitleYOffset(1.1);
  myStyle->SetMarkerStyle(8);
  myStyle->SetMarkerSize(0.4);
  myStyle->SetFuncWidth(1.);
  myStyle->SetPadTopMargin(0.1);
  myStyle->SetPadBottomMargin(0.12);
  myStyle->SetPadLeftMargin(0.15);
  myStyle->SetPadRightMargin(0.15);
  myStyle->SetTitleSize(0.06,"X");
  myStyle->SetTitleSize(0.06,"Y");
  myStyle->SetTitleSize(0.06,"Z");
  myStyle->SetTitleFont(52,"X");
  myStyle->SetTitleFont(52,"Y");
  myStyle->SetTitleFont(52,"Z");
  myStyle->SetLabelFont(42,"X");
  myStyle->SetLabelFont(42,"Y");
  myStyle->SetLabelFont(42,"Z");

  myStyle->SetErrorX(0.);
	
	gROOT->SetStyle("myStyle");

	/*
	const int ncolors= 2000;	
	const int Number= 12;
	int myTemperaturePalette[ncolors];
  double r[]    = {0.164, 0.150, 0.250, 0.450, 0.670, 0.880, 1.000, 1.000,1.000,0.970,0.850,0.650};
  double g[]    = {0.043, 0.306, 0.630, 0.853, 0.973, 1.000,1.000,0.880,0.679,0.430,0.150,0.000};
  double b[]    = {0.850,1.000,1.000,1.000,1.000,1.000,0.750,0.600,0.450,0.370,0.196,0.130};
  double stop[] = {0.,0.2,0.3,0.4,0.5,0.7,0.75,0.8,0.85,0.9,0.95,1};
  int FI = TColor::CreateGradientColorTable(Number, stop, r, g, b, ncolors);
  for (int i=0;i<ncolors;i++) {
		myTemperaturePalette[i] = FI+i;
	}
	*/

	/*
  int HotToColdPalette[ncolors];
  double r_HotToCold[Number];
  double g_HotToCold[Number];
  double b_HotToCold[Number];
  //double stop_HotToCold[] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.5,0.6,0.7,0.8,1};
	double stop_HotToCold[] = {0.,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1};

  for (int i=0;i<Number;i++) {
   r_HotToCold[i]= r[Number-i-1];
   g_HotToCold[i]= g[Number-i-1];
   b_HotToCold[i]= b[Number-i-1];
  }
  int FI_HotToCold = TColor::CreateGradientColorTable(Number, stop_HotToCold, r_HotToCold, g_HotToCold, b_HotToCold, ncolors);
  for (int i=0;i<ncolors;i++) {
  	HotToColdPalette[i] = FI_HotToCold + i;
  }
	*/


	gStyle->SetNumberContours(999);
	//gStyle->SetPalette(ncolors,myTemperaturePalette);

	//gStyle->SetPalette(53);//Black Body
	//gStyle->SetPalette(52);//Black & White
	gStyle->SetPalette(55);//Rainbow


}//close SetGraphicsStyle()

