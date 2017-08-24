#include <SLICSegmenter.h>
#include <SLICUtils.h>
#include <OutlierDetector.h>
#include <HClust.h>

#include <TColor.h>
#include <TMath.h>
#include <TStyle.h>
#include <TColor.h>
#include <TPad.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TText.h>
#include <TPolyLine.h>
#include <Math/QuantFuncMathCore.h>

#include <mathop.h>

#include <chrono>
using namespace std;
using namespace std::chrono;

#define SMALL_NUMBER 0.00000000001;

SLICSegmenter::SLICSegmenter() {

	fLUT= 0;
	fImg= 0;
	fInputImg= 0;
	fLaplImg= 0;
	fEdgeFilterImg= 0;
	fKirschEdgeFilterImg= 0;
	fRegions.clear();
	fMergedRegions.clear();
	fUseLogNormalizedImage= true;
	fSPMergingRatio= 0.3;
	fSPMergingRegularization= 0.5;
	fMinMergedSP= 1;
	fSPMergingDistThreshold= 0.1;
	fPixelRatioCut= 0.5;
	fSignificantSPRatio= 0.5;
	
	fSPMergingMaxDissRatio= 1.15;	
	fSPMergingMaxDissRatio_2ndNeighbor= 1.15;
	fSPSize= 20;
	fSPRegularization= 0.01;

	fSPMergingIncludeSpatialPars= false;
	fSPMergingUseAdaptingDistThreshold= false;
	fSPMergingAdaptingDistThresholdScale= 500;

	fSPMergingAggloMethod= 2;
	fSPMergingMinClustSize= 3;
	fSPMergingMaxHeightQ= 0.95;
	fSPMergingDeepSplitLevel= 1;
		
	fSPMergingUseCurvature= true;
	fSPMergingEdgeModel= 1;
	fUse2ndNeighborsInSPMerging= true;

	fSaliencyImg= 0;
	fPreMergingSegmentedImg= 0;

	//Saliency options
	fSaliencyThresholdFactor= 2;
	fBkgSaliencyThresholdFactor= 1;
	fSaliencyImgThresholdFactor= 1;
	fSaliencyMinReso= 20;
	fSaliencyMaxReso= 60;
	fSaliencyResoStepSize= 10;
	fSaliencyUseRobustPars= true;
	fSaliencyUseBkgMap= true;
	fSaliencyUseNoiseMap= true;
	fSaliencyUseCurvatureMap= false;
	fSaliencyNNFactor= 0.1;
	fSaliencyFilterThresholdFactor= 0.8;
	fSaliencyNormalizationMode= 2;

	fBoxSizeX= 90;
	fBoxSizeY= 90;
	fGridSizeX= 18;
	fGridSizeY= 18;

	fNPixMin= 3000;

	//Chan-Vese 
	fCVSegmentation= 0;
	fCVTimeStep= 0.1;
	fCVWindowSize= 1;	
	fCVLambda1Par= 1;
	fCVLambda2Par= 2;
	fCVMuPar= 0.5;
	fCVNuPar= 0;
	fCVPPar= 1;

}//close constructor


SLICSegmenter::~SLICSegmenter() {
	
	cout<<"SLICSegmenter::~SLICSegmenter(): INFO: Deleting regions..."<<endl;
	for(unsigned int i=0;i<fRegions.size();i++){
		if(fRegions[i]){
			delete fRegions[i];
			fRegions[i]= 0;
		}
	}	
	fRegions.clear();

	for(unsigned int i=0;i<fMergedRegions.size();i++){
		if(fMergedRegions[i]){
			delete fMergedRegions[i];
			fMergedRegions[i]= 0;
		}
	}	
	fMergedRegions.clear();
	
	cout<<"SLICSegmenter::~SLICSegmenter(): INFO: Deleting filter images..."<<endl;
	if(fLUT) {
		delete fLUT;
		fLUT= 0;
	}
	if(fLaplImg) {
		//fLaplImg->Delete();
		delete fLaplImg;
		fLaplImg= 0;
	}
	if(fEdgeFilterImg) {
		//fEdgeFilterImg->Delete();
		delete fEdgeFilterImg;
		fEdgeFilterImg= 0;
	}
	if(fKirschEdgeFilterImg) {
		//	fKirschEdgeFilterImg->Delete();
		delete fKirschEdgeFilterImg;
		fKirschEdgeFilterImg= 0;
	}
	
	cout<<"SLICSegmenter::~SLICSegmenter(): INFO: Deleting norm image..."<<endl;
	if(fImg) fImg->Delete();
	
	fPixelClusterIds.clear();

	cout<<"SLICSegmenter::~SLICSegmenter(): INFO: Deleting Chan-Vese..."<<endl;
	if(fCVSegmentation){
		delete fCVSegmentation;
		fCVSegmentation= 0;
	}

	
}//close constructor


int SLICSegmenter::Init(Img* inputImage){

	//## Clear previous data (if any), and re-initialize it
  fPixelClusterIds.clear();
	fRegions.clear();
	fMergedRegions.clear();
	if(fImg) fImg->Delete();
	if(fLaplImg) fLaplImg->Delete();
	if(fEdgeFilterImg) fEdgeFilterImg->Delete();
	if(fKirschEdgeFilterImg) fKirschEdgeFilterImg->Delete();
	if(fLUT) fLUT->Delete();

	//## Check input image
	if(!inputImage) {
		cerr<<"SLICSegmenter::Init(): ERROR: Null input image ptr!"<<endl;
		return -1;
	}
	fInputImg= inputImage;
	
	//## Compute image stats
	if(!fInputImg->HasStats()){
		if(fInputImg->ComputeStats(true,false,true)<0) {
			cerr<<"SLICSegmenter::Init(): ERROR: Cannot get input image stats!"<<endl;
			return -1;
		}
	}
	
	
	//## Normalize image to range [1-256] and define color LUT for pixel mapping	
	double normmin= 1;
	double normmax= 256;
	if(fUseLogNormalizedImage) fImg= inputImage->GetLogNormalizedImage(normmin,normmax,true);
	else fImg= inputImage->GetNormalizedImage(normmin,normmax,true);

	int nlevels= 255;
	double zmin = fImg->GetMinimum();
	double zmax = fImg->GetMaximum();
	if(!fLUT) {
		//Find if histo with the same name exist
		//bool alreadyExist= (gROOT->FindObject("LUT")!=nullptr);
		fLUT= new TH1D("","",nlevels,zmin,zmax);//do not give a name!
	}
	
	//## Compute laplacian images
	fLaplImg= fInputImg->GetLaplacianImage(true);
	if(!fLaplImg->HasStats()){
		if(fLaplImg->ComputeStats(true,false,true)<0) {
			cerr<<"SLICSegmenter::Init(): ERROR: Cannot get input laplacian image stats!"<<endl;
			return -1;
		}
	}
	
	//## Compute edge image
	fKirschEdgeFilterImg= fInputImg->GetKirschImage();

	if(fSPMergingEdgeModel==1){//Kirsch edge model
		fEdgeFilterImg= fInputImg->GetKirschImage();	
	}
	else if(fSPMergingEdgeModel==2){
		if(!fCVSegmentation) fCVSegmentation= new ChanVeseSegmentation;
		int status= fCVSegmentation->RunSegmentation(fInputImg,fCVTimeStep,fCVWindowSize,fCVLambda1Par,fCVLambda2Par,fCVMuPar,fCVNuPar,fCVPPar);
		if(status<0){
			cerr<<"SLICSegmenter::Init(): ERROR: ChanVese Segmentation failed!"<<endl;
			return -1;
		}
		fEdgeFilterImg= fCVSegmentation->GetContourImage();
	}//close else if
	else{
		cerr<<"SLICSegmenter::Init(): ERROR: Invalid edge model selected!"<<endl;
		return -1;
	}

	if(!fEdgeFilterImg->HasStats()){
		if(fEdgeFilterImg->ComputeStats(true,false,true)<0) {
			cerr<<"SLICSegmenter::Init(): ERROR: Cannot get input edge filtered image stats!"<<endl;
			return -1;
		}
	}

	if(!fKirschEdgeFilterImg->HasStats()){
		if(fKirschEdgeFilterImg->ComputeStats(true,false,true)<0) {
			cerr<<"SLICSegmenter::Init(): ERROR: Cannot get input Kirsch edge filtered image stats!"<<endl;
			return -1;
		}
	}
	
	return 0;

}//close Init()


int SLICSegmenter::RunSegmentation(Img* inputImg,int regionSize,double regularization,int minRegionSize,bool mergeRegions){

	//## Check input image
	if(!inputImg) {
		cerr<<"SLICSegmenter::RunSegmentation(): ERROR: Null ptr to input image given!"<<endl;
		return -1;
	}

	//## Initialize data for segmentation tasks (LUT, vectors/lists, ...)
	//## Transform image for segmentation purposes
	cout<<"SLICSegmenter::RunSegmentation(): INFO: Initializing segmentation task..."<<endl;
	if( Init(inputImg)<0 ){
		cerr<<"SLICSegmenter::RunSegmentation(): ERROR: Init failed!"<<endl;	
		return -1;
	}

	//## Perform SLIC segmentation 
	cout<<"SLICSegmenter::RunSegmentation(): INFO: Generate superpixels..."<<endl;
	fSPSize= regionSize;
	fSPRegularization= regularization;
	fSPMinArea= minRegionSize;
	if(SuperpixelGenerator(regionSize, regularization, minRegionSize)<0){
		cerr<<"SLICSegmenter::RunSegmentation(): ERROR: Superpixel generation failed!"<<endl;
		return -1;
	}
	
	//## Merge the superpixels?
	if(mergeRegions) {
		if(MultiStepSPMerger()<0) {
			cerr<<"SLICSegmenter::RunSegmentation(): ERROR: Hierarchical superpixel merging failed!"<<endl;
			return -1;
		}
	}

	return 0;

}//close RunSegmentation()



int SLICSegmenter::SuperpixelGenerator(int regionSize,double regularization, int minRegionSize){

	if(!fImg){
		cerr<<"SLICSegmenter::SuperpixelGenerator(): INFO: Null ptr to image!"<<endl;
		return -1;
	}
	
	//## Init data
	int Nx= fImg->GetNbinsX();
	int Ny= fImg->GetNbinsY();
	int numChannels= 1;
	std::vector< std::vector<bool> > istaken;
	
	float* image = new float[Nx*Ny];//input image
	
	vl_uint32* segmentation = new vl_uint32[Nx*Ny];//vector where the segmentation is stored

	vl_index i, x, y, u, v, k, region;
  vl_uindex iter;
  vl_size const numRegionsX = (vl_size) ceil( (double)Nx/regionSize) ;
  vl_size const numRegionsY = (vl_size) ceil( (double)Ny/regionSize) ;
  vl_size const numRegions = numRegionsX * numRegionsY ;
  vl_size const numPixels = Nx * Ny ;
  float* centers;
  float* edgeMap;
  float previousEnergy = VL_INFINITY_F;
  float startingEnergy= 0;
  vl_uint32* masses ;
  vl_size const maxNumIterations = 100;

  assert(segmentation);
  assert(image);
  assert(Nx >= 1);
  assert(Ny >= 1);
  assert(numChannels >= 1);
  assert(regionSize >= 1);
  assert(regularization >= 0);

#define atimage(x,y,k) image[(x)+(y)*Nx+(k)*Nx*Ny]
#define atEdgeMap(x,y) edgeMap[(x)+(y)*Nx]

  edgeMap = (float*)vl_calloc(numPixels, sizeof(float)) ;
  masses = (vl_uint32*)vl_malloc(sizeof(vl_uint32) * numPixels) ;
  centers = (float*)vl_malloc(sizeof(float) * (2 + numChannels) * numRegions) ;
	
	//## Fill input image
	for (int i=0;i<Nx; i++) { 
		std::vector<long int> cr;
		std::vector<bool> nb;

		for (int j=0;j<Ny; j++) {		
			double w= fImg->GetBinContent(i+1,j+1);
			image[i + j*Nx]= w;
			cr.push_back(-1);
			nb.push_back(false);
		}//end loop Ny
		fPixelClusterIds.push_back(cr);
		istaken.push_back(nb);
	}//end loop Nx

	//## Generate the superpixels
	//## The algorithm will store the final segmentation in a one-dimensional array.        
	
  //## Compute edge map (gradient strength)
	//cout<<"SLICSegmenter::RunSLICSegmentation(): INFO: Computing edge map..."<<endl;
  for (k = 0 ; k < (signed)numChannels ; ++k) {
    for (y = 1 ; y < (signed)Ny-1 ; ++y) {
      for (x = 1 ; x < (signed)Nx-1 ; ++x) {
        float a = atimage(x-1,y,k) ;
        float b = atimage(x+1,y,k) ;
        float c = atimage(x,y+1,k) ;
        float d = atimage(x,y-1,k) ;
        atEdgeMap(x,y) += (a - b)  * (a - b) + (c - d) * (c - d) ;
      }
    }
  }

  //## Initialize K-means centers
  i = 0;
	cout<<"SLICSegmenter::SuperpixelGenerator(): INFO: Initializing K-means centers..."<<endl;
 
  for (v = 0 ; v < (signed)numRegionsY ; ++v) {
    for (u = 0 ; u < (signed)numRegionsX ; ++u) {
      vl_index xp;
      vl_index yp;
      vl_index centerx = 0;
      vl_index centery = 0;
      float minEdgeValue = VL_INFINITY_F;

      x = (vl_index) vl_round_d(regionSize * (u + 0.5));
      y = (vl_index) vl_round_d(regionSize * (v + 0.5));

      x = VL_MAX(VL_MIN(x, (signed)Nx-1),0) ;
      y = VL_MAX(VL_MIN(y, (signed)Ny-1),0) ;

      // search in a 3x3 neighbourhood the smallest edge response
      for (yp = VL_MAX(0, y-1) ; yp <= VL_MIN((signed)Ny-1, y+1) ; ++ yp) {
        for (xp = VL_MAX(0, x-1) ; xp <= VL_MIN((signed)Nx-1, x+1) ; ++ xp) {
          float thisEdgeValue = atEdgeMap(xp,yp) ;
          if (thisEdgeValue < minEdgeValue) {
            minEdgeValue = thisEdgeValue ;
            centerx = xp ;
            centery = yp ;
          }
        }
      }

      // initialize the new center at this location
      centers[i++] = (float) centerx ;
      centers[i++] = (float) centery ;
      for (k  = 0 ; k < (signed)numChannels ; ++k) {
        centers[i++] = atimage(centerx,centery,k) ;
      }
    }
  }

  //## Run k-means iterations
	cout<<"SLICSegmenter::SuperpixelGenerator(): INFO: Running "<<maxNumIterations<<" iterations..."<<endl;
 
  for (iter = 0 ; iter < maxNumIterations ; ++iter) {
		
    float factor = regularization / (regionSize * regionSize) ;
    float energy = 0 ;

    // assign pixels to centers
    for (y = 0 ; y < (signed)Ny ; ++y) {
      for (x = 0 ; x < (signed)Nx ; ++x) {
        vl_index u = floor((double)x / regionSize - 0.5) ;
        vl_index v = floor((double)y / regionSize - 0.5) ;
        vl_index up, vp ;
        float minDistance = VL_INFINITY_F ;

        for (vp = VL_MAX(0, v) ; vp <= VL_MIN((signed)numRegionsY-1, v+1) ; ++vp) {
          for (up = VL_MAX(0, u) ; up <= VL_MIN((signed)numRegionsX-1, u+1) ; ++up) {
            vl_index region = up  + vp * numRegionsX ;
            float centerx = centers[(2 + numChannels) * region + 0]  ;
            float centery = centers[(2 + numChannels) * region + 1] ;
            float spatial = (x - centerx) * (x - centerx) + (y - centery) * (y - centery) ;
            float appearance = 0 ;
            float distance ;
            for (k = 0 ; k < (signed)numChannels ; ++k) {
              float centerz = centers[(2 + numChannels) * region + k + 2]  ;
              float z = atimage(x,y,k) ;
              appearance += (z - centerz) * (z - centerz) ;
            }
            distance = appearance + factor * spatial ;
            if (minDistance > distance) {
              minDistance = distance ;
              segmentation[x + y * Nx] = (vl_uint32)region ;
            }
          }
        }
        energy += minDistance ;
      }
    }

    // check energy termination conditions
    if (iter == 0) {
      startingEnergy = energy ;
    } 
		else {
      if ((previousEnergy - energy) < 1e-5 * (startingEnergy - energy)) {
        break ;
      }
    }
    previousEnergy = energy ;

    // recompute centers
    memset(masses, 0, sizeof(vl_uint32) * Nx * Ny) ;
    memset(centers, 0, sizeof(float) * (2 + numChannels) * numRegions) ;

    for (y = 0 ; y < (signed)Ny ; ++y) {
      for (x = 0 ; x < (signed)Nx ; ++x) {
        vl_index pixel = x + y * Nx;
        vl_index region = segmentation[pixel] ;
        masses[region]++;
        centers[region * (2 + numChannels) + 0] += x;
        centers[region * (2 + numChannels) + 1] += y;
        for (k = 0 ; k < (signed)numChannels ; ++k) {
          centers[region * (2 + numChannels) + k + 2] += atimage(x,y,k) ;
        }
      }
    }

    for (region = 0 ; region < (signed)numRegions ; ++region) {
      float mass = VL_MAX(masses[region], 1e-8) ;
			for (i = (2 + numChannels)*region; i<(signed)(2 + numChannels)*(region + 1); ++i) {
        centers[i]/= mass;
			}
    }

  }//end loop iterations


	//Free stuff
  vl_free(masses);
  vl_free(centers);
  vl_free(edgeMap);

	//## Eliminate small regions
  cout<<"SLICSegmenter::SuperpixelGenerator(): INFO: Eliminating small regions..."<<endl;
  vl_uint32 * cleaned = (vl_uint32*)vl_calloc(numPixels, sizeof(vl_uint32));
  vl_uindex * segment = (vl_uindex*)vl_malloc(sizeof(vl_uindex) * numPixels);
  vl_size segmentSize;
  vl_uint32 label;
  vl_uint32 cleanedLabel;
  vl_size numExpanded;
  vl_index const dx [] = {+1, -1,  0,  0} ;
  vl_index const dy [] = { 0,  0, +1, -1} ;
  vl_index direction;
  vl_index pixel;

  for (pixel = 0 ; pixel < (signed)numPixels ; ++pixel) {
		if (cleaned[pixel]) continue;
		
    label = segmentation[pixel];
    numExpanded = 0 ;
    segmentSize = 0 ;
    segment[segmentSize++] = pixel;

		// Find an adjacent label, for possible use later
    // Find cleanedLabel as the label of an already cleaned region neighbour of this pixel
		cleanedLabel = label + 1;
    cleaned[pixel] = label + 1;
    x = pixel % Nx;
    y = pixel / Nx;
    for (direction = 0 ; direction < 4 ; ++direction) {
    	vl_index xp = x + dx[direction] ;
      vl_index yp = y + dy[direction] ;
      vl_index neighbor = xp + yp * Nx ;
      if (0<=xp && xp<(signed)Nx && 0<=yp && yp<(signed)Ny && cleaned[neighbor]) {
      	cleanedLabel = cleaned[neighbor] ;
      }
    }//end loop direction

    // Expand the segment	
		while (numExpanded < segmentSize) {
			vl_index open = segment[numExpanded++] ;
      x = open % Nx ;
      y = open / Nx ;
      for (direction = 0 ; direction < 4 ; ++direction) {
      	vl_index xp = x + dx[direction] ;
        vl_index yp = y + dy[direction] ;
        vl_index neighbor = xp + yp * Nx ;
        if (0<=xp && xp<(signed)Nx && 0<=yp && yp<(signed)Ny && cleaned[neighbor]==0 && segmentation[neighbor]==label) {
        	cleaned[neighbor] = label + 1 ;
          segment[segmentSize++] = neighbor ;
        }//close if
      }//end loop for
    }//end while loop

    //Change label to cleanedLabel if the segment is too small
		if (segmentSize < minRegionSize) {
			while (segmentSize > 0) {
      	cleaned[segment[--segmentSize]] = cleanedLabel ;
      }//end while loop
    }//close if 

  }//end pixel loop

	//Restore base 0 indexing of the regions
  for (pixel=0; pixel<(signed)numPixels; ++pixel) cleaned[pixel]--;

  memcpy(segmentation, cleaned, numPixels * sizeof(vl_uint32)) ;
	
	//## Delete data
	//cout<<"SLICSegmenter::RunSLICSegmentation(): INFO: Delete allocated data..."<<endl;
	vl_free(cleaned);
  vl_free(segment);
	vl_free(image);

	//## Fill segmentation labels and compute stats/contours/...
	Img::StatsData* ImgStats= fInputImg->GetPixelStats();
	if(!ImgStats) return -1;
	Img::StatsData* EdgeFilterImgStats= fEdgeFilterImg->GetPixelStats();
	if(!EdgeFilterImgStats) return -1;
	Img::StatsData* LaplImgStats= fLaplImg->GetPixelStats();
	if(!LaplImgStats) return -1;

	for(unsigned int k=0;k<numRegions;k++) {
		fRegions.push_back( new Region() );
		fRegions[k]->fId= -1;
		fRegions[k]->fImageSizeX= Nx;
		fRegions[k]->fImageSizeY= Ny;
		fRegions[k]->fImageMinX= fInputImg->GetXaxis()->GetXmin();
		fRegions[k]->fImageMaxX= fInputImg->GetXaxis()->GetXmax();
		fRegions[k]->fImageMinY= fInputImg->GetYaxis()->GetXmin();
		fRegions[k]->fImageMaxY= fInputImg->GetYaxis()->GetXmax();
		fRegions[k]->fImageMinS= ImgStats->min;
		fRegions[k]->fImageMaxS= ImgStats->max;
		fRegions[k]->fImageMinScurv= LaplImgStats->min;
		fRegions[k]->fImageMaxScurv= LaplImgStats->max;
		fRegions[k]->fImageMinSedge= EdgeFilterImgStats->min;
		fRegions[k]->fImageMaxSedge= EdgeFilterImgStats->max;
		fRegions[k]->fImageRMS= ImgStats->rms;			
	}
	
	float colors[3];
  for (int ix=0; ix<Nx; ix++) {
		double binX= fImg->GetXaxis()->GetBinCenter(ix+1);
  	for (int iy=0;iy<Ny; iy++) {
			double binY= fImg->GetYaxis()->GetBinCenter(iy+1);
			int id= ix + iy * Nx;
			int pixelId= fImg->GetBin(ix+1,iy+1);
			double S= fInputImg->GetBinContent(ix+1,iy+1);
			double binContent= fImg->GetBinContent(ix+1,iy+1);//already in log scale is log mapping is activated	
			
			int LUTBinId= fLUT->FindBin(binContent);
			int colorId= TColor::GetColorPalette(LUTBinId); 
			gROOT->GetColor(colorId)->GetRGB(*colors,*(colors+1),*(colors+2));//colors are in range [0,1]		
			for(int kk=0;kk<3;kk++) colors[kk]*= 255;//convert colors to 0-255 
			
			int pixelLabel= (int)segmentation[id];

			double S_curv= fLaplImg->GetBinContent(ix+1,iy+1);
			double S_edge= fEdgeFilterImg->GetBinContent(ix+1,iy+1);
    	fPixelClusterIds[ix][iy] = pixelLabel;
			
			if(pixelLabel<0 || pixelLabel>=(int)numRegions){
				cout<<"SLICSegmenter::SuperpixelGenerator(): WARNING: Skip this pixel label: "<<pixelLabel<<endl;
				continue;
			}

			//Create and fill a pixel	
			Region::Pixel aPixel;
			aPixel.id= pixelId;
			aPixel.S= S;
			aPixel.color= TVector3(colors[0],colors[1],colors[2]);	
			aPixel.x= binX;
			aPixel.y= binY;
			aPixel.ix= ix;
			aPixel.iy= iy;
			aPixel.isOnEdge= false;
			aPixel.distanceToEdge= 1.e+99;
			aPixel.S_curv= S_curv;
			aPixel.S_edge= S_edge;
		
			//Add pixel to region
			fRegions[pixelLabel]->fId= pixelLabel;
			if(S!=0) fRegions[pixelLabel]->AddPixel(aPixel);
			
    }//end loop Ny
  }//end loop Nx

	//cout<<"SLICSegmenter::RunSLICSegmentation(): INFO: Freeing segmentation..."<<endl;
	free(segmentation);
	

	//## Remove regions without pixels (important in case of empty image zones)
	std::vector<int> regionsToBeDeleted;
	std::map<int,int> goodRegionIdMap;
	for(unsigned int k=0;k<fRegions.size();k++) {	
		int regionId= fRegions[k]->fId;
		int nPix= fRegions[k]->fNPix;
		goodRegionIdMap.insert( std::pair<int,int>(regionId,1) );
		if(nPix<=0) {
			regionsToBeDeleted.push_back(k);
			goodRegionIdMap[regionId]= -1;
		}
	}
	for (int ix=0; ix<Nx; ix++) {
		for (int iy=0;iy<Ny; iy++) {
			int labelId= fPixelClusterIds[ix][iy];
			int corrFactor= goodRegionIdMap[labelId];
			fPixelClusterIds[ix][iy]*= corrFactor;
		}
	}
	Utils::DeleteItems(fRegions, regionsToBeDeleted);
	//##########################################################################

	cout<<"SLICSegmenter::SuperpixelGenerator(): INFO: nRegions after cleanup: "<<fRegions.size()<<endl;

	//## Compute region parameters
	for(unsigned int k=0;k<fRegions.size();k++) {	
		int regionId= fRegions[k]->fId;
		int nPix= fRegions[k]->fNPix;
		if(nPix<=0) continue;
		
		//Compute stats params
		fRegions[k]->ComputeParameters(false,true,false);

		double Cx= fRegions[k]->fX0;
		double Cy= fRegions[k]->fY0;
		double Cz= fRegions[k]->fMean;
	
		//Compute RGB corresponding to mean
		double regionMean= fRegions[k]->fMean;
		double regionMean_lg= log10(regionMean);
		int LUTBinId= fLUT->FindBin(regionMean);
		if(fUseLogNormalizedImage) LUTBinId= fLUT->FindBin(regionMean_lg);//use log mapping
		int colorId= TColor::GetColorPalette(LUTBinId); 
		gROOT->GetColor(colorId)->GetRGB(*colors,*(colors+1),*(colors+2));//colors are in range [0,1]		
		for(int kk=0;kk<3;kk++) colors[kk]*= 255;//convert colors to 0-255 

		//Set region color
		(fRegions[k]->fColor).SetXYZ(colors[0],colors[1],colors[2]);
			
		//Print regions
		//fRegions[k]->Dump();

		//cout<<"SLICSegmenter::SuperpixelGenerator(): INFO: Region no. "<<fRegions[k]->fId<<": N="<<nPix<<" C("<<Cx<<","<<Cy<<","<<Cz<<") color("<<fRegions[k]->fColor.X()<<","<<fRegions[k]->fColor.Y()<<","<<fRegions[k]->fColor.Z()<<"), connection(";
		//for(unsigned int l=0;l<(fRegions[k]->fNeighbourRegions).size();l++) cout<<(fRegions[k]->fNeighbourRegions)[l]<<",";
		//cout<<")"<<endl;
	}//end loop regions

	cout<<"SLICSegmenter::SuperpixelGenerator(): INFO: End superpixel generation..."<<endl;
	
	/*
	//## Get segmented image
	cout<<"SLICSegmenter::SuperpixelGenerator(): INFO: Get segmented image..."<<endl;
	Img* segmentedImg= SLICUtils::GetSegmentedImage(fInputImg,fRegions,false,false,true);
	if(segmentedImg) segmentedImg->SetNameTitle("initSegmentedImg","initSegmentedImg");

	//## Get text labels
	std::vector<TText*> regionTextList= SLICUtils::GetRegionTextLabels(fRegions);

	//## Compute region contour info
	cout<<"SLICSegmenter::SuperpixelGenerator(): INFO: Get contour data..."<<endl;
	SLICUtils::SLICContourData* contourData= SLICUtils::ComputeBoundaryContours(fInputImg,fPixelClusterIds,fRegions);

	//## Compute similarity matrix
	cout<<"SLICSegmenter::SuperpixelGenerator(): INFO: Get similarity data..."<<endl;
	SLICUtils::SLICSimilarityData* similarityData= SLICUtils::ComputeDissimilarityMatrix(fEdgeFilterImg,contourData,fRegions,fSPMergingRegularization,fSPMergingIncludeSpatialPars);
	std::vector<double> saliencyList= similarityData->saliencyList;
	

	SLICUtils::SLICSimilarityData* dissimilarityData= SLICUtils::ComputeRegionSimilarity(fEdgeFilterImg,contourData,fRegions,fSPMergingRegularization,fSPMergingIncludeSpatialPars);
	fAbsDissMedian= dissimilarityData->Dmedian;
 	fAbsDissMedianRMS= dissimilarityData->Dmedianrms;
 	if(dissimilarityData){
		delete dissimilarityData;
		dissimilarityData= 0;
	}
	*/
	
	//## Compute saliency map	
	/*
	//fSaliencyImg=	SLICUtils::GetSaliencyMap(fInputImg,fRegions,100);
	//fSaliencyImg=	SLICUtils::GetSaliencyMap(fInputImg,similarityData->SaliencyDissimilarityMatrix,fRegions,100);
	TString imgName= Form("%s_saliency",fInputImg->GetName());
	fSaliencyImg= (Img*)fInputImg->Clone(imgName);
	fSaliencyImg->SetNameTitle(imgName,imgName);
	fSaliencyImg->Reset();
	for(unsigned int i=0;i<saliencyList.size();i++){
		double Saliency= saliencyList[i];
		int nPixelsInRegion= (int)fRegions[i]->fNPix;
		cout<<"Region no. "<<i<<" (id="<<fRegions[i]->fId<<") nPix="<<nPixelsInRegion<<", Saliency="<<Saliency<<endl;
		for(int j=0;j<nPixelsInRegion;j++){//loop on pixels inside region
			int thisPixelId= (fRegions[i]->fPixelCollection)[j].id;
			fSaliencyImg->SetBinContent(thisPixelId,Saliency);
		}//end loop pixels in region	
	}//end loop regions
	*/
	//fSaliencyImg= fInputImg->GetSaliencyMap(fRegions,100);


	/*
	fSaliencyImg= fInputImg->GetSaliencyMap(fRegions,100);
	

	TCanvas* SaliencyPlot= new TCanvas("SaliencyPlot","SaliencyPlot");
	SaliencyPlot->cd();
	if(fSaliencyImg) fSaliencyImg->Draw("COLZ");

	//## Tag significative & salient regions
	cout<<"SLICSegmenter::SuperpixelGenerator(): INFO: Tagging superpixels..."<<endl;
	//SLICUtils::TagRegions(fInputImg,fRegions,fSaliencyThresholdFactor,fBkgSaliencyThresholdFactor,0.975,false,false);
	SLICUtils::TagRegions(fInputImg,fSaliencyImg,fRegions,fSaliencyThresholdFactor,fBkgSaliencyThresholdFactor,0.975,false,false);
	*/

	/*
	cout<<"SLICSegmenter::SuperpixelGenerator(): INFO: Drawing results..."<<endl;
	TCanvas* Plot= new TCanvas("Plot","Plot");
	Plot->cd();
	if(segmentedImg) segmentedImg->Draw("COLZ");
	*/
	/*
	for(unsigned int i=0;i<regionTextList.size();i++) {	
		if(regionTextList[i]) regionTextList[i]->Draw("same");
	}
	if(contourData && (contourData->contour) ){
		TGraph* contourGraph= (contourData->contour)->GetGraph();
		if(contourGraph) {
			contourGraph->SetLineColor(kRed);
			contourGraph->SetMarkerStyle(1);
			contourGraph->SetMarkerColor(kRed);
			contourGraph->Draw("P");
		}
	}
	*/

	return 0;

}//close SLICSegmenter::SuperpixelGenerator()




int SLICSegmenter::MultiStepSPMerger(){

	//## Check regions
	int nRegions= (int)fRegions.size();
	if(nRegions<=0) {
		cerr<<"SLICSegmenter::MultiStepSPMerger(): WARN: No regions available, nothing to be merged!"<<endl;
		return 0;
	}
	
	//## Init region list to be passed to algorithm
	cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: Copying initial superpixel partition to tmp partition..."<<endl;
	std::vector<Region*> regions;
	Region* aRegion= 0;
	for(int i=0;i<fRegions.size();i++){
		aRegion= new Region;
		*aRegion= *(fRegions[i]);
		regions.push_back(aRegion);
	}//end loop regions
	//regions.assign(fRegions.begin(),fRegions.end());

	std::map<int,int> regionIdMap_initialSegm;	
	for(int k=0;k<fRegions.size();k++){
		int regionId= fRegions[k]->fId;
		regionIdMap_initialSegm.insert( std::pair<int,int>(regionId,k) );
	}//end loop regions


	//## Compute saliency map 
	cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: Computing saliency map: reso pars("<<fSaliencyMinReso<<","<<fSaliencyMaxReso<<","<<fSaliencyResoStepSize<<"), knn="<<fSaliencyNNFactor<<", thr="<<fSaliencyFilterThresholdFactor<<" ..."<<endl;
	//fSaliencyImg= fInputImg->GetMultiResoSaliencyMap(fSaliencyMinReso,fSaliencyMaxReso,fSaliencyResoStepSize,fSPRegularization,fSPMinArea,fSaliencyNNFactor,fSaliencyUseRobustPars,fSPMergingUseCurvature,fSaliencyFilterThresholdFactor,fSaliencyUseCurvatureMap,fSaliencyUseBkgMap,fSaliencyUseNoiseMap,fSaliencyNormalizationMode,fSaliencyThresholdFactor,fSaliencyImgThresholdFactor);
	fSaliencyImg= fInputImg->GetMultiResoSaliencyMap(fSaliencyMinReso,fSaliencyMaxReso,fSaliencyResoStepSize,fSPRegularization,fSPMinArea,fSaliencyNNFactor,fSaliencyUseRobustPars,fSaliencyDissExpFalloffPar,fSaliencySpatialDistRegPar,fSaliencyFilterThresholdFactor,fSaliencyUseCurvatureMap,fSaliencyUseBkgMap,fSaliencyUseNoiseMap,fSaliencyNormalizationMode,fSaliencyThresholdFactor,fSaliencyImgThresholdFactor);
	

	//## Compute saliency stats
	if(!fSaliencyImg->ComputeStats(true,false,true)<0){
		cerr<<"SLICSegmenter::MultiStepSPMerger(): ERROR: Saliency stats computing failed!"<<endl;
		return -1;
	}	

	//## Compute saliency local background
	bool useTwoPass= false;
	fSaliencyImg->ComputeBkg(Img::eMedianBkg);
	//fSaliencyImg->ComputeLocalBkg(Img::eSuperpixelBkg,Img::eMedianBkg,fBoxSizeX,fBoxSizeY,fGridSizeX,fGridSizeY,fSPSize,fSPRegularization,fSPMinArea,useTwoPass);

	
	
	//## Tag regions
	cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: Tagging superpixels..."<<endl;
	SLICUtils::TagRegions_last(fInputImg,fSaliencyImg,regions,fSaliencyThresholdFactor,fBkgSaliencyThresholdFactor,fNPixMin);
	for(int i=0;i<fRegions.size();i++){
		fRegions[i]->fTag= regions[i]->fTag;
	}//end loop regions
	

	//## COmpute the pre-merging segmented image
	//## 3 channel: 1: bkg, 2: unknown, 3: signal
	fPreMergingSegmentedImg= (Img*)fSaliencyImg->Clone("PreMergingSegmentedImg");
	fPreMergingSegmentedImg->Reset();

	for(unsigned int i=0;i<fRegions.size();i++){//loop on regions
		int regionId= fRegions[i]->fId;	
		int regionTag= fRegions[i]->fTag;
		long int nPixelsInRegion= regions[i]->fNPix;

		double binContent= 0;
		if(regionTag==Region::eBkgTag) binContent= 1;
		else if(regionTag==Region::eSignalTag) binContent= 3;
		else if(regionTag==Region::eUntagged) binContent= 2;
		
		for(long int j=0;j<nPixelsInRegion;j++){//loop on pixels inside region
			int thisPixelId= (fRegions[i]->fPixelCollection)[j].id;
			fPreMergingSegmentedImg->SetBinContent(thisPixelId,binContent);
		}//end loop pixels in region	
	}//end loop regions

	//## Create the mapping of regionId and vector index
	cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: Create the mapping of regionId and region list index..."<<endl;
	std::map<int,int> regionIdMap;
	std::map<int,int> regionIdMap_top;	
	std::vector<std::vector<int>> mergedRegionList;
	
	for(int k=0;k<nRegions;k++){
		int regionId= regions[k]->fId;
		regionIdMap.insert( std::pair<int,int>(regionId,k) );
		regionIdMap_top.insert( std::pair<int,int>(regionId,k) );
		mergedRegionList.push_back( std::vector<int>() );
	}//end loop regions
	
	//Init pixel labels: copy initial pixel labels to tmp list
	std::vector< std::vector<long int> > labels;
	for(unsigned int i=0;i<fPixelClusterIds.size();i++) {
		labels.push_back( std::vector<long int>() );
		labels[i].assign(fPixelClusterIds[i].begin(),fPixelClusterIds[i].end());

		fMergedPixelClusterIds.push_back(  std::vector<long int>() );
		fMergedPixelClusterIds[i].assign(fPixelClusterIds[i].begin(),fPixelClusterIds[i].end());
	}


	bool stopAlgo= false;
	int stageCounter= 0;
	int nTotMergedRegions= 0;
	int nBkg= 0;
	int nSig= 0;
	int nBoh= 0;
	int nBkg_allRegions= 0;
	int nSig_allRegions= 0;
	int nBoh_allRegions= 0;
	for(unsigned int i=0;i<regions.size();i++) {
		int tag= regions[i]->fTag;
		if(tag==Region::eBkgTag) nBkg++;
		else if(tag==Region::eSignalTag) nSig++;
		else if(tag==Region::eUntagged) nBoh++;
	}//end loop regions
	cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: Initial NR="<<regions.size()<<": (nBkg,nSig,nBoh)=("<<nBkg<<","<<nSig<<","<<nBoh<<")"<<endl;
			

	if(nSig>0){
	while(!stopAlgo){

		nTotMergedRegions= 0;
		cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: Stage no. "<<stageCounter<<", nTotMergedRegions="<<nTotMergedRegions<<endl;
		stageCounter++;
			
		//## Compute region tagging stats
		nBkg= 0;
		nSig= 0;
		nBoh= 0;	
		for(unsigned int i=0;i<regions.size();i++) {
			int tag= regions[i]->fTag;
			if(tag==Region::eBkgTag) nBkg++;
			else if(tag==Region::eSignalTag) nSig++;
			else if(tag==Region::eUntagged) nBoh++;
		}//end loop regions
		cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: NR="<<regions.size()<<": (nBkg,nSig,nBoh)=("<<nBkg<<","<<nSig<<","<<nBoh<<")"<<endl;
			
		//=== DEBUG ====
		nBkg_allRegions= 0;
		nSig_allRegions= 0;
		nBoh_allRegions= 0;	
		for(unsigned int i=0;i<fRegions.size();i++) {
			int tag= fRegions[i]->fTag;
			if(tag==Region::eBkgTag) nBkg_allRegions++;
			else if(tag==Region::eSignalTag) nSig_allRegions++;
			else if(tag==Region::eUntagged) nBoh_allRegions++;
		}//end loop regions
		cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: NR="<<fRegions.size()<<": (nBkg_allRegions,nSig_allRegions,nBoh_allRegions)=("<<nBkg_allRegions<<","<<nSig_allRegions<<","<<nBoh_allRegions<<")"<<endl;
		//=== DEBUG ====
		
		//## == STAGE 1==
		//## Merge non-tagged regions with bkg-tagged regions
		bool is1stStageEnd= false;
		int nMergedRegions_1stStage= 0;

		if(nBkg>0 && nBoh>0) {//start 1st stage if there are untagged regions present
		
			cout<<"==============="<<endl; 
			cout<<"=== STAGE 1 ==="<<endl;
			cout<<"==============="<<endl;
			cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: Start 1st stage..."<<endl;			
			nMergedRegions_1stStage= 0;
			int nregions_preStage1= (int)regions.size();
			SPMaxSimilarityMerger(regions,labels,Region::eBkgTag,Region::eUntagged);
			int nregions_postStage1= (int)regions.size();
			nMergedRegions_1stStage= nregions_preStage1-nregions_postStage1;
			nTotMergedRegions+= nMergedRegions_1stStage;
		}//close if	

		
		//## ReCompute region tagging stats
		nBkg= 0;
		nSig= 0;
		nBoh= 0;
		for(unsigned int i=0;i<regions.size();i++) {
			int tag= regions[i]->fTag;
			if(tag==Region::eBkgTag) nBkg++;
			else if(tag==Region::eSignalTag) nSig++;
			else if(tag==Region::eUntagged) nBoh++;
		}//end loop regions
		cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: After 1st stage: NR="<<regions.size()<<": (nBkg,nSig,nBoh)=("<<nBkg<<","<<nSig<<","<<nBoh<<")"<<endl;
		
		//=== DEBUG ====
		nBkg_allRegions= 0;
		nSig_allRegions= 0;
		nBoh_allRegions= 0;	
		for(unsigned int i=0;i<fRegions.size();i++) {
			int tag= fRegions[i]->fTag;
			if(tag==Region::eBkgTag) nBkg_allRegions++;
			else if(tag==Region::eSignalTag) nSig_allRegions++;
			else if(tag==Region::eUntagged) nBoh_allRegions++;
		}//end loop regions
		cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: After 1st stage: NR="<<fRegions.size()<<": (nBkg_allRegions,nSig_allRegions,nBoh_allRegions)=("<<nBkg_allRegions<<","<<nSig_allRegions<<","<<nBoh_allRegions<<")"<<endl;
		//=== DEBUG ====
		
		//## == STAGE 2 ==
		//## Merge non-tagged regions adaptively
		bool is2ndStageEnd= false;
		int nMergedRegions_2ndStage= 0;

		if(nBoh>0){//start 2nd stage if there are untagged regions available

			cout<<"==============="<<endl; 
			cout<<"=== STAGE 2 ==="<<endl;
			cout<<"==============="<<endl;
			cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: Start 2nd stage..."<<endl;
			nMergedRegions_2ndStage= 0;

			int nregions_preStage2= (int)regions.size();
			SPMaxSimilarityMerger(regions,labels,Region::eUntagged,Region::eUntagged);
			int nregions_postStage2= (int)regions.size();
			nMergedRegions_2ndStage= nregions_preStage2-nregions_postStage2;
			nTotMergedRegions+= nMergedRegions_2ndStage;	
		}//close if


		//## Check end condition
		nBkg= 0;
		nSig= 0;
		nBoh= 0;
		for(unsigned int i=0;i<regions.size();i++) {
			int tag= regions[i]->fTag;
			if(tag==Region::eBkgTag) nBkg++;
			else if(tag==Region::eSignalTag) nSig++;
			else if(tag==Region::eUntagged) nBoh++;
		}//end loop regions
		cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: After 2nd stage: NR="<<regions.size()<<": (nBkg,nSig,nBoh)=("<<nBkg<<","<<nSig<<","<<nBoh<<")"<<endl;

		//=== DEBUG ====
		nBkg_allRegions= 0;
		nSig_allRegions= 0;
		nBoh_allRegions= 0;	
		for(unsigned int i=0;i<fRegions.size();i++) {
			int tag= fRegions[i]->fTag;
			if(tag==Region::eBkgTag) nBkg_allRegions++;
			else if(tag==Region::eSignalTag) nSig_allRegions++;
			else if(tag==Region::eUntagged) nBoh_allRegions++;
		}//end loop regions
		cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: After 2nd stage: NR="<<fRegions.size()<<": (nBkg_allRegions,nSig_allRegions,nBoh_allRegions)=("<<nBkg_allRegions<<","<<nSig_allRegions<<","<<nBoh_allRegions<<")"<<endl;
		//=== DEBUG ====
		

		if(nBoh<=0) {
			cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: NR="<<regions.size()<<": No untagged regions left, algorithm end!"<<endl;
			stopAlgo= true;
		}
		if(nTotMergedRegions==0){
			cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: NR="<<regions.size()<<": No regions merged in all stage, mark remaining as signal and end algorithm!"<<endl;
			stopAlgo= true;
			for(unsigned int i=0;i<regions.size();i++) {
				int tag= regions[i]->fTag;
				if(tag==Region::eUntagged) regions[i]->fTag= Region::eSignalTag;
				//if(tag==Region::eUntagged) regions[i]->fTag= Region::eBkgTag;
			}//end loop regions
		}//close if
		
	}//end loop algo
	}//close if
	else{
		cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: No signal regions have been tagged...no segmentation will be performed!"<<endl;
		//Tag all untagged regions as bkg
		for(unsigned int i=0;i<regions.size();i++) {
			int tag= regions[i]->fTag;
			if(tag==Region::eUntagged) regions[i]->fTag= Region::eBkgTag;
		}//end loop regions	
	}//close else

	/*	
	//## Check if untagged regions are occluded (i.e. not surrounded by any bkg region)
	//## Occluded untagged regions are not assigned to the signal!
	SLICUtils::TagOccludedRegions(fInputImg,labels,regions,Region::eUntagged,Region::eBkgTag);
	for(unsigned int i=0;i<regions.size();i++) {
		int tag= regions[i]->fTag;
		bool isOccluded= regions[i]->fIsOccluded;
		if(tag==Region::eUntagged) {
			if(isOccluded) regions[i]->fTag= Region::eBkgTag;
			else regions[i]->fTag= Region::eSignalTag;
		}
	}//end loop regions
	*/

	
	//## Change tag in original segmented regions
	for(unsigned int i=0;i<fRegions.size();i++) {
		fRegions[i]->fTag= Region::eBkgTag;
	}//end loop regions

	for(unsigned int i=0;i<regions.size();i++) {
		int regionTag= regions[i]->fTag;
		int regionId= regions[i]->fId;
		int regionIndex_initialSegm= regionIdMap_initialSegm[regionId]; 
		fRegions[regionIndex_initialSegm]->fTag= regionTag;
		
		std::vector<int> subregionIds= regions[i]->fSubRegionIds;
		for(unsigned int j=0;j<subregionIds.size();j++){
			int subregionId= subregionIds[j];
			int subregionIndex_initialSegm= regionIdMap_initialSegm[subregionId]; 
			fRegions[subregionIndex_initialSegm]->fTag= regionTag;
		}	
	}//end loop regions

	/*
	for(unsigned int i=0;i<regions.size();i++) {
		int regionTag= regions[i]->fTag;
		int regionId= regions[i]->fId;
		int regionIndex_initialSegm= regionIdMap_initialSegm[regionId]; 
		fRegions[regionIndex_initialSegm]->fTag= regionTag;
	}//end loop regions

	for(unsigned int i=0;i<fRegions.size();i++) {
		int regionTag= fRegions[i]->fTag;
		if(regionTag==Region::eUntagged) fRegions[i]->fTag= Region::eSignalTag;
		cout<<"-> Region no. "<<i<<" id="<<fRegions[i]->fId<<", tag="<<fRegions[i]->fTag<<endl;
	}
	*/
	//## Draw results before signal merging
	Img* maxSimilarityMergedImg= SLICUtils::GetSegmentedImage(fInputImg,regions,Region::eSignalTag,true,false);
	if(maxSimilarityMergedImg) maxSimilarityMergedImg->SetNameTitle("maxSimilarityMergedImg","maxSimilarityMergedImg");
	std::vector<TText*> regionTextList= SLICUtils::GetRegionTextLabels(regions,Region::eSignalTag);
		
	TCanvas* maxSimilarityMergedImgPlot= new TCanvas("maxSimilarityMergedImgPlot","maxSimilarityMergedImgPlot");
	maxSimilarityMergedImgPlot->cd();
	maxSimilarityMergedImg->Draw("COLZ");
	/*
	for(unsigned int i=0;i<regionTextList.size();i++) {	
		if(regionTextList[i]) regionTextList[i]->Draw("same");
	}
	*/

	Img* maxSimilarityMergedImg2= SLICUtils::GetSegmentedImage(fInputImg,fRegions,Region::eSignalTag,true,false);
	if(maxSimilarityMergedImg2) maxSimilarityMergedImg->SetNameTitle("maxSimilarityMergedImg2","maxSimilarityMergedImg2");
	std::vector<TText*> regionTextList2= SLICUtils::GetRegionTextLabels(fRegions,Region::eSignalTag);
		
	TCanvas* maxSimilarityMergedImgPlot2= new TCanvas("maxSimilarityMergedImgPlot2","maxSimilarityMergedImgPlot2");
	maxSimilarityMergedImgPlot2->cd();
	maxSimilarityMergedImg2->Draw("COLZ");
	for(unsigned int i=0;i<regionTextList2.size();i++) {	
		if(regionTextList2[i]) regionTextList2[i]->Draw("same");
	}

	//## Hierarchically merge signal regions
	if(regions.size()>2 && nSig>0){
		bool includeSpatialPars= true;
		int edgeModel= 1;
		int mergingStatus= SPHierarchicalMerger(fRegions,fPixelClusterIds,Region::eSignalTag,Region::eSignalTag,includeSpatialPars,edgeModel);
		if(mergingStatus<0){
			cerr<<"SLICSegmenter::MultiStepSPMerger(): ERROR: Merging of signal regions failed!"<<endl;
			return -1;
		}
	}
	cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: Merged signal regions (NR="<<fMergedRegions.size()<<")"<<endl;
	


	//## Add background regions together
	int bkgLabel= -1;
	Region* bkgRegion= new Region;
	bkgRegion->fId= 0;
	std::vector<Region*> signalPlusMergedBkgRegions;
	bool isLabelSet= false;
	int nBkgRegions= 0;
	int nSignalRegions= 0;
	std::map<int,int> relabelMap;

	//for(unsigned int i=0;i<regions.size();i++) {
	for(unsigned int i=0;i<fRegions.size();i++) {
		//Region* thisRegion= regions[i];
		Region* thisRegion= fRegions[i];

		//int tag= regions[i]->fTag;
		//int regionId= regions[i]->fId;
		int tag= thisRegion->fTag;
		int regionId= thisRegion->fId;
		relabelMap.insert( std::pair<int,int>(regionId,regionId) );

		if(tag==Region::eBkgTag) {
			if(!isLabelSet){
				bkgLabel= regionId;
				isLabelSet= true;
			}
			relabelMap[regionId]= bkgLabel;
		
			bkgRegion->AddRegion(thisRegion);
			bkgRegion->fImageSizeX= thisRegion->fImageSizeX;
			bkgRegion->fImageSizeY= thisRegion->fImageSizeY;
			bkgRegion->fImageMinX= thisRegion->fImageMinX;
			bkgRegion->fImageMaxX= thisRegion->fImageMaxX;
			bkgRegion->fImageMinY= thisRegion->fImageMinY;
			bkgRegion->fImageMaxY= thisRegion->fImageMaxY;
			bkgRegion->fImageMinS= thisRegion->fImageMinS;
			bkgRegion->fImageMaxS= thisRegion->fImageMaxS;
			bkgRegion->fImageMinScurv= thisRegion->fImageMinScurv;
			bkgRegion->fImageMaxScurv= thisRegion->fImageMaxScurv;
			bkgRegion->fImageMinSedge= thisRegion->fImageMinSedge;
			bkgRegion->fImageMaxSedge= thisRegion->fImageMaxSedge;
			bkgRegion->fImageRMS= thisRegion->fImageRMS;
			bkgRegion->fTag= Region::eBkgTag;
			bkgRegion->fIsSalient= false;	
			bkgRegion->fIsSignificative= false;
			nBkgRegions++;
		}//close if bkg regions
		else if(tag==Region::eSignalTag){
			thisRegion->fIsSignificative= true;
			thisRegion->fIsSalient= true;
			fMergedRegions.push_back(thisRegion);	
			signalPlusMergedBkgRegions.push_back(thisRegion);
			nSignalRegions++;
		}
		else{
			cerr<<"SLICSegmenter::MultiStepSPMerger(): WARN: Untagged region present (CHECK!!!!)"<<endl;
		}
	}//end loop regions
	bkgRegion->fId= bkgLabel;
	fMergedRegions.push_back(bkgRegion);
	signalPlusMergedBkgRegions.push_back(bkgRegion);

	
	for(unsigned int i=0;i<fMergedRegions.size();i++){
		fMergedRegions[i]->ComputeParameters(false,false,true);
	}//end loop merged regions
	
	/*
	for(unsigned int i=0;i<labels.size();i++) {
		for(unsigned int j=0;j<labels[i].size();j++) {
			int oldLabel= labels[i][j];
			int newLabel= relabelMap[oldLabel];
			labels[i][j]= newLabel;
			fMergedPixelClusterIds[i][j]= labels[i][j];
		}
	}
	*/
	for(unsigned int i=0;i<fPixelClusterIds.size();i++) {
		for(unsigned int j=0;j<fPixelClusterIds[i].size();j++) {
			int oldLabel= fPixelClusterIds[i][j];
			int newLabel= relabelMap[oldLabel];
			labels[i][j]= newLabel;
			fMergedPixelClusterIds[i][j]= labels[i][j];
		}
	}

	cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: Get segmented image before merging..."<<endl;
	bool drawOnlySignificant= true;
	bool drawOnlySalient= false;
	bool normalize= false; 
	bool binarize= true;
	Img* segmentedSignalPlusMergedBkgImg= SLICUtils::GetSegmentedImage(fInputImg,signalPlusMergedBkgRegions,drawOnlySignificant,drawOnlySalient,normalize,binarize);
	segmentedSignalPlusMergedBkgImg->SetNameTitle("segmentedSignalPlusMergedBkgImg","segmentedSignalPlusMergedBkgImg");

	TCanvas* SignalPlusMergedBkgPlot= new TCanvas("SignalPlusMergedBkg","SignalPlusMergedBkg");
	SignalPlusMergedBkgPlot->cd();
	segmentedSignalPlusMergedBkgImg->Draw("COLZ");

	/*
	//if(fMergedRegions.size()>2){
	if(regions.size()>2){
		//int mergingStatus= SPHierarchicalMerger(fMergedRegions,fMergedPixelClusterIds,false,Region::eSignalTag,Region::eSignalTag);
		int mergingStatus= SPHierarchicalMerger(regions,labels,false,Region::eSignalTag,Region::eSignalTag);

		if(mergingStatus<0){
			cerr<<"SLICSegmenter::MultiStepSPMerger(): ERROR: Merging of signal regions failed!"<<endl;
			return -1;
		}
	}

	fMergedRegions.push_back(bkgRegion);
	for(unsigned int i=0;i<regions.size();i++){
		int tag= regions[i]->fTag;
		if(tag==Region::eSignalTag){
			regions[i]->fIsSignificative= true;
			regions[i]->fIsSalient= true;
			fMergedRegions.push_back(regions[i]);
		}
	}//end loop regions	
	
	for(unsigned int i=0;i<labels.size();i++) {
		for(unsigned int j=0;j<labels[i].size();j++) {
			fMergedPixelClusterIds[i][j]= labels[i][j];
		}
	}
	cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: Merged signal regions (NR="<<fMergedRegions.size()<<")"<<endl;
	
	*/
	

	//## Drawing results
	cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: Get segmented image..."<<endl;
	Img* segmentedImg= SLICUtils::GetSegmentedImage(fInputImg,fMergedRegions,false,false,true);
	if(segmentedImg) segmentedImg->SetNameTitle("finalSegmentedImg","finalSegmentedImg");

	cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: Drawing results..."<<endl;
	TCanvas* FinalPlot= new TCanvas("FinalPlot","FinalPlot");
	FinalPlot->cd();
	if(segmentedImg) segmentedImg->Draw("COLZ");

	cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: Drawing saliency results..."<<endl;
	TCanvas* SaliencyPlot= new TCanvas("SaliencyPlot","SaliencyPlot");
	SaliencyPlot->cd();
	if(fSaliencyImg) fSaliencyImg->Draw("COLZ");
	
	return 0;

}//close MultiStepSPMerger()


int SLICSegmenter::SPMaxSimilarityMerger(std::vector<Region*>& regions,std::vector< std::vector<long int> >& labels,int mergerTag,int mergedTag){

	//## Check regions
	int nRegions= (int)regions.size();
	if(nRegions<=0) {
		cerr<<"SLICSegmenter::SPMaxSimilarityMerger(): WARN: No regions available, nothing to be merged!"<<endl;
		return 0;
	}
	
	//## Create the mapping of regionId and vector index
	cout<<"SLICSegmenter::SPMaxSimilarityMerger(): INFO: Create the mapping of regionId and region list index..."<<endl;
	std::map<int,int> regionIdMap;
	std::map<int,int> regionIdMap_top;
	std::map<int,int> regionIdMap_initialSegm;	
	std::vector< std::vector<int> > mergedRegionList;
	
	for(int k=0;k<nRegions;k++){
		int regionId= regions[k]->fId;
		regionIdMap.insert( std::pair<int,int>(regionId,k) );
		regionIdMap_top.insert( std::pair<int,int>(regionId,k) );
		mergedRegionList.push_back( std::vector<int>() );
	}//end loop regions

	for(int k=0;k<fRegions.size();k++){
		int regionId= fRegions[k]->fId;
		regionIdMap_initialSegm.insert( std::pair<int,int>(regionId,k) );
	}//end loop regions
	
	//## Compute region contour info (neighbors, ...)
	cout<<"SLICSegmenter::SPMaxSimilarityMerger(): INFO: Finding region neighbors (NR="<<regions.size()<<") ..."<<endl;
	SLICUtils::SLICContourData* contourData= SLICUtils::ComputeBoundaryContours(fInputImg,labels,regions);
	SLICUtils::SLICConnectedRegions connectedRegionIds= contourData->connectedRegionIds;


	
	//## Decomment if you want to use 2nd neighbors in max similarity merging
	//Fill list of 2-nd neighbors	
	if(fUse2ndNeighborsInSPMerging){
		std::vector< std::vector<int> > neighborIndexList_2nd;
		for(int i=0; i<regions.size(); i++) {
			int regionId= regions[i]->fId;
			neighborIndexList_2nd.push_back( std::vector<int>() );		
		
			for(unsigned int j=0;j<connectedRegionIds[i].size();j++){//loop over 1st neighbors
				int neighborIndex= connectedRegionIds[i][j];
				for(unsigned int t=0;t<connectedRegionIds[neighborIndex].size();t++){//loop over 2nd neighbors
					int neighborIndex_2nd= connectedRegionIds[neighborIndex][t];
				
					std::vector<int>::iterator it= std::find(connectedRegionIds[i].begin(),connectedRegionIds[i].end(),neighborIndex_2nd);
					std::vector<int>::iterator it2= std::find(neighborIndexList_2nd[i].begin(),neighborIndexList_2nd[i].end(),neighborIndex_2nd);
		
					if( it==connectedRegionIds[i].end() && (it2==neighborIndexList_2nd[i].end() || neighborIndexList_2nd[i].size()==0) ){
						neighborIndexList_2nd[i].push_back(neighborIndex_2nd);
					}
				}//end loop 2nd-neighbors
			}//end loop 1st-neighbors

			connectedRegionIds[i].insert(connectedRegionIds[i].end(),neighborIndexList_2nd[i].begin(),neighborIndexList_2nd[i].end());
		}//end loop regions
	}//close if use 2nd neighbors
	

	//## Compute similarity matrix	
	cout<<"SLICSegmenter::SPMaxSimilarityMerger(): INFO: Compute region similarity matrix..."<<endl;
	SLICUtils::SLICSimilarityData* similarityData= SLICUtils::ComputeRegionSimilarity(fEdgeFilterImg,contourData,regions,fSPMergingRegularization,fSPMergingIncludeSpatialPars,mergedTag);
	TMatrixD* AdjacencyMatrix= similarityData->AdjacencyMatrix;
	TMatrixD* DissimilarityMatrix= similarityData->DissimilarityMatrix; 
	TMatrixD* AbsDissimilarityMatrix= similarityData->AbsDissimilarityMatrix; 
	std::vector< std::vector<int> > DissSortIndexes= similarityData->DissimilaritySortIndexMatrix;
		

	//## Compute page rank of segments and sort
	cout<<"SLICSegmenter::SPMaxSimilarityMerger(): INFO: Compute page rank ..."<<endl;	
	std::vector<double> ranks= Utils::ComputePageRank(AdjacencyMatrix->T());//pass transpose of adjacency matrix
	if(ranks.size()<=0){
		cerr<<"SLICSegmenter::SPMaxSimilarityMerger(): WARN: PageRank failed, cannot perform region merging!"<<endl;
		return -1;
	}
	std::vector<size_t> sort_index;//sorting index
	std::vector<double> ranks_sorted;
	Utils::sort_descending(ranks,ranks_sorted,sort_index);
		
	//## Loop over sorted ranks and select regions to be merged
	cout<<"SLICSegmenter::SPMaxSimilarityMerger(): INFO: Start merging loop ..."<<endl;	
	int nMergedRegions= 0;
	int maxRegionsToMerge= std::round(regions.size()*fSPMergingRatio);
	
	std::vector< std::vector<int> > regionsToBeMerged;
	for(unsigned int i=0;i<regions.size();i++) regionsToBeMerged.push_back( std::vector<int>() );
	std::vector<int> regionsToBeDeleted;
	std::vector<int> regionsIdToBeDeleted;

	for(unsigned int i=0;i<sort_index.size();i++){
		size_t index= sort_index[i];//region index in regions list
		double thisRank= ranks[index];//region rank
		int regionId= regions[index]->fId;
		int regionTag= regions[index]->fTag;
		int mapIndex= regionIdMap[regionId];
		bool isSignificative= regions[index]->fIsSignificative;
		if(regionTag!=mergerTag && mergerTag!=-1) {
			cout<<"SLICSegmenter::SPMaxSimilarityMerger(): INFO: Skip region no. "<<i<<" (id="<<regionId<<", tag="<<regionTag<<")..."<<endl;
			continue;
		}
			
		//Check if this seed was not already merged by a previous (best ranked) region
		std::vector<int>::iterator seedfinderIt= std::find(regionsToBeDeleted.begin(),regionsToBeDeleted.end(),index);
		if(seedfinderIt!=regionsToBeDeleted.end()){
			cout<<"SLICSegmenter::SPMaxSimilarityMerger(): INFO: Seed ranked region (id="<<regionId<<") was already selected for merging in a previous node, skip this merging!"<<endl;
			continue;
		}

	
		//Loop over untagged neighbors and find best merging
		int nGoodNeighbors= 0;

		for(unsigned int j=0;j<connectedRegionIds[index].size();j++){//loop over 1st-neighbors
			int neighborIndex= connectedRegionIds[index][j];
			int neighborId= regions[neighborIndex]->fId;
			int neighborTag= regions[neighborIndex]->fTag;
			//if(neighborTag!=Region::eUntagged) {//look at untagged neighbors only
			//	continue;
			//}

			if(mergerTag==Region::eBkgTag && neighborTag==mergerTag){//for bkg-untagged merging look at neighbors with tag different than bkg
				continue;
			}
			if(mergerTag==Region::eUntagged && neighborTag!=mergedTag){//for untagged-untagged merging look at neighbors with 'untagged' tag
				continue;
			}
			nGoodNeighbors++;
			
			double Delta_ij= (*DissimilarityMatrix)(index,neighborIndex);

			//Loop over neighbors of this neighbors (2nd neighbors)
			int closerIndex= -1;
			double closerDiss= 1.e+99;
			for(unsigned int k=0;k<connectedRegionIds[neighborIndex].size();k++){
				int neighborIndex_2nd= connectedRegionIds[neighborIndex][k];
				int neighborId_2nd= regions[neighborIndex_2nd]->fId;
						
				double Delta_jk= (*DissimilarityMatrix)(neighborIndex,neighborIndex_2nd);
				if(Delta_jk<closerDiss && Delta_jk>0){	
					closerDiss= Delta_jk;	
					closerIndex= neighborIndex_2nd;
				}
			}//end loop 2nd neighbors

			//If diss Dij is the maximum among Djk(s) merge ij!
			if(closerIndex!=index) continue;


			//Check if this closer region was not already selected to be merged to another node previously	
			std::vector<int>::iterator finderIt= std::find(regionsToBeDeleted.begin(),regionsToBeDeleted.end(),neighborIndex);				
			if(finderIt!=regionsToBeDeleted.end()){	
				cout<<"SLICSegmenter::SPMaxSimilarityMerger(): INFO: Closer neighbor no. "<<neighborIndex<<" (id="<<neighborId<<") was already selected for merging in a previous node, skip this merging!"<<endl;
				continue;
			}

			//Check if this neighbor was not previously selected as a merger
			if(regionsToBeMerged[neighborIndex].size()>0){
				cout<<"SLICSegmenter::SPMaxSimilarityMerger(): INFO: Closer neighbor no. "<<neighborIndex<<" (id="<<neighborId<<") was previously selected as a merger, skip this merging!"<<endl;
				continue;
			}

			//Merging selected!
			cout<<"SLICSegmenter::SPMaxSimilarityMerger(): INFO: New merging ("<<index<<"-"<<neighborIndex<<"), Delta_ij="<<Delta_ij<<" closerDiss="<<closerDiss<<endl;
			int regionIndex_A= regionIdMap_top[regionId];
			int regionIndex_B= regionIdMap_top[neighborId];		

			regionsToBeMerged[index].push_back(neighborIndex);
			regionsToBeDeleted.push_back(neighborIndex);
			regionsIdToBeDeleted.push_back(neighborId);		
			mergedRegionList[regionIndex_A].push_back(regionIndex_B);
			nMergedRegions++;

		}//end loop neighbors
				
	}//end loop ranks

	//## Delete contour data
	if(contourData){
		delete contourData;
		contourData= 0;
	}
	if(similarityData){
		delete similarityData;
		similarityData= 0;
	}
		
	//## Merge regions and update region and label list
	if(nMergedRegions>0){
		cout<<"SLICSegmenter::SPMaxSimilarityMerger(): INFO: Merge the selected regions..."<<endl;
		std::map<int,int> newLabelMap;

		for(unsigned int k=0;k<regions.size();k++){
			int regionId= regions[k]->fId;	
			int regionTag= regions[k]->fTag;
			newLabelMap.insert( std::pair<int,int>(regionId,regionId) );
			for(unsigned int j=0;j<regionsToBeMerged[k].size();j++){
				int mergedRegionIndex= regionsToBeMerged[k][j];
				int mergedRegionId= regions[mergedRegionIndex]->fId;
				int mergedRegionTag= regions[mergedRegionIndex]->fTag;
				newLabelMap[mergedRegionId]= regionId;

				//Change tag in original segmentation
				int mergedRegionIndex_originalSegm= regionIdMap_initialSegm[mergedRegionId];
				int mergedRegionId_originalSegm= fRegions[mergedRegionIndex_originalSegm]->fId;
				//fRegions[mergedRegionIndex_originalSegm]->fTag= regions[k]->fTag;

				cout<<"SLICSegmenter::SPMaxSimilarityMerger(): INFO: Region no. "<<k<<" (id="<<regionId<<",tag="<<regionTag<<") : merging region id="<<mergedRegionId<<", tag="<<mergedRegionTag<<" (index="<<mergedRegionIndex<<",orig index="<<mergedRegionIndex_originalSegm<<", orig id="<<mergedRegionId_originalSegm<<")"<<endl;		
				regions[k]->AddRegion(regions[mergedRegionIndex]);
				regions[k]->AddSubRegionId(mergedRegionId);
			
			}//end loop regions to be merged

			//Change tag in original segmentation
			int regionIndex_originalSegm= regionIdMap_initialSegm[regionId];
			//fRegions[regionIndex_originalSegm]->fTag= regionTag;
			
		}//end loop regions
	
		//## Delete aggregated region from region list and index map
		cout<<"SLICSegmenter::SPMaxSimilarityMerger(): INFO: Deleting regions aggregated in this step from the main list..."<<endl;
		Utils::DeleteItems(regions, regionsToBeDeleted);
		for(size_t k=0;k<regionsIdToBeDeleted.size();k++) regionIdMap.erase(regionsIdToBeDeleted[k]);

		//## Update map and recompute parameters & contours (will be used by nearest neighbors search)
		cout<<"SLICSegmenter::SPMaxSimilarityMerger(): INFO: Updating region parameters & contours..."<<endl;
		for(unsigned int k=0;k<regions.size();k++){
			regions[k]->ComputeParameters(false,false,true);
			int regionId= regions[k]->fId;
			regionIdMap[regionId]= k;
		}//end loop regions

		//## Update pixel labels
		cout<<"SLICSegmenter::MultiStepSPMerger(): INFO: Updating pixel labels..."<<endl;	
		for(unsigned int i=0;i<labels.size();i++) {
			for(unsigned int j=0;j<labels[i].size();j++) {
				int oldLabel= labels[i][j];
				int newLabel= newLabelMap[oldLabel];
				labels[i][j]= newLabel;
			}
		}

	}//close if merge regions
		
	cout<<"SLICSegmenter::SPMaxSimilarityMerger(): INFO: "<<nMergedRegions<<"/"<<regions.size()<<" regions merged at this stage..."<<endl;

	return 0;

}//close SPMaxSimilarityMerger()



int SLICSegmenter::SPHierarchicalMerger(std::vector<Region*>& inputRegions,std::vector< std::vector<long int> >& inputLabels,int mergerTag,int mergedTag,bool includeSpatialPars,int edgeModel){

	//## Check regions
	int nRegions= (int)inputRegions.size();
	if(nRegions<=0) {
		cerr<<"SLICSegmenter::SPHierarchicalMerger(): WARN: No regions available, nothing to be merged!"<<endl;
		return 0;
	}
	double Nx= fInputImg->GetNbinsX();
	double Ny= fInputImg->GetNbinsY();
	double N= Nx*Ny;

	//## Init region list to be passed to algorithm
	cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Copying initial superpixel partition to tmp partition..."<<endl;
	std::vector<Region*> regions;
	regions.assign(inputRegions.begin(),inputRegions.end());
	
	//## Create the mapping of regionId and vector index
	cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Create the mapping of regionId and region list index..."<<endl;
	std::map<int,int> regionIdMap;
	std::map<int,int> regionIdMap_top;	
	std::vector<std::vector<int>> mergedRegionList;
	
	for(int k=0;k<nRegions;k++){
		int regionId= regions[k]->fId;
		regionIdMap.insert( std::pair<int,int>(regionId,k) );
		regionIdMap_top.insert( std::pair<int,int>(regionId,k) );
		mergedRegionList.push_back( std::vector<int>() );
	}//end loop regions
	
	//Init pixel labels: copy initial pixel labels to tmp list
	std::vector< std::vector<long int> > labels;
	for(unsigned int i=0;i<inputLabels.size();i++) {
		labels.push_back( std::vector<long int>() );
		labels[i].assign(inputLabels[i].begin(),inputLabels[i].end());
	}

	//## Run a hierarchical clustering till all segments are merged in one
	int hierarchyLevel= 0;
	double AbsDissMedian= 0;
	double AbsDissMedianRMS= 0;
	double AbsDissMin= 0;
	double AbsDissMax= 0;
	double AbsDissMedian0= 0;
	double AbsDissMedianRMS0= 0;
	int nMergedRegionsInHierarchyLevel= 0;
	fSPMergingInfo.clear();	
	cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Starting hierarchical merging..."<<endl;	

	while(regions.size()>fMinMergedSP){
		nMergedRegionsInHierarchyLevel= 0;
		
		//Merge info for this hierarchy
		SPMergingInfo spMergeInfo;
		spMergeInfo.levelId= hierarchyLevel;
		spMergeInfo.NR= (int)regions.size();
		spMergeInfo.MSE= 0;
		spMergeInfo.DissMin= 0;
		spMergeInfo.DissMax= 0;
		spMergeInfo.DissMean= 0;
		spMergeInfo.DissRMS= 0;
		
		//## Compute region contour info (neighbors, ...)
		cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Finding region neighbors at hierarchy level "<<hierarchyLevel<<", NR="<<regions.size()<<" ..."<<endl;
		SLICUtils::SLICContourData* contourData= SLICUtils::ComputeBoundaryContours(fInputImg,labels,regions);
		SLICUtils::SLICConnectedRegions connectedRegionIds= contourData->connectedRegionIds;
		std::vector< std::vector<int> > neighborIndexList_1st;
		std::vector< std::vector<int> > neighborIndexList_2nd;

		for(int i=0; i<regions.size(); i++) {
			int regionId= regions[i]->fId;
			//Fill list of 1st and 2-nd neighbors	
			neighborIndexList_1st.push_back( std::vector<int>() );
			neighborIndexList_2nd.push_back( std::vector<int>() );
		
			for(unsigned int j=0;j<connectedRegionIds[i].size();j++){//loop over 1st neighbors
				int neighborIndex= connectedRegionIds[i][j];
				int neighborTag= regions[neighborIndex]->fTag;
				int neighborId= regions[neighborIndex]->fId;
				if(mergedTag!=-1 && neighborTag==mergedTag){
					neighborIndexList_1st[i].push_back(neighborIndex);
				}

				for(unsigned int t=0;t<connectedRegionIds[neighborIndex].size();t++){//loop over 2nd neighbors
					int neighborIndex_2nd= connectedRegionIds[neighborIndex][t];
					int neighborTag_2nd= regions[neighborIndex_2nd]->fTag;
					int neighborId_2nd= regions[neighborIndex_2nd]->fId;
					if(mergedTag!=-1 && neighborTag_2nd==mergedTag && neighborId_2nd!=regionId){
						std::vector<int>::iterator it= std::find(connectedRegionIds[i].begin(),connectedRegionIds[i].end(),neighborIndex_2nd);
						std::vector<int>::iterator it2= std::find(neighborIndexList_2nd[i].begin(),neighborIndexList_2nd[i].end(),neighborIndex_2nd);
		
						if( it==connectedRegionIds[i].end() && (it2==neighborIndexList_2nd[i].end() || neighborIndexList_2nd[i].size()==0) ){
							neighborIndexList_2nd[i].push_back(neighborIndex_2nd);
						}	
					}//close if
				}//end loop 2nd-neighbors
			}//end loop 1st-neighbors
		}//end loop regions

		//## Compute similarity matrix	
		cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Compute region similarity at hierarchy level "<<hierarchyLevel<<" ..."<<endl;

		Img* edgeImg= fEdgeFilterImg;
		if(edgeModel==1) edgeImg= fKirschEdgeFilterImg;
		SLICUtils::SLICSimilarityData* similarityData= SLICUtils::ComputeRegionSimilarity(edgeImg,contourData,regions,fSPMergingRegularization,includeSpatialPars,mergedTag);
		TMatrixD* AdjacencyMatrix= similarityData->AdjacencyMatrix;
		TMatrixD* DissimilarityMatrix= similarityData->DissimilarityMatrix; 
		TMatrixD* AbsDissimilarityMatrix= similarityData->AbsDissimilarityMatrix; 
		TMatrixD* NeighborMatrix= similarityData->NeighborMatrix; 
		std::vector< std::vector<int> > DissSortIndexes= similarityData->DissimilaritySortIndexMatrix;
		AbsDissMedian= similarityData->Dmedian;
 		AbsDissMedianRMS= similarityData->Dmedianrms;
 		AbsDissMin= similarityData->Dmin;
		AbsDissMax= similarityData->Dmax;

		if(hierarchyLevel==0) {
			AbsDissMedian0= AbsDissMedian;
			AbsDissMedianRMS0= AbsDissMedianRMS;
		}
		
		double MSE= 0.;
		for(unsigned int k=0;k<regions.size();k++){		
			double M2= regions[k]->fM2;
			MSE+= M2;
		}//end loop regions
		MSE/= N;
		spMergeInfo.MSE= MSE;

		//## Compute page rank of segments and sort
		cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Compute page rank at hierarchy level "<<hierarchyLevel<<" ..."<<endl;	
		std::vector<double> ranks= Utils::ComputePageRank(AdjacencyMatrix->T());//pass transpose of adjacency matrix
		if(ranks.size()<=0){
			cerr<<"SLICSegmenter::SPHierarchicalMerger(): WARN: PageRank failed, cannot perform region merging!"<<endl;
			return -1;
		}
		std::vector<size_t> sort_index;//sorting index
		std::vector<double> ranks_sorted;
		Utils::sort_descending(ranks,ranks_sorted,sort_index);
		
		
		//## Count max number of mergeable regions
		int nMaxMergeableRegions= 0;
		for(unsigned int i=0;i<sort_index.size();i++){
			size_t index= sort_index[i];//region index in regions list
			int regionTag= regions[index]->fTag;
			if(regionTag!=mergerTag && mergerTag!=-1) continue;
			nMaxMergeableRegions++;
		}//end loop regions
		//int maxRegionsToMerge= std::round(regions.size()*fSPMergingRatio);
		int maxRegionsToMerge= std::round(nMaxMergeableRegions*fSPMergingRatio);

		//## Loop over sorted ranks and select regions to be merged
		int nMergedRegions= 0;
		int nMergeableRegions= 0;
		std::vector< std::vector<int> > regionsToBeMerged;
		for(unsigned int i=0;i<regions.size();i++) regionsToBeMerged.push_back( std::vector<int>() );
		std::vector<int> regionsToBeDeleted;
		std::vector<int> regionsIdToBeDeleted;

		for(unsigned int i=0;i<sort_index.size();i++){
			size_t index= sort_index[i];//region index in regions list
			double thisRank= ranks[index];//region rank
			int regionId= regions[index]->fId;
			int regionTag= regions[index]->fTag;
			int mapIndex= regionIdMap[regionId];
			bool isSignificative= regions[index]->fIsSignificative;
			if(regionTag!=mergerTag && mergerTag!=-1) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Skip ranked seed region (id="<<regionId<<",tag="<<regionTag<<") as tag is different from mergerTag!"<<endl;
				continue;
			}

			//Stop merging above nmerge threshold						
			if(nMergeableRegions+1>maxRegionsToMerge) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Maximum number of regions mergeable reached for this hierarchy level ("<<nMergeableRegions+1<<">"<<maxRegionsToMerge<<")"<<endl;
				break;
			}
			
			//Chech if the seed region has any neighbors (according to the selected merged tag)
			int NN_1st= (int)neighborIndexList_1st[index].size();
			int NN_2nd= (int)neighborIndexList_2nd[index].size();	
			if(NN_1st<=0 && NN_2nd<=0){
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Skip ranked seed region (id="<<regionId<<",tag="<<regionTag<<") as no neighbors are available!"<<endl;
				continue;
			}

			//Get closest neighbor
			if(DissSortIndexes[index].size()<=0) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Skip seed ranked region (id="<<regionId<<") as no neighbors are available!"<<endl;
				continue;
			}
			int closerNeighborIndex= DissSortIndexes[index][0];
			int closerNeighborId= regions[closerNeighborIndex]->fId;
			int closerNeighborTag= regions[closerNeighborIndex]->fTag;
			bool closerNeighborIsSignificative= regions[closerNeighborIndex]->fIsSignificative;
			
			if(closerNeighborId==regionId) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Skip neighbor (id="<<closerNeighborId<<") as equal to seed region!"<<endl;
				continue;//skip same region
			}
			if(closerNeighborTag!=mergedTag && mergedTag!=-1) {	
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Skip ranked seed region (id="<<regionId<<",tag="<<regionTag<<") as neighbor (id="<<closerNeighborId<<",tag="<<closerNeighborTag<<") tag is different from mergedTag!"<<endl;
				continue;
			}
			nMergeableRegions++;
			
			//Check if this seed was not already merged by a previous (best ranked) region
			std::vector<int>::iterator seedfinderIt= std::find(regionsIdToBeDeleted.begin(),regionsIdToBeDeleted.end(),regionId);
			if(seedfinderIt!=regionsIdToBeDeleted.end()){
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Seed ranked region (id="<<regionId<<") was already selected for merging in a previous node, skip this merging!"<<endl;
				continue;
			}


			double Delta_ij= (*DissimilarityMatrix)(index,closerNeighborIndex);//include edgeness
			double Delta_ji= (*DissimilarityMatrix)(closerNeighborIndex,index);//include edgeness
			double AbsDelta_ij= (*AbsDissimilarityMatrix)(index,closerNeighborIndex);//without edgeness
			double AbsDelta_ji= (*AbsDissimilarityMatrix)(closerNeighborIndex,index);//without edgeness
			int closerNeighborness= (*NeighborMatrix)(index,closerNeighborIndex);

			
			cout<<"Region no. "<<index<<" (id="<<regionId<<") bestNN="<<closerNeighborId<<" 1stNN(";
			for(unsigned int l=0;l<neighborIndexList_1st[index].size();l++){
				int nnIndex= neighborIndexList_1st[index][l];
				int nnId= regions[nnIndex]->fId;
				cout<<nnId<<",";
			}
			cout<<") 2ndNN(";
			for(unsigned int l=0;l<neighborIndexList_2nd[index].size();l++){
				int nnIndex= neighborIndexList_2nd[index][l];
				int nnId= regions[nnIndex]->fId;
				cout<<nnId<<",";
			}
			cout<<")"<<endl;
			cout<<") SortNN(";
			for(unsigned int l=0;l<DissSortIndexes[index].size();l++){
				int nnIndex= DissSortIndexes[index][l];
				int nnId= regions[nnIndex]->fId;
				cout<<nnId<<",";
			}
			cout<<")"<<endl;
			

			//Get neighbors of the closer region
			std::vector<int> closerNeighborNNIndexes;
			closerNeighborNNIndexes.assign(
				DissSortIndexes[closerNeighborIndex].begin(),
				DissSortIndexes[closerNeighborIndex].begin() + ceil(DissSortIndexes[closerNeighborIndex].size()/2.)
			);

			//Find current region index among best neighbors of closer
			std::vector<int>::iterator it= std::find(closerNeighborNNIndexes.begin(),closerNeighborNNIndexes.end(),index);
			if(it==closerNeighborNNIndexes.end() ) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Neighbor region id="<<closerNeighborId<<" rejected for merging with region "<<regionId<<" since this is not among the closest neighbors of selected neighbor!"<<endl;
				continue;
			}

			
			//## Check similarity difference
			//## 1st neighbors are merged (flood-fill)
			//## 2nd neighbors are merged if their similarities are not too different
			if(closerNeighborness==1 && (Delta_ji>Delta_ij*fSPMergingMaxDissRatio) ) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Neighbor region id="<<closerNeighborId<<" (1st neighbor) rejected for merging with region "<<regionId<<"(Delta_ji="<<Delta_ji<<">"<<fSPMergingMaxDissRatio<<"*Delta_ij="<<Delta_ij<<endl;
				continue;
			}
			else if(closerNeighborness==2 && (Delta_ji>Delta_ij*fSPMergingMaxDissRatio_2ndNeighbor) ) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Neighbor region id="<<closerNeighborId<<" (2nd neighbor) rejected for merging with region "<<regionId<<"(Delta_ji="<<Delta_ji<<">"<<fSPMergingMaxDissRatio_2ndNeighbor<<"*Delta_ij="<<Delta_ij<<endl;
				continue;
			}
			

			/*
			if(closerNeighborness==1 && (AbsDelta_ij>fSPMergingMaxDissRatio) ) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Neighbor region id="<<closerNeighborId<<" (1st neighbor) rejected for merging with region "<<regionId<<"(AbsDelta_ij="<<AbsDelta_ij<<">"<<fSPMergingMaxDissRatio<<")"<<endl;
				continue;
			}
			else if(closerNeighborness==2 && (AbsDelta_ij>fSPMergingMaxDissRatio_2ndNeighbor) ) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Neighbor region id="<<closerNeighborId<<" (2nd neighbor) rejected for merging with region "<<regionId<<"(AbsDelta_ij="<<AbsDelta_ij<<">"<<fSPMergingMaxDissRatio_2ndNeighbor<<")"<<endl;
				continue;
			}
			*/


			//Check if this closer region was not already selected as a seed merger previously
			if( regionsToBeMerged[closerNeighborIndex].size()>0 ){	
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Closer neighbor (id="<<closerNeighborId<<") selected for merging was before selected as a primary merger, skip this merging!"<<endl;
				continue;
			}

			//Check if this closer region was not already selected to be merged to another node previously	
			std::vector<int>::iterator finderIt= std::find(regionsIdToBeDeleted.begin(),regionsIdToBeDeleted.end(),closerNeighborId);
			if(finderIt!=regionsIdToBeDeleted.end()){
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Closer neighbor (id="<<closerNeighborId<<") was already selected for merging in a previous node, skip this merging!"<<endl;
				continue;
			}

			
			//## Apply dissimilarity threshold
			//double DissThreshold= fabs(fAbsDissMedian + fSPMergingDistThreshold*fAbsDissMedianRMS);
			double DissThreshold= fSPMergingDistThreshold*AbsDissMedian0;
			//if(fSPMergingUseAdaptingDistThreshold) DissThreshold+= hierarchyLevel*(AbsDissMax-AbsDissMin)/fSPMergingAdaptingDistThresholdScale;
			
			if(closerNeighborness==2 && AbsDelta_ij>DissThreshold){
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Region id="<<regionId<<": closer neighbor (id="<<closerNeighborId<<") cannot be merged as dissimilarity is too large (AbsDiss="<<AbsDelta_ij<<", Diss="<<Delta_ij<<">"<<DissThreshold<<")"<<endl;
				continue;
			}
			//cout<<"Region id="<<regionId<<": closer neighbor (id="<<closerNeighborId<<") selected for merging (Delta_ij="<<Delta_ij<<", AbsDelta_ij="<<AbsDelta_ij<<", DissThreshold="<<DissThreshold<<")"<<endl;
			
			cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Region (id="<<regionId<<") merging neighbor (id="<<closerNeighborId<<"): closerNeighborness="<<closerNeighborness<<", Delta_ij="<<Delta_ij<<", Delta_ji="<<Delta_ji<<" ratio="<<Delta_ji/Delta_ij<<endl;

			
			int regionIndex_A= regionIdMap_top[regionId];
			int regionIndex_B= regionIdMap_top[closerNeighborId];
			
			regionsToBeMerged[index].push_back(closerNeighborIndex);
			regionsToBeDeleted.push_back(closerNeighborIndex);
			regionsIdToBeDeleted.push_back(closerNeighborId);		
			mergedRegionList[regionIndex_A].push_back(regionIndex_B);
	
			nMergedRegions++;

			/*
			//Stop merging above nmerge threshold			
			int maxRegionsToMerge= std::round(regions.size()*fSPMergingRatio);
			if(nMergedRegions>=maxRegionsToMerge) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Maximum number of regions mergeable reached for this hierarchy level!"<<endl;
				break;
			}
			*/
		}//end loop ranks

		//## Add merging info
		spMergeInfo.DissMin= AbsDissMin;
		spMergeInfo.DissMax= AbsDissMax;
		spMergeInfo.DissMedian= AbsDissMedian;
		spMergeInfo.DissMedianRMS= AbsDissMedianRMS;
		spMergeInfo.DissMedian0= AbsDissMedian0;
		spMergeInfo.DissMedianRMS0= AbsDissMedianRMS0;
		fSPMergingInfo.push_back(spMergeInfo);
		//cout<<"SLICSegmenter::HierarchicalMergeRegions(): INFO: Merge info @ level="<<spMergeInfo.levelId<<": NR="<<spMergeInfo.NR<<", Diss min/max="<<spMergeInfo.DissMin<<"/"<<spMergeInfo.DissMax<<", Diss0 median/medianrms="<<spMergeInfo.DissMedian0<<"/"<<spMergeInfo.DissMedianRMS0<<", MSE="<<MSE<<endl;

		//## Delete contour data
		if(contourData){
			delete contourData;
			contourData= 0;
		}
		if(similarityData){
			delete similarityData;
			similarityData= 0;
		}
		
	
		//## Merge regions
		cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Merge the selected regions..."<<endl;
		std::map<int,int> newLabelMap;

		for(unsigned int k=0;k<regions.size();k++){
			int regionId= regions[k]->fId;
			newLabelMap.insert( std::pair<int,int>(regionId,regionId) );
			for(unsigned int j=0;j<regionsToBeMerged[k].size();j++){
				int mergedRegionIndex= regionsToBeMerged[k][j];
				int mergedRegionId= regions[mergedRegionIndex]->fId;
				newLabelMap[mergedRegionId]= regionId;
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Region no. "<<k<<" (id="<<regionId<<") : merging region id="<<mergedRegionId<<" (index="<<mergedRegionIndex<<")"<<endl;		
				regions[k]->AddRegion(regions[mergedRegionIndex]);
			}//end loop regions to be merged
		}//end loop regions
	
		//## Delete aggregated region from region list and index map
		cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Deleting regions aggregated in this step from the main list..."<<endl;
		Utils::DeleteItems(regions, regionsToBeDeleted);
		for(size_t k=0;k<regionsIdToBeDeleted.size();k++) regionIdMap.erase(regionsIdToBeDeleted[k]);

		//## Update map and recompute parameters & contours (will be used by nearest neighbors search)
		cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Updating region parameters & contours..."<<endl;
		for(unsigned int k=0;k<regions.size();k++){
			regions[k]->ComputeParameters(false,false,true);
			int regionId= regions[k]->fId;
			regionIdMap[regionId]= k;
		}//end loop regions
	
		//## Update pixel labels
		cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Updating pixel labels..."<<endl;	
		for(unsigned int i=0;i<labels.size();i++) {
			for(unsigned int j=0;j<labels[i].size();j++) {
				int oldLabel= labels[i][j];
				int newLabel= newLabelMap[oldLabel];
				labels[i][j]= newLabel;
			}
		}
	
		nMergedRegionsInHierarchyLevel= nMergedRegions;
		hierarchyLevel++;

		if(nMergedRegionsInHierarchyLevel==0){
			cerr<<"SLICSegmenter::SPHierarchicalMerger(): WARN: No regions merged in this stage, exit to avoid stuck!"<<endl;	
			break;
		}
		cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: "<<nMergedRegionsInHierarchyLevel<<"/"<<regions.size()<<" regions aggregated at this level hierarchy..."<<endl;

	}//end while loop


	cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: "<<hierarchyLevel<<" hierarchy levels aggregated: N="<<regions.size()<<" regions left"<<endl;


	//## Final tag of composite regions according to a majority criterion, that is region tagged as bkg if the majority of sub-regions are non-significative
	cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Tagging significant regions..."<<endl;
	for(int i=0;i<regions.size();i++){
		int regionId= regions[i]->fId;
		int regionTag= regions[i]->fTag;
		int regionIndex_top= regionIdMap_top[regionId];//region index at the top of the hierarchy
		std::vector<int> subregionIds= regions[i]->fSubRegionIds;

		bool isSignificant= inputRegions[regionIndex_top]->fIsSignificative;
		bool isSalient= inputRegions[regionIndex_top]->fIsSalient;
		double salientFraction= 0;
		int nSubRegions= 1+subregionIds.size();
		if(isSalient) salientFraction++;
		bool isAnySignificant= false;
		if(isSignificant) isAnySignificant= true;

		for(unsigned int j=0;j<subregionIds.size();j++){
			int subRegionId= subregionIds[j];
			int subRegionIndex_top= regionIdMap_top[subRegionId];
			//bool subIsSignificant= fRegions[subRegionIndex_top]->fIsSignificative;
			//bool subIsSalient= fRegions[subRegionIndex_top]->fIsSalient;	
			bool subIsSignificant= inputRegions[subRegionIndex_top]->fIsSignificative;
			bool subIsSalient= inputRegions[subRegionIndex_top]->fIsSalient;
			if(subIsSalient) salientFraction++;
			if(subIsSignificant) isAnySignificant= true;
		}
		salientFraction/= (double)nSubRegions;
			
		if( (salientFraction>fSignificantSPRatio) || isAnySignificant ){
			cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Region id="<<regionId<<" tagged as significant as (significanceRatio="<<salientFraction<<">"<<fSignificantSPRatio<<")..."<<endl;
			regions[i]->fIsSignificative= true;
		}
		else {
			regions[i]->fIsSignificative= false;
		}	
	}//end loop regions
		
	//## Tag background?
	if(fUsePixelRatioCut){
		//Set region with maximum number of pixels
		for(int i=0;i<regions.size();i++){
			int regionId= regions[i]->fId;	
			int npix= regions[i]->fNPix;
			double pixelRatio= (double)npix/(double)N;
			if(pixelRatio>fPixelRatioCut) {
				cout<<"SLICSegmenter::SPHierarchicalMerger(): INFO: Region id="<<regionId<<" has too many pixels (ratio="<<pixelRatio<<">"<<fPixelRatioCut<<") and will be tagged as bkg..."<<endl;
				regions[i]->fIsSignificative= false;
			}
		}
	}//close if
	
	//## Copy final results to global data (fRegions, fPixelClusterIds)
	inputRegions.clear();
	inputRegions.assign(regions.begin(),regions.end());
	for(unsigned int i=0;i<labels.size();i++) {
		for(unsigned int j=0;j<labels[i].size();j++) {
			inputLabels[i][j]= labels[i][j];
		}
	}
	
	//## Get segmented image
	cout<<"SLICSegmenter::SuperpixelMerger(): INFO: Get segmented image..."<<endl;
	Img* finalMergedImg= SLICUtils::GetSegmentedImage(fInputImg,inputRegions,false,false,true);
	finalMergedImg->SetNameTitle("finalMergedImg","finalMergedImg");
	std::vector<TText*> regionTextList= SLICUtils::GetRegionTextLabels(inputRegions);

	//## Compute region contour info
	cout<<"SLICSegmenter::SuperpixelMerger(): INFO: Compute contour data..."<<endl;
	SLICUtils::SLICContourData* contourData= SLICUtils::ComputeBoundaryContours(fInputImg,inputLabels,inputRegions);

	//## Compute similarity matrix
	cout<<"SLICSegmenter::SuperpixelMerger(): INFO: Compute similarity data..."<<endl;
	SLICUtils::SLICSimilarityData* similarityData= SLICUtils::ComputeRegionSimilarity(fEdgeFilterImg,contourData,inputRegions,fSPMergingRegularization,fSPMergingIncludeSpatialPars);

	cout<<"SLICSegmenter::SuperpixelMerger(): INFO: Drawing results..."<<endl;
	TCanvas* SignalMergedPlot= new TCanvas("SignalMergedPlot","SignalMergedPlot");
	SignalMergedPlot->cd();
	finalMergedImg->Draw("COLZ");
	//for(unsigned int i=0;i<regionTextList.size();i++) regionTextList[i]->Draw("same");
	if(contourData){
		TGraph* contourGraph= (contourData->contour)->GetGraph();
		if(contourGraph) {
			contourGraph->SetLineColor(kRed);
			contourGraph->SetMarkerStyle(1);
			contourGraph->SetMarkerColor(kRed);
			contourGraph->Draw("P");
		}
	}
	
	
	return 0;

}//close SPHierarchicalMerger()



