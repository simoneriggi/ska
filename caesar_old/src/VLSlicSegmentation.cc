
#include <VLSlicSegmentation.h>
#include <OutlierDetector.h>

//#include "Img.h"

/*
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
*/

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

VLSlicSegmentation::VLSlicSegmentation() {

	fContours= 0;
	fLUT= 0;

	fImg= 0;
	fInputImg= 0;
	fImgStats= 0;
	fBackgroundImg= 0;
	fLaplImg= 0;
	fLaplImgStats= 0;
	fEdgeFilterImg= 0;
	fEdgeFilterImgStats= 0;
	
	fRegions.clear();

	fUseLogNormalizedImage= true;
	fSPMergingRatio= 0.3;
	fSPMergingRegularization= 0.5;
	fUse2ndNeighborsInSPMerging= true;
	fMinMergedSP= 1;
	fSPMergingDistThreshold= 0.1;
	fUseRobustParamsInSPMerging= false;

	fPixelRatioCut= 0.5;
	fTagSignificativeSP= false;
	fSPTaggingMethod= eSaliencyThresholding;
	fSignificantSPRatio= 0.5;
	fSaliencyThresholdFactor= 2;
	fSaliencyImg= 0;
	fSumSaliencyImg= 0;
	fSPMergingUseAdaptingDistThreshold= false;
	fSPMergingAdaptingDistThresholdScale= 500;

}//close constructor


VLSlicSegmentation::~VLSlicSegmentation() {
		
	/*
	while (!fRegions.empty()){
		cout<<"VLSlicSegmentation::~VLSlicSegmentation(): INFO: Deleting regions ..."<<endl;
  	Region* thisRegion = fRegions.back();
    fRegions.pop_back();
    if(thisRegion) delete thisRegion;
  }
	fRegions.clear();	
	*/


	cout<<"VLSlicSegmentation::~VLSlicSegmentation(): INFO: Deleting regions..."<<endl;
	for(unsigned int i=0;i<fRegions.size();i++){
		if(fRegions[i]){
			delete fRegions[i];
			fRegions[i]= 0;
		}
	}	
	fRegions.clear();
	
	cout<<"VLSlicSegmentation::~VLSlicSegmentation(): INFO: Deleting filter images..."<<endl;
	if(fLUT) fLUT->Delete();
	if(fLaplImg) fLaplImg->Delete();
	if(fEdgeFilterImg) fEdgeFilterImg->Delete();
	//if(fInputImg) fInputImg->Delete();
	cout<<"VLSlicSegmentation::~VLSlicSegmentation(): INFO: Deleting norm image..."<<endl;
	if(fImg) fImg->Delete();
	//cout<<"VLSlicSegmentation::~VLSlicSegmentation(): INFO: Deleting stats objects..."<<endl;
	//if(fImgStats) delete fImgStats;
	//if(fLaplImgStats) delete fLaplImgStats;
	//if(fEdgeFilterImgStats) delete fEdgeFilterImgStats;

	cout<<"VLSlicSegmentation::~VLSlicSegmentation(): INFO: Deleting bkg image..."<<endl;
	if(fBackgroundImg) fBackgroundImg->Delete();

	fPixelClusterIds.clear();
	fConnectedClusterList.clear();

	cout<<"VLSlicSegmentation::~VLSlicSegmentation(): INFO: Deleting contour..."<<endl;
	if(fContours) fContours->Delete();
	cout<<"VLSlicSegmentation::~VLSlicSegmentation(): INFO: End delete..."<<endl;
	
	
}//close constructor


int VLSlicSegmentation::Init(Img* inputImage){

	//## Clear previous data (if any), and re-initialize it
  fPixelClusterIds.clear();
	fRegions.clear();
	fConnectedClusterList.clear();
	if(fContours) fContours->Delete();
	if(fImg) fImg->Delete();
	//if(fImgStats) delete fImgStats;
	if(fLaplImg) fLaplImg->Delete();
	//if(fLaplImgStats) delete fLaplImgStats;
	if(fEdgeFilterImg) fEdgeFilterImg->Delete();
	//if(fEdgeFilterImgStats) delete fEdgeFilterImgStats;
	//if(fInputImg) fInputImg->Delete();
	if(fLUT) fLUT->Delete();

	//## Check input image
	if(!inputImage) {
		cerr<<"VLSlicSegmentation::Init(): ERROR: Null input image ptr!"<<endl;
		return -1;
	}
	fInputImg= inputImage;
	
	//## Compute image stats
	if(!fInputImg->HasStats()){
		if(fInputImg->ComputeStats(true,false,true)<0) {
			cerr<<"VLSlicSegmentation::Init(): ERROR: Cannot get input image stats!"<<endl;
			return -1;
		}
	}
	fImgStats= fInputImg->GetPixelStats();
	if(!fImgStats) return -1;
	
	//## Normalize image to range [1-256] and define color LUT for pixel mapping	
	double normmin= 1;
	double normmax= 256;
	if(fUseLogNormalizedImage) fImg= inputImage->GetLogNormalizedImage(normmin,normmax,true);
	else fImg= inputImage->GetNormalizedImage(normmin,normmax,true);

	int nlevels= 255;
	double zmin = fImg->GetMinimum();
	double zmax = fImg->GetMaximum();
	if(!fLUT) fLUT= new TH1D("LUT","LUT",nlevels,zmin,zmax);
	
	//## Compute laplacian images
	fLaplImg= fInputImg->GetLaplacianImage(true);
	if(!fLaplImg->HasStats()){
		if(fLaplImg->ComputeStats(true,false,true)<0) {
			cerr<<"VLSlicSegmentation::Init(): ERROR: Cannot get input laplacian image stats!"<<endl;
			return -1;
		}
	}
	fLaplImgStats= fLaplImg->GetPixelStats();
	if(!fLaplImgStats) return -1;

	//## Compute edge filtered images
	fEdgeFilterImg= fInputImg->GetKirschImage();
	if(!fEdgeFilterImg->HasStats()){
		if(fEdgeFilterImg->ComputeStats(true,false,true)<0) {
			cerr<<"VLSlicSegmentation::Init(): ERROR: Cannot get input edge filtered image stats!"<<endl;
			return -1;
		}
	}
	fEdgeFilterImgStats= fEdgeFilterImg->GetPixelStats();
	if(!fEdgeFilterImgStats) return -1;

	return 0;

}//close Init()


int VLSlicSegmentation::RunSegmentation(Img* inputImg,int regionSize,double regularization,int minRegionSize,bool mergeRegions,int mergeAlgoType){

	//## Check input image
	if(!inputImg) {
		cerr<<"VLSlicSegmentation::RunSegmentation(): ERROR: Null ptr to input image given!"<<endl;
		return -1;
	}

	//## Initialize data for segmentation tasks (LUT, vectors/lists, ...)
	//## Transform image for segmentation purposes
	cout<<"VLSlicSegmentation::RunSegmentation(): INFO: Initializing segmentation task..."<<endl;
	if( Init(inputImg)<0 ){
		cerr<<"VLSlicSegmentation::RunSegmentation(): ERROR: Init failed!"<<endl;	
		return -1;
	}

	//## Perform SLIC segmentation 
	cout<<"VLSlicSegmentation::RunSegmentation(): INFO: Generate superpixels..."<<endl;
	int status= RunSLICSegmentation(regionSize, regularization, minRegionSize);
	if(status<0){
		cerr<<"VLSlicSegmentation::RunSegmentation(): ERROR: Superpixel generation failed!"<<endl;
		return -1;
	}
	
	
	
	if(mergeRegions) {
		//## Merge regions with Hierarchical clustering
		if(mergeAlgoType==eHIER){
			int mergingStatus= HierarchicalMergeRegions();
			if(mergingStatus<0) {
				cerr<<"VLSlicSegmentation::RunSegmentation(): ERROR: Hierarchical region merging failed!"<<endl;
				return -1;
			}
		}
		//## Merge regions with DBSCAN
		else if(mergeAlgoType==eDBSCAN){	
			MergeRegions();
		}
		else{
			cerr<<"VLSlicSegmentation::RunSegmentation(): ERROR: Invalid merge algorithm selected!"<<endl;
			return -1;
		}
	}//close if merging
	

	return 0;

}//close RunSegmentation()



int VLSlicSegmentation::RunSLICSegmentation(int regionSize,double regularization, int minRegionSize){

	if(!fImg){
		cerr<<"VLSlicSegmentation::RunSLICSegmentation(): INFO: Null ptr to image!"<<endl;
		return -1;
	}
	
	//## Init data
	int Nx= fImg->GetNbinsX();
	int Ny= fImg->GetNbinsY();
	int numChannels= 1;//3;
	float colors[3];
	std::vector< std::vector<bool> > istaken;
	
	//float* image = new float[Nx*Ny*numChannels];//input image (when color mapping is used)
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

			/*
			if(w<=0) {
				for(int k=0;k<numChannels;k++){
					colors[k]= 255;//use white encoding for empty pixels
					image[i + j*Nx + Ny*Nx*k]= colors[k];
				}
			}//close if
			else{
				double lgw= log10(w);
				int LUTBinId= fLUT->FindBin(lgw);//use log mapping
				
				int colorId= TColor::GetColorPalette(LUTBinId); 
				gROOT->GetColor(colorId)->GetRGB(*colors,*(colors+1),*(colors+2));//colors are in range [0,1]
			
				for(int k=0;k<numChannels;k++){
					colors[k]*= 255;//convert colors to 0-255 
					image[i + j*Nx + Ny*Nx*k]= colors[k];
				}
 			}//close else
			*/

			cr.push_back(-1);
			nb.push_back(false);

		}//end loop Ny
		fPixelClusterIds.push_back(cr);
		istaken.push_back(nb);
	}//end loop Nx

	//## Generate the superpixels
	//## The algorithm will store the final segmentation in a one-dimensional array.        
	
  //## Compute edge map (gradient strength)
	//cout<<"VLSlicSegmentation::RunSLICSegmentation(): INFO: Computing edge map..."<<endl;
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
	cout<<"VLSlicSegmentation::RunSLICSegmentation(): INFO: Initializing K-means centers..."<<endl;
 
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

  // run k-means iterations
	cout<<"VLSlicSegmentation::RunSLICSegmentation(): INFO: Running "<<maxNumIterations<<" iterations..."<<endl;
 
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


	/*
	//## Print info
	for (region = 0 ; region < (signed)numRegions ; ++region) {
  	cout<<"VLSlicSegmentation::RunSLICSegmentation(): INFO: Region no. "<<region+1<<", npix="<<masses[region]<<" center(";
    for (i = (2 + numChannels)*region; i<(signed)(2 + numChannels)*(region + 1); ++i) {
    	cout<< centers[i]<<",";
    }
		cout<<")"<<endl;
  }
	*/

  vl_free(masses);
  vl_free(centers);
  vl_free(edgeMap);

	//## Eliminate small regions
  cout<<"VLSlicSegmentation::RunSLICSegmentation(): INFO: Eliminating small regions..."<<endl;
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
	//cout<<"VLSlicSegmentation::RunSLICSegmentation(): INFO: Delete allocated data..."<<endl;
	vl_free(cleaned);
  vl_free(segment);
	vl_free(image);

	//## Fill segmentation labels and compute stats/contours/...
	//cout<<"VLSlicSegmentation::RunSLICSegmentation(): INFO: Init pixel labels"<<endl;
	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};
	std::vector<TVector3> contours;
	contours.clear();
	contours.resize(0);

	for(unsigned int k=0;k<numRegions;k++) {
		fConnectedClusterList.push_back( std::vector<int>() );
		fRegions.push_back( new Region() );
		fRegions[k]->fId= -1;
		fRegions[k]->fImageSizeX= Nx;
		fRegions[k]->fImageSizeY= Ny;
		fRegions[k]->fImageMinX= fInputImg->GetXaxis()->GetXmin();
		fRegions[k]->fImageMaxX= fInputImg->GetXaxis()->GetXmax();
		fRegions[k]->fImageMinY= fInputImg->GetYaxis()->GetXmin();
		fRegions[k]->fImageMaxY= fInputImg->GetYaxis()->GetXmax();
		fRegions[k]->fImageMinS= fImgStats->min;
		fRegions[k]->fImageMaxS= fImgStats->max;
		fRegions[k]->fImageMinScurv= fLaplImgStats->min;
		fRegions[k]->fImageMaxScurv= fLaplImgStats->max;
		fRegions[k]->fImageMinSedge= fEdgeFilterImgStats->min;
		fRegions[k]->fImageMaxSedge= fEdgeFilterImgStats->max;
		fRegions[k]->fImageRMS= fImgStats->rms;			
	}
	//cout<<"VLSlicSegmentation::RunSLICSegmentation(): INFO: Filling pixel labels"<<endl;

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
				cout<<"VLSlicSegmentation::RunSLICSegmentation(): WARNING: Skip this pixel label: "<<pixelLabel<<endl;
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
			
			// Compare the pixel to its 8 neighbours to fill contour
			int nr_p = 0;
			std::vector<int> connectedRegionIds(numRegions,0);
			for (int k=0; k<8; k++) {
      	int ix_neighbor = ix + dx8[k];
				int iy_neighbor = iy + dy8[k];
				int id_neighbor= ix_neighbor + iy_neighbor*Nx;
        
        if (ix_neighbor>=0 && ix_neighbor<Nx && iy_neighbor>=0 && iy_neighbor<Ny) {
					int pixelLabel_neighbor= (int)segmentation[id_neighbor];
        
        	//if (istaken[ix_neighbor][iy_neighbor]==false && pixelLabel_neighbor!=pixelLabel && pixelLabel_neighbor>=0 ) {	
					if (pixelLabel_neighbor!=pixelLabel && pixelLabel_neighbor>=0 ) {	
          	nr_p++;
						connectedRegionIds[pixelLabel_neighbor]++;
          }
        }
      }//end loop neighbours
            
      //## Add the pixel to the contour list if desired
			//cout<<"Add the pixel to the contour list if desired?"<<endl;
     	if (nr_p >= 2) {
				contours.push_back(TVector3(ix,iy,binContent));	
				istaken[ix][iy] = true;

				for(unsigned int l=0;l<connectedRegionIds.size();l++){
					if( connectedRegionIds[l]>=2){ //4-connectivity 
						std::vector<int>::iterator it = std::find (fConnectedClusterList[pixelLabel].begin(), fConnectedClusterList[pixelLabel].end(), l);
						if( fConnectedClusterList[pixelLabel].empty() || it==fConnectedClusterList[pixelLabel].end() ) {
							fConnectedClusterList[pixelLabel].push_back(l);
							(fRegions[pixelLabel]->fNeighbourRegions).push_back(l);
						}//close if
					}//close if	
				}//end loop connected regions
     	}//close if
			
    }//end loop Ny
  }//end loop Nx

	//cout<<"VLSlicSegmentation::RunSLICSegmentation(): INFO: Freeing segmentation..."<<endl;
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

	cout<<"VLSlicSegmentation::RunSLICSegmentation(): INFO: nRegions after cleanup: "<<fRegions.size()<<endl;
	fRegionText= 0;
	fRegionTextList.clear();
	
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
			
		//fRegions[k]->Dump();

		cout<<"VLSlicSegmentation::RunSLICSegmentation(): INFO: Region no. "<<fRegions[k]->fId<<": N="<<nPix<<" C("<<Cx<<","<<Cy<<","<<Cz<<") color("<<fRegions[k]->fColor.X()<<","<<fRegions[k]->fColor.Y()<<","<<fRegions[k]->fColor.Z()<<"), connection(";
		for(unsigned int l=0;l<(fRegions[k]->fNeighbourRegions).size();l++) cout<<(fRegions[k]->fNeighbourRegions)[l]<<",";
		cout<<")"<<endl;

		//double textX= Cx + fImg->GetXaxis()->GetXmin();
		//double textY= Cy + fImg->GetYaxis()->GetXmin();
		double textX= Cx;
		double textY= Cy;

		fRegionText= new TText(textX,textY,Form("%d",regionId));
		fRegionText->SetTextSize(0.015);
		fRegionText->SetTextColor(kBlack);
		fRegionTextList.push_back(fRegionText);
	}//end loop regions

	
	
	//## Fill contours
	int nContourPoints= (int)contours.size();	
	//cout<<"VLSlicSegmentation::RunSLICSegmentation(): INFO: nContourPoints="<<nContourPoints<<endl;
	
	TString contourGraphName= Form("%s-SLICContours",std::string(fImg->GetName()).c_str());
	if(!fContours) fContours= new TGraph(nContourPoints);
  fContours->SetNameTitle(contourGraphName,contourGraphName);

  for (int i=0; i<nContourPoints; i++) {
		double cx= fImg->GetXaxis()->GetBinCenter(contours[i].X()+1);
		double cy= fImg->GetYaxis()->GetBinCenter(contours[i].Y()+1);
		fContours->SetPoint(i,cx, cy);
  }//end loop

	fContours->SetMarkerSize(1);
	fContours->SetMarkerStyle(1);
	fContours->SetLineColor(kRed);
	fContours->SetMarkerColor(kRed);	
	


	//## Get cluster colored image
	TString imgName= Form("%s-SLICColored",std::string(fImg->GetName()).c_str());
	Img* colored_img= GetClusterColoredImage(fImg,fRegions);
	colored_img->SetNameTitle(imgName,imgName);

	TString canvasName= Form("%s-SLICColoredPlot",std::string(fImg->GetName()).c_str());
	TCanvas* SlicColoredPlot= new TCanvas(canvasName,canvasName);
	SlicColoredPlot->cd();
	
	colored_img->Draw("COLZ");
	//for(unsigned int k=0;k<fRegionTextList.size();k++) fRegionTextList[k]->Draw("same");

	fContours->Draw("Psame");
	//fClusterCenterGraph->Draw("Psame");

	SlicColoredPlot->Update();
	

	return 0;

}//close VLSlicSegmentation::RunSLICSegmentation()



void VLSlicSegmentation::MergeRegions(int minNpts,double epsColor){

	//## Init clustered and visited
	int Nx= fImg->GetNbinsX();
	int Ny= fImg->GetNbinsY();
	int nRegions= (int)fRegions.size(); 
	std::vector< std::vector<Region*> > clusters;
	std::vector<bool> clustered(nRegions,false);
	std::vector<bool> visited(nRegions,false);
	std::vector<int> noise;
	std::vector<int> neighborPts;
	std::vector<int> neighborPts_;
	std::vector<int> regionLabels;
	int clusterCounter= 0;

	float colors[3]= {0,0,0};

	for(int i=0;i<nRegions;i++){
		regionLabels.push_back(-1);		
	}

	//clusters.push_back(std::vector<Pixel>()); //will stay empty?
	
	//## For each unvisited region center P
	cout<<"VLSlicSegmentation::MergeRegions(): INFO: Start region merging (nRegions="<<nRegions<<") with params: minNpts="<<minNpts<<", epsColor="<<epsColor<<endl;
	for(int i=0; i<nRegions; i++) {
		if(visited[i]) continue;

		visited[i] = true;
		int regionId= fRegions[i]->fId;
      
		//## Find Neighbors of region id in color space
		neighborPts = findNeighbors(regionId,epsColor);
		//cout<<"VLSlicSegmentation::MergeRegions(): INFO: "<<neighborPts.size()<<" neighbors for region no. "<<i<<" ..."<<endl;

		//## If the number of neighbors is too small mark as noise and skip point
    if(neighborPts.size() < minNpts){
    	noise.push_back(i);
			continue;
		}

		//## Create a cluster and add point P
		//cout<<"VLSlicSegmentation::MergeRegions(): INFO: Merged region detected including regions ("<<i<<", ";
		//for(unsigned l=0;l<neighborPts.size();l++) cout<<neighborPts[l]<<", ";	
		//cout<<")"<<endl;
			
		clustered[i] = true;//add point to cluster
		clusters.push_back( std::vector<Region*>() );
		clusters[clusterCounter].push_back(fRegions[i]);
		regionLabels[regionId]= clusterCounter;
		
      
		//## For each point P' in neighborPts
    for(unsigned int j=0; j<neighborPts.size(); j++){
			int neighborId= neighborPts[j];

      if(!visited[neighborId]) {//if P' is not visited mark as visited!
      	visited[neighborId] = true;
        neighborPts_ = findNeighbors(neighborId,epsColor);

				//for(unsigned l=0;l<neighborPts_.size();l++) cout<<neighborPts_[l]<<", ";	
				//cout<<")"<<endl;

        if((int)neighborPts_.size() >= minNpts) {
        	neighborPts.insert(neighborPts.end(),neighborPts_.begin(),neighborPts_.end());
        }
      }//close if
                
			// if P' is not yet a member of any cluster add P' to cluster c
      if(!clustered[neighborId]) {
				clusters[clusterCounter].push_back(fRegions[neighborId]);
				clustered[neighborId] = true;
				regionLabels[neighborId]= clusterCounter;//label region
			}
    }//end for loop neighborPts
     
		clusterCounter++;
    
	}//close loop on points

	cout<<"VLSlicSegmentation::MergeRegions(): INFO: "<<clusters.size()<<" region clusters found!"<<endl;
		
	if(clusters.size()<=0) {
		cout<<"VLSlicSegmentation::MergeRegions(): INFO: Merging completed without merged regions..."<<endl;
		return;
	}
		
	//## Relabel all pixels
	int nLabels= clusters.size();
	for(int i=0;i<nRegions;i++){
		if(regionLabels[i]==-1){//unlabeled, assign a new label starting from first available
			regionLabels[i]= nLabels;
			nLabels++;			
		}
		//cout<<"Region no. "<<i<<", new label="<<regionLabels[i]<<endl;
	}
	cout<<"VLSlicSegmentation::MergeRegions(): INFO: "<<nLabels<<" region clusters after region merging..."<<endl;
	
	//Reset previous vectors and re-label all pixels	
	for(int i=0;i<fRegions.size();i++){
		if(fRegions[i]) {
			delete fRegions[i];
			fRegions[i]= 0;
		}
	}
	fRegions.clear();
	fRegions.resize(0);

	for(int k=0;k<nLabels;k++){
		fRegions.push_back( new Region() );
		fRegions[k]->fId= -1;
		fRegions[k]->fImageSizeX= Nx;
		fRegions[k]->fImageSizeY= Ny;
		fRegions[k]->fImageMinS= fImgStats->min;
		fRegions[k]->fImageMaxS= fImgStats->max;
		fRegions[k]->fImageRMS= fImgStats->rms;	
	}//end loop new regions

	for(int i=0;i<Nx;i++){
		double binX= fImg->GetXaxis()->GetBinCenter(i+1);
		for(int j=0;j<Ny;j++){	
			double binY= fImg->GetYaxis()->GetBinCenter(j+1);
			int binId= fImg->GetBin(i+1,j+1);
			double S= fInputImg->GetBinContent(i+1,j+1);
			double binContent= fImg->GetBinContent(i+1,j+1);//already in log scale if log mapping is activated
			//double binContent_lg= log10(binContent);
			//int LUTBinId= fLUT->FindBin(binContent_lg);//use log mapping
			int LUTBinId= fLUT->FindBin(binContent);//use log mapping
			int colorId= TColor::GetColorPalette(LUTBinId); 
			gROOT->GetColor(colorId)->GetRGB(*colors,*(colors+1),*(colors+2));//colors are in range [0,1]		
			for(int kk=0;kk<3;kk++) colors[kk]*= 255;//convert colors to 0-255 

			double S_curv= fLaplImg->GetBinContent(i+1,j+1);
			double S_edge= fEdgeFilterImg->GetBinContent(i+1,j+1);

			int oldRegionId= fPixelClusterIds[i][j];
			int newRegionId= regionLabels[oldRegionId];

			//Create a pixel
			Region::Pixel aPixel;
			aPixel.id= binId;
			aPixel.S= S;
			aPixel.color= TVector3(colors[0],colors[1],colors[2]);	
			//aPixel.x= i;
			//aPixel.y= j;
			aPixel.x= binX;
			aPixel.y= binY;
			aPixel.ix= i;
			aPixel.iy= j;
			aPixel.isOnEdge= false;
			aPixel.distanceToEdge= 1.e+99;
			aPixel.S_curv= S_curv;
			aPixel.S_edge= S_edge;
			fPixelClusterIds[i][j]= newRegionId;
			
			fRegions[newRegionId]->fId= newRegionId;
			if(S!=0) fRegions[newRegionId]->AddPixel(aPixel);
				
		}//end loop y
	}//end loop x

	
	fRegionTextList.clear();
	fRegionText= 0;

	for(int k=0;k<nLabels;k++){
		int nPix= fRegions[k]->fNPix;
		int regionId= fRegions[k]->fId;
		if(nPix<=0) continue;

		fRegions[k]->ComputeParameters();

		double Cx= fRegions[k]->fX0;
		double Cy= fRegions[k]->fY0;
		double Cz= fRegions[k]->fMean;

		double regionMean= fRegions[k]->fMean;
		double regionMean_lg= log10(regionMean);
		int LUTBinId= fLUT->FindBin(regionMean);//use log mapping
		if(fUseLogNormalizedImage) LUTBinId= fLUT->FindBin(regionMean_lg);//use log mapping
		int colorId= TColor::GetColorPalette(LUTBinId); 
		gROOT->GetColor(colorId)->GetRGB(*colors,*(colors+1),*(colors+2));//colors are in range [0,1]		
		for(int kk=0;kk<3;kk++) colors[kk]*= 255;//convert colors to 0-255 
		(fRegions[k]->fColor).SetXYZ(colors[0],colors[1],colors[2]);

		cout<<"VLSlicSegmentation::MergeRegions(): INFO: Region id "<<regionId<<": N="<<nPix<<" C("<<Cx<<","<<Cy<<","<<Cz<<") color("<<(fRegions[k]->fColor).X()<<","<<(fRegions[k]->fColor).Y()<<","<<(fRegions[k]->fColor).Z()<<"), connection(";
		cout<<endl;
			
		//double textX= Cx + fImg->GetXaxis()->GetXmin();
		//double textY= Cy + fImg->GetYaxis()->GetXmin();
		double textX= Cx;
		double textY= Cy;
		fRegionText= new TText(textX,textY,Form("%d",regionId));
		fRegionText->SetTextSize(0.015);
		fRegionText->SetTextColor(kBlack);
		fRegionTextList.push_back(fRegionText);	
	}//end loop new regions	
	
	
	
	//## Recompute contours
	TGraph* clusterContours= ComputeClusterContours(fImg,fPixelClusterIds,fRegions);
	
	
	//## Recolor image with cluster means	
	TString imgName= Form("%s-MergedRegionImg",std::string(fImg->GetName()).c_str());
	Img* ColoredImg= GetClusterColoredImage(fImg,fRegions);
	ColoredImg->SetNameTitle(imgName,imgName);

	TString canvasName= Form("%s-MergedRegionImgPlot",std::string(fImg->GetName()).c_str());
	TCanvas* MergeRegionColoredPlot= new TCanvas(canvasName,canvasName);
	MergeRegionColoredPlot->cd();
		
	if(ColoredImg) ColoredImg->Draw("COLZ");
	//for(unsigned int k=0;k<fRegionTextList.size();k++) fRegionTextList[k]->Draw("same");
	if(clusterContours) clusterContours->Draw("Psame");
	

}//close VLSlicSegmentation::MergeRegions()


Region* VLSlicSegmentation::FindSimpleBackgroundRegion(int nPixThreshold){

	//## Check if there are regions
	if(fRegions.size()<=0) return 0;

	Region* bkgRegion= new Region;
	int nPix_max= -1.e+99;
	int BkgRegionId= -1;
	bool isFound= false;

	for(unsigned int i=0;i<fRegions.size();i++){
		int nPix= fRegions[i]->fNPix;
		if(nPix<=0) continue;
		int RegionId= fRegions[i]->fId;
		if(nPix>nPixThreshold && nPix>nPix_max){
			BkgRegionId= RegionId;
			*bkgRegion= *fRegions[i];
			isFound= true;
		}
	}//end loop regions
	
	if(!isFound){
		cerr<<"VLSlicSegmentation::FindSimpleBackgroundRegion(): WARNING: Cannot find the background region given the current criteria, returning null ptr!"<<endl;
		return 0;
	}
	
	cout<<"VLSlicSegmentation::FindSimpleBackgroundRegion(): INFO: Background region id: "<<BkgRegionId<<endl;

	return bkgRegion;

}//close VLSlicSegmentation::FindSimpleBackgroundRegion()



Region* VLSlicSegmentation::FindBackgroundRegion(double CL,bool includeSpatialPar,bool includeCurvPar,bool useOnlyPosExcess,bool useRobustPars){

	//## Check if there are regions
	if(fRegions.size()<=0) return 0;
	
	Region* bkgRegion= new Region;
	int BkgRegionId= -1;
	
	//## Find significant regions
	ApplyMahalanobisThresholding(fRegions,CL,includeSpatialPar,includeCurvPar,useOnlyPosExcess,useRobustPars);

	//## Fill background region
	int nMergedRegions= 0;
	for(unsigned int i=0;i<fRegions.size();i++){
		int nPix= fRegions[i]->fNPix;
		int RegionId= fRegions[i]->fId;
		bool isSignificant= fRegions[i]->fIsSignificative;
		cout<<"Region id="<<RegionId<<": isSignificant? "<<isSignificant<<endl;			

		if(!isSignificant){//Add this region to the background
			if(nMergedRegions==0) *bkgRegion= *fRegions[i];
			else bkgRegion->AddRegion(fRegions[i],true);
			nMergedRegions++;
		}
	}//end loop regions
		

	//## Compute parameters of background region
	if(nMergedRegions<=0){
		cerr<<"VLSlicSegmentation::FindBackgroundRegion(): WARN: Cannot find the background region given the current criteria, returning null ptr!"<<endl;
		return 0;
	}
	cout<<"VLSlicSegmentation::FindBackgroundRegion(): INFO: nRegions in bkg="<<nMergedRegions<<endl;
	bool computeContours= false;
	bool computeRobustStats= useRobustPars;
	bool forceRecomputing= true;
	bkgRegion->ComputeParameters(computeContours,computeRobustStats,forceRecomputing);

	//## Fill background image
	if(!fBackgroundImg) fBackgroundImg= (Img*)fImg->Clone("BackgroundImg");
	fBackgroundImg->Reset();

	std::vector<Region::Pixel> bkgPixCollection= bkgRegion->fPixelCollection;
	for(unsigned int i=0;i<bkgPixCollection.size();i++){
		int gBinId= bkgPixCollection[i].id;
		double S= fInputImg->GetBinContent(gBinId);
		fBackgroundImg->SetBinContent(gBinId,S);	
	}

	return bkgRegion;

}//close FindBackgroundRegion()




//std::pair<TMatrixD,TMatrixD> VLSlicSegmentation::ComputeRegionSimilarity(std::vector<Region*> regions,std::map<int,int> mapping,double beta){
VLSlicSegmentation::SLICSimilarityData* VLSlicSegmentation::ComputeRegionSimilarity(std::vector<Region*> regions,std::map<int,int> mapping,double beta){

	int nRegions= (int)regions.size();	

	//## Init data
	SLICSimilarityData* SimilarityData= new SLICSimilarityData;
	(SimilarityData->DissimilarityMatrix)= new TMatrixD(nRegions,nRegions);
	(SimilarityData->DissimilarityMatrix)->Zero();

	(SimilarityData->SimilarityMatrix)= new TMatrixD(nRegions,nRegions);
	(SimilarityData->SimilarityMatrix)->Zero();

	(SimilarityData->AbsDissimilarityMatrix)= new TMatrixD(nRegions,nRegions);
	(SimilarityData->AbsDissimilarityMatrix)->Zero();
	/*
	TMatrixD DissimilarityMatrix(nRegions,nRegions);
	DissimilarityMatrix.Zero();

	TMatrixD SimilarityMatrix(nRegions,nRegions);
	SimilarityMatrix.Zero();
	*/
	cout<<"VLSlicSegmentation::ComputeRegionSimilarity(): INFO: Start region merging (nRegions="<<nRegions<<")"<<endl;

	//## Find min & max distance for normalization	
	double Dmin_norm= 0;
	double Dmax_norm= 1;
	double Dmin= 1.e+99;
	double Dmax= -1.e+99;
	double Emin_norm= 0;
	double Emax_norm= 1;
	double Emin= 1.e+99;
	double Emax= -1.e+99;
	std::vector<double> DList;

	for(int i=0; i<nRegions; i++) {
		Region* thisRegion= regions[i];
		std::vector<Region::NeighborInfo> neighbours= thisRegion->fNeighbourRegionInfo;
		for(unsigned int j=0;j<neighbours.size();j++){
			Region::NeighborInfo info= neighbours[j];	
			int order= info.order;
			double E= info.Edgeness;
			double D= info.D;
			DList.push_back(D);
			if(D<Dmin) Dmin= D;
			if(D>Dmax) Dmax= D;
			if(E<Emin && order==1) Emin= E;
			if(E>Emax && order==1) Emax= E;
		}
	}//end loop regions

	double Dmedian= Utils::GetMedian(DList,true);
	double Dmad= Utils::GetMAD(DList,Dmedian);
	double Dmedianrms= Dmad*1.4826;
	SimilarityData->Dmedian= Dmedian;
	SimilarityData->Dmedianrms= Dmedianrms;
	SimilarityData->Dmin= Dmin;
	SimilarityData->Dmax= Dmax;
	SimilarityData->DList= DList;
	cout<<"VLSlicSegmentation::ComputeRegionSimilarity(): INFO: Dmedian="<<Dmedian<<" Dmedianrms="<<Dmedianrms<<" Dmin/Dmax: "<<Dmin<<"/"<<Dmax<<", Emin/Emax="<<Emin<<"/"<<Emax<<endl;

	//## Fill similarity matrix
	for(int i=0; i<nRegions; i++) {
		Region* thisRegion= regions[i];
		int thisRegionId= thisRegion->fId;
		std::vector<Region::NeighborInfo> neighbours= thisRegion->fNeighbourRegionInfo;

		for(unsigned int j=0;j<neighbours.size();j++){
			Region::NeighborInfo info= neighbours[j];
			double E= info.Edgeness;
			double E_norm= Emin_norm + (Emax_norm-Emin_norm)*(E-Emin)/(Emax-Emin);
			if(info.order>1) E_norm= Emax_norm;//fixed penalty for 2nd order neighbors
			
			double D= info.D;
			double D_norm= Dmin_norm + (Dmax_norm-Dmin_norm)*(D-Dmin)/(Dmax-Dmin);
			
			//== NORM DIST ==
			double dist= D_norm;
			if(nRegions>2 && Emax>Emin) dist+= beta*E_norm;
			//if(nRegions>2 && Emax>Emin){
			//	dist= (1-beta)*D_norm + beta*E_norm;
			//}			

			
			dist+= SMALL_NUMBER;//to avoid dividing by zero!

			int neighbourId= info.id;
			int neighbourIndex= mapping[neighbourId];
			//DissimilarityMatrix(i,neighbourIndex)= dist;
			//SimilarityMatrix(i,neighbourIndex)= 1./dist;
			(*(SimilarityData->DissimilarityMatrix))(i,neighbourIndex)= dist;
			(*(SimilarityData->SimilarityMatrix))(i,neighbourIndex)= 1./dist;
			(*(SimilarityData->AbsDissimilarityMatrix))(i,neighbourIndex)= D;
		}//end loop neighbors
	}//end loop regions
	

	//Normalize similarity matrix by rows
	/*
	for(int i=0;i<SimilarityMatrix.GetNrows();i++){
		double sum= 0;	
		for(int j=0;j<SimilarityMatrix.GetNcols();j++) sum+= SimilarityMatrix(i,j);
		if(sum!=0) {
			for(int j=0;j<SimilarityMatrix.GetNcols();j++) SimilarityMatrix(i,j)/= sum;
		}
	}//end loop rows
	*/
	for(int i=0;i<(SimilarityData->SimilarityMatrix)->GetNrows();i++){
		double sum= 0;	
		for(int j=0;j<(SimilarityData->SimilarityMatrix)->GetNcols();j++) sum+= (*(SimilarityData->SimilarityMatrix))(i,j);
		if(sum!=0) {
			for(int j=0;j<(SimilarityData->SimilarityMatrix)->GetNcols();j++) (*(SimilarityData->SimilarityMatrix))(i,j)/= sum;
		}
	}//end loop rows

	
	//cout<<"== Similarity Matrix =="<<endl;
	//SimilarityMatrix.Print();
	//cout<<"======================="<<endl;

	//cout<<"== Dissimilarity Matrix =="<<endl;
	//DissimilarityMatrix.Print();
	//cout<<"======================="<<endl;
	
	//std::vector<TMatrixD> MatrixList;
	//MatrixList.push_back(SimilarityMatrix);
	//MatrixList.push_back(DissimilarityMatrix);

	//return std::make_pair(SimilarityMatrix,DissimilarityMatrix);

	return SimilarityData;

}//close ComputeSimilarity()


int VLSlicSegmentation::TagRegions(std::vector<Region*> regions,double saliencyThresholdFactor,double CL,bool includeSpatialPar,bool includeCurvPar){

	cout<<"VLSlicSegmentation::TagRegions(): INFO: Tag regions..."<<endl;
	int nRegions= (int)regions.size();
	
	//## Compute saliency
	double saliencyAvg= 0;
	for(int i=0;i<nRegions;i++){
		double sum= 0;
		for(int j=0;j<nRegions;j++){
			if(i==j) continue;
			double diss= regions[i]->GetDissimilarity(regions[j],false,true);
			sum+= diss;
		}//end loop regions
		double S= 0;
		if(sum!=0) S= 1.-exp(-sum/(double)(nRegions));
		regions[i]->fSaliency= S;
		saliencyAvg+= S;
	}//end loop regions
	if(saliencyAvg>0) saliencyAvg/= (double)nRegions;


	//## Compute robust mean/covariance/Mahalanobis distance
	int nPars= 2;
	if(includeSpatialPar) nPars+= 2;
	if(includeCurvPar) nPars+= 2;
	TMatrixD* data_matrix= new TMatrixD(nRegions,nPars);
	for(int i=0;i<nRegions;i++){
		TVectorD regionPars= fRegions[i]->GetParamVector(includeSpatialPar,includeCurvPar);//only color params (no centroids!)
		for(int j=0;j<regionPars.GetNoElements();j++)	{
			(*data_matrix)(i,j)= regionPars(j);	
		}//end loop par dim
	}//end loop regions
		
	//## Perform outlier analysis and retrieve results
	OutlierDetector outlierDetector;
	if( outlierDetector.FindOutliers(data_matrix)<0 ) {
		cerr<<"VLSlicSegmentation::TagRegions(): ERROR: Outlier detection failed!"<<endl;
		return -1;
	}

	std::vector<int> data_flags= outlierDetector.GetDataFlags();
	std::vector<int> outliers_ids= outlierDetector.GetOutliers();
	std::vector<double> MHDs= outlierDetector.GetDistances();
	TMatrixD* robustCov= outlierDetector.GetRobustCov();
	TVectorD* robustMean= outlierDetector.GetRobustMean();
	double cutoffValue= outlierDetector.GetCutoff();
		
	//## Tag significative regions according to Mahalanobis distance
	//## CL: 0.975 (conventional) 0.997 (3sigma) 0.954 (2sigma), 0.683 (1sigma)
	//## Tag salient regions
	double chi2quantile= ROOT::Math::chisquared_quantile(CL,nPars);
	double outlierCut= chi2quantile;
	double saliencyThr= saliencyThresholdFactor*saliencyAvg;
	cout<<"VLSlicSegmentation::TagRegions(): INFO: <saliency>="<<saliencyAvg<<" saliencyThr="<<saliencyThr<<endl;

	for(int i=0;i<nRegions;i++){
		int regionId= regions[i]->fId;
		int outlierFlag= data_flags[i];	
		double MHD= MHDs[i];
		bool isOutlier= (MHD>outlierCut);
		bool isPositiveOutlier= ( (*data_matrix)(i,0)>(*robustMean)(0) );
		double S= regions[i]->fSaliency;
		bool isSalient= (S>saliencyThr);
		regions[i]->fMahalanobisDistance= MHD;
		
		if(isSalient && isPositiveOutlier) regions[i]->fIsSalient= true;
		else regions[i]->fIsSalient= false;
		if(isOutlier && isPositiveOutlier) regions[i]->fIsSignificative= true;
		else regions[i]->fIsSignificative= false;

		cout<<"VLSlicSegmentation::TagRegions(): INFO: Region id="<<regionId<<": MHD="<<MHD<<" (cut="<<outlierCut<<",cutoffValue="<<cutoffValue<<") isOutlier?"<<isOutlier<<" isPositiveOutlier? "<<isPositiveOutlier<<" outlierFlag="<<outlierFlag<<" isSalient?"<<isSalient<<endl;		
	}//end loop regions

	return 0;

}//close TagRegions()


Img* VLSlicSegmentation::ComputeSaliencyMap(std::vector<Region*> regions,int knn){

	//## Compute saliency
	cout<<"VLSlicSegmentation::ApplySaliencyThresholding(): INFO: Applying saliency thresholding..."<<endl;
	int nRegions= (int)regions.size();
	if(nRegions<=0 || !fInputImg) return 0;
	
	Img* saliencyImg= (Img*)fInputImg->Clone("saliencyImg");
	saliencyImg->Reset();

	for(int i=0;i<nRegions;i++){

		double Saliency= 0;

		//Find the closest k-neighbors in apperance
		if(knn<nRegions){
			std::vector<double> distList;
			for(int j=0;j<nRegions;j++){
				//if(i==j) continue;
				std::pair<double,double> dists= regions[i]->GetDistance(regions[j],false,true);
				double dist_appearance= dists.first;
				distList.push_back(dist_appearance);
			}//end loop regions

			std::vector<double> sorted;
			std::vector<size_t> sort_index;//sorting index
			Utils::sort( distList,sorted,sort_index);

			//Compute saliency over k-th neighbors
			double sum= 0;
			for(int k=0;k<knn;k++){
				size_t index= sort_index[k];
				double diss= regions[i]->GetDissimilarity(regions[index],false,true);
				sum+= diss;
			}	
			if(sum!=0) Saliency= 1.-exp(-sum/(double)(knn));
		}//close if knn
		else {

			double sum= 0;
			for(int j=0;j<nRegions;j++){
				if(i==j) continue;
				double diss= regions[i]->GetDissimilarity(regions[j],false,true);
				sum+= diss;
			}//end loop regions
			if(sum!=0) Saliency= 1.-exp(-sum/(double)(nRegions-1));

		}//close else

		//Fill image
		int nPixelsInRegion= (int)regions[i]->fNPix;
		for(int j=0;j<nPixelsInRegion;j++){//loop on pixels inside region
			int thisPixelId= (regions[i]->fPixelCollection)[j].id;
			saliencyImg->SetBinContent(thisPixelId,Saliency);
		}//end loop pixels in region	

	}//end loop regions

	return saliencyImg;

}//close GetSaliencyMap()




int VLSlicSegmentation::ApplySaliencyThresholding(std::vector<Region*> regions,double saliencyThresholdFactor){

	//## Compute saliency
	cout<<"VLSlicSegmentation::ApplySaliencyThresholding(): INFO: Applying saliency thresholding..."<<endl;
	int nRegions= (int)regions.size();
	
	double saliencyAvg= 0;
	for(int i=0;i<nRegions;i++){
		double sum= 0;
		for(int j=0;j<nRegions;j++){
			if(i==j) continue;
			double diss= regions[i]->GetDissimilarity(regions[j],false,true);
			sum+= diss;
		}//end loop regions
		double S= 0;
		if(sum!=0) S= 1.-exp(-sum/(double)(nRegions));
		regions[i]->fSaliency= S;
		saliencyAvg+= S;
	}//end loop regions

	if(saliencyAvg>0) saliencyAvg/= (double)nRegions;
		
	//## Tag regions as significative if above saliency threshold
	double thr= saliencyThresholdFactor*saliencyAvg;
	cout<<"VLSlicSegmentation::ApplySaliencyThresholding(): INFO: <saliency>="<<saliencyAvg<<" thr="<<thr<<endl;
	for(int i=0;i<nRegions;i++){
		int regionId= regions[i]->fId;	
		double S= regions[i]->fSaliency;
		if(S>thr) {
			cerr<<"VLSlicSegmentation::ApplySaliencyThresholding(): INFO: Region id="<<regionId<<" tagged as significant (Saliency="<<S<<">"<<thr<<")..."<<endl;
			regions[i]->fIsSalient= true;
		}
		else regions[i]->fIsSalient= false;
	}//end loop regions

	return 0;

}//close ApplySaliencyThresholding()


int VLSlicSegmentation::ApplyMahalanobisThresholding(std::vector<Region*> regions,double CL,bool includeSpatialPar,bool includeCurvPar,bool useOnlyPosExcess, bool useRobustPars){

	//## Init parameter data matrix
	int nRegions= (int)regions.size();
	int nPars= 2;
	if(includeSpatialPar) nPars+= 2;
	if(includeCurvPar) nPars+= 2;
	TMatrixD* data_matrix= new TMatrixD(nRegions,nPars);
	for(int i=0;i<nRegions;i++){
		TVectorD regionPars= fRegions[i]->GetParamVector(includeSpatialPar,includeCurvPar,useRobustPars);//only color params (no centroids!)
		for(int j=0;j<regionPars.GetNoElements();j++)	{
			(*data_matrix)(i,j)= regionPars(j);	
		}//end loop par dim
	}//end loop regions
		
	//## Perform outlier analysis and retrieve results
	OutlierDetector outlierDetector;
	if( outlierDetector.FindOutliers(data_matrix)<0 ) {
		cerr<<"VLSlicSegmentation::ApplyMahalanobisThresholding(): ERROR: Outlier detection failed!"<<endl;
		return -1;
	}

	std::vector<int> data_flags= outlierDetector.GetDataFlags();
	std::vector<int> outliers_ids= outlierDetector.GetOutliers();
	std::vector<double> MHDs= outlierDetector.GetDistances();
	TMatrixD* robustCov= outlierDetector.GetRobustCov();
	TVectorD* robustMean= outlierDetector.GetRobustMean();
	double cutoffValue= outlierDetector.GetCutoff();
		
	//## Tag significative regions according to Mahalanobis distance
	//## CL: 0.975 (conventional) 0.997 (3sigma) 0.954 (2sigma), 0.683 (1sigma)
	double chi2quantile= ROOT::Math::chisquared_quantile(CL,nPars);
	double outlierCut= chi2quantile;
	
	for(int i=0;i<nRegions;i++){
		int regionId= regions[i]->fId;
		int outlierFlag= data_flags[i];	
		double MHD= MHDs[i];
		bool isOutlier= (MHD>outlierCut);
		bool isPositiveOutlier= ( (*data_matrix)(i,0)>(*robustMean)(0) );
		bool isSelectedOutlier= (useOnlyPosExcess && isOutlier && isPositiveOutlier) || (!useOnlyPosExcess && isOutlier);

		cout<<"VLSlicSegmentation::ApplyMahalanobisThresholding(): INFO: Region id="<<regionId<<": MHD="<<MHD<<" (cut="<<outlierCut<<",cutoffValue="<<cutoffValue<<") isOutlier?"<<isOutlier<<" isPositiveOutlier? "<<isPositiveOutlier<<" outlierFlag="<<outlierFlag<<endl;

		regions[i]->fMahalanobisDistance= MHD;
	
		if(isSelectedOutlier){
			regions[i]->fIsSignificative= true;
			cout<<"VLSlicSegmentation::ApplyMahalanobisThresholding(): INFO: Region id="<<regionId<<" tagged as significative!"<<endl;
		}
		else{
			regions[i]->fIsSignificative= false;
		}
	}//end loop regions

	return 0;

}//close ApplyMahalanobisThresholding()

int VLSlicSegmentation::HierarchicalMergeRegions(){

	//Init timers
	std::chrono::high_resolution_clock::time_point tstart = high_resolution_clock::now();
	double elapseTime_tot= 0;
	double elapseTime_init= 0;
	double elapseTime_findNeighbors= 0;
	double elapseTime_computeSimMatrix= 0;
	double elapseTime_computePageRank= 0;
	double elapseTime_selectMergingRegions= 0;
	double elapseTime_copyRegions= 0;
	
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
  
	int nRegions= (int)fRegions.size();
	if(nRegions<=0) return -1;

	//## Set algo options
	int nMinRegions= fMinMergedSP;
	double mergeRatio= fSPMergingRatio;
	bool use2ndOrderNeighbors= fUse2ndNeighborsInSPMerging;
	double HThreshold= fSPMergingDistThreshold;
	bool useRobustParams= fUseRobustParamsInSPMerging;
	double regFactor= fSPMergingRegularization;
	  
	
	//## Select significative regions
	if(fTagSignificativeSP){
		TagRegions(fRegions,fSaliencyThresholdFactor,0.975,false,false);
	}//close if

	//## Set region list to be passed to algorithm
	cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Copying initial superpixel partition to tmp partition..."<<endl;
	std::vector<Region*> regions;
	regions.assign(fRegions.begin(),fRegions.end());
	fSPMergingInfo.clear();
	fSaliencyImg= ComputeSaliencyMap(regions);
	fSaliencyImg->SetNameTitle("saliencyMap","saliencyMap");
	
	fSumSaliencyImg= (Img*)fSaliencyImg->Clone("saliencySumMap");
	fSumSaliencyImg->SetNameTitle("saliencySumMap","saliencySumMap");

	//Init norm consts
	double NormMin= 0;
	double NormMax= 1;
	double A= NormMin - (NormMax-NormMin)*regions[0]->fImageMinS/(regions[0]->fImageMaxS-regions[0]->fImageMinS);
	double B= (NormMax-NormMin)/(regions[0]->fImageMaxS-regions[0]->fImageMinS);
	double Acurv= NormMin - (NormMax-NormMin)*regions[0]->fImageMinScurv/(regions[0]->fImageMaxScurv-regions[0]->fImageMinScurv);	
	double Bcurv= (NormMax-NormMin)/(regions[0]->fImageMaxScurv-regions[0]->fImageMinScurv);
	double Nx= fInputImg->GetNbinsX();
	double Ny= fInputImg->GetNbinsY();
	double N= Nx*Ny;

	//Create the mapping of regionId and vector index
	//Compute the initial MSE
	cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Create the mapping of regionId and region list index..."<<endl;
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
	}


	high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto dt_init = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
 	cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: init dt(usec)="<<dt_init<<endl;
	
	//## Run a hierarchical clustering till all segments are merged in one
	double d= 0.85;
	double tol= 1.e-4;
	int hierarchyLevel= 0;
	double AbsDissMedian= 0;
	double AbsDissMedianRMS= 0;
	double AbsDissMin= 0;
	double AbsDissMax= 0;
	double AbsDissMedian0= 0;
	double AbsDissMedianRMS0= 0;
	int nMergedRegionsInHierarchyLevel= 0;
	cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Starting hierarchical merging..."<<endl;	

	while(regions.size()>nMinRegions){
		nMergedRegionsInHierarchyLevel= 0;

		//Merge info for this hierarchy
		SPMergingInfo spMergeInfo;
		spMergeInfo.levelId= hierarchyLevel;
		spMergeInfo.NR= (int)regions.size();
		spMergeInfo.MSE= 0;
		spMergeInfo.DissMin= 0;
		spMergeInfo.DissMax= 0;
	
		//## Find neighbors for initial segmentation
		cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Finding region neighbors at hierarchy level "<<hierarchyLevel<<", NR="<<regions.size()<<" ..."<<endl;
		t1 = high_resolution_clock::now();
		findNeighbors_v2(regions,regionIdMap,labels,use2ndOrderNeighbors,useRobustParams);

		double MSE= 0.;
		for(unsigned int k=0;k<regions.size();k++){		
			regions[k]->SortNeighborInfo(true,true,regFactor);//ascending order	& tot dist (D+Edgeness)	
			double M2= regions[k]->fM2;
			MSE+= M2;
		}//end loop regions
		MSE/= N;
		spMergeInfo.MSE= MSE;

		t2 = high_resolution_clock::now();
		auto dt_findNeighbor= std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
		elapseTime_findNeighbors+= dt_findNeighbor;
		//cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: findNeighbors dt(usec)="<<dt_findNeighbor<<endl;
		
		//## Compute similarity matrix of current segmented regions
		cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Compute region similarity at hierarchy level "<<hierarchyLevel<<" ..."<<endl;
		t1 = high_resolution_clock::now();	
		
		//std::pair<TMatrixD,TMatrixD> WList= ComputeRegionSimilarity(regions,regionIdMap,regFactor);
		//TMatrixD AdjacencyMatrix= WList.first;
		//TMatrixD DissimilarityMatrix= WList.second;
		SLICSimilarityData* similarityData= ComputeRegionSimilarity(regions,regionIdMap,regFactor);
		TMatrixD* AdjacencyMatrix= similarityData->SimilarityMatrix;
		TMatrixD* DissimilarityMatrix= similarityData->DissimilarityMatrix;
		TMatrixD* AbsDissimilarityMatrix= similarityData->AbsDissimilarityMatrix;
		std::vector<double> DList= similarityData->DList;
		AbsDissMedian= similarityData->Dmedian;
 		AbsDissMedianRMS= similarityData->Dmedianrms;
 		AbsDissMin= similarityData->Dmin;
		AbsDissMax= similarityData->Dmax;
		if(hierarchyLevel==0){
			AbsDissMedian0= AbsDissMedian;
			AbsDissMedianRMS0= AbsDissMedianRMS;
		}
		std::vector<double> MergedDList;
	
		t2 = high_resolution_clock::now();
		auto dt_computeSimMatrix= std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
		elapseTime_computeSimMatrix+= dt_computeSimMatrix;
		//cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: computeSimMatrix dt(usec)="<<dt_computeSimMatrix<<endl;
		

		//## Compute page rank of segments
		//cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Compute page rank at hierarchy level "<<hierarchyLevel<<" ..."<<endl;	
		t1 = high_resolution_clock::now();
		//std::vector<double> ranks= ComputePageRank(AdjacencyMatrix.T(),d,tol);//pass transpose of similarity matrix
		std::vector<double> ranks= ComputePageRank(AdjacencyMatrix->T(),d,tol);//pass transpose of similarity matrix
		if(ranks.size()<=0){
			cerr<<"VLSlicSegmentation::HierarchicalMergeRegions(): WARN: PageRank failed, cannot perform region merging!"<<endl;
			return -1;
		}
		t2 = high_resolution_clock::now();
		auto dt_computePageRank= std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
		elapseTime_computePageRank+= dt_computePageRank;
		//cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: computePageRank dt(usec)="<<dt_computePageRank<<endl;
		
		//Original code computes a weighted mean between previous rank vector and current vector (NOT EXPLAINED IN THE PAPER!)
		// p_new= 0.7 p_old + 0.3 p
		

		t1 = high_resolution_clock::now();
		
		//Sort ranks
		std::vector<size_t> sort_index;//sorting index
		std::vector<double> ranks_sorted;
		Utils::sort_descending(ranks,ranks_sorted,sort_index);
		
		//Loop over sorted ranks and select regions to be merged
		int nMergedRegions= 0;
		std::vector< std::vector<int> > regionsToBeMerged;
		for(unsigned int i=0;i<regions.size();i++) regionsToBeMerged.push_back( std::vector<int>() );
		std::vector<int> regionsToBeDeleted;
		std::vector<int> regionsIdToBeDeleted;

		double Diss_min= 1.e+99;
		double Diss_max= -1.e+99;
		
		for(unsigned int i=0;i<sort_index.size();i++){
			size_t index= sort_index[i];//region index in regions list
			double thisRank= ranks[index];//region rank
			int regionId= regions[index]->fId;
			int mapIndex= regionIdMap[regionId];
			//cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Region index="<<index<<" mapIndex="<<mapIndex<<" regionId="<<regionId<<" rank="<<thisRank<<endl;

			//Check if this seed was not already merged by a previous (best ranked) region
			std::vector<int>::iterator seedfinderIt= std::find(regionsIdToBeDeleted.begin(),regionsIdToBeDeleted.end(),regionId);
			if(seedfinderIt!=regionsIdToBeDeleted.end()){
				cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Seed ranked region (id="<<regionId<<") was already selected for merging in a previous node, skip this merging!"<<endl;
				continue;
			}

			//Loop over neighbours to decide which one to merge given the similarity matrix
			std::vector<Region::NeighborInfo> neighborsInfo= regions[index]->fNeighbourRegionInfo;//CHECK!!!!!!!
			if(neighborsInfo.size()<=0) continue;
			
			//cout<<"Region id="<<regionId<<" neighbors info(";
			//for(unsigned int ss=0;ss<neighborsInfo.size();ss++) cout<<neighborsInfo[ss].D+regFactor*neighborsInfo[ss].Edgeness<<", ";
			//cout<<endl;

			int closerNeighborId= neighborsInfo[0].id;
			int closerNeighborIndex= regionIdMap[closerNeighborId];	
			double closerNeighborDist= neighborsInfo[0].D;	
			double closerNeighborEdgeness= neighborsInfo[0].Edgeness;
			double closerNeighborEntropy= neighborsInfo[0].H;
			//double Delta_ij= DissimilarityMatrix(index,closerNeighborIndex);//should be equal to closerNeighborDist without edgeness
			//double Delta_ji= DissimilarityMatrix(closerNeighborIndex,index);
			double Delta_ij= (*DissimilarityMatrix)(index,closerNeighborIndex);//should be equal to closerNeighborDist without edgeness
			double Delta_ji= (*DissimilarityMatrix)(closerNeighborIndex,index);
			double AbsDelta_ij= (*AbsDissimilarityMatrix)(index,closerNeighborIndex);//should be equal to closerNeighborDist without edgeness
			double AbsDelta_ji= (*AbsDissimilarityMatrix)(closerNeighborIndex,index);
			
			//Get neighbors of the closer region
			std::vector<Region::NeighborInfo> closerNeighborsNNs= regions[closerNeighborIndex]->fNeighbourRegionInfo;
			std::vector<Region::NeighborInfo> closerNeighborsGoodNNs;
			closerNeighborsGoodNNs.assign(closerNeighborsNNs.begin(),closerNeighborsNNs.begin()+ceil(closerNeighborsNNs.size()/2.));
				
			//Find current region among best neighbors of closer
			std::vector<Region::NeighborInfo>::iterator it = std::find_if(closerNeighborsGoodNNs.begin(), closerNeighborsGoodNNs.end(), Region::FindNeighborId(regionId));

			if(it==closerNeighborsGoodNNs.end()) {
				cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Neighbor region id="<<closerNeighborId<<" rejected for merging with region "<<regionId<<" as cannot be found within its closer neighbors..."<<endl;
				continue;
			}

			double dissTolFactor= 1.15; 
			if(Delta_ji>Delta_ij*dissTolFactor ) {
				cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Neighbor region id="<<closerNeighborId<<" rejected for merging with region "<<regionId<<"(Delta_ji="<<Delta_ji<<">1.15*Delta_ij="<<Delta_ij<<endl;
				continue;
			}
			
			
			//Check if this closer region was not already selected as a seed merger previously
			if( regionsToBeMerged[closerNeighborIndex].size()>0 ){	
				cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Closer neighbor (id="<<closerNeighborId<<") selected for merging was before selected as a primary merger, skip this merging!"<<endl;
				continue;
			}
			//Check if this closer region was not already selected to be merged to another node previously	
			std::vector<int>::iterator finderIt= std::find(regionsIdToBeDeleted.begin(),regionsIdToBeDeleted.end(),closerNeighborId);
			if(finderIt!=regionsIdToBeDeleted.end()){
				cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Closer neighbor (id="<<closerNeighborId<<") was already selected for merging in a previous node, skip this merging!"<<endl;
				continue;
			}

			int regionIndex_A= regionIdMap_top[regionId];
			int regionIndex_B= regionIdMap_top[closerNeighborId];
			
			//Similarity distance threshold
			double DissThreshold= fabs(AbsDissMedian0 + fSPMergingDistThreshold*AbsDissMedianRMS0);
			//double DissThreshold= fabs(AbsDissMedian + fSPMergingDistThreshold*AbsDissMedianRMS);
			if(fSPMergingUseAdaptingDistThreshold) DissThreshold+= hierarchyLevel*(AbsDissMax-AbsDissMin)/fSPMergingAdaptingDistThresholdScale;
			//DissThreshold= 1.e+99;//DEBUG!!!			

			//if(Delta_ij>HThreshold){
			if(AbsDelta_ij>DissThreshold){
				//cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Region id="<<regionId<<": closer neighbor (id="<<closerNeighborId<<") cannot be merged as dissimilarity is too large (Diss="<<Delta_ij<<">"<<HThreshold<<")"<<endl;
				cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Region id="<<regionId<<": closer neighbor (id="<<closerNeighborId<<") cannot be merged as dissimilarity is too large (AbsDiss="<<AbsDelta_ij<<">"<<DissThreshold<<")"<<endl;
				continue;
			}
		
			cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Neighbor region id="<<closerNeighborId<<" selected to be merged to region "<<regionId<<" (Delta_ij="<<Delta_ij<<", Delta_ji="<<Delta_ji<<" AbsDelta_ij="<<AbsDelta_ij<<", closerNeighborDist="<<closerNeighborDist<<")"<<endl;	
	
			//Store min & max dissimilarities at this hierarchy
			if(Delta_ij<Diss_min) Diss_min= Delta_ij;
			if(Delta_ij>Diss_max) Diss_max= Delta_ij;
			

			MergedDList.push_back(AbsDelta_ij);

			regionsToBeMerged[mapIndex].push_back(closerNeighborId);
			regionsToBeDeleted.push_back(closerNeighborIndex);
			regionsIdToBeDeleted.push_back(closerNeighborId);		
			mergedRegionList[regionIndex_A].push_back(regionIndex_B);
	
			nMergedRegions++;

			//Stop merging above nmerge threshold			
			int maxRegionsToMerge= std::round(regions.size()*mergeRatio);
			if(nMergedRegions>=maxRegionsToMerge) {
				cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Maximum number of regions mergeable reached for this hierarchy level!"<<endl;
				break;
			}
		}//end loop ranks
		

		if(nMergedRegions==0){
			cerr<<"VLSlicSegmentation::HierarchicalMergeRegions(): WARN: No regions merged in this stage, exit to avoid stuck!"<<endl;	
			break;
		}

		//Add merging info
		spMergeInfo.DissMin= AbsDissMin;
		spMergeInfo.DissMax= AbsDissMax;
		spMergeInfo.DissMedian= AbsDissMedian;
		spMergeInfo.DissMedianRMS= AbsDissMedianRMS;
		spMergeInfo.DissMedian0= AbsDissMedian0;
		spMergeInfo.DissMedianRMS0= AbsDissMedianRMS0;
		spMergeInfo.DissList= DList;

		std::sort(MergedDList.begin(),MergedDList.end());
		spMergeInfo.MergedDissList= MergedDList;
		spMergeInfo.MergedDissMin= MergedDList[0];
		spMergeInfo.MergedDissMax= MergedDList[MergedDList.size()-1];
		
		fSPMergingInfo.push_back(spMergeInfo);
		cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Merge info @ level="<<spMergeInfo.levelId<<": NR="<<spMergeInfo.NR<<", Diss min/max="<<spMergeInfo.DissMin<<"/"<<spMergeInfo.DissMax<<", DissMedian="<<spMergeInfo.DissMedian<<", DissMedianRMS="<<spMergeInfo.DissMedianRMS<<", MSE="<<MSE<<endl;

		t2 = high_resolution_clock::now();
		auto dt_selectMergingRegions= std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
		elapseTime_selectMergingRegions+= dt_selectMergingRegions;
		//cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: selectMergingRegions dt(usec)="<<dt_selectMergingRegions<<endl;
	

		//## Merge regions
		t1 = high_resolution_clock::now();
		cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Merge the selected regions..."<<endl;
		std::map<int,int> newLabelMap;

		for(unsigned int k=0;k<regions.size();k++){
			int regionId= regions[k]->fId;
			newLabelMap.insert( std::pair<int,int>(regionId,regionId) );
			for(unsigned int j=0;j<regionsToBeMerged[k].size();j++){
				int mergedRegionId= regionsToBeMerged[k][j];
				int mergedRegionIndex= regionIdMap[mergedRegionId];
				newLabelMap[mergedRegionId]= regionId;
				cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Region no. "<<k<<" (id="<<regionId<<") : merging region id="<<mergedRegionId<<" (index="<<mergedRegionIndex<<")"<<endl;
				regions[k]->AddRegion(regions[mergedRegionIndex]);
			}//end loop regions to be merged
		}//end loop regions
	
		//## Delete aggregated region from region list and index map
		cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Deleting regions aggregated in this step from the main list..."<<endl;
		Utils::DeleteItems(regions, regionsToBeDeleted);
		for(size_t k=0;k<regionsIdToBeDeleted.size();k++) regionIdMap.erase(regionsIdToBeDeleted[k]);

		//## Update map and recompute parameters & contours (will be used by nearest neighbors search)
		cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Updating region parameters & contours..."<<endl;
		//cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Surviving regions (";
		for(unsigned int k=0;k<regions.size();k++){
			regions[k]->ComputeParameters(false,useRobustParams,false);
			//regions[k]->Dump();
			int regionId= regions[k]->fId;
			regionIdMap[regionId]= k;
			//cout<<"(id="<<regionId<<",pos="<<k<<"), ";
		}//end loop regions
		//cout<<")"<<endl;

		//## Update pixel labels
		cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Updating pixel labels..."<<endl;	
		for(unsigned int i=0;i<labels.size();i++) {
			for(unsigned int j=0;j<labels[i].size();j++) {
				int oldLabel= labels[i][j];
				int newLabel= newLabelMap[oldLabel];
				labels[i][j]= newLabel;
			}
		}

		//## Update saliency map
		Img* thisLevelSaliencyMap= ComputeSaliencyMap(regions);
		if(thisLevelSaliencyMap){
			fSumSaliencyImg->Add(thisLevelSaliencyMap);
			thisLevelSaliencyMap->Delete();
		}

		t2 = high_resolution_clock::now();
		auto dt_copyRegions= std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
		elapseTime_copyRegions+= dt_copyRegions;
		//cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: copyRegions dt(usec)="<<dt_copyRegions<<endl;
	
		nMergedRegionsInHierarchyLevel= nMergedRegions;
		hierarchyLevel++;

		/*
		if(nMergedRegionsInHierarchyLevel==0){
			cerr<<"VLSlicSegmentation::HierarchicalMergeRegions(): WARN: No regions merged in this stage, exit to avoid stuck!"<<endl;	
			break;
		}
		*/
		cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: "<<nMergedRegionsInHierarchyLevel<<"/"<<regions.size()<<" regions aggregated at this level hierarchy..."<<endl;

	}//end while loop


	cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: "<<hierarchyLevel<<" hierarchy levels aggregated: N="<<regions.size()<<" regions left"<<endl;


	if(fTagSignificativeSP){
		cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Tagging significant regions..."<<endl;

		/*
		/#### DEBUG ####
		for(int i=0;i<regions.size();i++){
			int regionId= regions[i]->fId;
			std::vector<int> subregionIds= regions[i]->fSubRegionIds;
			cout<<"Region id="<<regionId<<"(";
			for(unsigned int j=0;j<subregionIds.size();j++){
				cout<<subregionIds[j]<<",";
			}
			cout<<")"<<endl;	
		}
		
		cout<<"----"<<endl;
		for(unsigned int i=0;i<mergedRegionList.size();i++){
			int regionId= fRegions[i]->fId;
			cout<<"Region id="<<regionId<<"(";
			for(unsigned int j=0;j<mergedRegionList[i].size();j++){
				int subRegionIndex= mergedRegionList[i][j];
				int subRegionId= fRegions[subRegionIndex]->fId;
				cout<<subRegionId<<",";
			}
			cout<<")"<<endl;
		}
		//#######################
		*/

		//## Final tag of composite regions according to a majority criterion, that is region tagged as bkg if the majority of sub-regions are non-significative
		for(int i=0;i<regions.size();i++){
			int regionId= regions[i]->fId;
			int regionIndex_top= regionIdMap_top[regionId];//region index at the top of the hierarchy
			std::vector<int> subregionIds= regions[i]->fSubRegionIds;

			bool isSignificant= fRegions[regionIndex_top]->fIsSignificative;
			bool isSalient= fRegions[regionIndex_top]->fIsSalient;
			double salientFraction= 0;
			int nSubRegions= 1+subregionIds.size();
			if(isSalient) salientFraction++;
			bool isAnySignificant= false;
			if(isSignificant) isAnySignificant= true;

			//cout<<"Region id="<<regionId<<"(";
			for(unsigned int j=0;j<subregionIds.size();j++){
				int subRegionId= subregionIds[j];
				int subRegionIndex_top= regionIdMap_top[subRegionId];
				bool subIsSignificant= fRegions[subRegionIndex_top]->fIsSignificative;
				bool subIsSalient= fRegions[subRegionIndex_top]->fIsSalient;	
				if(subIsSalient) salientFraction++;
				if(subIsSignificant) isAnySignificant= true;
				//cout<<subRegionId<<",";
			}
			//cout<<")"<<endl;
			salientFraction/= (double)nSubRegions;
			
			if( (salientFraction>fSignificantSPRatio) || isAnySignificant ){
				cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Region id="<<regionId<<" tagged as significant as (significanceRatio="<<salientFraction<<">"<<fSignificantSPRatio<<")..."<<endl;
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
					cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Region id="<<regionId<<" has too many pixels (ratio="<<pixelRatio<<">"<<fPixelRatioCut<<") and will be tagged as bkg..."<<endl;
					regions[i]->fIsSignificative= false;
				}
			}
		}
	}//close if fTagSignificativeSP

	//## Copy final results to global data (fRegions, fPixelClusterIds)
	fRegions.clear();
	fRegions.assign(regions.begin(),regions.end());
	for(unsigned int i=0;i<labels.size();i++) {
		for(unsigned int j=0;j<labels[i].size();j++) {
			fPixelClusterIds[i][j]= labels[i][j];
			//cout<<"Pixel("<<i<<","<<j<<") label="<<fPixelClusterIds[i][j]<<endl;
		}
	}
	
	fRegionTextList.clear();
	fRegionText= 0;
	float colors[3];
	for(unsigned int k=0;k<fRegions.size();k++){
		fRegions[k]->Dump();
		int nPix= fRegions[k]->fNPix;
		int regionId= fRegions[k]->fId;
		bool isSignificative= fRegions[k]->fIsSignificative;
		if(nPix<=0) continue;

		fRegions[k]->ComputeParameters(false,useRobustParams,true);

		double Cx= fRegions[k]->fX0;
		double Cy= fRegions[k]->fY0;
		double Cz= fRegions[k]->fMean;

		double regionMean= fRegions[k]->fMean;
		double regionMean_lg= log10(regionMean);
		int LUTBinId= fLUT->FindBin(regionMean);//use log mapping
		if(fUseLogNormalizedImage) LUTBinId= fLUT->FindBin(regionMean_lg);//use log mapping
		int colorId= TColor::GetColorPalette(LUTBinId); 
		gROOT->GetColor(colorId)->GetRGB(*colors,*(colors+1),*(colors+2));//colors are in range [0,1]		
		for(int kk=0;kk<3;kk++) colors[kk]*= 255;//convert colors to 0-255 
		(fRegions[k]->fColor).SetXYZ(colors[0],colors[1],colors[2]);

		//cout<<"VLSlicSegmentation::HierarchicalMergeRegions(): INFO: Region id "<<regionId<<", isSignificative? "<<isSignificative<<": N="<<nPix<<" C("<<Cx<<","<<Cy<<","<<Cz<<") color("<<(fRegions[k]->fColor).X()<<","<<(fRegions[k]->fColor).Y()<<","<<(fRegions[k]->fColor).Z()<<"), connection(";
		//cout<<endl;
			
		double textX= Cx;
		double textY= Cy;
		if(isSignificative){
			fRegionText= new TText(textX,textY,Form("%d",regionId));
			fRegionText->SetTextSize(0.015);
			fRegionText->SetTextColor(kBlack);
			fRegionTextList.push_back(fRegionText);		
		}
	}//end loop new regions	
	
	
	
	//## Recompute contours
	TString contourGraphName= Form("%s-HierarchMergedRegionContours",std::string(fImg->GetName()).c_str());
	TGraph* clusterContours= ComputeClusterContours(fImg,fPixelClusterIds,fRegions);
	clusterContours->SetNameTitle(contourGraphName,contourGraphName);

	//## Recolor image with cluster means	
	TString imgName= Form("%s-HierarchMergedRegionImg",std::string(fImg->GetName()).c_str());
	//Img* ColoredImg= GetClusterColoredImage(fImg,fRegions,selectSignificativeRegions);
	Img* ColoredImg= GetClusterColoredImage(fImg,fRegions,fTagSignificativeSP);
	ColoredImg->SetNameTitle(imgName,imgName);

	TString canvasName= Form("%s-HierarchMergedRegionImgPlot",std::string(fImg->GetName()).c_str());
	TCanvas* MergeRegionColoredPlot= new TCanvas(canvasName,canvasName);
	MergeRegionColoredPlot->cd();
		
	if(ColoredImg) ColoredImg->Draw("COLZ");
	for(unsigned int k=0;k<fRegionTextList.size();k++) fRegionTextList[k]->Draw("same");
	if(clusterContours) clusterContours->Draw("Psame");
		

	return 0;

}//close HierarchicalMergeRegions()

std::vector<double> VLSlicSegmentation::ComputePageRank(TMatrixD M,double d,double tol){

	bool hasConverged= true;
	int maxIterToStop= 100;
	int N= M.GetNcols();	
	
	//Initialize rank with uniform prob (0.5)
	TVectorD v(N);
	TVectorD v_last(N);
	TMatrixD UnoMatrix(N,N);	
	UnoMatrix.Zero();
	for(int i=0;i<N;i++) {
		v(i)= 0.5;
		v_last(i)= 1.e+99;//TMath::Infinity();
		for(int j=0;j<N;j++) UnoMatrix(i,j)= 1;
	}
	
	double norm2= v.Norm2Sqr();
	v*= 1./sqrt(norm2);
	
	//Compute Mhat= (d x M) + (1-d)/N*UnoMatrix(N,N)
	TMatrixD M_hat(N,N);
	M_hat.Zero();
	M_hat= d*M + ((1.-d)/(double)(N))*UnoMatrix;
	
	//cout<<"== M_hat =="<<endl;	
	//M_hat.Print();
	//cout<<"==========="<<endl;

	int iter= 0;
	int nIterToConverge= 0;
	double diff= (v-v_last).Norm2Sqr();
	//cout<<"--> Iter no. "<<iter<<" diff="<<diff<<" tol="<<tol<<endl;
	//v.Print();

	//while( (v-v_last).Norm2Sqr()>tol ){
	while( diff>tol ){
	//while (diff > tol && iter < maxIterToStop) {
  	v_last= v;
    v = M_hat*v;
		double normFactor= sqrt(v.Norm2Sqr());
    v*= 1./normFactor;
		
		diff= sqrt((v-v_last).Norm2Sqr());
		//cout<<"--> Iter no. "<<iter<<" diff="<<diff<<" tol="<<tol<<endl;
		//v.Print();

		iter++;
		nIterToConverge= iter;

    if(iter >= maxIterToStop){
    	hasConverged= false;
			break;
		}
	}//end loop 
	
	std::vector<double> ranks;
	ranks.clear();

	if(hasConverged){
		cout<<"== GOOGLE RANKS =="<<endl;
		cout<<"nIter="<<nIterToConverge<<" p(";
		for(int k=0;k<v.GetNoElements();k++) {
			ranks.push_back(v[k]);
			cout<<v[k]<<",";
		}
		cout<<")"<<endl;
		cout<<"=================="<<endl;
	}
	else{
		cerr<<"ComputePageRank(): WARN: Page rank did not converge!"<<endl;
	}

	return ranks;

}//close VLSlicSegmentation::ComputePageRank()


void VLSlicSegmentation::findNeighbors_v2(std::vector<Region*> regions,std::map<int,int> mapping,std::vector< std::vector<long int> > labels,bool add2ndOrderNeighbors,bool useRobustParams) {
  	
	std::chrono::high_resolution_clock::time_point t1, t2, tstart, tstop;
	tstart = high_resolution_clock::now();
		
	double elapseTime_computeParams= 0;
	double elapseTime_computeEdges= 0;
	double elapseTime_add1stNeighbors= 0;
	double elapseTime_add2ndNeighbors= 0;
	double elapseTime_sortNeighbors= 0;
	
	//## Loop over region list to find neighbors
	int Nx= fInputImg->GetNbinsX();
	int Ny= fInputImg->GetNbinsY();
	int nRegions= (int)regions.size();
	cout<<"VLSlicSegmentation::findNeighbors(): INFO: Finding neighbors (NR="<<nRegions<<")"<<endl;
	
	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};
	
	// Initialize the contour vector and the matrix detailing whether a pixel is already taken to be a contour
	vec2db istaken;
	vec2db isNeighborLinkTaken;
	Contour* contour= 0;
	std::vector<Contour*> contourCollection;
	std::vector< std::vector<int> > connectedRegionIds;
	std::vector< std::map<int,std::vector<int>> > connectedRegionSharedBoundaryPixelIds;
	for(int k=0;k<nRegions;k++){
		contour= new Contour;
		contourCollection.push_back(contour);
		connectedRegionIds.push_back( std::vector<int>() );
		connectedRegionSharedBoundaryPixelIds.push_back( std::map<int,std::vector<int>>() );
		isNeighborLinkTaken.push_back( std::vector<bool>() );
		for(int s=0;s<nRegions;s++) isNeighborLinkTaken[k].push_back(false);
	}
	
	for (int i=0; i<Nx; i++) { 
  	std::vector<bool> nb;
    for (int j = 0; j < Ny; j++) {
   		nb.push_back(false);
    }
    istaken.push_back(nb);
  }//end loop Nx
  

	//## Clear existing neighbor lists and compute parameters for all regions
	for(int k=0;k<nRegions;k++){
		//Clear list of neighbor info for this region
		(regions[k]->fNeighbourRegions).clear();
		(regions[k]->fNeighbourRegionInfo).clear();

		//Compute region parameters & contour info	
		t1 = high_resolution_clock::now();
		regions[k]->ComputeParameters(false,useRobustParams,false);
		t2 = high_resolution_clock::now();
		auto dt_computeParams= std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
		elapseTime_computeParams+= dt_computeParams;
  }//end loop regions


  //## Go through all the pixels and identify contours and neighbors
	t1 = high_resolution_clock::now();
  for (int i=0; i<Nx; i++) {
  	for (int j=0; j<Ny; j++) {
			int regionId= labels[i][j];
			if(regionId<0) continue;
			int regionIndex= mapping[regionId];			
   		int nr_p = 0;
			std::map<int,int> connectedRegionCounterMap;
			connectedRegionCounterMap.clear();
			std::map<int,std::vector<int>> connectedRegionBoundaryPixelListMap;
			connectedRegionBoundaryPixelListMap.clear();
            
      // Compare the pixel to its 8 neighbours
      for (int k=0; k<8; k++) {
				int x = i + dx8[k];
				int y = j + dy8[k];
				         
        if (x >= 0 && x < Nx && y >= 0 && y < Ny) {
					int regionId_neighbor= labels[x][y];
       	
        	//if (!istaken[x][y] && regionId!=regionId_neighbor && regionId_neighbor>=0) {	
					if (regionId!=regionId_neighbor && regionId_neighbor>=0) {	
          	nr_p++;
						connectedRegionCounterMap[regionId_neighbor]++;
						int gBin_neighbor= fInputImg->GetBin(x+1,y+1);
						connectedRegionBoundaryPixelListMap[regionId_neighbor];
						connectedRegionBoundaryPixelListMap[regionId_neighbor].push_back(gBin_neighbor);
          }
        }
      }//end loop neighbours
            
      // Add the pixel to the contour list of corresponding region
     	if (nr_p>= 2) {
				int gBin= fInputImg->GetBin(i+1,j+1);
				double cx= fInputImg->GetXaxis()->GetBinCenter(i+1);
				double cy= fInputImg->GetYaxis()->GetBinCenter(j+1);
				double Sedge= fEdgeFilterImg->GetBinContent(i+1,j+1);
				contourCollection[regionIndex]->AddPoint(cv::Point2f(cx,cy));
        istaken[i][j] = true;

				/*
				for (std::map<int,int>::iterator counterIt=connectedRegionCounterMap.begin(); counterIt!=connectedRegionCounterMap.end(); ++counterIt){
					int neighborId= counterIt->first;
					int counts= counterIt->second;	
					int neighbourIndex= mapping[neighborId];
					if(counts>=2) {	
						std::vector<int>::iterator it = std::find (connectedRegionIds[regionIndex].begin(), connectedRegionIds[regionIndex].end(), neighborId);
						if( connectedRegionIds[regionIndex].empty() || it==connectedRegionIds[regionIndex].end() ) {
							connectedRegionIds[regionIndex].push_back(neighborId);
						}//close if
					}//close if counts>2
				}//end loop map counter iterator
				*/

				//Add neighbor ids and boundary pixel ids
				std::map<int,std::vector<int>>::iterator counterListIt = connectedRegionBoundaryPixelListMap.begin();
				for (counterListIt=connectedRegionBoundaryPixelListMap.begin(); counterListIt!=connectedRegionBoundaryPixelListMap.end(); ++counterListIt){
					int neighborId= counterListIt->first;
					std::vector<int> neighborPixIds= counterListIt->second;
					int counts= (int)neighborPixIds.size();	
					int neighbourIndex= mapping[neighborId];

					if(counts>=2) {
						std::vector<int>::iterator it = std::find (connectedRegionIds[regionIndex].begin(), connectedRegionIds[regionIndex].end(), neighborId);
						if( connectedRegionIds[regionIndex].empty() || it==connectedRegionIds[regionIndex].end() ) {
							connectedRegionIds[regionIndex].push_back(neighborId);
						}//close if

						(connectedRegionSharedBoundaryPixelIds[regionIndex])[neighborId];//add connection with neighbor if not existing in the map
						((connectedRegionSharedBoundaryPixelIds[regionIndex])[neighborId]).push_back(gBin);//add this contour pixel
						for(int t=0;t<counts;t++){
							int gBin_neighbor= neighborPixIds[t];
							it = std::find ( ((connectedRegionSharedBoundaryPixelIds[regionIndex])[neighborId]).begin(), ((connectedRegionSharedBoundaryPixelIds[regionIndex])[neighborId]).end(), gBin_neighbor);
							if( ((connectedRegionSharedBoundaryPixelIds[regionIndex])[neighborId]).empty() || it==((connectedRegionSharedBoundaryPixelIds[regionIndex])[neighborId]).end() ) {
								((connectedRegionSharedBoundaryPixelIds[regionIndex])[neighborId]).push_back(gBin_neighbor);
							}//close if
						}
					}//close if counts>2
				}//end loop map counter iterator

      }//close if nr_p>2

    }//end loop image Ny
  }//end loop image Nx
	t2 = high_resolution_clock::now();
	auto dt_computeEdges= std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
	elapseTime_computeEdges+= dt_computeEdges;

	/*
	bool hasBoundaryImg= false;
	if(gDirectory->FindObject("BoundaryImg")) hasBoundaryImg= true;
	if(!hasBoundaryImg){
		Img* boundaryImg= (Img*)fInputImg->Clone("BoundaryImg");
		boundaryImg->Reset();
		for(int k=0;k<nRegions;k++){
			//#### PRINT FOR DEBUG ####
			std::map<int,std::vector<int>> sharedPixelsMap= connectedRegionSharedBoundaryPixelIds[k];
			cout<<"VLSlicSegmentation::findNeighbors(): INFO: Region id="<<regions[k]->fId<<": "<<connectedRegionIds[k].size()<<" neighbors to be added at 1st-order (";
		
			for(unsigned int l=0;l<connectedRegionIds[k].size();l++){
				int neighborId= connectedRegionIds[k][l];
				int neighbourIndex= mapping[neighborId];
				cout<<neighborId<<"{";
				std::vector<int> sharedPixelIds= sharedPixelsMap[neighborId];

				for(unsigned int t=0;t<sharedPixelIds.size();t++) {
					int gBin= sharedPixelIds[t];
					cout<<gBin<<",";
					boundaryImg->SetBinContent(gBin,1);
				}
				cout<<"},";
			}//end loop neighbors
			cout<<")"<<endl;
		}//end loop regions
		TCanvas* BoundaryPlot= new TCanvas("BoundaryPlot","BoundaryPlot");
		BoundaryPlot->cd();
		boundaryImg->Draw("COLZ");
	}//close if
	*/

	double NormMin= 0;
	double NormMax= 1;
	double EdgeImgMin= fEdgeFilterImgStats->min;
	double EdgeImgMax= fEdgeFilterImgStats->max;
	
	//Add 1st neighbors to the list
	for(int k=0;k<nRegions;k++){
		t1 = high_resolution_clock::now();
	
		int regionId= regions[k]->fId;
		
		regions[k]->fContour= contourCollection[k];	
		regions[k]->fNeighbourRegions= connectedRegionIds[k];
		
		for(unsigned int l=0;l<connectedRegionIds[k].size();l++){
			int neighborId= connectedRegionIds[k][l];
			int neighbourIndex= mapping[neighborId];

			//Check if this link has been already computed
			if(isNeighborLinkTaken[k][neighbourIndex] && isNeighborLinkTaken[neighbourIndex][k]) continue;

			//Compute asymm distance between regions
			std::pair<double,double> dists= regions[k]->GetAsymmDistance(regions[neighbourIndex],useRobustParams);

			//Compute edgeness between regions
			double Edgeness= 0;	
			double Edgeness_norm= 0;		
			std::vector<int> sharedPixelIds= (connectedRegionSharedBoundaryPixelIds[k])[neighborId];
			int nBoundaryPixels= (int)sharedPixelIds.size();
			for(int t=0;t<nBoundaryPixels;t++) {
				int gBin= sharedPixelIds[t];
				double S_edge= fEdgeFilterImg->GetBinContent(gBin);
				double S_edge_norm= NormMin + (NormMax-NormMin)*(S_edge-EdgeImgMin)/(EdgeImgMax-EdgeImgMin);
				Edgeness+= S_edge;
				Edgeness_norm+= S_edge_norm;
			}	
			if(nBoundaryPixels>0) {
				Edgeness/= (double)nBoundaryPixels; 
				Edgeness_norm/= (double)nBoundaryPixels; 
			}


			/*
			//Compute Entropy pars
			std::vector<double> entropyPars= regions[k]->GetEntropyDistance(regions[neighbourIndex]);
			double H1= entropyPars[0];
			double H2= entropyPars[1];
			double CE12= entropyPars[2];
			double CE21= entropyPars[3];
			double entropyCost= 1 - 0.5*(H1/CE12 + H2/CE21);//nearly 0 when the two regions are similar (H1=CE12, H2=CE21)
			cout<<"VLSlicSegmentation::findNeighbors(): INFO: Region link ("<<regionId<<","<<neighborId<<") H1="<<H1<<" CE12="<<CE12<<" H2="<<H2<<" CE21="<<CE21<<" HCost="<<entropyCost<<endl;
			*/
			double entropyCost= 0;

			//Fill neighbor info struct
			Region::NeighborInfo info;
			info.id= neighborId;	
			info.D= dists.first;
			info.Edgeness= Edgeness_norm;//Edgeness;
			info.order= 1;
			info.H= entropyCost;

			Region::NeighborInfo info_neighbor;
			info_neighbor.id= regionId;	
			info_neighbor.D= dists.second;
			info_neighbor.Edgeness= Edgeness_norm;//Edgeness;
			info_neighbor.order= 1;
			info_neighbor.H= entropyCost;
			
			//regions[k]->AddNeighborInfo(regions[neighbourIndex],1,useRobustParams);
			regions[k]->AddNeighborInfo(info);
			regions[neighbourIndex]->AddNeighborInfo(info_neighbor);
			isNeighborLinkTaken[k][neighbourIndex]= true;
			isNeighborLinkTaken[neighbourIndex][k]= true;
			
		}//end loop neighbors
		
		t2 = high_resolution_clock::now();
		auto dt_add1stNeighbors= std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
		elapseTime_add1stNeighbors+= dt_add1stNeighbors;
	}//end loop regions



	//Add 2nd order neighbors
	if(add2ndOrderNeighbors){
		t1 = high_resolution_clock::now();
		
		std::vector< std::vector<int> > NeighborIdsToBeAdded;

		for(int k=0;k<nRegions;k++){
			NeighborIdsToBeAdded.push_back( std::vector<int>() );
			int regionId= regions[k]->fId;
			
			//Get list of 1st neighbors
			std::vector<Region::NeighborInfo> neighbors= regions[k]->fNeighbourRegionInfo;

			for(unsigned int s=0;s<neighbors.size();s++){
				int neighborId= neighbors[s].id;
				int neighborIndex= mapping[neighborId];
				
				//Get list of 2nd neighbors
				std::vector<Region::NeighborInfo> neighbors_2nd= regions[neighborIndex]->fNeighbourRegionInfo;
				for(unsigned int t=0;t<neighbors_2nd.size();t++){
					int neighborId_2nd= neighbors_2nd[t].id;
					int neighborIndex_2nd= mapping[neighborId_2nd];
					//cout<<"neighborId="<<neighborId<<": adding 2nd neighbor id "<<neighborId_2nd<<"("<<regions[neighborIndex_2nd]->fId<<") to region "<<regionId<<"?"<<endl;
					
					//Find if neighbors was already added in the lists
					std::vector<int>::iterator it= std::find(NeighborIdsToBeAdded[k].begin(),NeighborIdsToBeAdded[k].end(),neighborId_2nd);
					
					if(neighborId_2nd==regionId || it!=NeighborIdsToBeAdded[k].end()) continue; 
					NeighborIdsToBeAdded[k].push_back(neighborId_2nd);
				}//end loop 2nd neighbors
			}//end loop 1st neighbors
		}//end loop regions

		std::chrono::high_resolution_clock::time_point t3= high_resolution_clock::now();
		auto dt_add2ndNeighbors_firstStep= std::chrono::duration_cast<std::chrono::microseconds>(t3-t1).count();

		//Add 2nd neighbors at the end
		for(int k=0;k<nRegions;k++){	
			int regionId= regions[k]->fId;

			//cout<<"VLSlicSegmentation::findNeighbors(): INFO: Region id="<<regions[k]->fId<<": "<<NeighborIdsToBeAdded[k].size()<<" neighbors to be added at 2nd-order..."<<endl;
			for(unsigned int j=0;j<NeighborIdsToBeAdded[k].size();j++) {
				int neighborId_2nd= NeighborIdsToBeAdded[k][j];
				int neighborIndex_2nd= mapping[neighborId_2nd];
				
				//Check if this link has been already computed
				if(isNeighborLinkTaken[k][neighborIndex_2nd] && isNeighborLinkTaken[neighborIndex_2nd][k]) continue;

				//Compute asymm distance between regions
				std::pair<double,double> dists= regions[k]->GetAsymmDistance(regions[neighborIndex_2nd],useRobustParams);

				//Compute Entropy pars
				double entropyCost= 0;
				/*
				std::vector<double> entropyPars= regions[k]->GetEntropyDistance(regions[neighborIndex_2nd]);
				double H1= entropyPars[0];
				double H2= entropyPars[1];
				double CE12= entropyPars[2];
				double CE21= entropyPars[3];
				double entropyCost= 1 - 0.5*(H1/CE12 + H2/CE21);//nearly 0 when the two regions are similar (H1=CE12, H2=CE21)
				cout<<"VLSlicSegmentation::findNeighbors(): INFO: Region link ("<<regionId<<","<<neighborId_2nd<<") HCost="<<entropyCost<<endl;
				*/

				//Fill neighbor info struct
				Region::NeighborInfo info;
				info.id= neighborId_2nd;	
				info.D= dists.first;
				info.Edgeness= NormMax;//1
				info.order= 2;
				info.H= entropyCost;

				Region::NeighborInfo info_neighbor;
				info_neighbor.id= regionId;	
				info_neighbor.D= dists.second;
				info_neighbor.Edgeness= NormMax;//1
				info_neighbor.order= 1;
				info_neighbor.H= entropyCost;
			
				//regions[k]->AddNeighborInfo(regions[neighborIndex_2nd],2,useRobustParams);

				regions[k]->AddNeighborInfo(info);
				regions[neighborIndex_2nd]->AddNeighborInfo(info_neighbor);
				isNeighborLinkTaken[k][neighborIndex_2nd]= true;
				isNeighborLinkTaken[neighborIndex_2nd][k]= true;

			}//end loop 2nd order neighbors
		}//end loop regions
		t2 = high_resolution_clock::now();
		auto dt_add2ndNeighbors= std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
		elapseTime_add2ndNeighbors+= dt_add2ndNeighbors;
		
		auto dt_add2ndNeighbors_secondStep= std::chrono::duration_cast<std::chrono::microseconds>(t2-t3).count();
		//cout<<"VLSlicSegmentation::findNeighbors(): INFO: t_add2ndNeighbors="<<dt_add2ndNeighbors_firstStep<<", "<<dt_add2ndNeighbors_secondStep<<endl;
	}//close if add2ndOrderNeighbors
	

	
	//## Sort region neighbors
	t1 = high_resolution_clock::now();	
	/*
	for(int k=0;k<nRegions;k++){		
		regions[k]->SortNeighborInfo(true,false);//ascending order	
		//regions[k]->SortNeighborInfo(true,true);//ascending order	& tot dist (D+Edgeness)
	}//end loop regions	
	*/
	t2 = high_resolution_clock::now();
	auto dt_sortNeighbors= std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
	elapseTime_sortNeighbors+= dt_sortNeighbors;
	
	tstop = high_resolution_clock::now();
	auto dt_tot= std::chrono::duration_cast<std::chrono::microseconds>(tstop-tstart).count();
				
	//cout<<"VLSlicSegmentation::findNeighbors(): INFO: t_params="<<elapseTime_computeParams/dt_tot*100.<<" t_edges="<<elapseTime_computeEdges/dt_tot*100.<<" t_add1stNeighbors="<<elapseTime_add1stNeighbors/dt_tot*100.<<" t_add2ndNeighbors="<<elapseTime_add2ndNeighbors/dt_tot*100.<<" t_sort="<<elapseTime_sortNeighbors/dt_tot*100.<<endl;
	

}//close VLSlicSegmentation::findNeighbors_v2()






void VLSlicSegmentation::findNeighbors(std::vector<Region*> regions,std::map<int,int> mapping,std::vector< std::vector<long int> > labels,bool add2ndOrderNeighbors,bool useRobustParams){

	std::chrono::high_resolution_clock::time_point t1, t2, tstart, tstop;
	tstart = high_resolution_clock::now();
		
	double elapseTime_computeParams= 0;
	double elapseTime_computeEdges= 0;
	double elapseTime_add1stNeighbors= 0;
	double elapseTime_add2ndNeighbors= 0;
	double elapseTime_sortNeighbors= 0;
	
	//## Loop over region list to find neighbors
	int nRegions= (int)regions.size();
	cout<<"VLSlicSegmentation::findNeighbors(): INFO: Finding neighbors (NR="<<nRegions<<")"<<endl;
		
	std::vector< std::vector<int> > connectedRegionIds;
	connectedRegionIds.clear();		

	for(int k=0;k<nRegions;k++){
		//Clear list of neighbor info for this region
		(regions[k]->fNeighbourRegions).clear();
		(regions[k]->fNeighbourRegionInfo).clear();

		//Compute region parameters & contour info	
		t1 = high_resolution_clock::now();
		regions[k]->ComputeParameters(true,useRobustParams,false);
		t2 = high_resolution_clock::now();
		auto dt_computeParams= std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
		elapseTime_computeParams+= dt_computeParams;

		//Get points around contour
		t1 = high_resolution_clock::now();
		std::vector<Region::EdgeInfo> edges= regions[k]->fAroundContourPoints;	
		int binX, binY, binZ;
		std::map<int,int> connectedRegionCounterMap;
		connectedRegionCounterMap.clear();

		connectedRegionIds.push_back( std::vector<int>() );
		
		for(unsigned l=0;l<edges.size();l++){
			int gBin= fImg->FindBin(edges[l].x,edges[l].y);
			fImg->GetBinXYZ(gBin,binX,binY,binZ);
			if( fImg->IsBinOverflow(gBin) || fImg->IsBinUnderflow(gBin) ) continue;
			int ix= binX-1;
			int iy= binY-1;
			int pixelLabel= labels[ix][iy];
			double S_curv= fLaplImg->GetBinContent(gBin);
			double S_edge= fEdgeFilterImg->GetBinContent(gBin);
				
			//Increment region counter in map (if label does not exist, a new label is created and init to zero)
			if(pixelLabel>=0) {
				connectedRegionCounterMap[pixelLabel]++;
				((regions[k]->fAroundContourPoints)[l]).labelId= pixelLabel;//add region label to edge info
				((regions[k]->fAroundContourPoints)[l]).S_edge= S_edge;
				((regions[k]->fAroundContourPoints)[l]).S_curv= S_curv;
			}
		}//end loop edge points

		//Add 1st neighbors to the list for this region
		for (std::map<int,int>::iterator counterIt=connectedRegionCounterMap.begin(); counterIt!=connectedRegionCounterMap.end(); ++counterIt){
			int neighborId= counterIt->first;
			int counts= counterIt->second;	
			int neighbourIndex= mapping[neighborId];
			if(counts>=2) {	
				connectedRegionIds[k].push_back(neighborId);
			}
		}//end loop map counter iterator
		
		t2 = high_resolution_clock::now();
		auto dt_computeEdges= std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
		elapseTime_computeEdges+= dt_computeEdges;
		
	}//end loop regions

	cout<<"VLSlicSegmentation::findNeighbors(): INFO: <t_params>="<<elapseTime_computeParams/(double)(nRegions)<<endl;


	//Add 1st neighbors to the list
	for(int k=0;k<nRegions;k++){
		t1 = high_resolution_clock::now();
	
		for(unsigned int l=0;l<connectedRegionIds[k].size();l++){
			int neighborId= connectedRegionIds[k][l];
			int neighbourIndex= mapping[neighborId];
			regions[k]->AddNeighborInfo(regions[neighbourIndex],1,useRobustParams);
		}//end loop neighbors

		t2 = high_resolution_clock::now();
		auto dt_add1stNeighbors= std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
		elapseTime_add1stNeighbors+= dt_add1stNeighbors;
	}//end loop regions
		


	/*
	for(int k=0;k<nRegions;k++){
		int regionId= regions[k]->fId;

		//Clear list of neighbor info for this region
		(regions[k]->fNeighbourRegions).clear();
		(regions[k]->fNeighbourRegionInfo).clear();
	
		//Compute region parameters & contour info	
		t1 = high_resolution_clock::now();
		regions[k]->ComputeParameters(true);
		t2 = high_resolution_clock::now();
		auto dt_computeParams= std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
		elapseTime_computeParams+= dt_computeParams;
		
		//Get points around contour
		t1 = high_resolution_clock::now();
		//std::vector<cv::Point2f> edges= regions[k]->fAroundContourPoints;
		std::vector<Region::EdgeInfo> edges= regions[k]->fAroundContourPoints;	
		int binX, binY, binZ;
		std::map<int,int> connectedRegionCounterMap;
		connectedRegionCounterMap.clear();

		for(unsigned l=0;l<edges.size();l++){
			int gBin= fImg->FindBin(edges[l].x,edges[l].y);
			fImg->GetBinXYZ(gBin,binX,binY,binZ);
			if( fImg->IsBinOverflow(gBin) || fImg->IsBinUnderflow(gBin) ) continue;
			int ix= binX-1;
			int iy= binY-1;
			int pixelLabel= labels[ix][iy];
			double S_curv= fLaplImg->GetBinContent(gBin);
			double S_edge= fEdgeFilterImg->GetBinContent(gBin);
				
			//Increment region counter in map (if label does not exist, a new label is created and init to zero)
			if(pixelLabel>=0) {
				connectedRegionCounterMap[pixelLabel]++;
				((regions[k]->fAroundContourPoints)[l]).labelId= pixelLabel;//add region label to edge info
				((regions[k]->fAroundContourPoints)[l]).S_edge= S_edge;
				((regions[k]->fAroundContourPoints)[l]).S_curv= S_curv;
			}
		}//end loop edge points

		for (std::map<int,int>::iterator counterIt=connectedRegionCounterMap.begin(); counterIt!=connectedRegionCounterMap.end(); ++counterIt){
			int neighborId= counterIt->first;
			int counts= counterIt->second;	
			int neighbourIndex= mapping[neighborId];
			if(counts>=2) {
				//regions[k]->AddNeighbor(neighborId);
				regions[k]->AddNeighborInfo(regions[neighbourIndex]);
			}
		}//end loop map counter iterator
		
		t2 = high_resolution_clock::now();
		auto dt_computeEdges= std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
		elapseTime_computeEdges+= dt_computeEdges;
		

		//Sort neighbor info
		//regions[k]->SortNeighborInfo(true);//ascending order

		//cout<<"Region id="<<regionId<<" Neighbors(";
		//for(unsigned int l=0;l<(regions[k]->fNeighbourRegionInfo).size();l++) {
		//	Region::NeighborInfo info= (regions[k]->fNeighbourRegionInfo)[l];
		//	cout<<"(id="<<info.id<<",D="<<info.D<<"), ";
		//}
		//cout<<")"<<endl;		
	}//end loop regions
	*/


	//Add 2nd order neighbors
	if(add2ndOrderNeighbors){
		t1 = high_resolution_clock::now();
		
		std::vector< std::vector<int> > NeighborIdsToBeAdded;

		for(int k=0;k<nRegions;k++){
			NeighborIdsToBeAdded.push_back( std::vector<int>() );
			int regionId= regions[k]->fId;
			
			//Get list of 1st neighbors
			std::vector<Region::NeighborInfo> neighbors= regions[k]->fNeighbourRegionInfo;

			for(unsigned int s=0;s<neighbors.size();s++){
				int neighborId= neighbors[s].id;
				int neighborIndex= mapping[neighborId];
				
				//Get list of 2nd neighbors
				std::vector<Region::NeighborInfo> neighbors_2nd= regions[neighborIndex]->fNeighbourRegionInfo;
				for(unsigned int t=0;t<neighbors_2nd.size();t++){
					int neighborId_2nd= neighbors_2nd[t].id;
					int neighborIndex_2nd= mapping[neighborId_2nd];
					//cout<<"neighborId="<<neighborId<<": adding 2nd neighbor id "<<neighborId_2nd<<"("<<regions[neighborIndex_2nd]->fId<<") to region "<<regionId<<"?"<<endl;
					
					//Find if neighbors was already added in the lists
					std::vector<int>::iterator it= std::find(NeighborIdsToBeAdded[k].begin(),NeighborIdsToBeAdded[k].end(),neighborId_2nd);
					
					if(neighborId_2nd==regionId || it!=NeighborIdsToBeAdded[k].end()) continue; 
					NeighborIdsToBeAdded[k].push_back(neighborId_2nd);
				}//end loop 2nd neighbors
			}//end loop 1st neighbors
		}//end loop regions

		std::chrono::high_resolution_clock::time_point t3= high_resolution_clock::now();
		auto dt_add2ndNeighbors_firstStep= std::chrono::duration_cast<std::chrono::microseconds>(t3-t1).count();

		//Add 2nd neighbors at the end
		for(int k=0;k<nRegions;k++){
			cout<<"VLSlicSegmentation::findNeighbors(): INFO: Region id="<<regions[k]->fId<<": "<<NeighborIdsToBeAdded[k].size()<<" neighbors to be added as 2nd level"<<endl;
			for(unsigned int j=0;j<NeighborIdsToBeAdded[k].size();j++) {
				int neighborId_2nd= NeighborIdsToBeAdded[k][j];
				int neighborIndex_2nd= mapping[neighborId_2nd];
				regions[k]->AddNeighborInfo(regions[neighborIndex_2nd],2,useRobustParams);
			}
		}
		t2 = high_resolution_clock::now();
		auto dt_add2ndNeighbors= std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
		elapseTime_add2ndNeighbors+= dt_add2ndNeighbors;
		
		auto dt_add2ndNeighbors_secondStep= std::chrono::duration_cast<std::chrono::microseconds>(t2-t3).count();
		cout<<"VLSlicSegmentation::findNeighbors(): INFO: t_add2ndNeighbors="<<dt_add2ndNeighbors_firstStep<<", "<<dt_add2ndNeighbors_secondStep<<endl;
	}//close if add2ndOrderNeighbors
	

	//## Sort region neighbors
	t1 = high_resolution_clock::now();
		
	for(int k=0;k<nRegions;k++){		
		regions[k]->SortNeighborInfo(true);//ascending order
		/*
		int regionId= regions[k]->fId;
		cout<<"Region id="<<regionId<<" Neighbors(";
		for(unsigned int l=0;l<(regions[k]->fNeighbourRegionInfo).size();l++) {
			Region::NeighborInfo info= (regions[k]->fNeighbourRegionInfo)[l];
			cout<<"(id="<<info.id<<",D="<<info.D<<"), ";
		}
		cout<<")"<<endl;
		*/		
	}//end loop regions
	t2 = high_resolution_clock::now();
	auto dt_sortNeighbors= std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
	elapseTime_sortNeighbors+= dt_sortNeighbors;
	
	tstop = high_resolution_clock::now();
	auto dt_tot= std::chrono::duration_cast<std::chrono::microseconds>(tstop-tstart).count();
				
	cout<<"VLSlicSegmentation::findNeighbors(): INFO: t_params="<<elapseTime_computeParams/dt_tot*100.<<" t_edges="<<elapseTime_computeEdges/dt_tot*100.<<" t_add1stNeighbors="<<elapseTime_add1stNeighbors/dt_tot*100.<<" t_add2ndNeighbors="<<elapseTime_add2ndNeighbors/dt_tot*100.<<" t_sort="<<elapseTime_sortNeighbors/dt_tot*100.<<endl;
		
}//close findNeighbors()


std::vector<int> VLSlicSegmentation::findNeighbors(int regionId, double eps){

	std::vector<int> list_of_neighbours;	
	list_of_neighbours.clear();
	list_of_neighbours.resize(0);

	//## Check given regionId
	if( regionId<0 || regionId>=(int)fRegions.size() || fRegions.size()<=0 ) {
		cout<<"VLSlicSegmentation::findNeighbors(): ERROR: Invalid region id given ("<<regionId<<") or empty region list (size="<<fRegions.size()<<"), returning empty list!"<<endl;
		return list_of_neighbours;
	}

	double eps2= eps*eps;

	Region* thisRegion= fRegions[regionId];
	std::vector<int> neighbours= thisRegion->fNeighbourRegions;
	int nConnectedRegions= (int)neighbours.size();

	for(int j=0;j<nConnectedRegions;j++){
		int neighbourId= neighbours[j];
		Region* thisNeighbourRegion= fRegions[neighbourId];
					
		//double dist2= thisRegion->GetDistance(thisNeighbourRegion);
		std::pair<double,double> dists= thisRegion->GetDistance(thisNeighbourRegion);
		double dist= dists.first;
		double dist2= dist*dist;	
		//cout<<"Region no. "<<regionId<<": 1st neighbour dist("<<regionId<<","<<neighbourId<<")="<<sqrt(dist2)<<endl;
		
		//if(dist<=eps) list_of_neighbours.push_back(neighbourId);
		if(dist2<=eps2 && dist2>0.0f && thisNeighbourRegion->fNPix>0) list_of_neighbours.push_back(neighbourId);
						
		/*
		//Loop over 2nd neighbours
		std::vector<int> neighbours_2nd= thisNeighbourRegion->fNeighbourRegions;

		for(unsigned int k=0;k<neighbours_2nd.size();k++){
			int neighbourId_2nd= neighbours_2nd[k];
			Region* thisNeighbourRegion_2nd= fRegions[neighbourId_2nd];
			//Region mergedRegion_2nd= *(thisRegion);
			//mergedRegion_2nd.AddRegion(thisNeighbourRegion_2nd);
			//mergedRegion_2nd.ComputeParameters();

			//Compute distance between current and merged region 	
			//dist= thisRegion->GetDistance(&mergedRegion_2nd);
			//similarity= 1./dist_2nd;

			dist2= thisRegion->GetDistance(thisNeighbourRegion_2nd);
			cout<<"Region no. "<<regionId<<": 2nd neighbour dist("<<regionId<<","<<regionId<<"+"<<neighbourId_2nd<<")="<<sqrt(dist2)<<endl;

			//if(dist<=eps) list_of_neighbours.push_back(neighbourId_2nd);
			if(dist2<=eps2 && dist2>0.0f && thisNeighbourRegion_2nd->fNPix>0) list_of_neighbours.push_back(neighbourId_2nd);
		
		}//end loop 2nd neighbours
		*/
	
	}//end loop neighbours
	
	return list_of_neighbours;

}//close VLSlicSegmentation::findNeighbors()



Img* VLSlicSegmentation::GetClusterColoredImage(Img* image,std::vector<Region*> regions,bool drawOnlySignificative,bool colorWithSaliency) {

	//## Check input image
	if(!image) return 0;
	
	int Nx= image->GetNbinsX();
	int Ny= image->GetNbinsY();
	int nRegions= (int)regions.size();
	
	if(nRegions<=0){
		cerr<<"VLSlicSegmentation::GetClusterColoredImage(): ERROR: Empty region list given...nothing will be done!"<<endl;
		return 0;
	}

	Img* colored_img= (Img*)image->Clone("colored_img");
	colored_img->Reset();

	double NormMin= image->GetMinimum();
	double NormMax= image->GetMaximum();
	double A= NormMin - (NormMax-NormMin)*fImgStats->min/(fImgStats->max-fImgStats->min);
	double B= (NormMax-NormMin)/(fImgStats->max-fImgStats->min);

  for(int i=0;i<nRegions;i++){//loop on regions
		int regionId= regions[i]->fId;
		bool isSignificative= regions[i]->fIsSignificative;
		
		if(drawOnlySignificative && !isSignificative) continue;
		int nPixelsInRegion= (int)regions[i]->fNPix;

		double Mean= regions[i]->fMean;
		double Mean_norm= A+B*Mean; 
		double LgMean= log10(Mean);
		double LgMean_norm= NormMin + log10(Mean_norm/NormMin)/log10(NormMax/NormMin) * (NormMax-NormMin);
		double saliency= regions[i]->fSaliency;

		double colorValue= Mean_norm;
		if(colorWithSaliency) colorValue= saliency;

		//int LUTBinId= fLUT->FindBin(Mean_norm);	
		//if(fUseLogNormalizedImage) LUTBinId= fLUT->FindBin(LgMean_norm);//use log mapping				
		
		cout<<"VLSlicSegmentation::GetClusterColoredImage(): INFO: Region no. "<<regionId<<" N="<<nPixelsInRegion<<" Mean="<<Mean<<" Mean_norm="<<Mean_norm<<endl;
		for(int j=0;j<nPixelsInRegion;j++){//loop on pixels inside region
			int thisPixelId= (regions[i]->fPixelCollection)[j].id;
			colored_img->SetBinContent(thisPixelId,colorValue);
		}//end loop pixels in region	
	}//end loop regions
	
	//colored_img->SetMinimum(NormMin);
	//colored_img->SetMaximum(NormMax);


	return colored_img;

}//close GetClusterColoredImage()


TGraph* VLSlicSegmentation::ComputeClusterContours(Img* image,std::vector< std::vector<long int> > pixelLabels, std::vector<Region*> regions) {
  	
	if(!image) return 0;

	int Nx= image->GetNbinsX();
	int Ny= image->GetNbinsY();
	int nRegions= (int)regions.size();
	
	cout<<"VLSlicSegmentation::GetClusterContours(): INFO: Computing contours from NR="<<regions.size()<<" regions..."<<endl;
	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};
	
	// Initialize the contour vector and the matrix detailing whether a pixel is already taken to be a contour
	//std::vector<TVector3> contours;
	
	TGraph* contourGraph= new TGraph;
	contourGraph->SetNameTitle("RegionContourGraph","RegionContourGraph");
	contourGraph->SetMarkerSize(1);
	contourGraph->SetMarkerStyle(1);
	contourGraph->SetLineColor(kWhite);
	contourGraph->SetMarkerColor(kWhite);	
	contourGraph->Set(0);
	int nContourPts= 0;

	vec2db istaken;
	
	for (int i=0; i<Nx; i++) { 
  	std::vector<bool> nb;
    for (int j = 0; j < Ny; j++) {
   		nb.push_back(false);
    }
    istaken.push_back(nb);
  }//end loop Nx
    
  // Go through all the pixels
  for (int i=0; i<Nx; i++) {
  	for (int j=0; j<Ny; j++) {
			long int pixelLabel= pixelLabels[i][j];
			if(pixelLabel<0) continue;
   		int nr_p = 0;
            
      // Compare the pixel to its 8 neighbours
      for (int k=0; k<8; k++) {
				int x = i + dx8[k];
				int y = j + dy8[k];
                
        if (x >= 0 && x < Nx && y >= 0 && y < Ny) {
					long int pixelLabel_neighbour= pixelLabels[x][y];
      	
        	if (!istaken[x][y] && pixelLabel!=pixelLabel_neighbour && pixelLabel_neighbour>=0) {	
          	nr_p++;
          }
        }
      }//end loop neighbours
            
      // Add the pixel to the contour list if desired
     	if (nr_p >= 2) {
				double col= image->GetBinContent(i+1,j+1);
				double cx= image->GetXaxis()->GetBinCenter(i+1);
				double cy= image->GetYaxis()->GetBinCenter(j+1);
      	contourGraph->SetPoint(nContourPts,cx,cy);
				//cout<<"VLSlicSegmentation::GetClusterContours(): INFO: Npt="<<nContourPts<<": Pix("<<i<<","<<j<<") ("<<cx<<","<<cy<<") label="<<pixelLabel<<" n="<<nr_p<<endl;
				nContourPts++;
        istaken[i][j] = true;
      }
    }//end lop image Ny
  }//end loop image Nx
    
	return contourGraph;

}//close VLSlicSegmentation::GetClusterContours()




