/**
* @file ChanVeseSegmentation.cc
* @class ChanVeseSegmentation
* @brief ChanVeseSegmentation
*
* @author S. Riggi
* @date 15/06/2015
*/

#include <ChanVeseSegmentation.h>

#include <TColor.h>
#include <TMath.h>
#include <TStyle.h>
#include <TColor.h>
#include <TPad.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TText.h>
#include <TPolyLine.h>
#include <TMatrixD.h>
#include <TRandom.h>

//#include <mathop.h>
#include <iostream>

using namespace std;
using std::cout;
using std::endl;
using std::atan;
using std::sqrt;
using std::abs;
using std::pow;
using std::min;
using std::max;


ClassImp(ChanVeseSegmentation)

ChanVeseSegmentation::ChanVeseSegmentation() {

	fSegmentedImg= 0;	
	fContourImg= 0;
	fContour= 0;
	fPlot= 0;

}//close constructor


ChanVeseSegmentation::~ChanVeseSegmentation() {

	/*
	if(fSegmentedImg) {
		delete fSegmentedImg;
		fSegmentedImg= 0;
	}
	*/
	if(fContour) {
		delete fContour;
		fContour= 0;
	}
	/*
	if(fContourImg) {
		delete fContourImg;
		fContourImg= 0;
	}
	*/
	if(fPlot) {	
		fPlot->Delete();
	}

}//close destructor


int ChanVeseSegmentation::RunSegmentation(Img* img,double dt,double h,double lambda1,double lambda2,double mu,double nu,double p,double initContourRadius){

	//Check input image
	if(!img){
		return -1;
	}

	//## Set up parameters
  struct CVsetup* pCVinputs = new struct CVsetup;
  pCVinputs->dt = dt;
  pCVinputs->h = h;
  pCVinputs->lambda1 = lambda1;
  pCVinputs->lambda2 = lambda2;
  pCVinputs->mu = mu;
  pCVinputs->nu = nu;
  pCVinputs->p = p;
	cout<<"== CV PARS =="<<endl;
	cout<<"lambda1="<<lambda1<<", lambda2="<<lambda2<<" mu="<<mu<<" nu="<<nu<<" p="<<p<<endl;
	cout<<"============="<<endl;

	//## Normalize image
	double norm_min= 0;
	double norm_max= 255;
	Img* img_norm= img->GetNormalizedImage(norm_min,norm_max);

	//## Convert image to matrix
	TMatrixD* imgMatrix= img_norm->GetMatrix();
	if(!imgMatrix) return -1;
	int nRows= imgMatrix->GetNrows();
	int nCols= imgMatrix->GetNcols();

	//## Init images
	TMatrixD* phi0 = new TMatrixD(nRows,nCols);
	phi0->Zero();
  TMatrixD* phi = new TMatrixD(nRows,nCols);
	phi->Zero();
  TMatrixD* edges = new TMatrixD(nRows,nCols);
	edges->Zero();
	TMatrixD* imgMatrixOut = new TMatrixD(nRows,nCols);
	imgMatrixOut->Zero();
  
	
	//## Set up initial circular contour for a 256x256 image
  double x;
  double y;	
	//double rowCenter= gRandom->Uniform(2*initContourRadius,nRows-2*initContourRadius);
	//double colCenter= gRandom->Uniform(2*initContourRadius,nCols-2*initContourRadius);
	double rowCenter= nRows/2.0;
	double colCenter= nCols/2.0;
	double mean= imgMatrix->Sum()/((double)(nRows)*(double)(nCols));
	cout<<"ChanVeseSegmentation::RunSegmentation(): INFO: Init contour center("<<rowCenter<<","<<colCenter<<")"<<endl;
	
	for (int i=0; i<nRows; ++i) {
    for (int j=0; j<nCols; ++j) {
			double x= double(i) - rowCenter;
			double y= double(j) - colCenter;

			double I= (*imgMatrix)(i,j);
			//double w = sqrt( pow(i-rowCenter,2) + pow(j-colCenter,2) ) - initContourRadius;
			double w= 900.0/(900.0 + x*x + y*y ) - initContourRadius;
			//double w = sqrt( x*x + y*y ) - initContourRadius;
			(*phi0)(i,j)= w;

			//(*phi0)(i,j)= w;
			//if(I>mean) (*phi0)(i,j)= I-mean;
			//else (*phi0)(i,j)= 0;
      //double rowId = double(i) - nRows/2.0;
      //double colId = double(j) - nCols/2.0;
			//double w= 900.0/(900.0 + rowId*rowId + colId*colId) - 0.5;
			//(*phi0)(i,j)= w;
      //double ix= colId;
			//double iy= rowId;
			//x= img->GetXaxis()->GetBinCenter(ix+1);
			//y= img->GetYaxis()->GetBinCenter(iy+1);
			//cout<<"(rowId,colId)=("<<rowId<<","<<colId<<") w="<<w<<" (x,y)=("<<x<<","<<y<<")"<<endl;
    }
  }

	fInitPlot= new TCanvas("InitContour","Init Contour");
	fInitPlot->cd();
	phi0->Draw("COLZ");
	
	cout << "ChanVeseSegmentation::RunSegmentation(): INFO: Performing segmentation...\n";
  CVSegmentation(imgMatrix,phi0, &phi, pCVinputs);

  cout << "ChanVeseSegmentation::RunSegmentation(): INFO: Getting zero crossings...\n";
  ZeroCrossings(phi,&edges,norm_max,norm_min);
  
  (*imgMatrixOut)= (*imgMatrix);
  
  
	//## Fill contour
	cout<<"ChanVeseSegmentation::RunSegmentation(): INFO: Getting contour..."<<endl;
	if(!fContour) fContour= new Contour;
	fContour->Reset();
	
	if(!fContourImg) fContourImg= (Img*)img->Clone("CVContourImg");
	fContourImg->Reset();

  for (int i = 0; i < imgMatrix->GetNrows(); ++i) {
		int iy= i;
		double y= img->GetYaxis()->GetBinCenter(iy+1);
    for (int j = 0; j < imgMatrix->GetNcols(); ++j) {
			int ix= j;
			double x= img->GetXaxis()->GetBinCenter(ix+1);
			
      if ( (*edges)(i,j) == norm_max) {
				fContour->AddPoint( cv::Point2f(x,y) );
				fContourImg->SetBinContent(ix+1,iy+1,1);
			}
    }//end loop 
  }//end loop
  
	//## Fill segmented image
	cout<<"ChanVeseSegmentation::RunSegmentation(): INFO: Filling segmented image..."<<endl;
	
	if(!fSegmentedImg) fSegmentedImg= (Img*)img->Clone("CVSegmImg");
	fSegmentedImg->Reset();
	
	for (int i = 0; i < nRows; ++i) {
		int iy= i;
		double y= img->GetYaxis()->GetBinCenter(iy+1);
    for (int j = 0; j < nCols; ++j) {
			int ix= j;
			double x= img->GetXaxis()->GetBinCenter(ix+1);
			
			if ( (*phi)(i,j) >= 0) fSegmentedImg->SetBinContent(ix+1,iy+1,norm_max);
      
      //if ( (*phi)(i,j) >= 0) fSegmentedImg->SetBinContent(ix+1,iy+1,0);
      //else fSegmentedImg->SetBinContent(ix+1,iy+1,255);
    }//end loop
  }//end loop

	fPlot= new TCanvas("SegmentedImage","Segmented Image");
	fPlot->cd();
	fSegmentedImg->Draw("COLZ");

	return 0;

}//close RunSegmentation()



//---------------------------------------------------------------------------//
// function GetRegionAverages
// Compute c1 and c2 as used in the Chan-Vese segmentation algorithm.
// c1 and c2 are given by 
//         c1 = integral(u0*H(phi))dxdy/integral(H(phi))dxdy
//         c2 = integral(u0*(1-H(phi))dxdy/integral(1-H(phi))dxdy
//
// If epsilon == 0, we define H as the usual Heaviside function. Then  c1 is 
// the average of the image pixels over the set where phi is >= 0, and c2 is
// the average over {phi < 0}.  
// If epsilon > 0, we use a smoothed version of the Heaviside function with
// parameter epsilon.
//void ChanVeseSegmentation::GetRegionAverages(Image<unsigned char>* img,Image<double>* phi,double epsilon, double &c1, double &c2) {
void ChanVeseSegmentation::GetRegionAverages(TMatrixD* img,TMatrixD* phi,double epsilon, double &c1, double &c2) {
  
	// Non-smoothed calculation
  if (0 == epsilon) {
    int n1 = 0;
    int n2 = 0;
    double Iplus = 0;
    double Iminus = 0;

    //for (unsigned int i = 0; i < img->nRows(); ++i) {
		//	for (unsigned int j = 0; j < img->nCols(); ++j) {
		for (int i=0; i<img->GetNrows(); ++i) {
			for (int j=0; j<img->GetNcols();++j) {
				double w_phi= (*phi)(i,j);
				double w= (*img)(i,j);
        if (w_phi>= 0){
          ++n1;
          Iplus += w;
        }
        else {
          ++n2;
          Iminus += w;
        }
      }//end loop bins X
    }//end loop bins Y
    c1 = Iplus/double(n1);
    c2 = Iminus/double(n2);
  }//close if non smoothed
  // Smoothed calculation
  else {
    double num1 = 0;
    double den1 = 0;
    double num2 = 0;
    double den2 = 0;
    double H_phi;
    for (int i = 0; i < phi->GetNrows(); ++i){
      for (int j = 0; j < phi->GetNcols(); ++j) {
        // Compute value of H_eps(phi) where H_eps is a mollified Heavyside function
				double w_phi= (*phi)(i,j);
				double w= (*img)(i,j);
        H_phi = .5*(1+(2/TMath::Pi())*atan(w_phi/epsilon));
        num1 += w*H_phi;
        den1 += H_phi;
        num2 += w*(1-H_phi);
        den2 += 1-H_phi;
      }//end loop bins Y
    }//end loop bins X
    
    c1 = num1/den1;
    c2 = num2/den2;
  }//close else

}//close ChanVeseSegmentation::GetRegionAverages()

//---------------------------------------------------------------------------//
// function ReinitPhi
// Reinitialize phi to the signed distance function to its zero contour
//void ChanVeseSegmentation::ReinitPhi(Image<double>* phiIn,Image<double>** psiOut,double dt,double h,unsigned int numIts) {
void ChanVeseSegmentation::ReinitPhi(TMatrixD* phiIn,TMatrixD** psiOut,double dt,double h,unsigned int numIts) {

  if (*psiOut == 0){
    //(*psiOut) = new Image<double>(phiIn->nRows(),phiIn->nCols());
		(*psiOut) = new TMatrixD(phiIn->GetNrows(),phiIn->GetNcols());
	}
  //else if ((*psiOut)->nRows() != phiIn->nRows() || (*psiOut)->nCols() != phiIn->nCols() ){
	else if ((*psiOut)->GetNrows() != phiIn->GetNrows() || (*psiOut)->GetNcols() != phiIn->GetNcols() ){
		(*psiOut)->Delete();
		(*psiOut) = new TMatrixD(phiIn->GetNrows(),phiIn->GetNcols());
    //(*psiOut)->Allocate(phiIn->nRows(),phiIn->nCols());
	}
  
  //(*psiOut)->CopyFrom(phiIn);
	(**psiOut) = *phiIn;


  double a;
  double b;
  double c;
  double d;
  double x;
  double G;
  
  bool fStop = false;
  double Q;
  unsigned int M;
  //Image<double>* psiOld = new Image<double>(phiIn->nRows(),phiIn->nCols());
  TMatrixD* psiOld = new TMatrixD(phiIn->GetNrows(),phiIn->GetNcols());
  
  for (unsigned int k = 0; k < numIts && fStop == false; ++k) {
    //psiOld->CopyFrom(*psiOut);
		*psiOld= (**psiOut);
    for (int i = 1; i < phiIn->GetNrows()-1; ++i) {
      for (int j = 1; j < phiIn->GetNcols()-1; ++j){
				double phiIn_ij= (*phiIn)(i,j);
				double phiIn_iminus1_j= (*phiIn)(i-1,j);
				double phiIn_iplus1_j= (*phiIn)(i+1,j);
				double phiIn_i_jminus1= (*phiIn)(i,j-1);
				double phiIn_i_jplus1= (*phiIn)(i,j+1);

        a = (phiIn_ij - phiIn_iminus1_j)/h;
        b = (phiIn_iplus1_j - phiIn_ij)/h;
        c = (phiIn_ij - phiIn_i_jminus1)/h;
        d = (phiIn_i_jplus1 - phiIn_ij)/h;
        
        if (phiIn_ij > 0)
          G = sqrt(max(max(a,0.0)*max(a,0.0),min(b,0.0)*min(b,0.0))
                 + max(max(c,0.0)*max(c,0.0),min(d,0.0)*min(d,0.0))) - 1.0;
        else if (phiIn_ij < 0)
          G = sqrt(max(min(a,0.0)*min(a,0.0),max(b,0.0)*max(b,0.0))
                 + max(min(c,0.0)*min(c,0.0),max(d,0.0)*max(d,0.0))) - 1.0;
        else
          G = 0;
        
        x = (phiIn_ij >= 0)?(1.0):(-1.0);
				(**psiOut)(i,j)= (**psiOut)(i,j) - dt*x*G;
        //(*psiOut)->data()[i][j] = (*psiOut)->data()[i][j] - dt*x*G;
      }
    }
    
    // Check stopping condition
    Q = 0.0;
    M = 0.0;
    for (int i = 0; i < phiIn->GetNrows(); ++i) {
      for (int j = 0; j < phiIn->GetNcols(); ++j) {
				double w_old= (*psiOld)(i,j);	
				double w_phi= (**psiOut)(i,j);
        if (abs(w_old) <= h) {
        	++M;
          Q += abs(w_old - w_phi);
        }
      }
    }
    if (M != 0)
      Q = Q/((double)M);
    else
      Q = 0.0;
    
    if (Q < dt*h*h)
    {
      fStop = true;
      cout << "ChanVeseSegmentation::ReinitPhi(): INFO: Stopping condition reached at " << k+1 << " iterations; Q = " << Q << endl;
    }
    else
    {
      //cout << "Iteration " << k << ", Q = " << Q << " > " << dt*h*h << endl;
    }
  }

}//close ChanVeseSegmentation::ReinitPhi()

//---------------------------------------------------------------------------//
// function GetChanVeseCoefficients
// Compute coefficients needed in Chan-Vese segmentation algorithm given current 
// level set function
void ChanVeseSegmentation::GetChanVeseCoefficients(TMatrixD* phi,
                             struct CVsetup* pCVinputs,
                             unsigned int i,
                             unsigned int j,
                             double L,
                             double& F1,
                             double& F2,
                             double& F3,
                             double& F4,
                             double& F,
                             double& deltaPhi)
{
  // factor to avoid division by zero
  double eps = 0.000001;
  double h = pCVinputs->h;
  double dt = pCVinputs->dt;
  double mu = pCVinputs->mu;
  unsigned int p = pCVinputs->p;
  
	double phi_ij= (*phi)(i,j);	
	double phi_iminus1_j= (*phi)(i-1,j);
	double phi_iplus1_j= (*phi)(i+1,j);
	double phi_i_jminus1= (*phi)(i,j-1);
	double phi_i_jplus1= (*phi)(i,j+1);
	double phi_iminus1_jminus1= (*phi)(i-1,j-1);
	double phi_iplus1_jminus1= (*phi)(i+1,j-1);
	double phi_iminus1_jplus1= (*phi)(i-1,j+1);
	

  double C1 = 1/sqrt(eps + pow((phi_iplus1_j - phi_ij),2)
                         + pow((phi_i_jplus1 - phi_i_jminus1),2)/4.0);
  double C2 = 1/sqrt(eps + pow((phi_ij - phi_iminus1_j),2)
                         + pow((phi_iminus1_jplus1 - phi_iminus1_jminus1),2)/4.0);
  double C3 = 1/sqrt(eps + pow((phi_iplus1_j - phi_iminus1_j),2)/4.0
                         + pow((phi_i_jplus1- phi_ij),2));
  double C4 = 1/sqrt(eps + pow((phi_iplus1_jminus1 - phi_iminus1_jminus1),2)/4.0
                         + pow((phi_ij - phi_i_jminus1),2));

  deltaPhi = h/(TMath::Pi()*(h*h + (phi_ij)*(phi_ij)));
  
  double Factor = dt*deltaPhi*mu*(double(p)*pow(L,p-1));
  F = h/(h+Factor*(C1+C2+C3+C4));
  Factor = Factor/(h+Factor*(C1+C2+C3+C4));
  
  F1 = Factor*C1;
  F2 = Factor*C2;
  F3 = Factor*C3;
  F4 = Factor*C4;
}

//---------------------------------------------------------------------------//
// Main segmentation algorithm.  Segment a grayscale image into foreground and
// background regions, given an initial contour defined by the level set function
// phi.  Based on the algorithm described in the paper
// "Active Contours Without Edges" by Chan & Vese.
//void ChanVeseSegmentation::ChanVeseSegmentation(Image<unsigned char>* img,Image<double>* phi0,Image<double>** phi,struct CVsetup* pCVinputs) {
void ChanVeseSegmentation::CVSegmentation(TMatrixD* img,TMatrixD* phi0,TMatrixD** phi,struct CVsetup* pCVinputs) {
  
	double P_ij;
  double deltaPhi;
  double F1;
  double F2; 
  double F3;
  double F4;
  double F;
  double L;
  double c1;
  double c2;
  
  // Segmentation parameters
  double h = pCVinputs->h;
  double dt = pCVinputs->dt;
  double nu = pCVinputs->nu;
  double lambda1 = pCVinputs->lambda1;
  double lambda2 = pCVinputs->lambda2;
  unsigned int p = pCVinputs->p;
  
  // Variables to evaluate stopping condition
  bool fStop = false;
  double Q;
  unsigned int M;
  TMatrixD* phiOld = new TMatrixD(img->GetNrows(),img->GetNcols());
  
  // Initialize phi
	if (*phi == 0){
    (*phi) = new TMatrixD(img->GetNrows(),img->GetNcols());
	}
 	else if ((*phi)->GetNrows() != phi0->GetNrows() || (*phi)->GetNcols() != phi0->GetNcols() ){
		(*phi)->Delete();
		(*phi) = new TMatrixD(phi0->GetNrows(),phi0->GetNcols());
	}
  
  (**phi) = *phi0;


	/*
  if (0 == phi){
    *phi = new Image<double>(img->nRows(),img->nCols());
		
	}
  else if ((*phi)->nRows() != phi0->nRows() || (*phi)->nCols() != phi0->nCols()){
    (*phi)->Allocate(phi0->nRows(),phi0->nCols());
	}
  
  (*phi)->CopyFrom(phi0);
  */
  
  // Main loop
  for (unsigned int k = 0; k < 5 && fStop == false; ++k)
  {
    //phiOld->CopyFrom(*phi);
    *phiOld= (**phi);
    

    // Compute region averages for current level set function
    // Main segmentation algorithm
    GetRegionAverages(img, *phi, h, c1, c2);

    // Inner loop...
    for (unsigned int l = 0; l < 5; ++l) {

      // Compute length of contour if p > 1
      if (1 == p) L = 1.0;
      else L = 1.0; // fix this!!
     
      
      // Loop through all interior image pixels
      for (int i = 1; i < img->GetNrows()-1; ++i) {
        for (int j = 1; j < img->GetNcols()-1; ++j) {
          // Compute coefficients needed in update
          GetChanVeseCoefficients(*phi,pCVinputs,i, j,L,F1,F2,F3,F4,F,deltaPhi);

					double phi_ij= (**phi)(i,j);
					double w= (*img)(i,j);
					P_ij = phi_ij - dt*deltaPhi*(nu + lambda1*pow(w-c1,2) - lambda2*pow(w-c2,2));

          //P_ij = (*phi)->data()[i][j] - dt*deltaPhi*(nu + lambda1*pow(w-c1,2) - lambda2*pow(w-c2,2));

          // Update level set function
					double phi_iplus1_j= (**phi)(i+1,j);
					double phi_iminus1_j= (**phi)(i-1,j);
					double phi_i_jplus1= (**phi)(i,j+1);
					double phi_i_jminus1= (**phi)(i,j-1);
					(**phi)(i,j) = F1*phi_iplus1_j + F2*phi_iminus1_j + F3*phi_i_jplus1 + F4*phi_i_jminus1 + F*P_ij;

          //(*phi)->data()[i][j] = F1*(*phi)->data()[i+1][j] + F2*(*phi)->data()[i-1][j] + F3*(*phi)->data()[i][j+1] + F4*(*phi)->data()[i][j-1] + F*P_ij;
        }
      }
      
      // Update border values of phi by reflection
      for (int i = 0; i < img->GetNrows(); ++i) {
				(**phi)(i,0)= (**phi)(i,1);
				(**phi)(i,img->GetNcols()-1)= (**phi)(i,img->GetNcols()-2);
				 
        //(*phi)->data()[i][0] = (*phi)->data()[i][1];
        //(*phi)->data()[i][img->nCols()-1] = (*phi)->data()[i][img->nCols()-2];
      }
      for (int j = 0; j < img->GetNcols(); ++j) {
				(**phi)(0,j)= (**phi)(1,j);
				(**phi)(img->GetNrows()-1,j)= (**phi)(img->GetNrows()-2,j);

        //(*phi)->data()[0][j] = (*phi)->data()[1][j];
        //(*phi)->data()[img->nRows()-1][j] = (*phi)->data()[img->nRows()-2][j];
      }
 
      // Reinitialize phi to the signed distance function to its zero contour
			int nIterations= 200;//100
      ReinitPhi(*phi, phi, 0.1, h, nIterations);
    }
    
    // Check stopping condition
    //...
  }
}//close CVSegmentation()



void ChanVeseSegmentation::ZeroCrossings(TMatrixD* imageIn,TMatrixD** edges,double fg,double bg) {

  // Allocate output image if necessary
  if (0 == (*edges)){
		(*edges) = new TMatrixD(imageIn->GetNrows(),imageIn->GetNcols());
  }
  else {
    if ((*edges)->GetNrows() != imageIn->GetNrows() || (*edges)->GetNcols() != imageIn->GetNcols()) {
			(*edges)->Delete();
			(*edges) = new TMatrixD(imageIn->GetNrows(),imageIn->GetNcols());
    } 
  }
  (**edges)= bg;

  for (int i = 0; i < imageIn->GetNrows(); ++i) {
    for (int j = 0; j < imageIn->GetNcols(); ++j) {
      // Currently only checking interior pixels to avoid bounds checking
      if (i > 0 && i < (imageIn->GetNrows()-1) && j > 0 && j < (imageIn->GetNcols()-1)) {

				double w_ij= (*imageIn)(i,j);
				double w_iminus1_jminus1= (*imageIn)(i-1,j-1);
				double w_iplus1_jplus1= (*imageIn)(i+1,j+1);
				double w_iminus1_j= (*imageIn)(i-1,j);
				double w_iplus1_j= (*imageIn)(i+1,j);
				double w_iplus1_jminus1= (*imageIn)(i+1,j-1);
				double w_iminus1_jplus1= (*imageIn)(i-1,j+1);
				double w_i_jminus1= (*imageIn)(i,j-1);
				double w_i_jplus1= (*imageIn)(i,j+1);

        if (0 == w_ij) {
          if (0 != w_iminus1_jminus1
           || 0 != w_iminus1_j
           || 0 != w_iminus1_jplus1
           || 0 != w_i_jminus1
           || 0 != w_i_jplus1
           || 0 != w_iplus1_jminus1
           || 0 != w_iplus1_j
           || 0 != w_iplus1_jplus1)
           {
             (**edges)(i,j) = fg;
           }
        }//close if
        else {
          if (abs(w_ij) < abs(w_iminus1_jminus1) && (w_ij>0) != (w_iminus1_jminus1>0))
             (**edges)(i,j) = fg;
     			else if (abs(w_ij) < abs(w_iminus1_j) && (w_ij>0) != (w_iminus1_j>0))
             (**edges)(i,j) = fg;
     			else if (abs(w_ij) < abs(w_iminus1_jplus1) && (w_ij>0) != (w_iminus1_jplus1>0))
             (**edges)(i,j) = fg;
     			else if (abs(w_ij) < abs(w_i_jminus1) && (w_ij>0) != (w_i_jminus1>0))
             (**edges)(i,j) = fg;
     			else if (abs(w_ij) < abs(w_i_jplus1) && (w_ij>0) != (w_i_jplus1>0))
             (**edges)(i,j) = fg;
     			else if (abs(w_ij) < abs(w_iplus1_jminus1) && (w_ij>0) != (w_iplus1_jminus1>0))
             (**edges)(i,j) = fg;
     			else if (abs(w_ij) < abs(w_iplus1_j) && (w_ij>0) != (w_iplus1_j>0))
             (**edges)(i,j) = fg;
     			else if (abs(w_ij) < abs(w_iplus1_jplus1) && (w_ij>0) != (w_iplus1_jplus1>0))
             (**edges)(i,j) = fg;
        }
      }
    }
  }
}//close ZeroCrossing()




