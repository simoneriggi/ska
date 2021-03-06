/**
* @file Contour.cc
* @class Contour
* @brief Contour
*
* Class representing image contour with methods for morphological parameter extraction 
* @author S. Riggi
* @date 11/07/2015
*/

#include <Contour.h>
#include <Utils.h>

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
#include <TEllipse.h>
#include <TPolyLine.h>
#include <TPaveText.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMatrixDEigen.h>

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

using namespace std;
//using namespace cv;


ClassImp(Contour)

Contour::Contour(){

	fPoints.clear();	
	HasParameters= false;	
	HasEllipseFit= false;
	EllipseMajAxis= 0;
	EllipseMinAxis= 0;
	EllipseRotAngle= 0;
	EllipseFitRedChi2= 1.e+99;
	EllipseAreaRatio= 0;
	Area= 0.;
	Perymeter= 0.;
	IsConvexContour= false;
	CircularityRatio= -999;
	BoundingBoxMaj= -999;
	BoundingBoxMin= -999;
	BoundingBoxAngle= -999;
	Elongation= -999;
	Rectangularity= -999;
	Roundness= -999;
	Eccentricity= -999;
	TiltAngle= -999;
	//Moments= ...
	for(int k=0;k<7;k++) HuMoments[k]= 0;
	for(int k=0;k<4;k++) BoundingBoxVertex[k]= cv::Point2f(0,0);	
	BoundingBoxCenter= cv::Point2f(0,0);
	Centroid= cv::Point2f(0,0);
	FDs.clear();
	ModFDs.clear();


}//close costructor


Contour::~Contour(){

}//close destructor


TGraph* Contour::GetGraph(){

	//Check number of contour pts
	int nContourPts= fPoints.size();
	if(nContourPts<=0) {
		cerr<<"Contour::GetGraph(): WARN: No contour points available (did you fill the contour?), returning null ptr graph!"<<endl;
		return 0;
	}
	//Fill contour graph
	TGraph* ContourGraph= new TGraph(nContourPts+1);
	for(int i=0;i<nContourPts;i++){
		//cout<<"Contour::GetGraph(): Pnt no. "<<i<<" (x,y)= ("<<fPoints[i].x<<","<<fPoints[i].y<<")"<<endl;
		ContourGraph->SetPoint(i,fPoints[i].x,fPoints[i].y);
	}
	ContourGraph->SetPoint(nContourPts,fPoints[0].x,fPoints[0].y);//Add another point as OpenCV does not close the contour

	return ContourGraph;

}//close Contour::GetGraph()


TPolyLine* Contour::GetBoundingBoxLine(){

	if(!HasPoints() || !HasParameters){
		cerr<<"Contour::GetBoundingBoxLine(): ERROR: No contour points and/or parameters available (did you fill the contour and compute its parameters?)!"<<endl;
		return 0;
	}

	TPolyLine* boundingRectangle= new TPolyLine(5);	
	for(int k=0;k<4;k++) boundingRectangle->SetPoint(k,BoundingBoxVertex[k].x,BoundingBoxVertex[k].y);
	boundingRectangle->SetPoint(4,BoundingBoxVertex[0].x,BoundingBoxVertex[0].y);	
	return boundingRectangle;

}//close Contour::GetBoundingBoxLine()

TPaveText* Contour::GetParamInfoBox(){

	TPaveText* infoBox = new TPaveText(0.15,0.7,0.75,0.85,"NDC");
	infoBox->AddText(Form("Eccentricity: %1.2f",Eccentricity));
	infoBox->AddText(Form("Elongation: %1.2f",Elongation));
	infoBox->AddText(Form("TiltAngle(deg): %1.2f",TiltAngle));
	infoBox->AddText(Form("Rectangularity: %1.2f",Rectangularity));
	infoBox->AddText(Form("Roundness: %1.2f, CircRatio: %1.2f",Roundness,CircularityRatio));
	infoBox->AddText(Form("BoundingBox: (%1.2f,%1.2f,%1.2f)",BoundingBoxMaj,BoundingBoxMin,BoundingBoxAngle));
	infoBox->AddText(Form("HuMoments: (%1.2g,%1.2g,%1.2g,%1.2g,%1.2g,%1.2g,%1.2g)",HuMoments[0],HuMoments[1],HuMoments[2],HuMoments[3],HuMoments[4],HuMoments[5],HuMoments[6]));
	if(HasEllipseFit) infoBox->AddText(Form("(x0,y0,a,b,Theta): (%1.2f,%1.2f,%1.2f,%1.2f,%1.2f), EllipseAreaRatio: %1.2f",EllipseCenter.x,EllipseCenter.y,EllipseMajAxis/2,EllipseMinAxis/2,EllipseRotAngle,EllipseAreaRatio));
	infoBox->SetTextAlign(12);
	infoBox->SetTextSize(0.02);
	infoBox->SetTextFont(52);
	infoBox->SetFillColor(0);
	infoBox->SetBorderSize(1);

	return infoBox;

}//close GetParamInfoBox()

int Contour::ComputeParameters(){

	int nContourPts= fPoints.size();
	if(nContourPts<=0) {
		cerr<<"Contour::ComputeParameters(): WARN: No contour points available (did you fill the contour?)!"<<endl;
		return -1;
	}
	if(nContourPts<4) {
		cerr<<"Contour::ComputeParameters(): WARN: Too few contour points available (n="<<nContourPts<<") to get any reliable parameter estimate!"<<endl;
		return -1;
	}

	try{
		//## Compute Area
		ComputeArea();
	
		//## Compute perimeter
		ComputePerymeter();

		//## Compute circularity ratio (require area & perimeter for that!)
		ComputeCircularityRatio();

		//## Compute bounding box and related params
		ComputeBoundingBox();
		ComputeElongation();//require bounding box 
		ComputeRectangularity();//require bounding box & area
		ComputeRoundness();//require bounding box & area

		//## Compute fitted ellipse
		ComputeFittedEllipse();

		//## Compute moments and HuMoments
		ComputeMoments();
		ComputeHuMoments();//require moment for HuMoments!
	
		//## Compute Eccentricity 
		ComputeEccentricity();//require moments!

		//## Compute Fourier descriptor
		//ComputeFourierDescriptors();
		//ComputeCentroidDistanceFD();

		//## Compute average bending energy
		//ComputeBendingEnergy();		

		HasParameters= true;
	}//close try block
	catch(cv::Exception ex){//something goes wrong!
  	cout<<"Contour::ComputeParameters(): ERROR: Computing contour parameters failed with status: "<<ex.msg <<endl;
		HasParameters= false;
		return -1;
  }		

	return 0;

}//close Contour::ComputeParameters()

void Contour::ComputeArea(){	
	Area= cv::contourArea(fPoints,false);
}//close ComputeArea()

void Contour::ComputePerymeter(){	
	bool isContourClosed= true;//assuming contours are always closed
	Perymeter= cv::arcLength(fPoints, isContourClosed);
}//close ComputePerymeter()

void Contour::ComputeCircularityRatio(){
	CircularityRatio= 4*TMath::Pi()*Area/pow(Perymeter,2);
}//close ComputeCircularityRatio()

void Contour::ComputeBoundingBox(){
	//cv::Rect BoundingRect= cv::boundingRect(fPoints);//bounding box
	cv::RotatedRect MinBoundingRect= cv::minAreaRect(fPoints);//rotated bounding box
	MinBoundingRect.points(BoundingBoxVertex);//bounding box vertexes
	double MinBoundingRect_height= MinBoundingRect.size.height;
	double MinBoundingRect_width= MinBoundingRect.size.width;
	BoundingBoxCenter= MinBoundingRect.center;
	BoundingBoxMaj= std::max(MinBoundingRect_height,MinBoundingRect_width);
	BoundingBoxMin= std::min(MinBoundingRect_height,MinBoundingRect_width);
	BoundingBoxAngle= MinBoundingRect.angle;//counterclockwise in degree 
	if(MinBoundingRect.size.width < MinBoundingRect.size.height){
    BoundingBoxAngle+= 180;
  }
	else{
		BoundingBoxAngle+= 90;  
	}

}//close ComputeBoundingBox()

void Contour::ComputeElongation(){
	Elongation= 1.-BoundingBoxMin/BoundingBoxMaj;
}
	
void Contour::ComputeRectangularity(){	
	Rectangularity= Area/(BoundingBoxMaj*BoundingBoxMin);
}
		
void Contour::ComputeRoundness(){
	Roundness= 4.*Area/(TMath::Pi()*pow(BoundingBoxMaj,2));
}

TEllipse* Contour::GetFittedEllipse(){

	if(!HasEllipseFit) return 0;
	TEllipse* ellipse= new TEllipse;
	ellipse->SetX1(EllipseCenter.x);
	ellipse->SetY1(EllipseCenter.y);
	//ellipse->SetTheta(-EllipseRotAngle);
	ellipse->SetTheta(EllipseRotAngle);
	ellipse->SetR1(EllipseMajAxis/2.);
	ellipse->SetR2(EllipseMinAxis/2.);	
	return ellipse;

}//close GetFittedEllipse()

void Contour::ComputeFittedEllipse(){
	
	HasEllipseFit= true;//checked in internal routines and swithed to false if the fit fails
	//Fit an ellipse to contour points
	
	
	// Fit the ellipse
	TGraph* contourGraph= GetGraph();
  TVectorD conic = EllipseFitter(contourGraph);
  TVectorD ellipseParams = ConicToParametric(conic);

	if(!HasEllipseFit){
		cerr<<"Contour::ComputeFittedEllipse(): WARN: Ellipse fit failed!"<<endl;
		if(contourGraph) contourGraph->Delete();	
		return;
	}

	//ellipse[0] = x0; // ellipse's "x" center
  //ellipse[1] = y0; // ellipse's "y" center
  //ellipse[2] = a; // ellipse's "semimajor" axis along "x"
  //ellipse[3] = b; // ellipse's "semiminor" axis along "y"
  //ellipse[4] = theta; // ellipse's axes rotation angle (in degrees)
	EllipseCenter= cv::Point2f(ellipseParams[0],ellipseParams[1]);
	//EllipseMajAxis= 2*std::max(ellipseParams[2],ellipseParams[3]);
	//EllipseMinAxis= 2*std::min(ellipseParams[2],ellipseParams[3]);
	EllipseMajAxis= 2*ellipseParams[2];
	EllipseMinAxis= 2*ellipseParams[3];
	EllipseRotAngle= ellipseParams[4];

	//Compute the area ratio
	double EllipseArea= TMath::Pi()*ellipseParams[2]*ellipseParams[3];
	EllipseAreaRatio= Area/EllipseArea;
	
	//Compute the fit residual
	double chi2= EllipseFitChi2(contourGraph,ellipseParams);
	double ndf= contourGraph->GetN()-ellipseParams.GetNoElements();
	EllipseFitRedChi2= chi2/ndf;

	
	cout<<"*** ELLIPSE FIT ***"<<endl;
	cout<<"(x0,y0)=("<<ellipseParams[0]<<","<<ellipseParams[1]<<")"<<endl;
	cout<<"a="<<ellipseParams[2]<<" b="<<ellipseParams[3]<<endl;
	cout<<"theta(deg)="<<ellipseParams[4]<<endl;
	cout<<"chi2="<<chi2<<" redchi2="<<EllipseFitRedChi2<<endl;
	cout<<"*******************"<<endl;
	
	/*
	//OpenCV method
	std::vector<cv::Point2f> hull;
	cv::convexHull(fPoints, hull);
	//cv::RotatedRect fittedEllipse = cv::fitEllipse(cv::Mat(fPoints));
	cv::RotatedRect fittedEllipse = cv::fitEllipse(hull);
	EllipseCenter= fittedEllipse.center;
	EllipseMajAxis= std::max(fittedEllipse.size.width,fittedEllipse.size.height);
	EllipseMinAxis= std::min(fittedEllipse.size.width,fittedEllipse.size.height);
	EllipseRotAngle= fittedEllipse.angle;//in degrees
	EllipseArea= TMath::Pi()*fittedEllipse.size.width/2*fittedEllipse.size.height/2;
	EllipseAreaRatio= Area/EllipseArea;
	cout<<"*** OPENCV ELLIPSE FIT ***"<<endl;
	cout<<"(x0,y0)=("<<EllipseCenter.x<<","<<EllipseCenter.y<<")"<<endl;
	cout<<"a="<<EllipseMajAxis/2<<" b="<<EllipseMinAxis/2<<endl;
	cout<<"theta(deg)="<<EllipseRotAngle<<endl;
	cout<<"*******************"<<endl;	
	*/
	if(contourGraph) contourGraph->Delete();
	
}//close ComputeFittedEllipse()

 
TVectorD Contour::EllipseFitter(TGraph* contourGraph){

  TVectorD ellipse;
  if (!contourGraph) {
		HasEllipseFit= false;
		return ellipse; // just a precaution
	}

	//Check number of points
	int N = contourGraph->GetN();
  if(N<6) {
		HasEllipseFit= false;
		return ellipse;
	}
  
	int i= 0;
  double tmp= 0;
  double xmin, xmax, ymin, ymax, X0, Y0;
  contourGraph->ComputeRange(xmin, ymin, xmax, ymax);

	#if 1 // 0 or 1 
	  X0 = (xmax + xmin) / 2.0;
  	Y0 = (ymax + ymin) / 2.0;
	#else // 0 or 1 
  	X0 = Y0 = 0.0;
	#endif // 0 or 1
  
  TMatrixD D1(N, 3); // quadratic part of the design matrix
  TMatrixD D2(N, 3); // linear part of the design matrix
  
  for(i=0;i<N;i++) {
    double x = (contourGraph->GetX())[i] - X0;
    double y = (contourGraph->GetY())[i] - Y0;
    D1[i][0] = x*x;
    D1[i][1] = x*y;
    D1[i][2] = y*y;
    D2[i][0] = x;
    D2[i][1] = y;
    D2[i][2] = 1.0;
  }
  
  // Quadratic part of the scatter matrix
  TMatrixD S1(TMatrixD::kAtA, D1);

  // Combined part of the scatter matrix
  TMatrixD S2(D1, TMatrixD::kTransposeMult, D2);

  // Linear part of the scatter matrix
  TMatrixD S3(TMatrixD::kAtA, D2);
  S3.Invert(&tmp); 
	S3*= -1.0;
  if (tmp == 0.0) {
    cout << "Contour::EllipseFitter(): WARN: Linear part of the scatter matrix is singular!" << endl;
		HasEllipseFit= false;
    return ellipse;
  }
  // For getting a2 from a1
  TMatrixD T(S3, TMatrixD::kMultTranspose, S2);

  // Reduced scatter matrix
  TMatrixD M(S2, TMatrixD::kMult, T); M += S1;

  // Premultiply by inv(C1)
  for (i=0;i<3;i++) {
    tmp = M[0][i]/2.0;
    M[0][i] = M[2][i]/2.0;
    M[2][i] = tmp;
    M[1][i]*= -1.0;
  }

  // Solve eigensystem
  TMatrixDEigen eig(M); // note: eigenvectors are not normalized
  const TMatrixD &evec = eig.GetEigenVectors();
  // const TVectorD &eval = eig.GetEigenValuesRe();
  if ((eig.GetEigenValuesIm()).Norm2Sqr() != 0.0) {
    cout << "Contour::EllipseFitter(): WARN: Eigenvalues have nonzero imaginary parts!" << endl;
		HasEllipseFit= false;
    return ellipse;
  }

  // Evaluate a’Ca (in order to find the eigenvector for min. pos. eigenvalue)
  for (i=0; i<3; i++) {
    tmp = 4.0 * evec[0][i] * evec[2][i] - evec[1][i] * evec[1][i];
    if (tmp > 0.0) break;
  }
  if (i > 2) {
    cout << "Contour::EllipseFitter(): WARN: No min. pos. eigenvalue found!" << endl;
    // i = 2;
		HasEllipseFit= false;
    return ellipse;
  }

  // Eigenvector for min. pos. eigenvalue
  TVectorD a1(TMatrixDColumn_const(evec, i));
  tmp = a1.Norm2Sqr();
  if (tmp > 0.0) {
    a1*= 1.0/sqrt(tmp); // normalize this eigenvector
  } 
	else {
    cout << "Contour::EllipseFitter(): WARN: Eigenvector for min. pos. eigenvalue is NULL!" <<endl;
		HasEllipseFit= false;
    return ellipse;
  }
  TVectorD a2(T*a1);
  
  // Ellipse coefficients
  ellipse.ResizeTo(8);
  ellipse[0] = X0; // "X0"
  ellipse[1] = Y0; // "Y0"
  ellipse[2] = a1[0]; // "A"
  ellipse[3] = a1[1]; // "B"
  ellipse[4] = a1[2]; // "C"
  ellipse[5] = a2[0]; // "D"
  ellipse[6] = a2[1]; // "E"
  ellipse[7] = a2[2]; // "F"
  
  return ellipse;
}//close EllipseFitter()

TVectorD Contour::ConicToParametric(const TVectorD &conic) {
  
	TVectorD ellipse;  
  if (conic.GetNrows() != 8) {
    cout << "Contour::ConicToParametric(): ERROR: Improper input vector length!" << endl;
		HasEllipseFit= false;
    return ellipse;
  }
  
  double a, b, theta;
  double x0 = conic[0]; // = X0
  double y0 = conic[1]; // = Y0
  
  // http://mathworld.wolfram.com/Ellipse.html
  double A = conic[2];
  double B = conic[3] / 2.0;
  double C = conic[4];
  double D = conic[5] / 2.0;
  double F = conic[6] / 2.0;
  double G = conic[7];
  
  double J = B * B - A * C;
  double Delta = A * F * F + C * D * D + J * G - 2.0 * B * D * F;
  double I = - (A + C);
  
  // http://mathworld.wolfram.com/QuadraticCurve.html
  if (!( (Delta != 0.0) && (J < 0.0) && (I != 0.0) && (Delta / I < 0.0) )) {
    cout << "Contour::ConicToParametric(): ERROR: Ellipse (real) specific constraints not met!" << endl;
		HasEllipseFit= false;
    return ellipse;
  }
  
  x0 += (C * D - B * F) / J;
  y0 += (A * F - B * D) / J;
  
  double tmp = sqrt((A - C) * (A - C) + 4.0 * B * B);
  a = sqrt(2.0 * Delta / J / (I + tmp));
  b = sqrt(2.0 * Delta / J / (I - tmp));
  
  theta = 0.0;
  if (B != 0.0) {
    tmp = (A - C) / 2.0 / B;
    theta = -45.0 * (std::atan(tmp) / TMath::PiOver2());
    if (tmp < 0.0) { theta -= 45.0; } else { theta += 45.0; }
    if (A > C) theta += 90.0;
  } 
	else if (A > C) theta = 90.0;
  
  // try to keep "a" > "b"
  if (a < b) { tmp = a; a = b; b = tmp; theta -= 90.0; }
  // try to keep "theta" = -45 ... 135 degrees
  if (theta < -45.0) theta += 180.0;
  if (theta > 135.0) theta -= 180.0;
  
  // ellipse coefficients
  ellipse.ResizeTo(5);
  ellipse[0] = x0; // ellipse's "x" center
  ellipse[1] = y0; // ellipse's "y" center
  ellipse[2] = a; // ellipse's "semimajor" axis along "x"
  ellipse[3] = b; // ellipse's "semiminor" axis along "y"
  ellipse[4] = theta; // ellipse's axes rotation angle (in degrees)
  
  return ellipse;

}//close ConicToParametric()




void Contour::ComputeMoments(){
	// spatial moments: m00, m10, m01, m20, m11, m02, m30, m21, m12, m03
  // central moments: mu20, mu11, mu02, mu30, mu21, mu12, mu03
  // central normalized moments: nu20, nu11, nu02, nu30, nu21, nu12, nu03
	Moments = cv::moments(fPoints);
	Centroid.x = Moments.m10/Moments.m00; 
	Centroid.y = Moments.m01/Moments.m00;
}

void Contour::ComputeHuMoments(){
	cv::HuMoments(Moments, HuMoments);
}


void Contour::ComputeEccentricity(){
	//Compute covariance matrix and its eigenvectors
	double Cxx= Moments.mu20/Moments.m00;
	double Cyy= Moments.mu02/Moments.m00;
	double Cxy= Moments.mu11/Moments.m00;
	double delta= sqrt(4*Cxy*Cxy+(Cxx-Cyy)*(Cxx-Cyy));
	double lambda1= ((Cxx+Cyy) + delta)/2.; 
	double lambda2= ((Cxx+Cyy) - delta)/2.;	
	Eccentricity= sqrt(1-lambda2/lambda1);
	TiltAngle= 0.5*atan(2.*Cxy/(Cxx-Cyy))*TMath::RadToDeg();
}


void Contour::ComputeCentroidDistanceFD(){

	//Reset lists
	CentroidDistanceModFDs.clear();
	int N= (int)fPoints.size();
	int n= N;//computing all Fourier descriptors (truncate to 10 in case)

	//Compute centroid distance function r(t)= sqrt( (x-xc)^2 + (y-yc)^2 )
	std::vector<double> r;
	double Xc= Centroid.x;
	double Yc= Centroid.y;
	for(int i=0;i<N;i++){
		double X= fPoints[i].x;
		double Y= fPoints[i].y;
		double dist= sqrt( pow(X-Xc,2) + pow(Y-Yc,2) );
		r.push_back(dist);
	}//end loop contour points
	
	//Compute the Discrete Fourier Transform of r	
	std::vector< std::complex<double> > Fn;
	for(int i=0; i<n; i++) {//loop over n
		Fn.push_back( std::complex<double>(0.,0.) );
		//int s= -floor(N/2.) + i;
		int s= i;

		for(int j=0;j<N;j++) {//loop over data size
			int k= j;
			double arg= 2.*TMath::Pi()*s*k/N;
			std::complex<double> prod= std::polar(1.,-arg);
			Fn[i]+= r[j]*prod;
		}//end loop data size
		Fn[i]/= N;
	}//end loop n

	for(unsigned int i=0; i<Fn.size(); i++) {
		//Normalize FD dividing by first FD
		Fn[i]/= Fn[1];

		//Compute modulus
		double FDMod= std::abs(Fn[i]); 
		CentroidDistanceModFDs.push_back(FDMod);
	}//end loop n
		
}//close Contour::ComputeCentroidDistanceFD()



void Contour::ComputeFourierDescriptors(){

	//Reset lists
	FDs.clear();
	ModFDs.clear();
	int N= (int)fPoints.size();
	int n= N;//computing all Fourier descriptors (truncate to 10 in case)

	//Put contour point in a vector of complex numbers U
	std::vector< std::complex<double> > U= GetComplexPointRepresentation(true);//shift in centroid for convenience
		
	//Compute the Discrete Fourier Transform of contour	complex points	
	std::vector< std::complex<double> > Fn= Utils::DFTShifted(U,n);
		
	//Compute the Fourier descriptors 
	//- translational invariance ==> set F(0)=0
	//- scale invariance ==> set F(k)= F(k)/|F(1)|
	//- rotational invariance ==>?
	//- invariance against shifting of the initial point ==> set F(k)= F(k) exp(-i phi(1) k) 
	double phi_1= std::arg(Fn[1+floor(N/2.)]);//in radians
	double phi_minus1= std::arg(Fn[-1+floor(N/2.)]);//in radians
	double mod_1= std::abs(Fn[1+floor(N/2.)]);
	double scaleInvariance= mod_1;
	
	for(unsigned int k=0;k<Fn.size();k++){
		int index= -floor(N/2.) + k;
		//std::complex<double> shiftInvariance= std::polar(1.,-phi_1*index);
		//if(index==0) Fn[k]= std::complex<double>(0.,0.);
		//else Fn[k]*= shiftInvariance/scaleInvariance;

		double Fn_mod= std::abs(Fn[k]);
		double Fn_phase= std::arg(Fn[k]);
		
		std::complex<double> Fn_norm= std::complex<double>(0.,0.);
		if(index!=0) Fn_norm= std::polar(Fn_mod/mod_1,Fn_phase-0.5*(phi_minus1+phi_1)+index*0.5*(phi_minus1-phi_1) );
		 
		//if(index==0) Fn[k]= std::complex<double>(0.,0.);
		//else Fn[k]/= scaleInvariance;
		//double FDMod= std::abs(Fn[k]);
		double FDMod= std::abs(Fn_norm);
	
		FDs.push_back(Fn_norm);
		ModFDs.push_back(FDMod);
		cout<<"Contour::ComputeFourierDescriptors(): INFO: FD no. "<<k<<" scale="<<index<<" FD="<<real(Fn[k])<<" + i "<<imag(Fn[k])<<" |FD|="<<FDMod<<endl;
	}//end loop fourier coeff

}//close ComputeFourierDescriptors()


void Contour::ComputeBendingEnergy(){

	int N= (int)fPoints.size();
	
	//Put contour point in a vector of complex numbers ut
	std::vector< std::complex<double> > ut= GetComplexPointRepresentation(true);//translate points in centroid coordinate system
			
	//Compute ut'= IDST(i x 2pi x s x Us)	
	std::vector< std::complex<double> > Us= Utils::DFT(ut,N);
	std::vector< std::complex<double> > Us_firstDeriv; 
	for(int i=0;i<N;i++){
		int s= -floor(N/2.) + i;
		double arg= 2*TMath::Pi()*s;
		std::complex<double> z(0,arg); 
		std::complex<double> thisUs_firstDeriv= z*Us[i];	
		Us_firstDeriv.push_back(thisUs_firstDeriv);
	}//end loop points

	std::vector< std::complex<double> > ut_firstDeriv= Utils::IDFT(Us_firstDeriv,N);
	double L= 0;
	for(unsigned int i=0;i<ut_firstDeriv.size();i++) L+= std::abs(ut_firstDeriv[i]);	
	L*= 2.*TMath::Pi()/N;
		

	//Compute average bending energy
	const int NSCALES= 18;
	double SigmaScales[]= {
		1.e-7,5.e-7,1.e-6,5.e-6,1.e-5,5.e-5,1.e-4,5.e-4,1.e-3,5.e-3,1.e-2,5.e-2,0.1,0.5,1,2,5,10
	};
	BendingEnergies.clear();

	for(int k=0;k<NSCALES;k++){
		double smoothPar= SigmaScales[k];

		//Compute curvature at each contour point		
		std::vector<double> Curvature= Utils::GetContourCurvature(ut,smoothPar);

		//Compute average bending energy = L^2/N x sum(curv^2)
		double BendingEnergy= 0;
		for(unsigned int i=0;i<Curvature.size();i++){
			BendingEnergy+= pow(Curvature[i],2);
		}//end loop points
		BendingEnergy*= pow(L,2)/N;
		//BendingEnergy/= (double)N;
		BendingEnergies.push_back(BendingEnergy);
		cout<<"Contour::ComputeBendingEnergy(): INFO: Scale no. "<<k<<"="<<smoothPar<<" BE="<<BendingEnergy<<endl;
	}//end loop scales

}//close ComputeBendingEnergy()


