/**
* @file SourceFitter.cc
* @class SourceFitter
* @brief SourceFitter
*
* Class to fit a source image with a mixture of gaussian/skew normal/skew-t bivariate functions
* @author S. Riggi
* @date 01/09/2015
*/

#include <SourceFitter.h>

#include <ConfigParser.h>
#include <ImgFITSReader.h>
#include <Img.h>
#include <Contour.h>
#include <Source.h>

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

ClassImp(SourceFitter)


SourceFitter::SourceFitter(){


}//close costructor


SourceFitter::~SourceFitter(){

}//close destructor


int SourceFitter::FitSource(Source* aSource,double curvThr,int componentMinNPix){

	if(!aSource){
		cerr<<"SourceFitter::FitSource(): WARN: Null ptr to source given!"<<endl;
		return -1;
	}

	//## Get source pars
	double Xmin= aSource->fXmin;
	double Xmax= aSource->fXmax;
	double Ymin= aSource->fYmin;
	double Ymax= aSource->fYmax;

	//## Get source and curvature image
	Img* sourceImg= aSource->GetImage(Source::eFluxMap);
	Img* curvImg= aSource->GetImage(Source::eFluxCurvMap);
	
	//## Find peaks in curvature map
	int tol= 1;
	TGraph* peaks= curvImg->FindPeaks(tol);
	if(!peaks) {
		cout<<"SourceFitter::FitSource(): INFO: No peaks detected...no deblending performed!"<<endl;
		sourceImg->Delete();
		curvImg->Delete();
		return -1;
	}
	int nPeaks= peaks->GetN();
	cout<<"SourceFitter::FitSource(): INFO: #"<<nPeaks<<" peaks detected!"<<endl;

	//## Find blobs image by thresholding the curvature map and then by flood-fill
	double fgValue= 1;
	Img* blobImg= curvImg->GetBinarized(curvThr);
	for(int k=0;k<nPeaks;k++){
		double x, y;
		peaks->GetPoint(k,x,y);
		int gBin= blobImg->FindBin(x,y);
		if(blobImg->IsBinOverflow(gBin) || blobImg->IsBinUnderflow(gBin) ) continue;
		blobImg->AddBinContent(gBin,1);
	}	
	blobImg->FindCompactSource(blobImg,fgValue+1,fgValue,componentMinNPix,false,false);
	std::vector<Source*> blobs= blobImg->GetSources();
		
	int nComponents= (int)blobs.size();
	if(nComponents<=0){
		cout<<"SourceFitter::FitSource(): INFO: No blobs detected...no deblending performed!"<<endl;
		sourceImg->Delete();
		curvImg->Delete();
		blobImg->Delete();
		return -1;	
	}
	cout<<"SourceFitter::FitSource(): INFO: #"<<nComponents<<" components detected!"<<endl;

	//## Get blob parameters
	Img* blobComponentImg= 0;
	std::vector<Img*> blobComponentImgList;
	std::vector< std::vector<double> > ComponentFitPars;
	std::vector< std::vector<double> > ComponentFitParsErr;
	std::string parNamePrefix[6]= {"A","meanX","meanY","sigmaX","sigmaY","theta"};
	TF2* componentFitFcn= new TF2("componentFitFcn",SourceFitter::Gaus2DGeomFcn,Xmin,Xmax,Ymin,Ymax,6);	
	for(int j=0;j<6;j++) componentFitFcn->SetParName(j,parNamePrefix[j].c_str());
	componentFitFcn->SetParLimits(6,0,180);
	componentFitFcn->SetNpx(1000);
	componentFitFcn->SetNpy(1000);

	TEllipse* ellipse= 0;
	std::vector<TEllipse*> ellipseList;
	double CLLevel= 0.683;//0.9;
	double quantile= ROOT::Math::chisquared_quantile(CLLevel,2);
	Source::FitInfo* fitInfo= 0;
	std::vector<Source::FitInfo*> fitInfoList;

	for(int k=0;k<nComponents;k++){
		//Get source mask
		cout<<"SourceFitter::FitSource(): INFO: Resetting sources in source image..."<<endl;
		sourceImg->ResetSources();
		cout<<"SourceFitter::FitSource(): INFO: Adding source no. "<<k<<" in source image..."<<endl;
		sourceImg->AddSource(blobs[k]);
		TString sourceName= Form("%s-Debl%d",sourceImg->GetName(),k);
		blobComponentImg= sourceImg->GetSourceMap(false);
		blobComponentImg->SetNameTitle(sourceName,sourceName);
		blobComponentImgList.push_back(blobComponentImg);
		
		double Smax= blobComponentImg->GetMaximum();
		double Smin= blobComponentImg->GetMinimum();
		double Sx= blobComponentImg->GetMean(1);
		double Sy= blobComponentImg->GetMean(2);
		double sigmaX= blobComponentImg->GetRMS(1);
		double sigmaY= blobComponentImg->GetRMS(2);
		double covXY= blobComponentImg->GetCovariance(1,2);
		double rho= covXY/(sigmaX*sigmaY);

		double theta= 0.5*atan((2*covXY)/(sigmaX*sigmaX-sigmaY*sigmaY));
		double ellMajAxis= sqrt(pow(sigmaX,2)*pow(sigmaY,2)*(1-rho*rho)/(pow(sigmaY*cos(theta),2)-2*covXY*cos(theta)*sin(theta)+pow(sigmaX*sin(theta),2)));
		double ellMinAxis= sqrt(pow(sigmaX,2)*pow(sigmaY,2)*(1-rho*rho)/(pow(sigmaY*sin(theta),2)+2*covXY*cos(theta)*sin(theta)+pow(sigmaX*cos(theta),2)));
		theta*= TMath::RadToDeg();

		cout<<"SourceFitter::FitSource(): INFO: Component no. "<<k<<": Smax="<<Smax<<" Mean("<<Sx<<","<<Sy<<"), CovMatrix("<<sigmaX<<","<<covXY<<","<<sigmaY<<") Ellipse("<<ellMajAxis<<","<<ellMinAxis<<", theta="<<theta<<")"<<endl;

		//Fit component
		componentFitFcn->SetParameters(Smax,Sx,Sy,sigmaX,sigmaY,theta);
		blobComponentImg->Fit(componentFitFcn,"RN");

		const double* fitPars= componentFitFcn->GetParameters();
		const double* fitParsErr= componentFitFcn->GetParErrors();
		ComponentFitPars.push_back( std::vector<double>() );
		ComponentFitParsErr.push_back( std::vector<double>() );
		cout<<"SourceFitter::FitSource(): INFO: Component no. "<<k<<": fitPars(";
		for(int l=0;l<componentFitFcn->GetNpar();l++){
			double parValue= *(fitPars+l);
			double parErrValue= *(fitParsErr+l);
			ComponentFitPars[k].push_back(parValue);	
			ComponentFitParsErr[k].push_back(parErrValue);
			cout<<parValue<<"+- "<<parErrValue<<",";
		}
		cout<<")"<<endl;
			
		double A_fit= ComponentFitPars[k][0];
		double Cx_fit= ComponentFitPars[k][1];
		double Cy_fit= ComponentFitPars[k][2];
		double sigmaX_fit= ComponentFitPars[k][3];
		double sigmaY_fit= ComponentFitPars[k][4];
		double theta_fit= ComponentFitPars[k][5];
		double AErr_fit= ComponentFitParsErr[k][0];
		double CxErr_fit= ComponentFitParsErr[k][1];
		double CyErr_fit= ComponentFitParsErr[k][2];
		double sigmaXErr_fit= ComponentFitParsErr[k][3];
		double sigmaYErr_fit= ComponentFitParsErr[k][4];
		double thetaErr_fit= ComponentFitParsErr[k][5];
		double covXY_fit= 0.5*tan(2*theta_fit*TMath::DegToRad())*(pow(sigmaX_fit,2)-pow(sigmaY_fit,2));
		double rho_fit= covXY_fit/(sigmaX_fit*sigmaY_fit);
		double ellMajAxis_fit= sqrt(pow(sigmaX_fit,2)*pow(sigmaY_fit,2)*(1-rho_fit*rho_fit)/(pow(sigmaY_fit*cos(theta_fit),2)-2*covXY_fit*cos(theta_fit)*sin(theta_fit)+pow(sigmaX_fit*sin(theta_fit),2)));
		double ellMinAxis_fit= sqrt(pow(sigmaX_fit,2)*pow(sigmaY_fit,2)*(1-rho_fit*rho_fit)/(pow(sigmaY_fit*sin(theta_fit),2)+2*covXY_fit*cos(theta_fit)*sin(theta_fit)+pow(sigmaX_fit*cos(theta_fit),2)));
		
		cout<<"SourceFitter::FitSource(): INFO: Component no. "<<k<<": Smax="<<Smax<<" Mean("<<Cx_fit<<","<<Cy_fit<<"), CovMatrix("<<sigmaX_fit<<","<<covXY_fit<<","<<sigmaY_fit<<") Ellipse("<<ellMajAxis_fit<<","<<ellMinAxis_fit<<", theta="<<theta_fit<<")"<<endl;
		
		ellipse= new TEllipse(Cx_fit,Cy_fit,ellMajAxis_fit,ellMinAxis_fit,0,360,theta_fit);
		ellipse->SetLineWidth(2);
		ellipse->SetFillColor(0);
		ellipse->SetFillStyle(0);
		ellipseList.push_back(ellipse);

		fitInfo= new Source::FitInfo;
		fitInfo->Cx= Cx_fit;
		fitInfo->Cy= Cy_fit;
		fitInfo->sigmaX= sigmaX_fit;
		fitInfo->sigmaY= sigmaY_fit;
		fitInfo->theta= theta_fit;
		fitInfo->CxErr= CxErr_fit;
		fitInfo->CyErr= CyErr_fit;
		fitInfo->sigmaXErr= sigmaXErr_fit;
		fitInfo->sigmaYErr= sigmaYErr_fit;
		fitInfo->thetaErr= thetaErr_fit;
		fitInfo->ellMaj= ellMajAxis_fit;
		fitInfo->ellMin= ellMinAxis_fit;
		fitInfoList.push_back(fitInfo);

		//if(blobComponentImg) blobComponentImg->Delete();
	}//end loop blobs


	//## Draw debug
	cout<<"SourceFitter::FitSource(): INFO: Draw debug plot..."<<endl;
	TF2* drawFcn= 0;
	std::vector<TF2*> drawFcnList;
	TCanvas* drawPlot;
	std::vector<TCanvas*> drawPlotList;	

	for(int k=0;k<ComponentFitPars.size();k++){
		TString canvasName= Form("drawPlot%d",k);
		drawPlot= new TCanvas(canvasName,canvasName);
		drawPlotList.push_back(drawPlot);

		drawPlotList[k]->cd();
		blobComponentImgList[k]->Draw("COLZ");
		ellipseList[k]->Draw("lsame");
		
		TString fcnName= Form("drawFcn%d",k+1);
		drawFcn= new TF2(fcnName,SourceFitter::Gaus2DGeomFcn,Xmin,Xmax,Ymin,Ymax,6);
		drawFcn->SetNpx(200);
		drawFcn->SetNpy(200);
		drawFcn->SetContour(100);
		for(int j=0;j<ComponentFitPars[k].size();j++) drawFcn->SetParameter(j,ComponentFitPars[k][j]);
		drawFcn->Draw("CONT1Z same");
		
		drawFcnList.push_back(drawFcn);	
	}//end loop components

	//## Fit component mixture
	if(nComponents>1){
		//Init mixture fit
		int parCounter= 0;
		TF2* mixtureFitFcn= new TF2("mixtureFitFcn",SourceFitter::Gaus2DMixtureGeomFcn,Xmin,Xmax,Ymin,Ymax,1+nComponents*6);
		mixtureFitFcn->SetNpx(1000);
		mixtureFitFcn->SetNpy(1000);
		mixtureFitFcn->SetParName(0,"nComponents");
		mixtureFitFcn->FixParameter(0,nComponents);
		parCounter++;
		for(int k=0;k<nComponents;k++) {
			for(int j=0;j<6;j++) {
				TString parName= Form("%s_%d",parNamePrefix[j].c_str(),k+1);
				mixtureFitFcn->SetParName(parCounter,parName);
				if(j==5) mixtureFitFcn->SetParLimits(parCounter,0,180);
				mixtureFitFcn->SetParameter(parCounter,ComponentFitPars[k][j]);
				parCounter++;
			}
		}

		//Fit mixture
		sourceImg->Fit(mixtureFitFcn,"RN");

		//Retrieve fit pars	
		parCounter= 1;
		TEllipse* ellipseFinal= 0;
		std::vector<TEllipse*> ellipseFinalList;
		for(int k=0;k<nComponents;k++){
			for(int j=0;j<6;j++) {
				double parValue= mixtureFitFcn->GetParameter(parCounter);
				double parValueErr= mixtureFitFcn->GetParameter(parCounter);
				ComponentFitPars[k][j]= parValue;
				ComponentFitParsErr[k][j]= parValueErr;
				parCounter++;
			}//end loop pars

			//Store pars
			double A_fit= ComponentFitPars[k][0];
			double Cx_fit= ComponentFitPars[k][1];
			double Cy_fit= ComponentFitPars[k][2];
			double sigmaX_fit= ComponentFitPars[k][3];
			double sigmaY_fit= ComponentFitPars[k][4];
			double theta_fit= ComponentFitPars[k][5];
			double AErr_fit= ComponentFitParsErr[k][0];
			double CxErr_fit= ComponentFitParsErr[k][1];
			double CyErr_fit= ComponentFitParsErr[k][2];
			double sigmaXErr_fit= ComponentFitParsErr[k][3];
			double sigmaYErr_fit= ComponentFitParsErr[k][4];
			double thetaErr_fit= ComponentFitParsErr[k][5];
			double covXY_fit= 0.5*tan(2*theta_fit*TMath::DegToRad())*(pow(sigmaX_fit,2)-pow(sigmaY_fit,2));
			double rho_fit= covXY_fit/(sigmaX_fit*sigmaY_fit);
			double ellMajAxis_fit= sqrt(pow(sigmaX_fit,2)*pow(sigmaY_fit,2)*(1-rho_fit*rho_fit)/(pow(sigmaY_fit*cos(theta_fit),2)-2*covXY_fit*cos(theta_fit)*sin(theta_fit)+pow(sigmaX_fit*sin(theta_fit),2)));
			double ellMinAxis_fit= sqrt(pow(sigmaX_fit,2)*pow(sigmaY_fit,2)*(1-rho_fit*rho_fit)/(pow(sigmaY_fit*sin(theta_fit),2)+2*covXY_fit*cos(theta_fit)*sin(theta_fit)+pow(sigmaX_fit*cos(theta_fit),2)));
			double Chi2= mixtureFitFcn->GetChisquare();
			double NDF= mixtureFitFcn->GetNDF();
			double NFreePars= mixtureFitFcn->GetNumberFreeParameters();
			
			ellipseFinal= new TEllipse(Cx_fit,Cy_fit,ellMajAxis_fit,ellMinAxis_fit,0,360,theta_fit);
			ellipseFinal->SetLineWidth(2);
			ellipseFinal->SetFillColor(0);
			ellipseFinal->SetFillStyle(0);
			ellipseFinalList.push_back(ellipseFinal);

			fitInfoList[k]->Cx= Cx_fit;
			fitInfoList[k]->Cy= Cy_fit;
			fitInfoList[k]->sigmaX= sigmaX_fit;
			fitInfoList[k]->sigmaY= sigmaY_fit;
			fitInfoList[k]->theta= theta_fit;
			fitInfoList[k]->CxErr= CxErr_fit;
			fitInfoList[k]->CyErr= CyErr_fit;
			fitInfoList[k]->sigmaXErr= sigmaXErr_fit;
			fitInfoList[k]->sigmaYErr= sigmaYErr_fit;
			fitInfoList[k]->thetaErr= thetaErr_fit;
			fitInfoList[k]->ellMaj= ellMajAxis_fit;
			fitInfoList[k]->ellMin= ellMinAxis_fit;

			aSource->AddFitInfo(fitInfoList[k]);
			aSource->fFitChi2= Chi2;
			aSource->fFitNDF= NDF;
			aSource->fFitNFreePars= NFreePars; 
		}//end loop components

		TCanvas* Plot= new TCanvas("Plot","Plot");
		Plot->cd();
		sourceImg->Draw("COLZ");
		mixtureFitFcn->Draw("CONT1Zsame");
		for(int k=0;k<nComponents;k++) ellipseFinalList[k]->Draw("lsame");
		//Clearup
		//if(mixtureFitFcn) mixtureFitFcn->Delete();
	}//close if

	//## Clearup
	//if(sourceImg) sourceImg->Delete();
	if(curvImg)	curvImg->Delete();
	if(blobImg)	blobImg->Delete();
	if(componentFitFcn) componentFitFcn->Delete();

	return 0;

}//close FitSource()


double SourceFitter::Gaus2DMixtureGeomFcn(double* x, double* p){

	int par_counter= 0;
	int nComponents= p[par_counter];
	par_counter++;
	
	int nPars= 6;
	double pdfPars[nPars];
	double fcnValue= 0.;
	
	for(int k=0;k<nComponents;k++){
		for(int j=0;j<nPars;j++){
			pdfPars[j]= p[par_counter];
			par_counter++;
		}		
		fcnValue+= Gaus2DGeomFcn(x,pdfPars);
	}//end loop mixtures

	return fcnValue;

}//close SourceFitter::Gaus2DMixtureGeomFcn()


double SourceFitter::Gaus2DGeomFcn(double* x, double* par){

	double X= x[0];
	double Y= x[1];
	
	double A= par[0];
	double X0= par[1];
	double Y0= par[2];
	double sigmaX= par[3];
	double sigmaY= par[4];
	double theta= par[5];
	theta*= TMath::DegToRad();

	double sigmaX2= sigmaX*sigmaX;
	double sigmaY2= sigmaY*sigmaY;

	double cosTheta= cos(theta);
	double cos2Theta= cos(2*theta);
	double cosTheta2= cosTheta*cosTheta;
	double sinTheta= sin(theta);
	double sin2Theta= sin(2*theta);
	double sinTheta2= sinTheta*sinTheta;

	double a= cosTheta2/(2*sigmaX2) + sinTheta2/(2*sigmaY2);
	//double b= -sin2Theta/(4*sigmaX2) + sin2Theta/(4*sigmaY2);
	double b= sin2Theta/(4*sigmaX2) - sin2Theta/(4*sigmaY2);
	double c= sinTheta2/(2*sigmaX2) + cosTheta2/(2*sigmaY2);

	double argX= a*(X-X0)*(X-X0);
	double argY= c*(Y-Y0)*(Y-Y0);
	double argXY= 2*b*(X-X0)*(Y-Y0);
	double z= argX+argXY+argY;

	double fcn= A*exp(-z);

	return fcn;

}//close SourceFitter::Gaus2DFcn()

/*
double SourceFitter::Gaus2DMixtureFcn(double* x, double* p){

	int par_counter= 0;
	double norm= p[par_counter];
	par_counter++;

	double w[fNComponents];
	for(int k=0;k<fNComponents;k++){
		w[k]= p[par_counter];
		par_counter++;	
	}
	
	double fcnValue= 0.;
	double pdfPars[6];
	pdfPars[0]= 1;

	for(int k=0;k<fNComponents;k++){
		pdfPars[1]= p[par_counter];
		par_counter++;
		pdfPars[2]= p[par_counter];
		par_counter++;
		pdfPars[3]= p[par_counter];
		par_counter++;		
		pdfPars[4]= p[par_counter];
		par_counter++;
		pdfPars[5]= p[par_counter];
		par_counter++;

		fcnValue+= w[k]*Gaus2DFcn(x,pdfPars);
	}//end loop mixtures

	fcnValue*= norm;

	return fcnValue;

}//close SourceFitter::Gaus2DMixtureFcn()


double SourceFitter::Gaus2DFcn(double* x, double* par){

	double X= x[0];
	double Y= x[1];
	
	double norm= par[0];
	double meanX= par[1];
	double meanY= par[2];
	double sigmaX= par[3];
	double sigmaY= par[4];
	double sigmaXY= par[5];
	double corrXY= sigmaXY/(sigmaX*sigmaY);
	double r= 1-corrXY*corrXY; 

	double normFactor= 1./(2*TMath::Pi()*sigmaX*sigmaY*sqrt(r));
	double argX= pow((X-meanX)/sigmaX,2);
	double argY= pow((Y-meanY)/sigmaY,2);
	double argXY= 2*corrXY*(X-meanX)*(Y-meanY)/(sigmaX*sigmaY);
	double z= argX-argXY+argY;

	double fcn= norm* exp(-z/(2*r));

	return fcn;

}//close SourceFitter::Gaus2DFcn()
*/




