#include <TROOT.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TTree.h>
#include <TFile.h>
#include <TStyle.h>
#include <TH2D.h>

#include <MathUtils.h>
#include <AstroUtils.h>

using namespace Caesar;

//Vars
std::vector<double> LgFluxBins= {
	-5,-4.5,-4,-3.5,-3,-2.75,-2.5,-2.25,-2,-1.75,-1.5,-1.25,-1,-0.5,0,1
};

std::map<int,int> typeToHistoIdMap;


int DrawImagingPerformances(std::string fileName){

	//## Set logging level
	LoggerManager::Instance().CreateConsoleLogger("INFO","logger","System.out");
	
	//## Read file
	TFile* inputFile= new TFile(fileName.c_str(),"READ");
	if(!inputFile){
		ERROR_LOG("Failed to read file "<<fileName<<"!");
		return -1;		
	}

	//## Get tree
	TTree* data= (TTree*)inputFile->Get("data");
	if(!data){
		ERROR_LOG("Failed to read tree from file "<<fileName<<"!");
		return -1;
	}
	
	//# Set tree branches
	std::string name;
	int type;
	int simtype;
	double simmaxscale;
	double S;
	double S_true;
	double S_rec;
	double beamArea;

	//data->SetBranchAddress("name",&name);
	data->SetBranchAddress("type",&type);
	data->SetBranchAddress("simtype",&simtype);
	data->SetBranchAddress("simmaxscale",&simmaxscale);
	data->SetBranchAddress("S",&S);
	data->SetBranchAddress("S_true",&S_true);
	data->SetBranchAddress("S_rec",&S_rec);
	data->SetBranchAddress("beamArea",&beamArea);

	//## Define histos
	int nBins= (int)(LgFluxBins.size()-1);
	INFO_LOG("nBins="<<nBins);
	TH1D* dummyHisto= new TH1D("dummyHisto","dummyHisto",nBins,LgFluxBins.data());
	
	TGraph* fluxPoints_compact= new TGraph;
	TGraphAsymmErrors* fluxAccuracy_compact= new TGraphAsymmErrors;
	int nPoints_compact= 0;
		
	//## Init data vector
	std::vector<std::vector<double>> fluxList_compact;
	for(int i=0;i<nBins;i++){
		fluxList_compact.push_back( std::vector<double>() );
	}

	//## Loop over data and fill histos
	INFO_LOG("#"<<data->GetEntries()<<" sources to be read...");
	for(int i=0;i<data->GetEntries();i++){
		data->GetEntry(i);

		//double fluxDensity_true= S_true/beamArea;
		double fluxDensity_true= S/beamArea;
		double lgFlux_true= log10(fluxDensity_true);
		double fluxRatio= S_rec/S;

		int gBin= dummyHisto->FindBin(lgFlux_true);

		if(type==Source::eCompact || type==Source::ePointLike){
			fluxList_compact[gBin-1].push_back(fluxRatio);
			fluxPoints_compact->SetPoint(nPoints_compact,lgFlux_true,fluxRatio);	
			nPoints_compact++;
		}
		
	}//end loop sources

	//===============================================
	//==          FILL STATS HISTO
	//===============================================	
	//Compute data list stats
	nPoints_compact= 0;
	for(int i=0;i<nBins;i++){
		double x= LgFluxBins[i] + 0.5*(LgFluxBins[i+1]-LgFluxBins[i]);

		//Compute flux accuracy
		if(fluxList_compact[i].size()>=5){
			Caesar::BoxStats<double> stats= StatsUtils::ComputeBoxStats(fluxList_compact[i]);
			double ymin= *(std::min_element(fluxList_compact[i].begin(),fluxList_compact[i].end()));
			double ymax= *(std::max_element(fluxList_compact[i].begin(),fluxList_compact[i].end()));
			double y= stats.median;
			double yerr_low= y-stats.Q1;
			double yerr_up= stats.Q3-y;
			fluxAccuracy_compact->SetPoint(nPoints_compact,x,y);
			fluxAccuracy_compact->SetPointError(nPoints_compact,0,0,yerr_low,yerr_up);
			nPoints_compact++;
			INFO_LOG("Flux bin "<<i<<": x="<<x<<", y="<<y<<" (min/max="<<ymin<<"/"<<ymax<<")");
		}

	}//end loop bins

	//===============================================
	//==          DRAW POSITIONAL ACCURACY
	//===============================================	
	gROOT->SetStyle("myStyle2");
	double fluxRatio_min= -10;
	double fluxRatio_max= 10;

	TCanvas* FluxAccuracyPlot= new TCanvas("FluxAccuracyPlot","FluxAccuracyPlot");
	FluxAccuracyPlot->cd();

	TH2D* FluxAccuracyPlotBkg= new TH2D("FluxAccuracyPlotBkg","",100,LgFluxBins[0]-0.5,LgFluxBins[LgFluxBins.size()-1]+0.5,100,fluxRatio_min,fluxRatio_max);
	FluxAccuracyPlotBkg->GetXaxis()->SetTitle("lg(flux/Jy)");
	FluxAccuracyPlotBkg->GetYaxis()->SetTitle("flux_{rec}/flux");
	FluxAccuracyPlotBkg->SetStats(0);
	FluxAccuracyPlotBkg->Draw();

	//Compact sources
	fluxAccuracy_compact->SetMarkerSize(1.3);
	fluxAccuracy_compact->SetMarkerStyle(8);
	fluxAccuracy_compact->SetMarkerColor(kBlack);
	fluxAccuracy_compact->SetLineColor(kBlack);
	fluxAccuracy_compact->Draw("ep same");

	fluxPoints_compact->SetMarkerSize(1);
	fluxPoints_compact->SetMarkerStyle(1);
	fluxPoints_compact->SetMarkerColor(kBlack);
	fluxPoints_compact->SetLineColor(kBlack);
	fluxPoints_compact->Draw("P same");

	//refernce line
	TLine* refLine_fluxAccuracy= new TLine(LgFluxBins[0],1,LgFluxBins[LgFluxBins.size()-1],1);
	refLine_fluxAccuracy->SetLineColor(kBlack);
	refLine_fluxAccuracy->SetLineStyle(kDashed);
	refLine_fluxAccuracy->Draw("same");

	return 0;

}//close DrawImagingPerformances()
