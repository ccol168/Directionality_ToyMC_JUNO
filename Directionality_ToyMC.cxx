#include <fstream>
#include <string>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <iomanip>

#include <TApplication.h>
#include <TVector3.h>
#include <TFile.h>
#include <TChain.h>
#include <TH1.h>
#include "TMath.h"
#include <TF1.h>
#include "TH1D.h"
#include "TH1.h"
#include "TH2F.h"
#include "TH3F.h"
#include <TROOT.h>
#include "TStyle.h"
#include <TCanvas.h>
#include <TRandom3.h>
#include <TSpectrum.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TVirtualFitter.h>

#include <TRandom3.h>
#include <TFractionFitter.h>
#include <vector>

#include <TFitResult.h>
#include <TMatrixDSym.h>
#include "TMinuit.h"

#include <TString.h>

#ifndef ROOT_TH1D
#endif

#define PI 3.141592653589793
#include "TMath.h"
#include "Math/Util.h"

using namespace std;

int N = 12;
TH1D** hPDF = new TH1D*[N];
TH1D** hPDF_Gen = new TH1D*[N];

TH1D** hPDF_Nig = new TH1D*[N];
TH1D** hPDF_Sub = new TH1D*[N];

TH1D *PseudoDataset_Nig;
TH1D *PseudoDataset_Day;
TH1D *PseudoDataset_Diff;

int Beg, End;

double density = 859.0; // [kg/m3]
double radius_FV = 14.0;
double Mass_FV = 4./3.*PI*pow(radius_FV,3)*density; // [m3]
double *Threshold = new double[N];
double *Events = new double[N];
double *Events_Nig = new double[N];
double *Events_Day = new double[N];

double *InjRate = new double[N];
double *InjRate_Nig = new double[N];
double *InjRate_Day = new double[N];
double *InjRate_AfterSampling = new double[N];
double *Fix = new double[N]();
double *Sigma = new double[N]();

double days;
double Exposure_LT;
double Efficiency;
double Exposure_LT_Nig;
double Exposure_LT_Day;
double Exposure_LT_Diff;

double GetFloatPrecision(double value, double precision){
    return (floor((value * pow(10, precision) + 0.5)) / pow(10, precision));
}

inline bool exists_test0 (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

void CartesianToSpherical(double & r,double & theta,double & phi,double x,double y,double z ){
	if ( r != 0.0 ){
		theta = TMath::ACos( z / r );
		phi = TMath::ATan2( y, x );
		}
	else
	theta = phi = 0.0;
}

void SphericalToCartesian(double & x, double & y, double & z, double r, double theta, double phi){
	if ( r < 0.0 ){
		throw "Negative radius in sphericalToCartesian()";
		}
	x = r * TMath::Sin( theta ) * TMath::Cos( phi );
	y = r * TMath::Sin( theta ) * TMath::Sin( phi );
	z = r * TMath::Cos( theta );
}


double Directionality_ToyMC(string Configuration_Text, string Output_Rootfile, string Output_Text){

	ifstream ReadCfgFile;
	ReadCfgFile.open(Configuration_Text.c_str());

	int LY = 0;
	cout << "#############" << endl;
	cout << "Cfg file: " << Configuration_Text.c_str() << endl;
	
	ReadCfgFile >> LY;
	cout << LY << endl;

	gRandom = new TRandom3(0);
	gRandom->SetSeed(0);
	
	int Photons = LY;
	
	std::vector<vector<double>> Scintillation_Cartesian;
	std::vector<vector<double>> InteractionVertex;
	std::vector<vector<double>> Scintillation_Cartesian_atPMTs;
	
	std::vector<vector<double>> Scintillation_Spherical_atPMTs;
		
	TH1F *h_Photon_Direction_x = new TH1F("h_Photon_Direction_x","h_Photon_Direction_x",500,-1,1);
	TH1F *h_Photon_Direction_y = new TH1F("h_Photon_Direction_y","h_Photon_Direction_y",500,-1,1);
	TH1F *h_Photon_Direction_z = new TH1F("h_Photon_Direction_z","h_Photon_Direction_z",500,-1,1);
	TH1F *h_Photon_Direction_r = new TH1F("h_Photon_Direction_r","h_Photon_Direction_r",500,-1,1);	

	double RefractionIndex = 1.5;
	double TravelledDistance = 17.5;
	double AbsorptionLength = 80.;
	double QE = 0.3;
	double AbsorptionProbability = (1-TMath::Exp(-TravelledDistance/AbsorptionLength));
	double ReemissionProbability = 0.75;
	
	double SurvivingProbability = (1-AbsorptionProbability*(1-ReemissionProbability))*QE;
	
	cout << "AbsorptionProbability = " << AbsorptionProbability << endl;
	cout << "SurvivingProbability = " << SurvivingProbability << endl;
			
	for(int i=0; i<Photons; i++){
		double xx = -1+2.*gRandom->TRandom::Uniform(1);
		double yy = -1+2.*gRandom->TRandom::Uniform(1);
		double zz = -1+2.*gRandom->TRandom::Uniform(1);
		double xx_at_PMTs, yy_at_PMTs, zz_at_PMTs;
		double theta;
		double phi;
		double norm = sqrt(xx*xx+yy*yy+zz*zz);
		//xx = xx/norm;
		//yy = yy/norm;
		//zz = zz/norm;
		Scintillation_Cartesian.push_back({xx,yy,zz});
		InteractionVertex.push_back({0,0,0});
		h_Photon_Direction_x->Fill(xx);
		h_Photon_Direction_y->Fill(yy);
		h_Photon_Direction_z->Fill(zz);
		
		CartesianToSpherical(norm,theta,phi,xx,yy,zz);
		
		cout << norm << " " << theta << " " << phi << " " << xx << " " << yy << " " << zz << endl;
		
		Scintillation_Spherical_atPMTs.push_back({TravelledDistance,theta,phi});
		
		SphericalToCartesian(xx_at_PMTs,yy_at_PMTs,zz_at_PMTs,TravelledDistance,theta,phi);
		
		Scintillation_Cartesian_atPMTs.push_back({xx_at_PMTs,yy_at_PMTs,zz_at_PMTs});
		
		cout << i << "\t" << xx_at_PMTs << "\t" << yy_at_PMTs << "\t" << zz_at_PMTs << endl;
						
	}

	TFile *foutput = new TFile (Output_Rootfile.c_str(), "RECREATE");
	foutput->cd();
	h_Photon_Direction_x->Write();
	h_Photon_Direction_y->Write();
	h_Photon_Direction_z->Write();
	
	foutput->Close();

	ofstream WriteOutputText;
	WriteOutputText.open(Output_Text.c_str(),ios::app);		   // APPENDING text output to Output_Text text file
	WriteOutputText.close();

	cout << "#############" << endl;
	
	return 0;

}

int main(int argc, char** argv) {
        string macro = argv[0];

        if(argc!=4) {
                cout << "\n     USAGE:  "<< macro << " Configuration_Text  Output_Rootfile  Output_Text \n" << endl;
                return 1;
        }

        string Configuration_Text = argv[1];
        string Output_Rootfile = argv[2];
	      string Output_Text = argv[3];

        Directionality_ToyMC(Configuration_Text, Output_Rootfile, Output_Text);
        return 0;
}
