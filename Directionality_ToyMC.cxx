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

int PMTNumber = 17611;

using namespace std;

double GetFloatPrecision(double value, double precision){
    return (floor((value * pow(10, precision) + 0.5)) / pow(10, precision));
}

inline bool exists_test0 (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

void CartesianToSpherical(double & r,double & theta,double & phi,double x,double y,double z ){
	
	r = sqrt(x*x+y*y+z*z);
	if ( r != 0.0 ){
		phi = TMath::ACos( z / r );
		theta = TMath::ATan2( y, x );
		}
	else
	theta = phi = 0.0;
}

void SphericalToCartesian(double & x, double & y, double & z, double r, double theta, double phi){
	if ( r < 0.0 ){
		throw "Negative radius in sphericalToCartesian()";
		}
	x = r * TMath::Sin( phi ) * TMath::Cos( theta );
	y = r * TMath::Sin( phi ) * TMath::Sin( theta );
	z = r * TMath::Cos( phi );
}

double Distance(double x1, double y1, double z1, double x2, double y2, double z2){
	return sqrt( pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2));
	}
	
double DistanceOnASphere(double r, double theta1, double phi1, double theta2, double phi2){
	return TMath::ACos(TMath::Cos(phi1)*TMath::Cos(phi2) + (TMath::Sin(phi1)*TMath::Sin(phi2))*TMath::Cos(theta1-theta2));
	}

double ClosestPMTIndex(double x_Event,double y_Event,double z_Event, vector<vector<double>> & PMT_Position_Spherical){

	double MinDistance=50000;
	double Distance_Temp;
	int Index=50000;
	int Closest=Index;
	double r_Event, theta_Event, phi_Event;
	double x_PMT,y_PMT,z_PMT, r_PMT, theta_PMT, phi_PMT;
	
	CartesianToSpherical(r_Event, theta_Event, phi_Event,x_Event, y_Event, z_Event);
		
	for(int PMT=0;PMT<PMTNumber;PMT++){
		r_PMT = PMT_Position_Spherical[PMT][0];
		theta_PMT = PMT_Position_Spherical[PMT][1];
		phi_PMT = PMT_Position_Spherical[PMT][2];
		//CartesianToSpherical(r_PMT, theta_PMT, phi_PMT,x_PMT,y_PMT,z_PMT);
		Distance_Temp = DistanceOnASphere(r_PMT, theta_PMT, phi_PMT,theta_Event, phi_Event);
		//Distance_Temp = r_PMT*TMath::ACos(TMath::Cos(phi_PMT)*TMath::Cos(phi_Event) + (TMath::Sin(phi_PMT)*TMath::Sin(phi_Event))*TMath::Cos(theta_PMT-theta_Event));	
		//cout << PMT <<  " r " << r_PMT << "  th " <<   theta_PMT << "  phi " << phi_PMT <<  endl;
		
		if(Distance_Temp<MinDistance){
			Closest = PMT;
			MinDistance = Distance_Temp;
			}
		
		//cout << theta_PMT << "  "  << phi_PMT << "    "  << Index << "   " << Distance_Temp << "  /   " << Closest << "   " << MinDistance << endl;	// DO NOT REMOVE
			
	}
		
	return Closest;

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
	std::vector<vector<double>> PMT_Position_Spherical;		
		
	TH1F *h_Photon_Direction_x = new TH1F("h_Photon_Direction_x","h_Photon_Direction_x",500,-1,1);
	TH1F *h_Photon_Direction_y = new TH1F("h_Photon_Direction_y","h_Photon_Direction_y",500,-1,1);
	TH1F *h_Photon_Direction_z = new TH1F("h_Photon_Direction_z","h_Photon_Direction_z",500,-1,1);
	TH1F *h_Photon_Direction_r = new TH1F("h_Photon_Direction_r","h_Photon_Direction_r",500,0,22000);
	TH1F *h_Photon_Direction_theta = new TH1F("h_Photon_Direction_theta","h_Photon_Direction_theta",500,0,2*PI);	
	TH1F *h_Photon_Direction_phi = new TH1F("h_Photon_Direction_phi","h_Photon_Direction_phi",500,0,PI);			
	TH1F *h_ClosestIndex = new TH1F("h_ClosestIndex","h_ClosestIndex",1000,0,20000);	
	
	double RefractionIndex = 1.5;
	double TravelledDistance = 19.0;
	double AbsorptionLength = 80.;
	double QE = 0.3;
	double AbsorptionProbability = (1-TMath::Exp(-TravelledDistance/AbsorptionLength));
	double ReemissionProbability = 0.75;
	
	double SurvivingProbability = (1-AbsorptionProbability*(1-ReemissionProbability))*QE;
	
	cout << "AbsorptionProbability = " << AbsorptionProbability << endl;
	cout << "SurvivingProbability = " << SurvivingProbability << endl;
	
	// LOAD PMTs position
	ifstream ReadPMTPosition;
	ReadPMTPosition.open("PMTPos_CD_LPMT.csv");
	double blank;
	int Index;
	double x_PMT,y_PMT,z_PMT, r_PMT, theta_PMT, phi_PMT;
		
	for(int PMT=0;PMT<PMTNumber;PMT++){		// TO BE CHANGED this is extremely not efficient since we are loading and reading the PMTs position for each photon
		ReadPMTPosition >> Index;
		ReadPMTPosition >> x_PMT;
		ReadPMTPosition >> y_PMT;
		ReadPMTPosition >> z_PMT;
		ReadPMTPosition >> blank >> blank;
		CartesianToSpherical(r_PMT, theta_PMT, phi_PMT,x_PMT,y_PMT,z_PMT);
		PMT_Position_Spherical.push_back({r_PMT,theta_PMT,phi_PMT});
		//cout << "x " << x_PMT << "  y  " <<   y_PMT << "  z  " << z_PMT <<  endl;			
	}	
			
	for(int iPh=0; iPh<Photons; iPh++){

		double xx_at_PMTs, yy_at_PMTs, zz_at_PMTs;
		
		//Generate unit vector over a sphere
		double rr = 1;
		double theta = gRandom->TRandom::Uniform(2*PI);
		double phi = TMath::ACos(-1.+2.*gRandom->TRandom::Uniform(0,1));

		//Scintillation_Cartesian.push_back({xx,yy,zz});
		InteractionVertex.push_back({0,0,0});
		double xx,yy,zz;
		SphericalToCartesian(xx,yy,zz,rr,theta,phi);
		
		h_Photon_Direction_x->Fill(xx);
		h_Photon_Direction_y->Fill(yy);
		h_Photon_Direction_z->Fill(zz);
		h_Photon_Direction_r->Fill(rr);
		h_Photon_Direction_theta->Fill(theta);
		h_Photon_Direction_phi->Fill(phi);
		
		Scintillation_Spherical_atPMTs.push_back({TravelledDistance,theta,phi});
		SphericalToCartesian(xx_at_PMTs,yy_at_PMTs,zz_at_PMTs,TravelledDistance,theta,phi);	
		Scintillation_Cartesian_atPMTs.push_back({xx_at_PMTs,yy_at_PMTs,zz_at_PMTs});
		
		//cout << iPh << "\t" << xx_at_PMTs << "\t" << yy_at_PMTs << "\t" << zz_at_PMTs << endl;

		//int IndexExample = ClosestPMTIndex(Scintillation_Cartesian_atPMTs[iPh][0],Scintillation_Cartesian_atPMTs[iPh][1],Scintillation_Cartesian_atPMTs[iPh][2]);

		int IndexExample = ClosestPMTIndex(Scintillation_Cartesian_atPMTs[iPh][0],Scintillation_Cartesian_atPMTs[iPh][1],Scintillation_Cartesian_atPMTs[iPh][2],PMT_Position_Spherical);
			
		cout << "PHOTON " << iPh << " : CLOSEST INDEX IS = " << IndexExample << endl;
	
		h_ClosestIndex->Fill(IndexExample);
						
	}
	

	TFile *foutput = new TFile (Output_Rootfile.c_str(), "RECREATE");
	foutput->cd();
	h_Photon_Direction_x->Write();
	h_Photon_Direction_y->Write();
	h_Photon_Direction_z->Write();
	h_Photon_Direction_r->Write();
	h_Photon_Direction_theta->Write();
	h_Photon_Direction_phi->Write();
	h_ClosestIndex->Write();
	
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
