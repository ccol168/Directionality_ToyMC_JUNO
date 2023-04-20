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
double PMTRadius = 0.25;
double JUNORadius = 19.0;
double FV;
int LY;
double ChScRatio;
int NEvents;
int TotalPhotons = 0;

//Useful values

double n = 1.55 ; //refraction index
double c = 299792458 ; // m/s
double m_e = 0.51099895; //MeV    electron mass
double Be7_energy = 0.862; //MeV    energy of a 7Be neutrino
//double Event_Energy = 0.5; //MeV
double G_F = 1.1663787*pow(10,-11); //MeV^-2  Fermi constant
double sin2_thetaW = 0.22290; //Weinberg angle

double max_eEnergy = 2*pow(Be7_energy,2)/(m_e+2*Be7_energy); //MeV  maximum electron energy from a Be7-neutrino scattering

double RefractionIndex = 1.5;
double TravelledDistance = 19.0;
double AbsorptionLength = 80.;
double QE = 0.3;
double AbsorptionProbability = (1-TMath::Exp(-TravelledDistance/AbsorptionLength));
double ReemissionProbability = 0.75;
	
double SurvivingProbability = (1-AbsorptionProbability*(1-ReemissionProbability))*QE;

double El_Direction_x_t,El_Direction_y_t,El_Direction_z_t,phi_t,theta_t,Closest_PMT_t,Start_Time_t,Arr_Time_t,Electron_Energy_t,Neutrino_Energy_t,type_t;
double Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t, Ph_r_AtPMT_t, Ph_theta_AtPMT_t, Ph_phi_AtPMT_t, TravelledDistance_t;
double Int_Vertex_x_t,Int_Vertex_y_t,Int_Vertex_z_t;
bool Type_t, IsFirst_t;

double Min_Distance_t;
bool Hit_t;

double MaxPMTDistance; //maximum distance in phi between two PMTs

using namespace std;

struct Tuple {
	double x;
	double y;
};

//to keep phi and theta in their intended range
double Pbc_phi (double in) {
    if (in<0) return in + M_PI;
    else if (in>M_PI) return in - M_PI;
    else return in;
}

double Pbc_theta (double in) {
    if (in<0) return in + 2*M_PI;
    else if (in>2*M_PI) return in - 2*M_PI;
    else return in;
}

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
		theta = TMath::ATan2( y, x ) + PI;
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

double cross_section (double T) {
	double gl = 0.5+sin2_thetaW;
	double gr = sin2_thetaW;

	return ((2*(pow(G_F,2))*m_e)/Be7_energy)*(pow(gl,2) + pow(gr,2)*pow(1-T/Be7_energy,2) - gl*gr*m_e*T/pow(Be7_energy,2));
}

double CalculateEventEnergy () {
	double EventEnergy, test;
	bool flag = false;

	while (flag == false) {
		flag = false;
		EventEnergy = gRandom -> Uniform(0.,max_eEnergy);
		test = gRandom -> Uniform(0.,cross_section(0.)); //the maximum cross section is at T=0
		if (test < cross_section(EventEnergy)) {
			flag = true;
		}
	}

	return EventEnergy;
}

double ClosestPMTIndex (double x_Event,double y_Event,double z_Event, vector<vector<double>> & PMT_Position_Spherical) {

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
		Distance_Temp = DistanceOnASphere(r_PMT, theta_PMT, phi_PMT,theta_Event, phi_Event);
		//cout << PMT <<  " r " << r_PMT << "  th " <<   theta_PMT << "  phi " << phi_PMT <<  endl;
		
		if(Distance_Temp<MinDistance){
			Closest = PMT;
			MinDistance = Distance_Temp;
		} 
		if (MinDistance*JUNORadius <= PMTRadius) {

			Min_Distance_t = JUNORadius * MinDistance;

			return Closest;
		}

		//Code to speed up the process, should work but a check is due
		if (Distance_Temp > MaxPMTDistance*JUNORadius + 0.25) {
			PMT += 9;
		}
		//cout << theta_PMT << "  "  << phi_PMT << "    "  << Index << "   " << Distance_Temp << "  /   " << Closest << "   " << MinDistance << endl;	// DO NOT REMOVE
	}

	Min_Distance_t = JUNORadius * MinDistance;
	return -1;
}

//outputs the position (x,y,z) on the sphere for a photon generated in (j,k,l) parallel to the vector (a,b,c)

void MovePhoton (double & x,double& y , double& z, double j, double k, double l, double a, double b , double c ) {

	double sroot,num,denom,t;

	denom = 2*(a*a + b*b + c*c);
	sroot = pow(2*j*a + 2*k*b + 2*l*c,2) - 4*(a*a + b*b + c*c)*(j*j + k*k + l*l - JUNORadius*JUNORadius);
	num= -(2*j*a + 2*k*b + 2*l*c) + sqrt(sroot);

	t = num/denom;
	x = j+a*t;
	y = k+b*t;
	z = l+c*t;

	return;
}

Tuple Generate_Cone (double theta_0, double phi_0, double angle) {

	Tuple out;
    double phi_out, theta_out;
    double theta, phi, cos_th1, cos_ph1, sin_th1, sin_ph_sin_th, sin_ph_cos_th;

	//generate random values on a cone centered in 0,0
    theta = gRandom->TRandom::Uniform(2*PI); 
    phi = angle; 

	//some exceptions
    if (phi_0 == 0) {

        phi_out = phi;
        theta_out = theta;
        

    } else if (phi_0 == M_PI) {

        phi_out = PI - phi;
        theta_out = theta;

	// traslate the cone on the original system
    } else {

        cos_ph1 = -sin(phi_0)*cos(theta)*sin(phi) + cos(phi_0)*cos(phi);
        sin_ph_cos_th = cos(phi_0)*cos(theta_0)*cos(theta)*sin(phi) - sin(theta_0)*sin(theta)*sin(phi) + cos(theta_0)*sin(phi_0)*cos(phi);
        sin_ph_sin_th = sin(theta_0)*cos(phi_0)*cos(theta)*sin(phi) + cos(theta_0)*sin(theta)*sin(phi) + sin(theta_0)*sin(phi_0)*cos(phi);        

        phi_out = acos(cos_ph1);
        cos_th1 = sin_ph_cos_th/sin(phi_out);
        sin_th1 = sin_ph_sin_th/sin(phi_out);

        if (sin_th1 >= 0) {
            theta_out = acos(cos_th1) ;
        } else if (sin_th1 < 0) {
            theta_out = Pbc_theta(-acos(cos_th1));
        }

    }

	out.x = theta_out;
	out.y = phi_out;

	return out;
}

double GenerateScintStartTime () {

	double StartTime;

	double Check = gRandom -> Uniform(0,1);

	if (Check < 0.707) {
		StartTime = gRandom -> Exp(4.6*pow(10,-9));
	} else if (Check < 0.912) {
		StartTime = gRandom -> Exp(15.1*pow(10,-9));
	} else if (Check < 0.972) {
		StartTime = gRandom -> Exp(76.1*pow(10,-9));
	} else {
		StartTime = gRandom -> Exp(397*pow(10,-9));
	}

	return StartTime;
}

int CheckHit (ofstream& WriteOutputText, int SeenPhotons) {

	if (Closest_PMT_t != -1) {

		WriteOutputText << theta_t << "  " << phi_t << "   " << Closest_PMT_t << "  " << Start_Time_t << "  " << Arr_Time_t << "  " << Type_t << endl;
		Hit_t = 1;
		SeenPhotons++;

	} else {Hit_t = 0;}


	return SeenPhotons;
}

//Generator for all the photons

int GeneratePhotons (ofstream& WriteOutputText, TTree* t, vector<vector<double>> PMT_Position_Spherical, bool RandomPos) {

	double x_Int,y_Int,z_Int,r_Int,theta_Int,phi_Int;
	double theta_vers,phi_vers,trash;
	int SeenPhotons = 0;

	//Generate a random interaction vertex

	if (RandomPos == true ) {
		r_Int = gRandom -> TRandom::Uniform(FV);
		theta_Int = gRandom -> TRandom::Uniform(2*PI);
		phi_Int = TMath::ACos(-1.+2.*gRandom->TRandom::Uniform(0,1));
		SphericalToCartesian(x_Int,y_Int,z_Int,r_Int,theta_Int,phi_Int);
	} else {
		x_Int = 0;
		y_Int = 0;
		z_Int = 0;
	}

	//Randomize the energy of the event
	double Event_Energy = CalculateEventEnergy();

	double theta_e = acos((1+m_e/Be7_energy)*pow(Event_Energy/(Event_Energy+2*m_e),0.5)); //angle between the solar-nu and the electron scattered (assuming 7Be-nu)
	double beta_el = pow(1-(pow(m_e/(Event_Energy+m_e),2)),0.5) ; //beta of the electron generated
	double theta_Cher = acos(1/(beta_el*n)); //Cherenkov angle

	int Photons = LY*Event_Energy*SurvivingProbability;
	int CherenkovPhotons = ChScRatio*Photons;

	//cout<<Photons<<"  "<<CherenkovPhotons<<"  "<<Event_Energy<<endl;
	//cout<<max_eEnergy<<endl;

	//Generate SCINTILLATION Photons
	
	for(int iPh=0; iPh<Photons; iPh++){

		//Generate unit vector over a sphere
		theta_vers = gRandom->TRandom::Uniform(2*PI);
		phi_vers = TMath::ACos(-1.+2.*gRandom->TRandom::Uniform(0,1));
		SphericalToCartesian(El_Direction_x_t,El_Direction_y_t,El_Direction_z_t,1,theta_vers,phi_vers);

		MovePhoton(Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t,x_Int,y_Int,z_Int,El_Direction_x_t,El_Direction_y_t,El_Direction_z_t);

		Ph_r_AtPMT_t = JUNORadius; // NOT REALLY TRUE
		Ph_theta_AtPMT_t = theta_t; // TRUE ONLY WITHOUT ABSORPTION
		Ph_phi_AtPMT_t = phi_t; // TRUE ONLY WITHOUT ABSORPTION
		
		CartesianToSpherical(trash,theta_t,phi_t,Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t);

		Closest_PMT_t = ClosestPMTIndex(Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t,PMT_Position_Spherical);

		TravelledDistance_t = Distance(x_Int,y_Int,z_Int,Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t);
		Type_t = 0; 
		Electron_Energy_t = Event_Energy;
		Neutrino_Energy_t = Be7_energy;

		Int_Vertex_x_t = x_Int;
		Int_Vertex_y_t = y_Int;
		Int_Vertex_z_t = z_Int;
		
		Start_Time_t = GenerateScintStartTime();
		Arr_Time_t = Start_Time_t + TravelledDistance_t/(n*c);

		if (iPh == 0) IsFirst_t = true;
		else IsFirst_t = false;
		
		//cout << iPh << "\t" << xx_at_PMTs << "\t" << yy_at_PMTs << "\t" << zz_at_PMTs << endl;
			
		//cout << "PHOTON " << iPh << " : CLOSEST INDEX IS = " << IndexExample << endl;

		//Generate the photon only if it hits the PMT
		SeenPhotons = CheckHit(WriteOutputText,SeenPhotons);

		TotalPhotons++;

		t -> Fill();
		
		//cout << i << "    " << iPh % (Photons/10) << endl;
						
	}

	//Generate CHERENKOV Photons

	Tuple a;
	Tuple b;

	a = Generate_Cone(0.,0.,theta_e);

	for (int iPh=0; iPh<CherenkovPhotons; iPh++) {
		
		b = Generate_Cone(a.x,a.y,theta_Cher);

		theta_vers = b.x;
		phi_vers = b.y;
		SphericalToCartesian(El_Direction_x_t,El_Direction_y_t,El_Direction_z_t,1,theta_vers,phi_vers);

		MovePhoton(Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t,x_Int,y_Int,z_Int,El_Direction_x_t,El_Direction_y_t,El_Direction_z_t);

		Ph_r_AtPMT_t = JUNORadius; // NOT REALLY TRUE
		Ph_theta_AtPMT_t = theta_t; // TRUE ONLY WITHOUT ABSORPTION
		Ph_phi_AtPMT_t = phi_t; // TRUE ONLY WITHOUT ABSORPTION
		
		CartesianToSpherical(trash,theta_t,phi_t,Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t);

		TravelledDistance_t = Distance(x_Int,y_Int,z_Int,Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t);

		Closest_PMT_t = ClosestPMTIndex(Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t,PMT_Position_Spherical);

		Type_t = 1;
		Electron_Energy_t = Event_Energy;
		Neutrino_Energy_t = Be7_energy;

		Int_Vertex_x_t = 0.;
		Int_Vertex_y_t = 0.;
		Int_Vertex_z_t = 0.;

		Start_Time_t = 0.;
		Arr_Time_t = Start_Time_t + TravelledDistance_t/(n*c);

		IsFirst_t = false;

		//Generate the photon only if it hits the PMT
		CheckHit(WriteOutputText,SeenPhotons);

		TotalPhotons++;

		t -> Fill();
	}

	return SeenPhotons;
}


double Directionality_ToyMC(string Configuration_Text, string Output_Rootfile, string Output_Text) {

	ifstream file(Configuration_Text);
	vector<string> col1;
	vector<string> col2;
	string line;

	// ### Parsing
	cout << "######### Configuration #########" << endl;
	cout << "Cfg file: " << Configuration_Text.c_str() << endl;
	while (getline(file, line)) {
		string s1, s2;
		istringstream iss(line);
		if (!(iss >> s1 >> s2)) { break; } // error
			col1.push_back(s1);
			col2.push_back(s2);
			cout << setw(10) << s1 <<  "\t\t" << s2 << endl;
	}
	
	cout << "################################" << endl;

	// quite rough; to be changed
	istringstream iss(col2[0]);
	iss >> LY;
	istringstream iss1(col2[1]);
	iss1 >> ChScRatio;
	istringstream iss2(col2[2]);	
	iss2 >> NEvents;
	istringstream iss3(col2[3]);	
	iss3 >> FV;	
	// ### End parsing	
	
	gRandom = new TRandom3(0);
	gRandom->SetSeed(0);
	
	std::vector<vector<double>> PMT_Position_Spherical;		
		
	std::cout << "AbsorptionProbability = " << AbsorptionProbability << endl;
	cout << "SurvivingProbability = " << SurvivingProbability << endl;
	//cout<<"Neutrino-electron angle = "<<theta_e*180./M_PI<<" deg"<<endl;
	//cout<<"Cherenkov angle = "<<theta_Cher*180./M_PI<<" deg"<<endl;

	cout <<"Scintillation photons generated @ 1 MeV = " << int(LY*SurvivingProbability) << endl;
	cout <<"Number of events generated = " << NEvents << endl << endl;

	//Making the tree
	TTree *t = new TTree("t","New Tree");
	t->Branch("El_Direction_x", &El_Direction_x_t, "El_Direction_x/D");
	t->Branch("El_Direction_y", &El_Direction_y_t, "El_Direction_y/D");
	t->Branch("El_Direction_z", &El_Direction_z_t, "El_Direction_z/D");
	t->Branch("Ph_x_AtPMT", &Ph_x_AtPMT_t, "Ph_x_AtPMT/D");
	t->Branch("Ph_y_AtPMT", &Ph_y_AtPMT_t, "Ph_y_AtPMT/D");
	t->Branch("Ph_z_AtPMT", &Ph_z_AtPMT_t, "Ph_z_AtPMT/D");
	t->Branch("TravelledDistance", &TravelledDistance_t, "TravelledDistance/D");
	t->Branch("Ph_theta_AtPMT", &Ph_theta_AtPMT_t, "Ph_theta_AtPMT/D");
	t->Branch("Ph_phi_AtPMT", &Ph_phi_AtPMT_t, "Ph_phi_AtPMT/D");
	t->Branch("Closest_PMT", &Closest_PMT_t, "Closest_PMT/D");
	t->Branch("Start_Time", &Start_Time_t, "Start_Time/D");
	t->Branch("Arr_Time", &Arr_Time_t, "Arr_Time/D");
	t->Branch("Electron_Energy", &Electron_Energy_t, "Electron_Energy/D");
	t->Branch("Int_Vertex_x", &Int_Vertex_x_t, "Int_Vertex_x/D");
	t->Branch("Int_Vertex_y", &Int_Vertex_y_t, "Int_Vertex_y/D");
	t->Branch("Int_Vertex_z", &Int_Vertex_z_t, "Int_Vertex_z/D");
	t->Branch("Neutrino_Energy", &Neutrino_Energy_t, "Neutrino_Energy/D");
	t->Branch("Type", &Type_t, "Type/O");

	t->Branch("Min_Distance", &Min_Distance_t, "Min_Distance/D");
	t->Branch("Hit", &Hit_t, "Hit/O");
	t->Branch("IsFirst", &IsFirst_t, "IsFirst/O"); //first photon generated (useful if you want to see only the property of the scattered electron)

	// LOAD PMTs position
	ifstream ReadPMTPosition;
	ReadPMTPosition.open("PMTPos_CD_LPMT.csv");
	double blank;
	int Index;
	double x_PMT,y_PMT,z_PMT, r_PMT, theta_PMT, phi_PMT;

	for(int PMT=0;PMT<PMTNumber;PMT++){		
		ReadPMTPosition >> Index;
		ReadPMTPosition >> x_PMT;
		ReadPMTPosition >> y_PMT;
		ReadPMTPosition >> z_PMT;
		ReadPMTPosition >> blank >> blank;
		CartesianToSpherical(r_PMT, theta_PMT, phi_PMT,x_PMT,y_PMT,z_PMT);
		PMT_Position_Spherical.push_back({r_PMT,theta_PMT,phi_PMT});			
	}	

	// #################### PMT_Position_Sperical is already sorted in phi ####################################
	//sorting PMT_Position_Spherical in phi
	//std::sort(PMT_Position_Spherical.begin(),PMT_Position_Spherical.end(),[] (const std::vector<int>& a, const std::vector<int>& b) {return a[2] < b[2];});

	/*ofstream WriteTry;
	WriteTry.open("Prova.txt");
	cout << "##########################" << endl;
	cout << "Sorted PMT_Position_Spherical" << endl;
	for(int PMT=0;PMT<PMTNumber;PMT++){	
		WriteTry<<PMT_Position_Spherical[PMT][2]<<endl;
	}
	cout << "##########################" << endl;*/
	

	//Find MaxPMTDistance
	MaxPMTDistance = 0;
	double ProvPMTDistance;

	for(int PMT=0;PMT<PMTNumber-1;PMT++){		
		ProvPMTDistance = PMT_Position_Spherical[PMT+1][2] - PMT_Position_Spherical[PMT][2];
		if (ProvPMTDistance > MaxPMTDistance) {
			MaxPMTDistance = ProvPMTDistance;
		}	
	}	
	//cout << "MaxPMTDistance = " << MaxPMTDistance<<endl;

	ofstream WriteOutputText;
	WriteOutputText.open(Output_Text.c_str(),ios::app);	

	TFile *foutput = new TFile (Output_Rootfile.c_str(), "RECREATE");
	foutput->cd();

	int SeenPhotons = 0;

	for (int i=0; i<NEvents; i++) {
		SeenPhotons += GeneratePhotons(WriteOutputText, t, PMT_Position_Spherical, true);
		if (NEvents > 10) {  //to avoid floating point exceptions for NEvents < 10
			if (i % (NEvents/10) == 0 && i != 0 && NEvents > 10) { // check if the index is a multiple of tenth
			std::cout << i << "-th Event ; " << (i / (NEvents/10)) * 10 << "% of events simulated \n";
    		}
		}
		
	}

	t->Write();
	
	foutput->Close();

   // APPENDING text output to Output_Text text file
	WriteOutputText.close();

	cout << "Geometric coverage = " << double(SeenPhotons)/double(TotalPhotons) <<endl;
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
