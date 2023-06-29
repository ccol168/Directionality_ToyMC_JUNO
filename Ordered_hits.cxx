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

double mod (double x, double y, double z) {
    return sqrt(x*x + y*y + z*z);
}

double cos_alpha (double PhX, double PhY, double PhZ, double XVertex, double YVertex, double ZVertex) {
    double x = PhX - XVertex;
    double y = PhY - YVertex;
    double z = PhZ - ZVertex;

    return  z/(mod(x,y,z));
}

int main(int argc, char** argv) {
        
    if(argc!=4) {
            cout << "\n     USAGE:  Input_rootfile Output_rootfile Nth_hits \n" << endl;
            return 1;
    }

    string Input_rootfile = argv[1];
    string Output_Rootfile = argv[2];
    int Nth_hits = atoi(argv[3]);

    double Start_Time,PhX,PhY,PhZ,XVertex,YVertex,ZVertex;
    bool Type,Hit,IsFirst;
    int NEvent = 0;

    cout << "Reading " << Input_rootfile  << endl;
    cout << "Nth_hits = " << Nth_hits << endl << endl;

    TFile *f = new TFile (Input_rootfile.c_str());
    TTree* tree= (TTree*)f->Get("t");

    tree -> SetBranchAddress("Start_Time",&Start_Time);
    tree -> SetBranchAddress("Type",&Type);
    tree -> SetBranchAddress("Hit",&Hit);
    tree -> SetBranchAddress("IsFirst",&IsFirst);
    tree -> SetBranchAddress("NEvent",&NEvent);
    tree -> SetBranchAddress("Ph_x_AtPMT",&PhX);
    tree -> SetBranchAddress("Ph_y_AtPMT",&PhY);
    tree -> SetBranchAddress("Ph_z_AtPMT",&PhZ);
    tree -> SetBranchAddress("Int_Vertex_x",&XVertex);
    tree -> SetBranchAddress("Int_Vertex_y",&YVertex);
    tree -> SetBranchAddress("Int_Vertex_z",&ZVertex);

    int TotalPhotons = tree -> GetEntries();

    //cout << "Done" << endl;

    TH1F *ChScRatio = new TH1F("ChScRatio","ChScRatio",Nth_hits,0,Nth_hits);
    TH1F *Cherenkov_cos_alpha = new TH1F("Cherenkov_cos_alpha","Cherenkov_cos_alpha",60,-1,1);
    TH1F *Scint_cos_alpha = new TH1F("Scint_cos_alpha","Scint_cos_alpha",60,-1,1);

    //cout << "Done" << endl;
    //cout<<Nth_hits<<endl;
    TH1F *Nth_hit_cos_alpha[Nth_hits];

    for (int i=0;i<Nth_hits;i++) {
        Nth_hit_cos_alpha[i] = new TH1F(TString::Format("hit_%ith_cos_alpha",i),TString::Format("hit_%ith_cos_alpha",i),60,-1,1);
        //cout << "DONE " <<TString::Format("hit_%ith_cos_alpha",i) << endl;
    }

    //cout << "Done" << endl;

    int* FirstTen = new int[Nth_hits]; //contains the number of times a Cherenkov photon arrived in the n-th position
	double* FirstTenValues = new double[Nth_hits]; //contains the time at which the first n-th photon arrived (could be Cherenkov or scintillation)
	int* FirstTenPlaces = new int[Nth_hits]; //contains the place in the vector of the n-th photon to arrive (Cher or scint)

	for (int i=0; i<Nth_hits; i++) { //cleans the counter
		FirstTen[i] = 0;
	}

    //cout << "Done" << endl;

    int CurrentEvent = 0;

    for (int i=0;i<TotalPhotons;) {

        int j = 0;

        vector <double> val_cos_alpha;
        vector <bool> val_Type;
        vector <bool> val_Hit;

        for (int k=0; k<Nth_hits; k++) {

		    FirstTenValues[k] = 10000;
		    FirstTenPlaces[k] = -1; //useless, just to be sure that it is updating correctly into the cycle

            //if ((k==0 || k==Nth_hits-1) && i==0 )cout << "Done " << k << endl;
        }

        //cout << "In while " << CurrentEvent << "  " << NEvent << "  " << i << "  " <<TotalPhotons<< endl;
       
        while (CurrentEvent == NEvent && i<TotalPhotons)  {
    
            //cout << "Readtree " << i << endl;
            tree -> GetEntry(i);

            //cout<<"Done "<<i<<endl;

            //do what you do
            if (Type == 0 && Hit == 1) {
                Scint_cos_alpha -> Fill(cos_alpha(PhX,PhY,PhZ,XVertex,YVertex,ZVertex));
            } else if (Type==1 && Hit == 1) {
                Cherenkov_cos_alpha -> Fill(cos_alpha(PhX,PhY,PhZ,XVertex,YVertex,ZVertex));
            }

            val_cos_alpha.push_back(cos_alpha(PhX,PhY,PhZ,XVertex,YVertex,ZVertex));
            val_Type.push_back(Type);
            val_Hit.push_back(Hit);

            for (int k=0; k<Nth_hits; k++) {	//looks if a photon is in the first 10 to arrive and updates the arrays accordingly

					if (Start_Time < FirstTenValues[k] && Hit==1) {
					
						if (k != Nth_hits -1 ) {
							
							int* provvint = new int[Nth_hits-1];
							double* provvdouble = new double[Nth_hits-1];
							int cost =0;

							for (int h=k;h<Nth_hits-1;h++) {			
								provvint[cost] = FirstTenPlaces[h];
								provvdouble[cost] = FirstTenValues[h];
								cost ++;
								
							}
							cost=0;
							for (int h=k;h<Nth_hits-1;h++) {
								FirstTenPlaces[h+1] = provvint[cost];
								 FirstTenValues[h+1] = provvdouble[cost];
								cost ++;

							}
						}	
						

						FirstTenPlaces[k] = j;
						FirstTenValues[k] = Start_Time;
						break;
					}
				
				}
            
            i++;
            j++;

        } 

        //cout << "Out while" << endl;
        
        for (int k=0;k<Nth_hits;k++) { //updates the counter

			//cout << FirstTenValues[k]<<"  ";

			if ( val_Type[FirstTenPlaces[k]] == 1) {
				FirstTen[k]++;
			} 
            
            for (int h=0;h<Nth_hits;h++) {
                if (k == h) {
                    Nth_hit_cos_alpha[h] -> Fill(val_cos_alpha[FirstTenPlaces[k]]);
                }
            }

		}

        val_cos_alpha.clear();
        val_Hit.clear();
        val_Type.clear();


        //cout << "Event #"<<CurrentEvent<<" ended with " << j <<" photons, total number cycled " <<i<<"/"<<TotalPhotons<<endl;
        CurrentEvent++;
        
    }

    for (int k=0;k<Nth_hits;k++) {

		ChScRatio -> Fill(k,double(FirstTen[k])/(CurrentEvent+1));

    }

    cout << "Creating " << Output_Rootfile << endl; 

    TFile* out = new TFile(Output_Rootfile.c_str(),"recreate");

    ChScRatio -> Write();
    Cherenkov_cos_alpha -> Write();
    Scint_cos_alpha -> Write();
    ChScRatio -> Write();

    
    for (int i=0; i<Nth_hits;i++) {
        Nth_hit_cos_alpha[i] -> Write();
    }

    

    out -> Close();
    
    return 0;
}