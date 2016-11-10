/*Include{{{*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <cmath>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <errno.h>
#include <sstream>
//#include <iterator>
//from root
#include <TMath.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TMatrix.h>
#include <TFile.h>
#include <TTree.h>
#include <TStopwatch.h>
#include <TRandom2.h>
/*}}}*/
//#include "christy_bosted_inelastic.h"
#include "DIS_Lite.h"
using namespace::std;
using namespace::TMath;
static const double PI=3.14159265358979323846;
static const double ALPHA=1/137.036; //fine-structure constant
static const double Deg2Rad = PI/180.0;
static const double MeV2GeV = 0.001;
static const double AP=ALPHA/PI;
static const double ELECTRON_MASS=0.510998918*MeV2GeV; //GeV
static const double NQe = 1.6*1e-19;
static const double Avogadro = 6.02*1e+23;
static const double M2=ELECTRON_MASS*ELECTRON_MASS*MeV2GeV; //GeV
static const double AMU=931.494043*MeV2GeV; //GeV
static const double PROTON_MASS=938.272029*MeV2GeV; //GeV
static const double DELTA_MASS=1232.*MeV2GeV; //GeV
static const double PI_MASS=139.57*MeV2GeV; //GeV
static const double HBARC=197.326968*MeV2GeV; //hbar*c GeV fm
static const double NBARN_TO_GEV2= 1./3.894 * 1e-5; //nbarn-->GeV^-2
static const double GEV2_TO_NBARN= 3.894 * 1e5; //GeV^-2-->nbarn
static const double NBARN_TO_CM2=1e-33; //nbarn-->CM^2


void PhaseSpace(){
	double E0 = 11.0;
	double Ep, Theta, xbj, nu, W2, eta, Q2;
	double xs_p,xs_n,xs_d;

	const double Ep_Min = 0.5; //GeV
	const double Ep_Max = 11.0; //GeV
	const double Theta_Min = 12.0; //Deg
	const double Theta_Max = 60.0; //Deg
	const int Nevnt = 1000000;

    DIS* dis = new DIS();
 
	/*Define Target Info{{{*/
	int A = 0, Z = 0;
	double Beam_Current = 0.0;//A
	//Upstream Window:   0.010" x 0.069g/cm-2
	//Downstream Window: 0.011" x 0.075g/cm-2
	const double Window_ArialDensity =  ( 0.069+0.075 );//g/cm-2, Al-2219 Alloy
	const double Window_Length = 25.;//cm
	double Target_Density = 0.0;//g*cm-3
	double Target_Thick = 0.0;//cm

	//H3
	A = 3; Z = 1;
	Beam_Current = 20*1e-6;//25uA
	Target_Thick = Window_Length;
	Target_Density = 0.081/Target_Thick;//g*cm-3
	//Target_Density = 0.00324;//g*cm-3
	double H3_Window_Lumi= Beam_Current/NQe * Window_ArialDensity*Avogadro/27;//cm-2*s-1
	double H3_Target_Lumi= Beam_Current/NQe * Target_Thick*Target_Density*Avogadro/A;//cm-2*s-1

	//He3
	A = 3; Z = 2;
	Beam_Current = 25*1e-6;//25uA
	Target_Thick = Window_Length;
	Target_Density = 0.093/Target_Thick;//g*cm-3
	double He3_Window_Lumi= Beam_Current/NQe * Window_ArialDensity*Avogadro/27;//cm-2*s-1
	double He3_Target_Lumi= Beam_Current/NQe * Target_Thick*Target_Density*Avogadro/A;//cm-2*s-1
	/*}}}*/

	/*Define ROOT file and tree{{{*/
	TString filename = Form("./rate_marathon.root");
	TFile *file = new TFile(filename,"recreate");
	TTree* T = new TTree("T","A New Tree");

	T->Branch("E0",&E0, "E0/D");
	T->Branch("Ep",&Ep, "Ep/D");
	T->Branch("Theta",&Theta, "Theta/D");
	T->Branch("Q2",&Q2, "Q2/D");
	T->Branch("W2",&W2, "W2/D");
	T->Branch("xbj",&xbj, "xbj/D");
	T->Branch("nu",&nu, "nu/D");
	T->Branch("eta",&eta, "eta/D");
	T->Branch("He3_Target_Lumi",&He3_Target_Lumi, "He3_Target_Lumi/D");
	T->Branch("He3_Window_Lumi",&He3_Window_Lumi, "He3_Window_Lumi/D");
	T->Branch("H3_Target_Lumi", &H3_Target_Lumi,  "H3_Target_Lumi/D");
	T->Branch("H3_Window_Lumi", &H3_Window_Lumi,  "H3_Window_Lumi/D");
	T->Branch("xs_p",&xs_p, "xs_p/D");
	T->Branch("xs_n",&xs_n, "xs_n/D");
	T->Branch("xs_d",&xs_d, "xs_d/D");
	/*}}}*/

	TRandom2 *fRan = new TRandom2();
	fRan->SetSeed(0); 

	for(int i=0;i<Nevnt;i++){  
		Ep = Ep_Min + (Ep_Max - Ep_Min) * fRan->Rndm();
		Theta = Theta_Min + (Theta_Max - Theta_Min) * fRan->Rndm();
		nu = E0-Ep;
		Q2 = 4.0*E0*Ep*pow(sin(Theta*Deg2Rad/2.0),2);
		xbj = Q2 / (2.0*PROTON_MASS*nu);
		W2 = ( PROTON_MASS*PROTON_MASS + 2.0*PROTON_MASS*nu - Q2 );
		double tao = Q2/4.0/PROTON_MASS/PROTON_MASS;
		eta = 1./(1.0 + 2.0*(1.0+tao)*pow(tan(Theta*Deg2Rad/2.0), 2));
    
        dis->SetKin(E0, Ep, Theta*Deg2Rad);

        xs_p = dis->Sigma("Proton");
        xs_n = dis->Sigma("Neutron");
        xs_d = dis->Sigma("Deutron");

        if(isnan(xs_p)||xs_p<-1e-20)
			xs_p=1e-66;
		if(isnan(xs_n)||xs_n<-1e-20)
			xs_n=1e-66;
		if(isnan(xs_d)||xs_d<-1e-20)
			xs_d=1e-66;

        T->Fill();
	}
	file->cd(); T->Write(); file->Close();
}
