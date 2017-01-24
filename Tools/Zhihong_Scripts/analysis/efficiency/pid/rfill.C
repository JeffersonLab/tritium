#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

void rfill(){

  TFile *outroot = new TFile("Detection.root","recreate");
  TTree *det = new TTree("det","Detection Tree");
  
  ifstream runL("L_Det_Eff.dat");
  ifstream runR("R_Det_Eff.dat");

  Int_t runnum,tgt,left,right; 
  TString kin,target;
  Double_t mom, angle;
  Double_t cer_t1 = -1.0;  Double_t cer_t6 = -1.0;
  Double_t cer_t3 = -1.0;  Double_t cer_t7 = -1.0;
  Double_t calo_t1 = -1.0;  Double_t calo_t6 = -1.0;
  Double_t calo_t3 = -1.0;  Double_t calo_t7 = -1.0;

  det->Branch("runnum",&runnum,"runnum/I");
  det->Branch("mom",&mom,"mom/D");
  det->Branch("tgt",&tgt,"tgt/I");
  det->Branch("left",&left,"left/I");
  det->Branch("right",&right,"right/I");
  det->Branch("cer_t1",&cer_t1,"cer_t1/D");
  det->Branch("cer_t3",&cer_t3,"cer_t3/D");
  det->Branch("cer_t6",&cer_t6,"cer_t6/D");
  det->Branch("cer_t7",&cer_t7,"cer_t7/D");
  det->Branch("calo_t1",&calo_t1,"calo_t1/D");
  det->Branch("calo_t3",&calo_t3,"calo_t3/D");
  det->Branch("calo_t6",&calo_t6,"calo_t6/D");
  det->Branch("calo_t7",&calo_t7,"calo_t7/D");

  while(!(runL.eof())){
    runL >> runnum >> kin >> target >> mom >> angle
	 >> cer_t3 >> cer_t7 >> calo_t3 >> calo_t7;
    left=1; right=0; 

    cerr << "Working on run#" << runnum <<endl;

    if(target == "BeO")
      tgt = 9;
    else if(target == "C12")
      tgt = 12;
    else if(target == "Ca40")
      tgt = 40;
    else if(target == "Ca48")
      tgt = 48;
    else if(target == "D20")
      tgt = 2720;
    else if(target == "D10")
      tgt = 2710;
    else if(target == "H2")
      tgt = 2;
    else if(target == "He3")
      tgt = 3;
    else if(target == "He4")
      tgt = 4;
    else if(target == "Optics")
      tgt = 55;
    else if(target == "BCM")
      tgt = 66;
    else
      tgt =-10;
    
    det->Fill();
  }
  
   while(!(runR.eof())){
   
     runR >> runnum >> kin >> target >> mom >> angle
	  >> cer_t1 >> cer_t6 >> calo_t1 >> calo_t6;

    cerr << "Working on run#" << runnum <<endl;

    left=0; right=1;

    if(target == "Be0")
      tgt = 9;
    else if(target == "C12")
      tgt = 12;
    else if(target == "Ca40")
      tgt = 40;
    else if(target == "Ca48")
      tgt = 48;
    else if(target == "D20")
      tgt = 2720;
    else if(target == "D10")
      tgt = 2710;
    else if(target == "H2")
      tgt = 2;
    else if(target == "He3")
      tgt = 3;
    else if(target == "He4")
      tgt = 4;
    else if(target == "Optics")
      tgt = 55;
    else if(target == "BCM")
      tgt = 66;
    else
      tgt =-10;

    det->Fill();
  }
  
  outroot->cd();
  det->Write();
  outroot->Close();
  runL.close(); runR.close();
}
