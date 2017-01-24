////////////////////////////////////////
// Detection Efficiency Study
//     For x>2 experiment
// -- Zhihong Ye 09/01/2011
////////////////////////////////////////

#include <iostream>
#include <iomanip>
using namespace std;

void L_Det_Sample()
{
  
  gStyle->SetOptStat(0);
 
  //  TString gDir = "/GoFlex/Rootfiles/";
  TString gDir = "/GoFlex/Sample/";
 
  TString GeneralCut = "L.tr.n==1 && abs(L.tr.x)<0.75 && abs(L.tr.y)<0.55 && abs(L.tr.th)<0.15 && abs(L.tr.ph)<0.045";
  TString TriggerCut3 = "((DBB.evtypebits>>3)&1)";
  TString TriggerCut7 = "((DBB.evtypebits>>7)&1)";
  
  TString electron_cut = "(L.prl2.e+L.prl1.e)/L.gold.p/1000.0>0.98&&(L.prl2.e+L.prl1.e)/L.gold.p/1000.0<1.02";

  TString aCut3 = TriggerCut3+"&&" + GeneralCut;
  TString aCut7 = TriggerCut7+"&&" + GeneralCut;

  ifstream runlist("/w/halla-2/e08014/yez/HRS_Cali/runsampleL.txt");
  //  ifstream runlist("runL.txt");
  ofstream outfile("L_Det_Eff.dat");
    
  Int_t run_number, goodr,goodl, comment;
  Float_t Cer_T3,Cer_T7,Calo_T3,Calo_T7;

  TString cut_cer1,cut_cer2,cut_cer3,cut_cer4;
  Float_t evt_cer1,evt_cer2,evt_cer3,evt_cer4;
  TString cut_calo1,cut_calo2,cut_calo3,cut_calo4;
  Float_t evt_calo1,evt_calo2,evt_calo3,evt_calo4;

  TString kin,tgt;
  Float_t mom,angle;
  Int_t Rt1,Rt3,Rt6,Rt7;

  while(!runlist.eof())
    {
      runlist >> kin >> tgt >> run_number >> mom >> angle >> Rt1 >> Rt6 >> Rt3 >> Rt7;
      
      cerr << "Working on Run#"<<run_number<<endl;
      
      filename = gDir + Form("e08014_HRS_%d.root",run_number);
      
      TFile *f1 = new TFile(filename.Data());
      if (!(f1->IsZombie())){
	TTree *T;
	T = (TTree*)gDirectory->Get("T");
	
	/////////////////////////////////////
	// Cerenkov Detection Eff
	/////////////////////////////////////
	cut_cer1  = aCut3 + "&&"+electron_cut;
	evt_cer1 = T->GetEntries(TCut(cut_cer1));
	
	cut_cer2 = "L.cer.asum_c>1.0 &&" + cut_cer1;
	evt_cer2 = T->GetEntries(TCut(cut_cer2));
	
	Cer_T3 = evt_cer2/evt_cer1*100.;
	
	cut_cer3 = aCut7 + "&&"+electron_cut;
	evt_cer3 = T->GetEntries(TCut(cut_cer3));
	
	cut_cer4 = "L.cer.asum_c>10.0 &&" + cut_cer3;
	evt_cer4 = T->GetEntries(TCut(cut_cer4));
	
	Cer_T7 = evt_cer4/evt_cer3*100.;
	cerr << "Cerenkov Detection Efficiency is: T3-->"<< Cer_T3<<"%" <<", T7-->"<<Cer_T7<<"%"<<endl;
	
	/////////////////////////////////////
	// PRL1&PRL2 Detection Eff
	/////////////////////////////////////
	cut_calo1 = aCut3 + "&&L.cer.asum_c>400&&L.cer.asum_c<700.";
	evt_calo1 = T->GetEntries(TCut(cut_calo1));
	
	cut_calo2 ="(L.prl2.e+L.prl1.e)>10. &&" + cut_calo1;
	evt_calo2 = T->GetEntries(TCut(cut_calo2));
	Calo_T3 = evt_calo2/evt_calo1*100.;
	
	cut_calo3 = aCut7 + "&&L.cer.asum_c>600&&L.cer.asum_c<900.";
	evt_calo3 = T->GetEntries(TCut(cut_calo3));
	
	cut_calo4 ="(L.prl2.e+L.prl1.e)>1. &&" + cut_calo3;
	evt_calo4 = T->GetEntries(TCut(cut_calo4));
	
	Calo_T7 = evt_calo4/evt_calo3*100.;
	cerr << "PRL Detection Efficiency is: T3-->"<< Calo_T3<<"%" <<", T7-->"<<Calo_T7<<"%"<<endl;
	
	outfile << setiosflags(ios::left) << setw(10) << setiosflags(ios::fixed) << setprecision(1) << run_number << " ";
	outfile << setiosflags(ios::left) << setw(10) << setiosflags(ios::fixed) << setprecision(1) << kin << " ";
	outfile << setiosflags(ios::left) << setw(10) << setiosflags(ios::fixed) << setprecision(1) << tgt << " ";
	outfile << setiosflags(ios::left) << setw(10) << setiosflags(ios::fixed) << setprecision(1) << Rt3 << " ";
	outfile << setiosflags(ios::left) << setw(10) << setiosflags(ios::fixed) << setprecision(1) << Rt7 << " ";
	outfile << setiosflags(ios::left) << setw(12) << setiosflags(ios::fixed) << setprecision(4) << Cer_T3 << " ";
	outfile << setiosflags(ios::left) << setw(12) << setiosflags(ios::fixed) << setprecision(4) << Cer_T7 << " ";
	outfile << setiosflags(ios::left) << setw(12) << setiosflags(ios::fixed) << setprecision(4) << Calo_T3 << " ";
	outfile << setiosflags(ios::left) << setw(12) << setiosflags(ios::fixed) << setprecision(4) << Calo_T7 << " ";
	outfile << setiosflags(ios::left) << setw(10) << setiosflags(ios::fixed) << setprecision(1) << comment <<endl;
      
	f1->Close();
      }//If(!(f1->IsZombie()))
    }//while(!runslit.eof())
  runlist.close();
  outfile.close();
}
