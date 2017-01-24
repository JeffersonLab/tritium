////////////////////////////////////////
// Detection Efficiency Study
//     For x>2 experiment
// -- Zhihong Ye 09/01/2011
////////////////////////////////////////

#include <iostream>
#include <iomanip>
using namespace std;

void R_Det_Scan_2nd(){
  
  gStyle->SetOptStat(0);
 
  TString gDir = "/w/halla-2/e08014/Rootfiles/";
  
  TString GeneralCut = "R.tr.n==1 && abs(R.tr.x)<0.75 && abs(R.tr.y)<0.55 && abs(R.tr.th)<0.15 && abs(R.tr.ph)<0.045";
  // GeneralCut += "&& abs(R.tr.tg_ph)<0.04 && abs(R.tr.tg_y)<0.06 && abs(R.tr.tg_th)<0.06 && abs(R.tr.tg_dp)<0.05 && R.tr.p>2.0&&R.tr.p<4.0";
  TString TriggerCut1 = "((DBB.evtypebits>>1)&1)";
  TString TriggerCut6 = "((DBB.evtypebits>>6)&1)";
  
  TString electron_cut = "(R.sh.e+R.ps.e)/R.gold.p/1000.0>0.98&&(R.sh.e+R.ps.e)/R.gold.p/1000.0<1.02";
  
  TString aCut1 = TriggerCut1+"&&" + GeneralCut;
  TString aCut6 = TriggerCut6+"&&" + GeneralCut;
  
  ifstream runlist("R_2nd.txt");
  ofstream outfile("R_Det_Eff_2nd.dat");
  ofstream log("R_Det.Eff_2nd.log");

  Int_t run_number; TString goodr,goodl, comment;
  Float_t Cer_T1,Cer_T6,Calo_T1,Calo_T6;

  TString cut_cer1,cut_cer2,cut_cer3,cut_cer4;
  Float_t evt_cer1,evt_cer2,evt_cer3,evt_cer4;
  TString cut_calo1,cut_calo2,cut_calo3,cut_calo4;
  Float_t evt_calo1,evt_calo2,evt_calo3,evt_calo4;
  
  while(!runlist.eof()){
    
    runlist >> run_number;
    cerr << "Working on Run#"<<run_number<<endl;
    log << "Working on Run#"<<run_number<<endl;

    filename = gDir + Form("e08014_Rs_%d.root",run_number);
    
    TFile *f1 = new TFile(filename.Data());
    if (f1->IsZombie()){
      log << "----------> Bad Run#"<<run_number<<endl;
      cerr << "----------> Bad Run#"<<run_number<<endl;
      
    }
    
    else {
      
      TTree *T = (TTree*)gDirectory->Get("T");
      
      /////////////////////////////////////
      // Cerenkov Detection Eff
      /////////////////////////////////////
      cut_cer1  = aCut1 + "&&"+electron_cut;
      evt_cer1 = T->GetEntries(TCut(cut_cer1));
      
      cut_cer2 = "R.cer.asum_c>5.0 &&" + cut_cer1;
      evt_cer2 = T->GetEntries(TCut(cut_cer2));
      
      Cer_T1 = evt_cer2/evt_cer1*100.;
      
      cut_cer3 = aCut6 + "&&"+electron_cut;
      evt_cer3 = T->GetEntries(TCut(cut_cer3));
      
      cut_cer4 = "R.cer.asum_c>5.0 &&" + cut_cer3;
      evt_cer4 = T->GetEntries(TCut(cut_cer4));
      
      Cer_T6 = evt_cer4/evt_cer3*100.;
      cerr << "Cerenkov Detection Efficiency is: T1-->"<< Cer_T1<<"%" <<", T6-->"<<Cer_T6<<"%"<<endl;
      
      /////////////////////////////////////
      // Shower&Preshower Detection Eff
      /////////////////////////////////////
      cut_calo1 = aCut1 + "&&R.cer.asum_c>800&&R.cer.asum_c<1000.";
      evt_calo1 = T->GetEntries(TCut(cut_calo1));
      
      cut_calo2 ="(R.sh.e+R.ps.e)>10. &&" + cut_calo1;
      evt_calo2 = T->GetEntries(TCut(cut_calo2));
      Calo_T1 = evt_calo2/evt_calo1*100.;
      
      cut_calo3 = aCut6 + "&&R.cer.asum_c>800&&R.cer.asum_c<1000.";
      evt_calo3 = T->GetEntries(TCut(cut_calo3));
      
      cut_calo4 ="(R.sh.e+R.ps.e)>10. &&" + cut_calo3;
      evt_calo4 = T->GetEntries(TCut(cut_calo4));
      
      Calo_T6 = evt_calo4/evt_calo3*100.;
      cerr << "Shower Detection Efficiency is: T1-->"<< Calo_T1<<"%" <<", T6-->"<<Calo_T6<<"%"<<endl;
      outfile << setiosflags(ios::left) << setw(10) << setiosflags(ios::fixed) << setprecision(1) << run_number << " ";
      outfile << setiosflags(ios::left) << setw(12) << setiosflags(ios::fixed) << setprecision(4) << Cer_T1 << " ";
      outfile << setiosflags(ios::left) << setw(12) << setiosflags(ios::fixed) << setprecision(4) << Cer_T6 << " ";
      outfile << setiosflags(ios::left) << setw(12) << setiosflags(ios::fixed) << setprecision(4) << Calo_T1 << " ";
      outfile << setiosflags(ios::left) << setw(12) << setiosflags(ios::fixed) << setprecision(4) << Calo_T6 << endl;
	
	f1->Close();
    }//If(!(f1->IsZombie()))
      
  }//while(!runslit.eof())

  runlist.close();
  outfile.close();

}
