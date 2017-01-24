////////////////////////////////////////
// Scintillator Efficiency Study
// For x>2 experiment
// -- Zhihong Ye 07/08/2011
////////////////////////////////////////

#include <iostream>
#include <iomanip>
using namespace std;

void L_scin_eff(){

  Int_t  run_number;    TString filename;
  TString gDir = "/w/halla-2/e08014/Rootfiles_0721/";
  //TString gDir = "/w/halla-2/e08014/Rootfiles/";

  ofstream myfile; myfile.open("L_Scin_Eff_all.3.txt");
  ifstream runfile("runlist3.txt"); const int NUMR = 10000;
 
  TString GeneralCut = "abs(L.tr.x)<0.75 && abs(L.tr.y)<0.55 && abs(L.tr.th)<0.15 && abs(L.tr.ph)<0.045";
  TString TriggerCut3 = "((DBB.evtypebits>>3)&1)";
  TString TriggerCut4 = "((DBB.evtypebits>>4)&1)";
  TString TriggerCut7 = "((DBB.evtypebits>>7)&1)";
  
  TFile *geocut = new TFile(Form("L_Geo_Cut_%d.root",run_number), "UPDATE");
  
  TCutG* phi_y_cut = (TCutG*) gROOT->FindObject("phi_y_cut"); //looking for old cut definition
  TCutG* theta_dp_cut = (TCutG*) gROOT->FindObject("theta_dp_cut"); //looking for old cut definition
  
  if(!phi_y_cut || !theta_dp_cut){
    
    TCanvas *c1= new TCanvas("c1","c1",1000,600);
    c1->Divide(2,1); 
    TPad *bpad = gPad;

    //    runfile >> run_number;
    TString file1 = gDir + Form("e08014_less_%d.root",3590);
    TChain *T1 = new TChain("T");
    T1->Add(file1);

    cerr <<"========= Create Acceptance CutG ..."<<endl;
    c1->cd(1);
    T1->Draw("L.tr.tg_ph:L.tr.tg_y>>gc1(300,-0.05,0.05,300,-0.05,0.05)",GeneralCut,"colz");
    gc1->SetXTitle("L.tr.tg_y"); gc1->SetYTitle("L.tr.tg_ph");  c1->Update();
    phi_y_cut = (TCutG*) (bpad->WaitPrimitive("CUTG", "CutG"));c1->Update();
    phi_y_cut->SetName("phi_y_cut"); phi_y_cut->SetVarX("L.tr.tg_y"); phi_y_cut->SetVarY("L.tr.tg_ph");
    cerr <<"Phi vs Y Done!"<<endl;

    c1->cd(2);
    T1->Draw("L.tr.tg_th:L.tr.tg_dp>>gc2(300,-0.05,0.05,300,-0.15,0.15)",GeneralCut,"colz");
    gc2->SetXTitle("L.tr.tg_dp"); gc2->SetYTitle("L.tr.tg_th");  c1->Update()
    theta_dp_cut = (TCutG*) (bpad->WaitPrimitive("CUTG", "CutG"));c1->Update();
    theta_dp_cut->SetName("theta_dp_cut"); theta_dp_cut->SetVarX("L.tr.tg_dp"); theta_dp_cut->SetVarY("L.tr.tg_th");
    cerr <<"Theta vs Dp Done!"<<endl;

    geocut->cd();
    phi_y_cut->Write("", TObject::kOverwrite); // Overwrite old cut
    theta_dp_cut->Write("", TObject::kOverwrite); // Overwrite old cut

    T1->Delete();
  }
  
  TString electron_cut = "L.cer.asum_c>300.&&(L.prl1.e+L.prl2.e)/L.tr.p/1000.>0.6";
  
   TString cut_t3 = TriggerCut3+"&&" + GeneralCut + "&&phi_y_cut && theta_dp_cut "+ "&&"+electron_cut;
   TString cut_t4 = TriggerCut4+"&&" + GeneralCut + "&&phi_y_cut && theta_dp_cut "+ "&&"+electron_cut;  
   TString cut_t7 = TriggerCut7+"&&" + GeneralCut + "&&phi_y_cut && theta_dp_cut "+ "&&"+electron_cut;
  
  // TString cut_t3 = TriggerCut3+"&&" + GeneralCut +  "&&electron_cut";
  // TString cut_t4 = TriggerCut4+"&&" + GeneralCut +  "&&electron_cut";
  // TString cut_t7 = TriggerCut7+"&&" + GeneralCut +  "&&electron_cut";
  
//   TString cut_t3 = TriggerCut3+"&&" + GeneralCut;
//   TString cut_t4 = TriggerCut4+"&&" + GeneralCut;
//   TString cut_t7 = TriggerCut7+"&&" + GeneralCut;

  float rnum[NUMR]; 
  float Teff[NUMR];

  int jj;  float evt,evt0,evt3,evt4,evt7;
  int ps1=1,ps2=1, ps3=1,ps4=1,ps6=1,ps7=1;
  
  int R,L,comment;

  while(!(runfile.eof())){
    runfile >> run_number >> L >>R >> comment ; rnum[jj] = run_number;

    if(L ==1){
    
      filename = gDir + Form("e08014_L_%d.root",run_number);
      //   TChain *T = new TChain("T");
      //   T->Add(filename);
      TFile  *ff1 = new TFile(filename);
      TTree * T = (TTree *) ff1->GetObjectChecked("T", "TTree");
    
      TArrayI prescales = Run_Data->GetParameters()->GetPrescales();
    
      cout << "Run: " << rnum[jj] << endl;
      //for( int i = 0; i<prescales->GetSize(); ++i )   
      for( int i = 0; i<8; ++i ) 
	cout << "prescale[" << i+1 << "] = " << prescales[i] << ",   ";
      cout<<endl; 
    
      ps1 = prescales[0]; ps2 = prescales[1]; ps3 = prescales[2];
      ps4 = prescales[3]; ps6 = prescales[5]; ps7 = prescales[6];

      evt3 = T->GetEntries(cut_t3);
      evt4 = T->GetEntries(cut_t4);
      evt7 = T->GetEntries(cut_t7);
    
      evt = ps3*evt3 + ps7*evt7;
      evt0 =  ps3*evt3 + ps4*evt4 + ps7*evt7;

      Teff[jj] = 100.0*evt/evt0;
       
      cerr << "-- Scin Eff for Run"<< rnum[jj] <<" is "
	   << Teff[jj] <<"% --"<<endl;

      myfile << setiosflags(ios::left) << setw(5) << setiosflags(ios::fixed) << setprecision(1) << run_number << "   ";
      myfile << setiosflags(ios::left) << setw(7) << setiosflags(ios::fixed) << setprecision(4) << Teff[jj] << endl;;

      jj++; ff1->Close();
    }
  }
 
  const int TTT1 = jj;
  gStyle->SetOptStat(0);
  TCanvas *cc = new TCanvas("cc","Trigger 1",1200,800);
  // cc->Divide(1,2);cc->cd(1);
  TH2F *hh = new TH2F("hh","",TTT1,3500,4320,TTT1,-10.0,110.0);
  hh->GetXaxis()->SetTitle("Run Number");
  hh->GetXaxis()->CenterTitle();
  hh->GetYaxis()->SetTitle("HRS-L Scintillator Efficiency (%)");
  hh->GetYaxis()->CenterTitle();
  hh->Draw();
  
  TGraph *t3 = new TGraph(TTT1,rnum,Teff);
  t3->SetMarkerStyle(20);
  t3->SetMarkerColor(2);
  t3->Draw("P");
 
  TLegend *l1 = new TLegend(0.4,0.4,0.6,0.6);
  l1->AddEntry(t3,"#epsilon_{Scin} = #frac{N_{t3}+N_{t7}}{N_{t3}+N_{t4}+N_{t7}}  ","p");
  l1->SetTextFont(22);
  l1->Draw();
}


void plot(){

  ifstream file1("L_Scin_Eff_a.txt");
  
  const int NUMR = 2000;
  float rnum[NUMR] ; 
  float Teff[NUMR];

  int i=0,j=0;
  while(!(file1.eof())){
    file1 >> rnum[i] >> Teff[i];
    i++;
  }
 
 const int TTT1 = i;
  gStyle->SetOptStat(0);
  TCanvas *cc = new TCanvas("cc","Trigger 1",1200,800);
  // cc->Divide(1,2);cc->cd(1);
  TH2F *hh = new TH2F("hh","",TTT1,3500,4320,TTT1,99.0,100.1);
  hh->GetXaxis()->SetTitle("Run Number");
  hh->GetXaxis()->CenterTitle();
  hh->GetYaxis()->SetTitle("HRS-L Scintillator Efficiency (%)");
  hh->GetYaxis()->CenterTitle();
  hh->Draw();
  
  TGraph *t3 = new TGraph(TTT1,rnum,Teff);
  t3->SetMarkerStyle(20);
  t3->SetMarkerColor(2);
  t3->Draw("P");
 
  // TLegend *l1 = new TLegend(0.4,0.4,0.6,0.6);
  // l1->SetTextFont(22);
  // l1->AddEntry(t3,"#epsilon_{Scin} = #frac{ps3*N_{t3}+ps7*N_{t7}}{ps3*N_{t3}+ps4*N_{t4}+ps7*N_{t7}}  ","p");
  // l1->Draw();


}
