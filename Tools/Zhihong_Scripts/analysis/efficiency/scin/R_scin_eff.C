////////////////////////////////////////
// Scintillator Efficiency Study
// For x>2 experiment
// -- Zhihong Ye 07/08/2011
////////////////////////////////////////

#include <iostream>
#include <iomanip>
using namespace std;

void R_scin_eff(){

  Int_t  run_number;    TString filename;
  TString gDir = "/w/halla-2/e08014/Rootfiles_0721/";

  ofstream myfile; myfile.open("R_Scin_Eff_all.4.txt");
  ifstream runfile("runlist3.txt"); const int NUMR = 10000;
 
  TString GeneralCut = "abs(R.tr.x)<0.75 && abs(R.tr.y)<0.55 && abs(R.tr.th)<0.15 && abs(R.tr.ph)<0.045";
  TString TriggerCut1 = "((DBB.evtypebits>>1)&1)";
  TString TriggerCut2 = "((DBB.evtypebits>>2)&1)";
  TString TriggerCut6 = "((DBB.evtypebits>>6)&1)";
  
  TFile *geocut = new TFile(Form("R_Geo_Cut_%d.root",run_number), "UPDATE");
  
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
    T1->Draw("R.tr.tg_ph:R.tr.tg_y>>gc1(300,-0.05,0.05,300,-0.05,0.05)",GeneralCut,"colz");
    gc1->SetXTitle("R.tr.tg_y"); gc1->SetYTitle("R.tr.tg_ph");  c1->Update();
    phi_y_cut = (TCutG*) (bpad->WaitPrimitive("CUTG", "CutG"));c1->Update();
    phi_y_cut->SetName("phi_y_cut"); phi_y_cut->SetVarX("R.tr.tg_y"); phi_y_cut->SetVarY("R.tr.tg_ph");
    cerr <<"Phi vs Y Done!"<<endl;

    c1->cd(2);
    T1->Draw("R.tr.tg_th:R.tr.tg_dp>>gc2(300,-0.05,0.05,300,-0.15,0.15)",GeneralCut,"colz");
    gc2->SetXTitle("R.tr.tg_dp"); gc2->SetYTitle("R.tr.tg_th");  c1->Update()
    theta_dp_cut = (TCutG*) (bpad->WaitPrimitive("CUTG", "CutG"));c1->Update();
    theta_dp_cut->SetName("theta_dp_cut"); theta_dp_cut->SetVarX("R.tr.tg_dp"); theta_dp_cut->SetVarY("R.tr.tg_th");
    cerr <<"Theta vs Dp Done!"<<endl;

    geocut->cd();
    phi_y_cut->Write("", TObject::kOverwrite); // Overwrite old cut
    theta_dp_cut->Write("", TObject::kOverwrite); // Overwrite old cut

    T1->Delete();
  }
  
//   TCutG* electron_cut = (TCutG*) gROOT->FindObject("electron_cut"); //looking for old cut definition
   
//   if(!electron_cut ){
    
//     //    runfile >> run_number;
//     TString file2 = gDir + Form("e08014_less_%d.root",3590);
//     TChain *T2 = new TChain("T");
//     T2->Add(file2);

//     TCanvas *c2= new TCanvas("c2","c2",800,800);
//     //  c2->Divide(2,1); 
//     TPad *cpad = gPad;

//     cerr <<"========= Create Electron and Pion CutG On Calorimeter ..."<<endl;
//     c2->cd(1);
//     T2->Draw("R.sh.asum_c:R.ps.asum_c>>gc3(300,-100.,6000.0,300,-100.0,6000.0)",GeneralCut,"colz");c2->Update();
//     electron_cut = (TCutG*) (cpad->WaitPrimitive("CUTG", "CutG"));c2->Update();
//     electron_cut->SetName("electron_cut"); electron_cut->SetVarX("R.ps.asum_c"); electron_cut->SetVarY("R.sh.asum_c");
//     cerr <<"Electron Done!"<<endl;

//      geocut->cd();
//      electron_cut->Write("", TObject::kOverwrite); // Overwrite old cut

//      T2->Delete();
//   }
  TString electron_cut = "R.cer.asum_c>300.&&(R.sh.e+R.ps.e)/R.tr.p/1000.>0.6";
  
   TString cut_t1 = TriggerCut1+"&&" + GeneralCut + "&&phi_y_cut && theta_dp_cut "+ "&&"+electron_cut;
   TString cut_t2 = TriggerCut2+"&&" + GeneralCut + "&&phi_y_cut && theta_dp_cut "+ "&&"+electron_cut;  
   TString cut_t6 = TriggerCut6+"&&" + GeneralCut + "&&phi_y_cut && theta_dp_cut "+ "&&"+electron_cut;
  
  // TString cut_t1 = TriggerCut1+"&&" + GeneralCut +  "&&electron_cut";
  // TString cut_t2 = TriggerCut2+"&&" + GeneralCut +  "&&electron_cut";
  // TString cut_t6 = TriggerCut6+"&&" + GeneralCut +  "&&electron_cut";
  
//   TString cut_t1 = TriggerCut1+"&&" + GeneralCut;
//   TString cut_t2 = TriggerCut2+"&&" + GeneralCut;
//   TString cut_t6 = TriggerCut6+"&&" + GeneralCut;

  float rnum[NUMR]; 
  float Teff[NUMR];

  int jj;  float evt,evt0,evt1,evt2,evt6;
  int ps1=1,ps2=1, ps3=1,ps4=1,ps6=1,ps7=1;
  
  int R,L,comment;

  while(!(runfile.eof())){
    runfile >> run_number >> L >> R >> comment; rnum[jj] = run_number;
    
    if( R == 1){

      filename = gDir + Form("e08014_R_%d.root",run_number);
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
      
      evt1 = T->GetEntries(cut_t1);
      evt2 = T->GetEntries(cut_t2);
      evt6 = T->GetEntries(cut_t6);
      
      evt = ps1*evt1 + ps6*evt6;
      evt0 =  ps1*evt1 + ps2*evt2 + ps6*evt6;
      
      Teff[jj] = 100.0*evt/evt0;
      
      cerr << "Run"<<run_number<<"-- Scin Eff for Run"<< rnum[jj] <<" is "
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
  TH2F *hh = new TH2F("hh","",TTT1,3500,4300,TTT1,-10.0,110.0);
  hh->GetXaxis()->SetTitle("Run Number");
  hh->GetXaxis()->CenterTitle();
  hh->GetYaxis()->SetTitle("HRS-R Scintillator Efficiency (%)");
  hh->GetYaxis()->CenterTitle();
  hh->Draw();
  
  TGraph *t1 = new TGraph(TTT1,rnum,Teff);
  t1->SetMarkerStyle(20);
  t1->SetMarkerColor(2);
  t1->Draw("P");
 
  TLegend *l1 = new TLegend(0.4,0.4,0.6,0.6);
  l1->AddEntry(t1,"#epsilon_{scin} = frac{ps1*T1+ps6*T6}{ps1*T1+ps2*T2+ps6*T6}","p");
  l1->Draw();
}


void plot(){

  ifstream file1("R_Scin_Eff_a.txt");
  
  const int NUMR = 2000;
  float rnum[NUMR] ; 
  float T1eff[NUMR], T6eff[NUMR];

  int i=0,j=0;
  while(!(file1.eof())){
    file1 >> rnum[i] >> T1eff[i] >> T6eff[i];
    i++;
  }
 
  const int TTT1 = i;
  gStyle->SetOptStat(0);
  TCanvas *cc = new TCanvas("cc","Trigger 1",1200,800);
  // cc->Divide(1,2);cc->cd(1);
  TH2F *hh = new TH2F("hh","",TTT1,3500,4300,TTT1,99.5,100.1);
  hh->GetXaxis()->SetTitle("Run Number");
  hh->GetXaxis()->CenterTitle();
  hh->GetYaxis()->SetTitle("HRS-R Scintillator Efficiency (%)");
  hh->GetYaxis()->CenterTitle();
  hh->Draw();
  
  TGraph *t1 = new TGraph(TTT1,rnum,T1eff);
  t1->SetMarkerStyle(20);
  t1->SetMarkerColor(2);
  t1->Draw("P");
  
  TGraph *t6 = new TGraph(TTT1,rnum,T6eff);
  t6->SetMarkerStyle(21);
  t6->SetMarkerColor(4);
  t6->Draw("P");
  
  TLegend *l1 = new TLegend(0.4,0.4,0.6,0.6);
  //l1->AddEntry(t1,"ps1*T1/(ps1*T1+ps2*T2)","p");
  //l1->AddEntry(t6,"ps6*T6/(ps6*T6+ps2*T2)","p");
  l1->AddEntry(t1,"#epsilon_{Scin} = #frac{ps1*N_{t1}+ps6*N_{t6}}{ps1*N_{t1}+ps2*N_{t2}+ps6*N_{t6}}  ","p");
  l1->Draw();

}
