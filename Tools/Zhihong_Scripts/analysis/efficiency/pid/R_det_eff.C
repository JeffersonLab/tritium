////////////////////////////////////////
// Cerenkov Cut Efficiency Study
// For x>2 experiment
// -- Zhihong Ye 07/01/2011
////////////////////////////////////////

#include <iostream>
#include <iomanip>
using namespace std;

void R_det_eff(){
  
  gStyle->SetOptStat(0);
 
  //  TString gDir = "/GoFlex/Rootfiles/";
  TString gDir = "/work/halla/e08014/disk1/Rootfiles/";
  Int_t  run_number;    
  cerr <<"Run number: "; cin >> run_number;

  filename = gDir + Form("e08014_Rs_%d.root",run_number);
  TChain *T = new TChain("T");
  T->Add(filename);

  TString GeneralCut = "R.tr.n==1&& abs(R.tr.x)<0.75 && abs(R.tr.y)<0.55 && abs(R.tr.th)<0.15 && abs(R.tr.ph)<0.045";
  GeneralCut += "&& abs(R.tr.tg_ph)<0.04 && abs(R.tr.tg_y)<0.06 && abs(R.tr.tg_th)<0.06 && abs(R.tr.tg_dp)<0.05";
  TString TriggerCut = "((DBB.evtypebits>>6)&1)";
    
  TFile *geocut = new TFile(Form("Geo_Cut_%d.root",run_number), "UPDATE");
  //TCutG* phi_y_cut, *theta_dp_cut;
  
  TCutG* phi_y_cut = (TCutG*) gROOT->FindObject("phi_y_cut"); //looking for old cut definition
  TCutG* theta_dp_cut = (TCutG*) gROOT->FindObject("theta_dp_cut"); //looking for old cut definition
  
  if(!phi_y_cut || !theta_dp_cut){
    
    TCanvas *c1= new TCanvas("c1","c1",1000,600);
    c1->Divide(2,1); 
    TPad *dpad = gPad;
 
    cerr <<"========= Create Acceptance CutG ..."<<endl;
    c1->cd(1);
    T->Draw("R.tr.tg_ph:R.tr.tg_y>>gc1(300,-0.05,0.05,300,-0.05,0.05)",GeneralCut,"colz");c1->Update();
    geocut->cd();
    phi_y_cut = (TCutG*) (dpad->WaitPrimitive("CUTG", "CutG"));c1->Update();
    phi_y_cut->SetName("phi_y_cut"); phi_y_cut->SetVarX("R.tr.tg_y"); phi_y_cut->SetVarY("R.tr.tg_ph");
    c1->cd(2);
    T->Draw("R.tr.tg_th:R.tr.tg_dp>>gc2(300,-0.05,0.05,300,-0.15,0.15)",GeneralCut,"colz");c1->Update();
    geocut->cd();
    theta_dp_cut = (TCutG*) (dpad->WaitPrimitive("CUTG", "CutG"));c1->Update();
    theta_dp_cut->SetName("theta_dp_cut"); theta_dp_cut->SetVarX("R.tr.tg_dp"); theta_dp_cut->SetVarY("R.tr.tg_th");
    geocut->cd();
    
    phi_y_cut->Write("", TObject::kOverwrite); // Overwrite old cut
    theta_dp_cut->Write("", TObject::kOverwrite); // Overwrite old cut
  }
  
  TCutG* electron_cut = (TCutG*) gROOT->FindObject("electron_cut"); //looking for old cut definition
  TCutG* pion_cut = (TCutG*) gROOT->FindObject("pion_cut"); //looking for old cut definition
  
  
  if(!electron_cut || !pion_cut){
    TCanvas *c2= new TCanvas("c2","c2",1000,600);
    c2->Divide(2,1); 
    TPad *cpad = gPad;

    cerr <<"========= Create Electron and Pion CutG On Calorimeter ..."<<endl;
    c2->cd(1);
    T->Draw("R.sh.e/R.tr.p/1000.:R.ps.e/R.tr.p/1000.>>gc3(100,-0.1,1.0,100,-0.1,1.0)",GeneralCut,"colz");c2->Update();
    geocut->cd();
    electron_cut = (TCutG*) (cpad->WaitPrimitive("CUTG", "CutG"));c2->Update();
    electron_cut->SetName("electron_cut"); electron_cut->SetVarX("R.ps.e/R.tr.p/1000."); electron_cut->SetVarY("R.sh.e/R.tr.p/1000.");

    c2->cd(2);
    T->Draw("R.sh.e/R.tr.p/1000.:R.ps.e/R.tr.p/1000.>>gc4(100,-0.1,0.5,100,-0.1,0.5)",GeneralCut,"colz");c2->Update();
    geocut->cd();
    pion_cut = (TCutG*) (cpad->WaitPrimitive("CUTG", "CutG"));c2->Update();
    pion_cut->SetName("pion_cut"); pion_cut->SetVarX("R.ps.e/R.tr.p/1000."); pion_cut->SetVarY("R.sh.e/R.tr.p/1000.");

    geocut->cd();
    electron_cut->Write("", TObject::kOverwrite); // Overwrite old cut
    pion_cut->Write("", TObject::kOverwrite); // Overwrite old cut
  }
  
  TString cute1 = TriggerCut+"&&" + GeneralCut + "&&R.tr.p>2.0&&R.tr.p<6.0";
  TCanvas *c3= new TCanvas("c3","c3",1000,600);
  c3->Divide(2,1); 
  
  /////////////////////////////////////
  // Cerenkov Detection Eff
  /////////////////////////////////////

  c3->cd(1);
  T->Draw("R.ps.e/R.tr.p/1000.:R.sh.e/R.tr.p/1000.>>f1(100,-0.1,1.0,100,-0.1,1.0)",TCut(cute1));
  f1->SetTitle("PS vs SH"); f1->SetXTitle("R.sh.e/R.tr.p/1000."); f1->SetYTitle("R.ps.e/R.tr.p/1000.");  c3->Update();
  T->SetMarkerColor(4);
  TString cute2 = cute1 + "&&electron_cut";
  T->Draw("R.ps.e/R.tr.p/1000.:R.sh.e/R.tr.p/1000.>>f2(100,-0.1,1.0,100,-0.1,1.0)",TCut(cute2));  c3->Update();
  float evt1 = f2->GetEntries();
  f1->Draw();f2->Draw("same");

  T->SetMarkerColor(1); c3->cd(2);
  TString cute3 = "R.cer.asum_c>70.0 &&" + cute2;
  T->Draw("R.cer.asum_c>>f3(300,-100.,5000.0)",TCut(cute3));gPad->SetLogy(1);
  f3->SetTitle("Cerenkov Sum After Calorimeter Cut"); f3->SetXTitle("R.cer.asum_c"); c3->Update();
  float evt2 = f3->GetEntries();

  cerr << "Cerenkov Detection Efficiency is: "<< evt2/evt1*100.<<"%"<<endl;
  TString num1 = Form("Detection Efficiency: %3.2lf",100*evt2/evt1);
  num1+="%";
  TLatex  *tex1 = new TLatex(0.3,0.85,num1);
  tex1->SetNDC();
  tex1->SetTextSize(0.04);
  tex1->SetTextColor(3);
  tex1->Draw();
  c3->Update(); 

  /*
  /////////////////////////////////////
  // Shower&Preshower Detection Eff
  /////////////////////////////////////
  TCanvas *c4= new TCanvas("c4","c4",1000,600);
  c4->Divide(2,1); 
  c4->cd(1);

  T->Draw("R.cer.asum_c>>f4(600,-100.,5000.)",TCut(cute1)); gPad->SetLogy(1);
  f4->SetTitle("Cerenkov Sum"); f4->SetXTitle("R.cer.asum_c"); c3->Update();
  T->SetFillColor(4);
  
  TString cute4 = cute1 + "&&R.cer.asum_c>500&&R.cer.asum_c<700.";
  T->Draw("R.cer.asum_c>>f5(600,10.,5000.)",TCut(cute4));gPad->SetLogy(1);
  f5->SetTitle("Cerenkov Sum"); f5->SetXTitle("R.cer.asum_c"); c3->Update();
  float evt3 = f5->GetEntries();
  f4->Draw();f5->Draw("same");
  T->SetFillColor(1);
  
  c4->cd(2);gPad->SetLogy(0);
  TString cute5 = cute4 + "&&(R.sh.e/R.tr.p/1000.+R.ps.e/R.tr.p/1000.)>0.";
  T->Draw("R.ps.e/R.tr.p/1000.:R.sh.e/R.tr.p/1000.>>f6(300,-0.1,1.0,300,-0.1,1.0)",TCut(cute5),"colz"); 
  f6->SetTitle("PS vs SH after Cerenkove Cut"); f6->SetXTitle("R.sh.e/R.tr.p/1000."); f6->SetYTitle("R.ps.e/R.tr.p/1000.");  c3->Update();
  float evt4 = f6->GetEntries();

  cerr << "Shower&Preshower Detection Efficiency is: "<< evt4/evt3*100.<<"%"<<endl;
  TString num2 = Form("Detection Efficiency: %3.2lf",100*evt4/evt3);
  num2+="%";
  TLatex  *tex2 = new TLatex(0.3,0.85,num2);
  tex2->SetNDC();
  tex2->SetTextSize(0.04);
  tex2->SetTextColor(3);
  tex2->Draw();
  c3->Update(); 
  */
}
