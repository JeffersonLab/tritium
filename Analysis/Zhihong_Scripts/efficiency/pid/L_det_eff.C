////////////////////////////////////////
// Cerenkov Cut Efficiency Study
// For x>2 experiment
// -- Zhihong Ye 07/01/2011
////////////////////////////////////////

#include <iostream>
#include <iomanip>
using namespace std;

void L_det_eff(){
  
  gStyle->SetOptStat(0);
  
  TString gDir = "/w/halla-2/e08014/Rootfiles/";
  //  TString gDir = "/GoFlex/Sample/";
  Int_t  run_number = 3735;    
  // cerr <<"Run number: "; cin >> run_number;

  //  filename = gDir + Form("e08014_HRS_%d.root",run_number);
  filename = gDir + Form("e08014_Ls_%d.root",run_number);
  TChain *T = new TChain("T");
  T->Add(filename);

  TString GeneralCut = "L.tr.n==1 && abs(L.tr.x)<0.75 && abs(L.tr.y)<0.55 && abs(L.tr.th)<0.15 && abs(L.tr.ph)<0.045";
  TString TriggerCut = "((DBB.evtypebits>>7)&1)";
  
  TFile *geocut = new TFile(Form("L_Geo_Cut_%d.root",run_number), "UPDATE");
 
  
  TCutG* phi_y_cut = (TCutG*) gROOT->FindObject("phi_y_cut"); //looking for old cut definition
  TCutG* theta_dp_cut = (TCutG*) gROOT->FindObject("theta_dp_cut"); //looking for old cut definition
  
  if(!phi_y_cut || !theta_dp_cut){
    
    TCanvas *c1= new TCanvas("c1","c1",1000,600);
    c1->Divide(2,1); 
    TPad *dpad = gPad;
 
    cerr <<"========= Create Acceptance CutG ..."<<endl;
    c1->cd(1);
    T->Draw("L.tr.tg_ph:L.tr.tg_y>>gc1(300,-0.05,0.05,300,-0.05,0.05)",GeneralCut,"colz");c1->Update();
    geocut->cd();
    phi_y_cut = (TCutG*) (dpad->WaitPrimitive("CUTG", "CutG"));c1->Update();
    phi_y_cut->SetName("phi_y_cut"); phi_y_cut->SetVarX("L.tr.tg_y"); phi_y_cut->SetVarY("L.tr.tg_ph");
    c1->cd(2);
    T->Draw("L.tr.tg_th:L.tr.tg_dp>>gc2(300,-0.05,0.05,300,-0.15,0.15)",GeneralCut,"colz");c1->Update();
    geocut->cd();
    theta_dp_cut = (TCutG*) (dpad->WaitPrimitive("CUTG", "CutG"));c1->Update();
    theta_dp_cut->SetName("theta_dp_cut"); theta_dp_cut->SetVarX("L.tr.tg_dp"); theta_dp_cut->SetVarY("L.tr.tg_th");
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
    T->Draw("L.prl2.e:L.prl1.e>>gc3(300,-100.,3000.0,300,-100.0,3000.0)",GeneralCut,"colz");c2->Update();
    geocut->cd();
    electron_cut = (TCutG*) (cpad->WaitPrimitive("CUTG", "CutG"));c2->Update();
    electron_cut->SetName("electron_cut"); electron_cut->SetVarX("L.prl1.e"); electron_cut->SetVarY("L.prl2.e");

    c2->cd(2);
    T->Draw("L.prl2.e:L.prl1.e>>gc4(300,-100.,500.0,300,-100.0,500.0)",GeneralCut,"colz");c2->Update();
    geocut->cd();
    pion_cut = (TCutG*) (cpad->WaitPrimitive("CUTG", "CutG"));c2->Update();
    pion_cut->SetName("pion_cut"); pion_cut->SetVarX("L.prl1.e"); pion_cut->SetVarY("L.prl2.e");

    geocut->cd();
    electron_cut->Write("", TObject::kOverwrite); // Overwrite old cut
    pion_cut->Write("", TObject::kOverwrite); // Overwrite old cut
  }


  
  /////////////////////////////////////
  // Cerenkov Detection Eff
  /////////////////////////////////////
  TString cute1 = TriggerCut+"&&" + GeneralCut +"&&phi_y_cut&& theta_dp_cut ";

  TCanvas *c3= new TCanvas("c3","c3",1000,600);
  c3->Divide(2,1); 
  
  c3->cd(1);
  T->Draw("L.prl1.e:L.prl2.e>>f1(600,-100.,3000.,600,-100.,3000.)",TCut(cute1));
  f1->SetTitle("PRL1 vs PRL2"); f1->SetXTitle("L.prl2.e"); f1->SetYTitle("L.prl1.e");  c3->Update();
  T->SetMarkerColor(4);
  TString cute2 = cute1 + "&&electron_cut";
  T->Draw("L.prl1.e:L.prl2.e>>f2(600,-100.,3000.,600,-100.,3000.)",TCut(cute2));  c3->Update();
  float evt1 = f2->GetEntries();
  f1->Draw();f2->Draw("same");


  T->SetMarkerColor(1);
  
  c3->cd(2);
  TString cute3 = "L.cer.asum_c>50.0 &&" + cute2;
  T->Draw("L.cer.asum_c>>f3(300,-100.,5000.0)",TCut(cute3));gPad->SetLogy(1);
  f3->SetTitle("Cerenkov Sum After Calorimeter Cut"); f3->SetXTitle("L.cer.asum_c"); c3->Update();

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
  // Pion Rejectors Detection Eff
  /////////////////////////////////////
 
  TCanvas *c4= new TCanvas("c4","c4",1000,600);
  c4->Divide(2,1); 
  
  c4->cd(1);

  T->Draw("L.cer.asum_c>>f4(600,-100.,5000.)",TCut(cute1)); gPad->SetLogy(1);
  f4->SetTitle("Cerenkov Sum"); f4->SetXTitle("L.cer.asum_c"); c3->Update();
  T->SetFillColor(4);
  TString cute4 = cute1 + "&&L.cer.asum_c>400&&L.cer.asum_c<600.";
  //TString cute4 = "((DBB.evtypebits>>6)&1)&&L.cer.asum_c>600&&L.cer.asum_c<900.";
 
  T->Draw("L.cer.asum_c>>f5(600,10.,5000.)",TCut(cute4));gPad->SetLogy(1);
  f5->SetTitle("Cerenkov Sum"); f5->SetXTitle("L.cer.asum_c"); c3->Update();
  float evt3 = f5->GetEntries();
  f4->Draw();f5->Draw("same");
  T->SetFillColor(1);
  
  
  c4->cd(2);gPad->SetLogy(0);
  TString cute5 = cute4 + "&&(L.prl2.e+L.prl1.e)>0.";
  //  TString cute5 = cute4 + "&&(L.prl2.asum_c+L.prl1.asum_c)>0.";
    //TString cute5 = "((DBB.evtypebits>>6)&1)&&L.cer.asum_c>600&&L.cer.asum_c<900.&&(L.prl2.e+L.prl1.e)>0.";
  T->Draw("L.prl1.e:L.prl2.e>>f6(600,-100.,3000.,600,-100.,3000.)",TCut(cute5),"colz"); 
  f6->SetTitle("PRL1 vs PRL2 after Cerenkove Cut"); f6->SetXTitle("L.prl2.e"); f6->SetYTitle("L.prl1.e");  c3->Update();

  float evt4 = f6->GetEntries();

  cerr << "Pion Rejectors Detection Efficiency is: "<< evt4/evt3*100.<<"%"<<endl;

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
