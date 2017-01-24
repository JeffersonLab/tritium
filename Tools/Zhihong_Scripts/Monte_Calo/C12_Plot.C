// #include <iostream>
// #include <iomanip>
// using namespace std;
//void plot()
{

	gStyle->SetOptStat(0);

  TString hrs="L";TString target="c12";  TString kin="kin3.1_clean"; int run_number = 3670;

  TString GeneralCut = "L.tr.n==1";//&&abs(L.tr.x)<0.75 && abs(L.tr.y)<0.55 && abs(L.tr.th)<0.15 && abs(L.tr.ph)<0.045";
  TString TriggerCut3 = "((DBB.evtypebits>>3)&1)";
  TString TriggerCut7 = "((DBB.evtypebits>>7)&1)";
  TString Acc_real = "abs(ExTgtL.th)<0.04&&abs(ExTgtL.ph)<0.025&&abs(ExTgtL.y)<0.02&&abs(ExTgtL.dp)<0.05&&abs(RctPtL.z)<0.07";

  TString Cut_All = GeneralCut+"&&"+TriggerCut3+"&&"+Acc_real
    +"&&L.cer.asum_c>=50.&&L.prl2.e>=100.&&(L.prl1.e+L.prl2.e)/L.tr.p/1000.>=0.5";
  //TString Cut_All ="1";

  double scale;  int norm = 100000;
  TString filename=Form("/work/halla/e08014/disk1/Rootfiles/e08014_%d.root",run_number);
  TFile *f1 = new TFile(filename);
  
  TString filename2=target+"_"+hrs+"_"+kin+".root";
  TFile *f2 = new TFile(filename2);
  // TChain* T = new TChain("SAMC");
  // SAMC->Add(filename2);
  // TString filename2=target+"_"+hrs+"_"+kin+"_1.root";
  // SAMC->Add(filename2);
  // TString filename2=target+"_"+hrs+"_"+kin+"_2.root";
  // SAMC->Add(filename2);

  // .x cutL.C;
  TCanvas *cc = new TCanvas("cc","cc",1200,800);
  cc->Divide(2,2);
  
  f1->cd();
  cc->cd(1);
  T->Draw("L.tr.r_x>>h1(200,-0.5,0.5)",Cut_All);
  h1->SetLineColor(2); scale = norm/h1->Integral(); h1->Scale(scale); 
  h1->SetTitle("X_fp (weighted), Red->real data, Blue->Simulated data");
  h1->Draw();
  
  cc->cd(2);
  T->Draw("L.tr.r_y>>h2(200,-0.05,0.05)",Cut_All);
  h2->SetLineColor(2); scale = norm/h2->Integral(); h2->Scale(scale); 
  h2->SetTitle("Y_fp (weighted), Red->real data, Blue->Simulated data");
  h2->Draw();
  
  cc->cd(3);
  T->Draw("L.tr.r_th>>h3(200,-0.02,0.02)",Cut_All);
  h3->SetLineColor(2); scale = norm/h3->Integral(); h3->Scale(scale); 
  h3->SetTitle("Th_fp (weighted), Red->real data, Blue->Simulated data");
  h3->Draw();
  
  cc->cd(4);
  T->Draw("L.tr.r_ph>>h4(200,-0.05,0.05)",Cut_All);
  h4->SetLineColor(2); scale = norm/h4->Integral(); h4->Scale(scale);
  h4->SetTitle("Ph_fp (weighted), Red->real data, Blue->Simulated data");
  h4->Draw();
  
  f2->cd();
  hxfp->SetLineColor(4);
  scale = norm/hxfp->Integral();  hxfp->Scale(scale);  cc->cd(1); hxfp ->Draw();
  hyfp->SetLineColor(4);
  scale = norm/hyfp->Integral();  hyfp->Scale(scale);  cc->cd(2); hyfp ->Draw();
  hthfp->SetLineColor(4);
  scale = norm/hthfp->Integral(); hthfp->Scale(scale); cc->cd(3); hthfp->Draw();
  hphfp->SetLineColor(4);
  scale = norm/hphfp->Integral(); hphfp->Scale(scale); cc->cd(4); hphfp->Draw();

  f1->cd();
  cc->cd(1); h1->Draw("same");  cc->cd(2); h2->Draw("same");
  cc->cd(3); h3->Draw("same");  cc->cd(4); h4->Draw("same");
  
  TString pic=target+"_"+hrs+"_"+kin+Form("_%d_fp.png",run_number);
  cc->Print(pic);
  
  TCanvas *cf = new TCanvas("cf","cf",1200,800);
  cf->Divide(2,2);

  f1->cd();
  cf->cd(1);
  T->Draw("ExTgtL.dp>>j1(200,-0.06,0.06)",Cut_All);
  j1->SetLineColor(2); scale = norm/j1->Integral(); j1->Scale(scale);
  j1->SetTitle("Dp_tg (weighted), Red->real data, Blue->Simulated data");
  j1->Draw();
  
   cf->cd(2);
   T->Draw("ExTgtL.y>>j2(200,-0.07,0.07)",Cut_All);
   j2->SetLineColor(2); scale = norm/j2->Integral(); j2->Scale(scale); 
   j2->SetXTitle("y_{tg}");
   j2->GetXaxis()->CenterTitle(1);
   j2->GetXaxis()->SetTitleFont(32);
   j2->GetXaxis()->SetTitleSize(0.06);
   j2->SetTitle("Y_tg (weighted), Red->real data, Blue->Simulated data");
   j2->Draw();
/*
  cf->cd(2);
  T->Draw("RctPtL.z>>j2(200,-0.12,0.12)",Cut_All);
  j2->SetLineColor(2); scale = norm/j2->Integral(); j2->Scale(scale); 
  j2->SetXTitle("RctPtL.z");
  j2->GetXaxis()->CenterTitle(1);
  j2->GetXaxis()->SetTitleFont(32);
  j2->GetXaxis()->SetTitleSize(0.06);
  j2->SetTitle("RctPnt_Z (weighted), Red->real data, Blue->Simulated data");
  j2->Draw();
*/
  cf->cd(3);
  T->Draw("ExTgtL.th>>j3(200,-0.06,0.06)",Cut_All);
  j3->SetLineColor(2); scale = norm/j3->Integral(); j3->Scale(scale); 
  j3->SetTitle("Th_tg (weighted), Red->real data, Blue->Simulated data");
  j3->Draw();
  
  cf->cd(4);
  T->Draw("ExTgtL.ph>>j4(200,-0.045,0.045)",Cut_All);
  j4->SetLineColor(2); scale = norm/j4->Integral(); j4->Scale(scale); 
  j4->SetTitle("Ph_tg (weighted), Red->real data, Blue->Simulated data");
  j4->Draw();
  
  f2->cd();
  hdptg->SetLineColor(4);
  scale = norm/hdptg->Integral();  hdptg->Scale(scale);  cf->cd(1); hdptg ->Draw("");
 
  hytg->SetLineColor(4);
  scale = norm/hytg->Integral();  hytg->Scale(scale);  cf->cd(2); hytg ->Draw("");
  
//  hrpz->SetLineColor(4);
//  scale = norm/hrpz->Integral();  hrpz->Scale(scale);  cf->cd(2); hrpz ->Draw("");

  hthtg->SetLineColor(4);
  scale = norm/hthtg->Integral(); hthtg->Scale(scale); cf->cd(3); hthtg->Draw("");
  hphtg->SetLineColor(4);
  scale = norm/hphtg->Integral(); hphtg->Scale(scale); cf->cd(4); hphtg->Draw("");
 
  f1->cd();
  cf->cd(1); j1->Draw("same"); 
  cf->cd(2);
  j2->SetLineColor(2); scale = norm/j2->Integral(); j2->Scale(scale); 
  j2->Draw("");
  cf->cd(3); j3->Draw("same");  cf->cd(4); j4->Draw("same");
  
  f2->cd();
  cf->cd(2); hytg->Draw("same"); //hrpz ->Draw("same");

  TString pic=target+"_"+hrs+"_"+kin+Form("_%d_tg.png",run_number);
  cf->Print(pic);
}
