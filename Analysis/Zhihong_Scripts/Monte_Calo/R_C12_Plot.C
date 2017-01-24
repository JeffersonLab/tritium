// #include <iostream>
// #include <iomanip>
// using namespace std;
//void plot()
{

  TString hrs="R";TString target="c12";  TString kin="kin3.1.x1mm"; int run_number = 3668;

  TString GeneralCut = "R.tr.n==1";//&&abs(R.tr.x)<0.75 && abs(R.tr.y)<0.55 && abs(R.tr.th)<0.15 && abs(R.tr.ph)<0.045";
  TString TriggerCut3 = "((DBB.evtypebits>>1)&1)";
  TString TriggerCut7 = "((DBB.evtypebits>>6)&1)";
  TString Acc_real = "abs(ExTgtR.th)<0.04&&abs(ExTgtR.ph)<0.025&&abs(ExTgtR.y)<0.15&&abs(ExTgtR.dp)<0.03&&abs(RctPtR.z)<0.05";
  TString Acc_simu = "abs(th)<0.04&&abs(R.gold.ph)<0.03&&abs(R.gold.y)<0.003";

  TString Cut_All = GeneralCut+"&&"+TriggerCut3+"&&"+Acc_real
    +"&&R.cer.asum_c>=50.&&R.sh.e>=200.&&epR>=0.5";
 
  double scale;  int norm = 100000;
  //  TString filename="e08014_"+hrs+"_"+Form("%d",run_number)+"_"+target+"_"+kin+".root";
  TString filename=Form("/w/halla-2/e08014/Rootfiles/e08014_%d.root",run_number);
  TFile *f1 = new TFile(filename);
  
  TString filename2=target+"_"+hrs+"_"+kin+".root";
  
  TFile *f2 = new TFile(filename2);

  
  // .x cutL.C;
  TCanvas *cc = new TCanvas("cc","cc",1000,1000);
  cc->Divide(2,2);
  
  f1->cd();
  cc->cd(1);
  T->Draw("R.tr.r_x>>h1(200,-1.0,1.0)",Cut_All);
  h1->SetLineColor(2); scale = norm/h1->Integral(); h1->Scale(scale); 
  h1->SetTitle("X_fp (weighted), Red->real data, Blue->Simulated data");
  h1->Draw();
  
  cc->cd(2);
  T->Draw("R.tr.r_y>>h2(200,-0.1,0.1)",Cut_All);
  h2->SetLineColor(2); scale = norm/h2->Integral(); h2->Scale(scale); 
  h2->SetTitle("Y_fp (weighted), Red->real data, Blue->Simulated data");
  h2->Draw();
  
  cc->cd(3);
  T->Draw("R.tr.r_th>>h3(200,-0.1,0.1)",Cut_All);
  h3->SetLineColor(2); scale = norm/h3->Integral(); h3->Scale(scale); 
  h3->SetTitle("Th_fp (weighted), Red->real data, Blue->Simulated data");
  h3->Draw();
  
  cc->cd(4);
  T->Draw("R.tr.r_ph>>h4(200,-0.1,0.1)",Cut_All);
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
  
  
  TCanvas *cf = new TCanvas("cf","cf",1000,1000);
  cf->Divide(2,2);

  f1->cd();
  cf->cd(1);
  T->Draw("R.gold.dp>>j1(200,-0.1,0.1)",Cut_All);
  j1->SetLineColor(2); scale = norm/j1->Integral(); j1->Scale(scale);
  j1->SetTitle("Dp_tg (weighted), Red->real data, Blue->Simulated data");
  j1->Draw();
  
  cf->cd(2);
  T->Draw("R.gold.y>>j2(200,-0.01,0.01)",Cut_All);
  j2->SetLineColor(2); scale = norm/j2->Integral(); j2->Scale(scale); 
  j2->SetTitle("Y_tg (weighted), Red->real data, Blue->Simulated data");
  j2->Draw();

  // cf->cd(2);
  // T->Draw("RctPtR.z>>j2(200,-0.01,0.01)",Cut_All);
  // j2->SetLineColor(2); scale = norm/j2->Integral(); j2->Scale(scale); 
  // j2->SetTitle("RctPnt_Z (weighted), Red->real data, Blue->Simulated data");
  // j2->Draw();

  cf->cd(3);
  T->Draw("R.gold.th>>j3(200,-0.1,0.1)",Cut_All);
  j3->SetLineColor(2); scale = norm/j3->Integral(); j3->Scale(scale); 
  j3->SetTitle("Th_tg (weighted), Red->real data, Blue->Simulated data");
  j3->Draw();
  
  cf->cd(4);
  T->Draw("R.gold.ph>>j4(200,-0.1,0.1)",Cut_All);
  j4->SetLineColor(2); scale = norm/j4->Integral(); j4->Scale(scale); 
  j4->SetTitle("Ph_tg (weighted), Red->real data, Blue->Simulated data");
  j4->Draw();
  
  f2->cd();
  hdptg->SetLineColor(4);
  scale = norm/hdptg->Integral();  hdptg->Scale(scale);  cf->cd(1); hdptg ->Draw("");
  hytg->SetLineColor(4);
  scale = norm/hytg->Integral();  hytg->Scale(scale);  cf->cd(2); hytg ->Draw("");
 
  // hrpz->SetLineColor(4);
  // scale = norm/hrpz->Integral();  hrpz->Scale(scale);  cf->cd(2); hrpz ->Draw("");

  hthtg->SetLineColor(4);
  scale = norm/hthtg->Integral(); hthtg->Scale(scale); cf->cd(3); hthtg->Draw("");
  hphtg->SetLineColor(4);
  scale = norm/hphtg->Integral(); hphtg->Scale(scale); cf->cd(4); hphtg->Draw("");
 
  f1->cd();
  cf->cd(1); j1->Draw("same");  cf->cd(2); j2->Draw("same");
  cf->cd(3); j3->Draw("same");  cf->cd(4); j4->Draw("same");
  
  TString pic=target+"_"+hrs+"_"+kin+Form("_%d_tg.png",run_number);
  cf->Print(pic);
 
}
