
  /*C/C++ Includes{{{*/
  #include <stdio.h>
  #include <string>
  #include <iostream>
  #include <fstream>
  #include <vector>
  #include <algorithm>
  #include <map>
  #include <cmath>
  /*}}}*/
  /*ROOT Includes{{{*/
  #include <TSystem.h>
  #include <TString.h>
  #include <TStyle.h>
  #include <Riostream.h>
  #include "TObjString.h"
  #include <TNamed.h>
  #include <TPRegexp.h>
  #include <TObjArray.h>
  #include <TChain.h>
  #include <TMath.h>
  #include <TH1.h>
  #include <TH1F.h>
  #include <TH2F.h>
  #include <TFile.h>
  #include <TROOT.h>
  #include <TF1.h>
  #include <TGraph.h>
  #include <TGraphErrors.h>
  #include <TCanvas.h>
  #include <TDatime.h>
  #include <TError.h>
  #include <TVirtualFitter.h>
  #include <TSQLServer.h>
  #include <TSQLResult.h>
  #include <TSQLRow.h>
  #include <TCut.h>
  #include <TMultiGraph.h>
  #include <TCutG.h>
  #include <TLorentzVector.h>
  #include <TMath.h>
  #include <TLatex.h>
  #include <TRandom3.h>
  //#include <TMatrix.h>
  /*}}}*/
  using namespace std;

//void Check_Point()
int main()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  TString cuteL = "L.tr.n==1 && DBB.evtypebits>>3&1 && epL>0.5 && L.prl2.e>=100. && L.cer.asum_c>=50.";
  TString cuteR = "R.tr.n==1 && DBB.evtypebits>>1&1 && epR>0.5 && R.sh.e>=200. && R.cer.asum_c>=50.";
  //  TString cuteL = "1"; 
  //TString cuteR = "1";
  
  TString infile = "Run_Chart.txt";
  //TString infile = "C12.txt";
  ifstream runfile(infile);
  ofstream outfile("bad_run.txt");
  Float_t mom, angle;
  TString iL, iR, Target,Kin;

  TCanvas *c1 = new TCanvas("c1","c1",1000,500);
  c1->Divide(2,1);

  TString rootfile; Int_t runnum=0;
  Float_t peak,peakL, peakR, peakL_err,peakR_err; Int_t zpeak;

  TFile* ff = new TFile("Point_Ytg.root","recreate");
  //TFile* ff = new TFile("Point_Ytg.root","update");
  TTree* T1  = new TTree("T1","A Pointing Tree");

  T1->Branch("runnum",   &runnum,   "runnum/I");
  T1->Branch("angle",    &angle,    "angle/F");
  T1->Branch("mom",      &mom,      "mom/F");
  T1->Branch("peakL",    &peakL,    "peakL/F");
  T1->Branch("peakL_err",&peakL_err,"peakL_err/F");
  T1->Branch("peakR",    &peakR,    "peakR/F");
  T1->Branch("peakR_err",&peakR_err,"peakR_err/F");
   TString ROOT_DIR = "/work/halla/e08014/disk1/Rootfiles";
   int NR=0, NL=0;
  while(!(runfile.eof())){

	  runfile >> runnum >> Kin >> mom >> angle >> Target >> iL >> iR;

	  if(Target=="C12"){	
		  rootfile = Form("%s/e08014_%d.root",ROOT_DIR.Data(),runnum);

		  if(gSystem->FindFile(ROOT_DIR,rootfile)){
			  cerr << "Working on Run#"<<runnum<<endl;
			  TFile *file1 = new TFile(rootfile,"r"); 
			  TTree *T = (TTree*)gDirectory->Get("T");
			  if(!(file1->IsZombie())){

				  peakL = -1000.00;
				  peakL_err = -1000.00;
				  peakR = -1000.00;
				  peakR_err = -1000.00;
					  
				  TH1F *h1 = new TH1F("h1",Form("Run#%d Y_Tg_L",runnum),500,-0.05,0.05);
				  TH1F *h2 = new TH1F("h2",Form("Run#%d Y_Tg_R",runnum),500,-0.05,0.05);

				  //if(iL=="L"){
				  NL = (int)T->GetEntries(cuteL.Data());
				  if(NL>1e3){
					  c1->cd(1);
					  T->Draw("L.gold.y>>h1",TCut(cuteL));
					  zpeak = h1->GetMaximumBin();
					  peak = h1->GetBinCenter(zpeak); 
					  h1->Fit("gaus","Q","",peak-0.005,peak+0.005);
					  peakL = h1->GetFunction("gaus")->GetParameter(1);
					  peakL_err = h1->GetFunction("gaus")->GetParameter(2);

					 // tex1 = new TLatex(0.2,0.7,Form("L: Angle=%4.2f,Mom=%5.3f,Peak=%5.2fmm",angle,mom,peakL*1000));
					  //tex1->SetNDC();
					 // tex1->SetTextSize(0.04);
					 // tex1->SetTextColor(1);
					 // tex1->Draw();
					  c1->Update();
				  }

				  //if(iR=="R"){
				  NR = (int)T->GetEntries(cuteR.Data());
				  if(NR>1e3){
					  c1->cd(2);
					  T->Draw("R.gold.y>>h2",TCut(cuteR));
					  zpeak = h2->GetMaximumBin();
					  peak = h2->GetBinCenter(zpeak); 
					  h2->Fit("gaus","Q","",peak-0.005,peak+0.005);
					  peakR = h2->GetFunction("gaus")->GetParameter(1);
					  peakR_err = h2->GetFunction("gaus")->GetParameter(2);

					 // tex2 = new TLatex(0.2,0.7,Form("R: Angle=%4.2f,Mom=%5.3f,Peak=%5.2fmm",angle,mom,peakR*1000));
					 // tex2->SetNDC();
					 // tex2->SetTextSize(0.04);
					 // tex2->SetTextColor(1);
					 // tex2->Draw();
					  c1->Update();
				  }
				  
				  c1->Print(Form("./Ytg/Ytg_%d.pdf",runnum));  
				  delete h1;
				  delete h2;
				//  delete tex1;
				//  delete tex2;
				  file1->Close();
				  ff->cd(); T1->Fill();// h1->Write();	
			  }
			  else
				  outfile << runnum << endl;
		  }
		  else
			  outfile << runnum << endl;

	  }
  }
  ff->Write(); ff->Close();
}

