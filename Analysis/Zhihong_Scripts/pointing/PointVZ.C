
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

  //Beam Position for X>2 Experiment
const Double_t Beam_X = -2.668; //mm
const Double_t Beam_Y = 3.022; //mm

const Double_t L_HRS = 8.458*1000.; //mm
const Double_t DegToRad = 3.1415926 / 180.00;
const Double_t RadToDeg = 180.00 / 3.1415926;

//void Check_Point()
int main()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  TString cuteL = "L.tr.n==1 && DBB.evtypebits>>3&1 && epL>0.5 && L.prl2.e>=100. && L.cer.asum_c>=50.";
  TString cuteR = "R.tr.n==1 && DBB.evtypebits>>1&1 && epR>0.5 && R.sh.e>=200. && R.cer.asum_c>=50.";
  //  TString cuteL = "1"; 
  //TString cuteR = "1";
  
  //TString infile = "Run_Chart.txt";
  TString infile = "C12_N.txt";
  ifstream runfile(infile);
  ofstream outfile("bad_run.txt");
  Float_t mom, angle, AngleL, AngleR, AngleL_Err, AngleR_Err,peak;
  TString iL, iR, Target,Kin;
  TString rootfile; Int_t runnum=0, zpeak=0;
  float YtgR, YtgR_Err, YtgL, YtgL_Err, RpzR, RpzL, RpzR_Err, RpzL_Err, DL, DR;

  TFile* ff = new TFile("Point_VZ_C12.root","recreate");
  //TFile* ff = new TFile("Point_VZ.root","update");
  TTree* T  = new TTree("T","A Pointing Tree");

  T->Branch("runnum",   &runnum,   "runnum/I");
  T->Branch("angle",    &angle,    "angle/F");
  T->Branch("AngleL",   &AngleL,   "AngleL/F");
  T->Branch("AngleR",   &AngleR,   "AngleR/F");
  T->Branch("AngleL_Err",   &AngleL_Err,   "AngleL_Err/F");
  T->Branch("AngleR_Err",   &AngleR_Err,   "AngleR_Err/F");
  T->Branch("mom",      &mom,      "mom/F");
  T->Branch("DL",       &DL,       "DL/F");
  T->Branch("DR",       &DR,       "DR/F");
  T->Branch("YtgL",     &YtgL,     "YtgL/F");
  T->Branch("YtgL_Err", &YtgL_Err, "YtgL_Err/F");
  T->Branch("YtgR",     &YtgR,     "YtgR/F");
  T->Branch("YtgR_Err", &YtgR_Err, "YtgR_Err/F");
  T->Branch("RpzL",     &RpzL,     "RpzL/F");
  T->Branch("RpzL_Err", &RpzL_Err, "RpzL_Err/F");
  T->Branch("RpzR",     &RpzR,     "RpzR/F");
  T->Branch("RpzR_Err", &RpzR_Err, "RpzR_Err/F");

  TString ROOT_DIR = "/work/halla/e08014/disk1/Rootfiles";
  int NR=0, NL=0;
 
  TCanvas *c1 = new TCanvas("c1","c1",1000,500);
  c1->Divide(2,1);
  TCanvas *c2 = new TCanvas("c2","c2",1000,500);
  c2->Divide(2,1);

  while(!(runfile.eof())){

	  runfile >> runnum >> Kin >> mom >> angle >> Target >> iL >> iR;

	  //if(Target=="C12"){	
	  if(1){	
		  rootfile = Form("%s/e08014_%d.root",ROOT_DIR.Data(),runnum);

		  if(gSystem->FindFile(ROOT_DIR,rootfile)){
			  cerr << "Working on Run#"<<runnum<<endl;
			  TFile *file1 = new TFile(rootfile,"r"); 
			  TTree *T1 = (TTree*)gDirectory->Get("T");
			  if(!(file1->IsZombie())){

				  YtgL = -1000.00; YtgL_Err = -1000.00; YtgR = -1000.00; YtgR_Err = -1000.00;
				  RpzL = -1000.00; RpzL_Err = -1000.00; RpzR = -1000.00; RpzR_Err = -1000.00;
				  AngleR = -1000.00; AngleL = -1000.00;	 AngleR_Err = -1000.00; AngleL_Err = -1000.00;	 
                  DL = -1000.0; DR = -1000.0;

				  TH1F *h1 = new TH1F("h1",Form("Run#%d Ytg_L",runnum),500,-0.03,0.03);
				  TH1F *h2 = new TH1F("h2",Form("Run#%d VZ_L",runnum),500,-0.03,0.03);
				  TH1F *h3 = new TH1F("h3",Form("Run#%d Ytg_R",runnum),500,-0.03,0.03);
				  TH1F *h4 = new TH1F("h4",Form("Run#%d VZ_R",runnum),500,-0.03,0.03);
	
				  //if(iL=="L"){
				  NL = (int)T1->GetEntries(cuteL.Data());
				  if(NL>1e3){
					  c1->cd(1);
					  T1->Draw("L.gold.y>>h1",TCut(cuteL));
					  zpeak = h1->GetMaximumBin();
					  peak = h1->GetBinCenter(zpeak); 
					  h1->Fit("gaus","Q","",peak-0.005,peak+0.005);
					  YtgL = 1000. * h1->GetFunction("gaus")->GetParameter(1);
					  YtgL_Err = 1000. * h1->GetFunction("gaus")->GetParameter(2);

					  c1->cd(2);
					  T1->Draw("RctPtL.z>>h2",TCut(cuteL));
					  zpeak = h1->GetMaximumBin();
					  peak = h1->GetBinCenter(zpeak); 
					  h2->Fit("gaus","Q","",peak-0.005,peak+0.005);
					  RpzL = 1000. * h2->GetFunction("gaus")->GetParameter(1);
					  RpzL_Err = 1000. * h2->GetFunction("gaus")->GetParameter(2);

                      DL = Beam_X * cos(angle*DegToRad) - RpzL * sin(angle*DegToRad) - YtgL;   
                      AngleL = angle - DL/L_HRS*RadToDeg;
                      AngleL_Err = AngleL / DL * sqrt(pow(YtgL_Err,2) + pow(RpzL_Err*sin(angle*DegToRad),2));

					  c1->Update();
					  c1->Print(Form("./LHRS/VZ_%d.png",runnum));  
				  }

				  //if(iR=="R"){
				  NR = (int)T1->GetEntries(cuteR.Data());
				  if(NR>1e3){
					  c2->cd(1);
					  T1->Draw("R.gold.y>>h3",TCut(cuteR));
					  zpeak = h3->GetMaximumBin();
					  peak = h3->GetBinCenter(zpeak); 
					  h3->Fit("gaus","Q","",peak-0.005,peak+0.005);
					  YtgR = 1000. * h3->GetFunction("gaus")->GetParameter(1);
					  YtgR_Err = 1000. * h3->GetFunction("gaus")->GetParameter(2);

					  c2->cd(2);
					  T1->Draw("RctPtR.z>>h4",TCut(cuteR));
					  zpeak = h4->GetMaximumBin();
					  peak = h4->GetBinCenter(zpeak); 
					  h4->Fit("gaus","Q","",peak-0.005,peak+0.005);
					  RpzR = 1000. * h4->GetFunction("gaus")->GetParameter(1);
					  RpzR_Err = 1000. * h4->GetFunction("gaus")->GetParameter(2);

					  DR = Beam_X*cos(angle*DegToRad) - RpzR*sin(angle*DegToRad) - YtgR;   
                      AngleR = angle - DR/L_HRS*RadToDeg;
                      AngleR_Err = AngleR / DR * sqrt(pow(YtgR_Err,2) + pow(RpzR_Err*sin(angle*DegToRad),2));

					  c2->Update();
					  c2->Print(Form("./RHRS/VZ_%d.png",runnum));  
				  }
				  
				  delete h1;
				  delete h2;
				  delete h3;
				  delete h4;
				  file1->Close();
				  ff->cd(); T->Fill();// h1->Write();	
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

