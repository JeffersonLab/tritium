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
#include <TLatex.h>
#include <TVirtualFitter.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TCut.h>
#include <TMultiGraph.h>
#include <TCutG.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>
//#include <TMatrix.h>
/*}}}*/

void Fitting(){
	TString RootFile_Name;
	TString Target, Arm;
	gSystem->Exec("ls -l *.root");
	cerr<<" --- Target (H2,He3,He4) = "; cin >> Target;
	cerr<<" --- Arm (L, R) = "; cin >> Arm;
	RootFile_Name = Form("%s_Yield_%s_All.root",Target.Data(),Arm.Data());
	TFile *f1 = new TFile(RootFile_Name.Data(),"r");
	TTree *T  = (TTree*) gDirectory->Get("T");

	TString OutFile_Name = RootFile_Name;
	OutFile_Name.ReplaceAll(".root",".out");
	ofstream outfile(OutFile_Name);

	TString File_Name1 = RootFile_Name;
	TString File_Name2 = RootFile_Name;
	TString File_Name3 = RootFile_Name;
	TString File_Name4 = RootFile_Name;
	File_Name1.ReplaceAll(".root",".png");
	File_Name2.ReplaceAll(".root",".pdf");
	File_Name3.ReplaceAll(".root","_0.png");
	File_Name4.ReplaceAll(".root","_0.pdf");

	double VZ, I, Y, Y_Err;
	int size, runno, Bin;

	T->SetBranchAddress("runsize",&size);
	T->SetBranchAddress("Bin",&Bin);
	T->SetBranchAddress("RunNo",&runno);
	T->SetBranchAddress("I",&I);
	T->SetBranchAddress("Y",&Y);
	T->SetBranchAddress("Y_Err",&Y_Err);
	T->SetBranchAddress("VZ",&VZ);

	int N = T->GetEntries();
	T->GetEntry(0);
	const int Size = size;
	int temp=20;
	if(Target=="C12"||Target=="Ca40"||Target=="Ca48")
		temp=1;
	const int Bin_Size = temp;

	/*Read Root File{{{*/
	cerr<<" @@@ Reading Root File ..."<<endl;
	double vI[Bin_Size][Size];
	double vI_Err[Bin_Size][Size];
	double vY[Bin_Size][Size];
	double vY_Err[Bin_Size][Size];
	double vVZ[Bin_Size];
	double Y_Max=-1e32,Y_Min=1e32, X_Max=-1000,X_Min=1000;
	int k=0;
	for(int j=0;j<Bin_Size; j++){
		k=0;
		for(int i=0;i<N;i++){
			T->GetEntry(i);
			if(Y>-0.5 && Bin==j){
				vI[j][k]=I;
				vI_Err[j][k]=0.;
				vY[j][k]=Y;
				vY_Err[j][k]=Y_Err;

				vVZ[j] = VZ;
				if(X_Max<vI[j][k])
					X_Max = vI[j][k];
				if(X_Min>vI[j][k])
					X_Min = vI[j][k];
				if(Y_Max<vY[j][k])
					Y_Max = vY[j][k];
				if(Y_Min>vY[j][k])
					Y_Min = vY[j][k];
				// Y_Max = Y_Max?vY[j][k]:Y_Max<vY[j][k];
				// Y_Min = Y_Min?vY[j][k]:Y_Min>vY[j][k];
				// X_Max = X_Max?vI[j][k]:X_Max<vI[j][k];
				// X_Min = X_Min?vI[j][k]:X_Min>vI[j][k];

				cerr<<Form("--- Bin = %d, VZ=%6.4f, I=%6.4f, Y=%e",
						Bin, VZ, I, Y)<<endl;
				k++;
				if(k>Size)
					cerr<<"*** Something wrong?!"<<endl;
			}
		}
	}
	/*}}}*/

	/*Plot and Fit{{{*/
	cerr<<" @@@ Making Plots and Fitting ..."<<endl;
	cerr<<Form(" ---- X_Min=%f, X_Max=%f, Y_Min=%e, Y_Max=%e",X_Min,X_Max, Y_Min, Y_Max)<<endl;

	TCanvas *c3 = new TCanvas("c3","c3",1200,700);
	c3->cd();
	TH2F *h3 = new TH2F("h3","",300,X_Min-2.0,X_Max+2.0,300,0.0,Y_Max/Y_Min*1.12);
	h3->SetStats(kFALSE);
	h3->SetXTitle("I (uA)");
	h3->GetXaxis()->CenterTitle(1);
	h3->GetXaxis()->SetTitleFont(32);
	h3->GetXaxis()->SetTitleSize(0.06);
	h3->GetYaxis()->SetTitleOffset(0.8);
	h3->SetYTitle("Yield_{normalized}");
	h3->GetYaxis()->CenterTitle(1);
	h3->GetYaxis()->SetTitleFont(32);
	h3->GetYaxis()->SetTitleSize(0.06);
	h3->GetYaxis()->SetTitleOffset(0.8);

	gStyle->SetOptStat(1);
	gStyle->SetOptFit(1);
	double S[Bin_Size],S_Err[Bin_Size];
	double Y0[Bin_Size],Y0_Err[Bin_Size];
	double BF[Bin_Size],BF_Err[Bin_Size];
	double Zero[Bin_Size];

	for(int l=0; l<Bin_Size; l++){
		for(int m=0; m<Size; m++){
			vY[l][m]/=Y_Min;   	
			vY_Err[l][m]/=Y_Min;   	
		}
	}
	TLatex *tex;
	outfile<<Form("%12s %12s %12s %12s %12s", "VZ","Y0","Y0_Err","BF","BF_Err")<<endl; 
	Y_Max=-1e32,Y_Min=1e32;
	for(int l=0; l<Bin_Size; l++){
		h3->Draw();
		for(int m=0;m<Size;m++)
			cerr<<Form(" VZ=%6.4f, I = %6.3f, Y = %e", vVZ[l], vI[l][m],vY[l][m])<<endl;

		TGraphErrors *ex = new TGraphErrors(Size, vI[l], vY[l], vI_Err[l], vY_Err[l]);
		ex->SetMarkerStyle(25);
		ex->SetMarkerColor(4);
		ex->SetLineColor(4);
		ex->Draw("P");
		ex->Fit("pol1","IM","",X_Min,X_Max);

		if(l==0){//For C12 only
//	if(l==30){
			tex=new TLatex(0.25, 0.6, Form("%s Boiling Effect Fit at z_{react} = %4.2f cm on HRS-%s",Target.Data(), vVZ[l]*100, Arm.Data()));
			tex->SetNDC();
			tex->SetTextFont(32);
			tex->SetTextSize(0.05);
			tex->Draw();
			c3->Print(File_Name3.Data());
			c3->Print(File_Name4.Data());
		}
		Y0[l]   = ex->GetFunction("pol1")->GetParameter(0) ;
		Y0_Err[l]   = ex->GetFunction("pol1")->GetParError(0) ;
		S[l] = ex->GetFunction("pol1")->GetParameter(1);
		S_Err[l] = ex->GetFunction("pol1")->GetParError(1);
		BF[l] = S[l]/Y0[l];
		BF_Err[l] = BF[l] * sqrt( pow(Y0_Err[l]/Y0[l],2) + pow(S_Err[l]/S[l],2) );
		Zero[l] = 0.0;

		if(Y_Max<S[l])
			Y_Max = S[l];
		if(Y_Min>S[l])
			Y_Min = S[l];

		cerr<<Form("  VZ=%f, Slop=%f",vVZ[l],S[l])<<endl; 
		outfile<<Form("  %12.8f %12.8f %12.8f %12.8f %12.8f",vVZ[l],Y0[l], Y0_Err[l], BF[l],BF_Err[l])<<endl; 

		c3->Clear();
	}
	c3->Close();

	TCanvas *c4 = new TCanvas("c4","c4",1200,700);
	c4->cd();
	TH2F *h4 = new TH2F("h4","",300, -0.14, 0.14, 300, -0.011, 0.001);
	h4->SetStats(kFALSE);
	h4->SetXTitle("z_{react} (m)");
	h4->GetXaxis()->CenterTitle(1);
	h4->GetXaxis()->SetTitleFont(32);
	h4->GetXaxis()->SetTitleSize(0.06);
	h4->GetXaxis()->SetTitleOffset(0.7);
	h4->SetYTitle("BF");
	h4->GetYaxis()->CenterTitle(1);
	h4->GetYaxis()->SetTitleFont(32);
	h4->GetYaxis()->SetTitleSize(0.04);
	h4->GetYaxis()->SetTitleOffset(1.2);
	h4->Draw();

	TGraphErrors *ex1 = new TGraphErrors(Bin_Size, vVZ, BF, Zero, BF_Err);
	ex1->SetMarkerStyle(20);
	ex1->SetMarkerColor(4);
	ex1->SetLineColor(4);
	ex1->Draw("P");
	tex=new TLatex(0.3, 0.2, Form("%s Boiling Effect on HRS-%s",Target.Data(), Arm.Data()));
	tex->SetNDC();
	tex->SetTextFont(32);
	tex->SetTextSize(0.05);
	tex->Draw();

	c4->Print(File_Name1.Data());
	c4->Print(File_Name2.Data());
	/*}}}*/
}
