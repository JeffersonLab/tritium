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
#include <TLegend.h>
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
#include <TRandom3.h>
//#include <TMatrix.h>
/*}}}*/

/*Correction{{{*/
double Correction(const TString& aTarget, const double aVZ, const double aCurrent, const double aRho){
	TString Boiling_File=Form("%s_Boiling.out",aTarget.Data());
	const int tBin = 60; //Bin Number used in Boiling Study, all equal to 60
	const double step = 0.24/tBin; //+- 0.12 used on VZ cut
	double *tVZ = new double[tBin];
	double *tSlop = new double[tBin];

	ifstream infile(Boiling_File); 
	TString aCom; infile >> aCom >> aCom; //Reading the first line

	for(int i=0;i<tBin;i++){
		infile >> tVZ[i] >> tSlop[i];
	}
	infile.close();

	int bin = (int) ((aVZ + 0.12)/step);
	double slop = ( tSlop[bin+1] - (tSlop[bin+1]-tSlop[bin])/(tVZ[bin+1]-tVZ[bin]) * (tVZ[bin+1]-aVZ) );
	double rho = aRho * (1 + slop*aCurrent);

	//cerr<<Form("Bin=%d, Slop=%f, VZ=%f/%f, I=%f, Rho0=%f, Rho=%f",bin,slop,aVZ,tVZ[bin], aCurrent, aRho, rho)<<endl;

	delete tVZ; delete tSlop;
	return rho;
}

/*}}}*/

/*Correction-Inverse{{{*/
double Correction_Inverse(const TString& aTarget, const double aVZ, const double aCurrent, const double aRho){
	TString Boiling_File=Form("%s_Boiling.out",aTarget.Data());
	const int tBin = 60; //Bin Number used in Boiling Study, all equal to 60
	const double step = 0.24/tBin; //+- 0.12 used on VZ cut
	double *tVZ = new double[tBin];
	double *tSlop = new double[tBin];

	ifstream infile(Boiling_File); 
	TString aCom; infile >> aCom >> aCom; //Reading the first line

	for(int i=0;i<tBin;i++){
		infile >> tVZ[i] >> tSlop[i];
	}
	infile.close();

	int bin = (int) ((aVZ + 0.12)/step);
	double slop = ( tSlop[bin+1] - (tSlop[bin+1]-tSlop[bin])/(tVZ[bin+1]-tVZ[bin]) * (tVZ[bin+1]-aVZ) );
	double rho = aRho / (1 + slop*aCurrent);

	cerr<<Form("Bin=%d, Slop=%f, VZ=%f/%f, I=%f, Rho0=%f, Rho=%f",bin,slop,aVZ,tVZ[bin], aCurrent, aRho, rho)<<endl;

	delete tVZ; delete tSlop;
	return rho;
}

/*}}}*/

/*int main{{{*/
void Check( TString Target_Name){

	/*Define{{{*/
	gROOT->SetStyle("Plain");
	gStyle->SetFillColor(0);
	TString Input_File = ""; 
	Input_File = Form("%s_L_kin5.1_unf.rho",Target_Name.Data()); 
	ifstream infile1(Input_File);
	Input_File = Form("../%s_Yield_L_All.out",Target_Name.Data()); 
	ifstream infile4(Input_File);

	Input_File = Form("%s.rho",Target_Name.Data()); 
	ofstream outfile(Input_File);
	
	const int bin = 60;
	double I1=0.;
	if(Target_Name=="H2"){
		I1=45.0;
	}
	else if(Target_Name=="He3"){
		I1=120.0;
	}
	else if(Target_Name=="He4"){
		I1 = 95.0; 
	}
	else
		cerr<<"***ERROR, Unknow Target!!!"<<endl;

	double *VZ1 = new double[bin];
	double *VZf = new double[bin];
	double *Rho1 = new double[bin];
	double *Rho1_Err = new double[bin];
	double *Rho0 = new double[bin];
	double *Rho0_Err = new double[bin];
	double *Rho0f = new double[bin];
	double *Rho0f_Err = new double[bin];
	double *Rho1f = new double[bin];
	double *Rho1f_Err = new double[bin];
	double *Zero = new double[bin];
	double *BF = new double[bin];
	double *BF_Err = new double[bin];

	double rho1_sum=0., rho0_sum=0.;
    double rho0f_sum=0.,rho1f_sum=0.;
    TString acom; infile4>>acom>>acom>>acom>>acom>>acom;
	/*}}}*/

	double aa=0;
	for(int i=0;i<bin;i++){
		infile1 >> VZ1[i] >> Rho1[i];	Rho1_Err[1]=sqrt(Rho1[i]);
		infile4 >> VZf[i] >> Rho0f[i] >> Rho0f_Err[i] >> BF[i] >> BF_Err[i];	

		Zero[i]=0.0;
		Rho0[i] = Correction_Inverse(Target_Name.Data(), VZ1[i], I1, Rho1[i] ); 
		Rho0_Err[i] = Rho0[i]/Rho1[i]*Rho1_Err[i];

		Rho1f[i] = Correction(Target_Name.Data(), VZf[i], I1, Rho0f[i] ); 
		Rho1f_Err[i] = Rho1f[i]/Rho0f[i]*Rho0f_Err[i];

		if(VZ1[i]>-0.075 && VZ1[i]<0.075){
           rho0_sum+=Rho0[i];
           rho1_sum+=Rho1[i];
		}
		else{
			Rho0[i]=0.0;
			Rho1[i]=0.0;
			Rho0_Err[i]=0.0;
			Rho1_Err[i]=0.0;
		}

		if(VZf[i]>-0.075 && VZf[i]<0.075){
           rho0f_sum+=Rho0f[i];
           rho1f_sum+=Rho1f[i];
		}
		else{
			Rho0f[i]=0.0;
			Rho1f[i]=0.0;
			Rho0f_Err[i]=0.0;
			Rho1f_Err[i]=0.0;
		}
	}
	infile1.close(); 

	cerr<<Form("Rho=%f, Rho1=%f, Rhof=%f, Rho1f=%f", rho0_sum, rho1_sum, rho0f_sum, rho1f_sum)<<endl;

	outfile<<Form("%12s %12s %12s %12s %12s", "VZ","Y0","Y0_Err","BF","BF_Err")<<endl; 
	for(int i=0;i<bin;i++){
		Rho0[i]/=rho0_sum;
		Rho0_Err[i]/=rho0_sum;
		Rho1[i]/=rho1_sum;
		Rho1_Err[i]/=rho1_sum;

		Rho0f[i]/=rho0f_sum;
		Rho0f_Err[i]/=rho0f_sum;
		Rho1f[i]/=rho1f_sum;
		Rho1f_Err[i]/=rho1f_sum;
		
		outfile<<Form("  %12.8f %12.8f %12.8f %12.8f %12.8f",VZf[i],Rho0f[i], Rho0f_Err[i], BF[i],BF_Err[i])<<endl; 
		}	
    outfile.close();

	TCanvas *c1 = new TCanvas("c1","c1",1200,700);
	c1->cd();
	TH2F *h1 = new TH2F("h1","",300,-0.12,0.12,300,-0.01,0.04);
	h1->SetStats(kFALSE);
	h1->SetXTitle("z_{react} (m)");
	h1->GetXaxis()->CenterTitle(1);
	h1->GetXaxis()->SetTitleFont(32);
	h1->GetXaxis()->SetTitleSize(0.06);
	h1->SetYTitle("#rho_{norm}");
	h1->GetYaxis()->CenterTitle(1);
	h1->GetYaxis()->SetTitleFont(32);
	h1->GetYaxis()->SetTitleSize(0.06);
	h1->GetYaxis()->SetTitleOffset(0.8);
	h1->SetLineWidth(1.5);
	h1->Draw();

	TGraphErrors *ex = new TGraphErrors(bin, VZ1, Rho0, Zero, Rho0_Err);
	ex->SetMarkerStyle(20);
	ex->SetMarkerColor(1);
	ex->SetLineColor(1);
	ex->SetLineWidth(1.5);

	TGraphErrors *ex4 = new TGraphErrors(bin, VZf, Rho0f, Zero, Rho0f_Err);
	ex4->SetMarkerStyle(20);
	ex4->SetMarkerColor(6);
	ex4->SetLineColor(6);
	ex4->SetLineWidth(1.5);
	

	TGraphErrors *ex1 = new TGraphErrors(bin, VZ1, Rho1, Zero, Rho1_Err);
	ex1->SetMarkerStyle(21);
	ex1->SetMarkerColor(6);
	ex1->SetLineColor(2);
	ex1->SetLineWidth(1.5);

	TGraphErrors *ex41 = new TGraphErrors(bin, VZf, Rho1f, Zero, Rho1f_Err);
	ex41->SetMarkerStyle(28);
	ex41->SetMarkerColor(4);
	ex41->SetLineColor(2);
	ex41->SetLineWidth(1.5);
	
	TLegend *l1  = new TLegend(0.33,0.2,0.72,0.45,Form("%s Density Profile:",Target_Name.Data()));
	l1->SetTextSize(0.03);
	l1->SetTextFont(32);

	c1->Clear(); h1->Draw();
	ex->Draw("l");
	l1->AddEntry(ex, Form("#rho_{0} calculated from #rho(%d uA) ",(int)(I1)),"l");
	ex4->Draw("p");
	l1->AddEntry(ex4, Form("#rho_{0} from Boiling Fitting"),"p");

	ex1->Draw("l");
	l1->AddEntry(ex1, Form("#rho(%d uA)  from data",(int) (I1)),"l");
	ex41->Draw("p");
	l1->AddEntry(ex41, Form("#rho(%d uA) calculated from Boiling Fitting",(int) (I1)),"p");

	l1->Draw();
//	c1->Print(Form("%s_Check_%d.png",Target_Name.Data(), (int) (I1)));
//	c1->Print(Form("%s_Check_%d.pdf",Target_Name.Data(), (int) (I1)));

	delete VZ1;
	delete VZf;
	delete Rho1;
	delete Rho0;  
	delete Rho1_Err;
	delete Rho0_Err;  
	delete Rho0f;  
	delete Rho1f;  
	delete Rho0f_Err;  
	delete Rho1f_Err;  
	delete Zero;
    delete BF;
    delete BF_Err;	
}
/*}}}*/

