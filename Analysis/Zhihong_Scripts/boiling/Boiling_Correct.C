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

	cerr<<Form("Bin=%d, Slop=%f, VZ=%f/%f, I=%f, Rho0=%f, Rho=%f",bin,slop,aVZ,tVZ[bin], aCurrent, aRho, rho)<<endl;

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

	//   cerr<<Form("Bin=%d, Slop=%f, VZ=%f/%f, I=%f, Rho0=%f, Rho=%f",bin,slop,aVZ,tVZ[bin], aCurrent, aRho, rho)<<endl;

	delete tVZ; delete tSlop;
	return rho;
}

/*}}}*/

/*int main{{{*/
void Boiling_Correct( TString Target_Name){

	/*Define{{{*/
	gROOT->SetStyle("Plain");
	gStyle->SetFillColor(0);
	TString Input_File = ""; 
	Input_File = Form("%s_L_kin5.0_unf_min.rho",Target_Name.Data()); 
	ifstream infile1(Input_File);
	Input_File = Form("%s_L_kin5.0_unf_mid.rho",Target_Name.Data()); 
	ifstream infile2(Input_File);
	Input_File = Form("%s_L_kin5.0_unf_max.rho",Target_Name.Data()); 
	ifstream infile3(Input_File);
	Input_File = Form("../%s_Yield_L_All.out",Target_Name.Data()); 
	ifstream infile4(Input_File);

	Input_File = Form("%s.rho",Target_Name.Data()); 
	ofstream outfile(Input_File);
	
	const int bin = 60;
	double I1=0., I2=0., I3=0.;
	if(Target_Name=="H2"){
		I1 = 15.0; I2 = 30.0; I3=45.0;
	}
	else if(Target_Name=="He3"){
		I1 = 20.0; I2 = 60.0; I3=120.0;
	}
	else if(Target_Name=="He4"){
		I1 = 15.0; I2 = 60.0; I3=95.0;
	}
	else
		cerr<<"***ERROR, Unknow Target!!!"<<endl;

	double *VZ1 = new double[bin];
	double *VZ2 = new double[bin];
	double *VZ3 = new double[bin];
	double *VZf = new double[bin];
	double *Rho1 = new double[bin];
	double *Rho2 = new double[bin];
	double *Rho3 = new double[bin];
	double *Rho0C = new double[bin];
	double *Rho2C = new double[bin];
	double *Rho3C = new double[bin];
	double *Rho0f = new double[bin];
	double *Rho1f = new double[bin];
	double *Rho2f = new double[bin];
	double *Rho3f = new double[bin];
	double *Rho1_Err = new double[bin];
	double *Rho2_Err = new double[bin];
	double *Rho3_Err = new double[bin];
	double *Rho0C_Err = new double[bin];
	double *Rho2C_Err = new double[bin];
	double *Rho3C_Err = new double[bin];
	double *Rho0f_Err = new double[bin];
	double *Rho1f_Err = new double[bin];
	double *Rho2f_Err = new double[bin];
	double *Rho3f_Err = new double[bin];
	double *Zero = new double[bin];

	double rho1_sum=0., rho0_sum=0., rho0c_sum=0.;
    double rho0f_sum=0.,rho1f_sum=0.,rho2f_sum=0.,rho3f_sum=0.;
	double rho2_sum=0., rho2c_sum=0.;
	double rho3_sum=0., rho3c_sum=0.;
    TString acom; infile4>>acom>>acom>>acom>>acom>>acom;
	/*}}}*/

	double aa=0;
	for(int i=0;i<bin;i++){
		infile1 >> VZ1[i] >> Rho1[i];	Rho1_Err[1]=sqrt(Rho1[i]);
		infile2 >> VZ2[i] >> Rho2[i];	Rho2_Err[1]=sqrt(Rho2[i]);	
		infile3 >> VZ3[i] >> Rho3[i];	Rho3_Err[1]=sqrt(Rho3[i]);	
		infile4 >> VZf[i] >> Rho0f[i] >> Rho0f_Err[i] >> aa >> aa;	

		Zero[i]=0.0;
		Rho0C[i] = Correction_Inverse(Target_Name.Data(), VZ1[i], I1, Rho1[i] ); 
		Rho0C_Err[i] = Rho0C[i]/Rho1[i]*Rho1_Err[i];
		Rho2C[i] = Correction(Target_Name.Data(), VZ2[i], I2, Rho0C[i] ); 
		Rho2C_Err[i] = Rho2C[i]/Rho0C[i]*Rho0C_Err[i];
		Rho3C[i] = Correction(Target_Name.Data(), VZ3[i], I3, Rho0C[i] );
		Rho3C_Err[i] = Rho3C[i]/Rho0C[i]*Rho0C_Err[i];

		Rho1f[i] = Correction(Target_Name.Data(), VZf[i], I1, Rho0f[i] ); 
		Rho1f_Err[i] = Rho1f[i]/Rho0f[i]*Rho0f_Err[i];
		Rho2f[i] = Correction(Target_Name.Data(), VZf[i], I2, Rho0f[i] ); 
		Rho2f_Err[i] = Rho2f[i]/Rho0f[i]*Rho0f_Err[i];
		Rho3f[i] = Correction(Target_Name.Data(), VZf[i], I3, Rho0f[i] ); 
		Rho3f_Err[i] = Rho3f[i]/Rho0f[i]*Rho0f_Err[i];

		if(VZ1[i]>-0.075 && VZ1[i]<0.075){
           rho0_sum+=Rho0C[i];
		}

		if(VZf[i]>-0.075 && VZf[i]<0.075){
           rho0f_sum+=Rho0f[i];
           rho1f_sum+=Rho1f[i];
           rho2f_sum+=Rho2f[i];
           rho3f_sum+=Rho3f[i];
		}
		else{
			Rho0f[i]=0.0;
			Rho1f[i]=0.0;
			Rho2f[i]=0.0;
			Rho3f[i]=0.0;
		}
		rho1_sum+=Rho1[i]; rho0c_sum+=Rho0C[i];
		rho2_sum+=Rho2[i]; rho2c_sum+=Rho2C[i];
		rho3_sum+=Rho3[i]; rho3c_sum+=Rho3C[i];
	}
	infile1.close(); infile2.close(); infile3.close();

    double rho00[bin],rho0f[bin];
	outfile<<Form("%8s %12s %12s %12s %12s",
			"VZ1","Rho(0)",Form("Rho(%d)",(int)(I1)),Form("Rho(%d)",(int)(I2)),Form("Rho(%d)",(int)(I3)))<<endl;
	for(int i=0;i<bin;i++){
		rho00[i]=Rho0C[i]/rho0_sum;
		Rho2[i]*=rho0c_sum/rho2_sum;
		Rho2_Err[i]*=rho0c_sum/rho2_sum;
		Rho2C[i]*=(rho0c_sum/rho2_sum) * (rho2_sum/rho2c_sum);
		Rho2C_Err[i]*=(rho0c_sum/rho2_sum) * (rho2_sum/rho2c_sum);
		Rho3[i]*=rho0c_sum/rho3_sum;
		Rho3_Err[i]*=rho0c_sum/rho3_sum;
		Rho3C[i]*=(rho0c_sum/rho3_sum) * (rho3_sum/rho3c_sum);
		Rho3C_Err[i]*=(rho0c_sum/rho3_sum) * (rho3_sum/rho3c_sum);

		Rho0f[i]/=rho0f_sum;
		Rho1f[i]/=rho1f_sum;
		Rho2f[i]/=rho2f_sum;
		Rho3f[i]/=rho3f_sum;
//		outfile<<VZ1[i]<<" "<<Rho0f[i]<<" "<<Rho1f[i]<<" "<<Rho2f[i]<<" "<<Rho3f[i]<<endl;
		outfile<<Form("%8.4f %12.4e %12.4e %12.4e %12.4e",VZ1[i],Rho0f[i],Rho1f[i],Rho2f[i],Rho3f[i])<<endl;
		Rho0f[i]*=rho0c_sum;
		Rho1f[i]*=rho0c_sum;
		Rho2f[i]*=rho0c_sum;
		Rho3f[i]*=rho0c_sum;
		Rho0f_Err[i]*=rho0c_sum/rho0f_sum;
		Rho1f_Err[i]*=rho0c_sum/rho0f_sum;
		Rho2f_Err[i]*=rho0c_sum/rho0f_sum;
		Rho3f_Err[i]*=rho0c_sum/rho0f_sum;
		}	
    outfile.close();

	TCanvas *c1 = new TCanvas("c1","c1",1200,700);
	c1->cd();
	TH2F *h1 = new TH2F("h1","",300,-0.14,0.14,300,-0.01,1.5);
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

	TGraphErrors *ex = new TGraphErrors(bin, VZ1, Rho0C, Zero, Rho0C_Err);
	ex->SetMarkerStyle(20);
	ex->SetMarkerColor(1);
	ex->SetLineColor(1);
	ex->SetLineWidth(1.5);

	TGraphErrors *ex4 = new TGraphErrors(bin, VZf, Rho0f, Zero, Rho0f_Err);
	ex4->SetMarkerStyle(20);
	ex4->SetMarkerColor(6);
	ex4->SetLineColor(6);
	ex4->SetLineWidth(1.5);
	
	TGraphErrors *ex41 = new TGraphErrors(bin, VZf, Rho1f, Zero, Rho1f_Err);
	ex41->SetMarkerStyle(28);
	ex41->SetMarkerColor(4);
	ex41->SetLineColor(2);
	ex41->SetLineWidth(1.5);
	
	TGraphErrors *ex2 = new TGraphErrors(bin, VZ2, Rho2, Zero, Rho2_Err);
	ex2->SetMarkerStyle(21);
	ex2->SetMarkerColor(6);
	ex2->SetLineColor(2);
	ex2->SetLineWidth(1.5);
	TGraphErrors *ex22 = new TGraphErrors(bin, VZ2, Rho2C, Zero, Rho2C_Err);
	ex22->SetMarkerStyle(25);
	ex22->SetMarkerColor(2);
	ex22->SetLineColor(4);
	ex22->SetLineWidth(1.5);
	
	TGraphErrors *ex42 = new TGraphErrors(bin, VZf, Rho2f, Zero, Rho2f_Err);
	ex42->SetMarkerStyle(21);
	ex42->SetMarkerColor(4);
	ex42->SetLineColor(2);
	ex42->SetLineWidth(1.5);
	
	TGraphErrors *ex3 = new TGraphErrors(bin, VZ3, Rho3, Zero, Rho3_Err);
	ex3->SetMarkerStyle(22);
	ex3->SetMarkerColor(6);
	ex3->SetLineColor(2);
	ex3->SetLineWidth(1.5);
	TGraphErrors *ex33 = new TGraphErrors(bin, VZ3, Rho3C, Zero, Rho3C_Err);
	ex33->SetMarkerStyle(26);
	ex33->SetMarkerColor(2);
	ex33->SetLineColor(4);
	ex33->SetLineWidth(1.5);
	
	TGraphErrors *ex43 = new TGraphErrors(bin, VZf, Rho3f, Zero, Rho3f_Err);
	ex43->SetMarkerStyle(21);
	ex43->SetMarkerColor(4);
	ex43->SetLineColor(2);
	ex43->SetLineWidth(1.5);

	TLegend *l1  = new TLegend(0.33,0.2,0.72,0.45,Form("%s Density Profile:",Target_Name.Data()));
	l1->SetTextSize(0.03);
	l1->SetTextFont(32);

	if(1){
		c1->Clear(); h1->Draw();
		ex->Draw("l");
		l1->AddEntry(ex, Form("#rho_{0} calculated from #rho(%d uA) ",(int)(I1)),"l");
		ex4->Draw("p");
		l1->AddEntry(ex4, Form("#rho_{0} from Boiling Fitting"),"p");

		ex2->Draw("l");
		ex22->Draw("l");
		l1->AddEntry(ex2, Form("#rho(%d uA)  from data",(int) (I2)),"l");
		l1->AddEntry(ex22, Form("#rho(%d uA) calculated from #rho_{0}",(int) (I2)),"l");
		ex42->Draw("p");
		l1->AddEntry(ex42, Form("#rho(%d uA) calculated from Boiling Fitting",(int) (I2)),"p");

		l1->Draw();
		c1->Print(Form("%s_Check_%d.png",Target_Name.Data(), (int) (I2)));
		c1->Print(Form("%s_Check_%d.pdf",Target_Name.Data(), (int) (I2)));
	}

	TLegend *l2  = new TLegend(0.33,0.2,0.72,0.45,Form("%s Density Profile:",Target_Name.Data()));
	l2->SetTextSize(0.03);
	l2->SetTextFont(32);

	if(1){
		c1->Clear(); h1->Draw();
		ex->Draw("l");
		l2->AddEntry(ex, Form("#rho_{0} calculated from #rho(%d uA) ",(int)(I1)),"l");
		ex4->Draw("p");
		l2->AddEntry(ex4, Form("#rho_{0} from Boiling Fitting"),"p");

		ex3->Draw("l");
		ex33->Draw("l");
		l2->AddEntry(ex3, Form("#rho(%d uA) from data",(int) (I3)),"l");
		l2->AddEntry(ex33, Form("#rho(%d uA) calculated from #rho_{0}",(int) (I3)),"l");
		ex43->Draw("p");
		l2->AddEntry(ex43, Form("#rho(%d uA) calculated from Boiling Fitting",(int) (I3)),"p");

		l2->Draw();
		c1->Print(Form("%s_Check_%d.png",Target_Name.Data(), (int) (I3)));
		c1->Print(Form("%s_Check_%d.pdf",Target_Name.Data(), (int) (I3)));
	}

	delete VZ1;delete VZ2; delete VZ3;
	delete Rho1; delete Rho2; delete Rho3;
	delete Rho0C; delete Rho2C; delete Rho3C;
}
/*}}}*/

