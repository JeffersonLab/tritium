#include "XGT2_Main.h"
#include "/work/halla/e08014/disk1/yez/Ep_Cross_Section/SRC/XGT2_XEMC.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
/*Get_Target_Name{{{*/
TString gGet_Target_Name(const TString& aName){
	TString aTarget_Name="MM";

	if(aName.Contains("H2")||aName.Contains("D2"))
		aTarget_Name="H2";
	else if(aName.Contains("He3"))
		aTarget_Name="He3";
	else if(aName.Contains("He4"))
		aTarget_Name="He4";
	else if(aName.Contains("C12"))
		aTarget_Name="C12";
	else if(aName.Contains("Ca40"))
		aTarget_Name="Ca40";
	else if(aName.Contains("Ca48"))
		aTarget_Name="Ca48";
	else
		cerr<<"*** I don't know the target name!!!"<<endl;
	return aTarget_Name;
}
/*}}}*/

/*Get_Kin_Name{{{*/
TString gGet_Kin_Name(const TString& aKin){
	TString aKin_Name="MM";
	if(aKin.Contains("3.1"))
		aKin_Name="3.1";
	else if(aKin.Contains("3.2"))
		aKin_Name="3.2";
	else if(aKin.Contains("4.1"))
		aKin_Name="4.1";
	else if(aKin.Contains("4.2"))
		aKin_Name="4.2";
	else if(aKin.Contains("5.0") && !aKin.Contains("5.05"))
		aKin_Name="5.0";
	else if(aKin.Contains("5.05"))
		aKin_Name="5.05";
	else if(aKin.Contains("5.1"))
		aKin_Name="5.1";
	else if(aKin.Contains("5.2"))
		aKin_Name="5.2";
	else if(aKin.Contains("6.5"))
		aKin_Name="6.5";
	else
		cerr<<"*** I don't know the kinematics name!!!"<<endl;
	return aKin_Name;
}
/*}}}*/

/*Cuts on the edge of Ep{{{*/
// Cut_L_Min Cut_L_Max Cut_R_Min Cut_R_Max
const double E0 = 3.356;
const double Ep_Cut_21_1[4] = {2.80, 3.00,2.80,3.00}; 
const double Ep_Cut_21_2[4] = {2.94, 3.18,2.94,3.18}; 
const double Ep_Cut_23_1[4] = {2.74, 2.94,2.74,2.94}; 
const double Ep_Cut_23_2[4] = {2.92, 3.16,2.92,3.16};
const double Ep_Cut_QE_1[4] = {2.42, 2.58,2.42,2.58}; 
const double Ep_Cut_QE_2[4] = {2.56, 2.74,2.56,2.74}; 
const double Ep_Cut_25_1[4] = {2.70, 2.88,2.70,2.88};
const double Ep_Cut_25_2[4] = {2.88, 3.12,2.88,3.12}; 
const double Ep_Cut_28_1[4] = {2.74, 2.98,2.74,2.98}; 
/*}}}*/

void Merge(const TString& aInputFileName, double* aXbj, double* aR_EX, double* aR_EX_Err, double* aR_MC, double* aR_MC_Err, int* aN);
TString Target0,Target,Com0,Com;
/*Comment{{{*/
/*Format in aFileName:
  Target 1:
  line#1 --> The first  setting of each angle on  L-HRS, like L_Kin3.1 ...
  line#2 --> The second setting of each angle on  L-HRS, like L_Kin3.2 ...
  line#3 --> The first  setting of each angle on  R-HRS, like R_Kin3.1 ...
  line#4 --> The second setting of each angle on  R-HRS, like R_Kin3.2 ...
  Target 2:  
  line#5 --> The first  setting of each angle on  L-HRS, like L_Kin3.1 ...
  line#6 --> The second setting of each angle on  L-HRS, like L_Kin3.2 ...

  line#7 --> The first  setting of each angle on  R-HRS, like R_Kin3.1 ...
  line#8 --> The second setting of each angle on  R-HRS, like R_Kin3.2 ...
  */
/*}}}*/

/*main{{{*/
int main()
{
	/*ROOT Define{{{*/
	// Fill color
	gStyle->SetStatColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetPadColor(0);
	gStyle->SetFrameFillColor(0);

	// Border mode
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetPalette(1,0);

	// Font
	gStyle->SetTextFont(32);
	gStyle->SetLabelFont(32, "xyz");
	gStyle->SetTitleFont(22, "xyz");
	gStyle->SetStatFont(32);
	gStyle->SetTitleFillColor(0);

	// Opt
	//    gStyle->SetOptDate(0);
	//    gStyle->SetOptTitle(1);
	//    gStyle->SetOptStat(1);
	//    gStyle->SetOptFit(1);
	//    // Margin
	// //    gStyle->SetPadLeftMargin(0.1);
	// //    gStyle->SetPadRightMargin(0.1);
	// //    gStyle->SetPadTopMargin(0.1);
	// //    gStyle->SetPadBottomMargin(0.1);


	//void HallA_style() {
	gROOT->SetStyle("Plain");
	gStyle->SetPaperSize(TStyle::kUSLetter);
	gStyle->SetPaperSize(18,22);
	gStyle->SetOptFit(1111);
	gStyle->SetPalette(1);
	gStyle->SetNdivisions(505);

	gStyle->SetCanvasColor(0);
	gStyle->SetPadTopMargin(.05);
	gStyle->SetPadLeftMargin(.15);
	gStyle->SetPadRightMargin(.1);
	gStyle->SetPadBottomMargin(.15);
	gStyle->SetTitleYOffset(1.3);
	gStyle->SetLabelFont(32,"X");
	gStyle->SetLabelFont(32,"Y");

	// prepare gStyle to be useful
	//   1 = solid
	//   2 = long dash (30 10)
	gStyle->SetTitleFillColor(-1); 
	//   3 = dotted (4 8)
	//   4 = dash-dot (15 12 4 12)
	//   5 = short dash ( 15 15 )
	//   6 = dash-dot-dot 
	gStyle->SetLineStyleString(1,"[]");
	gStyle->SetLineStyleString(2,"[30 10]");
	gStyle->SetLineStyleString(3,"[4 8]");
	gStyle->SetLineStyleString(4,"[15 12 4 12]");
	gStyle->SetLineStyleString(5,"[15 15]");
	gStyle->SetLineStyleString(6,"[15 12 4 12 4 12]");
	gStyle->SetLabelSize(0.05,"X");
	gStyle->SetLabelSize(0.05,"Y");
	gStyle->SetNdivisions(505,"Y");
	gStyle->SetOptDate(0);
	gStyle->SetDateY(.98);
	gStyle->SetStripDecimals(kFALSE);

	gStyle->SetHistFillColor(0); 
	gStyle->SetFillColor(0); 
	gStyle->SetFillStyle(4000); 
	gStyle->SetStatStyle(0); 
	gStyle->SetTitleStyle(0); 
	gStyle->SetCanvasBorderSize(0); 
	gStyle->SetFrameBorderSize(0); 
	gStyle->SetLegendBorderSize(0); 
	gStyle->SetStatBorderSize(0); 
	gStyle->SetTitleBorderSize(0); 

	/*}}}*/

	cerr<<" 1 vs 2 (He3 or He4 or C12) for 1 = "; cin >> Target0;
	cerr<<" 1 vs 2 (H2  or He3 or He4) for 2 = "; cin >> Target;
	cerr<<" Com for 1 = "; cin >> Com0;
	cerr<<" Com for 2 = "; cin >> Com;

	int A1 = 0, A2 = 0, Z1 = 0, Z2=0;
	/*Get Target Info 1{{{*/
	if(Target0 == "H2"||Target == "D2"){
		A1 = 2; Z1 = 1;  
	}
	else if(Target0 == "He3"){
		A1 = 3; Z1 = 2;  
	}
	else if(Target0 == "He4"){
		A1 = 4; Z1 = 2; 
	}
	else if(Target0 == "Be"){
		A1 = 9; Z1 = 4; 
	}
	else if(Target0 == "C12"){
		A1 = 12; Z1 = 6; 
	}
	else if(Target0 == "Al"){
		A1 = 27; Z1 = 13; 
	}
	else if(Target0 == "Ca40"){
		A1 = 40; Z1 = 20;
	}
	else if(Target0 == "Ca48"){
		A1 = 48; Z1 = 20; 
	}
	else
		cerr<<" **** ERROR, I don't know the target!!!"<<endl;
	/*}}}*/
	/*Get Target Info 2{{{*/
	if(Target == "H2"||Target == "D2"){
		A2 = 2; Z2 = 1;  
	}
	else if(Target == "He3"){
		A2 = 3; Z2 = 2;  
	}
	else if(Target == "He4"){
		A2 = 4; Z2 = 2; 
	}
	else if(Target == "Be"){
		A2 = 9; Z2 = 4; 
	}
	else if(Target == "C12"){
		A2 = 12; Z2 = 6; 
	}
	else if(Target == "Al"){
		A2 = 27; Z2 = 13; 
	}
	else if(Target == "Ca40"){
		A2 = 40; Z2 = 20;
	}
	else if(Target == "Ca48"){
		A2 = 48; Z2 = 20; 
	}
	else
		cerr<<" **** ERROR, I don't know the target!!!"<<endl;
	/*}}}*/

	/*Get Ratio{{{*/
	//bool bMC = kTRUE;
	bool bMC = kFALSE;
	double Theta =21;
	double Zero[100] = {100*0.0};

	double xbj_21[100], R_EX_21[100], R_EX_Err_21[100], R_MC_21[100];
	int N_21 = 0;
	Merge(Form("./input/R21_%s_%s.rt",Target0.Data(),Target.Data()), xbj_21, R_EX_21, R_EX_Err_21, R_MC_21, &N_21);
	double Q2_21_D = 4.0*E0*Ep_Cut_21_1[0]*pow(TMath::Sin(21.*3.1415926/180.0*0.5),2);
	double Q2_21_U = 4.0*E0*Ep_Cut_21_2[3]*pow(TMath::Sin(21.*3.1415926/180.0*0.5),2);
	
	double xbj_23[100], R_EX_23[100], R_EX_Err_23[100], R_MC_23[100];
	int N_23 = 0;
	if(!(Target=="H2"))
		Merge(Form("./input/R23_%s_%s.rt",Target0.Data(),Target.Data()), xbj_23, R_EX_23, R_EX_Err_23, R_MC_23, &N_23);
	double Q2_23_D = 4.0*E0*Ep_Cut_23_1[0]*pow(TMath::Sin(23.*3.1415926/180.0*0.5),2);
	double Q2_23_U = 4.0*E0*Ep_Cut_23_2[3]*pow(TMath::Sin(23.*3.1415926/180.0*0.5),2);
	
	double xbj_QE[100], R_EX_QE[100], R_EX_Err_QE[100], R_MC_QE[100];
	int N_QE = 0;
	if(!(Target=="H2"))
		Merge(Form("./input/RQE_%s_%s.rt",Target0.Data(),Target.Data()), xbj_QE, R_EX_QE, R_EX_Err_QE, R_MC_QE, &N_QE);
	double Q2_25_D = 4.0*E0*Ep_Cut_QE_1[0]*pow(TMath::Sin(25.*3.1415926/180.0*0.5),2);
	
	double xbj_25[100], R_EX_25[100], R_EX_Err_25[100], R_MC_25[100];
	int N_25 = 0;
	Merge(Form("./input/R25_%s_%s.rt",Target0.Data(),Target.Data()), xbj_25, R_EX_25, R_EX_Err_25, R_MC_25, &N_25);
	double Q2_25_U = 4.0*E0*Ep_Cut_25_2[3]*pow(TMath::Sin(25.*3.1415926/180.0*0.5),2);

	double xbj_28[100], R_EX_28[100], R_EX_Err_28[100], R_MC_28[100];
	int N_28 = 0;
	if(Target0=="C12" && Target=="He3"){
		Merge(Form("./input/R28_%s_%s.rt",Target0.Data(),Target.Data()), xbj_28, R_EX_28, R_EX_Err_28, R_MC_28, &N_28);
	}
	double Q2_28_D = 4.0*E0*Ep_Cut_28_1[0]*pow(TMath::Sin(28.*3.1415926/180.0*0.5),2);
	double Q2_28_U = 4.0*E0*Ep_Cut_28_1[3]*pow(TMath::Sin(28.*3.1415926/180.0*0.5),2);
	
	//Add E02019 Data
	double xbj_NF[1000], R_EX_NF[1000], R_EX_Err_NF[1000],R_EX_Err_Up[1000],R_EX_Err_Down[1000];
	//Only fit one angle 18,22,26,32,40,50, 0->for all:
	//if(!(Target0=="C12"&&Target=="He4")){
	int N_NF = 0;  Double_t Theta_NF=18;

	if(Target=="He3" && !(Target0=="He4"))
		E02019_A_He3(Target0.Data(), Theta_NF,xbj_NF, R_EX_NF, R_EX_Err_Up,  R_EX_Err_Down, &N_NF);
	if(Target=="He3" && Target0=="He4")
		E02019_He4_He3(Target0.Data(), Theta_NF,xbj_NF, R_EX_NF, R_EX_Err_Up,  R_EX_Err_Down, &N_NF);
	
	if(0){
		Theta_NF = 18; //0->all
		E02019(Form("./input/RNF_%s_%s.rt",Target0.Data(),Target.Data()), Theta_NF,xbj_NF, R_EX_NF, R_EX_Err_NF, &N_NF);
		Theta_NF = 22; //0->all
		E02019(Form("./input/RNF_%s_%s.rt",Target0.Data(),Target.Data()), Theta_NF,xbj_NF, R_EX_NF, R_EX_Err_NF, &N_NF);
		Theta_NF = 26; //0->all
		E02019(Form("./input/RNF_%s_%s.rt",Target0.Data(),Target.Data()), Theta_NF,xbj_NF, R_EX_NF, R_EX_Err_NF, &N_NF);
	}
	cerr<<Form("@@@@@@@@@@ N_21 = %d, N_23 = %d, N_25 = %d, N_QE = %d, N_28 = %d, N_NF = %d", N_21,N_23,N_25,N_QE,N_28,N_NF)<<endl;

	//Add CLAS Data
	double xbj_CS[1000], R_EX_CS[1000], R_EX_Err_CS[1000];
	int N_CS = 0;  
	if(Target=="He3" && Target0=="He4")
		CLAS_He4_He3( xbj_CS, R_EX_CS, R_EX_Err_CS, &N_CS);
	cerr<<Form("@@@@@@@@@@ N_CLAS = %d", N_CS)<<endl;
/*}}}*/

	/*Plot {{{*/
	gStyle->SetOptStat(0);

	double RMAX = -1000,RMIN=0.4;
	for(int i=0; i<N_21; i++)
		RMAX = RMAX?R_EX_21[i]:RMAX<R_EX_21[i];
	for(int i=0; i<N_23; i++)
		RMAX = RMAX?R_EX_23[i]:RMAX<R_EX_23[i];
	for(int i=0; i<N_25; i++)
		RMAX = RMAX?R_EX_25[i]:RMAX<R_EX_25[i];
	if(Target=="He3"){
		for(int i=0; i<N_28; i++)
			RMAX = RMAX?R_EX_28[i]:RMAX<R_EX_28[i];
	}
	for(int i=0; i<N_21; i++)
		RMAX = RMAX?R_MC_21[i]:RMAX<R_MC_21[i];
	for(int i=0; i<N_23; i++)
		RMAX = RMAX?R_MC_23[i]:RMAX<R_MC_23[i];
	for(int i=0; i<N_25; i++)
		RMAX = RMAX?R_MC_25[i]:RMAX<R_MC_25[i];
	if(Target=="He3"){
		for(int i=0; i<N_28; i++)
			RMAX = RMAX?R_MC_28[i]:RMAX<R_MC_28[i];
	}

	double XMax = 0;
	if(Target=="He4" &&Target0=="C12")
		XMax=4.2;
	else
		XMax=3.2;
/*
	if(Target=="He3" &&Target0=="C12")
		RMAX = 20.0;
	if(Target=="He3" &&Target0=="He4")
		RMAX = 15.0;
	if(Target=="He4" &&Target0=="C12")
		RMAX = 7;
	if(Target=="Ca40" &&Target0=="Ca48")
		RMAX = 1.5;
*/

	TCanvas *c3 = new TCanvas("c3","c3",1200,700);
	c3->cd();
	TH2F *h3 = new TH2F("h3","",300,0.7,XMax,300,RMIN,RMAX);
	h3->SetStats(kFALSE);
	h3->SetXTitle("x_{bj}");
	h3->GetXaxis()->CenterTitle(1);
	h3->GetXaxis()->SetTitleFont(32);
	h3->GetXaxis()->SetTitleSize(0.06);
	h3->GetYaxis()->SetTitle(Form("#frac{Y_{%s}}{Y_{%s}} #times #frac{%d}{%d}",Target0.Data(),Target.Data(),A2,A1));
	h3->GetYaxis()->CenterTitle(1);
	h3->GetYaxis()->SetTitleFont(32);
	h3->GetYaxis()->SetTitleSize(0.04);
	h3->GetYaxis()->CenterTitle();
	gPad->SetLogy(0);
	h3->Draw();

	if(Target=="He3" &&Target0=="C12")
		RMAX = 10.0;
	if(Target=="He3" &&Target0=="He4")
		RMAX = 5.0;
	if(Target=="He4" &&Target0=="C12")
		RMAX = 4;
	if(Target=="Ca40" &&Target0=="Ca48"){
		RMAX = 1.3; RMIN=0.7;
	}

	TCanvas *c4 = new TCanvas("c4","c4",1200,700);
	c4->cd();
	TH2F *h4 = new TH2F("h4","",300,0.7,XMax,300,RMIN,RMAX);
	h4->SetStats(kFALSE);
	h4->SetXTitle("x_{bj}");
	h4->GetXaxis()->CenterTitle(1);
	h4->GetXaxis()->SetTitleFont(32);
	h4->GetXaxis()->SetTitleSize(0.06);
	h4->GetYaxis()->SetTitle(Form("#frac{Y_{%s}}{Y_{%s}} #times #frac{%d}{%d}",Target0.Data(),Target.Data(),A2,A1));
	h4->GetYaxis()->CenterTitle(1);
	h4->GetYaxis()->SetTitleFont(32);
	h4->GetYaxis()->SetTitleSize(0.04);
	h4->GetYaxis()->CenterTitle();
	gPad->SetLogy(0);
	h4->Draw();

	TLegend *l1 = new TLegend(0.25,0.60,0.50, 0.90,Form("Y_{%s} / Y_{%s} Ratio: ",Target0.Data(),Target.Data()));
	l1->SetTextSize(0.03);
	l1->SetTextFont(22);

	TGraphErrors *ex21 = new TGraphErrors(N_21,xbj_21, R_EX_21, Zero,R_EX_Err_21);
	ex21->SetMarkerStyle(24);
	ex21->SetMarkerColor(3);
	ex21->SetLineColor(3);
	Theta = 21.0;
	c3->cd();	ex21->Draw("P");
	c4->cd();	ex21->Draw("P");
	l1->AddEntry(ex21, Form("%2.0f^{o}: %4.3f<Q^{2}<%4.3f (GeV/c^{2})",Theta,Q2_21_D,Q2_21_U),"p");

	TGraphErrors *ex23 = new TGraphErrors(N_23,xbj_23, R_EX_23, Zero,R_EX_Err_23);
	ex23->SetMarkerStyle(25);
	ex23->SetMarkerColor(4);
	ex23->SetLineColor(4);
	Theta = 23.0;
	if(!(Target=="H2")){
		c3->cd();	ex23->Draw("P");
		c4->cd();	ex23->Draw("P");
		l1->AddEntry(ex23, Form("%2.0f^{o}: %4.3f<Q^{2}<%4.3f (GeV/c^{2})",Theta,Q2_23_D,Q2_23_U),"p");
	}

	TGraphErrors *exQE = new TGraphErrors(N_QE,xbj_QE, R_EX_QE, Zero,R_EX_Err_QE);
	exQE->SetMarkerStyle(26);
	exQE->SetMarkerColor(2);
	exQE->SetLineColor(2);
	c3->cd(); exQE->Draw("P");
	c4->cd(); exQE->Draw("P");

	TGraphErrors *ex25 = new TGraphErrors(N_25,xbj_25, R_EX_25, Zero,R_EX_Err_25);
	ex25->SetMarkerStyle(26);
	ex25->SetMarkerColor(2);
	ex25->SetLineColor(2);
	Theta = 25.0;
	c3->cd(); ex25->Draw("P");
	c4->cd(); ex25->Draw("P");
	l1->AddEntry(ex25, Form("%2.0f^{o}: %4.3f<Q^{2}<%4.3f (GeV/c^{2})",Theta,Q2_25_D,Q2_25_U),"p");

	TGraphErrors *ex28 = new TGraphErrors(N_28,xbj_28, R_EX_28, Zero,R_EX_Err_28);
	ex28->SetMarkerStyle(28);
	ex28->SetMarkerColor(6);
	ex28->SetLineColor(6);
	Theta = 28.0;
	if(Target0=="C12" && Target=="He3"){
		c3->cd(); 	ex28->Draw("P");
		c4->cd(); 	ex28->Draw("P");
		l1->AddEntry(ex28, Form("%2.0f^{o}: %4.3f<Q^{2}<%4.3f (GeV/c^{2})",Theta,Q2_28_D,Q2_28_U),"p");
	}


	if(bMC){	
		TGraphErrors *mc21 = new TGraphErrors(N_21,xbj_21, R_MC_21, Zero, R_MC_Err_21);
		//TGraphErrors *mc21 = new TGraphErrors(N_21,xbj_21, R_MC_21, Zero, Zero);
		mc21->SetLineWidth(2);
		mc21->SetLineColor(3);
		c3->cd(); mc21->Draw("l");
		c4->cd(); mc21->Draw("l");
		l1->AddEntry(mc21, Form("R_{MC} at %d^{0}",21),"l");

		TGraphErrors *mc23 = new TGraphErrors(N_23,xbj_23, R_MC_23, Zero, R_MC_Err_23);
		//TGraphErrors *mc23 = new TGraphErrors(N_23,xbj_23, R_MC_23, Zero, Zero);
		mc23->SetLineWidth(2);
		mc23->SetLineColor(4);
		c3->cd();	mc23->Draw("l");
		c4->cd();	mc23->Draw("l");
		l1->AddEntry(mc23, Form("R_{MC} at %d^{0}",23),"l");

		TGraphErrors *mcQE = new TGraphErrors(N_QE,xbj_QE, R_MC_QE, Zero, R_MC_Err_QE);
		//TGraphErrors *mcQE = new TGraphErrors(N_QE,xbj_QE, R_MC_QE, Zero, Zero);
		mcQE->SetLineWidth(2);
		mcQE->SetLineColor(2);
		c3->cd();	mcQE->Draw("l");
		c4->cd();	mcQE->Draw("l");

		TGraphErrors *mc25 = new TGraphErrors(N_25,xbj_25, R_MC_25, Zero, R_MC_Err_25);
		//TGraphErrors *mc25 = new TGraphErrors(N_25,xbj_25, R_MC_25, Zero, Zero);
		mc25->SetLineWidth(2);
		mc25->SetLineColor(2);
		c3->cd();	mc25->Draw("l");
		c4->cd();	mc25->Draw("l");
		l1->AddEntry(mc25, Form("R_{MC} at %d^{0}",25),"l");

		TGraphErrors *mc28 = new TGraphErrors(N_28,xbj_28, R_MC_28, Zero, R_MC_Err_28);
		//TGraphErrors *mc28 = new TGraphErrors(N_28,xbj_28, R_MC_28, Zero, Zero);
		mc28->SetLineWidth(2);
		mc28->SetLineColor(6);
		if(Target=="C12"&&Target=="He3"){
			c3->cd();	mc28->Draw("l");
			c4->cd();	mc28->Draw("l");
			l1->AddEntry(mc28, Form("R_{MC} at %d^{0}",28),"l");
		}
	}

	c3->cd(); l1->Draw();
	c4->cd(); l1->Draw();
	if(bMC){
		c3->Print(Form("%s_%s_Yield_Ratio_MC.png",Target0.Data(),Target.Data()));
		c3->Print(Form("%s_%s_Yield_Ratio_MC.pdf",Target0.Data(),Target.Data()));
	}
	else{
		c3->Print(Form("%s_%s_Yield_Ratio.png",Target0.Data(),Target.Data()));
		c3->Print(Form("%s_%s_Yield_Ratio.pdf",Target0.Data(),Target.Data()));
	}
	if(bMC){
		c4->Print(Form("%s_%s_Yield_Ratio_MC_Zoom.png",Target0.Data(),Target.Data()));
		c4->Print(Form("%s_%s_Yield_Ratio_MC_Zoom.pdf",Target0.Data(),Target.Data()));
	}
	else{
		c4->Print(Form("%s_%s_Yield_Ratio_Zoom.png",Target0.Data(),Target.Data()));
		c4->Print(Form("%s_%s_Yield_Ratio_Zoom.pdf",Target0.Data(),Target.Data()));
	}
/*}}}*/


	return 0;
}
/*}}}*/

/*gGet_Ep_Cut{{{*/
void gGet_Ep_Cut(const TString &aKinName, const TString & aArm, double* aEp_Min, double* aEp_Max){

	if(aKinName.Contains("Kin3.1")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_21_1[0];		
			aEp_Max[0] = Ep_Cut_21_1[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_21_1[2];		
			aEp_Max[0] = Ep_Cut_21_1[3];		
		}
	}
	else if(aKinName.Contains("Kin3.2")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_21_2[0];		
			aEp_Max[0] = Ep_Cut_21_2[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_21_2[2];		
			aEp_Max[0] = Ep_Cut_21_2[3];		
		}
	}

	else if(aKinName.Contains("Kin4.1")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_23_1[0];		
			aEp_Max[0] = Ep_Cut_23_1[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_23_1[2];		
			aEp_Max[0] = Ep_Cut_23_1[3];		
		}
	}
	else if(aKinName.Contains("Kin4.2")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_23_2[0];		
			aEp_Max[0] = Ep_Cut_23_2[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_23_2[2];		
			aEp_Max[0] = Ep_Cut_23_2[3];		
		}
	}

	else if(aKinName.Contains("Kin5.0") && ! aKinName.Contains("Kin5.05")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_QE_1[0];		
			aEp_Max[0] = Ep_Cut_QE_1[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_QE_1[2];		
			aEp_Max[0] = Ep_Cut_QE_1[3];		
		}
	}
	else if(aKinName.Contains("Kin5.05")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_QE_2[0];		
			aEp_Max[0] = Ep_Cut_QE_2[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_QE_2[2];		
			aEp_Max[0] = Ep_Cut_QE_2[3];		
		}
	}

	else if(aKinName.Contains("Kin5.1")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_25_1[0];		
			aEp_Max[0] = Ep_Cut_25_1[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_25_1[2];		
			aEp_Max[0] = Ep_Cut_25_1[3];		
		}
	}
	else if(aKinName.Contains("Kin5.2")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_25_2[0];		
			aEp_Max[0] = Ep_Cut_25_2[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_25_2[2];		
			aEp_Max[0] = Ep_Cut_25_2[3];		
		}
	}

	else if(aKinName.Contains("Kin6.5")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_28_1[0];		
			aEp_Max[0] = Ep_Cut_28_1[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_28_1[2];		
			aEp_Max[0] = Ep_Cut_28_1[3];		
		}
	}

	else{
		aEp_Min[0]=-1.0;
		aEp_Max[0]=10.0;
	}

}
/*}}}*/

/*Merge{{{*/
void Merge(const TString& aInputFileName, double* aXbj, double* aR_EX, double* aR_EX_Err, double* aR_MC, double* aR_MC_Err, int* aN){

	/*Pre-define{{{*/	
	int aline = 0;
	TString aInputFilePath = "/work/halla/e08014/disk1/yez/Ep_Cross_Section/results/";

	double scale = 1.0;
	const double He4_Norm = 1.0/0.864;

	//First Target:
	const TString aName_T1L1 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	const TString aName_T1L2 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	const TString aName_T1R1 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	const TString aName_T1R2 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	//Second Target:
	const TString aName_T2L1 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	const TString aName_T2L2 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	const TString aName_T2R1 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	const TString aName_T2R2 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	const double factor = gGet_InputFile_Var(aInputFileName, aline).Atof(); 

	double Theta = 21.00;
	TString Target = "C12";
	TString Kin = "3.1";
	double Ep_Min = 1000.0, Ep_Max = -1000.0;
	Double_t nu,Q2,xbj,angle,Ep,Livetime,Ntg,Ni,Nf_EX,Nf_SAMC,Nf_SAMC_Gen,BinCor;
	Double_t XS,XS_Stat_Err,XS_Sys_Err,XS_MC,Yield_MC,Yield_EX,Y_MC_Sys_Err,Y_EX_Sys_Err,Y_MC_Stat_Err,Y_EX_Stat_Err;
	Double_t Nf_AL, Yield_AL, Y_AL_Stat_Err, Y_AL_Sys_Err;
	TString comment;
	/*}}}*/

	/*First Target{{{*/
	cerr<<endl<<"------------------- First Target ----------------------------------------------------------"<<endl<<endl;

	/*Define{{{*/
	Double_t XS_T1L1[200], XS_Err_T1L1[200], Ep_T1L1[200], Xbj_T1L1[200],  Xbj_Err_T1L1[200],XS_MC_T1L1[200];
	Double_t XS_T1L2[200], XS_Err_T1L2[200], Ep_T1L2[200], Xbj_T1L2[200],  Xbj_Err_T1L2[200],XS_MC_T1L2[200];
	Double_t XS_T1R1[200], XS_Err_T1R1[200], Ep_T1R1[200], Xbj_T1R1[200],  Xbj_Err_T1R1[200],XS_MC_T1R1[200];
	Double_t XS_T1R2[200], XS_Err_T1R2[200], Ep_T1R2[200], Xbj_T1R2[200],  Xbj_Err_T1R2[200],XS_MC_T1R2[200];
	Double_t Y_EX_T1L1[200],  Y_EX_Err_T1L1[200],  Y_MC_T1L1[200],  Y_MC_Err_T1L1[200];
	Double_t Y_EX_T1L2[200],  Y_EX_Err_T1L2[200],  Y_MC_T1L2[200],  Y_MC_Err_T1L2[200];
	Double_t Y_EX_T1R1[200],  Y_EX_Err_T1R1[200],  Y_MC_T1R1[200],  Y_MC_Err_T1R1[200];
	Double_t Y_EX_T1R2[200],  Y_EX_Err_T1R2[200],  Y_MC_T1R2[200],  Y_MC_Err_T1R2[200];
		/*}}}*/

	/*First{{{*/
	TString infile_T1L1 = aInputFilePath + aName_T1L1 +"_"+ Com0;
	cerr<<"==== Reading 1st file from "<<infile_T1L1<<endl;
	ifstream run_T1L1; run_T1L1.open(infile_T1L1);

	Int_t Count_T1L1 =0;
	Target = gGet_Target_Name(infile_T1L1.Data());
	if(Target=="He4")
		scale = He4_Norm;
	else
		scale = 1.0;

	Kin = gGet_Kin_Name(infile_T1L1.Data());
	TString Kin0 = Kin;
	TString Arm_T1L1="X";
	if(aName_T1L1.Contains("L_"))
		Arm_T1L1="L";
	else if(aName_T1L1.Contains("R_"))
		Arm_T1L1="R";


	if(Kin=="MM"){
		cerr<<"1st: Kin (3.1,3.2,4.1,4.2,5.0,5.05,5.1,5.2,6.5) = "; cin >> Kin;
	}
	if(Kin=="3.1"||Kin=="3.2")
		Theta = 21.00; 
	else if(Kin=="4.1"||Kin=="4.2")
		Theta = 23.00;
	else if(Kin=="5.0"||Kin=="5.05"||Kin=="5.1"||Kin=="5.2")
		Theta =  25.00;
	else if(Kin=="6.5")
		Theta = 28.00;
	else
		cerr<<"*** ERROR, unknown kinematics setting!"<<endl;

	TString KinName_T1L1 = Form("Kin%s",Kin.Data());
	gGet_Ep_Cut(KinName_T1L1.Data(), Arm_T1L1.Data(),&Ep_Min, &Ep_Max);
	cerr<<Form("@@@  For %s, Ep_Min = %f, Ep_Max = %f", KinName_T1L1.Data(), Ep_Min,Ep_Max)<<endl;

	Get_XS* xs_T1L1 = new Get_XS();
	xs_T1L1->Load_Table(Target.Data(),KinName_T1L1.Data());

	comment.ReadLine(run_T1L1);
	while ( run_T1L1>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
	{
		if ( Ep>=Ep_Min && Ep<=Ep_Max && angle<180 &&xbj>0.5  && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
		{
			Q2 = 4.0 * E0 * Ep * pow(sin(Theta*3.1415926/180./2.0),2);
			xbj = Q2/2.0/0.938272/(E0-Ep);
			Ep_T1L1[Count_T1L1] = Ep; 
			XS_T1L1[Count_T1L1] = XS * scale; 
			XS_Err_T1L1[Count_T1L1] = sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2));

			Y_EX_T1L1[Count_T1L1] = Yield_EX; 
			Y_MC_T1L1[Count_T1L1] = Yield_MC; 
			Y_EX_Err_T1L1[Count_T1L1] = sqrt(pow(Y_EX_Sys_Err,2)+pow(Y_EX_Stat_Err,2));
			Y_MC_Err_T1L1[Count_T1L1] = sqrt(pow(Y_MC_Sys_Err,2)+pow(Y_MC_Stat_Err,2));
				//if(Target.Contains("Ca"))
			if(1)
				XS_Err_T1L1[Count_T1L1] = pow(XS_Stat_Err,1);
			Xbj_T1L1[Count_T1L1] = xbj;    Xbj_Err_T1L1[Count_T1L1] = 0.0;
			XS_MC_T1L1[Count_T1L1] = xs_T1L1->gGet_XS_Born(E0,Ep,Theta);
			Count_T1L1++;
		}
		if ( !run_T1L1.good() )
			break;
	}
	run_T1L1.close();
	delete xs_T1L1;
	cerr << "  Total XS Points for 1st = " << Count_T1L1 <<endl;
	/*}}}*/

	/*Second{{{*/
	TString KinName_T1L2 ="";
	Bool_t bMore_T1L2 = kTRUE;
	TString infile_T1L2 = aInputFilePath + aName_T1L2 +"_"+ Com0;
	cerr<<"==== Reading 2nd file from "<<infile_T1L2<<endl;
	TString Arm_T1L2="X";
	if(aName_T1L2 =="0")
		bMore_T1L2 = kFALSE;
	if(aName_T1L2.Contains("L_"))
		Arm_T1L2="L";
	else if(aName_T1L2.Contains("R_"))
		Arm_T1L2="R";

	Int_t Count_T1L2=0;
	if(bMore_T1L2){
		ifstream run_T1L2; run_T1L2.open(infile_T1L2);
		Target = gGet_Target_Name(infile_T1L2.Data());
		if(Target=="He4")
		scale = He4_Norm;
		else
			scale = 1.0;

	Kin = gGet_Kin_Name(infile_T1L2.Data());
		if(Kin=="MM"){
			cerr<<"2nd: Kin (3.1,3.2,4.1,4.2,5.0,5.05,5.1,5.2,6.5) = "; cin >> Kin;
		}
		if(Kin=="3.1"||Kin=="3.2")
			Theta = 21.00;
		else if(Kin=="4.1"||Kin=="4.2")
			Theta = 23.00;
		else if(Kin=="5.0"||Kin=="5.05"||Kin=="5.1"||Kin=="5.2")
			Theta =  25.00;
		else if(Kin=="6.5")
			Theta = 28.00;
		else
			cerr<<"*** ERROR, unknown kinematics setting!"<<endl;
		KinName_T1L2 = Form("Kin%s",Kin.Data());
		gGet_Ep_Cut(KinName_T1L2.Data(), Arm_T1L2.Data(),&Ep_Min, &Ep_Max);
		cerr<<Form("@@@  For %s, Ep_Min = %f, Ep_Max = %f", KinName_T1L2.Data(), Ep_Min,Ep_Max)<<endl;

		Get_XS* xs_T1L2 = new Get_XS();
		xs_T1L2->Load_Table(Target.Data(),KinName_T1L2.Data());

		comment.ReadLine(run_T1L2);
	while ( run_T1L2>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
		{
			if ( Ep>=Ep_Min && Ep<=Ep_Max && angle<180 &&xbj>0.5  && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
			{
				Q2 = 4.0 * E0 * Ep * pow(sin(Theta*3.1415926/180./2.0),2);
				xbj = Q2/2.0/0.938272/(E0-Ep);
				Ep_T1L2[Count_T1L2] = Ep; 
				XS_T1L2[Count_T1L2] = XS * scale;    
				XS_Err_T1L2[Count_T1L2] = sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2));

				Y_EX_T1L2[Count_T1L2] = Yield_EX; 
				Y_MC_T1L2[Count_T1L2] = Yield_MC; 
				Y_EX_Err_T1L2[Count_T1L2] = sqrt(pow(Y_EX_Sys_Err,2)+pow(Y_EX_Stat_Err,2));
				Y_MC_Err_T1L2[Count_T1L2] = sqrt(pow(Y_MC_Sys_Err,2)+pow(Y_MC_Stat_Err,2));
				//	if(Target.Contains("Ca"))
				if(1)
						XS_Err_T1L2[Count_T1L2] = pow(XS_Stat_Err,2);
				Xbj_T1L2[Count_T1L2] = xbj;    Xbj_Err_T1L2[Count_T1L2] = 0.0;
				XS_MC_T1L2[Count_T1L2] = xs_T1L2->gGet_XS_Born(E0,Ep,Theta);
				Count_T1L2++;
			}
			if ( !run_T1L2.good() )
				break;
		}
		run_T1L2.close();
		delete xs_T1L2;
		cerr << "  Total XS Points for 2nd  = " << Count_T1L2 <<endl;
	}
	/*}}}*/

	/*Third{{{*/
	TString KinName_T1R1 ="";
	Bool_t bMore_T1R1 = kTRUE;
	TString infile_T1R1 =  aInputFilePath +aName_T1R1 +"_"+ Com0;
	cerr<<"==== Reading 3rd file from "<<infile_T1R1<<endl;
	TString Arm_T1R1="X";
	if(aName_T1R1 =="0")
		bMore_T1R1 = kFALSE;
	if(aName_T1R1.Contains("L_"))
		Arm_T1R1="L";
	else if(aName_T1R1.Contains("R_"))
		Arm_T1R1="R";

	Int_t Count_T1R1=0;
	bMore_T1R1 = kFALSE;
	if(bMore_T1R1){
		ifstream run_T1R1; run_T1R1.open(infile_T1R1);
		Target = gGet_Target_Name(infile_T1R1.Data());
		if(Target=="He4")
		scale = He4_Norm;
		else
			scale = 1.0;

	Kin = gGet_Kin_Name(infile_T1R1.Data());
		if(Kin=="MM"){
			cerr<<"2nd: Kin (3.1,3.2,4.1,4.2,5.0,5.05,5.1,5.2,6.5) = "; cin >> Kin;
		}
		if(Kin=="3.1"||Kin=="3.2")
			Theta = 21.00;
		else if(Kin=="4.1"||Kin=="4.2")
			Theta = 23.00;
		else if(Kin=="5.0"||Kin=="5.05"||Kin=="5.1"||Kin=="5.2")
			Theta =  25.00;
		else if(Kin=="6.5")
			Theta = 28.00;
		else
			cerr<<"*** ERROR, unknown kinematics setting!"<<endl;
		KinName_T1R1 = Form("Kin%s",Kin.Data());
		gGet_Ep_Cut(KinName_T1R1.Data(), Arm_T1R1.Data(),&Ep_Min, &Ep_Max);
		cerr<<Form("@@@  For %s, Ep_Min = %f, Ep_Max = %f", KinName_T1R1.Data(), Ep_Min,Ep_Max)<<endl;

		Get_XS* xs_T1R1 = new Get_XS();
		xs_T1R1->Load_Table(Target.Data(),KinName_T1R1.Data());

		comment.ReadLine(run_T1R1);
	while ( run_T1R1>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
		{
			if ( Ep>=Ep_Min && Ep<=Ep_Max && angle<180 &&xbj>0.5  && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
			{
				Q2 = 4.0 * E0 * Ep * pow(sin(Theta*3.1415926/180./2.0),2);
				xbj = Q2/2.0/0.938272/(E0-Ep);
				Ep_T1R1[Count_T1R1] = Ep; 
				XS_T1R1[Count_T1R1] = XS * scale; 
				XS_Err_T1R1[Count_T1R1] = sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2));

				Y_EX_T1R1[Count_T1R1] = Yield_EX; 
				Y_MC_T1R1[Count_T1R1] = Yield_MC; 
				Y_EX_Err_T1R1[Count_T1R1] = sqrt(pow(Y_EX_Sys_Err,2)+pow(Y_EX_Stat_Err,2));
				Y_MC_Err_T1R1[Count_T1R1] = sqrt(pow(Y_MC_Sys_Err,2)+pow(Y_MC_Stat_Err,2));
				//	if(Target.Contains("Ca"))
				if(1)
					XS_Err_T1R1[Count_T1R1] = pow(XS_Stat_Err,1);
				Xbj_T1R1[Count_T1R1] = xbj;    Xbj_Err_T1R1[Count_T1R1] = 0.0;
				XS_MC_T1R1[Count_T1R1] = xs_T1R1->gGet_XS_Born(E0,Ep,Theta);
				Count_T1R1++;
			}
			if ( !run_T1R1.good() )
				break;
		}
		run_T1R1.close();
		delete xs_T1R1;
		cerr << "  Total XS Points for 2nd  = " << Count_T1R1 <<endl;
	}
	/*}}}*/

	/*Forth{{{*/
	TString KinName_T1R2 ="";
	Bool_t bMore_T1R2 = kTRUE;
	TString infile_T1R2 = aInputFilePath + aName_T1R2 +"_"+ Com0;
	cerr<<"==== Reading 4th file from "<<infile_T1R2<<endl;
	TString Arm_T1R2="X";
	if(aName_T1R2 =="0")
		bMore_T1R2 = kFALSE;
	if(aName_T1R2.Contains("L_"))
		Arm_T1R2="L";
	else if(aName_T1R2.Contains("R_"))
		Arm_T1R2="R";

	Int_t Count_T1R2=0;
	bMore_T1R2 = kFALSE;
	if(bMore_T1R2){
		ifstream run_T1R2; run_T1R2.open(infile_T1R2);
		Target = gGet_Target_Name(infile_T1R2.Data());
		if(Target=="He4")
		scale = He4_Norm;
		else
			scale = 1.0;

	Kin = gGet_Kin_Name(infile_T1R2.Data());
		if(Kin=="MM"){
			cerr<<"2nd: Kin (3.1,3.2,4.1,4.2,5.0,5.05,5.1,5.2,6.5) = "; cin >> Kin;
		}
		if(Kin=="3.1"||Kin=="3.2")
			Theta = 21.00;
		else if(Kin=="4.1"||Kin=="4.2")
			Theta = 23.00;
		else if(Kin=="5.0"||Kin=="5.05"||Kin=="5.1"||Kin=="5.2")
			Theta =  25.00;
		else if(Kin=="6.5")
			Theta = 28.00;
		else
			cerr<<"*** ERROR, unknown kinematics setting!"<<endl;
		KinName_T1R2 = Form("Kin%s",Kin.Data());
		gGet_Ep_Cut(KinName_T1R2.Data(), Arm_T1R2.Data(),&Ep_Min, &Ep_Max);
		cerr<<Form("@@@  For %s, Ep_Min = %f, Ep_Max = %f", KinName_T1R2.Data(), Ep_Min,Ep_Max)<<endl;

		Get_XS* xs_T1R2 = new Get_XS();
		xs_T1R2->Load_Table(Target.Data(),KinName_T1R2.Data());

		comment.ReadLine(run_T1R2);
	while ( run_T1R2>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
		{
			if ( Ep>=Ep_Min && Ep<=Ep_Max && angle<180 &&xbj>0.5  && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
			{
				Q2 = 4.0 * E0 * Ep * pow(sin(Theta*3.1415926/180./2.0),2);
				xbj = Q2/2.0/0.938272/(E0-Ep);
				Ep_T1R2[Count_T1R2] = Ep; 
				XS_T1R2[Count_T1R2] = XS * scale;
				XS_Err_T1R2[Count_T1R2] = sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2));

				Y_EX_T1R2[Count_T1R2] = Yield_EX; 
				Y_MC_T1R2[Count_T1R2] = Yield_MC; 
				Y_EX_Err_T1R2[Count_T1R2] = sqrt(pow(Y_EX_Sys_Err,2)+pow(Y_EX_Stat_Err,2));
				Y_MC_Err_T1R2[Count_T1R2] = sqrt(pow(Y_MC_Sys_Err,2)+pow(Y_MC_Stat_Err,2));
				//	if(Target.Contains("Ca"))
				if(1)
					XS_Err_T1R2[Count_T1R2] = pow(XS_Stat_Err,1);
				Xbj_T1R2[Count_T1R2] = xbj;    Xbj_Err_T1R2[Count_T1R2] = 0.0;
				XS_MC_T1R2[Count_T1R2] = xs_T1R2->gGet_XS_Born(E0,Ep,Theta);
				Count_T1R2++;
			}
			if ( !run_T1R2.good() )
				break;
		}
		run_T1R2.close();
		delete xs_T1R2;
		cerr << "  Total XS Points for 2nd  = " << Count_T1R2 <<endl;
	}
	/*}}}*/

	/*Merge{{{*/
	//Add First Two Input
	Double_t Ep_T1L[200],Xbj_T1L[200], XS_T1L[200], XS_Err_T1L[200], XS_MC_T1L[200];
	Double_t Y_EX_T1L[200],Y_EX_Err_T1L[200],Y_MC_T1L[200],Y_MC_Err_T1L[200];
	Double_t XS_Sum_T1L = 0.0, XS_Err_Sum_T1L = 0.0; Int_t fSum_T1L = 0;
	Double_t Y_EX_Sum_T1L = 0.0, Y_EX_Err_Sum_T1L = 0.0, Y_MC_Sum_T1L = 0.0, Y_MC_Err_Sum_T1L = 0.0;
	Int_t Count_T1L = 0, J_T1L = 0;
	for(int i=0;i<Count_T1L1;i++){
		XS_Sum_T1L = XS_T1L1[i];
		XS_Err_Sum_T1L = pow(XS_Err_T1L1[i]/XS_T1L1[i],2);

		Y_EX_Sum_T1L = Y_EX_T1L1[i];
		Y_MC_Sum_T1L = Y_MC_T1L1[i];
		Y_EX_Err_Sum_T1L = pow(Y_EX_Err_T1L1[i]/Y_EX_T1L1[i],2);
		Y_MC_Err_Sum_T1L = pow(Y_MC_Err_T1L1[i]/Y_MC_T1L1[i],2);
		fSum_T1L = 1;
		if(bMore_T1L2){
			for(int j=0;j<Count_T1L2;j++){
				//Find the verlapping bin for two setting, such as Kin3.1 and Kin3.2
				//if(abs(Xbj_T1L1[i]-Xbj_T1L2[j])< 0.0001){
				if(abs(Ep_T1L1[i]-Ep_T1L2[j])< 0.0001){
					XS_Sum_T1L += XS_T1L2[j];
					XS_Err_Sum_T1L +=	pow(XS_Err_T1L2[j]/XS_T1L2[j],2);

					Y_EX_Sum_T1L += Y_EX_T1L2[j];
					Y_MC_Sum_T1L += Y_MC_T1L2[j];
					Y_EX_Err_Sum_T1L += pow(Y_EX_Err_T1L2[j]/Y_EX_T1L2[j],2);
					Y_MC_Err_Sum_T1L += pow(Y_MC_Err_T1L2[j]/Y_MC_T1L2[j],2);
					fSum_T1L++; J_T1L = j;
					break; //No two xbj bin within one setting
				}
			}
		}
		Ep_T1L[i] = Ep_T1L1[i];
		Xbj_T1L[i] = Xbj_T1L1[i];
		XS_MC_T1L[i] = XS_MC_T1L1[i];
		XS_T1L[i] = XS_Sum_T1L / fSum_T1L;
		XS_Err_T1L[i] = XS_T1L[i] * sqrt(XS_Err_Sum_T1L); 
	
		Y_EX_T1L[i] = Y_EX_Sum_T1L / fSum_T1L;
		Y_MC_T1L[i] = Y_MC_Sum_T1L / fSum_T1L;
		Y_EX_Err_T1L[i] = Y_EX_T1L[i] * sqrt(Y_EX_Err_Sum_T1L); 
		Y_MC_Err_T1L[i] = Y_MC_T1L[i] * sqrt(Y_MC_Err_Sum_T1L); 
			
		Count_T1L++;
	}
	for(int k=J_T1L+1;k<Count_T1L2;k++){
		Ep_T1L[Count_T1L] = Ep_T1L2[k];
		Xbj_T1L[Count_T1L] = Xbj_T1L2[k];
		XS_T1L[Count_T1L] = XS_T1L2[k];
		XS_Err_T1L[Count_T1L] = XS_Err_T1L2[k];
		XS_MC_T1L[Count_T1L] = XS_MC_T1L2[k];

		Y_EX_T1L[Count_T1L] = Y_EX_T1L2[k];
		Y_MC_T1L[Count_T1L] = Y_MC_T1L2[k];
		Y_EX_Err_T1L[Count_T1L] = Y_EX_Err_T1L2[k];
		Y_MC_Err_T1L[Count_T1L] = Y_MC_Err_T1L2[k];

		Count_T1L++;	
	}

	//Add Last Two Input
	Double_t Ep_T1R[200], Xbj_T1R[200], XS_T1R[200], XS_Err_T1R[200], XS_MC_T1R[200];
	Double_t XS_Sum_T1R = 0.0, XS_Err_Sum_T1R = 0.0; Int_t fSum_T1R = 0;
	Double_t Y_EX_T1R[200],Y_EX_Err_T1R[200],Y_MC_T1R[200],Y_MC_Err_T1R[200];
	Double_t Y_EX_Sum_T1R = 0.0, Y_EX_Err_Sum_T1R = 0.0, Y_MC_Sum_T1R = 0.0, Y_MC_Err_Sum_T1R = 0.0;
	Int_t Count_T1R = 0, J_T1R = 0;
	if(bMore_T1R1){
		for(int i=0;i<Count_T1R1;i++){
			XS_Sum_T1R = XS_T1R1[i];
			XS_Err_Sum_T1R = pow(XS_Err_T1R1[i]/XS_T1R1[i],2);
		
			Y_EX_Sum_T1R = Y_EX_T1R1[i];
			Y_MC_Sum_T1R = Y_MC_T1R1[i];
			Y_EX_Err_Sum_T1R = pow(Y_EX_Err_T1R1[i]/Y_EX_T1R1[i],2);
			Y_MC_Err_Sum_T1R = pow(Y_MC_Err_T1R1[i]/Y_MC_T1R1[i],2);
			fSum_T1R = 1;
			if(bMore_T1R2){
				for(int j=0;j<Count_T1R2;j++){
					//Find the verlapping bin for two setting, such as Kin3.1 and Kin3.2
					//if(abs(Xbj_T1R1[i]-Xbj_T1R2[j])< 0.0001){
					if(abs(Ep_T1R1[i]-Ep_T1R2[j])< 0.0001){
						XS_Sum_T1R += XS_T1R2[j];
						XS_Err_Sum_T1R +=	pow(XS_Err_T1R2[j]/XS_T1R2[j],2);

						Y_EX_Sum_T1R += Y_EX_T1R2[j];
						Y_MC_Sum_T1R += Y_MC_T1R2[j];
						Y_EX_Err_Sum_T1R += pow(Y_EX_Err_T1R2[j]/Y_EX_T1R2[j],2);
						Y_MC_Err_Sum_T1R += pow(Y_MC_Err_T1R2[j]/Y_MC_T1R2[j],2);
						fSum_T1R++; J_T1R = j;
						break; //No two xbj bin within one setting
					}
				}
			}
			Ep_T1R[i] = Ep_T1R1[i];
			Xbj_T1R[i] = Xbj_T1R1[i];
			XS_MC_T1R[i] = XS_MC_T1R1[i];
			XS_T1R[i] = XS_Sum_T1R / fSum_T1R;
			XS_Err_T1R[i] = XS_T1R[i] * sqrt(XS_Err_Sum_T1R); 

			Y_EX_T1R[i] = Y_EX_Sum_T1R / fSum_T1R;
			Y_MC_T1R[i] = Y_MC_Sum_T1R / fSum_T1R;
			Y_EX_Err_T1R[i] = Y_EX_T1R[i] * sqrt(Y_EX_Err_Sum_T1R); 
			Y_MC_Err_T1R[i] = Y_MC_T1R[i] * sqrt(Y_MC_Err_Sum_T1R); 

			Count_T1R++;
		}
		for(int k=J_T1R+1;k<Count_T1R2;k++){
			Ep_T1R[Count_T1R] = Ep_T1R2[k];
			Xbj_T1R[Count_T1R] = Xbj_T1R2[k];
			XS_T1R[Count_T1R] = XS_T1R2[k];
			XS_MC_T1R[Count_T1R] = XS_MC_T1R2[k];
			XS_Err_T1R[Count_T1R] = XS_Err_T1R2[k];

			Y_EX_T1R[Count_T1R] = Y_EX_T1R2[k];
			Y_MC_T1R[Count_T1R] = Y_MC_T1R2[k];
			Y_EX_Err_T1R[Count_T1R] = Y_EX_Err_T1R2[k];
			Y_MC_Err_T1R[Count_T1R] = Y_MC_Err_T1R2[k];

			Count_T1R++;	
		}
	}

	vector<double> vEp1, vXbj1, vXbj_Err1, vXS1, vXS_Err1, vXS_MC1;
	vector<double> vY_EX1, vY_EX_Err1, vY_MC1, vY_MC_Err1;
	double XS_MC_T1, XS_T1,XS_Err_T1;
	double Y_EX_T1, Y_EX_Err_T1, Y_MC_T1, Y_MC_Err_T1;
	double X_Min1 = 100.0, X_Max1 = -100.0, Y_Min1 = 10000.0, Y_Max1 = -10000.0;
	int N1 = 0, J1R = 0;
	bool bFound_T1L = kFALSE;
	if(bMore_T1R1){
		for(int i=0;i<Count_T1L; i++){
			bFound_T1L = kFALSE;
			for(int j=0;j<Count_T1R; j++){
				//if( abs(Xbj_T1L[i] - Xbj_T1R[j])<0.0001 ){
				if( abs(Ep_T1L[i] - Ep_T1R[j])<0.0001 ){
					bFound_T1L = kTRUE;
					vEp1.push_back(Ep_T1L[i]);	
					vXbj1.push_back(Xbj_T1L[i]);	
					vXbj_Err1.push_back(0.0);	
					XS_T1 = (XS_T1L[i]+XS_T1R[j])/2.0;
					XS_Err_T1 = XS_T1 * sqrt( pow(XS_Err_T1L[i]/XS_T1L[i],2) + pow(XS_Err_T1R[j]/XS_T1R[j],2) );
					XS_MC_T1 = (XS_MC_T1L[i]+XS_MC_T1R[j])/2.0;
					vXS1.push_back(XS_T1);
					vXS_Err1.push_back(XS_Err_T1);
					vXS_MC1.push_back( XS_MC_T1 );

					Y_EX_T1 = (Y_EX_T1L[i]+Y_EX_T1R[j])/2.0;
					Y_MC_T1 = (Y_MC_T1L[i]+Y_MC_T1R[j])/2.0;
					Y_EX_Err_T1 = Y_EX_T1 * sqrt( pow(Y_EX_Err_T1L[i]/Y_EX_T1L[i],2) + pow(Y_EX_Err_T1R[j]/Y_EX_T1R[j],2) );
					Y_MC_Err_T1 = Y_MC_T1 * sqrt( pow(Y_MC_Err_T1L[i]/Y_MC_T1L[i],2) + pow(Y_MC_Err_T1R[j]/Y_MC_T1R[j],2) );
                    vY_EX1.push_back(Y_EX_T1);
                    vY_MC1.push_back(Y_MC_T1);
                    vY_EX_Err1.push_back(Y_EX_Err_T1);
                    vY_MC_Err1.push_back(Y_MC_Err_T1);

					N1++; J1R = j;

					if(Xbj_T1L[i] < X_Min1)
						X_Min1 = Xbj_T1L[i];
					if(Xbj_T1L[i] > X_Max1)
						X_Max1 = Xbj_T1L[i];
					if(XS_T1 < Y_Min1)
						Y_Min1 = XS_T1;
					if(XS_T1 > Y_Max1)
						Y_Max1 = XS_T1;
					if(XS_MC_T1 < Y_Min1)
						Y_Min1 = XS_MC_T1;
					if(XS_T1 > Y_Max1)
						Y_Max1 = XS_MC_T1;
				}
			}
			if(!bFound_T1L){
				vEp1.push_back(Ep_T1L[i]);	
				vXbj1.push_back(Xbj_T1L[i]);	
				vXbj_Err1.push_back(0.0);	
				XS_T1 = XS_T1L[i];
				XS_Err_T1 = XS_Err_T1L[i];
				XS_MC_T1 = XS_MC_T1L[i];
				vXS1.push_back(XS_T1);
				vXS_Err1.push_back(XS_Err_T1);
				vXS_MC1.push_back( XS_MC_T1 );

				Y_EX_T1 = Y_EX_T1L[i];
				Y_MC_T1 = Y_MC_T1L[i];
				Y_EX_Err_T1 = Y_EX_Err_T1L[i];
				Y_MC_Err_T1 = Y_MC_Err_T1L[i];
				vY_EX1.push_back(Y_EX_T1);
				vY_MC1.push_back(Y_MC_T1);
				vY_EX_Err1.push_back(Y_EX_Err_T1);
				vY_MC_Err1.push_back(Y_MC_Err_T1);

				N1++;

				if(Xbj_T1L[i] < X_Min1)
					X_Min1 = Xbj_T1L[i];
				if(Xbj_T1L[i] > X_Max1)
					X_Max1 = Xbj_T1L[i];
				if(XS_T1 < Y_Min1)
					Y_Min1 = XS_T1;
				if(XS_T1 > Y_Max1)
					Y_Max1 = XS_T1;
				if(XS_MC_T1 < Y_Min1)
					Y_Min1 = XS_MC_T1;
				if(XS_T1 > Y_Max1)
					Y_Max1 = XS_MC_T1;
			}
		}
		for(int k = J1R+1; k<Count_T1R; k++){
			vEp1.push_back(Ep_T1R[k]);	
			vXbj1.push_back(Xbj_T1R[k]);	
			vXbj_Err1.push_back(0.0);	
			XS_T1 = XS_T1R[k];
			XS_Err_T1 = XS_Err_T1R[k];
			XS_MC_T1 = XS_MC_T1R[k];
			vXS1.push_back(XS_T1);
			vXS_Err1.push_back(XS_Err_T1);
			vXS_MC1.push_back( XS_MC_T1 );

			Y_EX_T1 = Y_EX_T1R[k];
			Y_MC_T1 = Y_MC_T1R[k];
			Y_EX_Err_T1 = Y_EX_Err_T1R[k];
			Y_MC_Err_T1 = Y_MC_Err_T1R[k];
			vY_EX1.push_back(Y_EX_T1);
			vY_MC1.push_back(Y_MC_T1);
			vY_EX_Err1.push_back(Y_EX_Err_T1);
			vY_MC_Err1.push_back(Y_MC_Err_T1);

			N1++;

			if(Xbj_T1R[k] < X_Min1)
				X_Min1 = Xbj_T1R[k];
			if(Xbj_T1R[k] > X_Max1)
				X_Max1 = Xbj_T1R[k];
			if(XS_T1 < Y_Min1)
				Y_Min1 = XS_T1;
			if(XS_T1 > Y_Max1)
				Y_Max1 = XS_T1;
			if(XS_MC_T1 < Y_Min1)
				Y_Min1 = XS_MC_T1;
			if(XS_T1 > Y_Max1)
				Y_Max1 = XS_MC_T1;
		}
	}
	else{
		for(int i=0;i<Count_T1L; i++){
			vEp1.push_back(Ep_T1L[i]);	
			vXbj1.push_back(Xbj_T1L[i]);	
			vXbj_Err1.push_back(0.0);	
			XS_T1 = XS_T1L[i];
			XS_Err_T1 = XS_Err_T1L[i];
			XS_MC_T1 = XS_MC_T1L[i];
			vXS1.push_back(XS_T1);
			vXS_Err1.push_back(XS_Err_T1);
			vXS_MC1.push_back( XS_MC_T1 );

			Y_EX_T1 = Y_EX_T1L[i];
			Y_MC_T1 = Y_MC_T1L[i];
			Y_EX_Err_T1 = Y_EX_Err_T1L[i];
			Y_MC_Err_T1 = Y_MC_Err_T1L[i];
			vY_EX1.push_back(Y_EX_T1);
			vY_MC1.push_back(Y_MC_T1);
			vY_EX_Err1.push_back(Y_EX_Err_T1);
			vY_MC_Err1.push_back(Y_MC_Err_T1);

			N1++;

			if(Xbj_T1L[i] < X_Min1)
				X_Min1 = Xbj_T1L[i];
			if(Xbj_T1L[i] > X_Max1)
				X_Max1 = Xbj_T1L[i];
			if(XS_T1 < Y_Min1)
				Y_Min1 = XS_T1;
			if(XS_T1 > Y_Max1)
				Y_Max1 = XS_T1;
			if(XS_MC_T1 < Y_Min1)
				Y_Min1 = XS_MC_T1;
			if(XS_T1 > Y_Max1)
				Y_Max1 = XS_MC_T1;
		}
	}
	/*}}}*/

	/*}}}*/

	/*Second Target{{{*/
	cerr<<endl<<"------------------ Second Target ----------------------------------------------------------"<<endl<<endl;

	/*Define{{{*/
	Double_t XS_T2L1[200], XS_Err_T2L1[200], Ep_T2L1[200], Xbj_T2L1[200],  Xbj_Err_T2L1[200],XS_MC_T2L1[200];
	Double_t XS_T2L2[200], XS_Err_T2L2[200], Ep_T2L2[200], Xbj_T2L2[200],  Xbj_Err_T2L2[200],XS_MC_T2L2[200];
	Double_t XS_T2R1[200], XS_Err_T2R1[200], Ep_T2R1[200], Xbj_T2R1[200],  Xbj_Err_T2R1[200],XS_MC_T2R1[200];
	Double_t XS_T2R2[200], XS_Err_T2R2[200], Ep_T2R2[200], Xbj_T2R2[200],  Xbj_Err_T2R2[200],XS_MC_T2R2[200];
	Double_t Y_EX_T2L1[200],  Y_EX_Err_T2L1[200],  Y_MC_T2L1[200],  Y_MC_Err_T2L1[200];
	Double_t Y_EX_T2L2[200],  Y_EX_Err_T2L2[200],  Y_MC_T2L2[200],  Y_MC_Err_T2L2[200];
	Double_t Y_EX_T2R1[200],  Y_EX_Err_T2R1[200],  Y_MC_T2R1[200],  Y_MC_Err_T2R1[200];
	Double_t Y_EX_T2R2[200],  Y_EX_Err_T2R2[200],  Y_MC_T2R2[200],  Y_MC_Err_T2R2[200];
	/*}}}*/

	/*First{{{*/
	TString infile_T2L1 = aInputFilePath + aName_T2L1 +"_"+ Com;
	cerr<<"==== Reading 1st file from "<<infile_T2L1<<endl;
	ifstream run_T2L1; run_T2L1.open(infile_T2L1);

	Int_t Count_T2L1 =0;
	Target = gGet_Target_Name(infile_T2L1.Data());
	if(Target=="He4")
		scale = He4_Norm;
	else
		scale = 1.0;

	Kin = gGet_Kin_Name(infile_T2L1.Data());

	TString Arm_T2L1="X";
	if(aName_T2L1.Contains("L_"))
		Arm_T2L1="L";
	else if(aName_T2L1.Contains("R_"))
		Arm_T2L1="R";


	if(Kin=="MM"){
		cerr<<"1st: Kin (3.1,3.2,4.1,4.2,5.0,5.05,5.1,5.2,6.5) = "; cin >> Kin;
	}
	if(Kin=="3.1"||Kin=="3.2")
		Theta = 21.00; 
	else if(Kin=="4.1"||Kin=="4.2")
		Theta = 23.00;
	else if(Kin=="5.0"||Kin=="5.05"||Kin=="5.1"||Kin=="5.2")
		Theta =  25.00;
	else if(Kin=="6.5")
		Theta = 28.00;
	else
		cerr<<"*** ERROR, unknown kinematics setting!"<<endl;

	TString KinName_T2L1 = Form("Kin%s",Kin.Data());
	gGet_Ep_Cut(KinName_T2L1.Data(), Arm_T2L1.Data(),&Ep_Min, &Ep_Max);
	cerr<<Form("@@@  For %s, Ep_Min = %f, Ep_Max = %f", KinName_T2L1.Data(), Ep_Min,Ep_Max)<<endl;

	Get_XS* xs_T2L1 = new Get_XS();
	xs_T2L1->Load_Table(Target.Data(),KinName_T2L1.Data());

	comment.ReadLine(run_T2L1);
	while ( run_T2L1>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
	{
		if ( Ep>=Ep_Min && Ep<=Ep_Max && angle<180 &&xbj>0.5  && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
		{
			Q2 = 4.0 * E0 * Ep * pow(sin(Theta*3.1415926/180./2.0),2);
			xbj = Q2/2.0/0.938272/(E0-Ep);
			Ep_T2L1[Count_T2L1]=Ep;
			XS_T2L1[Count_T2L1] = XS * scale;
			XS_Err_T2L1[Count_T2L1] = sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2));

			Y_EX_T2L1[Count_T2L1] = Yield_EX; 
			Y_MC_T2L1[Count_T2L1] = Yield_MC; 
			Y_EX_Err_T2L1[Count_T2L1] = sqrt(pow(Y_EX_Sys_Err,2)+pow(Y_EX_Stat_Err,2));
			Y_MC_Err_T2L1[Count_T2L1] = sqrt(pow(Y_MC_Sys_Err,2)+pow(Y_MC_Stat_Err,2));
				//if(Target.Contains("Ca"))
			if(1)
				XS_Err_T2L1[Count_T2L1] = pow(XS_Stat_Err,1);
			Xbj_T2L1[Count_T2L1] = xbj;    Xbj_Err_T2L1[Count_T2L1] = 0.0;
			XS_MC_T2L1[Count_T2L1] = xs_T2L1->gGet_XS_Born(E0,Ep,Theta);
			Count_T2L1++;
		}
		if ( !run_T2L1.good() )
			break;
	}
	run_T2L1.close();
	delete xs_T2L1;
	cerr << "  Total XS Points for 1st = " << Count_T2L1 <<endl;
	/*}}}*/

	/*Second{{{*/
	TString KinName_T2L2 ="";
	Bool_t bMore_T2L2 = kTRUE;
	TString infile_T2L2 = aInputFilePath + aName_T2L2 +"_"+ Com;
	cerr<<"==== Reading 2nd file from "<<infile_T2L2<<endl;
	TString Arm_T2L2="X";
	if(aName_T2L2 =="0")
		bMore_T2L2 = kFALSE;
	if(aName_T2L2.Contains("L_"))
		Arm_T2L2="L";
	else if(aName_T2L2.Contains("R_"))
		Arm_T2L2="R";

	Int_t Count_T2L2=0;
	if(bMore_T2L2){
		ifstream run_T2L2; run_T2L2.open(infile_T2L2);
		Target = gGet_Target_Name(infile_T2L2.Data());
		if(Target=="He4")
		scale = He4_Norm;
		else
			scale = 1.0;

		Kin = gGet_Kin_Name(infile_T2L2.Data());
		if(Kin=="MM"){
			cerr<<"2nd: Kin (3.1,3.2,4.1,4.2,5.0,5.05,5.1,5.2,6.5) = "; cin >> Kin;
		}
		if(Kin=="3.1"||Kin=="3.2")
			Theta = 21.00;
		else if(Kin=="4.1"||Kin=="4.2")
			Theta = 23.00;
		else if(Kin=="5.0"||Kin=="5.05"||Kin=="5.1"||Kin=="5.2")
			Theta =  25.00;
		else if(Kin=="6.5")
			Theta = 28.00;
		else
			cerr<<"*** ERROR, unknown kinematics setting!"<<endl;
		KinName_T2L2 = Form("Kin%s",Kin.Data());
		gGet_Ep_Cut(KinName_T2L2.Data(), Arm_T2L2.Data(),&Ep_Min, &Ep_Max);
		cerr<<Form("@@@  For %s, Ep_Min = %f, Ep_Max = %f", KinName_T2L2.Data(), Ep_Min,Ep_Max)<<endl;

		Get_XS* xs_T2L2 = new Get_XS();
		xs_T2L2->Load_Table(Target.Data(),KinName_T2L2.Data());

		comment.ReadLine(run_T2L2);
	while ( run_T2L2>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
		{
			if ( Ep>=Ep_Min && Ep<=Ep_Max && angle<180 &&xbj>0.5  && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
			{
				Q2 = 4.0 * E0 * Ep * pow(sin(Theta*3.1415926/180./2.0),2);
				xbj = Q2/2.0/0.938272/(E0-Ep);
				Ep_T2L2[Count_T2L2]=Ep;
				XS_T2L2[Count_T2L2] = XS * scale; 
				XS_Err_T2L2[Count_T2L2] = sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2));
	
				Y_EX_T2L2[Count_T2L2] = Yield_EX; 
				Y_MC_T2L2[Count_T2L2] = Yield_MC; 
				Y_EX_Err_T2L2[Count_T2L2] = sqrt(pow(Y_EX_Sys_Err,2)+pow(Y_EX_Stat_Err,2));
				Y_MC_Err_T2L2[Count_T2L2] = sqrt(pow(Y_MC_Sys_Err,2)+pow(Y_MC_Stat_Err,2));
				//	if(Target.Contains("Ca"))
				if(1)
					XS_Err_T2L2[Count_T2L2] = pow(XS_Stat_Err,1);
				Xbj_T2L2[Count_T2L2] = xbj;    Xbj_Err_T2L2[Count_T2L2] = 0.0;
				XS_MC_T2L2[Count_T2L2] = xs_T2L2->gGet_XS_Born(E0,Ep,Theta);
				Count_T2L2++;
			}
			if ( !run_T2L2.good() )
				break;
		}
		run_T2L2.close();
		delete xs_T2L2;
		cerr << "  Total XS Points for 2nd  = " << Count_T2L2 <<endl;
	}
	/*}}}*/

	/*Third{{{*/
	TString KinName_T2R1 ="";
	Bool_t bMore_T2R1 = kTRUE;
	TString infile_T2R1 = aInputFilePath + aName_T2R1 +"_"+ Com;
	cerr<<"==== Reading 3rd file from "<<infile_T2R1<<endl;
	TString Arm_T2R1="X";
	if(aName_T2R1 =="0")
		bMore_T2R1 = kFALSE;
	if(aName_T2R1.Contains("L_"))
		Arm_T2R1="L";
	else if(aName_T2R1.Contains("R_"))
		Arm_T2R1="R";

	Int_t Count_T2R1=0;
		bMore_T2R1 = kFALSE;
	if(bMore_T2R1){
		ifstream run_T2R1; run_T2R1.open(infile_T2R1);
		Target = gGet_Target_Name(infile_T2R1.Data());
		if(Target=="He4")
		scale = He4_Norm;
		else
			scale = 1.0;

		Kin = gGet_Kin_Name(infile_T2R1.Data());
		if(Kin=="MM"){
			cerr<<"2nd: Kin (3.1,3.2,4.1,4.2,5.0,5.05,5.1,5.2,6.5) = "; cin >> Kin;
		}
		if(Kin=="3.1"||Kin=="3.2")
			Theta = 21.00;
		else if(Kin=="4.1"||Kin=="4.2")
			Theta = 23.00;
		else if(Kin=="5.0"||Kin=="5.05"||Kin=="5.1"||Kin=="5.2")
			Theta =  25.00;
		else if(Kin=="6.5")
			Theta = 28.00;
		else
			cerr<<"*** ERROR, unknown kinematics setting!"<<endl;
		KinName_T2R1 = Form("Kin%s",Kin.Data());
		gGet_Ep_Cut(KinName_T2R1.Data(), Arm_T2R1.Data(),&Ep_Min, &Ep_Max);
		cerr<<Form("@@@  For %s, Ep_Min = %f, Ep_Max = %f", KinName_T2R1.Data(), Ep_Min,Ep_Max)<<endl;

		Get_XS* xs_T2R1 = new Get_XS();
		xs_T2R1->Load_Table(Target.Data(),KinName_T2R1.Data());

		comment.ReadLine(run_T2R1);
	while ( run_T2R1>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
		{
			if ( Ep>=Ep_Min && Ep<=Ep_Max && angle<180 &&xbj>0.5  && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
			{
				Q2 = 4.0 * E0 * Ep * pow(sin(Theta*3.1415926/180./2.0),2);
				xbj = Q2/2.0/0.938272/(E0-Ep);
				Ep_T2R1[Count_T2R1]=Ep;
				XS_T2R1[Count_T2R1] = XS* scale;   
				XS_Err_T2R1[Count_T2R1] = sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2));

				Y_EX_T2R1[Count_T2R1] = Yield_EX; 
				Y_MC_T2R1[Count_T2R1] = Yield_MC; 
				Y_EX_Err_T2R1[Count_T2R1] = sqrt(pow(Y_EX_Sys_Err,2)+pow(Y_EX_Stat_Err,2));
				Y_MC_Err_T2R1[Count_T2R1] = sqrt(pow(Y_MC_Sys_Err,2)+pow(Y_MC_Stat_Err,2));
				//	if(Target.Contains("Ca"))
				if(1)
					XS_Err_T2R1[Count_T2R1] = pow(XS_Stat_Err,1);
				Xbj_T2R1[Count_T2R1] = xbj;    Xbj_Err_T2R1[Count_T2R1] = 0.0;
				XS_MC_T2R1[Count_T2R1] = xs_T2R1->gGet_XS_Born(E0,Ep,Theta);
				Count_T2R1++;
			}
			if ( !run_T2R1.good() )
				break;
		}
		run_T2R1.close();
		delete xs_T2R1;
		cerr << "  Total XS Points for 2nd  = " << Count_T2R1 <<endl;
	}
	/*}}}*/

	/*Forth{{{*/
	TString KinName_T2R2 ="";
	Bool_t bMore_T2R2 = kTRUE;
	TString infile_T2R2 = aInputFilePath + aName_T2R2 +"_"+ Com;
	cerr<<"==== Reading 4th file from "<<infile_T2R2<<endl;
	TString Arm_T2R2="X";
	if(aName_T2R2 =="0")
		bMore_T2R2 = kFALSE;
	if(aName_T2R2.Contains("L_"))
		Arm_T2R2="L";
	else if(aName_T2R2.Contains("R_"))
		Arm_T2R2="R";

	Int_t Count_T2R2=0;
		bMore_T2R2 = kFALSE;
	if(bMore_T2R2){
		ifstream run_T2R2; run_T2R2.open(infile_T2R2);
		Target = gGet_Target_Name(infile_T2R2.Data());
		if(Target=="He4")
		scale = He4_Norm;
		else
			scale = 1.0;

		Kin = gGet_Kin_Name(infile_T2R2.Data());
		if(Kin=="MM"){
			cerr<<"2nd: Kin (3.1,3.2,4.1,4.2,5.0,5.05,5.1,5.2,6.5) = "; cin >> Kin;
		}
		if(Kin=="3.1"||Kin=="3.2")
			Theta = 21.00;
		else if(Kin=="4.1"||Kin=="4.2")
			Theta = 23.00;
		else if(Kin=="5.0"||Kin=="5.05"||Kin=="5.1"||Kin=="5.2")
			Theta =  25.00;
		else if(Kin=="6.5")
			Theta = 28.00;
		else
			cerr<<"*** ERROR, unknown kinematics setting!"<<endl;
		KinName_T2R2 = Form("Kin%s",Kin.Data());
		gGet_Ep_Cut(KinName_T2R2.Data(), Arm_T2R2.Data(),&Ep_Min, &Ep_Max);
		cerr<<Form("@@@  For %s, Ep_Min = %f, Ep_Max = %f", KinName_T2R2.Data(), Ep_Min,Ep_Max)<<endl;

		Get_XS* xs_T2R2 = new Get_XS();
		xs_T2R2->Load_Table(Target.Data(),KinName_T2R2.Data());

		comment.ReadLine(run_T2R2);
	while ( run_T2R2>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
		{
			if ( Ep>=Ep_Min && Ep<=Ep_Max && angle<180 &&xbj>0.5  && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
			{
				Q2 = 4.0 * E0 * Ep * pow(sin(Theta*3.1415926/180./2.0),2);
				xbj = Q2/2.0/0.938272/(E0-Ep);
				Ep_T2R2[Count_T2R2]=Ep;
				XS_T2R2[Count_T2R2] = XS * scale;     
				XS_Err_T2R2[Count_T2R2] = sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,1));

				Y_EX_T2R2[Count_T2R2] = Yield_EX; 
				Y_MC_T2R2[Count_T2R2] = Yield_MC; 
				Y_EX_Err_T2R2[Count_T2R2] = sqrt(pow(Y_EX_Sys_Err,2)+pow(Y_EX_Stat_Err,2));
				Y_MC_Err_T2R2[Count_T2R2] = sqrt(pow(Y_MC_Sys_Err,2)+pow(Y_MC_Stat_Err,2));
				//	if(Target.Contains("Ca"))
				if(1)
					XS_Err_T2R2[Count_T2R2] = pow(XS_Stat_Err,1);
				Xbj_T2R2[Count_T2R2] = xbj;    Xbj_Err_T2R2[Count_T2R2] = 0.0;
				XS_MC_T2R2[Count_T2R2] = xs_T2R2->gGet_XS_Born(E0,Ep,Theta);
				Count_T2R2++;
			}
			if ( !run_T2R2.good() )
				break;
		}
		run_T2R2.close();
		delete xs_T2R2;
		cerr << "  Total XS Points for 2nd  = " << Count_T2R2 <<endl;
	}
	/*}}}*/

	/*Merge{{{*/
	//Add First Two Input
	Double_t Ep_T2L[200],Xbj_T2L[200], XS_T2L[200], XS_Err_T2L[200], XS_MC_T2L[200];
	Double_t Y_EX_T2L[200],Y_EX_Err_T2L[200],Y_MC_T2L[200],Y_MC_Err_T2L[200];
	Double_t XS_Sum_T2L = 0.0, XS_Err_Sum_T2L = 0.0; Int_t fSum_T2L = 0;
	Double_t Y_EX_Sum_T2L = 0.0, Y_EX_Err_Sum_T2L = 0.0, Y_MC_Sum_T2L = 0.0, Y_MC_Err_Sum_T2L = 0.0;
	Int_t Count_T2L = 0, J_T2L = 0;
	for(int i=0;i<Count_T2L1;i++){
		XS_Sum_T2L = XS_T2L1[i];
		XS_Err_Sum_T2L = pow(XS_Err_T2L1[i]/XS_T2L1[i],2);

		Y_EX_Sum_T2L = Y_EX_T2L1[i];
		Y_MC_Sum_T2L = Y_MC_T2L1[i];
		Y_EX_Err_Sum_T2L = pow(Y_EX_Err_T2L1[i]/Y_EX_T2L1[i],2);
		Y_MC_Err_Sum_T2L = pow(Y_MC_Err_T2L1[i]/Y_MC_T2L1[i],2);
		fSum_T2L = 1;
		if(bMore_T2L2){
			for(int j=0;j<Count_T2L2;j++){
				//Find the verlapping bin for two setting, such as Kin3.1 and Kin3.2
				//if(abs(Xbj_T2L1[i]-Xbj_T2L2[j])< 0.0001){
				if(abs(Ep_T2L1[i]-Ep_T2L2[j])< 0.0001){
					XS_Sum_T2L += XS_T2L2[j];
					XS_Err_Sum_T2L +=	pow(XS_Err_T2L2[j]/XS_T2L2[j],2);

					Y_EX_Sum_T2L += Y_EX_T2L2[j];
					Y_MC_Sum_T2L += Y_MC_T2L2[j];
					Y_EX_Err_Sum_T2L += pow(Y_EX_Err_T2L2[j]/Y_EX_T2L2[j],2);
					Y_MC_Err_Sum_T2L += pow(Y_MC_Err_T2L2[j]/Y_MC_T2L2[j],2);
					fSum_T2L++; J_T2L = j;
					break; //No two xbj bin within one setting
				}
			}
		}
		Ep_T2L[i] = Ep_T2L1[i];
		Xbj_T2L[i] = Xbj_T2L1[i];
		XS_MC_T2L[i] = XS_MC_T2L1[i];
		XS_T2L[i] = XS_Sum_T2L / fSum_T2L;
		XS_Err_T2L[i] = XS_T2L[i] * sqrt(XS_Err_Sum_T2L); 
	
		Y_EX_T2L[i] = Y_EX_Sum_T2L / fSum_T2L;
		Y_MC_T2L[i] = Y_MC_Sum_T2L / fSum_T2L;
		Y_EX_Err_T2L[i] = Y_EX_T2L[i] * sqrt(Y_EX_Err_Sum_T2L); 
		Y_MC_Err_T2L[i] = Y_MC_T2L[i] * sqrt(Y_MC_Err_Sum_T2L); 
			
		Count_T2L++;
	}
	for(int k=J_T2L+1;k<Count_T2L2;k++){
		Ep_T2L[Count_T2L] = Ep_T2L2[k];
		Xbj_T2L[Count_T2L] = Xbj_T2L2[k];
		XS_T2L[Count_T2L] = XS_T2L2[k];
		XS_Err_T2L[Count_T2L] = XS_Err_T2L2[k];
		XS_MC_T2L[Count_T2L] = XS_MC_T2L2[k];

		Y_EX_T2L[Count_T2L] = Y_EX_T2L2[k];
		Y_MC_T2L[Count_T2L] = Y_MC_T2L2[k];
		Y_EX_Err_T2L[Count_T2L] = Y_EX_Err_T2L2[k];
		Y_MC_Err_T2L[Count_T2L] = Y_MC_Err_T2L2[k];

		Count_T2L++;	
	}

	//Add Last Two Input
	Double_t Ep_T2R[200], Xbj_T2R[200], XS_T2R[200], XS_Err_T2R[200], XS_MC_T2R[200];
	Double_t XS_Sum_T2R = 0.0, XS_Err_Sum_T2R = 0.0; Int_t fSum_T2R = 0;
	Double_t Y_EX_T2R[200],Y_EX_Err_T2R[200],Y_MC_T2R[200],Y_MC_Err_T2R[200];
	Double_t Y_EX_Sum_T2R = 0.0, Y_EX_Err_Sum_T2R = 0.0, Y_MC_Sum_T2R = 0.0, Y_MC_Err_Sum_T2R = 0.0;
	Int_t Count_T2R = 0, J_T2R = 0;
	if(bMore_T2R1){
		for(int i=0;i<Count_T2R1;i++){
			XS_Sum_T2R = XS_T2R1[i];
			XS_Err_Sum_T2R = pow(XS_Err_T2R1[i]/XS_T2R1[i],2);
		
			Y_EX_Sum_T2R = Y_EX_T2R1[i];
			Y_MC_Sum_T2R = Y_MC_T2R1[i];
			Y_EX_Err_Sum_T2R = pow(Y_EX_Err_T2R1[i]/Y_EX_T2R1[i],2);
			Y_MC_Err_Sum_T2R = pow(Y_MC_Err_T2R1[i]/Y_MC_T2R1[i],2);
			fSum_T2R = 1;
			if(bMore_T2R2){
				for(int j=0;j<Count_T2R2;j++){
					//Find the verlapping bin for two setting, such as Kin3.1 and Kin3.2
					//if(abs(Xbj_T2R1[i]-Xbj_T2R2[j])< 0.0001){
					if(abs(Ep_T2R1[i]-Ep_T2R2[j])< 0.0001){
						XS_Sum_T2R += XS_T2R2[j];
						XS_Err_Sum_T2R +=	pow(XS_Err_T2R2[j]/XS_T2R2[j],2);

						Y_EX_Sum_T2R += Y_EX_T2R2[j];
						Y_MC_Sum_T2R += Y_MC_T2R2[j];
						Y_EX_Err_Sum_T2R += pow(Y_EX_Err_T2R2[j]/Y_EX_T2R2[j],2);
						Y_MC_Err_Sum_T2R += pow(Y_MC_Err_T2R2[j]/Y_MC_T2R2[j],2);
						fSum_T2R++; J_T2R = j;
						break; //No two xbj bin within one setting
					}
				}
			}
			Ep_T2R[i] = Ep_T2R1[i];
			Xbj_T2R[i] = Xbj_T2R1[i];
			XS_MC_T2R[i] = XS_MC_T2R1[i];
			XS_T2R[i] = XS_Sum_T2R / fSum_T2R;
			XS_Err_T2R[i] = XS_T2R[i] * sqrt(XS_Err_Sum_T2R); 

			Y_EX_T2R[i] = Y_EX_Sum_T2R / fSum_T2R;
			Y_MC_T2R[i] = Y_MC_Sum_T2R / fSum_T2R;
			Y_EX_Err_T2R[i] = Y_EX_T2R[i] * sqrt(Y_EX_Err_Sum_T2R); 
			Y_MC_Err_T2R[i] = Y_MC_T2R[i] * sqrt(Y_MC_Err_Sum_T2R); 

			Count_T2R++;
		}
		for(int k=J_T2R+1;k<Count_T2R2;k++){
			Ep_T2R[Count_T2R] = Ep_T2R2[k];
			Xbj_T2R[Count_T2R] = Xbj_T2R2[k];
			XS_T2R[Count_T2R] = XS_T2R2[k];
			XS_MC_T2R[Count_T2R] = XS_MC_T2R2[k];
			XS_Err_T2R[Count_T2R] = XS_Err_T2R2[k];

			Y_EX_T2R[Count_T2R] = Y_EX_T2R2[k];
			Y_MC_T2R[Count_T2R] = Y_MC_T2R2[k];
			Y_EX_Err_T2R[Count_T2R] = Y_EX_Err_T2R2[k];
			Y_MC_Err_T2R[Count_T2R] = Y_MC_Err_T2R2[k];

			Count_T2R++;	
		}
	}

	vector<double> vEp2, vXbj2, vXbj_Err2, vXS2, vXS_Err2, vXS_MC2;
	vector<double> vY_EX2, vY_EX_Err2, vY_MC2, vY_MC_Err2;
	double XS_MC_T2, XS_T2,XS_Err_T2;
	double Y_EX_T2, Y_EX_Err_T2, Y_MC_T2, Y_MC_Err_T2;
	double X_Min2 = 100.0, X_Max2 = -100.0, Y_Min2 = 10000.0, Y_Max2 = -10000.0;
	int N2 = 0, J2R = 0;
	bool bFound_T2L = kFALSE;
	if(bMore_T2R1){
		for(int i=0;i<Count_T2L; i++){
			bFound_T2L = kFALSE;
			for(int j=0;j<Count_T2R; j++){
				//if( abs(Xbj_T2L[i] - Xbj_T2R[j])<0.0001 ){
				if( abs(Ep_T2L[i] - Ep_T2R[j])<0.0001 ){
					bFound_T2L = kTRUE;
					vEp2.push_back(Ep_T2L[i]);	
					vXbj2.push_back(Xbj_T2L[i]);	
					vXbj_Err2.push_back(0.0);	
					XS_T2 = (XS_T2L[i]+XS_T2R[j])/2.0;
					XS_Err_T2 = XS_T2 * sqrt( pow(XS_Err_T2L[i]/XS_T2L[i],2) + pow(XS_Err_T2R[j]/XS_T2R[j],2) );
					XS_MC_T2 = (XS_MC_T2L[i]+XS_MC_T2R[j])/2.0;
					vXS2.push_back(XS_T2);
					vXS_Err2.push_back(XS_Err_T2);
					vXS_MC2.push_back( XS_MC_T2 );

					Y_EX_T2 = (Y_EX_T2L[i]+Y_EX_T2R[j])/2.0;
					Y_EX_Err_T2 = Y_EX_T2 * sqrt( pow(Y_EX_Err_T2L[i]/Y_EX_T2L[i],2) + pow(Y_EX_Err_T2R[j]/Y_EX_T2R[j],2) );
					Y_MC_T2 = (Y_MC_T2L[i]+Y_MC_T2R[j])/2.0;
					Y_MC_Err_T2 = Y_MC_T2 * sqrt( pow(Y_MC_Err_T2L[i]/Y_MC_T2L[i],2) + pow(Y_MC_Err_T2R[j]/Y_MC_T2R[j],2) );
                    vY_EX2.push_back(Y_EX_T2);
                    vY_MC2.push_back(Y_MC_T2);
                    vY_EX_Err2.push_back(Y_EX_Err_T2);
                    vY_MC_Err2.push_back(Y_MC_Err_T2);

					N2++; J2R = j;

					if(Xbj_T2L[i] < X_Min2)
						X_Min2 = Xbj_T2L[i];
					if(Xbj_T2L[i] > X_Max2)
						X_Max2 = Xbj_T2L[i];
					if(XS_T2 < Y_Min2)
						Y_Min2 = XS_T2;
					if(XS_T2 > Y_Max2)
						Y_Max2 = XS_T2;
					if(XS_MC_T2 < Y_Min2)
						Y_Min2 = XS_MC_T2;
					if(XS_T2 > Y_Max2)
						Y_Max2 = XS_MC_T2;
				}
			}
			if(!bFound_T2L){
				vEp2.push_back(Ep_T2L[i]);	
				vXbj2.push_back(Xbj_T2L[i]);	
				vXbj_Err2.push_back(0.0);	
				XS_T2 = XS_T2L[i];
				XS_Err_T2 = XS_Err_T2L[i];
				XS_MC_T2 = XS_MC_T2L[i];
				vXS2.push_back(XS_T2);
				vXS_Err2.push_back(XS_Err_T2);
				vXS_MC2.push_back( XS_MC_T2 );

				Y_EX_T2 = Y_EX_T2L[i];
				Y_MC_T2 = Y_MC_T2L[i];
				Y_EX_Err_T2 = Y_EX_Err_T2L[i];
				Y_MC_Err_T2 = Y_MC_Err_T2L[i];
				vY_EX2.push_back(Y_EX_T2);
				vY_MC2.push_back(Y_MC_T2);
				vY_EX_Err2.push_back(Y_EX_Err_T2);
				vY_MC_Err2.push_back(Y_MC_Err_T2);

				N2++;

				if(Xbj_T2L[i] < X_Min2)
					X_Min2 = Xbj_T2L[i];
				if(Xbj_T2L[i] > X_Max2)
					X_Max2 = Xbj_T2L[i];
				if(XS_T2 < Y_Min2)
					Y_Min2 = XS_T2;
				if(XS_T2 > Y_Max2)
					Y_Max2 = XS_T2;
				if(XS_MC_T2 < Y_Min2)
					Y_Min2 = XS_MC_T2;
				if(XS_T2 > Y_Max2)
					Y_Max2 = XS_MC_T2;
			}
		}
		for(int k = J2R+1; k<Count_T2R; k++){
			vEp2.push_back(Ep_T2R[k]);	
			vXbj2.push_back(Xbj_T2R[k]);	
			vXbj_Err2.push_back(0.0);	
			XS_T2 = XS_T2R[k];
			XS_Err_T2 = XS_Err_T2R[k];
			XS_MC_T2 = XS_MC_T2R[k];
			vXS2.push_back(XS_T2);
			vXS_Err2.push_back(XS_Err_T2);
			vXS_MC2.push_back( XS_MC_T2 );

			Y_EX_T2 = Y_EX_T2R[k];
			Y_MC_T2 = Y_MC_T2R[k];
			Y_EX_Err_T2 = Y_EX_Err_T2R[k];
			Y_MC_Err_T2 = Y_MC_Err_T2R[k];
			vY_EX2.push_back(Y_EX_T2);
			vY_MC2.push_back(Y_MC_T2);
			vY_EX_Err2.push_back(Y_EX_Err_T2);
			vY_MC_Err2.push_back(Y_MC_Err_T2);

			N2++;

			if(Xbj_T2R[k] < X_Min2)
				X_Min2 = Xbj_T2R[k];
			if(Xbj_T2R[k] > X_Max2)
				X_Max2 = Xbj_T2R[k];
			if(XS_T2 < Y_Min2)
				Y_Min2 = XS_T2;
			if(XS_T2 > Y_Max2)
				Y_Max2 = XS_T2;
			if(XS_MC_T2 < Y_Min2)
				Y_Min2 = XS_MC_T2;
			if(XS_T2 > Y_Max2)
				Y_Max2 = XS_MC_T2;
		}
	}
	else{
		for(int i=0;i<Count_T2L; i++){
			vEp2.push_back(Ep_T2L[i]);	
			vXbj2.push_back(Xbj_T2L[i]);	
			vXbj_Err2.push_back(0.0);	
			XS_T2 = XS_T2L[i];
			XS_Err_T2 = XS_Err_T2L[i];
			XS_MC_T2 = XS_MC_T2L[i];
			vXS2.push_back(XS_T2);
			vXS_Err2.push_back(XS_Err_T2);
			vXS_MC2.push_back( XS_MC_T2 );

			Y_EX_T2 = Y_EX_T2L[i];
			Y_MC_T2 = Y_MC_T2L[i];
			Y_EX_Err_T2 = Y_EX_Err_T2L[i];
			Y_MC_Err_T2 = Y_MC_Err_T2L[i];
			vY_EX2.push_back(Y_EX_T2);
			vY_MC2.push_back(Y_MC_T2);
			vY_EX_Err2.push_back(Y_EX_Err_T2);
			vY_MC_Err2.push_back(Y_MC_Err_T2);

			N2++;

			if(Xbj_T2L[i] < X_Min2)
				X_Min2 = Xbj_T2L[i];
			if(Xbj_T2L[i] > X_Max2)
				X_Max2 = Xbj_T2L[i];
			if(XS_T2 < Y_Min2)
				Y_Min2 = XS_T2;
			if(XS_T2 > Y_Max2)
				Y_Max2 = XS_T2;
			if(XS_MC_T2 < Y_Min2)
				Y_Min2 = XS_MC_T2;
			if(XS_T2 > Y_Max2)
				Y_Max2 = XS_MC_T2;
		}
	}
	/*}}}*/

	/*}}}*/

	/*Ratio{{{*/  
	ofstream outfile(Form("Y_%s_%s_%s.dat",Target0.Data(),Target.Data(),Kin0.Data()));
	outfile << Form("%7s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s",
		   	"Xbj","Ratio","Ratio_Err","Y_EX1","Y_EX_Err1","Y_EX2","Y_EX_Err2", "Y_MC1","Y_MC_Err1","Y_MC2","Y_MC_Err2" )<<endl;

	Int_t Count=0;
	Double_t x1,x2,ep1,ep2,ratio,ratio_mc,err, err_mc, max;
	vector<Double_t> vEp, vXbj, vXbj_Err,vR, vR_Err,vR_MC, vR_MC_Err;

	Double_t X_Min = 10, X_Max = -10, R_Min = 1000.0, R_Max = -1000.0;
	int i,j; max=-1.0;
	for( i=0;i<N1;i++){
		x1 = vXbj1[i];
		ep1 = vEp1[i];
		for( j=0;j<N2;j++){
			x2 = vXbj2[j];
			ep2 = vEp2[j];
			if(abs(ep1-ep2)<0.0001){
				ratio    = factor * vY_EX1[i]/vY_EX2[j];
				ratio_mc = factor * vY_MC1[i]/vY_MC2[j];
				err    = ratio    * sqrt( pow(vY_EX_Err2[j]/vY_EX2[j],2) + pow(vY_EX_Err1[i]/vY_EX1[i],2) );
				err_mc = ratio_mc * sqrt( pow(vY_MC_Err2[j]/vY_MC2[j],2) + pow(vY_MC_Err1[i]/vY_MC1[i],2) );
				if(ratio>1e-30 && ratio<100.&& (err/ratio)<0.5){
					vEp.push_back(ep1);
					vXbj.push_back(x1);
					vXbj_Err.push_back(0.0);
					vR.push_back(ratio);
					vR_Err.push_back(err);
					vR_MC.push_back(ratio_mc);
					vR_MC_Err.push_back(err_mc);
					cerr <<Form("-----> Found one: xbj1=%f, xbj2=%f, Y_EX1=%e, Y_EX2=%e,Ratio=%f ",x1,x2,vY_EX1[i],vY_EX2[j],ratio)<<endl;
					outfile << Form("%7.5f %14e %14e %14e %14e %14e %14e %14e %14e %14e %14e", 
							x1,ratio, err, vY_EX1[i], vY_EX_Err1[i], vY_EX2[j], vY_EX_Err2[j], vY_MC1[i], vY_MC_Err1[i], vY_MC2[j], vY_MC_Err2[j] )<<endl;

					if(R_Min > ratio)
						R_Min = ratio;
					if(R_Min > ratio_mc)
						R_Min = ratio_mc;
					if(R_Max < ratio)
						R_Max = ratio;
					//				if(R_Max < ratio_mc)
					//					R_Max = ratio_mc;
					if(X_Min > x1)
						X_Min = x1;
					if(X_Max < x1)
						X_Max = x1;
					Count++;
				}
			}
		}
	}

	//Return the value
	for(int k=0; k<Count; k++){
		aXbj[k] = vXbj[k];
		aR_EX[k] = vR[k];
		aR_EX_Err[k] = vR_Err[k];
		aR_MC[k] = vR_MC[k];
		aR_MC_Err[k] = vR_MC_Err[k];
	}
	aN[0] = Count;
	outfile.close();
	/*}}}*/
}
/*}}}*/

/*gGetTargetAZ{{{*/
void gGetTargetAZ(const TString& aTarget, int* aA, int* aZ){
   if(aTarget=="H2"||aTarget=="D2"){
     aA[0] = 2; aZ[0] = 1;
   }
   else if(aTarget=="He3"){
     aA[0] = 3; aZ[0] = 2;
   }
   else if(aTarget=="He4"){
     aA[0] = 4; aZ[0] = 2;
   }
   else if(aTarget=="C12"){
     aA[0] = 12; aZ[0] = 6;
   }
   else if(aTarget=="Ca40"){
     aA[0] = 40; aZ[0] = 20;
   }
   else if(aTarget=="Ca48"){
     aA[0] = 48; aZ[0] = 20;
   }
   else{
	   cerr<<" **** I don't know this target:  "<<aTarget.Data()<<endl;
	   aA[0]=0; aZ[0]=0;
   } 
}
/*}}}*/

