#include "XGT2_Main.h"
#include "/work/halla/e08014/disk1/yez/Ep_Cross_Section/SRC/XGT2_XEMC.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"

//const double He4_Norm = 1.0/0.864/0.864;
const double He4_Norm = 1.0;

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

/*Get He3 RadCorr{{{*/
inline double gGet_RadCorr_He3(const double aXbj){
  double aRadCor = 1-0.1*(aXbj-2.45)+0.025*pow((aXbj-2.45),2)-0.72*pow((aXbj-2.45),3);
  return 1.0/aRadCor;
}
/*}}}*/

void Merge(const TString& aInputFileName, double* aXbj, double* aR_EX, double* aR_EX_Err, double* aR_MC, int* aN);
void E02019(const TString& aInputFileName, double aTheta, double* aXbj, double* aR_EX, double* aR_EX_Err, int* aN);
void E02019_A_He3(const TString& aTarget, const double aTheta,double* aXbj, double* aR_EX, double* aR_EX_Err_Up,double* aR_EX_Err_Down, int* aN);
void E02019_He4_He3(const TString& aTarget, const double aTheta,double* aXbj, double* aR_EX, double* aR_EX_Err_Up,double* aR_EX_Err_Down, int* aN);
void CLAS_He4_He3( double* aXbj_CS, double* aR_EX_CS, double* aR_EX_Err_CS, int *aN_CS);
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
	gStyle->SetOptFit(1);
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

	/*Plot*/
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

	double XMax = 3.2;
	if(Target=="He4")
		XMax=3.2;
	if(Target=="He3")
		XMax=3.2;
	if(Target=="H2")
		XMax=2.2;

	if(Target=="He3" &&Target0=="C12")
		RMAX = 20.0;
	if(Target=="He3" &&Target0=="He4")
		RMAX = 12.0;
	if(Target=="He4" &&Target0=="C12")
		RMAX = 7;
	if(Target=="Ca40" &&Target0=="Ca48"){
		RMIN = 0.8;		RMAX = 1.1;
	}
	if(Target=="C12" &&Target0=="Ca40")
		RMAX = 1.5;
	if(Target=="C12" &&Target0=="Ca48")
		RMAX = 1.5;
	if(Target=="H2" &&Target0=="He3")
		RMAX = 3.0;
	if(Target=="H2" &&Target0=="He4")
		RMAX = 5.1;
	if(Target=="H2" &&Target0=="C12"){
		RMIN=0.1; RMAX = 7.5;
	}
	if(Target=="H2" &&Target0=="Ca40"){
		RMIN=0.1; RMAX = 9.5;
	}
	if(Target=="H2" &&Target0=="Ca48"){
		RMIN=0.1; RMAX = 9.5;
	}

	TCanvas *c3 = new TCanvas("c3","c3",1200,700);
	c3->cd();
	TH2F *h3 = new TH2F("h3","",300,0.7,XMax,300,RMIN,RMAX);
	h3->SetStats(kFALSE);
	h3->SetXTitle("x_{bj}");
	h3->GetXaxis()->CenterTitle(1);
	h3->GetXaxis()->SetTitleFont(32);
	h3->GetXaxis()->SetTitleSize(0.06);
	h3->GetYaxis()->SetTitle(Form("#frac{#sigma_{%s}}{#sigma_{%s}} #times #frac{%d}{%d}",Target0.Data(),Target.Data(),A2,A1));
	h3->GetYaxis()->CenterTitle(1);
	h3->GetYaxis()->SetTitleFont(32);
	h3->GetYaxis()->SetTitleSize(0.04);
	h3->GetYaxis()->CenterTitle();
	gPad->SetLogy(0);
	h3->Draw();

	if(Target=="He3" &&Target0=="C12")
		RMAX = 5.0;
	if(Target=="He3" &&Target0=="He4")
		RMAX = 4.0;
	if(Target=="He4" &&Target0=="C12")
		RMAX = 4;
	if(Target=="Ca40" &&Target0=="Ca48"){
		RMAX = 1.1; RMIN=0.8;
	}
	if(Target=="C12" &&Target0=="Ca40"){
		RMAX = 1.5; RMIN=0.7;
	}
	if(Target=="C12" &&Target0=="Ca48"){
		RMAX = 1.5; RMIN=0.7;
	}
	if(Target=="H2" &&Target0=="C12"){
		 RMIN=0.1;
	}

	TCanvas *c4 = new TCanvas("c4","c4",1200,700);
	c4->cd();
	TH2F *h4 = new TH2F("h4","",300,0.7,XMax,300,RMIN,RMAX);
	h4->SetStats(kFALSE);
	h4->SetXTitle("x_{bj}");
	h4->GetXaxis()->CenterTitle(1);
	h4->GetXaxis()->SetTitleFont(32);
	h4->GetXaxis()->SetTitleSize(0.08);
	h4->GetXaxis()->SetTitleOffset(0.7);
	h4->GetYaxis()->SetTitle(Form("#frac{#sigma_{%s}}{#sigma_{%s}} #times #frac{%d}{%d}",Target0.Data(),Target.Data(),A2,A1));
	h4->GetYaxis()->CenterTitle(1);
	h4->GetYaxis()->SetTitleFont(32);
	h4->GetYaxis()->SetTitleSize(0.07);
	h4->GetYaxis()->SetTitleOffset(0.6);
	h4->GetYaxis()->CenterTitle();
	gPad->SetLogy(0);
	h4->Draw();

	TLegend *l1 = new TLegend(0.17,0.64,0.50, 0.94,Form("#sigma_{%s} / #sigma_{%s} Ratio: ",Target0.Data(),Target.Data()));
	l1->SetTextSize(0.035);
	l1->SetTextFont(22);

	TGraphErrors *ex21 = new TGraphErrors(N_21,xbj_21, R_EX_21, Zero,R_EX_Err_21);
	ex21->SetMarkerStyle(24);
	ex21->SetMarkerSize(1.3);
	ex21->SetMarkerColor(2);
	ex21->SetLineColor(2);
	Theta = 21.0;
	c3->cd();	ex21->Draw("P");
	c4->cd();	ex21->Draw("P");
	l1->AddEntry(ex21, Form("%2.0f^{o}: %4.3f<Q^{2}<%4.3f (GeV/c^{2})",Theta,Q2_21_D,Q2_21_U),"p");
//	ex21->Fit("pol0","","",1.75, 2.0);

	TGraphErrors *ex23 = new TGraphErrors(N_23,xbj_23, R_EX_23, Zero,R_EX_Err_23);
	ex23->SetMarkerStyle(25);
	ex23->SetMarkerSize(1.3);
	ex23->SetMarkerColor(4);
	ex23->SetLineColor(4);
	Theta = 23.0;
	if(!(Target=="H2")){
		c3->cd();	ex23->Draw("P");
		c4->cd();	ex23->Draw("P");
		l1->AddEntry(ex23, Form("%2.0f^{o}: %4.3f<Q^{2}<%4.3f (GeV/c^{2})",Theta,Q2_23_D,Q2_23_U),"p");
	//	ex23->Fit("pol0","","",1.7, 2.0);
	}

	TGraphErrors *exQE = new TGraphErrors(N_QE,xbj_QE, R_EX_QE, Zero,R_EX_Err_QE);
	exQE->SetMarkerStyle(26);
	exQE->SetMarkerSize(1.3);
	exQE->SetMarkerColor(6);
	exQE->SetLineColor(6);
	if(!(Target=="H2")){
		c3->cd(); exQE->Draw("P");
		c4->cd(); exQE->Draw("P");
	}

	TGraphErrors *ex25 = new TGraphErrors(N_25,xbj_25, R_EX_25, Zero,R_EX_Err_25);
	ex25->SetMarkerStyle(26);
	ex25->SetMarkerSize(1.3);
	ex25->SetMarkerColor(6);
	ex25->SetLineColor(6);
	Theta = 25.0;
	c3->cd(); ex25->Draw("P");
	c4->cd(); ex25->Draw("P");
	l1->AddEntry(ex25, Form("%2.0f^{o}: %4.3f<Q^{2}<%4.3f (GeV/c^{2})",Theta,Q2_25_D,Q2_25_U),"p");
//	ex25->Fit("pol0","","",1.5, 1.8);

	TGraphErrors *ex28 = new TGraphErrors(N_28,xbj_28, R_EX_28, Zero,R_EX_Err_28);
	ex28->SetMarkerStyle(27);
	ex28->SetMarkerSize(1.7);
	ex28->SetMarkerColor(1);
	ex28->SetLineColor(1);
	Theta = 28.0;
	if(Target0=="C12" && Target=="He3"){
		c3->cd(); 	ex28->Draw("P");
		c4->cd(); 	ex28->Draw("P");
		//ex28->Fit("pol0","","",1.7, 2.0);
		l1->AddEntry(ex28, Form("%2.0f^{o}: %4.3f<Q^{2}<%4.3f (GeV/c^{2})",Theta,Q2_28_D,Q2_28_U),"p");
	}

	//if(!(Target0=="C12"&&Target=="He4")){
	if(Target=="He3"){
		//TGraphErrors *exNF = new TGraphErrors(N_NF,xbj_NF, R_EX_NF, Zero,R_EX_Err_NF);
		TGraphAsymmErrors *exNF = new TGraphAsymmErrors(N_NF,xbj_NF, R_EX_NF, Zero,Zero,R_EX_Err_Down,R_EX_Err_Up);
		exNF->SetMarkerStyle(30);
		exNF->SetMarkerSize(1.6);
		exNF->SetMarkerColor(8);
		exNF->SetLineColor(8);
		c3->cd(); exNF->Draw("P");
		c4->cd(); exNF->Draw("P");
		//l1->AddEntry(exNF,Form("E02019 Data (%2.0f^{0})",Theta_NF),"p");
		l1->AddEntry(exNF,Form("E02019 Data (18^{0},Q^{2}=2.7 GeV/c^{2})",Theta_NF),"p");
	}

	//CLAS Data
    cerr<<"---- Plotting CLAS Data"<<endl;
	if(Target0=="He4"&&Target=="He3"){
		TGraphErrors *exCS = new TGraphErrors(N_CS,xbj_CS, R_EX_CS, Zero,R_EX_Err_CS);
		exCS->SetMarkerStyle(20);
		exCS->SetMarkerSize(0.85);
		exCS->SetMarkerColor(7);
		exCS->SetLineColor(7);
		c3->cd(); exCS->Draw("P");
		c4->cd(); exCS->Draw("P");
		l1->AddEntry(exCS,Form("CLAS Data"),"p");
	}

    gStyle->SetOptFit(0);
    int N_ALL = N_21+N_23+N_25+N_28;
	double R_EX_ALL[N_ALL], R_EX_Err_ALL[N_ALL], xbj_ALL[N_ALL];
	int count = 0;
	for(int i=0;i<N_21;i++){
         R_EX_ALL[count] = R_EX_21[i];
         R_EX_Err_ALL[count] = R_EX_Err_21[i];
         xbj_ALL[count] = xbj_21[i];
		 count++;
	}
	for(int i=0;i<N_23;i++){
         R_EX_ALL[count] = R_EX_23[i];
         R_EX_Err_ALL[count] = R_EX_Err_23[i];
         xbj_ALL[count] = xbj_23[i];
		 count++;
	}
	for(int i=0;i<N_25;i++){
         R_EX_ALL[count] = R_EX_25[i];
         R_EX_Err_ALL[count] = R_EX_Err_25[i];
         xbj_ALL[count] = xbj_25[i];
		 count++;
	}
	for(int i=0;i<N_28;i++){
         R_EX_ALL[count] = R_EX_28[i];
         R_EX_Err_ALL[count] = R_EX_Err_28[i];
         xbj_ALL[count] = xbj_28[i];
		 count++;
	}

   	TCanvas *c5 = new TCanvas("c5","c5",1200,700);
	c5->cd();
	TH2F *h5 = new TH2F("h5","",300,0.7,XMax,300,RMIN,RMAX);
	h5->SetStats(kFALSE);
	h5->SetXTitle("x_{bj}");
	h5->GetXaxis()->CenterTitle(1);
	h5->GetXaxis()->SetTitleFont(32);
	h5->GetXaxis()->SetTitleSize(0.08);
	h5->GetXaxis()->SetTitleOffset(0.7);
	h5->GetYaxis()->SetTitle(Form("#frac{#sigma_{%s}}{#sigma_{%s}} #times #frac{%d}{%d}",Target0.Data(),Target.Data(),A2,A1));
	h5->GetYaxis()->CenterTitle(1);
	h5->GetYaxis()->SetTitleFont(32);
	h5->GetYaxis()->SetTitleSize(0.07);
	h5->GetYaxis()->SetTitleOffset(0.6);
	h5->GetYaxis()->CenterTitle();
	gPad->SetLogy(0);
	h5->Draw();

	TGraphErrors *exall = new TGraphErrors(N_ALL,xbj_ALL, R_EX_ALL, Zero,R_EX_Err_ALL);
	exall->SetMarkerStyle(24);
	exall->SetMarkerSize(1.3);
	exall->SetMarkerColor(2);
	exall->SetLineColor(2);
	c5->cd();	exall->Draw("P");
	exall->Fit("pol0","","",1.5, 1.8);

	if(bMC){	
		TGraphErrors *mc21 = new TGraphErrors(N_21,xbj_21, R_MC_21, Zero, Zero);
		mc21->SetLineWidth(2);
		mc21->SetLineColor(2);
		c3->cd(); mc21->Draw("l");
		c4->cd(); mc21->Draw("l");
		l1->AddEntry(mc21, Form("R_{MC} at %d^{0}",21),"l");

		TGraphErrors *mc23 = new TGraphErrors(N_23,xbj_23, R_MC_23, Zero, Zero);
		mc23->SetLineWidth(2);
		mc23->SetLineColor(4);
		c3->cd();	mc23->Draw("l");
		c4->cd();	mc23->Draw("l");
		l1->AddEntry(mc23, Form("R_{MC} at %d^{0}",23),"l");

		TGraphErrors *mcQE = new TGraphErrors(N_QE,xbj_QE, R_MC_QE, Zero, Zero);
		mcQE->SetLineWidth(2);
		mcQE->SetLineColor(6);
		c3->cd();	mcQE->Draw("l");
		c4->cd();	mcQE->Draw("l");
		l1->AddEntry(mcQE, Form("R_{MC} at %d^{0} (QE) ",25),"l");

		TGraphErrors *mc25 = new TGraphErrors(N_25,xbj_25, R_MC_25, Zero, Zero);
		mc25->SetLineWidth(2);
		mc25->SetLineColor(6);
		c3->cd();	mc25->Draw("l");
		c4->cd();	mc25->Draw("l");
		l1->AddEntry(mc25, Form("R_{MC} at %d^{0}",25),"l");

		TGraphErrors *mc28 = new TGraphErrors(N_28,xbj_28, R_MC_28, Zero, Zero);
		mc28->SetLineWidth(2);
		mc28->SetLineColor(3);
		if(Target=="C12"&&Target=="He3"){
			c3->cd();	mc28->Draw("l");
			c4->cd();	mc28->Draw("l");
			l1->AddEntry(mc28, Form("R_{MC} at %d^{0}",28),"l");
		}
	}

	c3->cd(); l1->Draw();
	c4->cd(); l1->Draw();
	if(bMC){
		c3->Print(Form("%s_%s_Y_Ratio_MC_%s_L.png",Target0.Data(),Target.Data(),Com0.Data()));
		c3->Print(Form("%s_%s_Y_Ratio_MC_%s_L.pdf",Target0.Data(),Target.Data(),Com0.Data()));
	}
	else{
		c3->Print(Form("%s_%s_Y_Ratio_%s_L.png",Target0.Data(),Target.Data(),Com0.Data()));
		c3->Print(Form("%s_%s_Y_Ratio_%s_L.pdf",Target0.Data(),Target.Data(),Com0.Data()));
	}
	if(bMC){
		c4->Print(Form("%s_%s_Y_Ratio_MC_Zoom_%s_L.png",Target0.Data(),Target.Data(),Com0.Data()));
		c4->Print(Form("%s_%s_Y_Ratio_MC_Zoom_%s_L.pdf",Target0.Data(),Target.Data(),Com0.Data()));
	}
	else{
		c4->Print(Form("%s_%s_Y_Ratio_Zoom_%s.png",Target0.Data(),Target.Data(),Com0.Data()));
		c4->Print(Form("%s_%s_Y_Ratio_Zoom_%s.pdf",Target0.Data(),Target.Data(),Com0.Data()));
	}

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
void Merge(const TString& aInputFileName, double* aXbj, double* aR_EX, double* aR_EX_Err, double* aR_MC, int* aN){

	/*Pre-define{{{*/	
	int aline = 0;
	TString aInputFilePath = "/work/halla/e08014/disk1/yez/Ep_Cross_Section/results/";

	double scale = 1.0;
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
	Int_t Count_T1L1 =0;
	Target = gGet_Target_Name(aName_T1L1.Data());
	if(Target=="He4")
		scale = He4_Norm;
	else
		scale = 1.0;

	Kin = gGet_Kin_Name(aName_T1L1.Data());

	TString infile_T1L1 = aInputFilePath+Target+"/"+Com0+"/"+aName_T1L1 +"_"+ Com0;
	cerr<<"==== Reading 1st file from "<<infile_T1L1<<endl;
	ifstream run_T1L1; run_T1L1.open(infile_T1L1);

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

			Y_EX_T1L1[Count_T1L1] = Yield_EX * scale; 
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

	Target = gGet_Target_Name(aName_T1L2.Data());
	if(Target=="He4")
		scale = He4_Norm;
	else
		scale = 1.0;
	Kin = gGet_Kin_Name(aName_T1L2.Data());
		
	TString infile_T1L2 = aInputFilePath+Target+"/"+Com0+"/"+aName_T1L2 +"_"+ Com0;
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

				Y_EX_T1L2[Count_T1L2] = Yield_EX * scale; 
				Y_MC_T1L2[Count_T1L2] = Yield_MC ; 
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

	Target = gGet_Target_Name(aName_T1R1.Data());
		if(Target=="He4")
		scale = He4_Norm;
		else
			scale = 1.0;
	Kin = gGet_Kin_Name(aName_T1R1.Data());
		
	TString infile_T1R1 = aInputFilePath+Target+"/"+Com0+"/"+aName_T1R1 +"_"+ Com0;
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

				Y_EX_T1R1[Count_T1R1] = Yield_EX * scale; 
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

	Target = gGet_Target_Name(aName_T1R2.Data());
	if(Target=="He4")
		scale = He4_Norm;
	else
		scale = 1.0;
	Kin = gGet_Kin_Name(aName_T1R2.Data());

	TString infile_T1R2 = aInputFilePath+Target+"/"+Com0+"/"+aName_T1R2 +"_"+ Com0;
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

				Y_EX_T1R2[Count_T1R2] = Yield_EX * scale; 
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
	Int_t Count_T2L1 =0;
	Target = gGet_Target_Name(aName_T2L1.Data());
	if(Target=="He4")
		scale = He4_Norm;
	else
		scale = 1.0;
	Kin = gGet_Kin_Name(aName_T2L1.Data());

	TString infile_T2L1 = aInputFilePath+Target+"/"+Com+"/"+aName_T2L1 +"_"+ Com;
	cerr<<"==== Reading 1st file from "<<infile_T2L1<<endl;
	ifstream run_T2L1; run_T2L1.open(infile_T2L1);

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

			Y_EX_T2L1[Count_T2L1] = Yield_EX * scale; 
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

	Target = gGet_Target_Name(aName_T2L2.Data());
	if(Target=="He4")
		scale = He4_Norm;
	else
		scale = 1.0;
	Kin = gGet_Kin_Name(aName_T2L2.Data());

	TString infile_T2L2 = aInputFilePath+Target+"/"+Com+"/"+aName_T2L2 +"_"+ Com;
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
	
				Y_EX_T2L2[Count_T2L2] = Yield_EX * scale; 
				Y_MC_T2L2[Count_T2L2] = Yield_MC ; 
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

	Target = gGet_Target_Name(aName_T2R1.Data());
	if(Target=="He4")
		scale = He4_Norm;
	else
		scale = 1.0;

	Kin = gGet_Kin_Name(aName_T2R1.Data());

	TString infile_T2R1 = aInputFilePath+Target+"/"+Com+"/"+aName_T2R1 +"_"+ Com;
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

				Y_EX_T2R1[Count_T2R1] = Yield_EX * scale; 
				Y_MC_T2R1[Count_T2R1] = Yield_MC ; 
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

	Target = gGet_Target_Name(aName_T2R2.Data());
	if(Target=="He4")
		scale = He4_Norm;
	else
		scale = 1.0;
	Kin = gGet_Kin_Name(aName_T2R2.Data());

	TString infile_T2R2 = aInputFilePath+Target+"/"+Com+"/"+aName_T2R2 +"_"+ Com;
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

				Y_EX_T2R2[Count_T2R2] = Yield_EX * scale; 
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
	}
	aN[0] = Count;
	outfile.close();
	/*}}}*/
}
/*}}}*/

/*E02-019 Data{{{*/
void E02019(const TString& aInputFileName, const double aTheta,double* aXbj, double* aR_EX, double* aR_EX_Err, int* aN){
	int aline = 0;
	TString aInputFilePath = "/work/halla/e08014/disk1/yez/Ep_Cross_Section/results/Ratio/";

	ofstream outfile(Form("E02019_XS_%s_%s.dat",Target0.Data(),Target.Data()));
	outfile << Form("%7s %14s %14s %14s %14s %14s %14s", "Xbj","Ratio","XS1","XS2","Ratio_Err","XS_Err1","XS_Err2" )<<endl;

	const TString aName_T1 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	const TString aName_T2 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	const double factor = gGet_InputFile_Var(aInputFileName, aline).Atof(); 
	cerr<<Form("==== Reading E02-019 Data from %s and %s ", aName_T1.Data(), aName_T2.Data())<<endl;

	TString Ascii_File1 = Form("%s/e02019/%s",aInputFilePath.Data(),aName_T1.Data());
	cerr<<Form(" === Reading E02019 data from %s", Ascii_File1.Data())<<endl;
	ifstream infile1(Ascii_File1);
	TString Ascii_File2 = Form("%s/e02019/%s",aInputFilePath.Data(),aName_T2.Data());
	cerr<<Form(" === Reading E02019 data from %s", Ascii_File2.Data())<<endl;
	ifstream infile2(Ascii_File2);

	Double_t xbj, Q2, nA, nTheta, nP0, nEp, nXS_Rad, nXS_Born,nXS_CC,nFrac_Err;  
	const double nE0 = 5.776;
	/*First Target:{{{*/
	vector<double> vXbj1, vXS1, vXS_Err1,vP1,vTheta1;
	int N1=0;
	while(!(infile1.eof())){
		infile1 >> nA >> nTheta >> nP0 >> nEp >> nXS_Rad >> nXS_Born >> nXS_CC >> nFrac_Err;
		if(abs(nTheta-aTheta)<1.0|| (aTheta-0.0)<1){
		//if(abs(nTheta-aTheta)<1.0|| (nTheta)<aTheta){
			Q2 = 4.0 * nE0 * nEp * pow(sin(nTheta*3.1415926/180./2.0),2);
			xbj = Q2/2.0/0.938272/(nE0-nEp);
			//	cerr<<Form("=== E0=%f, Ep=%f, Theta=%f",nE0,nEp,nTheta)<<endl;
			N1++;
			vXbj1.push_back(xbj);
			vXS1.push_back(nXS_Born);
			vXS_Err1.push_back(nXS_Born*nFrac_Err);
			vTheta1.push_back(nTheta);
			vP1.push_back(nP0);
		}
	}
	/*}}}*/

	/* Second  Target:{{{*/
	vector<double> vXbj2, vXS2, vXS_Err2,vP2,vTheta2;
	int N2=0;
	while(!(infile2.eof())){
		infile2 >> nA >> nTheta >> nP0 >> nEp >> nXS_Rad >> nXS_Born >> nXS_CC >> nFrac_Err;

		if(abs(nTheta-aTheta)<1.0||(aTheta-0.0)<1){
		//if(abs(nTheta-aTheta)<1.0|| (nTheta)<aTheta){
			Q2 = 4.0 * nE0 * nEp * pow(sin(nTheta*3.1415926/180./2.0),2);
			xbj = Q2/2.0/0.938272/(nE0-nEp);
			//	cerr<<Form("===  E0=%f, Ep=%f, Theta=%f",nE0,nEp,nTheta)<<endl;
			N2++;
			vXbj2.push_back(xbj);
			vXS2.push_back(nXS_Born);
			vXS_Err2.push_back(nXS_Born*nFrac_Err);
			vTheta2.push_back(nTheta);
			vP2.push_back(nP0);
		}
	}
	/*}}}*/

	/*Ratio{{{*/  
	Int_t Count=0;
	Double_t x1,x2,ratio,err;
	vector<Double_t> vXbj, vXbj_Err,vR, vR_Err;

	for(int i=0;i<N1;i++){
		x1 = vXbj1[i]; 
		for(int j=0;j<N2;j++){
			x2 = vXbj2[j];
			if(abs(x1-x2)<0.0001 && abs(vP1[i]-vP2[j])<0.1&& abs(vTheta1[i]-vTheta2[j])<1.0){
				ratio = factor*vXS1[i]/vXS2[j];
				err = ratio * sqrt( pow(vXS_Err2[j]/vXS2[j],2) + pow(vXS_Err1[i]/vXS1[i],2) );
				if(ratio>1e-30 && ratio<100.&&(err/ratio)<0.5){
					vXbj.push_back(x1);
					vXbj_Err.push_back(0.0);
					vR.push_back(ratio);
					vR_Err.push_back(err);
					cerr <<Form("-----> Found one: xbj1=%f, xbj2=%f, XS1=%f, XS2=%f,Ratio=%f ",x1,x2,vXS1[i],vXS2[j],factor*vXS1[i]/vXS2[j])<<endl;
					Count++;
					outfile << Form("%7.5f %14e %14e %14e %14e %14e %14e", x1,ratio, vXS1[i], vXS2[j],err, vXS_Err1[i], vXS_Err2[j] )<<endl;
				}
			}
		}
	}

	//Return the value
	for(int k=0; k<Count; k++){
		aXbj[aN[0]+k] = vXbj[k];
		aR_EX[aN[0]+k] = vR[k];
		aR_EX_Err[aN[0]+k] = vR_Err[k];
	}
	aN[0] += Count;
	cerr<<Form("******************** Count = %d, Factor =%f ",Count, factor)<<endl;
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

/*E02-019 A/He3{{{*/
void E02019_A_He3(const TString& aTarget, const double aTheta,double* aXbj, double* aR_EX, double* aR_EX_Err_Up,double* aR_EX_Err_Down, int* aN){
	int aline = 0;
	TString aInputFilePath = "/work/halla/e08014/disk1/yez/Ep_Cross_Section/results/Ratio_L/";
	TString aInputFile = aInputFilePath + "Nadia_A_He3.dat";
	ifstream outfile(aInputFile);
    Double_t nTheta,nXbj,nRatio,nErr1,nErr2;
	Double_t aErr_Up,aErr_Down;
    int aA,aZ,nA;
    gGetTargetAZ(aTarget.Data(), &aA,&aZ);
    
	int aNNN = 0;
	while(!outfile.eof()){
		outfile >> nA >> nTheta >> nXbj >> nRatio >> nErr1 >> nErr2;
		if(abs(nA-aA)<1.0&&(abs(aTheta-nTheta)<1.0 ||(aTheta-0.0)<1.0)&&nRatio>=1e-33){
			aXbj[aNNN] = nXbj;
			aR_EX[aNNN] = nRatio;
			aErr_Up = abs(nErr1 - nRatio);
			aErr_Down = abs(nRatio - nErr2);	
			aR_EX_Err_Up[aNNN] = aErr_Up;
			aR_EX_Err_Down[aNNN] = aErr_Down;
			aNNN++;
		}
	}
	aN[0]=aNNN;
	outfile.close();
}
/*}}}*/

/*E02-019 A/He3{{{*/
void E02019_He4_He3(const TString& aTarget, const double aTheta,double* aXbj, double* aR_EX, double* aR_EX_Err_Up,double* aR_EX_Err_Down, int* aN){
	int aline = 0;
	TString aInputFilePath = "/work/halla/e08014/disk1/yez/Ep_Cross_Section/results/Ratio_L/";
	TString aInputFile = aInputFilePath + "Nadia_A_He3.dat";
	ifstream outfile(aInputFile);
    Double_t nTheta,nXbj,nRatio,nErr,nErr1,nErr2;
	Double_t aErr_Up,aErr_Down;
    int aA,aZ,nA;
    gGetTargetAZ(aTarget.Data(), &aA,&aZ);
    double aXbj_Max = 1.55;
	int aNNN = 0;
	while(!outfile.eof()){
		outfile >> nA >> nTheta >> nXbj >> nRatio >> nErr1 >> nErr2;
		if(abs(nA-aA)<1.0&&(abs(aTheta-nTheta)<1.0 ||(aTheta-0.0)<1.0)&&nRatio>=1e-33&&nXbj<=aXbj_Max){
			aXbj[aNNN] = nXbj;
			aR_EX[aNNN] = nRatio;
			aErr_Up = abs(nErr1 - nRatio);
			aErr_Down = abs(nRatio - nErr2);	
			aR_EX_Err_Up[aNNN] = aErr_Up;
			aR_EX_Err_Down[aNNN] = aErr_Down;
			aNNN++;
		}
	}

	TString aInputFile1 = aInputFilePath + "Nadia_He4_He3.dat";
	ifstream outfile1(aInputFile1);
	while(!outfile1.eof()){
		outfile1 >> nXbj >> nRatio >> nErr;
			aXbj[aNNN] = nXbj;
			aR_EX[aNNN] = 1./nRatio;
			aErr_Up = 1./(nRatio-nErr);
			aErr_Down = 1./(nRatio + nErr);	
			aR_EX_Err_Up[aNNN] = aErr_Up - 1./nRatio;
			aR_EX_Err_Down[aNNN] = 1./nRatio - aErr_Down;
			aNNN++;
	}
	aN[0]=aNNN;
	outfile.close();
	outfile1.close();
}
/*}}}*/

/*CLAS_He4_He3{{{*/
void CLAS_He4_He3( double* aXbj_CS, double* aR_EX_CS, double* aR_EX_Err_CS, int *aN_CS){
	TString aInputFilePath = "/work/halla/e08014/disk1/yez/Ep_Cross_Section/results/Ratio_L/";
	TString aInputFile = aInputFilePath + "CLAS_He4_He3.dat";
	ifstream infile(aInputFile);

	double cXbj, cR, cR_Err;
	int aN = 0;

	while(!(infile.eof())){
      infile >> cXbj >> cR >> cR_Err;	
      aXbj_CS[aN] = cXbj;
	  aR_EX_CS[aN] = cR/1.15;
	  aR_EX_Err_CS[aN] = cR_Err/1.15;
	  aN++;	
	}
	aN_CS[0] = aN;
	infile.close();
}
/*}}}*/
