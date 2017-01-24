#include "/work/halla/e08014/disk1/yez/Xbj_Cross_Section/SRC/XGT2_Main.h"
#include "/work/halla/e08014/disk1/yez/Xbj_Cross_Section/SRC/XGT2_XEMC.h"
#include "TLegend.h"
#include "TMath.h"
#include "TGraph.h"
/*Get_Target_Name{{{*/
TString gGet_Target_Name(const TString& aName){
	TString aTarget_Name="MM";

	if(aName.Contains("H2"))
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
	else if(aKin.Contains("5.0")&&!aKin.Contains("5.05"))
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

void Merge(const TString& aInputFileName, double* aXbj, double* aEp, double* aXS_EX, double* aXS_EX_Err, double* aXS_MC, int* aN);
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
TString Com;
//Cuts on the edge of Ep
// Cut_L_Min Cut_L_Max Cut_R_Min Cut_R_Max
const double Ep_Cut_21_1[4] = {2.80, 3.00,2.80,3.00}; 
const double Ep_Cut_21_2[4] = {2.94, 3.18,2.94,3.18}; 
const double Ep_Cut_23_1[4] = {2.74, 2.94,2.74,2.94}; 
const double Ep_Cut_23_2[4] = {2.92, 3.16,2.92,3.16};
const double Ep_Cut_QE_1[4] = {2.42, 2.58,2.42,2.58}; 
const double Ep_Cut_QE_2[4] = {2.56, 2.74,2.56,2.74}; 
const double Ep_Cut_25_1[4] = {2.70, 2.88,2.70,2.88};
const double Ep_Cut_25_2[4] = {2.88, 3.12,2.88,3.12}; 
const double Ep_Cut_28_1[4] = {2.74, 2.98,2.74,2.98}; 
//int Ratio()
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

	TString Target; cerr<<" Which Target? (H2,He3,He4,C12,Ca40,Ca48) "; cin >> Target;

	cerr<<"--- Comment = "; cin>>Com;
	//bool bMC = kTRUE;
	//bool bMC = kFALSE;
	//double Theta =21;
	//double Zero[100] = {100*0.0};
    double E0 = 3.356;

	ofstream outfile(Form("XS_%s_%s.dat",Target.Data(),Com.Data()));
	outfile << Form("%4s %7s %4s %7s %7s %7s %16s %16s %16s", 
			"Target",  "Angle","E0","Ep","Nu","xbj","XS","XS_Err","XS_MC")
		<< endl;
	
	/*Merge{{{*/

	/*Degree 21 {{{*/
	double xbj_21[100], Ep_21[100], Nu_21[100], XS_EX_21[100], XS_EX_Err_21[100], XS_MC_21[100];
	int N_21 = 0;
	if(Target=="H2" ||Target=="He3" ||  Target=="He4" || Target=="C12"||Target=="Ca40"||Target=="Ca48"){
		Merge(Form("R21_%s.dat",Target.Data()), xbj_21, Ep_21, XS_EX_21, XS_EX_Err_21, XS_MC_21, &N_21);
		cerr<<"=== 21: Total Points = "<<N_21<<endl;
		for(int i=0;i<N_21;i++){
			Nu_21[i] = E0 - Ep_21[i];
			cerr    << Form("%4s %7.3f %7.3f %7.3f %7.3f %7.3f %16e %16e %16e",
				   	Target.Data(),21.00, E0, Ep_21[i], Nu_21[i], xbj_21[i],XS_EX_21[i],XS_EX_Err_21[i],XS_MC_21[i] )
				<< endl;
			outfile << Form("%4s %7.3f %7.3f %7.3f %7.3f %7.3f %16e %16e %16e", 
					Target.Data(),21.00, E0, Ep_21[i], Nu_21[i], xbj_21[i],XS_EX_21[i],XS_EX_Err_21[i],XS_MC_21[i] )
				<< endl;
		}
	}
	double Q2_21_D = 4.0*E0*Ep_Cut_21_1[0]*pow(TMath::Sin(21.*3.1415926/180.0*0.5),2);
	double Q2_21_U = 4.0*E0*Ep_Cut_21_2[3]*pow(TMath::Sin(21.*3.1415926/180.0*0.5),2);
	/*}}}*/

	/*Degree 23 {{{*/
	double xbj_23[100], Nu_23[100], Ep_23[100], XS_EX_23[100], XS_EX_Err_23[100], XS_MC_23[100];
	int N_23 = 0;
	if(Target=="He3" ||  Target=="He4" || Target=="C12"||Target=="Ca40"||Target=="Ca48"){
		Merge(Form("R23_%s.dat",Target.Data()), xbj_23, Ep_23, XS_EX_23, XS_EX_Err_23, XS_MC_23, &N_23);
		cerr<<"=== 23: Total Points = "<<N_23<<endl;
		for(int i=0;i<N_23;i++){
			Nu_23[i] = E0 - Ep_23[i];
			cerr    << Form("%4s %7.3f %7.3f %7.3f %7.3f %7.3f %16e %16e %16e", 
					Target.Data(),23.00, E0, Ep_23[i], Nu_23[i], xbj_23[i],XS_EX_23[i],XS_EX_Err_23[i],XS_MC_23[i] )
				<< endl;
			outfile << Form("%4s %7.3f %7.3f %7.3f %7.3f %7.3f %16e %16e %16e", 
					Target.Data(),23.00, E0, Ep_23[i], Nu_23[i], xbj_23[i],XS_EX_23[i],XS_EX_Err_23[i],XS_MC_23[i] )
				<< endl;
		}
	}
	double Q2_23_D = 4.0*E0*Ep_Cut_23_1[0]*pow(TMath::Sin(23.*3.1415926/180.0*0.5),2);
	double Q2_23_U = 4.0*E0*Ep_Cut_23_2[3]*pow(TMath::Sin(23.*3.1415926/180.0*0.5),2);
		/*}}}*/

	/*QE{{{*/
	double xbj_QE[100], Nu_QE[100], Ep_QE[100], XS_EX_QE[100], XS_EX_Err_QE[100], XS_MC_QE[100];
	int N_QE = 0;
	if(Target=="He3" ||  Target=="He4" || Target=="C12"||Target=="Ca40"||Target=="Ca48"){
		Merge(Form("RQE_%s.dat",Target.Data()), xbj_QE, Ep_QE, XS_EX_QE, XS_EX_Err_QE, XS_MC_QE, &N_QE);
		cerr<<"=== QE: Total Points = "<<N_QE<<endl;
		for(int i=0;i<N_QE;i++){
			Nu_QE[i] = E0 - Ep_QE[i];
			cerr    << Form("%4s %7.3f %7.3f %7.3f %7.3f %7.3f %16e %16e %16e", Target.Data(),25.00, E0, Ep_QE[i], Nu_QE[i], xbj_QE[i],XS_EX_QE[i],XS_EX_Err_QE[i],XS_MC_QE[i] )
				<< endl;

			outfile << Form("%4s %7.3f %7.3f %7.3f %7.3f %7.3f %16e %16e %16e", Target.Data(),25.00, E0, Ep_QE[i], Nu_QE[i], xbj_QE[i],XS_EX_QE[i],XS_EX_Err_QE[i],XS_MC_QE[i] )
				<< endl;
		}
	}
	double Q2_25_D = 4.0*E0*Ep_Cut_QE_1[0]*pow(TMath::Sin(25.*3.1415926/180.0*0.5),2);
	/*}}}*/

	/*Degree 25 {{{*/
	double xbj_25[100], Ep_25[100], Nu_25[100], XS_EX_25[100], XS_EX_Err_25[100], XS_MC_25[100];
	int N_25 = 0;
	if(Target=="H2" ||Target=="He3" ||  Target=="He4" || Target=="C12"||Target=="Ca40"||Target=="Ca48"){
		Merge(Form("R25_%s.dat",Target.Data()), xbj_25, Ep_25, XS_EX_25, XS_EX_Err_25, XS_MC_25, &N_25);
		cerr<<"=== 25: Total Points = "<<N_25<<endl;
		for(int i=0;i<N_25;i++){
			Nu_25[i] = E0 - Ep_25[i];
			cerr    << Form("%4s %7.3f %7.3f %7.3f %7.3f %7.3f %16e %16e %16e", 
					Target.Data(),25.00, E0, Ep_25[i], Nu_25[i], xbj_25[i],XS_EX_25[i],XS_EX_Err_25[i],XS_MC_25[i] )
				<< endl;
			outfile << Form("%4s %7.3f %7.3f %7.3f %7.3f %7.3f %16e %16e %16e", 
					Target.Data(),25.00, E0, Ep_25[i], Nu_25[i], xbj_25[i],XS_EX_25[i],XS_EX_Err_25[i],XS_MC_25[i] )
				<< endl;
		}
	}
	double Q2_25_U = 4.0*E0*Ep_Cut_25_2[3]*pow(TMath::Sin(25.*3.1415926/180.0*0.5),2);
	/*}}}*/

	/*Degree 28 {{{*/
	double xbj_28[100], Nu_28[100], Ep_28[100], XS_EX_28[100], XS_EX_Err_28[100], XS_MC_28[100];
	int N_28 = 0;
	if(Target=="He3" || Target=="C12"){
		Merge(Form("R28_%s.dat",Target.Data()), xbj_28, Ep_28, XS_EX_28, XS_EX_Err_28, XS_MC_28, &N_28);
		cerr<<"=== 28: Total Points = "<<N_28<<endl;
		for(int i=0;i<N_28;i++){
			Nu_28[i] = E0 - Ep_28[i];
			cerr    << Form("%4s %7.3f %7.3f %7.3f %7.3f %7.3f %16e %16e %16e", 
					Target.Data(),28.00, E0, Ep_28[i], Nu_28[i], xbj_28[i],XS_EX_28[i],XS_EX_Err_28[i],XS_MC_28[i] )
				<< endl;
			outfile << Form("%4s %7.3f %7.3f %7.3f %7.3f %7.3f %16e %16e %16e", 
					Target.Data(),28.00, E0, Ep_28[i], Nu_28[i], xbj_28[i],XS_EX_28[i],XS_EX_Err_28[i],XS_MC_28[i] )
				<< endl;
		}
	}
	double Q2_28_D = 4.0*E0*Ep_Cut_28_1[0]*pow(TMath::Sin(28.*3.1415926/180.0*0.5),2);
	double Q2_28_U = 4.0*E0*Ep_Cut_28_1[3]*pow(TMath::Sin(28.*3.1415926/180.0*0.5),2);
	/*}}}*/

	/*}}}*/

    /*Plot{{{*/
    /*Range{{{*/
	double Min = 1000, Max =-1000;
	for(int i=0; i<N_21;i++){
		Min = Min<XS_EX_21[i]?Min:XS_EX_21[i];
		Max = Max>XS_EX_21[i]?Max:XS_EX_21[i];
	}
	for(int i=0; i<N_23;i++){
		Min = Min<XS_EX_23[i]?Min:XS_EX_23[i];
		Max = Max>XS_EX_23[i]?Max:XS_EX_23[i];
	}
	for(int i=0; i<N_QE;i++){
		Min = Min<XS_EX_QE[i]?Min:XS_EX_QE[i];
		Max = Max>XS_EX_QE[i]?Max:XS_EX_QE[i];
	}
	for(int i=0; i<N_25;i++){
		Min = Min<XS_EX_25[i]?Min:XS_EX_25[i];
		Max = Max>XS_EX_25[i]?Max:XS_EX_25[i];
	}
	for(int i=0; i<N_28;i++){
		Min = Min<XS_EX_28[i]?Min:XS_EX_28[i];
		Max = Max>XS_EX_28[i]?Max:XS_EX_28[i];
	}
/*}}}*/

	/*Xbj{{{*/
	double xbj_err[200];
	for(int i=0;i<200;i++)
       xbj_err[200]=0;

	Double_t Xbj_Max = 2.2;
	Double_t Xbj_Min = 0.5;
	if(Target=="H2")
		Xbj_Max = 2.2;
	else if(Target=="He3"){
		Xbj_Max = 3.5;
		Max*=3.0;
	}
	else if(Target=="He4")
		Xbj_Max = 4.3;
	else
		Xbj_Max = 5.3;

	TCanvas *c1 = new TCanvas("c1","c1",1200,800);
	TH2F* xs = new TH2F("xs","",1000,Xbj_Min,Xbj_Max,1000,Min/2.1,Max*2.1);
	xs->SetStats(kFALSE);
	xs->SetXTitle("x_{bj}");
	xs->GetXaxis()->CenterTitle(1);
	xs->GetXaxis()->SetTitleFont(32);
	xs->GetXaxis()->SetTitleSize(0.06);
	xs->SetYTitle("#frac{d#sigma}{dE'd#Omega} (nb/MeV/sr)");
	xs->GetYaxis()->CenterTitle(1);
	xs->GetYaxis()->SetTitleFont(32);
	xs->GetYaxis()->SetTitleSize(0.04);
	xs->Draw(); 
	//gPad->SetLogy(1);
	gPad->SetLogy(1);

	TLegend *l1 = new TLegend(0.60,0.59,0.75, 0.94,Form("%s Target",Target.Data()));
	l1->SetBorderSize(0);
	l1->SetTextSize(0.03);
	l1->SetTextFont(32);

	TGraphErrors *pl_21 = new TGraphErrors(N_21,xbj_21,XS_EX_21,xbj_err,XS_EX_Err_21);
	pl_21-> SetMarkerSize(1.0);
	pl_21-> SetMarkerStyle(24);
	pl_21-> SetMarkerColor(2);
	pl_21-> GetXaxis()->SetTitle("x_{bj}");
	pl_21-> GetYaxis()->SetTitle("#sigma");
	pl_21->Draw("P");
	l1->AddEntry(pl_21, Form("21^{o}: %4.3f<Q^{2}<%4.3f (GeV/c^{2})",Q2_21_D,Q2_21_U),"p");

	TGraph *p1_21_mc = new TGraphErrors(N_21,xbj_21,XS_MC_21);
	p1_21_mc-> SetLineWidth(2);
	p1_21_mc-> SetLineStyle(2);
	//p1_21_mc-> SetLineColor(2);
	p1_21_mc-> SetLineColor(1);
	p1_21_mc-> GetXaxis()->SetTitle("x_{bj}");
	p1_21_mc-> GetYaxis()->SetTitle("#sigma");
	p1_21_mc->Draw("L");
	
	TGraphErrors *pl_23 = new TGraphErrors(N_23,xbj_23,XS_EX_23,xbj_err,XS_EX_Err_23);
	if(!(Target=="H2")){
		pl_23-> SetMarkerSize(1.0);
		pl_23-> SetMarkerStyle(25);
		pl_23-> SetMarkerColor(3);
		pl_23-> GetXaxis()->SetTitle("x_{bj}");
		pl_23-> GetYaxis()->SetTitle("#sigma");
		pl_23->Draw("P");
		l1->AddEntry(pl_23, Form("23^{o}: %4.3f<Q^{2}<%4.3f (GeV/c^{2})",Q2_23_D,Q2_23_U),"p");
		
		TGraph *p1_23_mc = new TGraphErrors(N_23,xbj_23,XS_MC_23);
		p1_23_mc-> SetLineWidth(2);
		p1_23_mc-> SetLineStyle(2);
		//p1_23_mc-> SetLineColor(3);
		p1_23_mc-> SetLineColor(1);
		p1_23_mc-> GetXaxis()->SetTitle("x_{bj}");
		p1_23_mc-> GetYaxis()->SetTitle("#sigma");
		p1_23_mc->Draw("L");
	}

	TGraphErrors *pl_QE = new TGraphErrors(N_QE,xbj_QE,XS_EX_QE,xbj_err,XS_EX_Err_QE);
	pl_QE-> SetMarkerSize(1.0);
	pl_QE-> SetMarkerStyle(26);
	pl_QE-> SetMarkerColor(4);
	pl_QE-> GetXaxis()->SetTitle("x_{bj}");
	pl_QE-> GetYaxis()->SetTitle("#sigma");
	pl_QE->Draw("P");

	TGraph *p1_QE_mc = new TGraphErrors(N_QE,xbj_QE,XS_MC_QE);
	p1_QE_mc-> SetLineWidth(2);
	p1_QE_mc-> SetLineStyle(2);
	//p1_QE_mc-> SetLineColor(4);
	p1_QE_mc-> SetLineColor(1);
	p1_QE_mc-> GetXaxis()->SetTitle("x_{bj}");
	p1_QE_mc-> GetYaxis()->SetTitle("#sigma");
	p1_QE_mc->Draw("L");

	TGraphErrors *pl_25 = new TGraphErrors(N_25,xbj_25,XS_EX_25,xbj_err,XS_EX_Err_25);
	pl_25-> SetMarkerSize(1.0);
	pl_25-> SetMarkerStyle(26);
	pl_25-> SetMarkerColor(4);
	pl_25-> GetXaxis()->SetTitle("x_{bj}");
	pl_25-> GetYaxis()->SetTitle("#sigma");
	pl_25->Draw("P");
	l1->AddEntry(pl_25, Form("25^{o}: %4.3f<Q^{2}<%4.3f (GeV/c^{2})",Q2_25_D,Q2_25_U),"p");

	TGraph *p1_25_mc = new TGraphErrors(N_25,xbj_25,XS_MC_25);
	p1_25_mc-> SetLineWidth(2);
	p1_25_mc-> SetLineStyle(2);
	//p1_25_mc-> SetLineColor(4);
	p1_25_mc-> SetLineColor(1);
	p1_25_mc-> GetXaxis()->SetTitle("x_{bj}");
	p1_25_mc-> GetYaxis()->SetTitle("#sigma");
	p1_25_mc->Draw("L");
	
	TGraphErrors *pl_28 = new TGraphErrors(N_28,xbj_28,XS_EX_28,xbj_err,XS_EX_Err_28);
	if(Target=="He3" || Target=="C12"){
		pl_28-> SetMarkerSize(1.0);
		pl_28-> SetMarkerStyle(28);
		pl_28-> SetMarkerColor(6);
		pl_28-> GetXaxis()->SetTitle("x_{bj}");
		pl_28-> GetYaxis()->SetTitle("#sigma");
		pl_28->Draw("P");
		l1->AddEntry(pl_28, Form("28^{o}: %4.3f<Q^{2}<%4.3f (GeV/c^{2})",Q2_28_D,Q2_28_U),"p");

		TGraph *p1_28_mc = new TGraphErrors(N_28,xbj_28,XS_MC_28);
		p1_28_mc-> SetLineWidth(2);
		p1_28_mc-> SetLineStyle(2);
		//p1_28_mc-> SetLineColor(6);
		p1_28_mc-> SetLineColor(1);
		p1_28_mc-> GetXaxis()->SetTitle("x_{bj}");
		p1_28_mc-> GetYaxis()->SetTitle("#sigma");
		p1_28_mc->Draw("L");
	}

	//l1->AddEntry(p1_21_mc, "Theory","l");
	l1->AddEntry(p1_21_mc, "Model","l");
	l1->Draw();
	c1->Print(Form("%s_XS_All_xbj_%s.pdf",Target.Data(),Com.Data()));
	c1->Print(Form("%s_XS_All_xbj_%s.png",Target.Data(),Com.Data()));
	/*}}}*/

	/*Ep{{{*/
	const double Ep_Min = 2.30;
	const double Ep_Max = 3.30;

	TCanvas *c2 = new TCanvas("c2","c2",1200,800);
	TH2F* ep = new TH2F("ep","",1000,E0-Ep_Max,E0-Ep_Min,1000,Min/2.1,Max*2.1);
	ep->SetStats(kFALSE);
	ep->SetXTitle("#nu");
	ep->GetXaxis()->CenterTitle(1);
	ep->GetXaxis()->SetTitleFont(32);
	ep->GetXaxis()->SetTitleSize(0.06);
	ep->SetYTitle("#frac{d#sigma}{dE'd#Omega} (nb/MeV/sr)");
	ep->GetYaxis()->CenterTitle(1);
	ep->GetYaxis()->SetTitleFont(32);
	ep->GetYaxis()->SetTitleSize(0.04);
	ep->Draw(); 
	//gPad->SetLogy(1);
	gPad->SetLogy(1);

	TLegend *l2 = new TLegend(0.60,0.19,0.75, 0.54,Form("%s Target",Target.Data()));
	l2->SetBorderSize(0);
	l2->SetTextSize(0.03);
	l2->SetTextFont(32);

	TGraphErrors *p2_21 = new TGraphErrors(N_21,Nu_21,XS_EX_21,xbj_err,XS_EX_Err_21);
	p2_21-> SetMarkerSize(1.0);
	p2_21-> SetMarkerStyle(24);
	p2_21-> SetMarkerColor(2);
	p2_21-> GetXaxis()->SetTitle("#nu");
	p2_21-> GetYaxis()->SetTitle("#sigma");
	p2_21->Draw("P");
	l2->AddEntry(p2_21, Form("21^{o}: %4.3f<Q^{2}<%4.3f (GeV/c^{2})",Q2_21_D,Q2_21_U),"p");

	TGraph *p2_21_mc = new TGraphErrors(N_21,Nu_21,XS_MC_21);
	p2_21_mc-> SetLineWidth(2);
	p2_21_mc-> SetLineStyle(2);
	//p2_21_mc-> SetLineColor(2);
	p2_21_mc-> SetLineColor(1);
	p2_21_mc-> GetXaxis()->SetTitle("#nu");
	p2_21_mc-> GetYaxis()->SetTitle("#sigma");
	p2_21_mc->Draw("L");
	
	TGraphErrors *p2_23 = new TGraphErrors(N_23,Nu_23,XS_EX_23,xbj_err,XS_EX_Err_23);
	if(!(Target=="H2")){
		p2_23-> SetMarkerSize(1.0);
		p2_23-> SetMarkerStyle(25);
		p2_23-> SetMarkerColor(3);
		p2_23-> GetXaxis()->SetTitle("#nu");
		p2_23-> GetYaxis()->SetTitle("#sigma");
		p2_23->Draw("P");
		l2->AddEntry(p2_23, Form("23^{o}: %4.3f<Q^{2}<%4.3f (GeV/c^{2})",Q2_23_D,Q2_23_U),"p");
	
		TGraph *p2_23_mc = new TGraphErrors(N_23,Nu_23,XS_MC_23);
		p2_23_mc-> SetLineWidth(2);
		p2_23_mc-> SetLineStyle(2);
		//p2_23_mc-> SetLineColor(3);
		p2_23_mc-> SetLineColor(1);
		p2_23_mc-> GetXaxis()->SetTitle("#nu");
		p2_23_mc-> GetYaxis()->SetTitle("#sigma");
		p2_23_mc->Draw("L");

}

	TGraphErrors *p2_QE = new TGraphErrors(N_QE,Nu_QE,XS_EX_QE,xbj_err,XS_EX_Err_QE);
	p2_QE-> SetMarkerSize(1.0);
	p2_QE-> SetMarkerStyle(26);
	p2_QE-> SetMarkerColor(4);
	p2_QE-> GetXaxis()->SetTitle("#nu");
	p2_QE-> GetYaxis()->SetTitle("#sigma");
	p2_QE->Draw("P");
	TGraphErrors *p2_25 = new TGraphErrors(N_25,Nu_25,XS_EX_25,xbj_err,XS_EX_Err_25);
	p2_25-> SetMarkerSize(1.0);
	p2_25-> SetMarkerStyle(26);
	p2_25-> SetMarkerColor(4);
	p2_25-> GetXaxis()->SetTitle("#nu}");
	p2_25-> GetYaxis()->SetTitle("#sigma");
	p2_25->Draw("P");
	l2->AddEntry(p2_25, Form("25^{o}: %4.3f<Q^{2}<%4.3f (GeV/c^{2})",Q2_25_D,Q2_25_U),"p");

	TGraph *p2_QE_mc = new TGraphErrors(N_QE,Nu_QE,XS_MC_QE);
	p2_QE_mc-> SetLineWidth(2);
	p2_QE_mc-> SetLineStyle(2);
	//p2_QE_mc-> SetLineColor(4);
	p2_QE_mc-> SetLineColor(1);
	p2_QE_mc-> GetXaxis()->SetTitle("#nu");
	p2_QE_mc-> GetYaxis()->SetTitle("#sigma");
	p2_QE_mc->Draw("L");

	TGraph *p2_25_mc = new TGraphErrors(N_25,Nu_25,XS_MC_25);
	p2_25_mc-> SetLineWidth(2);
	p2_25_mc-> SetLineStyle(2);
	p2_25_mc-> SetLineColor(1);
	//p2_25_mc-> SetLineColor(4);
	p2_25_mc-> GetXaxis()->SetTitle("#nu");
	p2_25_mc-> GetYaxis()->SetTitle("#sigma");
	p2_25_mc->Draw("L");

	TGraphErrors *p2_28 = new TGraphErrors(N_28,Nu_28,XS_EX_28,xbj_err,XS_EX_Err_28);
	if(Target=="He3" || Target=="C12"){
		p2_28-> SetMarkerSize(1.0);
		p2_28-> SetMarkerStyle(28);
		p2_28-> SetMarkerColor(6);
		p2_28-> GetXaxis()->SetTitle("#nu");
		p2_28-> GetYaxis()->SetTitle("#sigma");
		p2_28->Draw("P");
		l2->AddEntry(p2_28, Form("28^{o}: %4.3f<Q^{2}<%4.3f (GeV/c^{2})",Q2_28_D,Q2_28_U),"p");

		TGraph *p2_28_mc = new TGraphErrors(N_28,Nu_28,XS_MC_28);
		p2_28_mc-> SetLineWidth(2);
		p2_28_mc-> SetLineStyle(2);
		//p2_28_mc-> SetLineColor(6);
		p2_28_mc-> SetLineColor(1);
		p2_28_mc-> GetXaxis()->SetTitle("#nu");
		p2_28_mc-> GetYaxis()->SetTitle("#sigma");
		p2_28_mc->Draw("L");
	}

    //l2->AddEntry(p2_21_mc,"Theory","l");
    l2->AddEntry(p2_21_mc,"Model","l");

    l2->Draw();
    c2->Print(Form("%s_XS_All_%s.png",Target.Data(),Com.Data()));
	c2->Print(Form("%s_XS_All_%s.pdf",Target.Data(),Com.Data()));
 	/*}}}*/
	/*}}}*/

	outfile.close();

	return 0;
}
/*}}}*/

/*gGet_Ep_Cut{{{*/
void gGet_Ep_Cut(const TString &aKinName, const TString & aArm, double* aEp_Min, double* aEp_Max){

	if(aKinName=="Kin3.1"){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_21_1[0];		
			aEp_Max[0] = Ep_Cut_21_1[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_21_1[2];		
			aEp_Max[0] = Ep_Cut_21_1[3];		
		}
	}
	else if(aKinName=="Kin3.2"){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_21_2[0];		
			aEp_Max[0] = Ep_Cut_21_2[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_21_2[2];		
			aEp_Max[0] = Ep_Cut_21_2[3];		
		}
	}

	else if(aKinName=="Kin4.1"){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_23_1[0];		
			aEp_Max[0] = Ep_Cut_23_1[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_23_1[2];		
			aEp_Max[0] = Ep_Cut_23_1[3];		
		}
	}
	else if(aKinName=="Kin4.2"){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_23_2[0];		
			aEp_Max[0] = Ep_Cut_23_2[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_23_2[2];		
			aEp_Max[0] = Ep_Cut_23_2[3];		
		}
	}

	else if(aKinName=="Kin5.0"){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_QE_1[0];		
			aEp_Max[0] = Ep_Cut_QE_1[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_QE_1[2];		
			aEp_Max[0] = Ep_Cut_QE_1[3];		
		}
	}
	else if(aKinName=="Kin5.05"){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_QE_2[0];		
			aEp_Max[0] = Ep_Cut_QE_2[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_QE_2[2];		
			aEp_Max[0] = Ep_Cut_QE_2[3];		
		}
	}

	else if(aKinName=="Kin5.1"){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_25_1[0];		
			aEp_Max[0] = Ep_Cut_25_1[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_25_1[2];		
			aEp_Max[0] = Ep_Cut_25_1[3];		
		}
	}
	else if(aKinName=="Kin5.2"){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_25_2[0];		
			aEp_Max[0] = Ep_Cut_25_2[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_25_2[2];		
			aEp_Max[0] = Ep_Cut_25_2[3];		
		}
	}

	else if(aKinName=="Kin6.5"){
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

/*void XS_Misak(const TString& aTarget_Name, const double aTheta, const int aN, double* aEp, double* aXS_Model){{{*/
void XS_Misak(const TString& aTarget_Name, const double aTheta, const int aN, double* aEp, double* aXS_Model){

	TString aXS_Input_Name = " ";
	cerr<<Form("--- Working on %s at Theta = %d, reading from file %s", aTarget_Name.Data(), (int)(aTheta), aXS_Input_Name.Data())<<endl;

	if(aTarget_Name=="He3")
		aXS_Input_Name = Form("/work/halla/e08014/disk1/yez/inclusive/output/he3/he3_%d.incl", (int)(aTheta));	
	else if(aTarget_Name=="He4")
		aXS_Input_Name = Form("/work/halla/e08014/disk1/yez/inclusive/output/he4/he4_%d.incl", (int)(aTheta));	
	else if(aTarget_Name=="C12")
		aXS_Input_Name = Form("/work/halla/e08014/disk1/yez/inclusive/output/c12/c12_%d.incl", (int)(aTheta));	
	else if(aTarget_Name=="Ca40")
		aXS_Input_Name = Form("/work/halla/e08014/disk1/yez/inclusive/output/ca40/ca40_%d.incl", (int)(aTheta));	
	else if(aTarget_Name=="Ca48")
		aXS_Input_Name = Form("/work/halla/e08014/disk1/yez/inclusive/output/ca48/ca48_%d.incl", (int)(aTheta));	
	else if(aTarget_Name=="H2"||aTarget_Name=="D2") //D2 XS is calculated together with any other target
		aXS_Input_Name = Form("/work/halla/e08014/disk1/yez/inclusive/output/he3/he3_%d.incl", (int)(aTheta));	
	else
		cerr<<" --- in XS_Misak(), I don't know the target name!!!"<<endl;	

	cerr<<Form("--- Working on %s at Theta = %d, reading from file %s", aTarget_Name.Data(), (int)(aTheta), aXS_Input_Name.Data())<<endl;
	
	ifstream infile;
    double aEp_M, aXS_M, aDum, aH2_XS;
	for(int i=0;i<aN;i++){
		aXS_Model[i] = -1.000;
		infile.open(aXS_Input_Name.Data());
        while(!(infile.eof())){
			infile >> aDum >> aDum >> aEp_M >> aDum >> aXS_M >>  aH2_XS >> aDum >> aDum;

			if(fabs(aEp_M-aEp[i])<0.001){
			   if(aTarget_Name=="H2"||aTarget_Name=="D2")
				   aXS_Model[i] = aH2_XS;
			   else
				   aXS_Model[i] = aXS_M;
			   break;
		   }
		}
		infile.close();
		if(aXS_Model[i]<1e-16)
			cerr<<"--- What!!! I can not find the cross section at Ep = "<< aEp[i]<<endl;
	}
}
/*}}}*/

void Merge(const TString& aInputFileName, double* aXbj, double* aEp, double* aXS_EX, double* aXS_EX_Err, double* aXS_MC, int* aN){

	/*Pre-define{{{*/	
	int aline = 0;
	TString aInputFilePath = "/work/halla/e08014/disk1/yez/Xbj_Cross_Section/results/";

	//First Target:
	const TString aName_T1L1 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	const TString aName_T1L2 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	const TString aName_T1R1 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	const TString aName_T1R2 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	TString Target = "C12";
	Target = gGet_Target_Name(aName_T1L1.Data());

	const double E0 = 3.356;
	double Theta = 21.00;
	TString Kin = "3.1";
	Double_t nu,Q2,xbj,angle,Ep,Livetime,Ntg,Ni,Nf_EX,Nf_SAMC,Nf_SAMC_Gen,BinCor;
	Double_t XS,XS_MC,XS_Stat_Err,XS_Sys_Err,Yield_MC,Yield_EX,Y_MC_Sys_Err,Y_EX_Sys_Err,Y_MC_Stat_Err,Y_EX_Stat_Err;
	Double_t Nf_AL, Yield_AL, Y_AL_Stat_Err, Y_AL_Sys_Err;
	TString comment;
	double Ep_Min, Ep_Max;
    double scale = 1.0;
	const double He4_Norm =1.0/0.864;
	/*}}}*/

	/*Define{{{*/
	Double_t XS_T1L1[200], XS_Err_T1L1[200], Ep_T1L1[200], Xbj_T1L1[200],  Xbj_Err_T1L1[200],XS_MC_T1L1[200];
	Double_t XS_T1L2[200], XS_Err_T1L2[200], Ep_T1L2[200], Xbj_T1L2[200],  Xbj_Err_T1L2[200],XS_MC_T1L2[200];
	Double_t XS_T1R1[200], XS_Err_T1R1[200], Ep_T1R1[200], Xbj_T1R1[200],  Xbj_Err_T1R1[200],XS_MC_T1R1[200];
	Double_t XS_T1R2[200], XS_Err_T1R2[200], Ep_T1R2[200], Xbj_T1R2[200],  Xbj_Err_T1R2[200],XS_MC_T1R2[200];
	/*}}}*/

	/*First{{{*/
	TString infile_T1L1 = aInputFilePath +Target+"/"+Com+"/"+ aName_T1L1+"_"+Com;
	cerr<<"==== Reading 1st file from "<<infile_T1L1<<endl;
	ifstream run_T1L1; run_T1L1.open(infile_T1L1);

	Int_t Count_T1L1 =0;
	Target = gGet_Target_Name(infile_T1L1.Data());
	Kin = gGet_Kin_Name(infile_T1L1.Data());

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
	gGet_Ep_Cut(KinName_T1L1.Data(), Arm_T1L1.Data(), &Ep_Min, &Ep_Max);
	Get_XS* xs_T1L1 = new Get_XS();
	xs_T1L1->Load_Table(Target.Data(),KinName_T1L1.Data());

	comment.ReadLine(run_T1L1);
	while ( run_T1L1>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
	{
		if ( Ep>=Ep_Min && Ep<=Ep_Max && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
		{
			XS_T1L1[Count_T1L1] = XS; 
			Ep_T1L1[Count_T1L1] = Ep; 
			XS_Err_T1L1[Count_T1L1] = sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2));
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
	TString infile_T1L2 = aInputFilePath +Target+"/"+Com+"/"+ aName_T1L2+"_"+Com;
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
		Get_XS* xs_T1L2 = new Get_XS();
		xs_T1L2->Load_Table(Target.Data(),KinName_T1L2.Data());

		comment.ReadLine(run_T1L2);
	while ( run_T1L2>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
		{
			if ( Ep>=Ep_Min && Ep<=Ep_Max && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
			{
				XS_T1L2[Count_T1L2] = XS;    
				Ep_T1L2[Count_T1L2] = Ep; 
				XS_Err_T1L2[Count_T1L2] = sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2));
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
	TString infile_T1R1 = aInputFilePath +Target+"/"+Com+"/"+ aName_T1R1+"_"+Com;
	cerr<<"==== Reading 3rd file from "<<infile_T1R1<<endl;
	TString Arm_T1R1="X";
	if(aName_T1R1 =="0")
		bMore_T1R1 = kFALSE;
	if(aName_T1R1.Contains("L_"))
		Arm_T1R1="L";
	else if(aName_T1R1.Contains("R_"))
		Arm_T1R1="R";

	Int_t Count_T1R1=0;
	if(bMore_T1R1){
		ifstream run_T1R1; run_T1R1.open(infile_T1R1);
		Target = gGet_Target_Name(infile_T1R1.Data());
		Kin = gGet_Kin_Name(infile_T1R1.Data());
		if(Kin=="MM"){
			cerr<<"3rd: Kin (3.1,3.2,4.1,4.2,5.0,5.05,5.1,5.2,6.5) = "; cin >> Kin;
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
		Get_XS* xs_T1R1 = new Get_XS();
		xs_T1R1->Load_Table(Target.Data(),KinName_T1R1.Data());

		comment.ReadLine(run_T1R1);
	while ( run_T1R1>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
		{
			if ( Ep>=Ep_Min && Ep<=Ep_Max && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
			{
				XS_T1R1[Count_T1R1] = XS; 
				Ep_T1R1[Count_T1R1] = Ep; 
				XS_Err_T1R1[Count_T1R1] = sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2));
				Xbj_T1R1[Count_T1R1] = xbj;    Xbj_Err_T1R1[Count_T1R1] = 0.0;
				XS_MC_T1R1[Count_T1R1] = xs_T1R1->gGet_XS_Born(E0,Ep,Theta);
				Count_T1R1++;
			}
			if ( !run_T1R1.good() )
				break;
		}
		run_T1R1.close();
		delete xs_T1R1;
		cerr << "  Total XS Points for 3rd  = " << Count_T1R1 <<endl;
	}
	/*}}}*/

	/*Forth{{{*/
	TString KinName_T1R2 ="";
	Bool_t bMore_T1R2 = kTRUE;
	TString infile_T1R2 = aInputFilePath +Target+"/"+Com+"/"+ aName_T1R2+"_"+Com;
	cerr<<"==== Reading 4th file from "<<infile_T1R2<<endl;
	TString Arm_T1R2="X";
	if(aName_T1R2 =="0")
		bMore_T1R2 = kFALSE;
	if(aName_T1R2.Contains("L_"))
		Arm_T1R2="L";
	else if(aName_T1R2.Contains("R_"))
		Arm_T1R2="R";

	Int_t Count_T1R2=0;
	if(bMore_T1R2){
		ifstream run_T1R2; run_T1R2.open(infile_T1R2);
		Target = gGet_Target_Name(infile_T1R2.Data());
		Kin = gGet_Kin_Name(infile_T1R2.Data());
		if(Kin=="MM"){
			cerr<<"4th: Kin (3.1,3.2,4.1,4.2,5.0,5.05,5.1,5.2,6.5) = "; cin >> Kin;
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
		gGet_Ep_Cut(KinName_T1R2.Data(), Arm_T1R2.Data(), &Ep_Min, &Ep_Max);
		Get_XS* xs_T1R2 = new Get_XS();
		xs_T1R2->Load_Table(Target.Data(),KinName_T1R2.Data());

		comment.ReadLine(run_T1R2);
	while ( run_T1R2>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
		{
			if ( Ep>=Ep_Min && Ep<=Ep_Max && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
			{
				XS_T1R2[Count_T1R2] = XS;
				Ep_T1R2[Count_T1R2] = Ep; 
				XS_Err_T1R2[Count_T1R2] = sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2));
				Xbj_T1R2[Count_T1R2] = xbj;    Xbj_Err_T1R2[Count_T1R2] = 0.0;
				XS_MC_T1R2[Count_T1R2] = xs_T1R2->gGet_XS_Born(E0,Ep,Theta);
				Count_T1R2++;
			}
			if ( !run_T1R2.good() )
				break;
		}
		run_T1R2.close();
		delete xs_T1R2;
		cerr << "  Total XS Points for 4th  = " << Count_T1R2 <<endl;
	}
	/*}}}*/

	/*Merge{{{*/
	//Add First Two Input
	Double_t Xbj_T1L[200], Ep_T1L[200], XS_T1L[200], XS_Err_T1L[200], XS_MC_T1L[200];
	Double_t XS_Sum_T1L = 0.0, XS_Err_Sum_T1L = 0.0; Int_t fSum_T1L = 0;
	Int_t Count_T1L = 0, J_T1L = 0;
	for(int i=0;i<Count_T1L1;i++){
		XS_Sum_T1L = XS_T1L1[i];
		XS_Err_Sum_T1L = pow(XS_Err_T1L1[i]/XS_T1L1[i],2);
		fSum_T1L = 1;
		if(bMore_T1L2){
			for(int j=0;j<Count_T1L2;j++){
				//Find the verlapping bin for two setting, such as Kin3.1 and Kin3.2
				if(abs(Ep_T1L1[i]-Ep_T1L2[j])< 0.0001){
					XS_Sum_T1L += XS_T1L2[j];
					XS_Err_Sum_T1L +=	pow(XS_Err_T1L2[j]/XS_T1L2[j],2);
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
		Count_T1L++;
	}
	for(int k=J_T1L+1;k<Count_T1L2;k++){
		Ep_T1L[Count_T1L] = Ep_T1L2[k];
		Xbj_T1L[Count_T1L] = Xbj_T1L2[k];
		XS_T1L[Count_T1L] = XS_T1L2[k];
		XS_Err_T1L[Count_T1L] = XS_Err_T1L2[k];
		XS_MC_T1L[Count_T1L] = XS_MC_T1L2[k];
		Count_T1L++;	
	}
	//Add Last Two Input
	Double_t Xbj_T1R[200],Ep_T1R[200], XS_T1R[200], XS_Err_T1R[200], XS_MC_T1R[200];
	Double_t XS_Sum_T1R = 0.0, XS_Err_Sum_T1R = 0.0; Int_t fSum_T1R = 0;
	Int_t Count_T1R = 0, J_T1R = 0;
	if(bMore_T1R1){
		for(int i=0;i<Count_T1R1;i++){
			XS_Sum_T1R = XS_T1R1[i];
			XS_Err_Sum_T1R = pow(XS_Err_T1R1[i]/XS_T1R1[i],2);
			fSum_T1R = 1;
			if(bMore_T1R2){
				for(int j=0;j<Count_T1R2;j++){
					//Find the verlapping bin for two setting, such as Kin3.1 and Kin3.2
					if(abs(Ep_T1R1[i]-Ep_T1R2[j])< 0.0001){
						XS_Sum_T1R += XS_T1R2[j];
						XS_Err_Sum_T1R +=	pow(XS_Err_T1R2[j]/XS_T1R2[j],2);
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
			Count_T1R++;
		}
		for(int k=J_T1R+1;k<Count_T1R2;k++){
			Ep_T1R[Count_T1R] = Ep_T1R2[k];
			Xbj_T1R[Count_T1R] = Xbj_T1R2[k];
			XS_T1R[Count_T1R] = XS_T1R2[k];
			XS_Err_T1R[Count_T1R] = XS_Err_T1R2[k];
			XS_MC_T1L[Count_T1R] = XS_MC_T1R2[k];
			Count_T1R++;	
		}
	}
	
	vector<double> vXbj1, vEp1, vXbj_Err1,vXS1, vXS_Err1, vXS_MC1;
	double XS_MC_T1, XS_T1, XS_Err_T1;
	int N1 = 0, J1R= 0;
	bool bFound_T1L = kFALSE;
	if(bMore_T1R1){
		for(int i=0;i<Count_T1L; i++){
			bFound_T1L = kFALSE;	
			for(int j=0;j<Count_T1R; j++){
				if( abs(Ep_T1L[i] - Ep_T1R[j])<0.0001 ){
					vEp1.push_back(Ep_T1L[i]);	
					vXbj1.push_back(Xbj_T1L[i]);	
					vXbj_Err1.push_back(0.0);	
					bFound_T1L = kTRUE;	
					XS_T1 = (XS_T1L[i]+XS_T1R[j])/2.0;
					XS_Err_T1 = XS_T1 * sqrt( pow(XS_Err_T1L[i]/XS_T1L[i],2) + pow(XS_Err_T1R[j]/XS_T1R[j],2) );
					//XS_MC_T1 = 	(XS_MC_T1L[i]+XS_MC_T1R[j])/2.0;
					XS_MC_T1 = 	XS_MC_T1L[i];
					vXS1.push_back(XS_T1);
					vXS_Err1.push_back(XS_Err_T1);
					vXS_MC1.push_back( XS_MC_T1 );
					N1++; J1R = j;
				}
			}
			if(!bFound_T1L){
				vEp1.push_back(Ep_T1L[i]);	
				vXbj1.push_back(Xbj_T1L[i]);	
				vXbj_Err1.push_back(0.0);	
				XS_T1 = XS_T1L[i];
				XS_Err_T1 = XS_Err_T1L[i];
				XS_MC_T1 = 	XS_MC_T1L[i];
				vXS1.push_back(XS_T1);
				vXS_Err1.push_back(XS_Err_T1);
				vXS_MC1.push_back( XS_MC_T1 );
				N1++;
			}
		}

		for(int k =J1R+1; k<Count_T1R; k++){
			vEp1.push_back(Ep_T1R[k]);	
			vXbj1.push_back(Xbj_T1R[k]);	
			vXbj_Err1.push_back(0.0);	
			XS_T1 = XS_T1R[k];
			XS_Err_T1 = XS_Err_T1R[k];
			XS_MC_T1 = 	XS_MC_T1R[k];
			vXS1.push_back(XS_T1);
			vXS_Err1.push_back(XS_Err_T1);
			vXS_MC1.push_back( XS_MC_T1 );
			N1++;
		}
	}
	else{
		for(int i=0;i<Count_T1L; i++){
			vEp1.push_back(Ep_T1L[i]);	
			vXbj1.push_back(Xbj_T1L[i]);	
			vXbj_Err1.push_back(0.0);	
			XS_T1 = XS_T1L[i];
			XS_Err_T1 = XS_Err_T1L[i];
			XS_MC_T1 = 	XS_MC_T1L[i];
			vXS1.push_back(XS_T1);
			vXS_Err1.push_back(XS_Err_T1);
			vXS_MC1.push_back( XS_MC_T1 );
			N1++;
		}
	}

	/*}}}*/

	if(Target=="He4")
		scale = He4_Norm;
	else
		scale = 1.0;

	for(unsigned int i=0;i<vEp1.size();i++)
		aEp[i] = vEp1[i];
	for(unsigned int i=0;i<vXbj1.size();i++)
		aXbj[i] = vXbj1[i];
	for(unsigned int i=0;i<vXS1.size();i++){
		aXS_EX[i] = vXS1[i]*scale;
		aXS_EX_Err[i] = vXS_Err1[i]*scale;
	}

	for(unsigned int i=0;i<vXS_MC1.size();i++)
		aXS_MC[i] = vXS_MC1[i];
//	cerr<<" --- Reading XS from Model ..."<<endl;
//	XS_Misak(Target.Data(), Theta, N1, aEp, aXS_MC);
	for(int i=0;i<N1;i++)
		cerr<<Form("---- Ep = %e, XS_MC = %e", aEp[i], aXS_MC[i])<<endl;
	aN[0] = N1; 

}
