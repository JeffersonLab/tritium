#include "XGT2_Main.h"
#include "/work/halla/e08014/disk1/yez/Xbj_Cross_Section/SRC/XGT2_XEMC.h"
#include "TLegend.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
 
/*ReadIn{{{*/
void ReadIn( const TString& aInputFile, double* aXbj, double* aR, double* aR_Err,double* aXS1, double* aXS_Err1,double* aXS2, double* aXS_Err2, int *aN){

	if(aInputFile!="0"){
		cerr<<"--- Reading data from "<< aInputFile.Data()<<endl;
		ifstream infile(aInputFile.Data());

		double cXbj, cR, cR_Err, cXS1, cXS2, cXS_Err1, cXS_Err2;
		int cN = 0;
		//Xbj          Ratio            XS1            XS2      Ratio_Err        XS_Err1        XS_Err2

		TString aCom;
		infile >> aCom >> aCom >> aCom >> aCom >> aCom >> aCom >> aCom;

		while(!(infile.eof())){
			infile >> cXbj >> cR >> cXS1 >> cXS2 >> cR_Err >> cXS_Err1 >> cXS_Err2;	
			aXbj[cN] = cXbj;
			aR[cN] = cR;
			aR_Err[cN] = cR_Err;
			aXS1[cN] = cXS1;
			aXS_Err1[cN] = cXS_Err1;
			aXS2[cN] = cXS2;
			aXS_Err2[cN] = cXS_Err2;
			cerr<<Form("====  Xbj = %5.3f, R = %5.3f ", aXbj[cN], aR[cN])<<endl;
			cN++;	
		}
		aN[0] = cN;
		infile.close();
	}
}
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


//	TString Input_Name1 = "./left/XS_5.0_Ca48_Ca40_DpTh40_DpTh40.dat",  Title_Name1 = "R(Ca48/Ca40) at 25^{#circ} on HRS-L";
//	TString Input_Name2 = "./right/XS_5.0_Ca48_Ca40_DpTh40_DpTh40.dat", Title_Name2 = "R(Ca48/Ca40) at 25^{#circ} on HRS-R";
//	TString Input_Name3 = "./left/XS_5.1_Ca48_Ca40_DpTh40_DpTh40.dat",  Title_Name3 = "R(Ca48/Ca40) at 25^{#circ} on HRS-L";
//	TString Input_Name4 = "./right/XS_5.1_Ca48_Ca40_DpTh40_DpTh40.dat", Title_Name4 = "R(Ca48/Ca40) at 25^{#circ} on HRS-R";
//	TString Output_Name = "Ca48_Ca40_25";

//	TString Input_Name1 = "./left/XS_5.0_C12_He3_DpTh40_DpTh40.dat",  Title_Name1 = "R(C12/He3) at 25^{#circ} on HRS-L";
//	TString Input_Name2 = "./right/XS_5.0_C12_He3_DpTh40_DpTh40.dat", Title_Name2 = "R(C12/He3) at 25^{#circ} on HRS-R";
//	TString Input_Name3 = "./left/XS_5.1_C12_He3_DpTh40_DpTh40.dat",  Title_Name3 = "R(C12/He3) at 25^{#circ} on HRS-L";
//	TString Input_Name4 = "./right/XS_5.1_C12_He3_DpTh40_DpTh40.dat", Title_Name4 = "R(C12/He3) at 25^{#circ} on HRS-R";
//	TString Output_Name = "C12_He3_25";

	TString Input_Name1 = "./left/XS_5.0_He4_He3_DpTh40_DpTh40.dat",  Title_Name1 = "R(He4/He3) at 25^{#circ} on HRS-L";
	TString Input_Name2 = "./right/XS_5.0_He4_He3_DpTh40_DpTh40.dat", Title_Name2 = "R(He4/He3) at 25^{#circ} on HRS-R";
	TString Input_Name3 = "./left/XS_5.1_He4_He3_DpTh40_DpTh40.dat",  Title_Name3 = "R(He4/He3) at 25^{#circ} on HRS-L";
	TString Input_Name4 = "./right/XS_5.1_He4_He3_DpTh40_DpTh40.dat", Title_Name4 = "R(He4/He3) at 25^{#circ} on HRS-R";
	TString Output_Name = "He4_He3_25";


//    TString Input_Name1 = "./left/XS_4.2_He4_He3_DpTh40_DpTh40.dat",  Title_Name1 = "R(He4/He3) at 23^{#circ} on HRS-L";
//	TString Input_Name2 = "./right/XS_4.2_He4_He3_DpTh40_DpTh40.dat", Title_Name2 = "R(He4/He3) at 23^{#circ} on HRS-R";
//	TString Input_Name3 = "0", Title_Name3 = "0";
//	TString Input_Name4 = "0", Title_Name4 = "0";
//	TString Output_Name = "He4_He3_23";


	double R_EX_1[50], R_EX_2[50], R_EX_3[50], R_EX_4[50];
	double R_EX_Err_1[50], R_EX_Err_2[50], R_EX_Err_3[50], R_EX_Err_4[50];
	double xbj_1[50], xbj_2[50], xbj_3[50], xbj_4[50];
	double XS_EX_11[50], XS_EX_21[50], XS_EX_31[50], XS_EX_41[50];
	double XS_EX_12[50], XS_EX_22[50], XS_EX_32[50], XS_EX_42[50];
	double XS_EX_Err11[50], XS_EX_Err21[50], XS_EX_Err31[50], XS_EX_Err41[50];
	double XS_EX_Err12[50], XS_EX_Err22[50], XS_EX_Err32[50], XS_EX_Err42[50];
	int N_1, N_2, N_3, N_4;

	double Zero[50];
	for(int i=0;i<50;i++)
		Zero[i]=0.0;

    ReadIn( Input_Name1.Data(), xbj_1, R_EX_1, R_EX_Err_1, XS_EX_11, XS_EX_Err11, XS_EX_12, XS_EX_Err12, &N_1);
    ReadIn( Input_Name2.Data(), xbj_2, R_EX_2, R_EX_Err_2, XS_EX_21, XS_EX_Err21, XS_EX_22, XS_EX_Err22, &N_2);
    ReadIn( Input_Name3.Data(), xbj_3, R_EX_3, R_EX_Err_3, XS_EX_31, XS_EX_Err31, XS_EX_32, XS_EX_Err32, &N_3);
    ReadIn( Input_Name4.Data(), xbj_4, R_EX_4, R_EX_Err_4, XS_EX_41, XS_EX_Err41, XS_EX_42, XS_EX_Err42, &N_4);


	/*Plot{{{*/
	cerr<<endl<<"@@@@@@@@ Now starting to make some plots ..."<<endl;

	gStyle->SetOptStat(0);
	
	/*Range{{{*/
	double RMAX = -1000, RMIN = 1000;
	double XMAX = -1000, XMIN = 1000;
	for(int i=0; i<N_1; i++){
		if(RMAX<R_EX_1[i]) RMAX = R_EX_1[i];
		if(XMAX <xbj_1[i]) XMAX =  xbj_1[i];
	}
	for(int i=0; i<N_2; i++){
		if(RMAX<R_EX_2[i]) RMAX = R_EX_2[i];
		if(XMAX <xbj_2[i]) XMAX =  xbj_2[i];
	}
	if(Input_Name3!="0"){
		for(int i=0; i<N_3; i++){
			if(RMAX<R_EX_3[i]) RMAX = R_EX_3[i];
			if(XMAX <xbj_3[i]) XMAX =  xbj_3[i];
		}
	}
	if(Input_Name4!="0"){
		for(int i=0; i<N_4; i++){
			if(RMAX<R_EX_4[i]) RMAX = R_EX_4[i];
			if(XMAX <xbj_4[i]) XMAX =  xbj_4[i];
		}
	}
	for(int i=0; i<N_1; i++){
		if(RMIN>R_EX_1[i]) RMIN = R_EX_1[i];
		if(XMIN >xbj_1[i]) XMIN =  xbj_1[i];
	}
	for(int i=0; i<N_2; i++){
		if(RMIN>R_EX_2[i]) RMIN = R_EX_2[i];
		if(XMIN >xbj_2[i]) XMIN =  xbj_2[i];
	}
	if(Input_Name3!="0"){
		for(int i=0; i<N_3; i++){
			if(RMIN>R_EX_3[i]) RMIN = R_EX_3[i];
			if(XMIN >xbj_3[i]) XMIN =  xbj_3[i];
		}
	}
	if(Input_Name4!="0"){
		for(int i=0; i<N_4; i++){
			if(RMIN>R_EX_4[i]) RMIN = R_EX_4[i];
			if(XMIN >xbj_4[i]) XMIN =  xbj_4[i];
		}
	}
	if(RMAX>20)
		RMAX=20;

	cerr<<Form("--- X_Min = %f, X_Max = %f, R_Min = %f, R_Max = %f", XMIN, XMAX, RMIN, RMAX)<<endl;

	//	RMIN = 0.0;
	//	RMAX = 2.0;
	/*}}}*/

	TCanvas *c3 = new TCanvas("c3","c3",1200,700);
	c3->cd();
	TH2F *h3 = new TH2F("h3","",300,XMIN-0.2,XMAX+0.2,300,RMIN-0.0,RMAX+0.0);
	h3->SetStats(kFALSE);
	h3->SetXTitle("x");
	h3->GetXaxis()->CenterTitle(1);
	h3->GetXaxis()->SetTitleFont(22);
	h3->GetXaxis()->SetTitleSize(0.08);
	h3->GetXaxis()->SetTitleOffset(0.7);
	h3->GetYaxis()->SetTitle("Ratio");
	h3->GetYaxis()->CenterTitle(1);
	h3->GetYaxis()->SetTitleFont(22);
	h3->GetYaxis()->SetTitleSize(0.07);
	h3->GetYaxis()->SetTitleOffset(0.8);
	h3->GetYaxis()->CenterTitle();
	gPad->SetLogy(0);
	h3->Draw();

	TLegend *l1 = new TLegend(0.17,0.74,0.50, 0.94,"");
	l1->SetTextSize(0.043);
	l1->SetTextFont(32);

	cerr<<"--- Ploting 1st set of data"<<endl;

	TGraphErrors *ex1 = new TGraphErrors(N_1,xbj_1, R_EX_1, Zero,R_EX_Err_1);
	ex1->SetMarkerStyle(26);
	ex1->SetMarkerSize(1.5);
	ex1->SetMarkerColor(2);
	ex1->SetLineColor(2);
	c3->cd(); ex1->Draw("P");
	l1->AddEntry(ex1, Form("%s",Title_Name1.Data()),"p");

	cerr<<"--- Ploting 2nd set of data"<<endl;
	TGraphErrors *ex2 = new TGraphErrors(N_2,xbj_2, R_EX_2, Zero,R_EX_Err_2);
	ex2->SetMarkerStyle(25);
	ex2->SetMarkerSize(1.5);
	ex2->SetMarkerColor(4);
	ex2->SetLineColor(4);
	c3->cd(); ex2->Draw("P");
	l1->AddEntry(ex2, Form("%s",Title_Name2.Data()),"p");

	if(Input_Name3!="0"){
	cerr<<"--- Ploting 3rd set of data"<<endl;
		TGraphErrors *ex3 = new TGraphErrors(N_3,xbj_3, R_EX_3, Zero,R_EX_Err_3);
		ex3->SetMarkerStyle(24);
		ex3->SetMarkerSize(1.5);
		ex3->SetMarkerColor(1);
		ex3->SetLineColor(1);
		c3->cd(); ex3->Draw("P");
		l1->AddEntry(ex3, Form("%s",Title_Name3.Data()),"p");
	}

	if(Input_Name4!="0"){
	cerr<<"--- Ploting 4th set of data"<<endl;
		TGraphErrors *ex4 = new TGraphErrors(N_4,xbj_4, R_EX_4, Zero,R_EX_Err_4);
		ex4->SetMarkerStyle(27);
		ex4->SetMarkerSize(1.5);
		ex4->SetMarkerColor(6);
		ex4->SetLineColor(1);
		c3->cd(); ex4->Draw("P");
		l1->AddEntry(ex4, Form("%s",Title_Name4.Data()),"p");
	}

	c3->cd(); l1->Draw();

	c3->Print(Form("./plot/%s.png",Output_Name.Data()));

	/*}}}*/
	return 0;
}
/*}}}*/

