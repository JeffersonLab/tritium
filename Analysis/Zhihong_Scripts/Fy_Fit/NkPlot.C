/*Headers{{{*/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include <string>
#include <cstdio>

#include "math.h"
#include "stdio.h"
#include <stdlib.h>

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TF1.h>
#include <TH1.h>
#include <TAxis.h>
#include <TLine.h>
#include <TLegend.h>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TObject.h"
#include "TSystem.h"
#include <TVirtualFitter.h>
#include <TMinuit.h>
//#include <TFitterMinuit.h>
#include <TApplication.h>
#include <vector>
using namespace std;

/*}}}*/
#include "XEMC_Born.h"

double Y_Min = 1000;
double Y_Max =-1000;

/*func{{{*/
Double_t func(Double_t *x, Double_t *par)
{
	Double_t scale = 1000.0;
	// Scale the paramters, base on XEM Model
	Double_t y = scale*x[0];
	Double_t A = par[0];
	Double_t f0 = par[1]/scale;
	Double_t bigB = par[2]/scale;
	Double_t alpha = par[3];
	Double_t a = par[4]/scale;
	Double_t b = par[5]/scale;

	Double_t Fy = -111.11;

	if(A<=2.0)//F(y) Fucntion: for light nuclear
		Fy = ((f0-bigB)*(alpha*alpha*exp(-a*a*y*y))/(alpha*alpha+y*y)+bigB*exp(-b*abs(y)))*scale;
	else  // F(y) Fucntion: for heavy nuclear, A>2
		Fy = (f0-bigB)*(alpha*alpha*exp(-a*a*y*y))/(alpha*alpha+y*y)+bigB*exp(-b*b*y*y);

	return Fy;
}
/*}}}*/

/*Momentum Distribution n(k){{{*/
Double_t nk(Double_t *x, Double_t *par)
{
  Double_t scale = 1000.0;
 
  // Scale the paramters, base on XEM Model
  Double_t y = scale*x[0];
  Double_t A = par[0];
  Double_t f0 = par[1]/scale;
  Double_t bigB = par[2]/scale;
  Double_t alpha = par[3];
  Double_t a = par[4]/scale;
  Double_t b = par[5]/scale;
  Double_t pi = 3.1415926;
  Double_t fnk = -111.11;

  if(A<3.0){  
	  if(y<0)
		  fnk =  1.0/pi*(f0-bigB)
			  * alpha*alpha/(alpha*alpha+y*y)
			  * (a*a+1/(alpha*alpha+y*y))
			  * exp(-a*a*y*y)
			  - bigB*b*exp(-b*abs(y))/2.0/pi/y;
	  else if(y>0)
		  fnk =  1.0/pi*(f0-bigB)
			  * alpha*alpha/(alpha*alpha+y*y)
			  * (a*a+1/(alpha*alpha+y*y))
			  * exp(-a*a*y*y)
			  + bigB*b*exp(-b*abs(y))/2.0/pi/y;
  }
  else
	  fnk =  1.0/pi*(f0-bigB)
		  * alpha*alpha/(alpha*alpha+y*y)
		  * (a*a+1/(alpha*alpha+y*y))
		  * exp(-a*a*y*y)
		  +b*b*bigB*exp(-b*b*y*y)/pi;

  fnk*=scale;

  fnk*=1e8;

	if(fnk>Y_Max)
		Y_Max = fnk;
	if(fnk<Y_Min)
		Y_Min = fnk;

  return fnk;
}
/*}}}*/

int main(){
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
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);
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

	/*Nk Data{{{*/
	XEM_TGT* H2 = new XEM_TGT();  
	H2->GetValueAZ(2,1,1);
	H2->Print();
//	TF1 *FyH2=new TF1("FyH2",func,-1.0,0.2,6);
//	FyH2->SetParameters(2.0,H2->f0,H2->B, H2->alpha, H2->a, H2->b);
//	FyH2->SetLineColor(1);
	TF1 *NkH2=new TF1("NkH2",nk,0.0,0.92,6);
	NkH2->SetParameters(2.0,H2->f0,H2->B, H2->alpha, H2->a, H2->b);
	NkH2->SetLineColor(1);
	NkH2->SetLineStyle(1);

	XEM_TGT* He3 = new XEM_TGT();  
	He3->GetValueAZ(3,2,1);
	He3->Print();
	TF1 *NkHe3=new TF1("NkHe3",nk,0.0,0.92,6);
	NkHe3->SetParameters(3.0,He3->f0,He3->B, He3->alpha, He3->a, He3->b);
	NkHe3->SetLineColor(2);
	NkHe3->SetLineStyle(2);

	XEM_TGT* He4 = new XEM_TGT();  
	He4->GetValueAZ(4,2,1);
	He4->Print();
	TF1 *NkHe4=new TF1("NkHe4",nk,0.0,0.92,6);
	NkHe4->SetParameters(4.0,He4->f0,He4->B, He4->alpha, He4->a, He4->b);
	NkHe4->SetLineColor(3);
	NkHe4->SetLineStyle(3);

	XEM_TGT* C12 = new XEM_TGT();  
	C12->GetValueAZ(12,6,1);
	C12->Print();
	TF1 *NkC12=new TF1("NkC12",nk,0.0,0.92,6);
	NkC12->SetParameters(12.0,C12->f0,C12->B, C12->alpha, C12->a, C12->b);
	NkC12->SetLineColor(4);
	NkC12->SetLineStyle(5);

	XEM_TGT* Ca40 = new XEM_TGT();  
	Ca40->GetValueAZ(40,20,1);
	Ca40->Print();
	TF1 *NkCa40=new TF1("NkCa40",nk,0.0,0.92,6);
	NkCa40->SetParameters(40.0,Ca40->f0,Ca40->B, Ca40->alpha, Ca40->a, Ca40->b);
	NkCa40->SetLineColor(6);
	NkCa40->SetLineStyle(6);

	XEM_TGT* Ca48 = new XEM_TGT();  
	Ca48->GetValueAZ(48,20,1);
	Ca48->Print();
	TF1 *NkCa48=new TF1("NkCa48",nk,0.0,0.92,6);
	NkCa48->SetParameters(48.0,Ca48->f0,Ca48->B, Ca48->alpha, Ca48->a, Ca48->b);
	NkCa48->SetLineColor(7);
	NkCa48->SetLineStyle(9);

	TCanvas *c1 = new TCanvas("c1","c1",1200,800);
    gStyle->SetOptStat(0);
	gPad->SetLogy(1);
	//TH2F *h1 = new TH2F("h1","Momentum Distribution", 300, -0.1, 0.92, 300,Y_Min*0.95, Y_Max*1.2);
	TH2F *h1 = new TH2F("h1","Momentum Distribution", 300, -0.02, 0.922, 300,1.*1e-4, 1e3);
	h1->SetFillColor(0);
	h1->SetTitle("");
	h1->SetXTitle("k");
	h1->GetXaxis()->CenterTitle(1);
	h1->GetXaxis()->SetTitleFont(32);
	h1->GetXaxis()->SetTitleSize(0.08);
	h1->GetXaxis()->SetTitleOffset(0.7);
	h1->GetXaxis()->SetNdivisions(10,2,5);
	h1->SetYTitle("n(k)");
	h1->GetYaxis()->CenterTitle(1);
	h1->GetYaxis()->SetTitleFont(32);
	h1->GetYaxis()->SetTitleSize(0.07);
	h1->GetYaxis()->SetTitleOffset(0.8);
	h1->Draw();

	TLegend *t1 = new TLegend(0.5,0.5,0.85,0.85,Form("from E08-014 Data: "));
	t1->SetDefaults();
	t1->SetBorderSize(0);
	t1->SetTextSize(0.05);
	t1->SetTextFont(32);

//	NkH2->Draw("");
	NkH2->Draw("same");
	t1->AddEntry(NkH2,Form("^{2}H"),"l");
	NkHe3->Draw("same");
	t1->AddEntry(NkHe3,Form("^{3}He"),"l");
	NkHe4->Draw("same");
	t1->AddEntry(NkHe4,Form("^{4}He"),"l");
	NkC12->Draw("same");
	t1->AddEntry(NkC12,Form("^{12}C"),"l");
	NkCa40->Draw("same");
	t1->AddEntry(NkCa40,Form("^{40}Ca"),"l");
	NkCa48->Draw("same");
	t1->AddEntry(NkCa48,Form("^{48}Ca"),"l");

    t1->Draw();
    c1->Print("Mom_Dis_XGT2.png");	
    c1->Print("Mom_Dis_XGT2.pdf");	
	t1->Delete();
    /*}}}*/

	/*Compare with Theory{{{*/
	ifstream cdf("cdf_mom_dis.dat");
	TString acom;
	cdf >> acom >> acom >> acom >> acom >> acom >> acom;
    double cK[82], cH2[82], cHe3[82], cHe4[82], cC12[82],cCa40[82];
  
	int cN=0;
	while(!(cdf.eof())){
	   cdf >> cK[cN] >> cH2[cN] >> cHe3[cN] >> cHe4[cN] >> cC12[cN] >> cCa40[cN];
	   cN++;
	}
    cdf.close();
	cN--;
	cerr<<" N = "<<cN<<endl;
    
	TCanvas *c2 = new TCanvas("c2","c2",1200,800);
	gPad->SetLogy(1);
    TH2F *h2 = (TH2F*)h1->Clone();
	gPad->SetLogy(1);
    h2->Draw();	
	TLegend *t2 = new TLegend(0.5,0.5,0.85,0.85,Form("from Theory: "));
	t2->SetDefaults();
	t2->SetBorderSize(0);
	t2->SetTextSize(0.05);
	t2->SetTextFont(32);

	TGraph *gh2=new TGraph(cN,cK,cH2);
	gh2->SetLineStyle(1);
	gh2->SetLineWidth(2);
	gh2->SetLineColor(1);
	gh2->Draw("C");
	t2->AddEntry(gh2,Form("^{2}H"),"l");

	TGraph *ghe3=new TGraph(cN,cK,cHe3);
	ghe3->SetLineStyle(1);
	ghe3->SetLineWidth(2);
	ghe3->SetLineColor(2);
	ghe3->Draw("C");
	t2->AddEntry(ghe3,Form("^{3}He"),"l");

	TGraph *ghe4=new TGraph(cN,cK,cHe4);
	ghe4->SetLineStyle(1);
	ghe4->SetLineWidth(2);
	ghe4->SetLineColor(3);
	ghe4->Draw("C");
	t2->AddEntry(ghe4,Form("^{4}He"),"l");

	TGraph *gc12=new TGraph(cN,cK,cC12);
	gc12->SetLineStyle(1);
	gc12->SetLineWidth(2);
	gc12->SetLineColor(4);
	gc12->Draw("C");
	t2->AddEntry(gc12,Form("^{12}C"),"l");

	TGraph *gca40=new TGraph(cN,cK,cCa40);
	gca40->SetLineStyle(1);
	gca40->SetLineWidth(2);
	gca40->SetLineColor(6);
	gca40->Draw("C");
	t2->AddEntry(gc12,Form("^{40}Ca"),"l");

    t2->Draw();
    c2->Print("Mom_Dis_Theo.png");	
    c2->Print("Mom_Dis_Theo.pdf");	

	TCanvas *c3 = new TCanvas("c3","c3",1200,800);
    c3->Divide(3,2);
	gStyle->SetOptStat(0);
	gPad->SetLogy(1);
    TH2F *hH2 = (TH2F*)h1->Clone();
    TH2F *hHe3 = (TH2F*)h1->Clone();
    TH2F *hHe4 = (TH2F*)h1->Clone();
    TH2F *hC12 = (TH2F*)h1->Clone();
    TH2F *hCa40 = (TH2F*)h1->Clone();

	c3->cd(1);gPad->SetLogy(1); hH2->SetTitle("H2"); hH2->Draw("");
	NkH2->SetLineColor(2); NkH2->SetLineStyle(1); NkH2->Draw("same");
	gh2->SetLineColor(4); gh2->SetLineStyle(1); gh2->Draw("C");
	c3->cd(2);gPad->SetLogy(1); hHe3->SetTitle("He3");  hHe3->Draw("");
	NkHe3->SetLineColor(2); NkHe3->SetLineStyle(1); NkHe3->Draw("same");
	ghe3->SetLineColor(4); ghe3->SetLineStyle(1); ghe3->Draw("C");
	c3->cd(3);gPad->SetLogy(1); hHe4->SetTitle("He4"); hHe4->Draw("");
	NkHe4->SetLineColor(2); NkHe4->SetLineStyle(1); NkHe4->Draw("same");
	ghe4->SetLineColor(4); ghe4->SetLineStyle(1); ghe4->Draw("C");
	c3->cd(4);gPad->SetLogy(1); hC12->SetTitle("C12"); hC12->Draw("");
	NkC12->SetLineColor(2); NkC12->SetLineStyle(1); NkC12->Draw("same");
	gc12->SetLineColor(4); gc12->SetLineStyle(1); gc12->Draw("C");
	c3->cd(5);gPad->SetLogy(1); hCa40->SetTitle("Ca40"); hCa40->Draw("");
	NkCa40->SetLineColor(2); NkCa40->SetLineStyle(1); NkCa40->Draw("same");
	gca40->SetLineColor(4); gca40->SetLineStyle(1); gca40->Draw("C");


    c3->Print("Mom_Dis_Comp.png");	
	/*}}}*/

}

