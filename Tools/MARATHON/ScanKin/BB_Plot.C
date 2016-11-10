/*Include{{{*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <errno.h>
#include <sstream>
//#include <iterator>
//from root
#include <TMath.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TMatrix.h>
#include <TFile.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TChain.h>
#include <TTree.h>
#include <TStopwatch.h>
/*}}}*/
using namespace::std;
using namespace::TMath;

#include "../DIS_Lite.h"
int main(){
    int Theta[10]={39,41,43,45,47,49,51,53,55,57};
    double H2_Rate[10]  = {178.53, 132.57,  98.95, 74.23, 55.99, 42.40, 32.14, 24.55, 18.66, 14.21};
    double H3_Rate[10]  = {108.96,  80.65,  60.00, 44.86, 33.73, 25.45, 19.23, 14.64, 11.08,  8.41};
    double He3_Rate[10] = {178.38, 132.82,  99.41, 74.78, 56.58, 42.97, 32.67, 25.03, 19.08, 14.58};

    const double RC_Eff = 0.8*0.8;

    //for(int i=0;i<5;i++){
    for(int i=0;i<10;i++){
        int T0 = Theta[i];

        //TFile *file = new TFile(Form("BB_rate_marathon_E11_P2_T%d.root",T0),"r");
        //TTree* T = (TTree*) file->Get("T");

        TChain *T = new TChain("T");
        T->Add(Form("BB_rate_marathon_E11_P2_T%d_lowE.root",T0));

        gStyle->SetOptStat(0);
        T->SetLineStyle(1);
        T->SetLineColor(1);
        T->SetMarkerSize(1.2);
        T->SetMarkerColor(4);
        gStyle->SetPaintTextFormat("5.2e"); 

        TCanvas *c1=new TCanvas("c1","c1",1400,800);
        gPad->SetGridx(1);
        gPad->SetGridy(1);
        TH1F* h1 = new TH1F("h1","",18,0.25,0.97);
        h1->SetXTitle("x");
        h1->SetYTitle("Rate (N/Day)");
        h1->GetXaxis()->CenterTitle(1);
        h1->GetYaxis()->CenterTitle(1);
        h1->GetXaxis()->SetTitleSize(0.07);
        h1->GetYaxis()->SetTitleSize(0.06);
        h1->GetXaxis()->SetTitleOffset(0.6);
        h1->GetYaxis()->SetTitleOffset(0.6);
        h1->GetXaxis()->SetTitleFont(20);
        h1->GetYaxis()->SetTitleFont(20);

        double N_Total = 0.0;
        double Rate_Total = 0.0;

        T->Draw("xb>>h1","(xs_d)*(xb<1.&&W>sqrt(3.))","");
        N_Total = h1->Integral();
        cerr<<"H2: N =  "<<N_Total<<endl;
        Rate_Total = H2_Rate[i]/N_Total*24*60.*60.*RC_Eff;
        h1->Scale(Rate_Total);

        h1->SetMarkerSize(1.5);
        h1->SetMarkerColor(4);
        h1->Draw("TEXT45");

        T->Draw("xb",Form("(xs_d*%f)*(xb<1.&&W>sqrt(3.))",Rate_Total),"same");

        TLatex *t1 = new TLatex();
        t1->SetNDC();
        t1->SetTextFont(32);
        t1->DrawLatex(0.6,0.8,Form("^{2}H Rates at %d^{#circ} with BigBite",T0));

        TCanvas *c2=new TCanvas("c2","c2",1400,800);
        gPad->SetGridx(1);
        gPad->SetGridy(1);
        TH1F* h2 = new TH1F("h2","",18,0.25,0.97);
        h2->SetXTitle("x");
        h2->SetYTitle("Rate (N/Day)");
        h2->GetXaxis()->CenterTitle(1);
        h2->GetYaxis()->CenterTitle(1);
        h2->GetXaxis()->SetTitleSize(0.07);
        h2->GetYaxis()->SetTitleSize(0.06);
        h2->GetXaxis()->SetTitleOffset(0.6);
        h2->GetYaxis()->SetTitleOffset(0.6);
        h2->GetXaxis()->SetTitleFont(20);
        h2->GetYaxis()->SetTitleFont(20);


        T->Draw("xb>>h2","(xs_d+xs_n)*(xb<1.&&W>sqrt(3.))","");
        N_Total = h2->Integral();
        cerr<<"H3: N =  "<<N_Total<<endl;
        Rate_Total = H3_Rate[i]/N_Total*24*60.*60.*RC_Eff;
        h2->Scale(Rate_Total);

        h2->SetMarkerSize(1.5);
        h2->SetMarkerColor(4);
        h2->Draw("TEXT45");

        T->Draw("xb",Form("((xs_d+xs_n)*%f)*(xb<1.&&W>sqrt(3.))",Rate_Total),"same");

        TLatex *t2 = new TLatex();
        t2->SetNDC();
        t2->SetTextFont(32);
        t2->DrawLatex(0.6,0.8,Form("^{3}H Rates at %d^{#circ} with BigBite",T0));

        TCanvas *c3=new TCanvas("c3","c3",1400,800);
        gPad->SetGridx(1);
        gPad->SetGridy(1);
        TH1F* h3 = new TH1F("h3","",18,0.25,0.97);
        h3->SetXTitle("x");
        h3->SetYTitle("Rate (N/Day)");
        h3->GetXaxis()->CenterTitle(1);
        h3->GetYaxis()->CenterTitle(1);
        h3->GetXaxis()->SetTitleSize(0.07);
        h3->GetYaxis()->SetTitleSize(0.06);
        h3->GetXaxis()->SetTitleOffset(0.6);
        h3->GetYaxis()->SetTitleOffset(0.6);
        h3->GetXaxis()->SetTitleFont(20);
        h3->GetYaxis()->SetTitleFont(20);

        T->Draw("xb>>h3","(xs_d+xs_p)*(xb<1.&&W>sqrt(3.))","");
        N_Total = h3->Integral();
        cerr<<"He3: N =  "<<N_Total<<endl;
        Rate_Total = He3_Rate[i]/N_Total*24*60.*60.*RC_Eff;
        h3->Scale(Rate_Total);

        h3->SetMarkerSize(1.5);
        h3->SetMarkerColor(4);
        h3->Draw("TEXT45");

        T->Draw("xb",Form("((xs_d+xs_p)*%f)*(xb<1.&&W>sqrt(3.))",Rate_Total),"same");

        TLatex *t3 = new TLatex();
        t3->SetNDC();
        t3->SetTextFont(32);
        t3->DrawLatex(0.6,0.8,Form("^{3}He Rates at %d^{#circ} with BigBite",T0));

        c1->Print(Form("H2_BB_xbj_Rate_LowE_%d.png",T0));
        c1->Print(Form("H2_BB_xbj_Rate_LowE_%d.pdf",T0));
        c2->Print(Form("H3_BB_xbj_Rate_LowE_%d.png",T0));
        c2->Print(Form("H3_BB_xbj_Rate_LowE_%d.pdf",T0));
        c3->Print(Form("He3_BB_xbj_Rate_LowE_%d.png",T0));
        c3->Print(Form("He3_BB_xbj_Rate_LowE_%d.pdf",T0));
        
        T->Delete();
    }
}
