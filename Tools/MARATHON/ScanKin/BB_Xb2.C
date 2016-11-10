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
    int Theta[9]={39,41,43,45,47,49,51,53,55};
    double H2_Rate[9]  = {178.50, 132.81,  98.90, 74.33, 56.01, 42.46, 32.19, 24.54, 18.65};
    double H3_Rate[9]  = {108.93,  80.81,  59.98, 44.93, 33.74, 25.49, 19.26, 14.63, 11.07};
    double He3_Rate[9] = {178.34, 133.06,  99.36, 74.89, 56.59, 43.03, 32.72, 25.02, 19.07};

    const double RC_Eff = 0.8*0.8 * 24 * 60 *60;//80% RC, 80% efficienciy, from Hz to count per day
    const int xb_Bin = 18;
    const double xb_Min = 0.25;
    const double xb_Max = 0.97;
    const double xb_step = (xb_Max - xb_Min) / xb_Bin;

    double N_Total = 0.0;
    double Rate_Total = 0.0;

    bool b1 = true;
    //bool b1 = false;
    bool b2 = true;
    //bool b2 = false;
    bool b3 = true;
    //bool b3 = false;


    for(int i=0;i<9;i++){
        int T0 = Theta[i];

        TChain *T = new TChain("T"); 
        T->Add(Form("./results3/BB_rate_marathon_E11_P1_T%d_lowE.root",T0));
        T->Add(Form("./results3/BB_rate_marathon_E11_P2_T%d_lowE.root",T0));

        gStyle->SetOptStat(0);
        T->SetLineStyle(1);
        T->SetLineColor(1);
        T->SetMarkerSize(1.2);
        T->SetMarkerColor(3);
        gStyle->SetPaintTextFormat("4.2e"); 

        if(b1){
            TCanvas *c1=new TCanvas("c1","c1",1400,800);/*{{{*/
            gPad->SetGridx(1);
            gPad->SetGridy(1);
            TH1F* h1 = new TH1F("h1","",xb_Bin,xb_Min,xb_Max);
            h1->SetXTitle("xb");
            h1->SetYTitle("Rate (N/Day)");
            //h1->SetYTitle("Rate (Hz)");
            h1->GetXaxis()->CenterTitle(1);
            h1->GetYaxis()->CenterTitle(1);
            h1->GetXaxis()->SetTitleSize(0.07);
            h1->GetYaxis()->SetTitleSize(0.06);
            h1->GetXaxis()->SetTitleOffset(0.6);
            h1->GetYaxis()->SetTitleOffset(0.6);
            h1->GetXaxis()->SetTitleFont(20);
            h1->GetYaxis()->SetTitleFont(20);

            T->Draw("xb>>h1","(xs_d)*(xb<1.&&W>sqrt(3.))","");
            N_Total = h1->Integral();
            cerr<<"H2: N =  "<<N_Total<<endl;
            Rate_Total = H2_Rate[i]/N_Total*RC_Eff;
            h1->Scale(Rate_Total);

            c1->cd();
            //h1->SetMinimum(-0.5);
            h1->SetMarkerSize(1.5);
            h1->SetMarkerColor(4);
            gStyle->SetTextFont(32);
            h1->Draw("TEXT45");
            TH1F* h12 = (TH1F*) h1->Clone(); h12->SetLineColor(6);h12->Draw("same");

            TLatex *t1 = new TLatex();
            //t1->SetNDC();
            t1->SetTextFont(22);
            t1->DrawLatex(0.6,h1->GetMaximum()*0.9,Form("^{2}H Rates at %d^{#circ} with BigBite",T0));
            t1->SetTextColor(1);
            t1->SetTextFont(22);
            t1->DrawLatex(0.6,h1->GetMaximum()*0.7,Form("(BB Momentum > 0.6 GeV/c)"));

            t1->SetTextColor(4);
            t1->SetTextFont(22);
            t1->DrawLatex(0.3,h1->GetMaximum()*0.05,Form("Rate in each bin in Count/Day"));
            t1->SetTextColor(2);
            t1->SetTextFont(22);
            t1->DrawLatex(0.3,h1->GetMaximum()*0.15,Form("(Average W^{2} in GeV^{2})"));


            TLatex *t11 = new TLatex();
            t11->SetTextFont(32);
            t11->SetTextSize(0.025);
            t11->SetTextAngle(45);
            t11->SetTextColor(2);
            for(int i=0;i<xb_Bin;i++){
                double xb_c = h1->GetBinCenter(i+1);
                double y_c = h1->GetBinContent(i+1);
                T->Draw("W*W",Form("(xs_d)*(abs(xb-%f)<=%f&&W>sqrt(3.))",xb_c,xb_step),"goff");
                TH1F *htemp = (TH1F*) gROOT->FindObject("htemp");
                double W2_mean = htemp->GetMean();
                delete htemp;
                
                t11->DrawLatex(xb_c-0.015,y_c+0.04,Form("(%4.2f)",W2_mean));
                
                cout<<Form("H2:  xb = %4.2f, W = %f", 
                        xb_c, W2_mean)<<endl;
            }

            c1->Print(Form("H2_BB_xb_LowE_%d.png",T0));
            c1->Print(Form("H2_BB_xb_LowE_%d.pdf",T0));
            /*}}}*/
        }

        if(b2){
            TCanvas *c2=new TCanvas("c2","c2",1400,800);/*{{{*/
            gPad->SetGridx(1);
            gPad->SetGridy(1);
            TH1F* h2 = new TH1F("h2","",xb_Bin,xb_Min,xb_Max);
            h2->SetXTitle("xb");
            //h2->SetYTitle("Rate (Hz)");
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
            Rate_Total = H3_Rate[i]/N_Total*RC_Eff;
            h2->Scale(Rate_Total);

            c2->cd();
            //h2->SetMinimum(-0.5);
            h2->SetMarkerSize(1.5);
            h2->SetMarkerColor(4);
            gStyle->SetTextFont(32);
            h2->Draw("TEXT45");
            TH1F* h22 = (TH1F*) h2->Clone(); h22->SetLineColor(6);h22->Draw("same");

            TLatex *t2 = new TLatex();
            t2->SetTextFont(22);
            t2->DrawLatex(0.6,h2->GetMaximum()*0.9,Form("^{3}H Rates at %d^{#circ} with BigBite",T0));
            t2->SetTextColor(1);
            t2->SetTextFont(22);
            t2->DrawLatex(0.6,h2->GetMaximum()*0.7,Form("(BB Momentum > 0.6 GeV/c)"));

            t2->SetTextColor(4);
            t2->SetTextFont(22);
            t2->DrawLatex(0.3,h2->GetMaximum()*0.05,Form("Rate in each bin in Count/Day"));
            t2->SetTextColor(2);
            t2->SetTextFont(22);
            t2->DrawLatex(0.3,h2->GetMaximum()*0.15,Form("(Average W^{2} in GeV^{2})"));



            TLatex *t21 = new TLatex();
            t21->SetTextFont(32);
            t21->SetTextSize(0.025);
            t21->SetTextAngle(45);
            t21->SetTextColor(2);
            for(int i=0;i<xb_Bin;i++){
                double xb_c = h2->GetBinCenter(i+1);
                double y_c = h2->GetBinContent(i+1);
                T->Draw("W*W",Form("(xs_d+xs_n)*(abs(xb-%f)<=%f&&W>sqrt(3.))",xb_c,xb_step),"goff");
                TH1F *htemp = (TH1F*) gROOT->FindObject("htemp");
                double W2_mean = htemp->GetMean();
                delete htemp;
                
                t21->DrawLatex(xb_c-0.015,y_c+0.04,Form("(%4.2f)",W2_mean));
                
                cout<<Form("H2:  xb = %4.2f, W = %f", 
                        xb_c, W2_mean)<<endl;
            }

            c2->Print(Form("H3_BB_xb_LowE_%d.png",T0));
            c2->Print(Form("H3_BB_xb_LowE_%d.pdf",T0));
            /*}}}*/
        }

        if(b3){
            TCanvas *c3=new TCanvas("c3","c3",1400,800);/*{{{*/
            gPad->SetGridx(1);
            gPad->SetGridy(1);
            TH1F* h3 = new TH1F("h3","",xb_Bin,xb_Min,xb_Max);
            h3->SetXTitle("xb");
            //h3->SetYTitle("Rate (Hz)");
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
            Rate_Total = He3_Rate[i]/N_Total*RC_Eff;
            h3->Scale(Rate_Total);

            c3->cd();
            //h3->SetMinimum(-0.5);
            h3->SetMarkerSize(1.5);
            h3->SetMarkerColor(4);
            gStyle->SetTextFont(32);
            h3->Draw("TEXT45");
            TH1F* h32 = (TH1F*) h3->Clone(); h32->SetLineColor(6);h32->Draw("same");

            TLatex *t3 = new TLatex();
            t3->SetTextFont(22);
            t3->DrawLatex(0.6,h3->GetMaximum()*0.9,Form("^{3}He Rates at %d^{#circ} with BigBite",T0));
            t3->SetTextColor(1);
            t3->SetTextFont(22);
            t3->DrawLatex(0.6,h3->GetMaximum()*0.7,Form("(BB Momentum > 0.6 GeV/c)"));

            t3->SetTextColor(4);
            t3->SetTextFont(22);
            t3->DrawLatex(0.3,h3->GetMaximum()*0.05,Form("Rate in each bin in Count/Day"));
            t3->SetTextColor(2);
            t3->SetTextFont(22);
            t3->DrawLatex(0.3,h3->GetMaximum()*0.15,Form("(Average W^{2} in GeV^{2})"));



            TLatex *t31 = new TLatex();
            t31->SetTextFont(32);
            t31->SetTextSize(0.025);
            t31->SetTextAngle(45);
            t31->SetTextColor(2);
            for(int i=0;i<xb_Bin;i++){
                double xb_c = h3->GetBinCenter(i+1);
                double y_c = h3->GetBinContent(i+1);
                T->Draw("W*W",Form("(xs_d+xs_p)*(abs(xb-%f)<=%f&&W>sqrt(3.))",xb_c,xb_step),"goff");
                TH1F *htemp = (TH1F*) gROOT->FindObject("htemp");
                double W2_mean = htemp->GetMean();
                delete htemp;
                
                t31->DrawLatex(xb_c-0.015,y_c+0.04,Form("(%4.2f)",W2_mean));
                
                cout<<Form("H2:  xb = %4.2f, W = %f", 
                        xb_c, W2_mean)<<endl;
            }

            c3->Print(Form("He3_BB_xb_LowE_%d.png",T0));
            c3->Print(Form("He3_BB_xb_LowE_%d.pdf",T0));
            /*}}}*/
        }
        delete T;
    }
}
