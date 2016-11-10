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
const double Deg2Rad = 3.1415926/180.0;
int main(){
    //const double RC_Eff = 0.8*0.8*24*60.*60.;/*{{{*/
    const double RC_Eff = 0.8*0.8;

    const double Ep_Min = 0.95, Ep_Max = 2.15;
    const int Ep_Bin = 12;
    const double Ep[Ep_Bin] = {1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1};


    const double E0 = 11.0;
    ////47
    int T0 = 47;
    double xbj[Ep_Bin] = {0.37, 0.41, 0.46, 0.50, 0.54, 0.59, 0.63, 0.68, 0.73, 0.78, 0.83, 0.88};
    double Sig_Pim[Ep_Bin] = {362.60, 190.50, 98.30, 49.71, 24.57, 11.84, 5.525, 2.481, 1.060, 0.423, 0.1528, 0.04657};
    double Sig_Pi0[Ep_Bin] = {501.00, 259.30,132.10, 66.08, 32.39, 15.50, 7.203, 3.230, 1.383, 0.556, 0.2041, 0.06455};

    ////43
    //const int T0 = 43;
    //double xbj[Ep_Bin] = {0.31,0.35,0.39,0.42,0.46,0.50,0.54,0.58,0.62,0.66,0.70,0.74};
    //const double Sig_Pim[Ep_Bin] = {647.50, 366.90,205.20,113.20, 61.57, 32.97, 17.35, 8.94, 4.51, 2.21, 1.044, 0.4723};
     //const double Sig_Pi0[Ep_Bin] = {891.80, 497.50,274.40,149.60, 80.57, 42.77, 22.34, 11.46,5.748,2.809,1.329, 0.6037};

    ////with E'>0.6GeV
    //int Theta[10]={39,41,43,45,47,49,51,53,55,57};
    //double H2_Rate[10]  = {178.53, 132.57,  98.95, 74.23, 55.99, 42.40, 32.14, 24.55, 18.66, 14.21};
    //double H3_Rate[10]  = {108.96,  80.65,  60.00, 44.86, 33.73, 25.45, 19.23, 14.64, 11.08,  8.41};
    //double He3_Rate[10] = {178.38, 132.82,  99.41, 74.78, 56.58, 42.97, 32.67, 25.03, 19.08, 14.58};
    double H2_Rate = 55.99;//63.56;[>{{{<]
    double H3_Rate = 33.73;//37.59;
    double He3_Rate = 56.58;//65.24;
    double Window_H2_Rate =100.84;//73.54;//43degrees
    double Window_H3_Rate =100.84;//73.54;//43degrees
    double Window_He3_Rate =100.84*25/20;//73.54*25/20;//43degrees
    if(T0==43){
        H2_Rate = 98.95;
        H3_Rate = 60.00;
        He3_Rate = 99.41;
        Window_H2_Rate =178.04;//43degrees
        Window_H3_Rate =178.04;//43degrees
        Window_He3_Rate =178.04*25/20;//43degrees
    }

    double Sig_H2[Ep_Bin] = {Ep_Bin*0};
    double Sig_H3[Ep_Bin] = {Ep_Bin*0};
    double Sig_He3[Ep_Bin] = {Ep_Bin*0};
    double Sig_Al[Ep_Bin] = {Ep_Bin*0};

    double rate_pim_h2[Ep_Bin] = {Ep_Bin*0};
    double rate_pim_h3[Ep_Bin] = {Ep_Bin*0};
    double rate_pim_he3[Ep_Bin] = {Ep_Bin*0};

    double rate_h2[Ep_Bin] = {Ep_Bin*0};
    double rate_h3[Ep_Bin] = {Ep_Bin*0};
    double rate_he3[Ep_Bin] = {Ep_Bin*0};

    double r_pim2e_h2[Ep_Bin] = {Ep_Bin*0};
    double r_pim2e_h3[Ep_Bin] = {Ep_Bin*0};
    double r_pim2e_he3[Ep_Bin] = {Ep_Bin*0};

    DIS* dis = new DIS();
    for(int i=0;i<Ep_Bin;i++){
        dis->SetKin(E0, Ep[i], T0*Deg2Rad);
        Sig_H2[i]  = dis->Sigma(2,1); 
        Sig_H3[i]  =  dis->Sigma(3,1);
        Sig_He3[i] = dis->Sigma(3,2);
        Sig_Al[i]  =  dis->Sigma(27,13);

        cout<<Form("Ep = %4.2f, xbj=%4.2f,  XS_H2=%6.2e, XS_H3=%6.2e, XS_He3=%6.2e, XS_AL=%6.2e", Ep[i], dis->Xbj(), Sig_H2[i],Sig_H3[i],Sig_He3[i],Sig_Al[i])<<endl;
    }
    delete dis;
    /*}}}*/

    //TFile *file = new TFile(Form("BB_rate_marathon_E11_P2_T%d.root",T0),"r");/*{{{*/
    //TTree *T = (TTree*) file->Get("T");
    TChain *T = new TChain("T"); 
    T->Add(Form("BB_rate_marathon_E11_P2_T%d_lowE.root",T0));

    gStyle->SetOptStat(0);
    T->SetLineStyle(1);
    T->SetLineColor(1);
    T->SetMarkerSize(1.2);
    T->SetMarkerColor(3);
    gStyle->SetPaintTextFormat("5.2e"); 

    double N_Total = 0.0;
    double Rate_Total = 0.0;
    double Window_N_Total = 0.0;
    double Window_Rate_Total = 0.0;
    /*}}}*/

    bool b1 = true;
    //bool b1 = false;
    bool b2 = true;
    //bool b2 = false;
    bool b3 = true;
    //bool b3 = false;

    if(b1){
        TCanvas *c1=new TCanvas("c1","c1",1400,800);/*{{{*/
        gPad->SetGridx(1);
        gPad->SetGridy(1);
        TH1F* h1 = new TH1F("h1","",Ep_Bin,Ep_Min,Ep_Max);
        h1->SetXTitle("Ep");
        //h1->SetYTitle("Rate (N/Day)");
        h1->SetYTitle("Rate (Hz)");
        h1->GetXaxis()->CenterTitle(1);
        h1->GetYaxis()->CenterTitle(1);
        h1->GetXaxis()->SetTitleSize(0.07);
        h1->GetYaxis()->SetTitleSize(0.06);
        h1->GetXaxis()->SetTitleOffset(0.6);
        h1->GetYaxis()->SetTitleOffset(0.6);
        h1->GetXaxis()->SetTitleFont(20);
        h1->GetYaxis()->SetTitleFont(20);
        //TH1F* h11 = (TH1F*) h1->Clone();
        TH1F* h11 = new TH1F("h11","",Ep_Bin,Ep_Min,Ep_Max);

        T->Draw("Ep>>h1","(xs_d)*(xb<1.&&W>sqrt(3.))","");
        N_Total = h1->Integral();
        cerr<<"H2: N =  "<<N_Total<<endl;
        Rate_Total = H2_Rate/N_Total*RC_Eff;
        h1->Scale(Rate_Total);
        
        T->Draw("Ep>>h11","(13*xs_p+14*xs_n)*(xb<1.&&W>sqrt(3.))","goff");
        Window_N_Total = h11->Integral();
        Window_Rate_Total = Window_H2_Rate/Window_N_Total*RC_Eff;
        h11->Scale(Window_Rate_Total);

        double total_rate_h2=0.0;
        for(int i=0;i<Ep_Bin;i++){
            double rate1 = h1->GetBinContent(i+1);
            double rate2 = h11->GetBinContent(i+1);
            if(isnan(rate1)) rate1 = 0.0;
            if(isnan(rate2)) rate2 = 0.0;
            rate_pim_h2[i] = rate1*Sig_Pi0[i]/Sig_H2[i] + rate2* Sig_Pi0[i]/Sig_Al[i];
            rate_h2[i] = rate_pim_h2[i]+rate1+rate2;
            r_pim2e_h2[i] = rate_pim_h2[i]/(rate1+rate2);
            if(isnan(rate_pim_h2[i])) rate_pim_h2[i] = 0.0;
            if(isnan(r_pim2e_h2[i])) r_pim2e_h2[i] = 0.0;

            cout<<Form("H2:  Ep = %4.2f, rate1=%4.2f, rate2=%4.2f, Pi/e=%6.2f", 
                    Ep[i],rate1,rate2, r_pim2e_h2[i])<<endl;
            total_rate_h2 += rate_h2[i];
        }

        c1->cd();
        h1->SetMinimum(-1.);
        h1->SetMarkerSize(1.5);
        h1->SetMarkerColor(4);
        gStyle->SetTextFont(32);
        h1->Draw("TEXT45");
        TH1F* h12 = (TH1F*) h1->Clone(); h12->SetLineColor(6);h12->Draw("same");

        TLatex *t1 = new TLatex();
        //t1->SetNDC();
        t1->SetTextFont(22);
        t1->DrawLatex(1.6,h1->GetMaximum()*0.9,Form("^{2}H Rates at %d^{#circ} with BigBite",T0));
        t1->SetTextColor(1);
        t1->SetTextFont(22);
        t1->DrawLatex(1.4,h1->GetMaximum()*0.8,Form("Total Rate(#pi^{0},e^{-}from Target+Windows) = %5.2f KHz",total_rate_h2/1000.));
        t1->SetTextFont(22);
        t1->DrawLatex(1.6,h1->GetMaximum()*0.7,Form("(BB Momentum > 0.6 GeV/c)"));

        t1->SetTextColor(1);
        t1->SetTextFont(22);
         t1->DrawLatex(1.0,h1->GetMaximum()*0.21,Form("Total Rate"));
        t1->SetTextFont(32);
        t1->SetTextColor(2);
        t1->DrawLatex(1.0,h1->GetMaximum()*0.07,Form("#pi^{0}/e^{-} Ratio(Target+Windows)"));
        t1->SetTextColor(1);
        t1->SetTextColor(4);
        t1->DrawLatex(1.0,h1->GetMaximum()*0.14,Form("eDIS Rate from Target"));

        t1->SetTextSize(0.03);
        t1->SetTextAngle(45);
        t1->SetTextColor(1);
        for(int i=0;i<Ep_Bin;i++){
            t1->SetTextColor(1);
            t1->SetTextFont(22);
            t1->DrawLatex(h1->GetBinCenter(i+1)-0.03,h1->GetBinContent(i+1)*1.00,Form("%5.2e",rate_h2[i]));
            t1->SetTextColor(2);
            t1->SetTextFont(32);
            t1->DrawLatex(h1->GetBinCenter(i+1)+0.03,h1->GetBinContent(i+1)*1.00,Form("%5.2f",r_pim2e_h2[i]));
        }
        
        TLatex *t11 = new TLatex();
        t11->SetTextFont(32);
        t11->SetTextSize(0.03);
        t11->SetTextAngle(0);
        t11->SetTextColor(6);
        for(int i=0;i<Ep_Bin;i++){
            t11->DrawLatex(h1->GetBinCenter(i+1)-0.02,-0.5,Form("%4.2f",xbj[i]));
        }
  
        c1->Print(Form("H2_BB_Ep_LowE_%d.png",T0));
        c1->Print(Form("H2_BB_Ep_LowE_%d.pdf",T0));
        /*}}}*/
    }
    if(b2){
        TCanvas *c2=new TCanvas("c2","c2",1400,800);/*{{{*/
        gPad->SetGridx(1);
        gPad->SetGridy(1);
        TH1F* h2 = new TH1F("h2","",Ep_Bin,Ep_Min,Ep_Max);
        h2->SetXTitle("Ep");
        //h2->SetYTitle("Rate (N/Day)");
        h2->SetYTitle("Rate (Hz)");
        h2->GetXaxis()->CenterTitle(1);
        h2->GetYaxis()->CenterTitle(1);
        h2->GetXaxis()->SetTitleSize(0.07);
        h2->GetYaxis()->SetTitleSize(0.06);
        h2->GetXaxis()->SetTitleOffset(0.6);
        h2->GetYaxis()->SetTitleOffset(0.6);
        h2->GetXaxis()->SetTitleFont(20);
        h2->GetYaxis()->SetTitleFont(20);
        //TH1F* h21 = (TH1F*) h2->Clone();
        TH1F* h21 = new TH1F("h21","",Ep_Bin,Ep_Min,Ep_Max);

        T->Draw("Ep>>h2","(xs_d+xs_n)*(xb<1.&&W>sqrt(3.))","");
        N_Total = h2->Integral();
        cerr<<"H3: N =  "<<N_Total<<endl;
        Rate_Total = H3_Rate/N_Total*RC_Eff;
        h2->Scale(Rate_Total);

        T->Draw("Ep>>h21","(13*xs_p+14*xs_n)*(xb<1.&&W>sqrt(3.))","goff");
        Window_N_Total = h21->Integral();
        Window_Rate_Total = Window_H3_Rate/Window_N_Total*RC_Eff;
        h21->Scale(Window_Rate_Total);

        double total_rate_h3=0;
        for(int i=0;i<Ep_Bin;i++){
            double rate1 = h2->GetBinContent(i+1);
            double rate2 = h21->GetBinContent(i+1);
            if(isnan(rate1)) rate1 = 0.0;
            if(isnan(rate2)) rate2 = 0.0;
            rate_pim_h3[i] = rate1*Sig_Pi0[i]/Sig_H3[i] + rate2* Sig_Pi0[i]/Sig_Al[i];
            r_pim2e_h3[i] = rate_pim_h3[i]/(rate1+rate2);
            rate_h3[i] = rate_pim_h3[i]+rate1+rate2;
            if(isnan(rate_pim_h3[i])) rate_pim_h3[i] = 0.0;
            if(isnan(r_pim2e_h3[i])) r_pim2e_h3[i] = 0.0;

            cout<<Form("H3:  Ep = %4.2f, rate1=%4.2f, rate2=%4.2f, Pi/e=%6.2f", 
                    Ep[i],rate1,rate2, r_pim2e_h3[i])<<endl;
            total_rate_h3 += rate_h3[i];
        }

        h2->SetMinimum(-1.);
        h2->SetMarkerSize(1.5);
        h2->SetMarkerColor(4);
        gStyle->SetTextFont(32);
        h2->Draw("TEXT45");
        TH1F* h22 = (TH1F*) h2->Clone(); h22->SetLineColor(6);h22->Draw("same");

        TLatex *t2 = new TLatex();
        //t2->SetNDC();
        t2->SetTextFont(22);
        t2->DrawLatex(1.6,h2->GetMaximum()*0.9,Form("^{3}H Rates at %d^{#circ} with BigBite",T0));
        t2->DrawLatex(1.4,h2->GetMaximum()*0.8,Form("Total Rate(#pi^{0},e^{-}from Target+Windows) = %5.2f KHz",total_rate_h3/1000.));
        t2->SetTextFont(22);
        t2->DrawLatex(1.6,h2->GetMaximum()*0.7,Form("(BB Momentum > 0.6 GeV/c)"));

        t2->SetTextColor(1);
        t2->SetTextFont(22);
 
        t2->DrawLatex(1.0,h2->GetMaximum()*0.21,Form("Total Rate"));
        t2->SetTextColor(2);
        t2->SetTextFont(32);
        t2->DrawLatex(1.0,h2->GetMaximum()*0.07,Form("#pi^{0}/e^{-} Ratio(Target+Windows)"));
        t2->SetTextColor(1);
        t2->SetTextColor(4);
        t2->DrawLatex(1.0,h2->GetMaximum()*0.14,Form("eDIS Rate from Target"));

        t2->SetTextSize(0.03);
        t2->SetTextAngle(45);
        t2->SetTextColor(1);
        t2->SetTextFont(32);
        for(int i=0;i<Ep_Bin;i++){
            t2->SetTextColor(1);
            t2->SetTextFont(22);
            t2->DrawLatex(h2->GetBinCenter(i+1)-0.03,h2->GetBinContent(i+1)*1.00,Form("%5.2e",rate_h3[i]));
            t2->SetTextColor(2);
            t2->SetTextFont(32);
            t2->DrawLatex(h2->GetBinCenter(i+1)+0.03,h2->GetBinContent(i+1)*1.00,Form("(%5.2f)",r_pim2e_h3[i]));
        }
      
        TLatex *t21 = new TLatex();
        t21->SetTextFont(32);
        t21->SetTextSize(0.03);
        t21->SetTextAngle(0);
        t21->SetTextColor(6);
        for(int i=0;i<Ep_Bin;i++){
            t21->DrawLatex(h2->GetBinCenter(i+1)-0.02,-0.5,Form("%4.2f",xbj[i]));
        }
  

        c2->Print(Form("H3_BB_Ep_LowE_%d.png",T0));
        c2->Print(Form("H3_BB_Ep_LowE_%d.pdf",T0));
    }

    if(b3){
        TCanvas *c3=new TCanvas("c3","c3",1400,800);/*{{{*/
        gPad->SetGridx(1);
        gPad->SetGridy(1);
        TH1F* h3 = new TH1F("h3","",Ep_Bin,Ep_Min,Ep_Max);
        h3->SetXTitle("Ep");
        //h3->SetYTitle("Rate (N/Day)");
        h3->SetYTitle("Rate (Hz)");
        h3->GetXaxis()->CenterTitle(1);
        h3->GetYaxis()->CenterTitle(1);
        h3->GetXaxis()->SetTitleSize(0.07);
        h3->GetYaxis()->SetTitleSize(0.06);
        h3->GetXaxis()->SetTitleOffset(0.6);
        h3->GetYaxis()->SetTitleOffset(0.6);
        h3->GetXaxis()->SetTitleFont(20);
        h3->GetYaxis()->SetTitleFont(20);
        //TH1F* h31 = (TH1F*) h3->Clone();
        TH1F* h31 = new TH1F("h31","",Ep_Bin,Ep_Min,Ep_Max);

        T->Draw("Ep>>h3","(xs_d+xs_p)*(xb<1.&&W>sqrt(3.))","");
        N_Total = h3->Integral();
        cerr<<"He3: N =  "<<N_Total<<endl;
        Rate_Total = He3_Rate/N_Total*RC_Eff;
        h3->Scale(Rate_Total);

        T->Draw("Ep>>h31","(13*xs_p+14*xs_n)*(xb<1.&&W>sqrt(3.))","goff");
        Window_N_Total = h31->Integral();
        Window_Rate_Total = Window_He3_Rate/Window_N_Total*RC_Eff;
        h31->Scale(Window_Rate_Total);

        double total_rate_he3=0.0;
        for(int i=0;i<Ep_Bin;i++){
            double rate1 = h3->GetBinContent(i+1);
            double rate2 = h31->GetBinContent(i+1);
            if(isnan(rate1)) rate1 = 0.0;
            if(isnan(rate2)) rate2 = 0.0;
            rate_pim_he3[i] = rate1*Sig_Pi0[i]/Sig_He3[i] + rate2* Sig_Pi0[i]/Sig_Al[i];
            r_pim2e_he3[i] = rate_pim_he3[i]/(rate1+rate2);
            rate_he3[i] = rate_pim_he3[i]+rate1+rate2;
            if(isnan(rate_pim_he3[i])) rate_pim_he3[i] = 0.0;
            if(isnan(r_pim2e_he3[i])) r_pim2e_he3[i] = 0.0;

            cout<<Form("He3:  Ep = %4.2f, rate1=%4.2f, rate2=%4.2f, Pi/e=%6.2f", 
                    Ep[i],rate1,rate2, r_pim2e_he3[i])<<endl;
            total_rate_he3 += rate_he3[i];
        }

        h3->SetMinimum(-1.);
        h3->SetMarkerSize(1.5);
        h3->SetMarkerColor(4);
        gStyle->SetTextFont(32);
        h3->Draw("TEXT45");
        TH1F* h32 = (TH1F*) h3->Clone(); h32->SetLineColor(6);h32->Draw("same");

        TLatex *t3 = new TLatex();
        //t3->SetNDC();
        t3->SetTextFont(22);
        t3->DrawLatex(1.6,h3->GetMaximum()*0.9,Form("^{3}He Rates at %d^{#circ} with BigBite",T0));
        t3->SetTextColor(1);
        t3->DrawLatex(1.4,h3->GetMaximum()*0.8,Form("Total Rate(#pi^{0},e^{-}from Target+Windows) = %5.2f KHz",total_rate_he3/1000.));
        t3->SetTextFont(22);
        t3->DrawLatex(1.6,h3->GetMaximum()*0.7,Form("(BB Momentum > 0.6 GeV/c)"));

        t3->SetTextColor(1);
        t3->SetTextFont(22);
 
        t3->SetTextFont(22);
        t3->SetTextColor(1);
        t3->DrawLatex(1.0,h3->GetMaximum()*0.21,Form("Total Rate"));
        t3->SetTextColor(2);
        t3->SetTextFont(32);
        t3->DrawLatex(1.0,h3->GetMaximum()*0.07,Form("#pi^{0}/e^{-} Ratio"));
        t3->SetTextFont(32);
        t3->SetTextColor(4);
        t3->DrawLatex(1.0,h3->GetMaximum()*0.14,Form("eDIS Rate from Target(Target+Windows)"));

        t3->SetTextSize(0.03);
        t3->SetTextAngle(45);
        t3->SetTextFont(32);
        t3->SetTextColor(1);
        for(int i=0;i<Ep_Bin;i++){
            t3->SetTextColor(1);
            t3->SetTextFont(22);
            t3->DrawLatex(h3->GetBinCenter(i+1)-0.03,h3->GetBinContent(i+1)*1.00,Form("%5.2e",rate_he3[i]));
            t3->SetTextColor(2);
            t3->SetTextFont(32);
            t3->DrawLatex(h3->GetBinCenter(i+1)+0.03,h3->GetBinContent(i+1)*1.00,Form("%5.2f",r_pim2e_he3[i]));
        }

        TLatex *t31 = new TLatex();
        t31->SetTextFont(32);
        t31->SetTextSize(0.03);
        t31->SetTextAngle(0);
        t31->SetTextColor(6);
        for(int i=0;i<Ep_Bin;i++){
            t31->DrawLatex(h3->GetBinCenter(i+1)-0.02,-0.5,Form("%4.2f",xbj[i]));
        }    
        c3->Print(Form("He3_BB_Ep_LowE_%d.png",T0));
        c3->Print(Form("He3_BB_Ep_LowE_%d.pdf",T0));
    }
    }
