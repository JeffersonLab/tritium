#include "GetRate.h"
#include <TH2F.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include "DIS.h"

int main(){
    /*Define{{{*/
    double F2_allm = 0.0, F2_in09 = 0.0, F2_cj = 0.0;
    double F2 = 0.0, F2_Err=0.0; 
    double Q2 = 0.0;
    double xb = 0; 
    double staterr, syserr, dum;
    vector<double> vF2_allm;
    vector<double> vF2_in09;
    vector<double> vF2_cj;
    vector<double> vF2, vF2_Err; 
    vector<double> vXb; 
    vector<double> vQ2; /*}}}*/
    vector<double> vR_allm, vR_in09, vR_cj;
    vector<double> vRErr_allm, vRErr_in09, vRErr_cj;
 
    //TString Nucleon = "Proton";/*{{{*/
    //int A = 1, Z=1;
    //TString DataFileName = "F2p_Data_BCDMS_allE.dat";
    //TString DataFileName = "F2p_Data_SLAC.dat";

    /*TString Nucleon = "Neutron";*/
    //int A = 1, Z=0;
    ////TString DataFileName = "F2n_Data_EMC.dat";
    //TString DataFileName = "F2n_Data_WA25.dat";
    ////TString DataFileName = "F2n_Data_NA2.dat";


    TString Nucleon = "Deutron";
    int A = 2, Z=1;
    TString DataFileName = "F2d_Whitlow.dat";
    //TString DataFileName = "F2d_Data_CLAS.dat";
    //TString DataFileName = "F2d_Data_MNC.dat";
    //TString DataFileName = "F2d_Data_PPE.dat";
     /*}}}*/

    ifstream infile(DataFileName.Data());/*{{{*/
    TString RootFileName = DataFileName;
    RootFileName.ReplaceAll(".dat",".root");

    TFile *file = new TFile(RootFileName.Data(),"recreate");
    TTree* T = new TTree("T", "a new tree");
    T->Branch("xb", &xb, "xb/D");
    T->Branch("Q2", &Q2, "Q2/D");
    T->Branch("F2_in09", &F2_in09, "F2_in09/D");
    T->Branch("F2_cj", &F2_cj, "F2_cj/D");
    T->Branch("F2_allm", &F2_allm, "F2_allm/D");
    T->Branch("F2", &F2, "F2/D");
/*}}}*/
   
    DIS* dis = new DIS();
    dis->SetCJ(100);

    TString com;
    infile >> com >> com >> com >> com >> com;

    ofstream outfile("out_F2d_Whitlow.dat");
    outfile<<Form("%8s %8s %12s %12s %12s %12s %12s", "xb", "Q2", "F2", "staterr", "syserr", "F2_in09", "F2_allm")<<endl;

    while(infile >> xb >> Q2 >> F2 >> staterr >> syserr){/*{{{*/
        //infile >> xb >> Q2 >> F2 >> staterr >> syserr;

        if(Nucleon=="Deutron")
            F2 *=2;//from per-nucleon to per-nucleus

        //if(xb<0.1) continue;

        double W2 = pow(0.938272,2) + Q2/xb - Q2; 

        cout<<Form("For xb = %4.2f, Q2 = %5.3f, F2=%7.5f", xb, Q2,F2)<<endl;

        F2_Err = sqrt( staterr*staterr + syserr*syserr );

        dis->SetF1F2Kin(xb, Q2);
        F2_allm = dis->GetF2("F2ALLM",Nucleon.Data());
        F2_in09 = dis->GetF2("F1F2IN09",Nucleon.Data());
        //F2_cj = dis->GetF2("CJ",Nucleon.Data());

        /*dis->Bosted_F1F209(A, Z, Q2, W2, &dum, &F2_in09);
          F2_allm = dis->GetF2_F2ALLM(xb, Q2, Nucleon.Data());
        //F2_cj = dis->GetF2_CJ(xb, Q2, Nucleon.Data());*/

        T->Fill();
        cout<<Form("  IN09: F2 = %4.3f", F2_in09)<<endl;
        cout<<Form(" ALL2M: F2 = %4.3f", F2_allm)<<endl;
        //cout<<Form("    CJ: F2 = %4.3f", F2_cj)<<endl;
        //
        outfile<<Form("%8.2f %8.2f %12.4f %12.4f %12.4f %12.4f %12.4f", xb, Q2, F2, staterr, syserr, F2_in09, F2_allm)<<endl;

        //if((xb<0.5&&F2_Err/F2<0.1)||(xb>0.5&&xb<0.7&&F2_Err/F2<0.5)||(xb>0.7&&xb<0.8&&F2_Err/F2<1.9)||(xb>0.8&&xb<0.9&&F2_Err/F2<4.0)||(xb>0.9&&xb<1.0&&F2_Err/F2<10.0)){
        if(1){
        vXb.push_back(xb);
        vQ2.push_back(Q2);
        vF2.push_back(F2);
        vF2_Err.push_back(F2_Err);
        //vF2_cj.push_back(F2_cj);
        vF2_in09.push_back(F2_in09);
        vF2_allm.push_back(F2_allm);

        vR_in09.push_back(F2_in09/F2);
        vR_allm.push_back(F2_allm/F2);
        //vR_cj.push_back(F2_cj/F2);
        vRErr_in09.push_back(F2_in09/F2* F2_Err); 
        vRErr_allm.push_back(F2_allm/F2* F2_Err);
        //vRErr_cj.push_back(F2_cj/F2* F2_Err/F2);
        }
    }/*}}}*/
    file->cd(); T->Write(); file->Close();
    delete dis;
    infile.close();
    outfile.close();

    gStyle->SetOptStat(0);/*{{{*/
    TCanvas *c1 = new TCanvas("c1","c1",1200,800);
    c1->Divide(1,2);
    c1->cd(1);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    TH2F* h1 = new TH2F("h1","",100,0.0,1.0,100,-0.01,0.7);
    h1->SetXTitle("xb");
    h1->SetYTitle("F2");
    h1->GetXaxis()->SetTitleSize(0.09);
    h1->GetXaxis()->SetTitleOffset(0.7);
    h1->GetXaxis()->CenterTitle(1);
    h1->GetYaxis()->SetTitleSize(0.09);
    h1->GetYaxis()->SetTitleOffset(0.5);
     h1->GetYaxis()->CenterTitle(1);
    h1->GetXaxis()->SetLabelSize(0.08);
    h1->GetYaxis()->SetLabelSize(0.08);
     h1->Draw();

    TGraphErrors* gr1 = new TGraphErrors(vXb.size(), &vXb[0], &vF2[0], 0, &vF2_Err[0]);
    gr1->SetMarkerStyle(24);
    gr1->SetMarkerColor(1);
    gr1->Draw("p");

    TGraph* gr2 = new TGraph(vXb.size(), &vXb[0], &vF2_in09[0]);
    gr2->SetMarkerStyle(25);
    gr2->SetMarkerColor(2);
     gr2->Draw("p");

    TGraph* gr3 = new TGraph(vXb.size(), &vXb[0], &vF2_allm[0]);
    gr3->SetMarkerStyle(26);
    gr3->SetMarkerColor(4);
    gr3->Draw("p");


    //TGraph* gr4 = new TGraph(vXb.size(), &vXb[0], &vF2_cj[0]);
    //gr4->SetMarkerStyle(27);
    //gr4->SetMarkerColor(6);
     //gr4->Draw("p");

    TLegend *l1 = new TLegend(0.80,0.30,0.95,0.50,"F2","");
    l1->SetTextSize(0.04);
    l1->SetTextFont(32);
    l1->AddEntry(gr1,"Data","p");
    l1->AddEntry(gr2,"F1F2IN09","p");
    l1->AddEntry(gr3,"ALLM97","p");
    //l1->AddEntry(gr4,"CJ","p");
    l1->Draw();

    c1->cd(2);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    TH2F* h2 = new TH2F("h2","",100,0.0,1.0,100,0.0,2.);
    h2->SetXTitle("xb");
    h2->SetYTitle("F2_{model} / F2_{data}");
    h2->GetXaxis()->CenterTitle(1);
    h2->GetYaxis()->CenterTitle(1);
    h2->GetXaxis()->SetTitleSize(0.09);
    h2->GetXaxis()->SetTitleOffset(0.7);
    h2->GetYaxis()->SetTitleSize(0.09);
    h2->GetYaxis()->SetTitleOffset(0.5);
    h2->GetXaxis()->SetLabelSize(0.08);
    h2->GetYaxis()->SetLabelSize(0.08);
 h2->Draw();

    TGraphErrors* rr1 = new TGraphErrors(vXb.size(), &vXb[0], &vR_in09[0], 0, &vRErr_in09[0]);
    rr1->SetMarkerStyle(25);
    rr1->SetLineColor(2);
    rr1->SetMarkerColor(2);
    rr1->Draw("p");

    TGraphErrors* rr2 = new TGraphErrors(vXb.size(), &vXb[0], &vR_allm[0], 0, &vRErr_allm[0]);
    rr2->SetMarkerStyle(26);
    rr2->SetMarkerColor(4);
    rr2->SetLineColor(4);
    rr2->Draw("p");

    //TGraphErrors* rr3 = new TGraphErrors(vXb.size(), &vXb[0], &vR_cj[0], 0, &vRErr_cj[0]);
    //rr3->SetMarkerStyle(22);
    //rr3->SetMarkerColor(6);
    //rr3->Draw("p");

    TString PNGFileName = DataFileName;
    PNGFileName.ReplaceAll(".dat",".png");
    c1->Print(PNGFileName.Data());
    PNGFileName.ReplaceAll(".png",".pdf");
    c1->Print(PNGFileName.Data());
/*}}}*/
}
