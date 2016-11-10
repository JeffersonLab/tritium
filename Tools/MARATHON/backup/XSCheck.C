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
    double XS_allm = 0.0, XS_in09 = 0.0, XS_cj = 0.0;
    double XS = 0.0, XS_Err=0.0; 
    double Q2 = 0.0;
    double xb = 0; 
    double staterr, syserr, dum;
    vector<double> vXS_allm;
    vector<double> vXS_in09;
    vector<double> vXS_cj;
    vector<double> vXS, vXS_Err; 
    vector<double> vXb; 
    vector<double> vQ2; /*}}}*/
    vector<double> vR_allm, vR_in09, vR_cj;
    vector<double> vRErr_allm, vRErr_in09, vRErr_cj;
 

    TString Nucleon = "Deutron";
    int A = 2, Z=1;
    TString DataFileName = "XS_Whitlow.dat";

    double E0, Ep, Theta, epsilon, RC,W2,W2_data;
    int i,j;

    ifstream infile(DataFileName.Data());/*{{{*/
    TString RootFileName = DataFileName;
    RootFileName.ReplaceAll(".dat",".root");

    TFile *file = new TFile(RootFileName.Data(),"recreate");
    TTree* T = new TTree("T", "a new tree");
    T->Branch("E0", &E0, "E0/D");
    T->Branch("Ep", &Ep, "Ep/D");
    T->Branch("Theta", &Theta, "Theta/D");
    T->Branch("xb", &xb, "xb/D");
    T->Branch("Q2", &Q2, "Q2/D");
    T->Branch("XS_in09", &XS_in09, "XS_in09/D");
    T->Branch("XS_cj", &XS_cj, "XS_cj/D");
    T->Branch("XS_allm", &XS_allm, "XS_allm/D");
    T->Branch("XS", &XS, "XS/D");
/*}}}*/
   
    DIS* dis = new DIS();
    dis->SetCJ(100);

    TString com;
    infile >> com >> com >> com >> com >> com  
        >> com >> com >> com >> com >> com
        >> com >> com >> com;

    ofstream outfile("out_XS_Whitlow.dat");
    outfile<<Form("%8s %8s %8s %8s %8s %12s %12s %12s %12s %12s", "E0", "Ep","Theta","xb", "Q2", "XS", "staterr", "syserr", "XS_in09", "XS_allm")<<endl;

    double min = 10000000,max=-1000000;
    while(infile >>i>>j>>E0>>Ep>>Theta>> xb >> Q2 >> epsilon>>W2_data>>RC>>XS >> staterr >> syserr){/*{{{*/
        //infile >> xb >> Q2 >> XS >> staterr >> syserr;

        if(Nucleon=="Deutron")
            XS *=2/1000.;//from per-nucleon to per-nucleus, orginal unit was in pb/sr/GeV

        //if(xb<0.1) continue;

        W2 = pow(0.938272,2) + Q2/xb - Q2; 

        cout<<Form("%d %d, For E0=%4.2f, Ep=%4.2f, T=%4.2f, xb = %4.2f, Q2 = %5.3f, W2=%5.3f/%5.3f, XS=%7.5f",i,j,E0,Ep,Theta, xb, Q2,W2_data,W2,XS)<<endl;

        XS_Err = sqrt( staterr*staterr + syserr*syserr );

        //dis->SetKin(E0,Ep,Theta*3.1415926/180.);
        //XS_allm = dis->Sigma("F2ALLM",Nucleon.Data());
        //XS_in09 = dis->Sigma("F1F2IN09",Nucleon.Data());

        double SinSQ = pow( sin(Theta*3.1415926/180. * 0.5), 2);
        double CosSQ = pow( cos(Theta*3.1415926/180. * 0.5), 2);
        double Nu = E0 - Ep;

        dis->SetF1F2Kin(xb, Q2);
        double F2_allm = dis->GetF2("F2ALLM",Nucleon.Data());
        double F2_in09 = dis->GetF2("F1F2IN09",Nucleon.Data());
  
        XS_allm = 4.0 * pow(1/137.036, 2) * Ep*Ep / Q2/Q2 
                    * (1./Nu * CosSQ + 1/0.938272/xb * SinSQ) * F2_allm*GeV2_to_nbarn;
 
        XS_in09 = 4.0 * pow(1/137.036, 2) * Ep*Ep / Q2/Q2 
                    * (1./Nu * CosSQ + 1/0.938272/xb * SinSQ) * F2_in09*GeV2_to_nbarn;
 
 
        T->Fill();
        cout<<Form("  IN09: XS = %4.3f", XS_in09)<<endl;
        cout<<Form(" ALL2M: XS = %4.3f", XS_allm)<<endl;
        //cout<<Form("    CJ: XS = %4.3f", XS_cj)<<endl;
        //
        outfile<<Form("%8.2f %8.2f %12.4f %12.4f %12.4f %12.4f %12.4f", xb, Q2, XS, staterr, syserr, XS_in09, XS_allm)<<endl;

        //if((xb<0.5&&XS_Err/XS<0.1)||(xb>0.5&&xb<0.7&&XS_Err/XS<0.5)||(xb>0.7&&xb<0.8&&XS_Err/XS<1.9)||(xb>0.8&&xb<0.9&&XS_Err/XS<4.0)||(xb>0.9&&xb<1.0&&XS_Err/XS<10.0)){
        if(1){
        vXb.push_back(xb);
        vQ2.push_back(Q2);
        vXS.push_back(XS);
        vXS_Err.push_back(XS_Err);
        //vXS_cj.push_back(XS_cj);
        vXS_in09.push_back(XS_in09);
        vXS_allm.push_back(XS_allm);

        vR_in09.push_back(XS_in09/XS);
        vR_allm.push_back(XS_allm/XS);
        //vR_cj.push_back(XS_cj/XS);
        vRErr_in09.push_back(XS_in09/XS* XS_Err); 
        vRErr_allm.push_back(XS_allm/XS* XS_Err);
        //vRErr_cj.push_back(XS_cj/XS* XS_Err/XS);
          
        if(XS>max) max=XS;
        if(XS<min) min=XS;
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
    gPad->SetLogy(1);
    TH2F* h1 = new TH2F("h1","",100,0.0,1.0,100,min*0.5,max*2);
    h1->SetXTitle("xb");
    h1->SetYTitle("XS");
    h1->GetXaxis()->SetTitleSize(0.09);
    h1->GetXaxis()->SetTitleOffset(0.7);
    h1->GetXaxis()->CenterTitle(1);
    h1->GetYaxis()->SetTitleSize(0.09);
    h1->GetYaxis()->SetTitleOffset(0.5);
     h1->GetYaxis()->CenterTitle(1);
    h1->GetXaxis()->SetLabelSize(0.08);
    h1->GetYaxis()->SetLabelSize(0.08);
     h1->Draw();

    TGraphErrors* gr1 = new TGraphErrors(vXb.size(), &vXb[0], &vXS[0], 0, &vXS_Err[0]);
    gr1->SetMarkerStyle(24);
    gr1->SetMarkerColor(1);
    gr1->Draw("p");

    TGraph* gr2 = new TGraph(vXb.size(), &vXb[0], &vXS_in09[0]);
    gr2->SetMarkerStyle(25);
    gr2->SetMarkerColor(2);
     gr2->Draw("p");

    TGraph* gr3 = new TGraph(vXb.size(), &vXb[0], &vXS_allm[0]);
    gr3->SetMarkerStyle(26);
    gr3->SetMarkerColor(4);
    gr3->Draw("p");


    //TGraph* gr4 = new TGraph(vXb.size(), &vXb[0], &vXS_cj[0]);
    //gr4->SetMarkerStyle(27);
    //gr4->SetMarkerColor(6);
     //gr4->Draw("p");

    c1->cd(2);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    TH2F* h2 = new TH2F("h2","",100,0.0,1.0,100,0.0,2.);
    h2->SetXTitle("xb");
    h2->SetYTitle("XS_{model} / XS_{data}");
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

    TLegend *l1 = new TLegend(0.80,1.20,0.95,1.70,"XS","");
    l1->SetTextSize(0.04);
    l1->SetTextFont(32);
    l1->AddEntry(gr1,"Data","p");
    l1->AddEntry(gr2,"F1XSIN09","p");
    l1->AddEntry(gr3,"ALLM97","p");
    //l1->AddEntry(gr4,"CJ","p");
    l1->Draw();


    TString PNGFileName = DataFileName;
    PNGFileName.ReplaceAll(".dat",".png");
    c1->Print(PNGFileName.Data());
    PNGFileName.ReplaceAll(".png",".pdf");
    c1->Print(PNGFileName.Data());
/*}}}*/
}
