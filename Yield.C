//========================================================================
// This code is to calculate the Normalized Charge Yield
// By: Sheren Alsalmi
// Oct,23,2015
//=======================================================================
void Yield()
{
    TCanvas *c = new TCanvas();
    Double_t E1, E2, p, ratio,x, live1, live2, tr_eff , Y, Yin, Ynorm;
    Double_t Ai = 0.0001208, Af= -1.413; // from calibration
    Double_t Q, Av_Q, Q_tot, inI, L_counts, count0, countf, t, LT;
    Double_t PS = 100;
    inI = 11;
    tr_eff = 0.936; // track effeciency
    TFile *ff = TFile::Open("file.root");
    TBranch* counts = T->GetBranch("left_bcm_d3c"); // OR u1c
    TBranch* prl1 = T->GetBranch("L.prl1.e");
    TBranch* prl2 = T->GetBranch("L.prl2.e");
    TBranch* tr_p = T->GetBranch("L.gold.p");
    TBranch* tr_vz = T->GetBranch("L.tr.vz");
    TBranch* clk = T->GetBranch("left_clkcount");
    TBranch* recev = T->GetBranch("DL.scalLive1");
    TBranch* totev = T->GetBranch("DL.scalLive2");
    prl1->SetAddress(&E1);
    prl2->SetAddress(&E2);
    tr_p->SetAddress(&p);
    counts->SetAddress(&L_counts);
    clk->SetAddress(&t);
    recev->SetAddress(&live1);
    totev->SetAddress(&live2);

    Int_t nentries = T->GetEntries();
    
    cout<<"No. of total events="<<nentries<<endl;
    TCut pion_cut1="abs(L.prl1.e+L.prl2.e)>500";
    TCut pion_cut="abs((L.prl1.e)>250 && abs(L.prl2.e)>50)";
    TCut cer_cut="abs(L.cer.asum_c)>250";
    TCut acc_cut="(L.tr.vz)>-0.075 &&(L.tr.vz)<0.075";
 //   TCut clust_cut = "L.vdc.u1.nclust==1 && L.vdc.u2.nclust==1 &&L.vdc.v1.nclust==1 && L.vdc.v2.nclust==1";
     TCut clust_cut = "L.tr.n==1";
     TCut entry_cut =  "Entry$<270000";  // for some runs
     TCut sum_cut=  acc_cut && pion_cut && cer_cut && clust_cut;

 //============= Live Time  ==========================
    TCanvas *c1 = new TCanvas();
    c1->cd();
    T->Draw("DL.scalLive1/DL.scalLive2>>ht(100,0.96,0.97)");
    ht->SetTi
    LT = ht->GetMean();
    cout<<"Live Time = "<<LT<<endl;
//====================================================
    TCanvas *c3 = new TCanvas();
    T->Draw("L.tr.vz>>hvz(100,-0.4,0.4)");
    TCanvas *c5 = new TCanvas();
    c5->Divide(1,2);
    c5->cd(1);
    T->Draw("L.tr.vz>>hvz1(100,-0.2,0.2)");
    c5->cd(2);
    T->Draw("L.tr.vz>>hvz2(100,-0.2,0.2)", acc_cut);

//===================================================
    if (T==0) return;
    for (Int_t j=0; j<=nentries; j++){
        prl1->GetEntry(j);
        prl2->GetEntry(j);
        tr_p->GetEntry(j);
        counts->GetEntry(j);
        clk->GetEntry(j);
        ratio = 0.001*(E1+E2)/p;
        
        if(j==0){
            count0= L_counts;}
    
        if (j==nentries) {

            countf = L_counts;}
    }
    
    
    TCanvas *c2 = new TCanvas();
    c2->cd();
    
    T->Draw("L.tr.n>>hh(1000,0,10)", sum_cut);
    
    Int_t ngood = hh->GetEntries();
    cout<<"no. of good events="<<ngood<<endl;
    
//=== Cer_cut efficiency ===========
    TCanvas *c4 = new TCanvas();
    c4->Divide(1,2);
    c4->cd(1);
    T->Draw("L.cer.asum_c>>hcer1(100,0,3000)", pion_cut);
    c4->cd(2);
    T->Draw("L.cer.asum_c>>hcer2(100,0,3000)", cer_cut);
    Double_t cer_eff = (hcer2->Integral())/(hcer1->Integral());
    cout<<"Cerenkove efficiency = "<<cer_eff<<endl;
//==================================
    
    
    Q_tot = inI * t*1E-5;
    Double_t in_Ne = Q_tot/(1.602*1E-19);
    
    Q = Ai*(countf-count0) + (Af*t*1E-5);
    
    Double_t Iav = Q/(t*1E-5);
    
    Double_t Ne = Q/(1.602*1E-19);
    
    cout<<"Curent = "<<Iav<<endl;
    Y = (ngood*PS)/(Q*LT*tr_eff*cer_eff);
    Yin = 1808.86;;
    Double_t Ynorm = Y/Yin;
    cout<<"Charge Yield="<<Y<<"Normalized Charge Yield="<<Ynorm<<endl;

    


}
