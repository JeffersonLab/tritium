{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TString GeneralCut = "L.tr.n==1 && abs(L.tr.x)<0.9 && abs(L.tr.y)<0.5 && abs(L.tr.th)<0.8 && abs(L.tr.ph)<0.5&&L.tr.p>2.0&&L.tr.p<4.0";
  GeneralCut += "&& abs(ExTgtL.ph)<0.02 && abs(ExTgtL.th)<0.04 && abs(ExTgtL.dp)<0.06 && abs(RctPtL.z)<0.3";
  TString TriggerCut = "((DBB.evtypebits>>3)&1)&&DBB.edtpl[0]==0";
  TString electrons = "epL>=0.5&&L.prl2.e>=100.0&&L.cer.asum_c>=50.0";
  TString cute1 = TriggerCut+"&&" + GeneralCut + "&&"+electrons;

  TString infile; cerr<<"InPut File: "<<endl; cin >> infile;
  ifstream runfile(infile);

  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  TCanvas *c2 = new TCanvas("c2","c2",1200,800);
  c2->Divide(4,4);
 
  TString rootfile; Int_t runnum=0, count=1;
  Float_t kin; Int_t angle;

  TString comment;
  comment.ReadLine(runfile);
  while(!(runfile.eof())){
    
    runfile >> runnum >> kin >> angle;
    rootfile = Form("/w/halla-2/e08014/Rootfiles/e08014_%d.root",runnum);
    cerr << "Working on Run#"<<runnum<<endl;
    
    TFile *file1 = new TFile(rootfile,"r"); 
    TH1F *h1 = new TH1F("h1","C12 Target",200,-0.02,0.02);
    c2->cd(count);
    T->Draw("RctPtL.z>>h1",TCut(cute1));
    Int_t peak1 = h1->GetMaximumBin();
    Double_t scale1 = 10000/h1->Integral(peak1-5,peak1+5);
    h1->Scale(scale1); h1->SetXTitle("ReactPoint.Z");
    
    Float_t z_peak = h1->GetBinCenter(peak1); 
    h1->Fit("gaus","","",z_peak-0.01,z_peak+0.006);
    z_peak = gaus->GetParameter(1);
    cerr <<"Highest Beam Center ="<<z_peak<<endl;
    TString name1 = Form("Kin%2.1f,Run#%d",kin,runnum);
    
    TLatex *tex = new TLatex(0.3,0.8,name1);
    tex->SetNDC();
    tex->SetTextSize(0.08);
    tex->SetTextColor(1);
    tex->Draw();
    TLatex *tex1 = new TLatex(0.3,0.7,Form("Angle=%d,Peak=%5.2f",angle,z_peak*1000));
    tex1->SetNDC();
    tex1->SetTextSize(0.08);
    tex1->SetTextColor(1);
    tex1->Draw();
    
    c1->cd();
  
   if(count ==1)
      h1->Draw();
    else{
      if(count<10)
	h1->SetLineColor(count);
      else if(count==5||count-9==5)
	h1->SetLineColor(13);
      else
	h1->SetLineColor(count-9);
      h1->Draw("same");
    }
    count++;
    
    
  }
}
 //    Int_t peak1 = h1->GetMaximumBin();
 //  Double_t z_he3 = h1->GetBinCenter(peak1); cerr <<"He3 First Beam Center ="<<z_he3<<endl;
 //  Double_t scale1 = 10000/h1->Integral(peak1-5,peak1+5);
 //  h1->Scale(scale1); h1->SetXTitle("ReactPoint.Z");h1->Draw();

 //  c1->cd(2);
 //  TFile *he4 = new TFile(fHe4,"r"); 
 //  TH1F *h3 = new TH1F("h3","He4 Target",300,-0.2,0.2);
  
 //  T->Draw("rpl.z>>h3",TCut(cute1));
 //  peak1 = h3->GetMaximumBin();
 //  Double_t z_he4 = h3->GetBinCenter(peak1); cerr <<"He4 First Beam Center ="<<z_he4<<endl;
 //  scale1 = 10000/h3->Integral(peak1-5,peak1+5);
 //  h3->Scale(scale1); h3->SetXTitle("ReactPoint.Z");h3->Draw();

 //  c1->cd(3);
 //  TFile *dum = new TFile(fDummy,"r"); 
  
 //  if(z_he3>0.)
 //    TH1F *temp = new TH1F("temp","Dummy Target",300, 0.0,0.2);
 //  else
 //    TH1F *temp = new TH1F("temp","Dummy Target",300, -0.2,0.0);

 //  T->Draw("rpl.z>>temp",TCut(cute1));
 //  Int_t peak2 = temp->GetMaximumBin(); 
 //  Double_t z_temp =temp->GetBinCenter(peak2); cerr <<"Dummy First Beam Center ="<<z_temp<<endl;

 //  Float_t shift = (z_he3-z_temp)*1000.0;

 //  TH1F *h2 = new TH1F("h2",Form("Dummy Target with %3.1fmm shift",shift),300,-0.2,0.2);

 //  T->Draw(Form("rpl.z+%f>>h2",z_he3-z_temp),TCut(cute1));
 //  Int_t peak3 = h2->GetMaximumBin(); 
 //  Double_t z_dum =h2->GetBinCenter(peak3); cerr <<"Now Dummy First Beam Center ="<<z_dum<<endl;
 //  Double_t scale2 = 10000/h2->Integral(peak3-5,peak3+5);
 //  h2->Scale(scale2); h2->SetLineColor(2);h2->SetXTitle("ReactPoint.Z");h2->Draw();


 // c1->cd(4);
 // dum->cd();
 //  if(z_he4>0.0)
 //    TH1F *temp1 = new TH1F("temp1","Dummy Target",300, 0.0,0.2);
 //  else
 //    TH1F *temp1 = new TH1F("temp1","Dummy Target",300, -0.2,0.0);

 //  T->Draw("rpl.z>>temp1",TCut(cute1));
 //  peak2 = temp1->GetMaximumBin(); 
 //  Double_t z_temp =temp1->GetBinCenter(peak2); cerr <<"Dummy First Beam Center ="<<z_temp<<endl;

 //  Double_t shift1 = (z_he4-z_temp)*1000.0;

 //  TH1F *h4 = new TH1F("h4",Form("Dummy Target with %3.1fmm shift",shift1),300,-0.2,0.2);

 //  T->Draw(Form("rpl.z+%f>>h4",z_he4-z_temp),TCut(cute1));
 //  peak3 = h4->GetMaximumBin(); 
 //  z_dum =h4->GetBinCenter(peak3); cerr <<"Now Dummy First Beam Center ="<<z_dum<<endl;
 //  scale2 = 10000/h4->Integral(peak3-5,peak3+5);
 //  h4->Scale(scale2); h4->SetLineColor(2);h4->SetXTitle("ReactPoint.Z");h4->Draw();


 //  TCanvas *c2 = new TCanvas("c2","c2",1200,800);
 //  dum.cd(); h2->SetLineColor(4); h2->SetLineWidth(3); h2.Draw();
 //  he3.cd(); h1->SetFillColor(2); h1->SetLineWidth(2);h1.Draw("same");h2.Draw("same");
 

 //  TLatex *tex = new TLatex(0.35,0.8,Form("Kin3.1,Red->He3,Blue->Dummy (shfited - %3.1lfmm)",shift));
 //  //TLatex *tex = new TLatex(0.35,0.8,Form("Kin3.1,Red->H2,Blue->Dummy (shfited - %3.1lfmm)",shift));
 //  tex->SetNDC();
 //  tex->SetTextSize(0.03);tex->SetTextColor(3);
 //  tex->Draw();
 //  c2.Print("Dummy_He3_Kin3.1.png");
 //  //c2.Print("Dummy_H2_Kin3.1.png");
  
 //  TCanvas *c3 = new TCanvas("c3","c3",1200,800);
 //  dum.cd(); h4->SetLineColor(4); h4->SetLineWidth(3); h4.Draw();
 //  he4.cd(); h3->SetFillColor(2); h3->SetLineWidth(2);h3.Draw("same");h4.Draw("same");
  
 //  TLatex *tex1 = new TLatex(0.35,0.8,Form("Kin3.1,Red->He4,Blue->Dummy (shfited - %3.1lfmm)",shift1));
 //  tex1->SetNDC();
 //  tex1->SetTextSize(0.03);tex1->SetTextColor(3);
 //  tex1->Draw();
    //  c3.Print("Dummy_He4_Kin3.1.png");

