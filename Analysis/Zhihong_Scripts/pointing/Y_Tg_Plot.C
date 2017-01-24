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
    rootfile = Form("/work/halla/e08014/disk1/Rootfiles/e08014_%d.root",runnum);
    cerr << "Working on Run#"<<runnum<<endl;
    
    TFile *file1 = new TFile(rootfile,"r"); 
    TH1F *h1 = new TH1F("h1","Dummy Y_Tg",200,-0.0,0.1);
    c2->cd(count);
    T->Draw("ExTgtL.y>>h1",TCut(cute1));
    Int_t peak1 = h1->GetMaximumBin();
    Double_t scale1 = 10000/h1->Integral(peak1-5,peak1+5);
    h1->Scale(scale1); h1->SetXTitle("ExTgt.y");
    
    Float_t z_peak = h1->GetBinCenter(peak1); 
    h1->Fit("gaus","Q","",z_peak-0.01,z_peak+0.006);
    z_peak = gaus->GetParameter(1);
    cerr <<"Highest Beam Center ="<<z_peak<<endl;
    TString name1 = Form("Kin%2.1f,Run#%d",kin,runnum);
    
    TLatex *tex = new TLatex(0.3,0.8,name1);
    tex->SetNDC();
    tex->SetTextSize(0.08);
    tex->SetTextColor(1);
    tex->Draw();
    TLatex *tex1 = new TLatex(0.2,0.7,Form("Angle=%d,Peak=%5.2fmm",angle,z_peak*1000));
    tex1->SetNDC();
    tex1->SetTextSize(0.08);
    tex1->SetTextColor(1);
    tex1->Draw();
    
    c1->cd();
  
   if(count ==1)
      h1->Draw();
    else{
      if(count<10){
	h1->SetLineColor(count);
	gaus->SetLineColor(count);
      }
      else if(count==5||count-9==5){
	h1->SetLineColor(13);
      	gaus->SetLineColor(13);
      }
      else{
	h1->SetLineColor(count-9);
      	gaus->SetLineColor(count-9);
      }
      h1->Draw("same");
    }
    count++;
    
    
  }
}
 
