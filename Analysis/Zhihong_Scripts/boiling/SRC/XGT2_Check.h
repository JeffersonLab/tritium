  //////////////////////////
  // Check Data and Cuts
  //////////////////////////
  
  /*inline void gCheck_Data_Cuts(TTree *aT_Tree,const TString& aArm,const TString& aCuts,Bool_t aIsExtTgt=kFALSE){{{*/
  void gCheck_Data_Cuts(const vector<Int_t>& aRunNoChain,const TString aArm,const TString aCuts,Bool_t aIsExtTgt=kFALSE)
  {
    Int_t Main_Trigger = -1;
    if(aArm=="L"){
      Main_Trigger = 3;
    }
    else if(aArm=="R"){
      Main_Trigger = 1;
    }
    
    TCanvas *c1 = new TCanvas("c1",Form("PID for HRS-%s",aArm.Data()), 1200,800);
    c1->Divide(2,2);
    
    TCanvas *c2 = new TCanvas("c2",Form("Focal Plane for HRS-%s",aArm.Data()), 1200,800);
    c2->Divide(2,2);
    
    TCanvas *c3 = new TCanvas("c3",Form("Target Plane for HRS-%s",aArm.Data()), 1200,800);
    c3->Divide(2,2);
    
    
    for ( unsigned int i=0; i<aRunNoChain.size(); i++ )
    {
      Int_t aRunNo=aRunNoChain[i];
      TChain* aT_Tree=new TChain("T");  //Aussming the first run get the most correct EPICS value
      aT_Tree->Add(Form("%s/%s_%d.root",ROOTFILES_DIR.Data(),ROOTFILES_NAME.Data(),aRunNo));
      
      aT_Tree->SetLineColor(1); aT_Tree->SetLineWidth(1);
      
      c1->cd(1); c1->Update();
      TH1F *cer = new TH1F("cer","Cherenkov Sum", 300,-10.,4000.);
      gPad->SetLogy(1);
      aT_Tree->Draw(Form("%s.cer.asum_c>>cer",aArm.Data()));
      aT_Tree->SetLineColor(3); aT_Tree->SetLineWidth(3);
      aT_Tree->Draw(Form("%s.cer.asum_c",aArm.Data()),aCuts.Data(),"same");
      
      c1->cd(2);c1->Update();
      TH1F *ep = new TH1F("ep","Calorimeter E/P", 300,-0.3,1.5);
      gPad->SetLogy(1); aT_Tree->SetLineColor(1); aT_Tree->SetLineWidth(1);
      aT_Tree->Draw(Form("ep%s>>ep",aArm.Data()),Form("DBB.evtypebits>>%d&1",Main_Trigger));
      aT_Tree->SetLineColor(3); aT_Tree->SetLineWidth(3);
      aT_Tree->Draw(Form("ep%s",aArm.Data()),aCuts.Data(),"same");
      
      c1->cd(3);c1->Update();
      TH2F *calo = new TH2F("calo","Calorimeter 1.vs.2", 100,-10.0,3500.0,100,-10.0,3500.0);
      gPad->SetLogy(0); aT_Tree->SetMarkerColor(2);
      if(aArm=="L"){
	aT_Tree->Draw("L.prl2.e:L.prl1.e>>calo");
	aT_Tree->SetMarkerColor(4);
	aT_Tree->Draw("L.prl2.e:L.prl1.e",aCuts.Data(),"same");
      }
      else if(aArm=="R"){
	aT_Tree->Draw("R.sh.e:R.ps.e>>calo");
	aT_Tree->SetMarkerColor(4);
	aT_Tree->Draw("R.sh.e:R.ps.e",aCuts.Data(),"same");
      }
      
      c1->cd(4);c1->Update();
      TH2F *calo1 = new TH2F("calo1","Calorimeter 1.vs.2", 100,-10.0,3500.0,100,-10.0,3500.0);
      gPad->SetLogy(0); aT_Tree->SetMarkerColor(1);
      if(aArm=="L")
	aT_Tree->Draw("L.prl2.e:L.prl1.e>>calo1",aCuts.Data(),"colz");
      else if(aArm=="R")
	aT_Tree->Draw("R.sh.e:R.ps.e>>calo1",aCuts.Data(),"colz");
      
      c1->Print(Form("%d_%s_PID_Check.png",aRunNo,aArm.Data()));
      
      //====================================================
      aT_Tree->SetLineColor(1); aT_Tree->SetLineWidth(1);
      c2->Update();c2->cd(1);
      TH1F *xfp = new TH1F("xfp","X_Fp",300,-1.0,1.0);
      aT_Tree->Draw(Form("%s.tr.r_x>>xfp",aArm.Data()));
      aT_Tree->SetLineColor(4); aT_Tree->SetLineWidth(3);
      aT_Tree->Draw(Form("%s.tr.r_x",aArm.Data()),aCuts.Data(),"same");
      
      c2->cd(2);c2->Update();
      TH1F *yfp = new TH1F("yfp","Y_Fp",300,-0.1,0.1);
      aT_Tree->SetLineColor(1); aT_Tree->SetLineWidth(1);
      aT_Tree->Draw(Form("%s.tr.r_y>>yfp",aArm.Data()));
      aT_Tree->SetLineColor(4); aT_Tree->SetLineWidth(3);
      aT_Tree->Draw(Form("%s.tr.r_y",aArm.Data()),aCuts.Data(),"same");
      
      c2->cd(3);c2->Update();
      TH1F *thfp = new TH1F("thfp","Theta_Fp",300,-0.06,0.06);
      aT_Tree->SetLineColor(1); aT_Tree->SetLineWidth(1);
      aT_Tree->Draw(Form("%s.tr.r_th>>thfp",aArm.Data()));
      aT_Tree->SetLineColor(4); aT_Tree->SetLineWidth(3);
      aT_Tree->Draw(Form("%s.tr.r_th",aArm.Data()),aCuts.Data(),"same");
      
      c2->cd(4);  c2->Update();
      TH1F *phfp = new TH1F("phfp","Phi_Fp",300,-0.06,0.06);
      aT_Tree->SetLineColor(1); aT_Tree->SetLineWidth(1);
      aT_Tree->Draw(Form("%s.tr.r_ph>>phfp",aArm.Data()));
      aT_Tree->SetLineColor(4); aT_Tree->SetLineWidth(3);
      aT_Tree->Draw(Form("%s.tr.r_ph",aArm.Data()),aCuts.Data(),"same");
      
      c2->Print(Form("%d_%s_Fp_Check.png",aRunNo,aArm.Data()));
      
      //========================
      
      c3->cd(1); c3->Update();
      TH1F *dptg = new TH1F("dptg","Dp_Tg",300,-0.12,0.12);
      aT_Tree->Draw(Form("%s.gold.dp>>dptg",aArm.Data()));
      aT_Tree->SetLineColor(4); aT_Tree->SetLineWidth(3);
      aT_Tree->Draw(Form("%s.gold.dp",aArm.Data()),aCuts.Data(),"same");
      
      c3->cd(2);c3->Update();
      if(aIsExtTgt)
  	 TH1F *ytg = new TH1F("ytg","Y_Tg",300,-0.1,0.1);
      else
	 TH1F *ytg = new TH1F("ytg","Y_Tg",300,-0.015,0.015);
      
      aT_Tree->SetLineColor(1); aT_Tree->SetLineWidth(1);
      aT_Tree->Draw(Form("%s.gold.y>>ytg",aArm.Data()));
      aT_Tree->SetLineColor(4); aT_Tree->SetLineWidth(3);
      aT_Tree->Draw(Form("%s.gold.y",aArm.Data()),aCuts.Data(),"same");
      
      c3->cd(3); c3->Update();
      TH1F *thtg = new TH1F("thtg","Theta_Tg",300,-0.1,0.1);
      aT_Tree->SetLineColor(1); aT_Tree->SetLineWidth(1);
      aT_Tree->Draw(Form("%s.gold.th>>thtg",aArm.Data()));
      aT_Tree->SetLineColor(4); aT_Tree->SetLineWidth(3);
      aT_Tree->Draw(Form("%s.gold.th",aArm.Data()),aCuts.Data(),"same");
      
      c3->Update();c3->cd(4);  
      TH1F *phtg = new TH1F("phtg","Phi_Tg",300,-0.06,0.06);
      aT_Tree->SetLineColor(1); aT_Tree->SetLineWidth(1);
      aT_Tree->Draw(Form("%s.gold.ph>>phtg",aArm.Data()));
      aT_Tree->SetLineColor(4); aT_Tree->SetLineWidth(3);
      aT_Tree->Draw(Form("%s.gold.ph",aArm.Data()),aCuts.Data(),"same");
      
      c3->Print(Form("%d_%s_Tgt_Check.png",aRunNo,aArm.Data()));
      
      delete aT_Tree;
    }
  }
  /*}}}*/
  
