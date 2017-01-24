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

/*Make Plots{{{*/
inline void gMake_Plot(TTree *aT_EX,const TString& aEX_Cut,TTree *aT_MC,const Double_t aAccCut,const TString& aTarget,const TString& aArm,const TString& aKin,const Double_t aP0,const Double_t aAngle, const Double_t aNorm ){

	/*HRS-L{{{*/
	if(aArm=="L"){
		gStyle->SetOptStat(0);

		const double aE0 = 3.356;

		/*Ploting EX Data{{{*/
		/*Focal Plane{{{*/
		//X_fp
		TH1F *h_xfp = new TH1F("h_xfp","", 100, -1.15,1.15);
		h_xfp->SetXTitle("x_{fp}");
		h_xfp->SetLineColor(2); 
		h_xfp->SetLineWidth(1.1);
		h_xfp->GetXaxis()->SetTitleSize(0.08);
		h_xfp->GetXaxis()->SetTitleOffset(0.7);
		h_xfp->GetXaxis()->CenterTitle(1);
		T_EX->Project("h_xfp","L.tr.r_x",Cut_All);
		scale = norm/h_xfp->Integral(); h_xfp->Scale(scale);
		//Y_fp
		TH1F *h_yfp = new TH1F("h_yfp","", 100, -0.05,0.05);
		h_yfp->SetXTitle("y_{fp}");
		h_yfp->SetLineColor(2); 
		h_yfp->SetLineWidth(1.1);
		h_yfp->GetXaxis()->SetTitleSize(0.08);
		h_yfp->GetXaxis()->SetTitleOffset(0.7);
		h_yfp->GetXaxis()->CenterTitle(1);
		T_EX->Project("h_yfp","L.tr.r_y",Cut_All);
		scale = norm/h_yfp->Integral(); h_yfp->Scale(scale);
		//Th_fp
		TH1F *h_tfp = new TH1F("h_tfp","", 100, -0.04,0.04);
		h_tfp->SetXTitle("#theta_{fp}");
		h_tfp->SetLineColor(2); 
		h_tfp->SetLineWidth(1.1);
		h_tfp->GetXaxis()->SetTitleSize(0.08);
		h_tfp->GetXaxis()->SetTitleOffset(0.7);
		h_tfp->GetXaxis()->CenterTitle(1);
		T_EX->Project("h_tfp","L.tr.r_th",Cut_All);
		scale = norm/h_tfp->Integral(); h_tfp->Scale(scale);
		//Ph_fp
		TH1F *h_pfp = new TH1F("h_pfp","", 100, -0.06,0.06);
		h_pfp->SetXTitle("#phi_{fp}");
		h_pfp->SetLineColor(2); 
		h_pfp->SetLineWidth(1.1);
		h_pfp->GetXaxis()->SetTitleSize(0.08);
		h_pfp->GetXaxis()->SetTitleOffset(0.7);
		h_pfp->GetXaxis()->CenterTitle(1);
		T_EX->Project("h_pfp","L.tr.r_ph",Cut_All);
		scale = norm/h_pfp->Integral(); h_pfp->Scale(scale);
		/*}}}*/

		/*Target Plane{{{*/
		//DeltaP
		TH1F *h_dtg = new TH1F("h_dtg","", 100, -0.07,0.07);
		h_dtg->SetXTitle("#delta P");
		h_dtg->SetLineColor(2); 
		h_dtg->SetLineWidth(1.1);
		h_dtg->GetXaxis()->SetTitleSize(0.08);
		h_dtg->GetXaxis()->SetTitleOffset(0.7);
		h_dtg->GetXaxis()->CenterTitle(1);
		T_EX->Project("h_dtg","ExTgtL.dp",Cut_All);
		scale = norm/h_dtg->Integral(); h_dtg->Scale(scale);
		//Y_tg  
		TH1F *h_ytg = new TH1F("h_ytg","", 100, -0.04,0.04);
		h_ytg->SetXTitle("y_{tg}");
		h_ytg->SetLineColor(2); 
		h_ytg->SetLineWidth(1.1);
		h_ytg->GetXaxis()->SetTitleSize(0.08);
		h_ytg->GetXaxis()->SetTitleOffset(0.7);
		h_ytg->GetXaxis()->CenterTitle(1);
		T_EX->Project("h_ytg","ExTgtL.y",Cut_All);
		scale = norm/h_ytg->Integral(); h_ytg->Scale(scale);
		//Z_react
		TH1F *h_ztg = new TH1F("h_ztg","", 100, -0.12,0.12);
		h_ztg->SetXTitle("z_{react}");
		h_ztg->SetLineColor(2); 
		h_ztg->SetLineWidth(1.1);
		h_ztg->GetXaxis()->SetTitleSize(0.08);
		h_ztg->GetXaxis()->SetTitleOffset(0.7);
		h_ztg->GetXaxis()->CenterTitle(1);
		T_EX->Project("h_ztg","RctPtL.z",Cut_All);
		scale = norm/h_ztg->Integral(); h_ztg->Scale(scale);
		//Th_tg  
		TH1F *h_ttg = new TH1F("h_ttg","", 100, -0.07,0.07);
		h_ttg->SetXTitle("#theta_{tg}");
		h_ttg->SetLineColor(2); 
		h_ttg->SetLineWidth(1.1);
		h_ttg->GetXaxis()->SetTitleSize(0.08);
		h_ttg->GetXaxis()->SetTitleOffset(0.7);
		h_ttg->GetXaxis()->CenterTitle(1);
		T_EX->Project("h_ttg","ExTgtL.th",Cut_All);
		scale = norm/h_ttg->Integral(); h_ttg->Scale(scale);
		//Ph_tg  
		TH1F *h_ptg = new TH1F("h_ptg","", 100, -0.05,0.05);
		h_ptg->SetXTitle("#phi_{tg}");
		h_ptg->SetLineColor(2); 
		h_ptg->SetLineWidth(1.1);
		h_ptg->GetXaxis()->SetTitleSize(0.08);
		h_ptg->GetXaxis()->SetTitleOffset(0.7);
		h_ptg->GetXaxis()->CenterTitle(1);
		T_EX->Project("h_ptg",Form("(ExTgtL.ph-%f)",aL_Ph_Offset),Cut_All);
		scale = norm/h_ptg->Integral(); h_ptg->Scale(scale);

		/*Define Acceptance Variables{{{*/
		double aSinA = sin(aAngle*TMath::DegToRad());
		double aCosA = cos(aAngle*TMath::DegToRad());
		double aSinSQ = pow(sin(aAngle*TMath::DegToRad()/2.0),2); //Define Sin(Theta/2.)^2

		TString aEX_Dp= Form("ExTgt%s.dp",aArm.Data());
		TString aEX_Th = Form("ExTgt%s.th",aArm.Data());
		TString aEX_Ph = Form("(ExTgt%s.ph-%f)",aArm.Data(),aL_Ph_Offset);
		TString aEX_SA = Form("acos((%f-%s*%f)/sqrt(1.0+pow(%s,2)+pow(%s,2)))",aCosA,aEX_Ph.Data(),aSinA,aEX_Th.Data(),aEX_Ph.Data());
		TString aEX_SinSQ = Form("pow(sin(%s/2.0),2)",aEX_SA.Data()); //Define Sin(Theta/2.)^2

		TString aEX_Ep =Form("(%7.5f*(1.0+%s))",aP0,aEX_Dp.Data());
		TString aEX_Nu= Form("(%7.5f-%7.5f*(1.0+%s))",aE0,aP0,aEX_Dp.Data());
		TString aEX_Q2 = Form("4.0*%7.5f*(%7.5f*(1.0+%s))*%s",aE0,aP0,aEX_Dp.Data(),aEX_SinSQ.Data());
		TString aEX_Xbj=Form("%s/(2.0*%8.6f*%s)",aEX_Q2.Data(),PROTON_MASS,aEX_Nu.Data());
		/*}}}*/

		cerr<<"---Xbj = "<<aEX_Xbj.Data()<<endl;

		//xbj   
		double Xbj_Max = 4.0;
		double Xbj_Min = 0.8;
		if(Kin.Contains("5.0")){
			Xbj_Min = 0.5;
			Xbj_Max = 1.5;
		}
		if(Kin.Contains("5.05")){
			Xbj_Min = 0.8;
			Xbj_Max = 1.8;
		}
		if(Kin.Contains("5.0")){
			Xbj_Min = 0.5;
			Xbj_Max = 1.5;
		}
		if(Kin.Contains("3.1")||Kin.Contains("4.1")||Kin.Contains("5.1")){
			Xbj_Min = 1.0;
			Xbj_Max = 2.5;
		}
		if(Kin.Contains("3.2")||Kin.Contains("4.2")||Kin.Contains("5.2")||Kin.Contains("6.5")){
			Xbj_Min = 1.5;
			Xbj_Max = 4.0;
			if(Target=="He3")
				Xbj_Max = 3.0;
		}

		TH1F *h_xbj = new TH1F("h_xbj","", 400, Xbj_Min,Xbj_Max);
		h_xbj->SetXTitle("x_{bj}");
		h_xbj->SetLineColor(2); 
		h_xbj->SetLineWidth(1.1);
		h_xbj->GetXaxis()->SetTitleSize(0.08);
		h_xbj->GetXaxis()->SetTitleOffset(0.7);
		h_xbj->GetXaxis()->CenterTitle(1);
		//aEX_Xbj=Form("ExtEK%s.x_bj",aArm.Data());
		T_EX->Project("h_xbj",aEX_Xbj.Data(),Cut_All);
		scale = norm/h_xbj->Integral(); h_xbj->Scale(scale);
		/*}}}*/

		/*}}}*/

		/*Ploting MC Data{{{*/

		/*Focal Plane{{{*/
		//X_fp
		TH1F *m_xfp = new TH1F("m_xfp","", 100,-1.15,1.15);
		m_xfp->SetXTitle("x_{fp}");
		m_xfp->SetLineColor(4); 
		m_xfp->SetLineWidth(1.1);
		m_xfp->GetXaxis()->SetTitleSize(0.08);
		m_xfp->GetXaxis()->SetTitleOffset(0.7);
		m_xfp->GetXaxis()->CenterTitle(1);
		//Y_fp
		TH1F *m_yfp = new TH1F("m_yfp","", 100, -0.05,0.05);
		m_yfp->SetXTitle("y_{fp}");
		m_yfp->SetLineColor(4); 
		m_yfp->SetLineWidth(1.1);
		m_yfp->GetXaxis()->SetTitleSize(0.08);
		m_yfp->GetXaxis()->SetTitleOffset(0.7);
		m_yfp->GetXaxis()->CenterTitle(1);
		//Th_fp
		TH1F *m_tfp = new TH1F("m_tfp","", 100, -0.04,0.04);
		m_tfp->SetXTitle("#theta_{fp}");
		m_tfp->SetLineColor(4); 
		m_tfp->SetLineWidth(1.1);
		m_tfp->GetXaxis()->SetTitleSize(0.08);
		m_tfp->GetXaxis()->SetTitleOffset(0.7);
		m_tfp->GetXaxis()->CenterTitle(1);
		//Ph_fp
		TH1F *m_pfp = new TH1F("m_pfp","", 100, -0.06,0.06);
		m_pfp->SetXTitle("#phi_{fp}");
		m_pfp->SetLineColor(4); 
		m_pfp->SetLineWidth(1.1);
		m_pfp->GetXaxis()->SetTitleSize(0.08);
		m_pfp->GetXaxis()->SetTitleOffset(0.7);
		m_pfp->GetXaxis()->CenterTitle(1);
		/*}}}*/

		/*Target Plane{{{*/
		//DeltaP
		TH1F *m_dtg = new TH1F("m_dtg","", 100, -0.07,0.07);
		m_dtg->SetXTitle("#delta P");
		m_dtg->SetLineColor(4); 
		m_dtg->SetLineWidth(1.1);
		m_dtg->GetXaxis()->SetTitleSize(0.08);
		m_dtg->GetXaxis()->SetTitleOffset(0.7);
		m_dtg->GetXaxis()->CenterTitle(1);
		//Y_tg  
		TH1F *m_ytg = new TH1F("m_ytg","", 100, -0.04,0.04);
		m_ytg->SetXTitle("y_{tg}");
		m_ytg->SetLineColor(4); 
		m_ytg->SetLineWidth(1.1);
		m_ytg->GetXaxis()->SetTitleSize(0.08);
		m_ytg->GetXaxis()->SetTitleOffset(0.7);
		m_ytg->GetXaxis()->CenterTitle(1);
		//Z_react
		TH1F *m_ztg = new TH1F("m_ztg","", 100, -0.12,0.12);
		m_ztg->SetXTitle("z_{react}");
		m_ztg->SetLineColor(4); 
		m_ztg->SetLineWidth(1.1);
		m_ztg->GetXaxis()->SetTitleSize(0.08);
		m_ztg->GetXaxis()->SetTitleOffset(0.7);
		m_ztg->GetXaxis()->CenterTitle(1);
		//Th_tg  
		TH1F *m_ttg = new TH1F("m_ttg","", 100, -0.07,0.07);
		m_ttg->SetXTitle("#theta_{tg}");
		m_ttg->SetLineColor(4); 
		m_ttg->SetLineWidth(1.1);
		m_ttg->GetXaxis()->SetTitleSize(0.08);
		m_ttg->GetXaxis()->SetTitleOffset(0.7);
		m_ttg->GetXaxis()->CenterTitle(1);
		//Ph_tg  
		TH1F *m_ptg = new TH1F("m_ptg","", 100, -0.05,0.05);
		m_ptg->SetXTitle("#phi_{tg}");
		m_ptg->SetLineColor(4); 
		m_ptg->SetLineWidth(1.1);
		m_ptg->GetXaxis()->SetTitleSize(0.08);
		m_ptg->GetXaxis()->SetTitleOffset(0.7);
		m_ptg->GetXaxis()->CenterTitle(1);
		//xbj
		TH1F *m_xbj = new TH1F("m_xbj","", 400, Xbj_Min,Xbj_Max);
		m_xbj->SetXTitle("x_{bj}");
		m_xbj->SetLineColor(4); 
		m_xbj->SetLineWidth(1.1);
		m_xbj->GetXaxis()->SetTitleSize(0.08);
		m_xbj->GetXaxis()->SetTitleOffset(0.7);
		m_xbj->GetXaxis()->CenterTitle(1);
		/*}}}*/

		/*Fill Histograms{{{*/
		/*Set BranchAddress{{{*/
		int IsPassed,IsQualified;
		double x_fp,y_fp,th_fp,ph_fp,x_tg_rec, y_tg_rec,th_tg_rec,ph_tg_rec,dp_rec,reactz_rec,Angle_rec;
		T_MC->SetBranchAddress("IsPassed",&IsPassed); //if pass through the magnet, 1=true, 0=false(not pass)
		T_MC->SetBranchAddress("IsQualified",&IsQualified); //if rec var is in the range of gen, 1=true, 0=false(not
		T_MC->SetBranchAddress("x_fp",&x_fp    ); //
		T_MC->SetBranchAddress("y_fp",&y_fp    ); //
		T_MC->SetBranchAddress("th_fp",&th_fp    ); //
		T_MC->SetBranchAddress("ph_fp",&ph_fp    ); //
		T_MC->SetBranchAddress("x_tg_rec",&x_tg_rec    ); //cm
		T_MC->SetBranchAddress("y_tg_rec",&y_tg_rec    ); //cm
		T_MC->SetBranchAddress("th_tg_rec",&th_tg_rec    ); //
		T_MC->SetBranchAddress("ph_tg_rec",&ph_tg_rec    ); //
		T_MC->SetBranchAddress("dp_rec",&dp_rec    );
		T_MC->SetBranchAddress("Angle_rec",&Angle_rec    );
		T_MC->SetBranchAddress("reactz_rec",&reactz_rec    );

		double E_s, E_p, Theta,Xbj;
		T_MC->SetBranchAddress("E_s",&E_s);
		T_MC->SetBranchAddress("E_p",&E_p);
		T_MC->SetBranchAddress("Theta",&Theta);
		T_MC->SetBranchAddress("Xbj",&Xbj);
		/*}}}*/
		Int_t Evt = T_MC->GetEntries();

		/*Loop and Fill{{{*/
		//Define XS table
		Get_XS* xs = new Get_XS();
		xs->Load_Table(Target.Data(),Kin.Data());
		double cs_Center = xs->gGet_XS_Center();
		cerr<<"--> XS Center = "<< cs_Center<<endl;
		if(cs_Center<1e-33){
			cs_Center = 1.0;
			cerr<<"---> Reset cs_Centr = 1"<<endl;
		}
		double cs = 0.0;
		double weight=1.0;

		/*Target Density{{{*/
		TString Density_Profile = Form("/work/halla/e08014/disk1/yez/Monte_Calo/SAMC/Generator/%s.rho",Target.Data());
		TString dummy;
		double VZ[52], Rho[52], Rho1[52],Rho2[52],Rho3[52];
		ifstream density(Density_Profile.Data()); 
		density >> dummy >> dummy >> dummy >> dummy >> dummy;
		double Rho_Sum = 0;
		for(int i=0;i<52;i++){
			density >> VZ[i] >> Rho[i] >> Rho1[i] >> Rho2[i] >> Rho3[i];
			Rho_Sum += Rho3[i];
		}
		Rho_Sum /=52;
		const double VZ_Max = 0.1020, VZ_Min = -0.1020;
		const double step = (VZ_Max-VZ_Min)/51.;//0.004
		/*}}}*/

		int VZ_Bin = 0; double Y_Offset = 0.0;
		int Angle_temp = (int)(aAngle+0.5);
		Double_t Angle_mc = (double)(Angle_temp);
		Double_t y_cut_mc = fabs(0.07 * sin(aAngle*TMath::DegToRad()) );
		//Double_t y_cut_mc = fabs(0.07 * sin(Angle_mc*TMath::DegToRad()) );
		for(int j=0;j<Evt/10;j++){
			T_MC->GetEntry(j);
			y_tg_rec += Y_Offset;//Correct the offset
			if(fabs(x_fp)       < aAccCut[0] &&
					fabs(y_fp)       < aAccCut[1] &&
					fabs(th_fp)      < aAccCut[2] &&
					fabs(ph_fp)      < aAccCut[3] &&
					fabs(dp_rec)     < aAccCut[4] &&
					fabs(y_tg_rec)   < aAccCut[5] &&
					fabs(th_tg_rec)  < aAccCut[6] &&
					fabs(ph_tg_rec)  < aAccCut[7] &&
					fabs(reactz_rec) < aAccCut[8] &&
					IsPassed==1 ){

				double E_p_new = aP0 * (1+dp_rec);
				double Angle_new = acos( (aCosA-ph_tg_rec*aSinA)/sqrt(1.+pow(th_tg_rec,2)+pow(ph_tg_rec,2))); 
				double Q2_new = 4.0*E_s*E_p_new*pow( sin(Angle_new/2.), 2);
				double xbj_new = Q2_new/2./PROTON_MASS/(E_s-E_p_new);

				VZ_Bin = (reactz_rec-VZ_Min)/step;
				cs = xs->gGet_XS(E_s,E_p_new,abs(Angle_new*180./3.14));
				//aNorm = Ntg*Ne*(Omega * DeltaP)/N_Total_Sim
				weight = cs*Rho3[VZ_Bin]/Rho_Sum * aNorm; 

				m_xfp->Fill(x_fp,weight);
				m_yfp->Fill(y_fp,weight);
				m_tfp->Fill(th_fp,weight);
				m_pfp->Fill(ph_fp,weight);

				m_dtg->Fill(dp_rec,weight);
				m_ytg->Fill(y_tg_rec,weight);
				m_ztg->Fill(reactz_rec,weight);
				m_ttg->Fill(th_tg_rec,weight);
				m_ptg->Fill(ph_tg_rec,weight);
				m_xbj->Fill(xbj_new,weight);
			}

		}
		delete xs;
		/*}}}*/

		/*}}}*/

		/*}}}*/

		/*Plot{{{*/
		TH1F *r_xfp = (TH1F*) h_xfp->Clone("r_xfp"); r_xfp->Divide(m_xfp);
		TH1F *r_yfp = (TH1F*) h_yfp->Clone("r_yfp"); r_yfp->Divide(m_yfp);
		TH1F *r_tfp = (TH1F*) h_tfp->Clone("r_tfp"); r_tfp->Divide(m_tfp);
		TH1F *r_pfp = (TH1F*) h_pfp->Clone("r_pfp"); r_pfp->Divide(m_pfp);

		r_xfp->SetTitle(Form("x_{fp} Ratio (EX/MC) %s@%s (#%d)",     Target.Data(),Kin.Data(),run_number)); 
		r_xfp->SetTitleSize(2.3); r_xfp->SetLineColor(1);
		r_yfp->SetTitle(Form("y_{fp} Ratio (EX/MC) %s@%s (#%d)",     Target.Data(),Kin.Data(),run_number)); 
		r_yfp->SetTitleSize(2.3); r_yfp->SetLineColor(1);
		r_tfp->SetTitle(Form("#theta_{fp} Ratio (EX/MC) %s@%s (#%d)",Target.Data(),Kin.Data(),run_number)); 
		r_tfp->SetTitleSize(2.3); r_tfp->SetLineColor(1);
		r_pfp->SetTitle(Form("#phi_{fp} Ratio (EX/MC) %s@%s (#%d)",  Target.Data(),Kin.Data(),run_number)); 
		r_pfp->SetTitleSize(2.3); r_pfp->SetLineColor(1);    

		TH1F *r_dtg = (TH1F*) h_dtg->Clone("r_dtg"); r_dtg->Divide(m_dtg);
		TH1F *r_ytg = (TH1F*) h_ytg->Clone("r_ytg"); r_ytg->Divide(m_ytg);
		TH1F *r_ztg = (TH1F*) h_ztg->Clone("r_ztg"); r_ztg->Divide(m_ztg);
		TH1F *r_ttg = (TH1F*) h_ttg->Clone("r_ttg"); r_ttg->Divide(m_ttg);
		TH1F *r_ptg = (TH1F*) h_ptg->Clone("r_ptg"); r_ptg->Divide(m_ptg);

		r_dtg->SetTitle(Form("#deltaP Ratio (EX/MC) %s@%s (#%d)",    Target.Data(),Kin.Data(),run_number)); 
		r_dtg->SetTitleSize(2.3); r_dtg->SetLineColor(1);
		r_ztg->SetTitle(Form("ReactZ Ratio (EX/MC) %s@%s (#%d)",     Target.Data(),Kin.Data(),run_number));
		r_ztg->SetTitleSize(2.3); r_ztg->SetLineColor(1);
		r_ytg->SetTitle(Form("y_{tg} Ratio (EX/MC) %s@%s (#%d)",     Target.Data(),Kin.Data(),run_number));
		r_ytg->SetTitleSize(2.3); r_ytg->SetLineColor(1);
		r_ttg->SetTitle(Form("#theta_{tg} Ratio (EX/MC) %s@%s (#%d)",Target.Data(),Kin.Data(),run_number));
		r_ttg->SetTitleSize(2.3); r_ttg->SetLineColor(1);
		r_ptg->SetTitle(Form("#phi_{tg} Ratio (EX/MC) %s@%s (#%d)",  Target.Data(),Kin.Data(),run_number)); 
		r_ptg->SetTitleSize(2.3); r_ptg->SetLineColor(1);    

		TCanvas *c1 = new TCanvas("c1","c1 Focal Plane",800,1000);
		c1->Divide(2,4);

		gPad->SetLogy(0);
		c1->cd(1); h_xfp->Draw();
		c1->cd(3); h_yfp->Draw();
		c1->cd(5); h_tfp->Draw();
		c1->cd(7); h_pfp->Draw();
		c1->cd(1); m_xfp->Draw("same");  
		c1->cd(3); m_yfp->Draw("same"); 
		c1->cd(5); m_tfp->Draw("same"); 
		c1->cd(7); m_pfp->Draw("same"); 

		c1->cd(2);gPad->SetLogy(1);gPad->SetGridy(1);r_xfp->GetYaxis()->SetRangeUser(0.1,10.);r_xfp->Draw("");  
		c1->cd(4);gPad->SetLogy(1);gPad->SetGridy(1);r_yfp->GetYaxis()->SetRangeUser(0.1,10.);r_yfp->Draw(""); 
		c1->cd(6);gPad->SetLogy(1);gPad->SetGridy(1);r_tfp->GetYaxis()->SetRangeUser(0.1,10.);r_tfp->Draw(""); 
		c1->cd(8);gPad->SetLogy(1);gPad->SetGridy(1);r_pfp->GetYaxis()->SetRangeUser(0.1,10.);r_pfp->Draw(""); 


		TCanvas *c2 = new TCanvas("c2","c2 Target Plane",800,1000);
		c2->Divide(2,4);

		gPad->SetLogy(0);
		c2->cd(1); h_dtg->Draw();
		//c2->cd(3); h_ztg->Draw();
		c2->cd(3); h_ytg->Draw();
		c2->cd(5); h_ttg->Draw();
		c2->cd(7); h_ptg->Draw();

		c2->cd(1); m_dtg->Draw("same");  
		//c2->cd(3); m_ztg->Draw("same"); 
		c2->cd(3); m_ytg->Draw("same"); 
		c2->cd(5); m_ttg->Draw("same"); 
		c2->cd(7); m_ptg->Draw("same"); 

		c2->cd(2);gPad->SetLogy(1);gPad->SetGridy(1);r_dtg->GetYaxis()->SetRangeUser(0.1,10.); r_dtg->Draw("");  
		//c2->cd(4);gPad->SetLogy(1);gPad->SetGridy(1);r_ztg->GetYaxis()->SetRangeUser(0.1,10.); r_ztg->Draw("same"); 
		c2->cd(4);gPad->SetLogy(1);gPad->SetGridy(1);r_ytg->GetYaxis()->SetRangeUser(0.1,10.);r_ytg->Draw(""); 
		c2->cd(6);gPad->SetLogy(1);gPad->SetGridy(1);r_ttg->GetYaxis()->SetRangeUser(0.1,10.); r_ttg->Draw(""); 
		c2->cd(8);gPad->SetLogy(1);gPad->SetGridy(1);r_ptg->GetYaxis()->SetRangeUser(0.1,10.); r_ptg->Draw(""); 


		SAMC_File.ReplaceAll(".root","_");
		TString pic=Target+"_"+Kin+"_"+Form("%d_fp_NoCut.png",run_number);
		c1->Print(pic);
		pic=Target+"_"+Kin+"_"+Form("%d_fp_NoCut.pdf",run_number);
		c1->Print(pic);
		pic=Target+"_"+Kin+"_"+Form("%d_tg_NoCut.png",run_number);
		c2->Print(pic);
		pic=Target+"_"+Kin+"_"+Form("%d_tg_NoCut.pdf",run_number);
		c2->Print(pic);
		//TCanvas *c2 = new TCanvas("c2","c2",1200,700);
		//c2->Divide(2,2);

		TH1F *r_xbj = (TH1F*) h_xbj->Clone("r_xbj"); r_xbj->Divide(m_xbj);
		r_xbj->SetTitle(Form("x_{bj} Ratio (EX/MC) %s@%s (#%d)",  Target.Data(),Kin.Data(),run_number)); 
		r_xbj->SetXTitle("x_{bj}"); r_xbj->SetTitleSize(2.3); r_xbj->GetYaxis()->SetRangeUser(0.1,10.);r_xbj->SetLineColor(1);    

		TCanvas *c3 = new TCanvas("c3","c3 Focal Plane",900,900);
		c3->Divide(1,2);
		c3->cd(1); gPad->SetLogy(1); h_xbj->Draw(); m_xbj->Draw("same");
		c3->cd(2); gPad->SetLogy(1); gPad->SetGridy(1); r_xbj->Draw();

		pic=Target+"_"+Kin+"_"+Form("%d_xbj_NoCut.png",run_number);
		c3->Print(pic);
		pic=Target+"_"+Kin+"_"+Form("%d_xbj_NoCut.pdf",run_number);
		c3->Print(pic);

		/*}}}*/
	}
	/*}}}*/




}
/*}}}*/
