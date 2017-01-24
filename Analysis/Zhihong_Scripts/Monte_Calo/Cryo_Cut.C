/*
  Weight some variables by cross section.
  -- Zhihong Ye, 08/08/2011
*/
#include "/work/halla/e08014/disk1/yez/Cross_Section/SRC/XGT2_Main.h"
#include "/work/halla/e08014/disk1/yez/Cross_Section/SRC/XGT2_XEMC.h"

void plot(){
	gStyle->SetOptStat(0);

	int run_number = 4085; cerr<<"---- Run Number = "; cin >> run_number;
	TString SAMC_File = "";cerr<<"---- SAMC File = "; cin >> SAMC_File;
    TString Target = ""; cerr<<"---- Target = "; cin >> Target;
	TString Kin = ""; cerr<<"---- Kin = "; cin >> Kin;

	/*Ploting EX Data{{{*/
	TString GeneralCut = "L.tr.n==1";//&&abs(L.tr.x)<0.75 && abs(L.tr.y)<0.55 && abs(L.tr.th)<0.15 && abs(L.tr.ph)<0.045";
	TString TriggerCut3 = "((DBB.evtypebits>>3)&1)";
	TString TriggerCut7 = "((DBB.evtypebits>>7)&1)";
	Double_t y_cut = 0.170;
	TString Acc_real = Form("abs(ExTgtL.th)<0.03&&abs(ExTgtL.ph)<0.02&&abs(ExTgtL.y)<%f&&abs(ExTgtL.dp)<0.04&&abs(RctPtL.z)<0.070",y_cut);
	TString Cut_All = GeneralCut+"&&"+TriggerCut3+"&&"+Acc_real
		+"&&L.cer.asum_c>=50.&&L.prl2.e>=100.&&(L.prl1.e+L.prl2.e)/L.tr.p/1000.>=0.5";
	//TString Cut_All ="1";

	double scale;  int norm = 100000;

	TString filename=Form("/work/halla/e08014/disk1/Rootfiles/e08014_%d.root",run_number);
	TChain *T = new TChain("T");
	T->Add(filename);

	/*Focal Plane{{{*/
	//X_fp
	TH1F *h_xfp = new TH1F("h_xfp","", 100, -1.15,1.15);
	h_xfp->SetXTitle("x_{fp}");
	h_xfp->SetLineColor(2); 
	h_xfp->SetLineWidth(1.1);
	h_xfp->GetXaxis()->SetTitleSize(0.08);
	h_xfp->GetXaxis()->SetTitleOffset(0.7);
	h_xfp->GetXaxis()->CenterTitle(1);
	T->Project("h_xfp","L.tr.r_x",Cut_All);
	scale = norm/h_xfp->Integral(); h_xfp->Scale(scale);
	//Y_fp
	TH1F *h_yfp = new TH1F("h_yfp","", 100, -0.05,0.05);
	h_yfp->SetXTitle("y_{fp}");
	h_yfp->SetLineColor(2); 
	h_yfp->SetLineWidth(1.1);
	h_yfp->GetXaxis()->SetTitleSize(0.08);
	h_yfp->GetXaxis()->SetTitleOffset(0.7);
	h_yfp->GetXaxis()->CenterTitle(1);
	T->Project("h_yfp","L.tr.r_y",Cut_All);
	scale = norm/h_yfp->Integral(); h_yfp->Scale(scale);
	//Th_fp
	TH1F *h_tfp = new TH1F("h_tfp","", 100, -0.02,0.02);
	h_tfp->SetXTitle("#theta_{fp}");
	h_tfp->SetLineColor(2); 
	h_tfp->SetLineWidth(1.1);
	h_tfp->GetXaxis()->SetTitleSize(0.08);
	h_tfp->GetXaxis()->SetTitleOffset(0.7);
	h_tfp->GetXaxis()->CenterTitle(1);
	T->Project("h_tfp","L.tr.r_th",Cut_All);
	scale = norm/h_tfp->Integral(); h_tfp->Scale(scale);
	//Ph_fp
	TH1F *h_pfp = new TH1F("h_pfp","", 100, -0.05,0.05);
	h_pfp->SetXTitle("#phi_{fp}");
	h_pfp->SetLineColor(2); 
	h_pfp->SetLineWidth(1.1);
	h_pfp->GetXaxis()->SetTitleSize(0.08);
	h_pfp->GetXaxis()->SetTitleOffset(0.7);
	h_pfp->GetXaxis()->CenterTitle(1);
	T->Project("h_pfp","L.tr.r_ph",Cut_All);
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
	T->Project("h_dtg","ExTgtL.dp",Cut_All);
	scale = norm/h_dtg->Integral(); h_dtg->Scale(scale);
	//Y_tg  
	TH1F *h_ytg = new TH1F("h_ytg","", 100, -0.07,0.07);
	h_ytg->SetXTitle("y_{tg}");
	h_ytg->SetLineColor(2); 
	h_ytg->SetLineWidth(1.1);
	h_ytg->GetXaxis()->SetTitleSize(0.08);
	h_ytg->GetXaxis()->SetTitleOffset(0.7);
	h_ytg->GetXaxis()->CenterTitle(1);
	T->Project("h_ytg","ExTgtL.y",Cut_All);
	scale = norm/h_ytg->Integral(); h_ytg->Scale(scale);
	//Z_react
	TH1F *h_ztg = new TH1F("h_ztg","", 100, -0.12,0.12);
	h_ztg->SetXTitle("z_{react}");
	h_ztg->SetLineColor(2); 
	h_ztg->SetLineWidth(1.1);
	h_ztg->GetXaxis()->SetTitleSize(0.08);
	h_ztg->GetXaxis()->SetTitleOffset(0.7);
	h_ztg->GetXaxis()->CenterTitle(1);
	T->Project("h_ztg","RctPtL.z",Cut_All);
	scale = norm/h_ztg->Integral(); h_ztg->Scale(scale);
	//Th_tg  
	TH1F *h_ttg = new TH1F("h_ttg","", 100, -0.05,0.05);
	h_ttg->SetXTitle("#theta_{tg}");
	h_ttg->SetLineColor(2); 
	h_ttg->SetLineWidth(1.1);
	h_ttg->GetXaxis()->SetTitleSize(0.08);
	h_ttg->GetXaxis()->SetTitleOffset(0.7);
	h_ttg->GetXaxis()->CenterTitle(1);
	T->Project("h_ttg","ExTgtL.th",Cut_All);
	scale = norm/h_ttg->Integral(); h_ttg->Scale(scale);
	//Ph_tg  
	TH1F *h_ptg = new TH1F("h_ptg","", 100, -0.03,0.03);
	h_ptg->SetXTitle("#phi_{tg}");
	h_ptg->SetLineColor(2); 
	h_ptg->SetLineWidth(1.1);
	h_ptg->GetXaxis()->SetTitleSize(0.08);
	h_ptg->GetXaxis()->SetTitleOffset(0.7);
	h_ptg->GetXaxis()->CenterTitle(1);
	T->Project("h_ptg","ExTgtL.ph",Cut_All);
	scale = norm/h_ptg->Integral(); h_ptg->Scale(scale);
	/*}}}*/

   	/*}}}*/

	/*Ploting MC Data{{{*/
	 TChain* S = new TChain("SAMC");
	 S->Add(SAMC_File);

	TString Cut_MC = Form("abs(th_tg_rec)<0.03&&abs(ph_tg_rec)<0.02&&abs(y_tg_rec)<%f&&abs(dp_rec)<0.04&&abs(reactz_rec)<0.10",y_cut);

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
	TH1F *m_tfp = new TH1F("m_tfp","", 100, -0.02,0.02);
	m_tfp->SetXTitle("#theta_{fp}");
	m_tfp->SetLineColor(4); 
	m_tfp->SetLineWidth(1.1);
	m_tfp->GetXaxis()->SetTitleSize(0.08);
	m_tfp->GetXaxis()->SetTitleOffset(0.7);
	m_tfp->GetXaxis()->CenterTitle(1);
	//Ph_fp
	TH1F *m_pfp = new TH1F("m_pfp","", 100, -0.05,0.05);
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
	TH1F *m_ytg = new TH1F("m_ytg","", 100, -0.07,0.07);
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
	TH1F *m_ttg = new TH1F("m_ttg","", 100, -0.05,0.05);
	m_ttg->SetXTitle("#theta_{tg}");
	m_ttg->SetLineColor(4); 
	m_ttg->SetLineWidth(1.1);
	m_ttg->GetXaxis()->SetTitleSize(0.08);
	m_ttg->GetXaxis()->SetTitleOffset(0.7);
	m_ttg->GetXaxis()->CenterTitle(1);
	//Ph_tg  
	TH1F *m_ptg = new TH1F("m_ptg","", 100, -0.03,0.03);
	m_ptg->SetXTitle("#phi_{tg}");
	m_ptg->SetLineColor(4); 
	m_ptg->SetLineWidth(1.1);
	m_ptg->GetXaxis()->SetTitleSize(0.08);
	m_ptg->GetXaxis()->SetTitleOffset(0.7);
	m_ptg->GetXaxis()->CenterTitle(1);
	/*}}}*/

	/*Fill Histograms{{{*/
	/*Set BranchAddress{{{*/
    int IsPassed,IsQualified;
	double x_fp,y_fp,th_fp,ph_fp,x_tg_rec, y_tg_rec,th_tg_rec,ph_tg_rec,dp_rec,reactz_rec,Angle_rec;
	S->SetBranchAddress("IsPassed",&IsPassed); //if pass through the magnet, 1=true, 0=false(not pass)
	S->SetBranchAddress("IsQualified",&IsQualified); //if rec var is in the range of gen, 1=true, 0=false(not
	S->SetBranchAddress("x_fp",&x_fp    ); //
	S->SetBranchAddress("y_fp",&y_fp    ); //
	S->SetBranchAddress("th_fp",&th_fp    ); //
	S->SetBranchAddress("ph_fp",&ph_fp    ); //
	S->SetBranchAddress("x_tg_rec",&x_tg_rec    ); //cm
	S->SetBranchAddress("y_tg_rec",&y_tg_rec    ); //cm
	S->SetBranchAddress("th_tg_rec",&th_tg_rec    ); //
	S->SetBranchAddress("ph_tg_rec",&ph_tg_rec    ); //
	S->SetBranchAddress("dp_rec",&dp_rec    );
	S->SetBranchAddress("Angle_rec",&Angle_rec    );
	S->SetBranchAddress("reactz_rec",&reactz_rec    );

	double E_s, E_p, Theta,Xbj;
	S->SetBranchAddress("E_s",&E_s);
	S->SetBranchAddress("E_p",&E_p);
	S->SetBranchAddress("Theta",&Theta);
	S->SetBranchAddress("Xbj",&Xbj);
	/*}}}*/
	Int_t Evt = S->GetEntries();

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
	for(int j=0;j<Evt;j++){
		S->GetEntry(j);
		y_tg_rec += Y_Offset;//Correct the offset
		//    if(IsPassed){
		if(IsPassed
				&&abs(reactz_rec)<0.07
		    	&&abs(y_tg_rec+0.00)<(0.175)
				&&abs(th_tg_rec)<0.03
				&&abs(ph_tg_rec)<0.02
				&&abs(dp_rec)<0.04
				){	
			VZ_Bin = (reactz_rec-VZ_Min)/step;
			cs = xs->gGet_XS(E_s,E_p,abs(Angle_rec*180./3.14));
			weight = cs/cs_Center*Rho3[VZ_Bin]/Rho_Sum; 
   
			m_xfp->Fill(x_fp,weight);
			m_yfp->Fill(y_fp,weight);
			m_tfp->Fill(th_fp,weight);
			m_pfp->Fill(ph_fp,weight);

			m_dtg->Fill(dp_rec,weight);
			m_ytg->Fill(y_tg_rec,weight);
			m_ztg->Fill(reactz_rec,weight);
			m_ttg->Fill(th_tg_rec,weight);
			m_ptg->Fill(ph_tg_rec,weight);
		}

	}
	delete xs;
	/*}}}*/

	scale = norm/m_xfp->Integral(); m_xfp->Scale(scale);
	scale = norm/m_yfp->Integral(); m_yfp->Scale(scale);
	scale = norm/m_tfp->Integral(); m_tfp->Scale(scale);
	scale = norm/m_pfp->Integral(); m_pfp->Scale(scale);

	scale = norm/m_dtg->Integral(); m_dtg->Scale(scale);
	scale = norm/m_ytg->Integral(); m_ytg->Scale(scale);
	scale = norm/m_ztg->Integral(); m_ztg->Scale(scale);
	scale = norm/m_ttg->Integral(); m_ttg->Scale(scale);
	scale = norm/m_ptg->Integral(); m_ptg->Scale(scale);
	/*}}}*/

	/*}}}*/

	/*Plot{{{*/
	TH1F *r_xfp = (TH1F*) h_xfp->Clone("r_xfp"); r_xfp->Divide(m_xfp);
	TH1F *r_yfp = (TH1F*) h_yfp->Clone("r_yfp"); r_yfp->Divide(m_yfp);
	TH1F *r_tfp = (TH1F*) h_tfp->Clone("r_tfp"); r_tfp->Divide(m_tfp);
	TH1F *r_pfp = (TH1F*) h_pfp->Clone("r_pfp"); r_pfp->Divide(m_pfp);
    
	r_xfp->SetTitle("#x_{fp} Ratio (EX/MC)");     r_xfp->SetLineColor(1);
	r_yfp->SetTitle("y_{fp} Ratio (EX/MC)");      r_yfp->SetLineColor(1);
	r_tfp->SetTitle("#theta_{fp} Ratio (EX/MC)"); r_tfp->SetLineColor(1);
	r_pfp->SetTitle("#phi_{fp} Ratio (EX/MC)");   r_pfp->SetLineColor(1);


	TH1F *r_dtg = (TH1F*) h_dtg->Clone("r_dtg"); r_dtg->Divide(m_dtg);
	TH1F *r_ytg = (TH1F*) h_ytg->Clone("r_ytg"); r_ytg->Divide(m_ytg);
	TH1F *r_ztg = (TH1F*) h_ztg->Clone("r_ztg"); r_ztg->Divide(m_ztg);
	TH1F *r_ttg = (TH1F*) h_ttg->Clone("r_ttg"); r_ttg->Divide(m_ttg);
	TH1F *r_ptg = (TH1F*) h_ptg->Clone("r_ptg"); r_ptg->Divide(m_ptg);
    
	r_dtg->SetTitle("#deltaP Ratio (EX/MC)");     r_dtg->SetLineColor(1);
	r_ztg->SetTitle("ReactZ Ratio (EX/MC)");      r_ztg->SetLineColor(1);
	r_ytg->SetTitle("y_{tg} Ratio (EX/MC)");      r_ytg->SetLineColor(1);
	r_ttg->SetTitle("#theta_{tg} Ratio (EX/MC)"); r_ttg->SetLineColor(1);
	r_ptg->SetTitle("#phi_{tg} Ratio (EX/MC)");   r_ptg->SetLineColor(1);

	TCanvas *c1 = new TCanvas("c1","c1 Focal Plane",800,1000);
	c1->Divide(2,4);

	c1->cd(1); h_xfp->Draw();
	c1->cd(3); h_yfp->Draw();
	c1->cd(5); h_tfp->Draw();
	c1->cd(7); h_pfp->Draw();
	c1->cd(1); m_xfp->Draw("same");  
	c1->cd(3); m_yfp->Draw("same"); 
	c1->cd(5); m_tfp->Draw("same"); 
	c1->cd(7); m_pfp->Draw("same"); 

	c1->cd(2); r_xfp->Draw("");  
	c1->cd(4); r_yfp->Draw(""); 
	c1->cd(6); r_tfp->Draw(""); 
	c1->cd(8); r_pfp->Draw(""); 


	TCanvas *c2 = new TCanvas("c2","c2 Target Plane",800,1000);
	c2->Divide(2,4);

	c2->cd(1); h_dtg->Draw();
	c2->cd(3); h_ztg->Draw();
	//c2->cd(5); h_ytg->Draw();
	c2->cd(5); h_ttg->Draw();
	c2->cd(7); h_ptg->Draw();

	c2->cd(1); m_dtg->Draw("same");  
	c2->cd(3); m_ztg->Draw("same"); 
	//c2->cd(3); m_ytg->Draw("same"); 
	c2->cd(5); m_ttg->Draw("same"); 
	c2->cd(7); m_ptg->Draw("same"); 

	c2->cd(2); r_dtg->Draw("");  
	c2->cd(4); r_ztg->Draw("same"); 
	//c2->cd(4); r_ytg->Draw(""); 
	c2->cd(6); r_ttg->Draw(""); 
	c2->cd(8); r_ptg->Draw(""); 


	SAMC_File.ReplaceAll(".root","_");
	TString pic=Target+"_"+SAMC_File+Form("%d_fp_cut.png",run_number);
	c1->Print(pic);
	pic=Target+"_"+SAMC_File+Form("%d_tg_cut.png",run_number);
	c2->Print(pic);
	//TCanvas *c2 = new TCanvas("c2","c2",1200,700);
	//c2->Divide(2,2);
	/*}}}*/
}
