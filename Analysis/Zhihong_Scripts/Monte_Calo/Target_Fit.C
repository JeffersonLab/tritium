{
	/*Define{{{*/
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	TString hrs="L";
	TString GeneralCut = "L.tr.n==1";
	TString TriggerCut3 = "((DBB.evtypebits>>3)&1)";
	TString TriggerCut7 = "((DBB.evtypebits>>7)&1)";
	TString Acc_real = "abs(ExTgtL.th)<0.03&&abs(ExTgtL.ph)<0.02&&abs(ExTgtL.y)<0.08&&abs(ExTgtL.dp)<0.03&&abs(RctPtL.z)<0.075";
	TString Cut_All = GeneralCut+"&&"+TriggerCut3+"&&"+Acc_real
		+"&&L.cer.asum_c>=50.&&L.prl2.e>=100.&&(L.prl1.e+L.prl2.e)/L.tr.p/1000.>=0.5";
	double scale;  int norm = 100000;
	double vz=0, rho=0, rho_sum=0;
    /*}}}*/

	TCanvas *c2 = new TCanvas("c2","c2",1200,1000);
    c2->Divide(2,3);

	/*H2{{{*/
	c2->cd(1);
	TString kin="kin3.1_unf"; int runno_h2 = 3681;
	TString filename1=Form("/work/halla/e08014/disk1/Rootfiles/e08014_%d.root",runno_h2);
	TFile *f1 = new TFile(filename1);
	TFile *f11 = new TFile("h2_"+hrs+"_"+kin+".root");

	f1->cd();
	T->SetLineColor(2);
	T->Draw("RctPtL.z>>j1(200,-0.12,0.12)",Cut_All);
	j1->SetLineColor(2); scale = norm/j1->Integral(); j1->Scale(scale); 
	j1->SetXTitle("z_{react}");
	j1->GetXaxis()->CenterTitle(1);
	j1->GetXaxis()->SetTitleFont(32);
	j1->GetXaxis()->SetTitleSize(0.08);
	j1->GetXaxis()->SetTitleOffset(0.7);
	j1->SetTitle("");
	j1->Draw();

	TH1F *j1_mc = (TH1F*) f11->Get("hrpz");
	j1_mc->SetLineColor(4);
	scale = norm/j1_mc->Integral();  j1_mc->Scale(scale);  
	j1_mc ->Draw("same");

	TLatex *tex1 = new TLatex(0.3,0.3,"H2 Density Profile: Blue->MC, Red->EX");
    tex1->SetNDC();
	tex1->SetTextFont(32);
	tex1->Draw();

	c2->cd(2);
    TH1F *r1 = (TH1F*) j1->Clone("r1");
	r1->Divide(j1_mc) ;
	r1->SetLineColor(4);
	r1->Draw();
    r1->Fit("pol5","","",-0.075,0.075);

	TLatex *tex11 = new TLatex(0.4,0.3,"H2 Density Distribution");
    tex11->SetNDC();
	tex11->SetTextFont(32);
	tex11->Draw();
	
	ofstream out_h2("h2_"+hrs+"_"+kin+".rho");
	rho_sum=0.;
	for(int i=0;i<200;i++)
		rho_sum+=r1->GetBinContent(i);
	cerr<<" --- H2: Rho_Sum = "<<rho_sum<<endl;
	for(int i=0;i<200;i++){
		vz =r1->GetBinCenter(i);
		rho=r1->GetBinContent(i);
        out_h2 << " "<< vz<<" "<<rho<<endl;
	}

/*}}}*/

	/*He3{{{*/
	c2->cd(3);
	TString kin="kin3.1_unf"; int runno_he3 = 4085;
	TString filename2=Form("/work/halla/e08014/disk1/Rootfiles/e08014_%d.root",runno_he3);
	TFile *f2 = new TFile(filename2);
	TFile *f22 = new TFile("he3_"+hrs+"_"+kin+".root");

	f2->cd();
	T->SetLineColor(2);
	T->Draw("RctPtL.z>>j2(200,-0.12,0.12)",Cut_All);
	j2->SetLineColor(2); scale = norm/j2->Integral(); j2->Scale(scale); 
	j2->SetXTitle("z_{react}");
	j2->GetXaxis()->CenterTitle(1);
	j2->GetXaxis()->SetTitleFont(32);
	j2->GetXaxis()->SetTitleSize(0.08);
	j2->GetXaxis()->SetTitleOffset(0.7);
	j2->SetTitle("");
	j2->Draw();

	TH1F *j2_mc = (TH1F*) f22->Get("hrpz");
	j2_mc->SetLineColor(4);
	scale = norm/j2_mc->Integral();  j2_mc->Scale(scale);  
	j2_mc ->Draw("same");

	TLatex *tex2 = new TLatex(0.3,0.3,"He3 Density Profile: Blue->MC, Red->EX");
    tex2->SetNDC();
	tex2->SetTextFont(32);
	tex2->Draw();

	c2->cd(4);
    TH1F *r2 = (TH1F*) j2->Clone("r2");
	r2->Divide(j2_mc) ;
	r2->SetLineColor(4);
	r2->Draw();
    r2->Fit("pol5","","",-0.075,0.075);

	TLatex *tex22 = new TLatex(0.4,0.3,"He3 Density Distribution");
    tex22->SetNDC();
	tex22->SetTextFont(32);
	tex22->Draw();
	
	ofstream out_he3("he3_"+hrs+"_"+kin+".rho");
	rho_sum=0.;
	for(int i=0;i<200;i++)
		rho_sum+=r2->GetBinContent(i);
	cerr<<" --- He3: Rho_Sum = "<<rho_sum<<endl;
	for(int i=0;i<200;i++){
		vz =r2->GetBinCenter(i);
		rho=r2->GetBinContent(i);
        out_he3 << " "<< vz<<" "<<rho<<endl;
	}

	/*}}}*/
	
	/*He4{{{*/
	c2->cd(5);
	TString kin="kin3.1_unf"; int runno_he4 = 3676;
	TString filename3=Form("/work/halla/e08014/disk1/Rootfiles/e08014_%d.root",runno_he4);
	TFile *f3 = new TFile(filename3);
	TFile *f33 = new TFile("he4_"+hrs+"_"+kin+".root");

	f3->cd();
	T->SetLineColor(2);
	T->Draw("RctPtL.z>>j3(200,-0.12,0.12)",Cut_All);
	j3->SetLineColor(2); scale = norm/j3->Integral(); j3->Scale(scale); 
	j3->SetXTitle("z_{react}");
	j3->GetXaxis()->CenterTitle(1);
	j3->GetXaxis()->SetTitleFont(32);
	j3->GetXaxis()->SetTitleSize(0.08);
	j3->GetXaxis()->SetTitleOffset(0.7);
	j3->SetTitle("");
	j3->Draw();

	TH1F *j3_mc = (TH1F*) f33->Get("hrpz");
	j3_mc->SetLineColor(4);
	scale = norm/j3_mc->Integral();  j3_mc->Scale(scale);  
	j3_mc ->Draw("same");
	
	TLatex *tex3 = new TLatex(0.3,0.3,"He4 Density Profile: Blue->MC, Red->EX");
    tex3->SetNDC();
	tex3->SetTextFont(32);
	tex3->Draw();

	c2->cd(6);
    TH1F *r3 = (TH1F*) j3->Clone("r3");
	r3->Divide(j3_mc) ;
	r3->SetLineColor(4);
	r3->Draw();
    r3->Fit("pol5","","",-0.075,0.075);
	
	TLatex *tex33 = new TLatex(0.4,0.3,"He4 Density Distribution");
    tex33->SetNDC();
	tex33->SetTextFont(32);
	tex33->Draw();
	
	ofstream out_he4("he4_"+hrs+"_"+kin+".rho");
	rho_sum=0.;
	for(int i=0;i<200;i++)
		rho_sum+=r3->GetBinContent(i);
	cerr<<" --- He4: Rho_Sum = "<<rho_sum<<endl;
	for(int i=0;i<200;i++){
		vz =r3->GetBinCenter(i);
		rho=r3->GetBinContent(i);
        out_he4 << " "<< vz<<" "<<rho<<endl;
	}

	/*}}}*/
}
