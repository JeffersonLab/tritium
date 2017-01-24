void plot(){

	TString Target_Name; cerr<<"---Target (H2,He3,He4) = "; cin >> Target_Name;
	ifstream infileL(Form("%s_Yield_L_All.out",Target_Name.Data()));
	ifstream infileR(Form("%s_Yield_R_All.out",Target_Name.Data()));
 
	const int bin=20;
	double Zero[bin], VZL[bin], VZR[bin],Y0R[bin],Y0L[bin], BFR[bin],BFL[bin];
	double Y0R_Err[bin],Y0L_Err[bin],BFR_Err[bin],BFL_Err[bin];
	double offset = 0.00;
	if(Target_Name=="H2")
		offset = 0.012;
	else if(Target_Name=="He3")
		offset = 0.016;
	else if(Target_Name=="He4")
		offset = 0.016;
	else
		cerr<<"*** ERROR ***"<<endl;
	TString com;
	infileL >> com  >> com >> com >> com>> com;
	infileR >> com  >> com >> com >> com>> com;
	double Y_Max = -19988, Y_Min = 19999;
	double Y_SumL=0.0, Y_SumR=0.0;
	for(int i=0;i<bin;i++){
	   infileL >> VZL[i] >> Y0L[i] >> Y0L_Err[i] >> BFL[i] >> BFL_Err[i];
	   infileR >> VZR[i] >> Y0R[i] >> Y0R_Err[i] >> BFR[i] >> BFR_Err[i];
	   VZR[i]-=offset;
	   Zero[i]=0.0;

	   if(abs(VZL[i])<0.075)
		   Y_SumL+=Y0L[i];
	   else
		   Y0L[i]=0.0;
	   if(abs(VZR[i])<0.075)
		   Y_SumR+=Y0R[i];
	   else
		   Y0R[i]=0.0;
	}

	cerr<<Form("---Y_Sum_L=%f, Y_Sum_R=%f", Y_SumL, Y_SumR)<<endl;
	for(int i=0;i<bin;i++){
       Y0L[i]/=Y_SumL;
       Y0R[i]/=Y_SumR;

	   Y0L_Err[i]/=Y_SumL;
	   Y0R_Err[i]/=Y_SumR;

	   if(Y_Max<Y0L[i])
		   Y_Max = Y0L[i];
	   if(Y_Min>Y0L[i])
		   Y_Min = Y0L[i];
	   if(Y_Max<Y0R[i])
		   Y_Max = Y0R[i];
	   if(Y_Min>Y0R[i])
		   Y_Min = Y0R[i];
	}
	
	/*Plot BF{{{*/
	gROOT->SetStyle("Plain");
	gStyle->SetFillColor(0);
	cerr<<" @@@ Making Plots BF ..."<<endl;
	TCanvas *c1 = new TCanvas("c1","c1",1200,700);
	c1->cd();
	TH2F *h1 = new TH2F("h1","",300, -0.14, 0.14, 300, -0.015, 0.005);
	h1->SetStats(kFALSE);
	h1->SetXTitle("z_{react} (m)");
	h1->GetXaxis()->CenterTitle(1);
	h1->GetXaxis()->SetTitleFont(32);
	h1->GetXaxis()->SetTitleSize(0.08);
	h1->GetXaxis()->SetTitleOffset(0.7);
	h1->SetYTitle("BF");
	h1->GetYaxis()->CenterTitle(1);
	h1->GetYaxis()->SetTitleFont(32);
	h1->GetYaxis()->SetTitleSize(0.07);
	h1->GetYaxis()->SetTitleOffset(0.8);
	h1->Draw();

	TGraphErrors *exL = new TGraphErrors(bin, VZL, BFL, Zero, BFL_Err);
	exL->SetMarkerStyle(20);
	exL->SetMarkerColor(4);
	exL->SetMarkerSize(1.2);
	exL->SetLineColor(4);
	exL->Draw("P");
	TGraphErrors *exR = new TGraphErrors(bin, VZR, BFR, Zero, BFL_Err);
	exR->SetMarkerStyle(25);
	exR->SetMarkerSize(1.2);
	exR->SetMarkerColor(4);
	exR->SetLineColor(2);
	exR->Draw("P");

	TLegend *l1 = new TLegend(0.35,0.3,0.7,0.45,Form("%s Boiling Study: #rho(I)=#rho_{0}(1+BF#times I)",Target_Name.Data()));
	l1->SetTextFont(32);
	l1->SetTextSize(0.03);
	l1->SetShadowColor(0);
	l1->AddEntry(exL,"HRS-L","p");
	l1->AddEntry(exR,"HRS-R","p");
	l1->Draw();	

	/*}}}*/

	/*Plot Y0{{{*/
	gROOT->SetStyle("Plain");
	gStyle->SetFillColor(0);
	cerr<<" @@@ Making Plots Y0 ..."<<endl;
	TCanvas *c2 = new TCanvas("c2","c2",1200,700);
	c2->cd();
	TH2F *h2 = new TH2F("h2","",300, -0.14, 0.14, 300, Y_Min*0.8, Y_Max*1.15);
	h2->SetStats(kFALSE);
	h2->SetXTitle("z_{react} (m)");
	h2->GetXaxis()->CenterTitle(1);
	h2->GetXaxis()->SetTitleFont(32);
	h2->GetXaxis()->SetTitleSize(0.08);
	h2->GetXaxis()->SetTitleOffset(0.7);
	h2->SetYTitle("Y_{0}");
	h2->GetYaxis()->CenterTitle(1);
	h2->GetYaxis()->SetTitleFont(32);
	h2->GetYaxis()->SetTitleSize(0.07);
	h2->GetYaxis()->SetTitleOffset(0.8);
	h2->Draw();

	TGraphErrors *fxL = new TGraphErrors(bin, VZL, Y0L, Zero, Y0L_Err);
	fxL->SetMarkerStyle(20);
	fxL->SetMarkerColor(4);
	fxL->SetMarkerSize(1.2);
	fxL->SetLineColor(4);
	fxL->Draw("P");
	TGraphErrors *fxR = new TGraphErrors(bin, VZR, Y0R, Zero, Y0R_Err);
	fxR->SetMarkerStyle(25);
	fxR->SetMarkerSize(1.2);
	fxR->SetMarkerColor(4);
	fxR->SetLineColor(2);
	fxR->Draw("P");

	TLegend *l2 = new TLegend(0.35,0.6,0.6,0.75,Form("%s Boiling Study:",Target_Name.Data()));
	l2->SetTextFont(32);
	l2->SetTextSize(0.03);
	l2->SetShadowColor(0);
	l2->AddEntry(exL,"HRS-L","p");
	l2->AddEntry(exR,"HRS-R","p");
	l2->Draw();	

	/*}}}*/

	c1->Print(Form("%s_BF.png",Target_Name.Data()));
	c1->Print(Form("%s_BF.pdf",Target_Name.Data()));
	c2->Print(Form("%s_Const.png",Target_Name.Data()));
	c2->Print(Form("%s_Const.pdf",Target_Name.Data()));

}
