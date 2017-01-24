//void plot(){
{
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

	TString Target_Name; cerr<<"---Target (H2,He3,He4) = "; cin >> Target_Name;
	TString Kin_Name;    cerr<<"---Kin (3.1,3.2, ...) = ";  cin >> Kin_Name;
	Kin_Name = "Kin"+Kin_Name;
	ifstream infile(Form("%s_Yield_L_%s_highx.dat",Target_Name.Data(),Kin_Name.Data()));
	ifstream infileD(Form("Dummy20cm_Yield_L_%s_highx.dat", Kin_Name.Data()));

	/*Define{{{*/	
	const int bin=200;
	double Zero[bin], VZ[bin], VZD[bin],D[bin],Y[bin],S[bin], DE[bin],YE[bin], SE[bin];
	double Dp[bin], DEp[bin];
	double Ep[bin],EpE[bin], EpD[bin],EpED[bin];
	//For the Thickness difference of the Dummy target and the Cell
	const double Cor_Down=0.1863;
	const double Cor_Up=0.2448;
	double Y_Max = -19988, Y_Min = 19999;
	double Y_Max_Up = -19988, Y_Max_Down = 19999;
	double D_Max_Up = -19988, D_Max_Down = 19999;
	double S_Max = -19988, S_Min = 19999;
	double offset = 0.00;
	if(Target_Name=="H2")
		offset = 0.012;
	else if(Target_Name=="He3")
		offset = 0.016;
	else if(Target_Name=="He4")
		offset = 0.016;
	else
		cerr<<"*** ERROR ***"<<endl;
	/*}}}*/

	/*Read In{{{*/
	int i=0;
	while(!(infile.eof())){
		infile >> VZ[i] >> Y[i] >> YE[i];
		infileD >> VZD[i] >> Dp[i] >> DEp[i];

		VZ[i] += 0.01;
		VZD[i] += 0.01;
		 
		if(VZ[i]< 1e-16){
//			Dp[i] *= Cor_Up;
//			DEp[i] *= Cor_Up;
		if(Y_Max_Up<Y[i])
			Y_Max_Up = Y[i];
		if(D_Max_Up<Dp[i])
			D_Max_Up = Dp[i];
		}
		else{
//			Dp[i] *= Cor_Down;
//			DEp[i] *= Cor_Down;
			if(Y_Max_Down<Y[i])
				Y_Max_Down = Y[i];
			if(D_Max_Down<Dp[i])
				D_Max_Down = Dp[i];
		}

		if(Y_Max<Y[i])
			Y_Max = Y[i];
		if(Y_Min>Y[i])
			Y_Min = Y[i];
		i++;
	}
	const int N = i;

	for(int i=0;i<N;i++){
		if(VZ[i]< 1e-16){
			Dp[i] *= Y_Max_Up/D_Max_Up;
			DEp[i] *= Y_Max_Up/D_Max_Up;
		}
		else{
			Dp[i] *= Y_Max_Down/D_Max_Down;
			DEp[i] *= Y_Max_Down/D_Max_Down;
		}
		if(Y_Max<Dp[i])
			Y_Max = Dp[i];
		if(Y_Min>Dp[i])
			Y_Min = Dp[i];
	
	}
	/*}}}*/
	
	/*Sub{{{*/
	int Nk=0;
	if(Target_Name=="He3")
		Nk=4;
	else
		Nk=4;
	const int No = Nk;
	for(int j=No;j<N;j++){
	    D[j]=Dp[j-No];
	    DE[j]=DEp[j-No];
		S[j] = (Y[j] - D[j])/Y[j]*100.0;
		//SE[j] = sqrt(pow(YE[j],2) + pow(DE[j],2))/Y[j]*100.0;
		SE[j] = S[j]*YE[j]/Y[j];
		//SE[j] = S[j]*sqrt(pow(YE[j]/Y[j],2) + pow(DE[j]/D[j],2));
		Zero[j]=0.0;

		if(S_Max<S[j])
			S_Max = S[j];
		if(S_Min>S[j])
			S_Min = S[j];
	}
	/*}}}*/
	
	/*Plot Yield{{{*/
	gROOT->SetStyle("Plain");
	gStyle->SetFillColor(0);
	cerr<<" @@@ Making Plots Yield ..."<<endl;
	TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
    c1->Divide(1,2);
	c1->cd(1);
	gPad->SetLogy(1);
	TH2F *h1 = new TH2F("h1","",300, -0.14, 0.14, 300, Y_Max*1e-2, Y_Max*1.5);
	h1->SetStats(kFALSE);
	h1->SetXTitle("z_{react} (m)");
	h1->GetXaxis()->CenterTitle(1);
	h1->GetXaxis()->SetTitleFont(32);
	h1->GetXaxis()->SetTitleSize(0.08);
	h1->GetXaxis()->SetTitleOffset(0.7);
	h1->GetXaxis()->SetNdivisions(10,2,5);
	h1->SetYTitle("Y");
	h1->GetYaxis()->CenterTitle(1);
	h1->GetYaxis()->SetTitleFont(32);
	h1->GetYaxis()->SetTitleSize(0.07);
	h1->GetYaxis()->SetTitleOffset(0.8);
	h1->Draw();

	TGraphErrors *ex = new TGraphErrors(N, VZ, Y, Zero, YE);
	ex->SetMarkerStyle(20);
	ex->SetMarkerColor(2);
	ex->SetMarkerSize(1.2);
	ex->SetLineColor(4);
	ex->Draw("P");

	TGraphErrors *exD = new TGraphErrors(N, VZD, D, Zero, DE);
	exD->SetMarkerStyle(25);
	exD->SetMarkerSize(1.2);
	exD->SetMarkerColor(4);
	exD->SetLineColor(2);
	exD->Draw("P");

	TLegend *l1 = new TLegend(0.35,0.2,0.7,0.35,Form("Yield at %s on HRS-L:",Kin_Name.Data()));
	l1->SetTextFont(32);
	l1->SetTextSize(0.03);
	l1->SetShadowColor(0);
	l1->AddEntry(ex,Form("%s",Target_Name.Data()),"p");
	l1->AddEntry(exD,"Dummy-20cm (normalized)","p");
	l1->Draw();	

	/*}}}*/

	/*Plot Yield Substruction{{{*/
	gROOT->SetStyle("Plain");
	gStyle->SetFillColor(0);
	cerr<<" @@@ Making Plots Yield Substruction ..."<<endl;
//	TCanvas *c2 = new TCanvas("c2","c2",1200,700);
//  c2->cd();
	c1->cd(2);
	gPad->SetLogy(0);
	//TH2F *h2 = new TH2F("h2","",300, -0.14, 0.14, 300, 10.0, S_Max*1.15);
	TH2F *h2 = new TH2F("h2","",300, -0.14, 0.14, 300, 5.0, 101.0);
	h2->SetStats(kFALSE);
	h2->SetXTitle("z_{react} (m)");
	h2->GetXaxis()->CenterTitle(1);
	h2->GetXaxis()->SetTitleFont(32);
	h2->GetXaxis()->SetTitleSize(0.08);
	h2->GetXaxis()->SetTitleOffset(0.7);
	h2->GetXaxis()->SetNdivisions(10,1,5);
	h2->SetYTitle("#Delta Y(\%)");
	h2->GetYaxis()->CenterTitle(1);
	h2->GetYaxis()->SetTitleFont(32);
	h2->GetYaxis()->SetTitleSize(0.07);
	h2->GetYaxis()->SetTitleOffset(0.8);
	h2->GetYaxis()->SetNdivisions(10,2,5);
	h2->Draw();

	TGraphErrors *fx = new TGraphErrors(N, VZ, S, Zero, SE);
	fx->SetMarkerStyle(20);
	fx->SetMarkerColor(1);
	fx->SetMarkerSize(1.2);
	fx->SetLineColor(4);
	fx->Draw("P");

	TLegend *l2 = new TLegend(0.35,0.2,0.6,0.35,Form("Yield Substruction at %s:",Kin_Name.Data()));
	l2->SetTextFont(32);
	l2->SetTextSize(0.03);
	l2->SetShadowColor(0);
	l2->AddEntry(fx,Form("%s on HRS-L",Target_Name.Data()),"p");
	l2->Draw();	

	/*}}}*/
	
	/*Plot Yield Substruction -- Zoom{{{*/
/*	gROOT->SetStyle("Plain");
	gStyle->SetFillColor(0);
	cerr<<" @@@ Making Plots Yield Substruction ..."<<endl;
	TCanvas *c3 = new TCanvas("c3","c3",1200,700);
	c3->cd();
	gPad->SetLogy(0);
	TH2F *h3 = new TH2F("h3","",300, -0.08, 0.08, 300, 90.0, 100.5);
	h3->SetStats(kFALSE);
	h3->SetXTitle("z_{react} (m)");
	h3->GetXaxis()->CenterTitle(1);
	h3->GetXaxis()->SetTitleFont(32);
	h3->GetXaxis()->SetTitleSize(0.08);
	h3->GetXaxis()->SetTitleOffset(0.7);
	h3->GetXaxis()->SetNdivisions(10);
	h3->SetYTitle("#Delta Y(\%)");
	h3->GetYaxis()->CenterTitle(1);
	h3->GetYaxis()->SetTitleFont(32);
	h3->GetYaxis()->SetTitleSize(0.07);
	h3->GetYaxis()->SetTitleOffset(0.8);
	h3->GetYaxis()->SetNdivisions(10,2,5);
	h3->Draw();
    fx->Draw("P");
	l2->Draw();	
*/
	/*}}}*/

	c1->Print(Form("%s_Yield_%s_highx.png",Target_Name.Data(),Kin_Name.Data()));
	c1->Print(Form("%s_Yield_%s_highx.pdf",Target_Name.Data(),Kin_Name.Data()));
//	c2->Print(Form("%s_Sub.png",Target_Name.Data()));
	//c3->Print(Form("%s_Sub_Zoom.png",Target_Name.Data()));

}
