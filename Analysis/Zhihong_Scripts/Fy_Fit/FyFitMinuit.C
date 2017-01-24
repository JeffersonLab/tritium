///////////////2//////////////////////
//     F(y) Parameters Fitting      //
//     TMinuit Minimization         //
//    ---  Zhihong Ye 12/24/2012    //
//////////////////////////////////////
#include "FyMinuit.h"
//Maximium Value of Y to fit:
Double_t Ymax = 0.20,       Ymin = -2.0;
Double_t Ymax_XGT2 = 0.20, Ymin_XGT2 = -2.0;
//Double_t Ymax = 0.1; //For C12 only
//Only fit one angle 18,22,26,32,40,50, 0->for all:
const int Angle = 18; //E02019 Data
//Only fit one angle 21,23,25,28, 0->for all:
const int Angle_XGT2 = 0; 
Int_t N_XGT2 = 0,N_NF=0;
Int_t A,Z;
double yaxis_min = 1e-05;
const TString XS_Date = "./Feb09/";

/*Main{{{*/
Int_t main()
{ 
	gStyle->SetOptFit(1);  
	gStyle->SetOptStat(0);

	cerr<<"Target (H2,He3,He4,C12) = "; cin >> Target_Name;

	if(Target_Name=="H2"){
		Ymax = 0.01; 
		Ymin_XGT2 = -0.8;
	}
	if(Target_Name=="He3"){
		//Ymin_XGT2 = -0.85;
		Ymin_XGT2 = -1.20;
	}
	if(Target_Name=="He4"){
		//Ymin_XGT2 = -0.90;
		Ymin_XGT2 = -1.20;
	}
	if(Target_Name=="C12"||Target_Name=="Ca40"||Target_Name=="Ca48"){
		Ymin_XGT2 = -0.80;
	}
			if(Target_Name=="C12"){
		yaxis_min = 1e-04;
	}
	if(Target_Name=="H2"){
		yaxis_min = 1e-03;
	}

	/*Get Target Info{{{*/
	if(Target_Name == "H2"){
		A = 2; Z = 1;  
	}
	else if(Target_Name == "He3"){
		A = 3; Z = 2;  
	}
	else if(Target_Name == "He4"){
		A = 4; Z = 2; 
	}
	else if(Target_Name == "Be"){
		A = 9; Z = 4; 
	}
	else if(Target_Name == "C12"){
		A = 12; Z = 6; 
	}
	else if(Target_Name == "Al"){
		A = 27; Z = 13; 
	}
	else if(Target_Name == "Ca40"){
		A = 40; Z = 20;
	}
	else if(Target_Name == "Ca48"){
		A = 48; Z = 20; 
	}
	else if(Target_Name == "Cu"){
		A = 64; Z = 29;
	}
	else if(Target_Name == "Au"){
		A = 197; Z = 79; 
	}
	else
		cerr<<" **** ERROR, I don't know the target!!!"<<endl;
	/*}}}*/
	cerr<<Form("Initializing parameters for target %s, A=%d, Z=%d",
			Target_Name.Data(),A,Z)<<endl;
	Target->GetValueAZ(A,Z);
	Target->Print();

	LoadDataXGT2();
	//if(!(Target_Name.Contains("Ca")))
	//	LoadDataAscii();

	double fyPar_Old[5] = {Target->f0,Target->B, Target->alpha, Target->a, Target->b};
	double fyPar[5] = {Target->f0,Target->B, Target->alpha, Target->a, Target->b};

	int iChoice = 2; //1-> TMinuit, 2->Graphic Fitting
	cerr<<"Fitting Method (1-> TMinuit, 2->Graphic Fitting) = "; cin >> iChoice;

	/////////////////////////////
	//TMinuit Minimization
	/////////////////////////////
	if(iChoice==1)
		DoMinuit(fyPar);

	////////////////////////////
	//Graphic Fitting
	////////////////////////////
	if(iChoice==2)
		DoGraphFit(fyPar);
		//DoNkFit(fyPar);

	//DoRatio(fyPar_Old,fyPar,iChoice);
	DoRatio_Angle(fyPar_Old,fyPar,iChoice);

	//  output(fyPar,iChoice);
	return 0;
}
/*}}}*/

/*DoMunuit{{{*/
void DoMinuit(double *par)
{
	cout << "      *************************************************" << endl; 
	cout << "      *          Minimization Fitting for             *" << endl;
	cout << "      *              F(y) Function                    *" << endl;
	cout << "      *             Z. Ye 09/20/2010                  *" << endl;
	cout << "      *************************************************" << endl; 
	cout << endl;
	gSystem->Load("libMinuit");
	static const Int_t iParaNum=5;

	TMinuit *pMinuit = new TMinuit(iParaNum);  //initialize TMinuit with a maximum of iParaNum params
	pMinuit->SetFCN(myFcn);

	Double_t arglist[10];
	Int_t ierflg = 0;

	/// set print level
	/*  SET PRIntout  <level>
		Sets the print level, determining how much output will be
		produced. Allowed values and their meanings are displayed
		after a SHOw PRInt command, and are currently <level>=:
		[-1]  no output except from SHOW commands
		[0]  minimum output
		[1]  default value, normal output
		[2]  additional output giving intermediate results.
		[3]  maximum output, showing progress of minimizations. */
	arglist[0] = 1;
	pMinuit->mnexcm("SET PRIntout",arglist,1,ierflg);

	/*
	   SET NOWarnings
	   Supresses Minuit warning messages.
	   SET WARnings
	   Instructs Minuit to output warning messages when suspicious
	   conditions arise which may indicate unreliable results.
	   This is the default.
	   */
	arglist[0] = 1;
	pMinuit->mnexcm("SET NOWarnings",arglist,0,ierflg);

	//Set Name of Parameters
	TString namePar[iParaNum]={"f0","bigB","alpha","a","b"};

	// Double_t inPar[iParaNum]={Target->f0,   
	// 			    Target->B,   
	// 			    Target->alpha,
	// 			    Target->a,
	// 			    Target->b};	
	Double_t inPar[iParaNum]={par[0],par[1],par[2],par[3],par[4]};   

	//Set Stepsize of Parameters
	Double_t step[iParaNum]={ 0.00000001,0.00000001,0.00001,0.00000001,0.00000001};

	//Set Min of Parameters, value==0 for No-Bound
	// Double_t minVal[iParaNum]={ 0.0001,0.0001,100.0,0.0001,0.0001};
	Double_t minVal[iParaNum]={ 0.0,0.0,0.0,0.0,0.0};

	//Set Max of Parameters, value==0 for No-Bound
	// Double_t maxVal[iParaNum]={ 0.01,0.01,300.0,0.01,0.01};
	Double_t maxVal[iParaNum]={ 0.0,0.0,0.0,0.0,0.0};

	//Initializing Parameters
	for(int ii=0;ii<iParaNum;ii++)
	{
		pMinuit->DefineParameter(ii,namePar[ii], 
				inPar[ii], step[ii], minVal[ii], maxVal[ii]);
	}

	//Fix parameters
	//arglist[0] = 3;   pMinuit->mnexcm("FIX ", arglist ,1,ierflg);
	//arglist[0] = 3;   pMinuit->mnexcm("FIX ", arglist ,1,ierflg);

	arglist[0] = 1;
	pMinuit->mnexcm("SET ERR", arglist,1,ierflg);

	/*
	//SET LIMits  0->[parno]  1->[lolim]  2->[uplim]
	//	  arglist[0] = 0;	   //this means set all parameters with the same limit
	//	  arglist[0] = 1;	   //this means set 1st parameters with the specified limit
	//	  arglist[0] = 2;	   //this means set 2nd parameters with the sepecified limit
	arglist[0] = 0;
	arglist[1] = 0.;
	arglist[2] = 0.;
	pMinuit->mnexcm("SET LIMits",arglist,3,ierflg);
	*/

	//Set Iteration Times and Use Migrad to start minimization
	arglist[0] = 500;
	arglist[1] = 0.01;
	pMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

	//Get Results

	for(int ii = 0; ii < iParaNum; ii++)
	{
		pMinuit->GetParameter(ii,outPar[ii],err[ii]);    
		par[ii] = outPar[ii];
	}
	//Put the results into a file
	cerr<<"Before:"<<endl;
	//  DoPlot(inPar,"Before");
	output(inPar,1);

	cerr<<"After:"<<endl;
	// DoPlot(outPar,"After");
	output(outPar,1);

}
/*}}}*/

/*myFcn{{{*/
// Calculate F(y), then get the Chi-Square Difference between Expertment Data,
// TMinuite will find a set of parameters to make this Chi-Square minimized
// Chi_Square = (F(y)- F_exp)**2
/*
   Instead of calculating Fy, calculate XS directly and direclty compare with expeirment XS
   */
void myFcn(Int_t& npar, Double_t* deriv, Double_t& f, Double_t *par, Int_t flag)
{

	Double_t aDelta = 0.0, aSum=0.0;
	//Set Parameters
	Target->f0 = par[0]; 
	Target->B = par[1]; 
	Target->alpha = par[2]; 
	Target->a = par[3]; 
	Target->b = par[4]; 
	Double_t aFy_Fit = 0.0, aFy_EX = 0.0, aFy_EX_Err = 0.0, aY = 0.0;
	int aNN =0;
	aNN = vY.size();
	//	if(Target_Name=="H2")
	//		aNN = vY.size();
	//	else
	//	aNN=vY_XGT2.size();

	for ( unsigned int ii=0; ii< aNN; ii++ )
	{
		if(Target_Name=="H2"||Target_Name=="He3"||Target_Name=="He4"){
			//	if(0){
			aFy_EX = vFy[ii];
			aFy_EX_Err = vFyErr[ii];
			aY = vY[ii];
		}
		else{
			aFy_EX = vFy_XGT2[ii];
			aFy_EX_Err = vFyErr_XGT2[ii];
			aY = vY_XGT2[ii];
		}

		if(aY <= Ymax){
			aFy_Fit = gGet_Fy(aY, Target);
			aDelta = (aFy_EX - aFy_Fit) / aFy_EX_Err;
			aSum+=aDelta*aDelta;
		}
		}

		f =  sqrt(aSum)/vXbj.size();
		//f = aSum;

		//   cerr << " >>> Chisq = " << f <<endl;
	}
	/*}}}*/

//Scaling Function F(y), from N. Fomin's Thesis (143P)
/*DoGraphFit{{{*/
void DoGraphFit(double* par)
{
	int n=vY.size();

	if (yMin-0.1>yaxis_min)
		yaxis_min = yMin-0.1;

	TCanvas *c1=new TCanvas("c1","",1600,800);
	c1->Divide(2,1);  c1->cd(1); gPad->SetLogy(0);
	TH2F *h1 = new TH2F("h1",Form("%s: F(y) Fitting (Before)",Target_Name.Data()),1000, xMin-0.1,xMax+0.5,1000,yaxis_min,yMax+1.);
	h1->SetXTitle("y");
	h1->SetYTitle("F(y)");
	h1->Draw();
	TGraphErrors *gr=new TGraphErrors(n,&vY[0],&vFy[0],&vZero[0],&vFyErr[0]);
	gr->SetMarkerStyle(21);
	gr->SetMarkerColor(4);
	gr->Draw("P");

	TF1 *thisF=new TF1("thisF",func,xMin-0.001,xMax+0.001,5);
	thisF->SetParameters(par[0],par[1],par[2],par[3],par[4]);
	thisF->SetParNames("f0","B","alpha","a","b");
	thisF->SetLineColor(2);
	thisF->Draw("same");

	c1->cd(2);gPad->SetLogy(0);
	TH2F *h2 = new TH2F("h2",Form("%s: F(y) Fitting (After)",Target_Name.Data()),1000, xMin-0.1,xMax+0.5,1000,yaxis_min,yMax+1.);
	h2->SetXTitle("y");
	h2->SetYTitle("F(y)");
	h2->Draw();
	TGraphErrors *gr1=new TGraphErrors(n,&vY[0],&vFy[0],&vZero[0],&vFyErr[0]);
	gr1->SetMarkerStyle(21);
	gr1->SetMarkerColor(4);
	gr1->Draw("P");

	Double_t oldpar[5];
	for(int i=0;i<5;i++)
		oldpar[i] = par[i];

	TF1 *thisF1=new TF1("thisF1",func,xMin-0.0,Ymax,5);
	thisF1->SetParameters(par[0],par[1],par[2],par[3],par[4]);
	thisF1->FixParameter(2,par[2]);
	thisF1->SetParNames("f0","B","alpha","a","b");
	thisF1->SetLineColor(2);
	thisF1->Draw("same");
	//    gr1->Fit(thisF1,"","",xMin,xMax);
	gr1->Fit(thisF1,"","",xMin,Ymax_XGT2);

	for(int i=0;i<5;i++)
		par[i] = thisF1->GetParameter(i);

	cerr<<"Before:"<<endl;
	output(oldpar,2);
	cerr<<"After:"<<endl;
	output(par,2);

	c1->Print(Form("%s_%d_GraphFit.png",Target_Name.Data(),Angle));
	delete gr;
	delete gr1;
	delete h1;
	delete h2;
	delete thisF;
	delete thisF1;
	delete c1;
}
/*}}}*/

/*DoNkFit{{{*/
void DoNkFit(double* par)
{
	const int n=vY.size();
    double Fy_Log[n], FyErr_Log[n];

	double fy_max = -1000, fy_min = 1000;
	for(int i=0;i<n;i++){
	     Fy_Log[i] = log(vFy[i]);
		 FyErr_Log[i] = log(vFyErr[i]);

		 if(fy_min > Fy_Log[i])
			 fy_min = Fy_Log[i];
		 if(fy_max < Fy_Log[i])
			 fy_max = Fy_Log[i];
	}

	if (yMin-0.1>yaxis_min)
		yaxis_min = yMin-0.1;

	TCanvas *c1=new TCanvas("c1","",1000,800);
	c1->Divide(1,2);  c1->cd(1); gPad->SetLogy(0);
	TH2F *h1 = new TH2F("h1",Form("%s: F(y) Fitting (Before)",Target_Name.Data()),1000, xMin-0.1,xMax+0.1,1000,fy_min*0.9,fy_max*1.1);
	h1->SetXTitle("y");
	h1->SetYTitle("F(y)");
	h1->Draw();
	TGraphErrors *gr=new TGraphErrors(n,&vY[0],Fy_Log,&vZero[0],FyErr_Log);
	gr->SetMarkerStyle(21);
	gr->SetMarkerColor(4);
	gr->Draw("P");

    double k[n], nk[n], nk_err[n];
	nk[0]= -0.01; nk_err[0] = 0.0;
	double nk_max = -1000, nk_min = 1000;
	for(int i=1;i<n;i++){
		cerr<<Form(" ===== Y1=%8.4f, Y2=%8.4f, Fy1 = %e, Fy2 = %e", vY[i-1], vY[i], vFy[i-1], vFy[i])<<endl;
		gr->Fit("pol1","","",vY[i-1]*0.99,vY[i]*1.01);
        nk[i] = exp(gr->GetFunction("pol1")->GetParameter(1));  
        nk_err[i] = exp(gr->GetFunction("pol1")->GetParError(1)); 
	    k[i] = 0.5*(vY[i-1]+vY[i]);	

		if(nk_min > nk[i])
			nk_min = nk[i];
		if(nk_max < nk[i])
			nk_max = nk[i];
	}

	c1->cd(2);gPad->SetLogy(1);
	TH2F *h2 = new TH2F("h2",Form("%s: n(k) Fitting",Target_Name.Data()),1000,0.0,1.1,1000,nk_min*0.7,nk_max*1.5);
	h2->SetXTitle("k");
	h2->SetYTitle("n(k)");
	h2->Draw();
	TGraphErrors *gr1=new TGraphErrors(n-1,k,nk,&vZero[0],nk_err);
	gr1->SetMarkerStyle(21);
	gr1->SetMarkerColor(4);
	gr1->Draw("P");

	ofstream outfile(Form("%s_Fy_Nk.out",Target_Name.Data()));

	for(int i=0;i<n;i++)
            outfile <<Form("%8.6f  %12e %12e %8.6f  %12e  %12e", vY[i], vFy[i], vFyErr[i], k[i], nk[i], nk_err[i])<<endl;

	cerr<<"After:"<<endl;
	output(par,2);

	c1->Print(Form("%s_%d_GraphFit.png",Target_Name.Data(),Angle));
	delete gr;
	delete gr1;
	delete h1;
	delete h2;
	delete c1;
}
/*}}}*/

/*DoPlot{{{*/
void DoPlot(Double_t *par, const TString& kCom)
{
	int n=vY.size();

	if (yMin-0.1>yaxis_min)
		yaxis_min = yMin-0.1;

	TCanvas *c1=new TCanvas("c1","",1400,800);
	c1->cd(); gPad->SetLogy(1);
	TH2F *h1 = new TH2F("h1",Form("%s: F(y) Fitting (%s)",Target_Name.Data(),kCom.Data()),1000, xMin-0.1,xMax+0.1,1000,yaxis_min,yMax+10.);
	h1->SetXTitle("y");
	h1->SetYTitle("F(y)");
	h1->Draw();
	TGraphErrors *gr=new TGraphErrors(n,&vY[0],&vFy[0],&vZero[0],&vFyErr[0]);
	gr->SetMarkerStyle(21);
	gr->SetMarkerColor(4);
	gr->Draw("P");

	TF1 *thisF=new TF1("thisF",func,xMin-0.001,xMax+0.001,5);
	thisF->SetParameters(par[0],par[1],par[2],par[3],par[4]);
	// thisF->SetParameters(Target->f0,Target->B,Target->alpha,Target->a,Target->b);
	thisF->SetParNames("f0","B","alpha","a","b");
	thisF->SetLineColor(2);
	thisF->Draw("same");

	c1->Print(Form("%s_%d_%s.png",Target_Name.Data(), Angle, kCom.Data()));

	/*
	c1->cd(); gPad->SetLogy(0);
	TH2F *h2 = new TH2F("h2",Form("%s: y .vs. x_{bj} (%s)",Target_Name.Data(),kCom.Data()),1000, xbjMin-0.1,3.5,1000,xMin-0.1,xMax+0.1);
	h2->SetYTitle("x_{bj}");
	h2->SetXTitle("y");
	h2->Draw();
	TGraphErrors *grr=new TGraphErrors(n,&vY[0],&vXbj[0],&vZero[0],&vZero[0]);
	grr->SetMarkerStyle(21);
	grr->SetMarkerColor(4);
	grr->Draw("P");
	c1->Print(Form("%s_%d_%s_xy.png",Target_Name.Data(), Angle, kCom.Data()));
*/


	delete gr;  
	delete h1;
	delete thisF;
	delete c1;
}
/*}}}*/

/*DoRatio{{{*/
void DoRatio(double* par_old, double* par_new, const int aChoice){
	/*ROOT Define{{{*/
	// Fill color
	gStyle->SetStatColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetPadColor(0);
	gStyle->SetFrameFillColor(0);

	// Border mode
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetPalette(1,0);

	// Font
	gStyle->SetTextFont(32);
	gStyle->SetLabelFont(32, "xyz");
	gStyle->SetTitleFont(22, "xyz");
	gStyle->SetStatFont(32);
	gStyle->SetTitleFillColor(0);

	// Opt
	//    gStyle->SetOptDate(0);
	//    gStyle->SetOptTitle(1);
	//    gStyle->SetOptStat(1);
	//    gStyle->SetOptFit(1);
	//    // Margin
	// //    gStyle->SetPadLeftMargin(0.1);
	// //    gStyle->SetPadRightMargin(0.1);
	// //    gStyle->SetPadTopMargin(0.1);
	// //    gStyle->SetPadBottomMargin(0.1);
	//void HallA_style() {
	gROOT->SetStyle("Plain");
	gStyle->SetPaperSize(TStyle::kUSLetter);
	gStyle->SetPaperSize(18,22);
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	gStyle->SetNdivisions(505);

	gStyle->SetCanvasColor(0);
	gStyle->SetPadTopMargin(.05);
	gStyle->SetPadLeftMargin(.15);
	gStyle->SetPadRightMargin(.1);
	gStyle->SetPadBottomMargin(.15);
	gStyle->SetTitleYOffset(1.3);
	gStyle->SetLabelFont(32,"X");
	gStyle->SetLabelFont(32,"Y");

	// prepare gStyle to be useful
	//   1 = solid
	//   2 = long dash (30 10)
	gStyle->SetTitleFillColor(-1); 
	//   3 = dotted (4 8)
	//   4 = dash-dot (15 12 4 12)
	//   5 = short dash ( 15 15 )
	//   6 = dash-dot-dot 
	gStyle->SetLineStyleString(1,"[]");
	gStyle->SetLineStyleString(2,"[30 10]");
	gStyle->SetLineStyleString(3,"[4 8]");
	gStyle->SetLineStyleString(4,"[15 12 4 12]");
	gStyle->SetLineStyleString(5,"[15 15]");
	gStyle->SetLineStyleString(6,"[15 12 4 12 4 12]");
	gStyle->SetLabelSize(0.05,"X");
	gStyle->SetLabelSize(0.05,"Y");
	gStyle->SetNdivisions(505,"Y");
	gStyle->SetOptDate(0);
	gStyle->SetDateY(.98);
	gStyle->SetStripDecimals(kFALSE);

	gStyle->SetHistFillColor(0); 
	gStyle->SetFillColor(0); 
	gStyle->SetFillStyle(4000); 
	gStyle->SetStatStyle(0); 
	gStyle->SetTitleStyle(0); 
	gStyle->SetCanvasBorderSize(0); 
	gStyle->SetFrameBorderSize(0); 
	gStyle->SetLegendBorderSize(0); 
	gStyle->SetStatBorderSize(0); 
	gStyle->SetTitleBorderSize(0); 

	/*}}}*/

	int n = vY.size();
	TString aCom="";
	if(aChoice==1)
		aCom="Minuit";
	if(aChoice==2)
		aCom="Graph";


	TCanvas *c1 = new TCanvas("c1","C1 - Ratio", 1200,800);
	c1->cd(); gPad->SetLogy(1);
	TH2F *h1 = new TH2F("h1",Form("%s: F(y) Fitting (Red -> Old,Blue -> New)",Target_Name.Data()),1000, xMin-0.1,xMax+0.1,1000,yaxis_min,yMax+1.);
	h1->SetFillColor(0);
	h1->SetTitle("");
	h1->SetXTitle("y");
	h1->GetXaxis()->CenterTitle(1);
	h1->GetXaxis()->SetTitleFont(32);
	h1->GetXaxis()->SetTitleSize(0.05);
	h1->SetYTitle("F(y)");
	h1->GetYaxis()->CenterTitle(1);
	h1->GetYaxis()->SetTitleFont(32);
	h1->GetYaxis()->SetTitleSize(0.04);
	h1->Draw();

	TGraphErrors *gr=new TGraphErrors(N_NF,&vY_NF[0],&vFy_NF[0],&vZero[0],&vFyErr_NF[0]);
	gr->SetMarkerStyle(25);
	gr->SetMarkerColor(6);
	gr->Draw("P");

	TGraphErrors *gr1=new TGraphErrors(N_XGT2,&vY_XGT2[0],&vFy_XGT2[0],&vZero[0],&vFyErr_XGT2[0]);
	gr1->SetMarkerStyle(23);
	gr1->SetMarkerColor(1);
	gr1->Draw("P");

	TF1 *thisF1=new TF1("thisF1",func,xMin-0.001,xMax+0.001,5);
	thisF1->SetParameters(par_old[0],par_old[1],par_old[2],par_old[3],par_old[4]);
	thisF1->SetParNames("f0","B","alpha","a","b");
	thisF1->SetLineColor(2);
	if(!(Target_Name.Contains("Ca"))){
		thisF1->Draw("same");
	}
	// TH2F *h2 = new TH2F("h2",Form("%s: F(y) Fitting (After)",Target_Name.Data()),1000, xMin-0.1,xMax+0.1,1000,yaxis_min,yMax+1.);
	// h2->SetXTitle("y");
	// h2->SetYTitle("F(y)");
	// h2->Draw();

	TF1 *thisF2=new TF1("thisF2",func,xMin-0.001,xMax+0.001,5);
	thisF2->SetParameters(par_new[0],par_new[1],par_new[2],par_new[3],par_new[4]);
	thisF2->SetParNames("f0","B","alpha","a","b");
	thisF2->SetLineColor(4);
	thisF2->Draw("same");

	TLegend *t1 = new TLegend(0.55,0.2,0.9,0.45,Form("Target: %s",Target_Name.Data()));
	t1->SetDefaults();
	t1->SetBorderSize(0);
	t1->SetTextSize(0.03);
	t1->SetTextFont(32);
//	if(!(Target_Name.Contains("Ca")))
//		t1->AddEntry(gr,Form("E02-019 Data (%d^{0})",Angle),"p");
	t1->AddEntry(gr1,Form("E08-014 Data (all angles)"),"p");
	if(!(Target_Name.Contains("Ca")))
		t1->AddEntry(thisF1,Form("F(y) Fitting from E02-019"),"l");
	t1->AddEntry(thisF2,Form("F(y) Fitting from E08-014"),"l");
	t1->SetFillColor(3);
	t1->Draw();

	c1->Print(Form("%s_%d_%s.png",Target_Name.Data(), Angle, aCom.Data()));
	c1->Print(Form("%s_%d_%s.pdf",Target_Name.Data(), Angle, aCom.Data()));
}
/*}}}*/

/*DoRatio_Angle{{{*/
void DoRatio_Angle(double* par_old, double* par_new, const int aChoice){
	/*ROOT Define{{{*/
	// Fill color
	gStyle->SetStatColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetPadColor(0);
	gStyle->SetFrameFillColor(0);

	// Border mode
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetPalette(1,0);

	// Font
	gStyle->SetTextFont(32);
	gStyle->SetLabelFont(32, "xyz");
	gStyle->SetTitleFont(22, "xyz");
	gStyle->SetStatFont(32);
	gStyle->SetTitleFillColor(0);

	// Opt
	//    gStyle->SetOptDate(0);
	//    gStyle->SetOptTitle(1);
	//    gStyle->SetOptStat(1);
	//    gStyle->SetOptFit(1);
	//    // Margin
	// //    gStyle->SetPadLeftMargin(0.1);
	// //    gStyle->SetPadRightMargin(0.1);
	// //    gStyle->SetPadTopMargin(0.1);
	// //    gStyle->SetPadBottomMargin(0.1);
	//void HallA_style() {
	gROOT->SetStyle("Plain");
	gStyle->SetPaperSize(TStyle::kUSLetter);
	gStyle->SetPaperSize(18,22);
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	gStyle->SetNdivisions(505);

	gStyle->SetCanvasColor(0);
	gStyle->SetPadTopMargin(.05);
	gStyle->SetPadLeftMargin(.15);
	gStyle->SetPadRightMargin(.1);
	gStyle->SetPadBottomMargin(.15);
	gStyle->SetTitleYOffset(1.3);
	gStyle->SetLabelFont(32,"X");
	gStyle->SetLabelFont(32,"Y");

	// prepare gStyle to be useful
	//   1 = solid
	//   2 = long dash (30 10)
	gStyle->SetTitleFillColor(-1); 
	//   3 = dotted (4 8)
	//   4 = dash-dot (15 12 4 12)
	//   5 = short dash ( 15 15 )
	//   6 = dash-dot-dot 
	gStyle->SetLineStyleString(1,"[]");
	gStyle->SetLineStyleString(2,"[30 10]");
	gStyle->SetLineStyleString(3,"[4 8]");
	gStyle->SetLineStyleString(4,"[15 12 4 12]");
	gStyle->SetLineStyleString(5,"[15 15]");
	gStyle->SetLineStyleString(6,"[15 12 4 12 4 12]");
	gStyle->SetLabelSize(0.05,"X");
	gStyle->SetLabelSize(0.05,"Y");
	gStyle->SetNdivisions(505,"Y");
	gStyle->SetOptDate(0);
	gStyle->SetDateY(.98);
	gStyle->SetStripDecimals(kFALSE);

	gStyle->SetHistFillColor(0); 
	gStyle->SetFillColor(0); 
	gStyle->SetFillStyle(4000); 
	gStyle->SetStatStyle(0); 
	gStyle->SetTitleStyle(0); 
	gStyle->SetCanvasBorderSize(0); 
	gStyle->SetFrameBorderSize(0); 
	gStyle->SetLegendBorderSize(0); 
	gStyle->SetStatBorderSize(0); 
	gStyle->SetTitleBorderSize(0); 

	/*}}}*/

	const int n = vY.size();
	TString aCom="";
	if(aChoice==1)
		aCom="Minuit";
	if(aChoice==2)
		aCom="Graph";

    double Y_NF1[100], Fy_NF1[100], FyErr_NF1[100];
    double Y_NF2[100], Fy_NF2[100], FyErr_NF2[100];
	int N_NF1 = 0, N_NF2 =0;
	for(int i=0;i<N_NF;i++){
	   if(abs(vTheta_NF[i]-18)<0.5){
	       Y_NF1[N_NF1]  = vY_NF[i];
		   Fy_NF1[N_NF1] = vFy_NF[i];
		   FyErr_NF1[N_NF1]  = vFyErr_NF[i];
	       N_NF1++; 
	   }
	   if(abs(vTheta_NF[i]-22)<0.5){
	       Y_NF2[N_NF2]  = vY_NF[i];
		   Fy_NF2[N_NF2] = vFy_NF[i];
		   FyErr_NF2[N_NF2]  = vFyErr_NF[i];
	       N_NF2++; 
	   }
	}	

    double Y_21[100], Y_23[100], Y_25[100], Y_28[100];
    double Fy_21[100], Fy_23[100], Fy_25[100], Fy_28[100];
    double FyErr_21[100], FyErr_23[100], FyErr_25[100], FyErr_28[100];

	int N_21 = 0, N_23 =0, N_25=0, N_28 = 0;
	for(int i=0;i<N_XGT2;i++){
	   if(abs(vTheta[i]-21)<0.5){
	       Y_21[N_21]  = vY_XGT2[i];
		   Fy_21[N_21] = vFy_XGT2[i];
		   FyErr_21[N_21]  = vFyErr_XGT2[i];
	       N_21++; 
	   }
	   if(abs(vTheta[i]-23)<0.5){
	       Y_23[N_23]  = vY_XGT2[i];
		   Fy_23[N_23] = vFy_XGT2[i];
		   FyErr_23[N_23]  = vFyErr_XGT2[i];
	       N_23++; 
	   }
	   if(abs(vTheta[i]-25)<0.5){
	       Y_25[N_25]  = vY_XGT2[i];
		   Fy_25[N_25] = vFy_XGT2[i];
		   FyErr_25[N_25]  = vFyErr_XGT2[i];
	       N_25++; 
	   }
	   if(abs(vTheta[i]-28)<0.5){
	       Y_28[N_28]  = vY_XGT2[i];
		   Fy_28[N_28] = vFy_XGT2[i];
		   FyErr_28[N_28]  = vFyErr_XGT2[i];
	       N_28++; 
	   }
	}
	if(Target_Name.Contains("Ca"))
		yaxis_min =2*1e-4;

		double X_Min = -0.81;
	if(Target_Name.Contains("He"))
		X_Min = -1.01;
	if(Target_Name.Contains("H2"))
		X_Min = -0.91;
	if(Target_Name.Contains("C12")){
		yaxis_min =1.2*1e-5;
		X_Min = -0.91;
	}

	TCanvas *c1 = new TCanvas("c1","C1 - Ratio", 1200,800);
	c1->cd(); gPad->SetLogy(1);
	TH2F *h1 = new TH2F("h1",Form("%s: F(y) Fitting (Red -> Old,Blue -> New)",Target_Name.Data()),1000, X_Min,0.2,1000,yaxis_min,yMax+1.);
	h1->SetFillColor(0);
	h1->SetTitle("");
	h1->SetXTitle("y");
	h1->GetXaxis()->CenterTitle(1);
	h1->GetXaxis()->SetTitleFont(32);
	h1->GetXaxis()->SetTitleSize(0.08);
	h1->GetXaxis()->SetTitleSize(0.08);
	h1->GetXaxis()->SetNdivisions(10,2,10);
	h1->GetXaxis()->SetTitleOffset(0.7);
	h1->SetYTitle("F(y)");
	h1->GetYaxis()->CenterTitle(1);
	h1->GetYaxis()->SetTitleFont(32);
	h1->GetYaxis()->SetTitleSize(0.07);
	h1->GetYaxis()->SetTitleOffset(0.7);
	h1->Draw();

	TLegend *t1 = new TLegend(0.55,0.15,0.9,0.42,Form("Target: %s",Target_Name.Data()));
	t1->SetDefaults();
	t1->SetBorderSize(0);
	t1->SetTextSize(0.035);
	t1->SetTextFont(32);


	TGraphErrors *gr21=new TGraphErrors(N_21,&Y_21[0],&Fy_21[0],&vZero[0],&FyErr_21[0]);
	gr21->SetMarkerStyle(2);
	gr21->SetMarkerSize(1.5);
	gr21->SetMarkerColor(2);
	if(N_21>1){
		gr21->Draw("P");
		t1->AddEntry(gr21,Form("E08-014 Data (21^{o}, 1.3 GeV^{2})"),"p");
	}
	TGraphErrors *gr23=new TGraphErrors(N_23,&Y_23[0],&Fy_23[0],&vZero[0],&FyErr_23[0]);
	gr23->SetMarkerStyle(3);
	gr23->SetMarkerSize(1.5);
	gr23->SetMarkerColor(3);
	if(N_23>1){
		gr23->Draw("P");
		t1->AddEntry(gr23,Form("E08-014 Data (23^{o}, 1.6 GeV^{2})"),"p");
	}
	TGraphErrors *gr25=new TGraphErrors(N_25,&Y_25[0],&Fy_25[0],&vZero[0],&FyErr_25[0]);
	gr25->SetMarkerStyle(5);
	gr25->SetMarkerSize(1.5);
	gr25->SetMarkerColor(4);
	if(N_25>1){
		gr25->Draw("P");
		t1->AddEntry(gr25,Form("E08-014 Data (25^{o}, 1.7 GeV^{2})"),"p");
	}
	TGraphErrors *gr28=new TGraphErrors(N_28,&Y_28[0],&Fy_28[0],&vZero[0],&FyErr_28[0]);
	gr28->SetMarkerStyle(30);
	gr28->SetMarkerSize(1.5);
	gr28->SetMarkerColor(6);
	if(N_28>1){
		gr28->Draw("P");
		t1->AddEntry(gr28,Form("E08-014 Data (28^{o}, 2.2 GeV^{2})"),"p");
	}

	TGraphErrors *gr1=new TGraphErrors(N_NF1,&Y_NF1[0],&Fy_NF1[0],&vZero[0],&FyErr_NF1[0]);
	gr1->SetMarkerStyle(28);
	gr1->SetMarkerSize(1.5);
	gr1->SetMarkerColor(7);
    if(N_NF1>1){
		gr1->Draw("P");
		t1->AddEntry(gr1,Form("E02-019 Data (%d^{o}, 2.5 GeV^{2})",18),"p");
	}
	TGraphErrors *gr2=new TGraphErrors(N_NF2,&Y_NF2[0],&Fy_NF2[0],&vZero[0],&FyErr_NF2[0]);
	gr2->SetMarkerStyle(27);
	gr2->SetMarkerSize(1.5);
	gr2->SetMarkerColor(12);
    if(N_NF2>1){
		gr2->Draw("P");
		t1->AddEntry(gr2,Form("E02-019 Data (%d^{o}, 3.3 GeV^{2})",22),"p");
	}

	TF1 *thisF2=new TF1("thisF2",func,xMin-0.001,xMax+0.001,5);
	thisF2->SetParameters(par_new[0],par_new[1],par_new[2],par_new[3],par_new[4]);
	thisF2->SetParNames("f0","B","alpha","a","b");
	thisF2->SetLineColor(1);
	thisF2->SetLineStyle(9);
	thisF2->SetLineWidth(2.0);
	thisF2->Draw("same");
	t1->AddEntry(thisF2,Form("F(y) Fitting from E08-014"),"l");
	t1->SetFillColor(3);
	t1->Draw();

	TLatex* l2 = new TLatex();
	l2->SetNDC();
	l2->SetTextFont(32);
	l2->SetTextColor(4);
	l2->SetTextSize(0.03);
	l2->DrawLatex(0.2, 0.9,  Form("f0 = %10.4f (new) %10.4f (old)",     par_new[0]*1000, par_old[0]*1000 ));
	l2->DrawLatex(0.2, 0.87, Form("B  = %10.4f (new) %10.4f (old)",     par_new[1]*1000, par_old[1]*1000 ));
	l2->DrawLatex(0.2, 0.84, Form(" a  = %10.4f (new) %10.4f (old)",     par_new[3]*1000, par_old[3]*1000 ));
	l2->DrawLatex(0.2, 0.81, Form(" b  = %10.4f (new) %10.4f (old)",     par_new[4]*1000, par_old[4]*1000 ));
	l2->DrawLatex(0.2, 0.78, Form("#alpha = %10.4f (new) %10.4f (old)", par_new[2], par_old[2] ));

	c1->Print(Form("%s_XGT2_%s.png",Target_Name.Data(),  aCom.Data()));
	c1->Print(Form("%s_XGT2_%s.pdf",Target_Name.Data(),  aCom.Data()));
/*
	c1->cd(); gPad->SetLogy(0);
	//TH2F *h2 = new TH2F("h2",Form("%s: y .vs. x_{bj} (%s)",Target_Name.Data(),aCom.Data()),1000, xbjMin-0.1,xbjMax+0.1,1000,xMin-0.1,xMax+0.1);
	TH2F *h2 = new TH2F("h2",Form("%s: y .vs. x_{bj} (%s)",Target_Name.Data(),aCom.Data()),1000, xbjMin-0.1,3.5,1000,xMin-0.1,xMax+0.1);
	h2->SetXTitle("x_{bj}");
	h2->SetYTitle("y");
	h2->Draw();
	TGraphErrors *grr=new TGraphErrors(n,&vXbj[0],&vY[0],&vZero[0],&vZero[0]);
	grr->SetMarkerStyle(21);
	grr->SetMarkerColor(4);
	grr->Draw("P");


	c1->Print(Form("%s_%s_xy.png",Target_Name.Data(), aCom.Data()));
*/
}
/*}}}*/

/*func{{{*/
Double_t func(Double_t *x, Double_t *par)
{
	//Set Parameters
	XEM_TGT* kTarget = new XEM_TGT();
	kTarget->GetValueAZ(Target->A,Target->Z);
	Double_t kY = x[0];
	kTarget->f0 = par[0]; 
	kTarget->B = par[1]; 
	kTarget->alpha = par[2]; 
	kTarget->a = par[3]; 
	kTarget->b = par[4]; 

	double  kFy_Fit = gGet_Fy(kY, kTarget);
	delete kTarget;
	return kFy_Fit;
}
/*}}}*/

/*LoadData{{{*/
void LoadData(){
	Double_t xbj,Ep,Theta, xs_rad, xs_born,xs_dis,xs_qe, nXS_Rad, nXS_Born,nXS_CC,nFrac_Err;
	Int_t A,Z;

	ofstream outfile(Form("Fy_Table_%s.dat",Target_Name.Data()));
	outfile << Form("%4s %4s %6s %10s %10s %16s  %16s %10s %16s %16s", "A", "Z", "Theta","Ep", "xbj", "XS_QE", "XS_QE_Err", "Y", "Fy", "Fy_Err" )
		<<endl;

	TString Root_File = Form("%s_XEMC.root",Target_Name.Data());
	TChain *T = new TChain("T");
	T->Add(Root_File);
	int nTotal = T->GetEntries();
	T->SetBranchAddress("xbj", &xbj);
	T->SetBranchAddress("Ep",  &Ep);
	T->SetBranchAddress("Theta", &Theta);
	T->SetBranchAddress("xs_rad", &xs_rad);
	T->SetBranchAddress("xs_born",&xs_born);
	T->SetBranchAddress("xs_dis", &xs_dis);
	T->SetBranchAddress("xs_qe", &xs_qe);
	T->SetBranchAddress("nXS_Rad", &nXS_Rad);
	T->SetBranchAddress("nXS_Born", &nXS_Born);
	T->SetBranchAddress("nXS_CC", &nXS_CC);
	T->SetBranchAddress("nFrac_Err", &nFrac_Err);
	T->SetBranchAddress("A", &A);
	T->SetBranchAddress("Z", &Z);

	Double_t nXS_QE;
	XEM_VAR2* aFy = new XEM_VAR2();
	Double_t  aFy_EX = 0.0, aFy_EX_Err = 0.0, aY = 0.0;

	for(int i=0;i<nTotal; i++){
		T->GetEntry(i);
		if(i==0){
			cerr<<Form("Initializing parameters for target %s, A=%d, Z=%d",Target_Name.Data(),A,Z)<<endl;
			Target->GetValueAZ(A,Z);
			Target->Print();
		}

		if(abs(Angle-Theta)<1e-5 || (Angle-0.0)<1e-5){
			//  	   	    if(1){
			//F(y)
			nXS_QE = (nXS_Born - xs_dis);
			cerr<<Form("---: XS_Born = %e, XS_DIS = %e", nXS_Born, xs_dis)<<endl;
			//nXS_QE = xs_qe;
			//nXS_QE = (xs_born - xs_dis);
			aFy = gCal_Sig2Fy(E0, Ep, Theta * pi /180., nXS_QE, Target);
			aY = aFy->First;
			aFy_EX = aFy->Second;
			aFy_EX_Err =  aFy_EX * nFrac_Err;

			//    if(aY>-1e4 &&aY<=Ymax+0.01 ){
			//    if(aY>-1e4){
			if(aY>-1e4 && aY>Ymin){
				vEp.push_back(Ep);
				vTheta.push_back(Theta);
				vXbj.push_back(xbj);
				vZero.push_back(0.0);

				vXS.push_back(nXS_QE);
				vXSErr.push_back(nXS_QE*nFrac_Err);

				vFy.push_back(aFy_EX);
				vFyErr.push_back(aFy_EX_Err);
				vY.push_back(aY);
				outfile << Form("%4d %4d %6.2f %10.4f %10.4f %16.6e  %16.6e %10.4f %16.6e %16.6e", A, Z, Theta, Ep, xbj, nXS_QE, nXS_QE*nFrac_Err, aY, aFy_EX, aFy_EX_Err )
					<<endl;


				if(xbjMin >= xbj)
					xbjMin = xbj;
				if(xbjMax <=xbj)
					xbjMax = xbj;

				if(xMin >= aY)
					xMin = aY;
				if(xMax <=aY)
					xMax = aY;
				if(yMin >= aFy_EX)
					yMin = aFy_EX;
				if(yMax <= aFy_EX)
					yMax = aFy_EX;
			}
		}
		}
		yMin = 1e-30;

		delete aFy;
		T->Delete();
		}
		/*}}}*/

/*LoadDataAscii{{{*/
void LoadDataAscii(){
	TString Ascii_File = Form("data/%s.dat",Target_Name.Data());
	ifstream infile(Ascii_File);

	Double_t nA, nTheta, nP0, nEp, nXS_Rad, nXS_Born,nXS_CC,nFrac_Err, nXS_QE;  
	Double_t xs_dis = -1.,xbj,Q2;
	XEM_VAR2* aFy = new XEM_VAR2();
	Double_t  aFy_EX = 0.0, aFy_EX_Err = 0.0, aY = 0.0;
	int i=0;
	while(!(infile.eof())){
		infile >> nA >> nTheta >> nP0 >> nEp >> nXS_Rad >> nXS_Born >> nXS_CC >> nFrac_Err;

		Q2 = 4.0 * E0 * nEp * pow(sin(nTheta*3.1415926/180./2.0),2);
		xbj = Q2/2.0/0.938272/(E0-nEp);

		xs_dis = XEMC_Born(E0*GeVToMeV, nEp*GeVToMeV, nTheta, A, Z, 3); //3->For DIS_F1F2IN09
		//xs_dis = XEMC_Born(E0*GeVToMeV, nEp*GeVToMeV, nTheta, A, Z, 4); //4->For DIS_F1F2IN06
		// double xs_qe = XEMC_Born(E0*GeVToMeV, nEp*GeVToMeV, nTheta, A, Z, 2); //2->For QE_XEM  
		// double xs_born = XEMC_Born(E0*GeVToMeV, nEp*GeVToMeV, nTheta, A, Z, 1); //1->For DIS_F1F2IN09 + QE
		// nXS_QE = xs_qe;
		// nXS_QE = xs_born - xs_dis;
		nXS_QE = nXS_Born - xs_dis;
		
		
		//if(nTheta>(Angle+1))
		if(nTheta>(22+1))
			break;
		//if(abs(Angle-nTheta)<1e-5 || (Angle-0.0)<1e-5){
		if(abs(nTheta-18)<1e-5 || (nTheta-22)<1e-5){
		//if(abs(nTheta-18)<1e-5){
			//  	   	    if(1){
			//	cerr<<Form("E0=%f, Ep=%f, Theta=%f",E0,nEp,nTheta)<<endl;
			//F(y)
			aFy = gCal_Sig2Fy(E0, nEp, nTheta * pi /180., nXS_QE, Target);
			aY = aFy->First;
			aFy_EX = aFy->Second;
			aFy_EX_Err =  aFy_EX * nFrac_Err;

			//    if(aY>-1e4 &&aY<=Ymax+0.01 ){
			if(aY<0.15 && aY>=-0.8){
				vEp.push_back(nEp);
				vTheta.push_back(nTheta);
				vXbj.push_back(xbj);
				vZero.push_back(0.0);
				nXS_QE = (nXS_Born - xs_dis);
				vXS.push_back(nXS_QE);
				vXSErr.push_back(nXS_QE*nFrac_Err);

				vFy.push_back(aFy_EX);
				vFyErr.push_back(aFy_EX_Err);
				vY.push_back(aY);

				vFy_NF.push_back(aFy_EX);
				vFyErr_NF.push_back(aFy_EX_Err);
				vY_NF.push_back(aY);
				vTheta_NF.push_back(nTheta);

				if(xbjMin >= xbj)
					xbjMin = xbj;
				if(xbjMax <=xbj)
					xbjMax = xbj;

				if(xMin >= aY)
					xMin = aY;
				if(xMax <=aY)
					xMax = aY;
				if(yMin >= aFy_EX)
					yMin = aFy_EX;
				if(yMax <= aFy_EX)
					yMax = aFy_EX;
			}
		}
		i++;
		}
		N_NF = vY_NF.size();
		delete aFy;
		}
		/*}}}*/

/*LoadDataXGT2{{{*/
void LoadDataXGT2(){
	//TString Ascii_File = Form("XS_%s.dat",Target_Name.Data());
	TString Ascii_File = Form("./%s/XS_%s_L.dat",XS_Date.Data(), Target_Name.Data());
    //cerr<<"==== Reading XS file from "<<Ascii_File.Data()<<endl;
	ifstream infile(Ascii_File);

	TString nTarget_Name, nCom; 
	Double_t nTheta, nE0, nEp, nXbj, nNu, nXS_Born, nXS_Born_Err,nXS_MC,nXS_QE, nFrac_Err;  
	Double_t xs_dis = -1.,xbj;
	XEM_VAR2* aFy = new XEM_VAR2();
	Double_t  aFy_EX = 0.0, aFy_EX_Err = 0.0, aY = 0.0;
	int i=0;
	infile >> nCom >> nCom >> nCom >> nCom >> nCom >> nCom >> nCom >> nCom >> nCom;
	while(!(infile.eof())){
		infile >> nTarget_Name >> nTheta >> nE0 >> nEp >> nNu >> nXbj >> nXS_Born >> nXS_Born_Err >> nXS_MC;

		if(nTheta>27) continue;

		if(abs(Angle_XGT2-nTheta)<1e-5 || (Angle_XGT2-0.0)<1e-5){
			cerr<<Form("--- %s: E0=%f, Ep=%f,Theta=%f, xbj=%f, XS=%e",nTarget_Name.Data(), nE0,nEp, nTheta, nXbj, nXS_Born)<<endl;
			xbj = nXbj;
			nFrac_Err = nXS_Born_Err/nXS_Born;

			xs_dis = XEMC_Born(nE0*GeVToMeV, nEp*GeVToMeV, nTheta, A, Z, 3); //3->For DIS_F1F2IN09
			//xs_dis = XEMC_Born(nE0*GeVToMeV, nEp*GeVToMeV, nTheta, A, Z, 4); //4->For DIS_F1F2IN06
			//   xs_qe = XEMC_Born(nE0, nEp, nTheta, A, Z, 2); //2->For QE_XEM  
			//   xs_born = XEMC_Born(nE0, nEp, nTheta, A, Z, 1); //1->For DIS_F1F2IN09 + QE
			nXS_QE = nXS_Born - xs_dis;

			//		cerr<<Form("--- %s: XS_Born = %e, XS_DIS = %e", nTarget_Name.Data(), nXS_Born, xs_dis)<<endl;

			//F(y)
			aFy = gCal_Sig2Fy(nE0, nEp, nTheta * pi /180., nXS_QE, Target);
			aY = aFy->First;
			aFy_EX = aFy->Second;
			aFy_EX_Err =  aFy_EX * nFrac_Err;

			//if(aY>Ymin_XGT2 &&aY<=Ymax+0.01 ){
			if(aY>Ymin_XGT2 ){
			//	if(aY>-1e4){
				cerr<<Form("--- %s: E0=%f, Ep=%f,Theta=%f, xbj=%f, y=%f",nTarget_Name.Data(), nE0,nEp, nTheta, nXbj, aY)<<endl;

				vEp.push_back(nEp);
				vTheta.push_back(nTheta);
				vXbj.push_back(xbj);
				vZero.push_back(0.0);
				vXS.push_back(nXS_QE);
				vXSErr.push_back(nXS_QE*nFrac_Err);

				vFy.push_back(aFy_EX);
				vFyErr.push_back(aFy_EX_Err);
				vY.push_back(aY);
		
				vTheta_XGT2.push_back(nTheta);
				vFy_XGT2.push_back(aFy_EX);
				vFyErr_XGT2.push_back(aFy_EX_Err);
				vY_XGT2.push_back(aY);

				if(xbjMin >= xbj)
					xbjMin = xbj;
				if(xbjMax <=xbj)
					xbjMax = xbj;
	
				if(xMin >= aY)
					xMin = aY;
				if(xMax <=aY)
					xMax = aY;
				if(yMin >= aFy_EX)
					yMin = aFy_EX;
				if(yMax <= aFy_EX)
					yMax = aFy_EX;
			}
			}
			i++;
		}
		N_XGT2 = vY_XGT2.size();
		delete aFy;
		}
		/*}}}*/

/*output{{{*/
//Output the results into a file
void output(Double_t *outPar,const int aChoice)
{
	TString aCom="";
	if(aChoice==1)
		aCom="Minuit";
	if(aChoice==2)
		aCom="Graph";

	TString name =  Form("%s_%d_Par_%s.dat",Target_Name.Data(),Angle,aCom.Data());

	ofstream out(name,ios::out);

	cerr<<endl<<"The Fitting Result is: "<<endl;

	cerr <<setiosflags(ios::fixed)
		<<"f0    = "
		<<setw(14)<<setprecision(8)<<outPar[0]<<endl;
	cerr <<setiosflags(ios::fixed)
		<<"bigB  = "
		<<setw(14)<<setprecision(8)<<outPar[1]<<endl;
	cerr <<setiosflags(ios::fixed)
		<<"alpha = "
		<<setw(14)<<setprecision(8)<<outPar[2]<<endl;

	cerr <<setiosflags(ios::fixed)
		<<"a     = "
		<<setw(14)<<setprecision(8)<<outPar[3]<<endl;
	cerr <<setiosflags(ios::fixed)
		<<"b     = "
		<<setw(14)<<setprecision(8)<<outPar[4]<<endl;

	out <<"======================================================================================="<<endl;
	out <<setiosflags(ios::fixed);
	out <<"   F(y) = (f0-bigB)*(alpha*alpha*Exp(-a*a*y*y)/(alpha*alpha+y*y)) + bigB*Exp(-b*|y|) "<<endl;
	out <<"----------------------------------------------------------------------------------------"<<endl;
	out <<" Where for Angle =  "<<Angle<<": "<<endl;
	out <<setiosflags(ios::fixed)
		<<"         f0    = "
		<<setw(14)<<setprecision(8)<<outPar[0]<<endl;
	out <<setiosflags(ios::fixed)
		<<"         bigB  = "
		<<setw(14)<<setprecision(8)<<outPar[1]<<endl;
	out <<setiosflags(ios::fixed)
		<<"         alpha = "
		<<setw(14)<<setprecision(8)<<outPar[2]<<endl;

	out <<setiosflags(ios::fixed)
		<<"         a     = "
		<<setw(14)<<setprecision(8)<<outPar[3]<<endl;
	out <<setiosflags(ios::fixed)
		<<"         b     = "
		<<setw(14)<<setprecision(8)<<outPar[4]<<endl;
	out <<" Note: Some values have to be scaled: alpha = alpha/1000.0, y = y*1000.0."<<endl;
	out <<"========================================================================================"<<endl;

	cerr << Form("%s %10.4f %10.4f %10.4f %10.4f %10.4f",Target_Name.Data(), outPar[0]*1000., outPar[1]*1000, outPar[3]*1000, outPar[4]*1000, outPar[2])<<endl;
	out<<endl << Form("%s %10.4f %10.4f %10.4f %10.4f %10.4f",Target_Name.Data(), outPar[0]*1000., outPar[1]*1000, outPar[3]*1000, outPar[4]*1000, outPar[2])<<endl;
}
/*}}}*/
