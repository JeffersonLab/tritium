//////////////////////////////////////
//     F(y) Parameters Fitting      //
//     TMinuit Minimization         //
//    ---  Zhihong Ye 12/24/2012    //
//////////////////////////////////////
/*Headers{{{*/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include <string>
#include <cstdio>

#include "math.h"
#include "stdio.h"
#include <stdlib.h>

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TF1.h>
#include <TH1.h>
#include <TAxis.h>
#include <TLine.h>
#include <TLegend.h>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TObject.h"
#include "TSystem.h"
#include <TVirtualFitter.h>
#include <TMinuit.h>
//#include <TFitterMinuit.h>
#include <TApplication.h>
#include <vector>
using namespace std;

/*}}}*/

/*Main{{{*/
Int_t main()
{ 
	gStyle->SetOptFit(1);  
	gStyle->SetOptStat(0);

	double fyPar_Old[5] = {Target->f0,Target->B, Target->alpha, Target->a, Target->b};
	double fyPar[5] = {Target->f0,Target->B, Target->alpha, Target->a, Target->b};
	
	LoadDataXGT2();

	/////////////////////////////
	//TMinuit Minimization
	/////////////////////////////
	DoMinuit(fyPar);

	output(fyPar,iChoice);
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

	//	LoadData();
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
	if(Target_Name=="H2")
		aNN = vY.size();
	else
		aNN=vY_XGT2.size();

	for ( unsigned int ii=0; ii< aNN; ii++ )
	{
		if(Target_Name=="H2"||Target_Name=="He3"){
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

/*LoadDataXGT2{{{*/
void LoadDataXGT2(){
	//TString Ascii_File = Form("XS_%s.dat",Target_Name.Data());
	TString Ascii_File = Form("./%s/XS_%s.dat",XS_Date.Data(), Target_Name.Data());
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
        cerr<<Form("--- %s: E0=%f, Ep=%f,Theta=%f, xbj=%f, XS=%e",nTarget_Name.Data(), nE0,nEp, nTheta, nXbj, nXS_Born)<<endl;

		if(abs(Angle_XGT2-nTheta)<1e-5 || (Angle_XGT2-0.0)<1e-5){
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
				//if(aY>-1e4){
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
